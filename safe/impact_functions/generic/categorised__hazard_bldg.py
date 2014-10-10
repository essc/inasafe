# coding=utf-8
"""Impact of flood on roads."""
from qgis.core import (
    QgsRectangle,
    QgsFeatureRequest,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform
)

from safe.metadata import (
    hazard_all,
    unit_categorised,
    layer_raster_numeric,
    exposure_structure,
    unit_building_type_type,
    unit_building_generic,
    layer_vector_polygon,
    hazard_definition,
    exposure_definition
)
from safe.common.utilities import OrderedDict
from safe.impact_functions.core import FunctionProvider
from safe.impact_functions.core import get_hazard_layer, get_exposure_layer
from safe.impact_functions.core import get_question
from safe.common.tables import Table, TableRow
from safe.common.utilities import ugettext as tr
from safe.impact_functions.impact_function_metadata import (
    ImpactFunctionMetadata)
from safe.storage.vector import Vector
from safe.common.utilities import get_utm_epsg
from safe.common.exceptions import GetDataError
from safe.common.qgis_raster_tools import (
    clip_raster, polygonize_gdal)
from safe.common.qgis_vector_tools import (
    split_by_polygon_in_out,
    extent_to_geo_array,
    reproject_vector_layer)


class CategoricalRasterBldgFunction(FunctionProvider):
    # noinspection PyUnresolvedReferences
    """Simple experimental impact function for inundation.

    :author ESSC
    :rating 3
    :param requires category=='hazard' and \
                    layertype=='raster'
    :param requires category=='exposure' and \
                    subcategory=='structure' and \
                    layertype=='vector'
        """
    def __init__(self):
        """Constructor."""
        self.extent = None

    class Metadata(ImpactFunctionMetadata):
        """Metadata for CategoricalRasterBldgFunction

        .. versionadded:: 2.2

        We only need to re-implement get_metadata(), all other behaviours
        are inherited from the abstract base class.
        """

        @staticmethod
        def get_metadata():
            """Return metadata as a dictionary.

            This is a static method. You can use it to get the metadata in
            dictionary format for an impact function.

            :returns: A dictionary representing all the metadata for the
                concrete impact function.
            :rtype: dict
            """
            dict_meta = {
                'id': 'CategoricalRasterBldgFunction',
                'name': tr('Hazard in Raster Building Function'),
                'impact': tr('Be flooded in given thresholds'),
                'author': 'Dianne Bencito',
                'date_implemented': 'N/A',
                'overview': tr('N/A'),
                'categories': {
                    'hazard': {
                        'definition': hazard_definition,
                        'subcategory': [hazard_all],
                        'units': [unit_categorised],
                        'layer_constraints': [layer_raster_numeric]
                    },
                    'exposure': {
                        'definition': exposure_definition,
                        'subcategory': exposure_structure,
                        'units': [
                            unit_building_type_type,
                            unit_building_generic],
                        'layer_constraints': [layer_vector_polygon]
                    }
                }
            }
            return dict_meta

    title = tr('Be flooded in given thresholds (bldgs)')

    parameters = OrderedDict([
        # This field of impact layer marks inundated roads by '1' value
        ('target_field', 'flooded'),
        # This field of the exposure layer contains
        # information about road types
        ('bldg_type_field', 'TYPE'),
        ('min threshold [m]', 1.0),
        ('max threshold [m]', float('inf')),
        ('postprocessors', OrderedDict([('RoadType', {'on': True})]))
    ])

    def get_function_type(self):
        """Get type of the impact function.

        :returns:   'qgis2.0'
        """
        return 'qgis2.0'

    def set_extent(self, extent):
        """Set up the extent of area of interest.

        It is mandatory to call this method before running the analysis.

        :param extent: Extents mutator [xmin, ymin, xmax, ymax].
        :type extent: list
        """
        self.extent = extent

    def run(self, layers):
        """Experimental impact function.

        :param layers: List of layers expected to contain at least:
            H: Polygon layer of inundation areas
            E: Vector layer of roads
        :type layers: list

        :returns: A new line layer with inundated roads marked.
        :type: safe_layer
        """
        target_field = self.parameters['target_field']
        bldg_type_field = self.parameters['bldg_type_field']
        threshold_min = self.parameters['min threshold [m]']
        threshold_max = self.parameters['max threshold [m]']

        if threshold_min > threshold_max:
            message = tr(
                'The minimal threshold is greater then the maximal specified '
                'threshold. Please check the values.')
            raise GetDataError(message)

        # Extract data
        H = get_hazard_layer(layers)    # Flood
        E = get_exposure_layer(layers)  # Roads

        question = get_question(
            H.get_name(), E.get_name(), self)

        H = H.get_layer()
        E = E.get_layer()

        #reproject self.extent to the hazard projection
        hazard_crs = H.crs()
        hazard_authid = hazard_crs.authid()

        if hazard_authid == 'EPSG:4326':
            viewport_extent = self.extent
        else:
            geo_crs = QgsCoordinateReferenceSystem()
            geo_crs.createFromSrid(4326)
            viewport_extent = extent_to_geo_array(
                QgsRectangle(*self.extent), geo_crs, hazard_crs)

        print viewport_extent

        #Align raster extent and viewport
        #assuming they are both in the same projection
        raster_extent = H.dataProvider().extent()
        clip_xmin = raster_extent.xMinimum()
        clip_xmax = raster_extent.xMaximum()
        clip_ymin = raster_extent.yMinimum()
        clip_ymax = raster_extent.yMaximum()
        if (viewport_extent[0] > clip_xmin):
            clip_xmin = viewport_extent[0]
        if (viewport_extent[1] > clip_ymin):
            clip_ymin = viewport_extent[1]
        if (viewport_extent[2] < clip_xmax):
            clip_xmax = viewport_extent[2]
        if (viewport_extent[3] < clip_ymax):
            clip_ymax = viewport_extent[3]

        height = ((viewport_extent[3] - viewport_extent[1]) /
                  H.rasterUnitsPerPixelY())
        height = int(height)
        width = ((viewport_extent[2] - viewport_extent[0]) /
                 H.rasterUnitsPerPixelX())
        width = int(width)

        raster_extent = H.dataProvider().extent()
        xmin = raster_extent.xMinimum()
        xmax = raster_extent.xMaximum()
        ymin = raster_extent.yMinimum()
        ymax = raster_extent.yMaximum()

        x_delta = (xmax - xmin) / H.width()
        x = xmin
        for i in range(H.width()):
            if abs(x - clip_xmin) < x_delta:
                # We have found the aligned raster boundary
                break
            x += x_delta
            _ = i

        y_delta = (ymax - ymin) / H.height()
        y = ymin
        for i in range(H.width()):
            if abs(y - clip_ymin) < y_delta:
                # We have found the aligned raster boundary
                break
            y += y_delta
        clip_extent = [x, y, x + width * x_delta, y + height * y_delta]

        # Clip and polygonize
        small_raster = clip_raster(
            H, width, height, QgsRectangle(*clip_extent))
        (flooded_polygon_inside, flooded_polygon_outside) = polygonize_gdal(
            small_raster, threshold_min, threshold_max)

        # Filter geometry and data using the extent
        extent = QgsRectangle(*self.extent)
        request = QgsFeatureRequest()
        request.setFilterRect(extent)

        if flooded_polygon_inside is None:
            message = tr(
                'There are no objects in the hazard layer with "value">%s.'
                'Please check the value or use other extent.' % (
                    threshold_min, ))
            raise GetDataError(message)

        #reproject the flood polygons to exposure projection
        exposure_crs = E.crs()
        exposure_authid = exposure_crs.authid()

        if hazard_authid != exposure_authid:
            flooded_polygon_inside = reproject_vector_layer(
                flooded_polygon_inside, E.crs())
            flooded_polygon_outside = reproject_vector_layer(
                flooded_polygon_outside, E.crs())

        # Clip exposure by the extent
        #extent_as_polygon = QgsGeometry().fromRect(extent)
        #no need to clip since It is using a bbox request
        #vector_layer = clip_by_polygon(
        #    E,
        #    extent_as_polygon
        #)
        # Find inundated roads, mark them
        polygon_layer = split_by_polygon_in_out(
            E,
            flooded_polygon_inside,
            flooded_polygon_outside,
            target_field, 1, request)

        target_field_index = polygon_layer.dataProvider().\
            fieldNameIndex(target_field)

        # Generate simple impact report
        epsg = get_utm_epsg(self.extent[0], self.extent[1])
        output_crs = QgsCoordinateReferenceSystem(epsg)
        transform = QgsCoordinateTransform(E.crs(), output_crs)
        bldg_len = flooded_len = 0  # Length of roads
        bldg_by_type = dict()      # Length of flooded roads by types

        bldg_data = polygon_layer.getFeatures()
        bldg_type_field_index = polygon_layer.fieldNameIndex(bldg_type_field)
        for bldg in bldg_data:
            attributes = bldg.attributes()
            bldg_type = attributes[bldg_type_field_index]
            if bldg_type.__class__.__name__ == 'QPyNullVariant':
                bldg_type = tr('Other')
            geom = bldg.geometry()
            geom.transform(transform)
            length = geom.length()
            bldg_len += length

            if not bldg_type in bldg_by_type:
                bldg_by_type[bldg_type] = {'flooded': 0, 'total': 0}
            bldg_by_type[bldg_type]['total'] += length

            if attributes[target_field_index] == 1:
                flooded_len += length
                bldg_by_type[bldg_type]['flooded'] += length
        table_body = [
            question,
            TableRow([
                tr('RBuilding Type'),
                tr('Flooded in the threshold (m)'),
                tr('Total (m)')],
                header=True),
            TableRow([tr('All'), int(flooded_len), int(bldg_len)])
        ]
        table_body.append(TableRow(
            tr('Breakdown by building type'), header=True))
        for t, v in bldg_by_type.iteritems():
            table_body.append(
                TableRow([t, int(v['flooded']), int(v['total'])])
            )

        impact_summary = Table(table_body).toNewlineFreeString()
        map_title = tr('Buildings inundated')

        style_classes = [dict(label=tr('Not Inundated'), value=0,
                              colour='#1EFC7C', transparency=50, size=0.5),
                         dict(label=tr('Inundated'), value=1,
                              colour='#F31A1C', transparency=50, size=0.5)]
        style_info = dict(target_field=target_field,
                          style_classes=style_classes,
                          style_type='categorizedSymbol')

        # Convert QgsVectorLayer to inasafe layer and return it
        polygon_layer = Vector(
            data=polygon_layer,
            name=tr('Flooded buildings'),
            keywords={
                'impact_summary': impact_summary,
                'map_title': map_title,
                'target_field': target_field},
            style_info=style_info)
        return polygon_layer
