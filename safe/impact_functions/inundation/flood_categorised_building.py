# -*- coding: utf-8 -*-
"""Categorised raster hazard impacting vector layer of buildings
"""

__license__ = "GPL"
__copyright__ = 'Copyright 2014, Australia Indonesia Facility for '
__copyright__ += 'Disaster Reduction'
from safe.impact_functions.core import (FunctionProvider,
                                        get_hazard_layer,
                                        get_exposure_layer,
                                        get_question)
from safe.common.utilities import (ugettext as tr,
                                   format_int)
from safe.metadata import (
    hazard_all,
    layer_raster_numeric,
    exposure_structure,
    unit_building_type_type,
    unit_building_generic,
    hazard_definition,
    layer_vector_polygon,
    exposure_definition,
    unit_categorised)
from safe.storage.vector import Vector
from safe.common.tables import Table, TableRow
from safe.engine.interpolation import assign_hazard_values_to_exposure_data
from safe.common.utilities import OrderedDict
from safe.impact_functions.impact_function_metadata import (
    ImpactFunctionMetadata)
from numpy import round


#FIXME: need to normalise all raster data Ole/Kristy
class CategorisedHazardBuildingImpactFunction(FunctionProvider):
    """Impact plugin for categorising hazard impact on building data

    :author AIFDR
    :rating 2
    :param requires category=='hazard' and \
                    unit=='categorised' and \
                    layertype=='raster'

    :param requires category=='exposure' and \
                    subcategory=='structure' and \
                    layertype=='vector'
    """

    class Metadata(ImpactFunctionMetadata):
        """Metadata for Categorised Hazard Population Impact Function.

        .. versionadded:: 2.1

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
                'id': 'CategorisedHazardBuildingImpactFunction',
                'name': tr('Categorised Hazard Building Impact Function'),
                'impact': tr('Be impacted'),
                'author': 'AIFDR',
                'date_implemented': 'N/A',
                'overview': tr(
                    'To assess the impacts of categorized hazards in raster '
                    'format on building vector layer.'),
                'categories': {
                    'hazard': {
                        'definition': hazard_definition,
                        'subcategory': hazard_all,
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

    # Function documentation
    title = tr('Be impacted by each category')
    synopsis = tr(
        'To assess the impacts of categorized hazard in raster '
        'format on structure/building raster layer.')
    actions = tr(
        'Provide details about how many building would likely need '
        'to be affected for each category.')
    hazard_input = tr(
        'A hazard raster layer where each cell represents '
        'the category of the hazard. There should be 3 '
        'categories: 1, 2, and 3.')
    exposure_input = tr(
        'Vector polygon layer which can be extracted from OSM '
        'where each polygon represents the footprint of a building.')
    output = tr(
        'Map of structure exposed to high category and a table with '
        'number of structure in each category')
    detailed_description = tr(
        'This function will calculate how many buildings will be affected '
        'per each category for all categories in the hazard layer. '
        'Currently there should be 3 categories in the hazard layer. After '
        'that it will show the result and the total of buildings that '
        'will be affected for the hazard given.')
    limitation = tr('The number of categories is three.')
    statistics_type = 'class_count'
    statistics_classes = [0, 1, 2, 3]
    parameters = OrderedDict([
        ('Categorical thresholds', [1, 2, 3]),
        ('postprocessors', OrderedDict([
            ('AggregationCategorical', {'on': True})]))])

    def run(self, layers):
        """Impact plugin for hazard impact.

        Counts number of building exposed to each categorised hazard zones.

        :param layers: List of layers expected to contain.
                * hazard_layer: Hazard layer of volcano
                * exposure_layer: Vector layer of structure data on
                the same grid as hazard_layer

        :returns: Map of building exposed to volcanic hazard zones.
                  Table with number of buildings affected
        :rtype: dict
        """

        # Extract data
        H = get_hazard_layer(layers)    # Value
        E = get_exposure_layer(layers)  # Building locations

        question = get_question(H.get_name(),
                                E.get_name(),
                                self)

        # Interpolate hazard level to building locations
        H = assign_hazard_values_to_exposure_data(H, E,
                                                  attribute_name='hazard_lev',
                                                  mode='constant')

        # Extract relevant numerical data
        coordinates = H.get_geometry()
        category = H.get_data()
        N = len(category)

        # The 3 category
        high_t = self.parameters['Categorical thresholds'][2]
        medium_t = self.parameters['Categorical thresholds'][1]
        low_t = self.parameters['Categorical thresholds'][0]

        # List attributes to carry forward to result layer
        #attributes = E.get_data()

        # Calculate building impact according to guidelines
        count3 = 0
        count2 = 0
        count1 = 0
        count0 = 0
        building_impact = []
        for i in range(N):
            # Get category value
            val = float(category[i]['hazard_lev'])

            ## FIXME it would be good if the affected were words not numbers
            ## FIXME need to read hazard layer and see category or keyword
            val = int(round(val))
            if val == high_t:
                affected = 3
                count3 += 1
            elif val == medium_t:
                affected = 2
                count2 += 1
            elif val == low_t:
                affected = 1
                count1 += 1
            elif val == 0:
                affected = 0
                count0 += 1

            # Collect depth and calculated damage
            result_dict = {self.target_field: affected,
                           'CATEGORY': val}

            # Record result for this feature
            building_impact.append(result_dict)

        # Create impact report
        # Generate impact summary
            table_body = [question,
                          TableRow([tr('Hazard Level'),
                                    tr('Buildings Affected')],
                          header=True),
                          TableRow([tr('Buildings in High flood area'),
                                    format_int(count3)]),
                          TableRow([tr('Buildings in Medium flood area'),
                                    format_int(count2)]),
                          TableRow([tr('Buildings in Low flood area'),
                                    format_int(count1)]),
                          TableRow([tr('Buildings in No flood area'),
                                    format_int(count0)]),
                          TableRow([tr('All Buildings'), format_int(N)])]

        table_body.append(TableRow(tr('Action Checklist:'), header=True))
        table_body.append(TableRow(
            tr('Are the critical facilities still open?')))
        table_body.append(TableRow(
            tr('Which structures have warning capacity (eg. sirens, speakers, '
               'etc.)?')))
        table_body.append(TableRow(
            tr('Which buildings will be evacuation centres?')))
        table_body.append(TableRow(
            tr('Where will we locate the operations centre?')))
        table_body.append(TableRow(
            tr('Where will we locate warehouse and/or distribution centres?')))


        table_body.append(TableRow(tr('Notes'), header=True))
        table_body.append(tr('Categorised hazard has only 3'
                             ' classes, high, medium and low.'))

        impact_summary = Table(table_body).toNewlineFreeString()
        impact_table = impact_summary
        map_title = tr('Categorised hazard impact on buildings')

        #FIXME it would be great to do categorized rather than grduated
        # Create style

        style_classes = [dict(label=tr('Not Flooded'), value=0,
                              colour='#1EFC7C', transparency=0, size=1),
                         dict(label=tr('Low'), value=1,
                              colour='#EBF442', transparency=0, size=1),
                         dict(label=tr('Medium'), value=2,
                              colour='#F4A442', transparency=0, size=1),
                         dict(label=tr('High'), value=3,
                              colour='#F31A1C', transparency=0, size=1)]
        style_info = dict(target_field=self.target_field,
                          style_classes=style_classes,
                          style_type='categorizedSymbol')

        # Create vector layer and return
        name = 'Buildings Affected by Flooding'

        V = Vector(data=building_impact,
                   projection=E.get_projection(),
                   geometry=coordinates,
                   geometry_type=E.geometry_type,
                   keywords={'impact_summary': impact_summary,
                             'impact_table': impact_table,
                             'map_title': map_title,
                             'target_field': self.target_field,
                             'statistics_type': self.statistics_type,
                             'statistics_classes': self.statistics_classes},
                   name=name,
                   style_info=style_info)
        return V
