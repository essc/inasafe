# coding=utf-8
"""InaSAFE Disaster risk tool by Australian Aid - Flood Impact on OSM
Buildings

Contact : ole.moller.nielsen@gmail.com

.. note:: This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.

"""

from safe.metadata import (
    hazard_flood,
    layer_vector_polygon,
    layer_raster_numeric,
    exposure_structure,
    unit_building_type_type,
    hazard_definition,
    exposure_definition,
    unit_building_generic,
    unit_categorised)
from safe.common.utilities import OrderedDict
from safe.impact_functions.core import (
    FunctionProvider, get_hazard_layer, get_exposure_layer, get_question)
from safe.storage.vector import Vector
from safe.common.utilities import ugettext as tr, format_int
from safe.common.tables import Table, TableRow
from safe.engine.interpolation import assign_hazard_values_to_exposure_data
from safe.impact_functions.impact_function_metadata import (
    ImpactFunctionMetadata)
import logging
from numpy import round

LOGGER = logging.getLogger('InaSAFE')


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
                        'subcategory': hazard_flood,
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

    # parameters
    parameters = OrderedDict([
        ('Categorical thresholds', [1.0, 2.0, 3.0]),
        ('postprocessors', OrderedDict([('BuildingType', {'on': True})]))
    ])

    def run(self, layers):
        """Flood impact to buildings (e.g. from Open Street Map).

         :param layers: List of layers expected to contain.
                * my_hazard: Hazard layer of flood
                * my_exposure: Vector layer of structure data on
                the same grid as my_hazard
        """
                # The 3 category
        high_t = self.parameters['Categorical thresholds'][2]
        medium_t = self.parameters['Categorical thresholds'][1]
        low_t = self.parameters['Categorical thresholds'][0]

        # Extract data
        my_hazard = get_hazard_layer(layers)  # Depth
        my_exposure = get_exposure_layer(layers)  # Building locations

        question = get_question(
            my_hazard.get_name(),
            my_exposure.get_name(),
            self)

        # Determine attribute name for hazard levels
        if my_hazard.is_raster:
            mode = 'grid'
            hazard_attribute = 'depth'
        else:
            mode = 'regions'
            hazard_attribute = None

        # Interpolate hazard level to building locations
        I = assign_hazard_values_to_exposure_data(
            my_hazard, my_exposure, attribute_name=hazard_attribute)

        # Extract relevant exposure data
        attribute_names = I.get_attribute_names()
        attributes = I.get_data()
        N = len(I)
        # Calculate building impact
        count = 0
        count1 = 0
        count2 = 0
        count3 = 0
        buildings = {}
        affected_buildings = {}
        for i in range(N):
            # Get category value
            val = float(attributes[i][hazard_attribute])

            ## FIXME it would be good if the affected were words not numbers
            ## FIXME need to read hazard layer and see category or keyword
            val = int(round(val))
            if val == high_t:
                count3 += 1
            elif val == medium_t:
                count2 += 1
            elif val == low_t:
                count1 += 1
            elif val == 0:
                count += 1

            # Count affected buildings by usage type if available
            if 'type' in attribute_names:
                usage = attributes[i]['type']
            elif 'TYPE' in attribute_names:
                usage = attributes[i]['TYPE']
            else:
                usage = None
            if 'amenity' in attribute_names and (usage is None or usage == 0):
                usage = attributes[i]['amenity']
            if 'building_t' in attribute_names and (usage is None
                                                    or usage == 0):
                usage = attributes[i]['building_t']
            if 'office' in attribute_names and (usage is None or usage == 0):
                usage = attributes[i]['office']
            if 'tourism' in attribute_names and (usage is None or usage == 0):
                usage = attributes[i]['tourism']
            if 'leisure' in attribute_names and (usage is None or usage == 0):
                usage = attributes[i]['leisure']
            if 'building' in attribute_names and (usage is None or usage == 0):
                usage = attributes[i]['building']
                if usage == 'yes':
                    usage = 'building'

            if usage is not None and usage != 0:
                key = usage
            else:
                key = 'unknown'

            if key not in buildings:
                buildings[key] = 0
                affected_buildings[key] = 0

            # Count all buildings by type
            buildings[key] += 1
            if val is True:
                # Count affected buildings by type
                affected_buildings[key] += 1

            # Add calculated impact to existing attributes
            attributes[i][self.target_field] = int(val)

        # Lump small entries and 'unknown' into 'other' category
        for usage in buildings.keys():
            val = buildings[usage]
            if val < 25 or usage == 'unknown':
                if 'other' not in buildings:
                    buildings['other'] = 0
                    affected_buildings['other'] = 0

                buildings['other'] += val
                affected_buildings['other'] += affected_buildings[usage]
                del buildings[usage]
                del affected_buildings[usage]

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
                                    format_int(count)]),
                          TableRow([tr('All Buildings'), format_int(N)])]

        school_closed = 0
        hospital_closed = 0
        # Generate break down by building usage type is available
        list_type_attribute = [
            'TYPE', 'type', 'amenity', 'building_t', 'office',
            'tourism', 'leisure', 'building']
        intersect_type = set(attribute_names) & set(list_type_attribute)
        if len(intersect_type) > 0:
            # Make list of building types
            building_list = []
            for usage in buildings:
                building_type = usage.replace('_', ' ')

                # Lookup internationalised value if available
                building_type = tr(building_type)
                building_list.append([
                    building_type.capitalize(),
                    format_int(affected_buildings[usage]),
                    format_int(buildings[usage])])
                if building_type == 'school':
                    school_closed = affected_buildings[usage]
                if building_type == 'hospital':
                    hospital_closed = affected_buildings[usage]

            # Sort alphabetically
            building_list.sort()

            table_body.append(TableRow(tr('Breakdown by building type'),
                                       header=True))
            for row in building_list:
                s = TableRow(row)
                table_body.append(s)

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

        if school_closed > 0:
            table_body.append(TableRow(
                tr('Where will the students from the %s closed schools go to '
                   'study?') % format_int(school_closed)))

        if hospital_closed > 0:
            table_body.append(TableRow(
                tr('Where will the patients from the %s closed hospitals go '
                   'for treatment and how will we transport them?') %
                format_int(hospital_closed)))

        table_body.append(TableRow(tr('Notes'), header=True))
        table_body.append(tr('Map shows buildings affected in'
                             ' low, medium and flood areas.'))

        # Result
        impact_summary = Table(table_body).toNewlineFreeString()
        impact_table = impact_summary

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

        # For printing map purpose
        map_title = tr('Buildings affected by flooding')
        legend_units = tr('(inundated or not inundated)')
        legend_title = tr('Structure inundated status')

        # Create vector layer and return
        vector_layer = Vector(
            data=attributes,
            projection=I.get_projection(),
            geometry=I.get_geometry(),
            name=tr('Estimated buildings affected'),
            keywords={
                'impact_summary': impact_summary,
                'impact_table': impact_table,
                'target_field': self.target_field,
                'map_title': map_title,
                'legend_units': legend_units,
                'legend_title': legend_title,
                'buildings_total': N,
                'buildings_affected': count1 + count2 + count3},
            style_info=style_info)
        return vector_layer
