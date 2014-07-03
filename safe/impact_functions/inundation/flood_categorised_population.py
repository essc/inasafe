# coding=utf-8
import numpy
from safe.common.utilities import OrderedDict
from safe.defaults import get_defaults
from safe.impact_functions.core import (
    FunctionProvider,
    get_hazard_layer,
    get_exposure_layer,
    get_question,
    get_function_title,
    default_minimum_needs,
    evacuated_population_weekly_needs)
from safe.impact_functions.styles import flood_population_style as style_info
from safe.metadata import (
    hazard_flood,
    layer_raster_numeric,
    exposure_population,
    unit_people_per_pixel,
    hazard_definition,
    exposure_definition,
    unit_categorised)
from safe.storage.raster import Raster
from safe.common.utilities import (
    ugettext as tr,
    format_int,
    round_thousand,
    get_thousand_separator)
from safe.common.tables import Table, TableRow
from safe.impact_functions.impact_function_metadata import (
    ImpactFunctionMetadata)


class CategorisedFloodPopulationImpactFunction(FunctionProvider):
    """Plugin for impact of population as derived by categorised flood.

    :author ESSC
    :rating 2
    :param requires category=='hazard' and \
                    unit=='categorised' and \
                    layertype=='raster'

    :param requires category=='exposure' and \
                    subcategory=='population' and \
                    layertype=='raster'
    """

    class Metadata(ImpactFunctionMetadata):
        """Metadata for Categorised Flood Population Impact Function.

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
                'id': 'CategorisedFloodPopulationImpactFunction',
                'name': tr('Categorised Flood Population Impact Function'),
                'impact': tr('Be impacted by each category'),
                'author': 'AIFDR',
                'date_implemented': 'N/A',
                'overview': tr(
                    'To assess the impacts of categorized floods in raster '
                    'format on population raster layer.'),
                'categories': {
                    'hazard': {
                        'definition': hazard_definition,
                        'subcategory': hazard_flood,
                        'units': [unit_categorised],
                        'layer_constraints': [layer_raster_numeric]
                    },
                    'exposure': {
                        'definition': exposure_definition,
                        'subcategory': exposure_population,
                        'units': [unit_people_per_pixel],
                        'layer_constraints': [layer_raster_numeric]
                    }
                }
            }
            return dict_meta

    # Function documentation
    title = tr('Be affected by each flood category')
    synopsis = tr('To assess the impacts of categorized floods in raster '
                  'format on population raster layer.')
    actions = tr('Provide details about how many people would likely need '
                 'to be impacted for each category.')
    hazard_input = tr('A flood raster layer where each cell represents '
                      'the category of the flood. There should be 3 '
                      'categories: 1, 2, and 3.')
    exposure_input = tr('An exposure raster layer where each cell represent '
                        'population count.')
    output = tr('Map of population exposed to high category and a table with '
                'number of people in each category')
    detailed_description = tr(
        'This function will calculate how many people will be impacted '
        'per each category for all categories in the flood layer. '
        'Currently there should be 3 categories in the flood layer. After '
        'that it will show the result and the total amount of people that '
        'will be impacted for the hazard given.')
    limitation = tr('The number of categories is three.')

    # Configurable parameters
    defaults = get_defaults()
    parameters = OrderedDict([
        ('Categorical thresholds', [1, 2, 3]),
        ('postprocessors', OrderedDict([
            ('Gender', {'on': True}),
            ('Age', {
                'on': True,
                'params': OrderedDict([
                    ('youth_ratio', defaults['YOUTH_RATIO']),
                    ('adult_ratio', defaults['ADULT_RATIO']),
                    ('elderly_ratio', defaults['ELDERLY_RATIO'])])}),
            ('MinimumNeeds', {'on': True}),
        ])),
        ('minimum needs', default_minimum_needs())
    ])

    def run(self, layers):
        """Plugin for impact of population as derived by categorised flood.

        Input
          layers: List of layers expected to contain
              my_hazard: Raster layer of categorised flood
              my_exposure: Raster layer of population data

        Counts number of people exposed to each category of the flood

        Return
          Map of population exposed to high category
          Table with number of people in each category
        """

        # The 3 category
        high_t = self.parameters['Categorical thresholds'][2]
        medium_t = self.parameters['Categorical thresholds'][1]
        low_t = self.parameters['Categorical thresholds'][0]

        # Identify flood and exposure layers
        my_hazard = get_hazard_layer(layers)    # Categorised Hazard
        my_exposure = get_exposure_layer(layers)  # Population Raster

        question = get_question(my_hazard.get_name(),
                                my_exposure.get_name(),
                                self)

        # Extract data as numeric arrays
        C = my_hazard.get_data(nan=0.0)  # Category

        # Calculate impact as population exposed to each category
        P = my_exposure.get_data(nan=0.0, scaling=True)
        if high_t == 0:
            H = numpy.where(0, P, 0)
        else:
            H = numpy.where(C == high_t, P, 0)
        if medium_t == 0:
            M = numpy.where(0, P, 0)
        else:
            M = numpy.where(C == medium_t, P, 0)
        if low_t == 0:
            L = numpy.where(0, P, 0)
        else:
            L = numpy.where(C == low_t, P, 0)
        Z = numpy.where(C == 0, P, 0)
        O = numpy.where(((C > 0) & (C <= high_t)), P, 0)

        # Count totals
        total = int(numpy.sum(P))
        high = int(numpy.sum(H))
        medium = int(numpy.sum(M))
        low = int(numpy.sum(L))
        zero = int(numpy.sum(Z))
        total_impact = high + medium + low

        # Don't show digits less than a 1000
        total = round_thousand(total)
        total_impact = round_thousand(total_impact)
        high = round_thousand(high)
        medium = round_thousand(medium)
        low = round_thousand(low)
        zero = round_thousand(zero)

        # Calculate estimated minimum needs
        minimum_needs = self.parameters['minimum needs']

        tot_needs = evacuated_population_weekly_needs(total_impact,
                                                      minimum_needs)

        # Generate impact report for the pdf map
        table_body = [question,
                      TableRow([tr('Total Population affected by flooding '),
                                '%s' % format_int(total_impact)],
                               header=True),
                      TableRow([tr('Population within high flood areas '),
                                '%s' % format_int(high)],
                               header=True),
                      TableRow([tr('Population within medium flood areas '),
                                '%s' % format_int(medium)],
                               header=True),
                      TableRow([tr('Population within low flood areas'),
                                '%s' % format_int(low)],
                               header=True),
                      TableRow([tr('Population not affected by flooding'),
                                '%s' % format_int(zero)],
                               header=True),
                      TableRow(tr('Table below shows the weekly minimum '
                                  'needs for all evacuated people')),
                      TableRow([tr('Needs per week'), tr('Total')],
                               header=True),
                      [tr('Rice [kg]'),
                       format_int(tot_needs['rice'])],
                      [tr('Drinking Water [l]'),
                       format_int(tot_needs['drinking_water'])],
                      [tr('Clean Water [l]'),
                       format_int(tot_needs['water'])],
                      [tr('Family Kits'),
                       format_int(tot_needs['family_kits'])],
                      [tr('Toilets'),
                       format_int(tot_needs['toilets'])]]

        impact_table = Table(table_body).toNewlineFreeString()

        table_body.append(TableRow(tr('Action Checklist:'), header=True))
        table_body.append(TableRow(tr('How will warnings be disseminated?')))
        table_body.append(TableRow(tr('How will we reach stranded people?')))
        table_body.append(TableRow(tr('Do we have enough relief items?')))
        table_body.append(TableRow(tr('If yes, where are they located and how '
                                      'will we distribute them?')))
        table_body.append(TableRow(tr(
            'If no, where can we obtain additional relief items from and how '
            'will we transport them to here?')))

        # Extend impact report for on-screen display
        table_body.extend([TableRow(tr('Notes'), header=True),
                           tr('Map shows population density in high, medium '
                              'and low flood areas'),
                           tr('Total population: %s') % format_int(total)])
        impact_summary = Table(table_body).toNewlineFreeString()

        # Generate 8 equidistant classes across the range of flooded population
        # 8 is the number of classes in the predefined flood population style
        # as imported
        # noinspection PyTypeChecker
        classes = numpy.linspace(numpy.nanmin(O.flat[:]),
                                 numpy.nanmax(O.flat[:]), 8)

        # Modify labels in existing flood style to show quantities
        style_classes = style_info['style_classes']

        style_classes[0]['label'] = tr('Not Affected [%i people/cell]')\
                                    % classes[0]
        style_classes[1]['label'] = tr('Low [%i people/cell]') % classes[1]
        style_classes[2]['label'] = tr(' ') % classes[2]
        style_classes[3]['label'] = tr(' ') % classes[3]
        style_classes[4]['label'] = tr('Medium [%i people/cell]') % classes[4]
        style_classes[5]['label'] = tr(' ') % classes[5]
        style_classes[6]['label'] = tr(' ') % classes[6]
        style_classes[7]['label'] = tr('High [%i people/cell]') % classes[7]

        style_info['legend_title'] = tr('Population Density')

        # For printing map purpose
        map_title = tr('Population affected by flooding')
        legend_notes = tr(
            'Thousand separator is represented by %s' %
            get_thousand_separator())
        legend_units = tr('(people per cell)')
        legend_title = tr('Population density')

        # Create raster object and return
        R = Raster(O,
                   projection=my_hazard.get_projection(),
                   geotransform=my_hazard.get_geotransform(),
                   name=tr('Population might %s') % (
                       get_function_title(self).lower()),
                   keywords={'impact_summary': impact_summary,
                             'impact_table': impact_table,
                             'map_title': map_title,
                             'legend_notes': legend_notes,
                             'legend_units': legend_units,
                             'legend_title': legend_title},
                   style_info=style_info)
        return R
