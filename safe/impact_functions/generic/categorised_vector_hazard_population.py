import numpy
from third_party.odict import OrderedDict
from safe.defaults import get_defaults
from safe.impact_functions.core import (FunctionProvider,
                                        get_hazard_layer,
                                        get_exposure_layer,
                                        get_question,
                                        get_function_title)
from safe.impact_functions.styles import flood_population_style as style_info
from safe.storage.vector import Vector
from safe.storage.vector import convert_polygons_to_centroids
from safe.storage.utilities import calculate_polygon_centroid
from safe.common.polygon import is_inside_polygon
from safe.common.utilities import (ugettext as tr,
                                   format_int,
                                   round_thousand)
from safe.common.tables import Table, TableRow
from safe.engine.interpolation import assign_hazard_values_to_exposure_data


class CategorisedVectorHazardPopulationImpactFunction(FunctionProvider):
    """Plugin for impact of population as derived by categorised hazard

    :author ESSC
    :rating 2
    :param requires category=='hazard' and \
                    unit=='normalised' and \
                    layertype=='vector'

    :param requires category=='exposure' and \
                    subcategory=='population' and \
                    layertype=='vector'
    """
    # Function documentation
    title = tr('Be vectoring impacted')
    synopsis = tr('To assess the impacts of categorized hazards in vector '
                  'format on population vector layer.')
    actions = tr('Provide details about how many people would likely need '
                 'to be impacted for each category.')
    hazard_input = tr('A hazard vector layer where each cell represents '
                      'the category of the hazard. There should be 3 '
                      'categories: 1, 2, and 3.')
    exposure_input = tr('An exposure vector layer where each cell represent '
                        'population count.')
    output = tr('Map of population exposed to high category and a table with '
                'number of people in each category')
    detailed_description = \
        tr('This function will calculate how many people will be impacted '
           'per each category for all categories in the hazard layer. '
           'Currently there should be 3 categories in the hazard layer. After '
           'that it will show the result and the total amount of people that '
           'will be impacted for the hazard given.')
    limitation = tr('The number of categories is three.')
    target_field = 'Impact'

    # Configurable parameters
    defaults = get_defaults()
    parameters = OrderedDict([
        ('population field', 'pop'),
        ('hazard field', 'haz_level'),
        ('impact field', 'haz_level'),
        ('impact population count field', 'pop_impact'),
        ('categories', [0.33, 0.66, 1]),
        ('postprocessors', OrderedDict([
            ('Gender', {'on': True}),
            ('Age', {
                'on': True,
                'params': OrderedDict([
                    ('youth_ratio', defaults['YOUTH_RATIO']),
                    ('adult_ratio', defaults['ADULT_RATIO']),
                    ('elder_ratio', defaults['ELDER_RATIO'])])})]))])


    def run(self, layers):
        """Plugin for impact of population as derived by categorised hazard

        Input
          layers: List of layers expected to contain
              my_hazard: Vector layer of categorised hazard
              my_exposure: Vector layer of population data

        Counts number of people exposed to each category of the hazard

        Return
          Map of population exposed to high category
          Table with number of people in each category
        """

        # The 3 category
        high_t = 1
        medium_t = 0.66
        low_t = 0.34

        # Identify hazard and exposure layers
        my_hazard = get_hazard_layer(layers)    # Categorised Hazard
        my_exposure = get_exposure_layer(layers)  # Population Vector
        my_exposure_keywords = my_exposure.get_keywords()
        question = get_question(my_hazard.get_name(),
                                my_exposure.get_name(),
                                self)
        my_exposure = self.add_density(my_exposure)
        my_impact = self.deintersect_exposure(my_exposure, my_hazard)
        my_impact = self.assign_hazard_level(my_impact, my_hazard)
        print my_impact.get_data()
        my_impact = self.multiply_density_by_area(my_impact)
        my_impact_keywords = {'impact_summary': 'impact summary',
                     'impact_table': 'impact table',
                     'map_title': 'map title',
                     'target_field': self.parameters['impact field'],
                     'statistics_type': 'class_count',
                     'statistics_classes': self.parameters['categories']}
        my_impact_keywords.update(my_exposure_keywords)
        my_impact.keywords = my_impact_keywords
        return my_impact


    def add_density(self, exposure_layer):
        population_field = self.parameters['population field']

        # Get population data from layer
        if population_field in exposure_layer.get_attribute_names():
            D = []
            for att in exposure_layer.get_data():
                population = att[population_field]
                # FIXME (DB) area needs to be derived.
                # See hub.qgis.org/issues/9060
                density = population / float(att['area'])
                att['density'] = density
                D.append(att)
            exposure_layer.data = D

        else:
            raise RuntimeError(tr('No population field found'))
        return exposure_layer


    def deintersect_exposure(self, exposure_layer, hazard_layer):
        # FIXME DB: Need to use the _prepare_polygon layer
        impact_layer = exposure_layer
        return impact_layer


    def assign_hazard_level(self, impact_layer, hazard_layer):
        impact_centroids_geom = convert_polygons_to_centroids(
            impact_layer).get_geometry()
        impact_field = self.parameters['hazard field']
        impact_attr = impact_layer.get_data()
        hazard_field = self.parameters['hazard field']
        hazard_attr = hazard_layer.get_data()
        hazard_geom = hazard_layer.get_geometry()

        for hazard_index, hazard_poly in enumerate(hazard_geom):
            hazard_level = hazard_attr[hazard_index][hazard_field]
            for impact_index, impact_centroid in enumerate(
                    impact_centroids_geom):
                if is_inside_polygon(impact_centroid, hazard_poly):
                    try:
                        impact_attr[impact_index][impact_field]
                        raise RuntimeError(
                            tr('%s field already defined in impact layer') %
                            impact_field)
                    except KeyError:
                        impact_attr[impact_index][hazard_field] = hazard_level

        impact_layer.data = impact_attr
        return impact_layer

    def multiply_density_by_area(self, impact_layer):
        impact_data = impact_layer.get_data()
        impact_count_field = self.parameters['impact population count field']
        for index, geom in enumerate(impact_layer.get_geometry()):
            impact_attr = impact_data[index]
            impact_attr[impact_count_field] = (impact_attr['density'] *
                                               impact_attr['area'])

        return impact_layer
