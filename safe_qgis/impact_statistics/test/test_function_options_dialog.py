"""
InaSAFE Disaster risk assessment tool developed by AusAid and World Bank
- **GUI Test Cases.**

Contact : ole.moller.nielsen@gmail.com

.. note:: This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.

"""
__author__ = 'misugijunz@gmail.com'
__date__ = '15/10/2012'
__copyright__ = ('Copyright 2012, Australia Indonesia Facility for '
                 'Disaster Reduction')
import unittest
import sys
import os
import logging

# noinspection PyPackageRequirements
from PyQt4.QtGui import QLineEdit

from safe.common.testing import get_qgis_app
# In our tests, we need to have this line below before importing any other
# safe_qgis.__init__ to load all the configurations that we make for testing
QGIS_APP, CANVAS, IFACE, PARENT = get_qgis_app()

# Add PARENT directory to path to make test aware of other modules
pardir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../../..///'))
sys.path.append(pardir)

from safe_qgis.safe_interface import get_plugins
from third_party.odict import OrderedDict

from safe_qgis.impact_statistics.function_options_dialog import (
    FunctionOptionsDialog)

LOGGER = logging.getLogger('InaSAFE')


class FunctionOptionsDialogTest(unittest.TestCase):
    """Test the InaSAFE GUI for Configurable Impact Functions"""
    def test_build_form(self):
        """Test that we can build a form by passing it a function and params.
        """
        # noinspection PyUnresolvedReferences
        # pylint: disable=W0612
        from safe.engine.impact_functions_for_testing import \
            itb_fatality_model_configurable
        # pylint: enable=W0612
        function_id = 'I T B Fatality Function Configurable'
        function_list = get_plugins(function_id)
        assert len(function_list) == 1
        assert function_list[0].keys()[0] == function_id

        dialog = FunctionOptionsDialog(None)
        parameter = {
            'thresholds': [1.0],
            'postprocessors': {
                'Gender': {'on': True},
                'Age': {
                    'on': True,
                    'params': {
                        'youth_ratio': 0.263,
                        'elderly_ratio': 0.078,
                        'adult_ratio': 0.659}}}}

        dialog.build_form(parameter)

        assert dialog.tabWidget.count() == 2

        children = dialog.tabWidget.findChildren(QLineEdit)
        assert len(children) == 4

    def test_build_form_minimum_needs(self):
        """Test that we can build a form by passing it a function and params.
        """
        function_id = 'Flood Evacuation Function Vector Hazard'
        function_list = get_plugins(function_id)
        assert len(function_list) == 1
        assert function_list[0].keys()[0] == function_id

        dialog = FunctionOptionsDialog(None)
        parameters = {
            'thresholds': [1.0],
            'postprocessors': {
                'Gender': {'on': True},
                'Age': {
                    'on': True,
                    'params': {
                        'youth_ratio': 0.263,
                        'elderly_ratio': 0.078,
                        'adult_ratio': 0.659}}}}

        dialog.build_form(parameters)

        assert dialog.tabWidget.count() == 2

        children = dialog.tabWidget.findChildren(QLineEdit)
        assert len(children) == 4

    def test_build_widget(self):
        dialog = FunctionOptionsDialog(None)
        value = dialog.build_widget(dialog.configLayout, 'foo', [2.3])
        widget = dialog.findChild(QLineEdit)

        # initial value must be same with default
        assert value() == [2.3]

        # change to 5.9
        widget.setText('5.9')
        assert value() == [5.9]

        widget.setText('5.9, 70')
        assert value() == [5.9, 70]

        widget.setText('bar')
        try:
            value()
        except ValueError:
            ## expected to raises this exception
            pass
        else:
            raise Exception("Fail: must be raise an exception")

    def test_parse_input(self):
        function_input = {
            'thresholds': lambda: [1.0],
            'postprocessors': {
                'Gender': {'on': lambda: True},
                'Age': {
                    'on': lambda: True,
                    'params': {
                        'youth_ratio': lambda: 0.263,
                        'elderly_ratio': lambda: 0.078,
                        'adult_ratio': lambda: 0.659}}}}

        dialog = FunctionOptionsDialog(None)
        result = dialog.parse_input(function_input)
        print result
        expected = OrderedDict([
            ('thresholds', [1.0]),
            ('postprocessors', OrderedDict([
                ('Gender', OrderedDict([('on', True)])),
                ('Age', OrderedDict([
                    ('on', True),
                    ('params', OrderedDict([
                        ('elderly_ratio', 0.078),
                        ('youth_ratio', 0.263),
                        ('adult_ratio', 0.659)]))]))]))])
        # noinspection PyPep8Naming
        self.maxDiff = None
        self.assertDictEqual(result, expected)

if __name__ == '__main__':
    suite = unittest.makeSuite(FunctionOptionsDialogTest, 'test')
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite)
