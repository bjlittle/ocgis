import csv
import numpy as np

from ocgis import constants
from ocgis.test.base import attr
from ocgis.calc.library.index.duration import Duration, FrequencyDuration
from ocgis.exc import DefinitionValidationError
from ocgis.api.operations import OcgOperations
from ocgis.test.test_ocgis.test_calc.test_calc_general import AbstractCalcBase


class TestDuration(AbstractCalcBase):
    def test_duration(self):
        duration = Duration()

        # # three consecutive days over 3
        values = np.array([1, 2, 3, 3, 3, 1, 1], dtype=float)
        values = self.get_reshaped(values)
        ret = duration.calculate(values, 2, operation='gt', summary='max')
        self.assertEqual(3.0, ret.flatten()[0])

        # # no duration over the threshold
        values = np.array([1, 2, 1, 2, 1, 2, 1], dtype=float)
        values = self.get_reshaped(values)
        ret = duration.calculate(values, 2, operation='gt', summary='max')
        self.assertEqual(0., ret.flatten()[0])

        ## no duration over the threshold
        values = np.array([1, 2, 1, 2, 1, 2, 1], dtype=float)
        values = self.get_reshaped(values)
        ret = duration.calculate(values, 2, operation='gte', summary='max')
        self.assertEqual(1., ret.flatten()[0])

        ## average duration
        values = np.array([1, 5, 5, 2, 5, 5, 5], dtype=float)
        values = self.get_reshaped(values)
        ret = duration.calculate(values, 4, operation='gte', summary='mean')
        self.assertEqual(2.5, ret.flatten()[0])

        ## add some masked values
        values = np.array([1, 5, 5, 2, 5, 5, 5], dtype=float)
        mask = [0, 0, 0, 0, 0, 1, 0]
        values = np.ma.array(values, mask=mask)
        values = self.get_reshaped(values)
        ret = duration.calculate(values, 4, operation='gte', summary='max')
        self.assertEqual(2., ret.flatten()[0])

        ## test with an actual matrix
        values = np.array([1, 5, 5, 2, 5, 5, 5, 4, 4, 0, 2, 4, 4, 4, 3, 3, 5, 5, 6, 9], dtype=float)
        values = values.reshape(5, 2, 2)
        values = np.ma.array(values, mask=False)
        ret = duration.calculate(values, 4, operation='gte', summary='mean')
        self.assertNumpyAll(np.ma.array([4., 2., 1.5, 1.5], dtype=ret.dtype), ret.flatten())

    def test_standard_operations(self):
        ret = self.run_standard_operations(
            [{'func': 'duration', 'name': 'max_duration',
              'kwds': {'operation': 'gt', 'threshold': 2, 'summary': 'max'}}],
            capture=True)
        for cap in ret:
            reraise = True
            if isinstance(cap['exception'], DefinitionValidationError):
                if cap['parms']['calc_grouping'] in [['month'], 'all']:
                    reraise = False
            if reraise:
                raise (cap['exception'])


class TestFrequencyDuration(AbstractCalcBase):
    def test_constructor(self):
        FrequencyDuration()

    def test_calculate(self):
        fduration = FrequencyDuration()

        values = np.array([1, 2, 3, 3, 3, 1, 1, 3, 3, 3, 4, 4, 1, 4, 4, 1, 10, 10], dtype=float)
        values = self.get_reshaped(values)
        ret = fduration.calculate(values, threshold=2, operation='gt')
        self.assertEqual(ret.flatten()[0].dtype.names, ('duration', 'count'))
        self.assertNumpyAll(np.ma.array([2, 3, 5], dtype=np.int32), ret.flatten()[0]['duration'])
        self.assertNumpyAll(np.ma.array([2, 1, 1], dtype=np.int32), ret.flatten()[0]['count'])

        calc = [{'func': 'freq_duration', 'name': 'freq_duration', 'kwds': {'operation': 'gt', 'threshold': 280}}]
        ret = self.run_standard_operations(calc, capture=True, output_format=None)
        for dct in ret:
            if isinstance(dct['exception'], NotImplementedError) and dct['parms']['aggregate']:
                pass
            elif isinstance(dct['exception'], DefinitionValidationError):
                if dct['parms']['output_format'] == 'nc' or dct['parms']['calc_grouping'] == ['month']:
                    pass
            else:
                raise (dct['exception'])

    @attr('slow')
    def test_real_data_multiple_datasets(self):
        kwds = {'time_region': {'year': [1991], 'month': [7]}}
        rd_tasmax = self.test_data.get_rd('maurer_2010_concatenated_tasmax', kwds=kwds)
        rd_tasmin = self.test_data.get_rd('maurer_2010_concatenated_tasmin', kwds=kwds)

        ops = OcgOperations(dataset=[rd_tasmax, rd_tasmin],
                            output_format=constants.OUTPUT_FORMAT_CSV_SHAPEFILE,
                            calc=[{'name': 'Frequency Duration', 'func': 'freq_duration',
                                   'kwds': {'threshold': 25.0, 'operation': 'gte'}}],
                            calc_grouping=['month', 'year'],
                            geom='us_counties', select_ugid=[2778], aggregate=True,
                            calc_raw=False, spatial_operation='clip',
                            headers=['did', 'ugid', 'gid', 'year', 'month', 'day', 'variable', 'calc_key', 'value'], )
        ret = ops.execute()

        with open(ret, 'r') as f:
            reader = csv.DictReader(f)
            variables = [row['VARIABLE'] for row in reader]
        self.assertEqual(set(variables), set(['tasmax', 'tasmin']))

    def test_real_data(self):
        """Test calculations on real data."""

        rd = self.test_data.get_rd('maurer_2010_concatenated_tasmax', kwds={'time_region': {'year': [1991],
                                                                                            'month': [7]}})
        for output_format in [constants.OUTPUT_FORMAT_NUMPY, constants.OUTPUT_FORMAT_CSV_SHAPEFILE,
                              constants.OUTPUT_FORMAT_SHAPEFILE, constants.OUTPUT_FORMAT_CSV]:
            ops = OcgOperations(dataset=rd,
                                output_format=output_format, prefix=output_format,
                                calc=[{'name': 'Frequency Duration',
                                       'func': 'freq_duration',
                                       'kwds': {'threshold': 15.0, 'operation': 'gte'}}],
                                calc_grouping=['month', 'year'],
                                geom='us_counties', select_ugid=[2778], aggregate=True,
                                calc_raw=False, spatial_operation='clip',
                                headers=['did', 'ugid', 'gid', 'year', 'month', 'day', 'variable', 'calc_key',
                                         'value'],
                                melted=True)
            ret = ops.execute()

            if output_format == 'numpy':
                ref = ret[2778]['tasmax'].variables['Frequency Duration'].value
                self.assertEqual(ref.compressed()[0].shape, (2,))

            if output_format == constants.OUTPUT_FORMAT_CSV_SHAPEFILE:
                real = [{'COUNT': '1', 'UGID': '2778', 'DID': '1', 'CALC_KEY': 'freq_duration', 'MONTH': '7',
                         'DURATION': '7', 'GID': '2778', 'YEAR': '1991', 'VARIABLE': 'tasmax', 'DAY': '16'},
                        {'COUNT': '1', 'UGID': '2778', 'DID': '1', 'CALC_KEY': 'freq_duration', 'MONTH': '7',
                         'DURATION': '23', 'GID': '2778', 'YEAR': '1991', 'VARIABLE': 'tasmax', 'DAY': '16'}]
                with open(ret, 'r') as f:
                    reader = csv.DictReader(f)
                    rows = list(reader)
                for row, real_row in zip(rows, real):
                    self.assertDictEqual(row, real_row)
