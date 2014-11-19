import os
import shutil
from ocgis import OcgOperations
from ocgis.test.base import TestBase


class TestCombinatorial(TestBase):

    def iter_dataset(self):
        for as_request_dataset in [True, False]:
            for k in self.test_data.iterkeys():
                kwds = {}
                if k == 'cmip3_extraction':
                    dimension_map = {'R': 'projection', 'T': 'time', 'Y': 'latitude', 'X': 'longitude'}
                    kwds['dimension_map'] = dimension_map
                rd = self.test_data.get_rd(k, kwds=kwds)
                if as_request_dataset:
                    yield k, rd
                else:
                    yield k, rd.get()

    def test(self):
        for key, dataset in self.iter_dataset():
            print key
            ops = OcgOperations(dataset=dataset, output_format='nc', prefix='nc1')
            try:
                ret = ops.execute()
            except ValueError:
                # realization dimensions may not be written to netCDF yet
                if key == 'cmip3_extraction':
                    continue
                else:
                    raise
            else:
                try:
                    ops2 = OcgOperations(dataset={'uri': ret}, output_format='nc', prefix='nc2')
                    ret2 = ops2.execute()
                    self.assertNcEqual(ret, ret2, ignore_attributes={'global': ['history']})
                finally:
                    folder = os.path.split(ret)[0]
                    shutil.rmtree(folder)
