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
        import logbook
        log = logbook.log()

        for key, dataset in self.iter_dataset():
            print key
            ops = OcgOperations(dataset=dataset, output_format='nc', prefix='nc1')
            try:
                ret1 = ops.execute()
            except ValueError:
                # realization dimensions may not be written to netCDF yet
                if key == 'cmip3_extraction':
                    continue
                else:
                    raise
            else:
                try:
                    ops2 = OcgOperations(dataset={'uri': ret1}, output_format='nc', prefix='nc2')
                    ret2 = ops2.execute()
                    self.assertNcEqual(ret1, ret2, ignore_attributes={'global': ['history']})
                finally:
                    for path in [ret1, ret2]:
                        folder = os.path.split(path)[0]
                        shutil.rmtree(folder)
