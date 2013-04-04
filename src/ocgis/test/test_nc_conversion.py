import unittest
from ocgis.test.base import TestBase
import ocgis
import netCDF4 as nc


class Test(TestBase):

    def test_projection_writing(self):
        uri = '/usr/local/climate_data/daymet/tmax.nc'
        variable = 'tmax'
        rd = ocgis.RequestDataset(uri=uri,variable=variable)
        ops = ocgis.OcgOperations(dataset=rd,snippet=True,output_format='nc')
        ret = ops.execute()
        ds = nc.Dataset(ret)
        self.assertTrue('lambert_conformal_conic' in ds.variables)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()