from ocgis.api.operations import OcgOperations
from ocgis.conv.nc import NcConverter
import numpy as np
from ocgis.test.test_simple.test_simple import nc_scope
from ocgis.test.test_ocgis.test_conv.test_base import AbstractTestConverter
import ocgis
from ocgis import constants


class TestNcConverter(AbstractTestConverter):
    
    def test_fill_value_modified(self):
        ## test the fill value is appropriately copied if reset inside the field
        coll = self.get_spatial_collection()
        ref = coll[1]['tas'].variables['tas']
        ref._dtype = np.int32
        ref._value = ref.value.astype(np.int32)
        ref._fill_value = None
        ncconv = NcConverter([coll],self.current_dir_output,'ocgis_output')
        ret = ncconv.write()
        with nc_scope(ret) as ds:
            var = ds.variables['tas']
            self.assertEqual(var.dtype,np.dtype('int32'))
            self.assertEqual(var.shape,(1,1,1))
            self.assertEqual(var._FillValue,np.ma.array([],dtype=np.dtype('int32')).fill_value)
        
    def test_fill_value_copied(self):
        rd = self.test_data.get_rd('cancm4_tas')
        with nc_scope(rd.uri) as ds:
            fill_value_test = ds.variables['tas']._FillValue
        ops = ocgis.OcgOperations(dataset=rd,snippet=True,output_format='nc')
        ret = ops.execute()
        with nc_scope(ret) as ds:
            self.assertEqual(fill_value_test,ds.variables['tas']._FillValue)

    def test_get_file_format(self):
        # use a field as the input dataset
        coll = self.get_spatial_collection(field=self.get_field())
        conv = NcConverter([coll], self.current_dir_output, 'foo')
        file_format = conv._get_file_format_()
        self.assertEqual(file_format, constants.netCDF_default_data_model)

        # add operations with a field as the dataset
        ops = OcgOperations(dataset=coll[1]['foo'], output_format='nc')
        conv = NcConverter([coll], self.current_dir_output, 'foo', ops=ops)
        file_format = conv._get_file_format_()
        self.assertEqual(file_format, constants.netCDF_default_data_model)

        # add operations and use a request dataset
        coll = self.get_spatial_collection()
        rd = self.test_data.get_rd('cancm4_tas')
        ops = OcgOperations(dataset=rd, output_format='nc')
        conv = NcConverter([coll], self.current_dir_output, 'foo', ops=ops)
        file_format = conv._get_file_format_()
        with nc_scope(rd.uri) as ds:
            self.assertEqual(file_format, ds.file_format)

    def test_write_coll(self):
        # todo: test with level
        # todo: test with bounds
        # use a field as the input dataset
        coll = self.get_spatial_collection(field=self.get_field())
        conv = NcConverter([coll], self.current_dir_output, 'foo')

        with nc_scope(conv.path, 'w') as ds:
            conv._write_coll_(ds, coll)