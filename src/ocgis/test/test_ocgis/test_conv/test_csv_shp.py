import os
from copy import deepcopy

import fiona

from ocgis import constants
from ocgis.test.base import TestBase
from ocgis.api.operations import OcgOperations
from ocgis.util.shp_cabinet import ShpCabinetIterator


class Test(TestBase):
    def test_geometries_not_duplicated_with_equivalent_ugid(self):
        # if geometries are equivalent, they should not have duplicates in the output shapefile.
        rd = self.test_data.get_rd('cancm4_tas')
        rd2 = self.test_data.get_rd('cancm4_tasmax_2011')
        ops = OcgOperations(dataset=[rd, rd2], geom='state_boundaries', select_ugid=[16],
                            output_format=constants.OUTPUT_FORMAT_CSV_SHAPEFILE, snippet=True)
        ops.execute()

        path_shp = os.path.join(self.current_dir_output, ops.prefix, 'shp', ops.prefix + '_ugid.shp')
        with fiona.open(path_shp) as source:
            self.assertEqual(len(list(source)), 1)

    def test_geometries_different_ugid(self):
        # equivalent geometries with different ugid values should be included
        row = list(ShpCabinetIterator(key='state_boundaries', select_ugid=[16]))
        row.append(deepcopy(row[0]))
        row[1]['properties']['UGID'] = 17
        rd = self.test_data.get_rd('cancm4_tas')
        rd2 = self.test_data.get_rd('cancm4_tasmax_2011')
        ops = OcgOperations(dataset=[rd, rd2], geom=row, output_format=constants.OUTPUT_FORMAT_CSV_SHAPEFILE,
                            snippet=True)
        ops.execute()

        path_shp = os.path.join(self.current_dir_output, ops.prefix, 'shp', ops.prefix + '_ugid.shp')
        with fiona.open(path_shp) as source:
            self.assertEqual(len(list(source)), 2)
