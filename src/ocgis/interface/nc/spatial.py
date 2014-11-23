from itertools import izip
from ocgis import constants
from ocgis.util.helpers import get_formatted_slice
from ocgis.interface.base.dimension.spatial import SpatialGridDimension
from ocgis.interface.base.variable import AbstractSourcedVariable
import numpy as np


class NcSpatialGridDimension(AbstractSourcedVariable, SpatialGridDimension):

    def __init__(self, *args, **kwargs):

        AbstractSourcedVariable.__init__(self, kwargs.pop('data', None), kwargs.pop('src_idx', None))
        SpatialGridDimension.__init__(self, *args, **kwargs)

    def __getitem__(self, item):
        ret = SpatialGridDimension.__getitem__(self, item)
        if ret._src_idx is not None:
            slice_row, slice_col = get_formatted_slice(item, 2)
            src_idx = {}
            for key, slc in izip(['row', 'col'], [slice_row, slice_col]):
                src_idx[key] = np.atleast_1d(ret._src_idx[key][slc])
            ret._src_idx = src_idx
        return ret

    @staticmethod
    def _format_src_idx_(value):
        if value is not None:
            assert isinstance(value, dict)
            assert value['row'] is not None
            assert value['col'] is not None
        return value

    def _get_uid_(self):
        if self._src_idx is not None:
            shp = (self._src_idx['row'].shape[0], self._src_idx['col'].shape[0])
        else:
            shp = None
        return SpatialGridDimension._get_uid_(self, shp=shp)

    def _set_value_from_source_(self):
        try:
            self.value = SpatialGridDimension._get_value_(self)
        except AttributeError:
            if self.row is None or self.col is None:
                ds = self._data.driver.open()
                try:
                    slices = {k: get_formatted_slice(self._src_idx[k], 1) for k in self._src_idx.keys()}
                    slice_row = slices['row']
                    slice_col = slices['col']
                    value_row = np.atleast_2d(ds.variables[self.name_row][slice_row, slice_col])
                    value_col = np.atleast_2d(ds.variables[self.name_col][slice_row, slice_col])
                    fill = np.zeros([2]+list(value_row.shape), dtype=value_row.dtype)
                    try:
                        fill_value = value_row.fill_value
                    except AttributeError:
                        fill_value = None
                    fill = np.ma.array(fill, fill_value=fill_value, mask=False)
                    fill[0, :, :] = value_row
                    fill[1, :, :] = value_col
                    self.value = fill
                finally:
                    self._data.driver.close(ds)
            else:
                raise

    def _validate_(self):
        try:
            SpatialGridDimension._validate_(self)
        except ValueError:
            if self._data is None:
                msg = 'With no value representations (i.e. row, column, value), a data source is required.'
                raise ValueError(msg)