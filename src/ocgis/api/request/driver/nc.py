from copy import deepcopy
import logging
import netCDF4 as nc

import numpy as np

from ocgis import constants
from ocgis.api.request.driver.base import AbstractDriver
from ocgis.exc import ProjectionDoesNotMatch, VariableNotFoundError, DimensionNotFound, DimensionShapeError
from ocgis.interface.base.crs import CFCoordinateReferenceSystem
from ocgis.interface.base.dimension.spatial import SpatialGridDimension, SpatialDimension
from ocgis.interface.base.variable import VariableCollection, Variable
from ocgis.interface.metadata import NcMetadata
from ocgis.interface.nc.dimension import NcVectorDimension
from ocgis.interface.nc.field import NcField
from ocgis.interface.nc.temporal import NcTemporalDimension
from ocgis.util.helpers import itersubclasses
from ocgis.util.logging_ocgis import ocgis_lh


class DriverNetcdf(AbstractDriver):
    key = 'netCDF'

    def __init__(self, *args, **kwargs):
        AbstractDriver.__init__(self, *args, **kwargs)
        self._raw_metadata = None

    @property
    def raw_metadata(self):
        if self._raw_metadata is None:
            ds = self.open()
            try:
                self._raw_metadata = NcMetadata(ds)
            finally:
                self.close(ds)
        return self._raw_metadata

    def close(self, obj):
        obj.close()

    def get_crs(self):
        crs = None
        for potential in itersubclasses(CFCoordinateReferenceSystem):
            try:
                crs = potential.load_from_metadata(self.rd._variable[0], self.rd.source_metadata)
                break
            except ProjectionDoesNotMatch:
                continue
        return crs

    def get_dimensioned_variables(self):
        metadata = self.raw_metadata
        variables = metadata['variables'].keys()
        ret = []

        ## check each variable for appropriate dimensions.
        for variable in variables:
            try:
                dim_map = get_dimension_map(variable, metadata)
            except DimensionNotFound:
                ## if a dimension is not located, then it is not an appropriate variable for subsetting.
                continue
            missing_dimensions = []

            ## these dimensions are required for subsetting.
            for required_dimension in ['X', 'Y', 'T']:
                if dim_map[required_dimension] is None:
                    missing_dimensions.append(required_dimension)

            if len(missing_dimensions) > 0:
                ## if any of the required dimensions are missing, the variable is not appropriate for subsetting.
                continue
            else:
                ret.append(variable)

        return ret

    def get_source_metadata(self):
        metadata = self.raw_metadata

        try:
            var = metadata['variables'][self.rd._variable[0]]
        except KeyError:
            raise VariableNotFoundError(self.rd.uri, self.rd._variable[0])
        if self.rd.dimension_map is None:
            metadata['dim_map'] = get_dimension_map(var['name'], metadata)
        else:
            for k, v in self.rd.dimension_map.iteritems():
                try:
                    variable_name = metadata['variables'][v]['name']
                except KeyError:
                    variable_name = None
                self.rd.dimension_map[k] = {'variable': variable_name,
                                            'dimension': v,
                                            'pos': var['dimensions'].index(v)}
                metadata['dim_map'] = self.rd.dimension_map

        return metadata

    def open(self):
        try:
            ret = nc.Dataset(self.rd.uri, 'r')
        except TypeError:
            ret = nc.MFDataset(self.rd.uri)
        return ret

    def _get_vector_dimension_(self, k, v, source_metadata):
        """
        :param str k: The string name/key of the dimension to load.
        :param dict v: A series of keyword parameters to pass to the dimension class.
        :param dict source_metadata: The request dataset's metadata as returned from
         :attr:`ocgis.api.request.base.RequestDataset.source_metadata`.
        :returns: A vector dimension object linked to the source data.
        :rtype: :class:`ocgis.interface.base.dimension.base.VectorDimension`
        """

        # this is the string axis representation
        axis_value = v['axis']
        # pull the axis information out of the dimension map
        ref_axis = source_metadata['dim_map'].get(axis_value)
        # if the axis is not represented, fill it with none. this happens when a dataset does not have a vertical
        # level or projection axis for example.
        if ref_axis is None:
            fill = None
        else:
            ref_variable = source_metadata['variables'].get(ref_axis['variable'])

            # for data with a projection/realization axis there may be no associated variable.
            try:
                ref_variable['axis'] = ref_axis
            except TypeError:
                if axis_value == 'R' and ref_variable is None:
                    ref_variable = {'axis': ref_axis, 'name': ref_axis['dimension'], 'attrs': {}}

            if len(ref_variable['dimensions']) > 1:
                msg = 'Vector dimensions must be one-dimensional. "{0}" has dimensions "{1}"'.format(k, ref_variable['dimensions'])
                raise DimensionShapeError(msg)

            # extract the data length to use when creating the source index arrays.
            length = source_metadata['dimensions'][ref_axis['dimension']]['len']
            src_idx = np.arange(0, length, dtype=constants.np_int)

            # get the target data type for the dimension
            try:
                dtype = np.dtype(ref_variable['dtype'])
            # the realization dimension may not be a associated with a variable
            except KeyError:
                if k == 'realization' and ref_variable['axis']['variable'] is None:
                    dtype = None
                else:
                    raise

            # get the name of the dimension
            name = ref_variable['axis']['dimension']

            # assemble parameters for creating the dimension class then initialize the class.
            kwds = dict(name_uid=v['name_uid'], src_idx=src_idx, data=self.rd, meta=ref_variable, axis=axis_value,
                        name_value=ref_variable.get('name'), dtype=dtype, attrs=ref_variable['attrs'].copy(),
                        name=name, name_bounds=ref_variable['axis'].get('bounds'))

            # there may be additional parameters for each dimension.
            if v['adds'] is not None:
                try:
                    kwds.update(v['adds'](ref_variable['attrs']))
                # adds may not be a callable object. assume they are a dictionary.
                except TypeError:
                    kwds.update(v['adds'])

            # check for the name of the bounds dimension in the source metadata. loop through the dimension map,
            # look for a bounds variable, and choose the bounds dimension if possible
            name_bounds_suffix = self._get_name_bounds_suffix_(source_metadata)
            kwds['name_bounds_suffix'] = name_bounds_suffix

            # create instance of the dimension
            fill = v['cls'](**kwds)

        return fill

    def _get_field_(self, format_time=True, interpolate_spatial_bounds=False):
        """
        :param bool format_time:
        :param bool interpolate_spatial_bounds:
        :raises ValueError:
        """

        # reference the request dataset's source metadata
        source_metadata = self.rd.source_metadata

        def _get_temporal_adds_(ref_attrs):
            ## calendar should default to standard if it is not present and the
            ## t_calendar overload is not used.
            calendar = self.rd.t_calendar or ref_attrs.get('calendar', None) or 'standard'

            return {'units': self.rd.t_units or ref_attrs['units'], 'calendar': calendar, 'format_time': format_time}

        # this dictionary contains additional keyword arguments for the row and column dimensions.
        adds_row_col = {'interpolate_bounds': interpolate_spatial_bounds}

        # parameters for the loading loop
        to_load = {'temporal': {'cls': NcTemporalDimension, 'adds': _get_temporal_adds_, 'axis': 'T', 'name_uid': 'tid',
                                'name': 'time'},
                   'level': {'cls': NcVectorDimension, 'adds': None, 'axis': 'Z', 'name_uid': 'lid',
                             'name': 'level'},
                   'row': {'cls': NcVectorDimension, 'adds': adds_row_col, 'axis': 'Y', 'name_uid': 'yc_id',
                           'name': 'yc'},
                   'col': {'cls': NcVectorDimension, 'adds': adds_row_col, 'axis': 'X', 'name_uid': 'xc_id',
                           'name': 'xc'},
                   'realization': {'cls': NcVectorDimension, 'adds': None, 'axis': 'R', 'name_uid': 'rlz_id',
                                   'name_value': 'rlz'}}
        loaded = {}

        for k, v in to_load.iteritems():
            fill = self._get_vector_dimension_(k, v, source_metadata)
            loaded[k] = fill

        if not {'temporal', 'row', 'col'}.issubset(set([k for k, v in loaded.iteritems() if v is not None])):
            raise ValueError('Target variable must at least have temporal, row, and column dimensions.')

        grid = SpatialGridDimension(row=loaded['row'], col=loaded['col'])

        spatial = SpatialDimension(name_uid='gid', grid=grid, crs=self.rd.crs, abstraction=self.rd.s_abstraction)

        vc = VariableCollection()
        for vdict in self.rd:
            variable_meta = deepcopy(source_metadata['variables'][vdict['variable']])
            variable_units = vdict['units'] or variable_meta['attrs'].get('units')
            dtype = np.dtype(variable_meta['dtype'])
            fill_value = variable_meta['fill_value']
            variable = Variable(vdict['variable'], vdict['alias'], units=variable_units, meta=variable_meta,
                                data=self.rd, conform_units_to=vdict['conform_units_to'], dtype=dtype,
                                fill_value=fill_value, attrs=variable_meta['attrs'].copy())
            vc.add_variable(variable)

        ret = NcField(variables=vc, spatial=spatial, temporal=loaded['temporal'], level=loaded['level'],
                      realization=loaded['realization'], meta=source_metadata.copy(), uid=self.rd.did,
                      name=self.rd.name, attrs=source_metadata['dataset'].copy())

        ## apply any subset parameters after the field is loaded
        if self.rd.time_range is not None:
            ret = ret.get_between('temporal', min(self.rd.time_range), max(self.rd.time_range))
        if self.rd.time_region is not None:
            ret = ret.get_time_region(self.rd.time_region)
        if self.rd.level_range is not None:
            try:
                ret = ret.get_between('level', min(self.rd.level_range), max(self.rd.level_range))
            except AttributeError:
                ## there may be no level dimension
                if ret.level == None:
                    msg = ("A level subset was requested but the target dataset does not have a level dimension. The "
                           "dataset's alias is: {0}".format(self.rd.alias))
                    raise (ValueError(msg))
                else:
                    raise

        return ret

    @staticmethod
    def _get_name_bounds_suffix_(source_metadata):
        """
        :param dict source_metadata: Metadata dictionary as returned from :attr:`~ocgis.RequestDataset.source_metadata`.
        :returns: The name of the bounds suffix to use when creating dimensions. If no bounds are found in the source
         metadata return ``None``.
        :rtype: str or None
        """

        name_bounds_suffix = None
        for v2 in source_metadata['dim_map'].itervalues():
            # it is possible the dimension itself is none
            try:
                if v2 is not None and v2['bounds'] is not None:
                    name_bounds_suffix = source_metadata['variables'][v2['bounds']]['dimensions'][1]
                    break
            except KeyError:
                # bounds key is likely just not there
                if 'bounds' in v2:
                    raise
        return name_bounds_suffix


def get_axis(dimvar, dims, dim):
    try:
        axis = dimvar['attrs']['axis']
    except KeyError:
        ocgis_lh('Guessing dimension location with "axis" attribute missing for variable "{0}".'.format(dimvar['name']),
                 logger='nc.dataset',
                 level=logging.WARN,
                 check_duplicate=True)
        axis = guess_by_location(dims, dim)
    return axis


def get_dimension_map(variable, metadata):
    """
    :param str variable:
    :param dict metadata:
    """
    dims = metadata['variables'][variable]['dimensions']
    mp = dict.fromkeys(['T', 'Z', 'X', 'Y'])

    ## try to pull dimensions
    for dim in dims:
        dimvar = None
        try:
            dimvar = metadata['variables'][dim]
        except KeyError:
            ## search for variable with the matching dimension
            for key, value in metadata['variables'].iteritems():
                if len(value['dimensions']) == 1 and value['dimensions'][0] == dim:
                    dimvar = metadata['variables'][key]
                    break
        ## the dimension variable may not exist
        if dimvar is None:
            ocgis_lh(logger='request.nc', exc=DimensionNotFound(dim))
        axis = get_axis(dimvar, dims, dim)
        ## pull metadata information the variable and dimension names
        mp[axis] = {'variable': dimvar['name'], 'dimension': dim}
        try:
            mp[axis].update({'pos': dims.index(dimvar['name'])})
        except ValueError:
            ## variable name may differ from the dimension name
            mp[axis].update({'pos': dims.index(dim)})

    ## look for bounds variables
    # bounds_names = set(constants.name_bounds)
    for key, value in mp.iteritems():

        if value is None:
            # this occurs for such things as levels or realizations where the dimensions is not present. the value is
            # set to none and should not be processed.
            continue

        # if the dimension is found, search for the bounds by various approaches.

        # try to get the bounds attribute from the variable directly. if the attribute is not present in the metadata
        # dictionary, continue looking for other options.
        bounds_var = metadata['variables'][value['variable']]['attrs'].get('bounds')
        var = metadata['variables'][variable]

        if bounds_var is None:
            # if no attribute is found, try some other options...

            ## if no bounds variable is found for time, it may be a climatological.
            if key == 'T':
                try:
                    bounds_var = metadata['variables'][value['variable']]['attrs']['climatology']
                    ocgis_lh('Climatological bounds found for variable: {0}'.format(var['name']), logger='request.nc',
                             level=logging.INFO)
                ## climatology is not found on time axis
                except KeyError:
                    pass

        # the bounds variable was found, but the variable is not actually present in the output file
        if bounds_var not in metadata['variables']:
            msg = 'Bounds listed for variable "{0}" but the destination bounds variable "{1}" does not exist.'.\
                format(var['name'], bounds_var)
            ocgis_lh(msg, logger='nc.driver', level=logging.WARNING, check_duplicate=True)
            bounds_var = None

        try:
            assert(isinstance(bounds_var, basestring))
        except AssertionError:
            assert(bounds_var is None)

        value.update({'bounds': bounds_var})

    return mp


def guess_by_location(dims, target):
    mp = {3: {0: 'T', 1: 'Y', 2: 'X'},
          4: {0: 'T', 2: 'Y', 3: 'X', 1: 'Z'}}
    try:
        axis_map = mp[len(dims)]
    except KeyError:
        # if there an improper number of dimensions, then the variable does not have appropriate dimensions for
        # subsetting
        raise DimensionNotFound(target)
    return axis_map[dims.index(target)]
