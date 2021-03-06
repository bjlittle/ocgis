import os.path
import abc
import csv
import logging

from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.polygon import Polygon
import fiona

from ocgis.api.request.driver.vector import DriverVector
from ocgis import constants
from ocgis.interface.base.field import Field
from ocgis.conv.meta import MetaConverter
from ocgis.util.inspect import Inspect
from ocgis.util.logging_ocgis import ocgis_lh


class AbstractConverter(object):
    """
    Base converter object. Intended for subclassing.

    :param colls: A sequence of :class:`~ocgis.SpatialCollection` objects.
    :type colls: sequence of :class:`~ocgis.SpatialCollection`
    :param str outdir: Path to the output directory.
    :param str prefix: The string prepended to the output file or directory.
    :param ops: Optional operations definition. This is required for some converters.
    :type ops: :class:`~ocgis.OcgOperations`
    :param bool add_meta: If False, do not add a source and OCGIS metadata file.
    :param bool add_auxiliary_files: If False, do not create an output folder. Write only the target ouput file.
    :param bool overwrite: If True, attempt to overwrite any existing output files.
    """

    __metaclass__ = abc.ABCMeta
    _ext = None
    _add_did_file = True  # add a descriptor file for the request datasets
    _add_ugeom = False  # added user geometry in the output folder
    _add_ugeom_nest = True  # nest the user geometry in a shp folder
    _add_source_meta = True  # add a source metadata file
    _use_upper_keys = True  # if headers should be capitalized

    def __init__(self, colls, outdir=None, prefix=None, ops=None, add_meta=True, add_auxiliary_files=True,
                 overwrite=False):
        self.colls = colls
        self.ops = ops
        self.prefix = prefix
        self.outdir = outdir
        self.add_meta = add_meta
        self.add_auxiliary_files = add_auxiliary_files
        self.overwrite = overwrite
        self._log = ocgis_lh.get_logger('conv')

        if self._ext is None:
            self.path = self.outdir
        else:
            self.path = os.path.join(self.outdir, prefix + '.' + self._ext)
            if os.path.exists(self.path):
                if not self.overwrite:
                    msg = 'Output path exists "{0}" and must be removed before proceeding. Set "overwrite" argument or env.OVERWRITE to True to overwrite.'.format(
                        self.path)
                    raise IOError(msg)

        ocgis_lh('converter initialized', level=logging.DEBUG, logger=self._log)

    def get_headers(self, coll):
        """
        :type coll: :class:`ocgis.SpatialCollection`
        :returns: A list of headers from the first element return from the collection iterator.
        :rtype: [str, ...]
        """

        ret = self.get_iter_from_spatial_collection(coll)
        ret = ret.next()
        ret = ret[1].keys()
        return ret

    def get_iter_from_spatial_collection(self, coll):
        """
        :type coll: :class:`ocgis.SpatialCollection`
        :returns: A generator from the input collection.
        :rtype: generator
        """

        itr = coll.get_iter_dict(use_upper_keys=self._use_upper_keys)
        return itr

    def _build_(self, *args, **kwargs):
        raise NotImplementedError

    def _clean_outdir_(self):
        """
        Remove previous output file from outdir.
        """

    def _get_return_(self):
        return self.path

    def _write_coll_(self, f, coll):
        raise NotImplementedError

    def _finalize_(self, *args, **kwargs):
        raise NotImplementedError

    def _get_or_create_shp_folder_(self):
        path = os.path.join(self.outdir, 'shp')
        if not os.path.exists(path):
            os.mkdir(path)
        return path

    def _get_should_append_to_unique_geometry_store_(self, store, geom, ugid):
        """
        :param sequence store:
        :param :class:`shapely.Geometry` geom:
        :param int ugid:
        """

        ret = True
        test_all = []
        for row in store:
            test_geom = row['geom'].almost_equals(geom)
            test_ugid = row['ugid'] == ugid
            test_all.append(all([test_geom, test_ugid]))
        if any(test_all):
            ret = False
        return ret

    def write(self):
        ocgis_lh('starting write method', self._log, logging.DEBUG)

        unique_geometry_store = []

        # indicates if user geometries should be written to file
        write_ugeom = False

        try:
            build = True

            for coll in iter(self.colls):
                if build:

                    # write the user geometries if configured and there is one present on the incoming collection.
                    if self._add_ugeom and coll.geoms.values()[0] is not None:
                        write_ugeom = True

                    f = self._build_(coll)
                    if write_ugeom:
                        ugid_shp_name = self.prefix + '_ugid.shp'
                        ugid_csv_name = self.prefix + '_ugid.csv'

                        if self._add_ugeom_nest:
                            fiona_path = os.path.join(self._get_or_create_shp_folder_(), ugid_shp_name)
                        else:
                            fiona_path = os.path.join(self.outdir, ugid_shp_name)

                        if coll.meta is None:
                            # convert the collection properties to fiona properties
                            from fiona_ import AbstractFionaConverter

                            fiona_properties = {}
                            archetype_properties = coll.properties.values()[0]
                            for name in archetype_properties.dtype.names:
                                fiona_properties[name] = AbstractFionaConverter.get_field_type(
                                    type(archetype_properties[name][0]))

                            fiona_schema = {'geometry': 'MultiPolygon', 'properties': fiona_properties}
                            fiona_meta = {'schema': fiona_schema, 'driver': 'ESRI Shapefile'}
                        else:
                            fiona_meta = coll.meta

                        # always use the CRS from the collection. shapefile metadata will always be WGS84, but it may be
                        # overloaded in the operations.
                        fiona_meta['crs'] = coll.crs.value

                        # selection geometries will always come out as MultiPolygon regardless if they began as points.
                        # points are buffered during the subsetting process.
                        fiona_meta['schema']['geometry'] = 'MultiPolygon'

                        fiona_object = fiona.open(fiona_path, 'w', **fiona_meta)

                    build = False
                self._write_coll_(f, coll)
                if write_ugeom:
                    # write the overview geometries to disk
                    r_geom = coll.geoms.values()[0]
                    if isinstance(r_geom, Polygon):
                        r_geom = MultiPolygon([r_geom])
                    # see if this geometry is in the unique geometry store
                    should_append = self._get_should_append_to_unique_geometry_store_(
                        unique_geometry_store,
                        r_geom,
                        coll.properties.values()[0]['UGID'])
                    if should_append:
                        unique_geometry_store.append({'geom': r_geom, 'ugid': coll.properties.values()[0]['UGID']})

                        # if it is unique write the geometry to the output files
                        coll.write_ugeom(fobject=fiona_object)
        finally:

            # errors are masked if the processing failed and file objects, etc. were not properly created. if there are
            # UnboundLocalErrors pass them through to capture the error that lead to the objects not being created.

            try:
                try:
                    self._finalize_(f)
                except UnboundLocalError:
                    pass
            except Exception as e:
                # this the exception we want to log
                ocgis_lh(exc=e, logger=self._log)
            finally:
                if write_ugeom:
                    try:
                        fiona_object.close()
                    except UnboundLocalError:
                        pass

        # the metadata and dataset descriptor files may only be written if OCGIS operations are present.
        if self.ops is not None and self.add_auxiliary_files == True:
            # added OCGIS metadata output if requested.
            if self.add_meta:
                ocgis_lh('adding OCGIS metadata file', 'conv', logging.DEBUG)
                lines = MetaConverter(self.ops).write()
                out_path = os.path.join(self.outdir, self.prefix + '_' + MetaConverter._meta_filename)
                with open(out_path, 'w') as f:
                    f.write(lines)

            # add the dataset descriptor file if requested
            if self._add_did_file:
                ocgis_lh('writing dataset description (DID) file', 'conv', logging.DEBUG)
                from ocgis.conv.csv_ import OcgDialect

                headers = ['DID', 'VARIABLE', 'ALIAS', 'URI', 'STANDARD_NAME', 'UNITS', 'LONG_NAME']
                out_path = os.path.join(self.outdir, self.prefix + '_did.csv')
                with open(out_path, 'w') as f:
                    writer = csv.writer(f, dialect=OcgDialect)
                    writer.writerow(headers)
                    for rd in self.ops.dataset.itervalues():
                        try:
                            for d in rd:
                                row = [rd.did, d['variable'], d['alias'], rd.uri]
                                try:
                                    ref_variable = rd.source_metadata['variables'][d['variable']]['attrs']
                                except KeyError:
                                    if isinstance(rd.driver, DriverVector):
                                        # not be present in metadata
                                        ref_variable = {}
                                    else:
                                        raise
                                row.append(ref_variable.get('standard_name', None))
                                row.append(ref_variable.get('units', None))
                                row.append(ref_variable.get('long_name', None))
                                writer.writerow(row)
                        except NotImplementedError:
                            if isinstance(rd, Field):
                                for variable in rd.variables.itervalues():
                                    row = [rd.uid, variable.name, variable.alias, None,
                                           variable.attrs.get('standard_name'), variable.units,
                                           variable.attrs.get('long_name')]
                                    writer.writerow(row)
                            else:
                                raise

            # add source metadata if requested
            if self._add_source_meta:
                ocgis_lh('writing source metadata file', 'conv', logging.DEBUG)
                out_path = os.path.join(self.outdir, self.prefix + '_source_metadata.txt')
                to_write = []

                for rd in self.ops.dataset.iter_request_datasets():
                    try:
                        metadata = rd.source_metadata
                        ip = Inspect(meta=metadata, uri=rd.uri)
                        to_write += ip.get_report_no_variable()
                    except AttributeError:
                        if isinstance(rd.driver, DriverVector):
                            to_write += rd.driver._inspect_get_lines_()
                        else:
                            raise

                with open(out_path, 'w') as f:
                    f.writelines('\n'.join(to_write))

        # return the internal path unless overloaded by subclasses.
        ret = self._get_return_()

        return ret

    @classmethod
    def get_converter_map(cls):
        """
        :returns: A dictionary with keys corresponding to an output format's short name. Values correspond to the
         converter class.
        :rtype: dict
        """

        from ocgis.conv.fiona_ import ShpConverter, GeoJsonConverter
        from ocgis.conv.csv_ import CsvConverter, CsvShapefileConverter
        from ocgis.conv.numpy_ import NumpyConverter
        from ocgis.conv.nc import NcConverter, NcUgrid2DFlexibleMeshConverter

        mmap = {constants.OUTPUT_FORMAT_SHAPEFILE: ShpConverter,
                constants.OUTPUT_FORMAT_CSV: CsvConverter,
                constants.OUTPUT_FORMAT_CSV_SHAPEFILE: CsvShapefileConverter,
                constants.OUTPUT_FORMAT_NUMPY: NumpyConverter,
                constants.OUTPUT_FORMAT_GEOJSON: GeoJsonConverter,
                # 'shpidx':ShpIdxConverter,
                # 'keyed':KeyedConverter,
                constants.OUTPUT_FORMAT_NETCDF: NcConverter,
                constants.OUTPUT_FORMAT_METADATA: MetaConverter,
                constants.OUTPUT_FORMAT_NETCDF_UGRID_2D_FLEXIBLE_MESH: NcUgrid2DFlexibleMeshConverter}

        return mmap

    @classmethod
    def get_converter(cls, output_format):
        """
        Return the converter based on output extensions or key.

        output_format :: str

        returns

        AbstractConverter
        """

        return cls.get_converter_map()[output_format]

    @classmethod
    def validate_ops(cls, ops):
        """
        Validate an operations object.

        :param ops: The input operations object to validate.
        :type ops: :class:`ocgis.OcgOperations`
        :raises: DefinitionValidationError
        """


class AbstractTabularConverter(AbstractConverter):
    """
    .. note:: Accepts all parameters to :class:`~ocgis.conv.base.AbstractConverter`.

    :keyword bool melted: (``=False``) If ``True``, use a melted tabular output format with variable values collected in
     a single column.
    """

    def __init__(self, *args, **kwargs):
        self.melted = kwargs.pop('melted', None) or False
        super(AbstractTabularConverter, self).__init__(*args, **kwargs)

    def get_iter_from_spatial_collection(self, coll):
        """
        :type coll: :class:`ocgis.SpatialCollection`
        :returns: A generator from the input collection.
        :rtype: generator
        """

        itr = coll.get_iter_dict(use_upper_keys=self._use_upper_keys, melted=self.melted)
        return itr
