import os
import tempfile
from shapely import wkt
from shapely.geometry import Point, MultiPoint, MultiPolygon, Polygon
from ocgis.test.base import TestBase
from ocgis.util.helpers import write_geom_dict, get_iter
from ocgis.util.spatial.wrap import Wrapper
import numpy as np


class TestWrapper(TestBase):

    @property
    def axis_multipolygon(self):
        id1 = 'POLYGON((-10.787192 72.513808,10.051467 68.773536,10.051467 56.484070,-7.581244 57.552720,-7.581244 57.552720,-10.787192 72.513808))'
        id2 = 'POLYGON((7.914168 46.866228,23.409581 34.576763,13.257414 28.699192,6.311195 38.317035,6.311195 38.317035,6.311195 38.317035,7.914168 46.866228))'
        id3 = 'POLYGON((-14.527464 40.454333,-5.978271 26.027569,-20.405034 24.424596,-23.076657 34.576763,-14.527464 40.454333))'

        return MultiPolygon([wkt.loads(i) for i in [id1, id2, id3]])

    @property
    def axis_polygon(self):
        return self.axis_multipolygon[0]

    @property
    def axis_multipoint(self):
        return MultiPoint([i.centroid for i in self.axis_multipolygon])

    @property
    def point(self):
        return Point(-99.85, 41.54)

    @property
    def nebraska(self):
        wkt_str = 'POLYGON((-101.407393 40.001003,-102.051535 39.998918,-102.047545 40.342644,-102.047620 40.431077,-102.046031 40.697319,-102.046992 40.743130,-102.047739 40.998071,-102.621257 41.000214,-102.652271 40.998124,-103.382956 41.000316,-103.572316 40.999648,-104.051705 41.003211,-104.054012 41.388085,-104.055500 41.564222,-104.053615 41.698218,-104.053513 41.999815,-104.056219 42.614669,-104.056199 43.003062,-103.501464 42.998618,-103.005875 42.999354,-102.788384 42.995303,-102.086701 42.989887,-101.231737 42.986843,-100.198142 42.991095,-99.532790 42.992335,-99.253971 42.992389,-98.497651 42.991778,-98.457444 42.937160,-98.391204 42.920135,-98.310339 42.881794,-98.167826 42.839571,-98.144869 42.835794,-98.123117 42.820223,-98.121820 42.808360,-98.033140 42.769192,-97.995144 42.766812,-97.963558 42.773690,-97.929477 42.792324,-97.889941 42.831271,-97.888659 42.855807,-97.818643 42.866587,-97.797028 42.849597,-97.772186 42.846164,-97.725250 42.858008,-97.685752 42.836837,-97.634970 42.861285,-97.570654 42.847990,-97.506132 42.860136,-97.483159 42.857157,-97.457263 42.850443,-97.389306 42.867433,-97.311414 42.861771,-97.271457 42.850014,-97.243189 42.851826,-97.224443 42.841202,-97.211831 42.812573,-97.161422 42.798619,-97.130469 42.773923,-97.015139 42.759542,-96.979593 42.758313,-96.970003 42.752065,-96.977869 42.727308,-96.970773 42.721147,-96.908234 42.731699,-96.810140 42.704084,-96.810437 42.681341,-96.799344 42.670019,-96.722658 42.668592,-96.699060 42.657715,-96.694596 42.641163,-96.715273 42.621907,-96.714059 42.612302,-96.636672 42.550731,-96.629294 42.522693,-96.605467 42.507236,-96.584753 42.518287,-96.547215 42.520499,-96.494701 42.488459,-96.439394 42.489240,-96.396074 42.467401,-96.397890 42.441793,-96.417628 42.414777,-96.411761 42.380918,-96.424175 42.349279,-96.389781 42.328789,-96.368700 42.298023,-96.342881 42.282081,-96.332658 42.260307,-96.337708 42.229522,-96.363512 42.214042,-96.352165 42.168185,-96.285123 42.123452,-96.265483 42.048897,-96.238725 42.028438,-96.236093 42.001258,-96.202842 41.996615,-96.185217 41.980685,-96.147328 41.966254,-96.145870 41.924907,-96.159970 41.904151,-96.135623 41.862620,-96.076417 41.791469,-96.099321 41.752975,-96.099771 41.731563,-96.085557 41.704987,-96.122202 41.694913,-96.120264 41.684094,-96.099306 41.654680,-96.111307 41.599006,-96.080835 41.576000,-96.091936 41.563145,-96.085840 41.537522,-96.050172 41.524335,-96.004592 41.536663,-95.993965 41.528103,-95.996688 41.511517,-96.013451 41.492994,-96.006897 41.481954,-95.953185 41.472387,-95.935065 41.462381,-95.940056 41.394805,-95.942895 41.340077,-95.889107 41.301389,-95.897591 41.286863,-95.911202 41.308469,-95.930230 41.302056,-95.910981 41.225245,-95.922250 41.207854,-95.916100 41.194063,-95.859198 41.180537,-95.859801 41.166865,-95.876685 41.164202,-95.858274 41.109187,-95.878804 41.065871,-95.859539 41.035002,-95.860897 41.002650,-95.837603 40.974258,-95.836541 40.901108,-95.834396 40.870300,-95.846435 40.848332,-95.851790 40.792600,-95.876616 40.730436,-95.767999 40.643117,-95.757546 40.620904,-95.767479 40.589048,-95.763412 40.549707,-95.737036 40.532373,-95.692066 40.524129,-95.687413 40.561170,-95.675693 40.565835,-95.662944 40.558729,-95.658060 40.530332,-95.684970 40.512205,-95.695361 40.485338,-95.636817 40.396390,-95.634185 40.358800,-95.616201 40.346497,-95.617933 40.331418,-95.645553 40.322346,-95.646827 40.309109,-95.595532 40.309776,-95.547137 40.266215,-95.476822 40.226855,-95.466636 40.213255,-95.460952 40.173995,-95.422476 40.131743,-95.392813 40.115416,-95.384542 40.095362,-95.403784 40.080379,-95.413764 40.048111,-95.390532 40.043750,-95.371244 40.028751,-95.345067 40.024974,-95.308697 39.999407,-95.329701 39.992595,-95.780700 39.993489,-96.001253 39.995159,-96.240598 39.994503,-96.454038 39.994172,-96.801420 39.994476,-96.908287 39.996154,-97.361912 39.997380,-97.816589 39.999729,-97.929588 39.998452,-98.264165 39.998434,-98.504479 39.997129,-98.720632 39.998461,-99.064747 39.998338,-99.178201 39.999577,-99.627859 40.002987,-100.180910 40.000478,-100.191111 40.000585,-100.735049 39.999172,-100.754856 40.000198,-101.322148 40.001821,-101.407393 40.001003))'
        return self.get_buffered(wkt.loads(wkt_str))

    @property
    def iowa(self):
        wkt_str = 'POLYGON((-91.120132 40.705443,-91.129303 40.682189,-91.162644 40.656352,-91.215060 40.643859,-91.262211 40.639587,-91.375762 40.603480,-91.411271 40.573012,-91.413026 40.548034,-91.382255 40.528538,-91.374946 40.503697,-91.385551 40.447294,-91.372908 40.403032,-91.385909 40.392405,-91.418968 40.386919,-91.448747 40.371946,-91.477038 40.391012,-91.490314 40.390806,-91.500377 40.405160,-91.527691 40.410169,-91.529607 40.435086,-91.538846 40.441288,-91.533208 40.455441,-91.579383 40.463760,-91.586028 40.484519,-91.616860 40.504873,-91.622536 40.532903,-91.692081 40.551677,-91.689959 40.581202,-91.716976 40.593435,-91.741711 40.609784,-91.946370 40.608266,-92.193174 40.600088,-92.361513 40.599576,-92.646432 40.591462,-92.717815 40.589667,-93.100938 40.584347,-93.370271 40.580491,-93.562910 40.580813,-93.786303 40.578448,-94.018059 40.574022,-94.238392 40.570966,-94.485231 40.574205,-94.639876 40.575744,-94.920616 40.577218,-95.217428 40.581892,-95.382555 40.584334,-95.767479 40.589048,-95.757546 40.620904,-95.767999 40.643117,-95.876616 40.730436,-95.851790 40.792600,-95.846435 40.848332,-95.834396 40.870300,-95.836541 40.901108,-95.837603 40.974258,-95.860897 41.002650,-95.859539 41.035002,-95.878804 41.065871,-95.858274 41.109187,-95.876685 41.164202,-95.859801 41.166865,-95.859198 41.180537,-95.916100 41.194063,-95.922250 41.207854,-95.910981 41.225245,-95.930230 41.302056,-95.911202 41.308469,-95.897591 41.286863,-95.889107 41.301389,-95.942895 41.340077,-95.940056 41.394805,-95.935065 41.462381,-95.953185 41.472387,-96.006897 41.481954,-96.013451 41.492994,-95.996688 41.511517,-95.993965 41.528103,-96.004592 41.536663,-96.050172 41.524335,-96.085840 41.537522,-96.091936 41.563145,-96.080835 41.576000,-96.111307 41.599006,-96.099306 41.654680,-96.120264 41.684094,-96.122202 41.694913,-96.085557 41.704987,-96.099771 41.731563,-96.099321 41.752975,-96.076417 41.791469,-96.135623 41.862620,-96.159970 41.904151,-96.145870 41.924907,-96.147328 41.966254,-96.185217 41.980685,-96.202842 41.996615,-96.236093 42.001258,-96.238725 42.028438,-96.265483 42.048897,-96.285123 42.123452,-96.352165 42.168185,-96.363512 42.214042,-96.337708 42.229522,-96.332658 42.260307,-96.342881 42.282081,-96.368700 42.298023,-96.389781 42.328789,-96.424175 42.349279,-96.411761 42.380918,-96.417628 42.414777,-96.397890 42.441793,-96.396074 42.467401,-96.439394 42.489240,-96.480243 42.517130,-96.489337 42.564028,-96.500942 42.573885,-96.488498 42.580480,-96.512844 42.629755,-96.541165 42.662405,-96.563039 42.668513,-96.626540 42.708354,-96.640709 42.748603,-96.632980 42.776835,-96.600875 42.799558,-96.587645 42.835381,-96.573126 42.834347,-96.556211 42.846660,-96.537511 42.896906,-96.544263 42.913866,-96.514935 42.952382,-96.517148 42.986458,-96.499020 43.012050,-96.520010 43.051508,-96.479573 43.061884,-96.462094 43.075582,-96.460805 43.087872,-96.451505 43.126308,-96.473114 43.209082,-96.487245 43.217909,-96.558605 43.225489,-96.566991 43.239633,-96.559567 43.253263,-96.570722 43.263612,-96.579131 43.290074,-96.540563 43.307659,-96.522894 43.356966,-96.525053 43.384225,-96.557708 43.400727,-96.589113 43.435539,-96.583796 43.481920,-96.598315 43.499849,-96.460454 43.499718,-96.061039 43.498533,-95.866912 43.498944,-95.464775 43.499541,-95.396558 43.500334,-94.920464 43.499371,-94.859839 43.500030,-94.455238 43.498102,-94.246787 43.498948,-93.973950 43.500298,-93.653699 43.500762,-93.500830 43.500488,-93.054380 43.501457,-93.027211 43.501278,-92.558008 43.500259,-92.453169 43.499462,-92.077532 43.499153,-91.730366 43.499571,-91.611099 43.500626,-91.223566 43.500808,-91.235903 43.464684,-91.210916 43.424051,-91.198243 43.370513,-91.177048 43.353946,-91.078498 43.313297,-91.066428 43.280683,-91.069052 43.257898,-91.161354 43.147576,-91.168571 43.082888,-91.159752 43.081182,-91.152214 43.001316,-91.139121 42.925893,-91.093428 42.871440,-91.082030 42.783365,-91.066168 42.744913,-90.999182 42.707058,-90.919409 42.680677,-90.892545 42.678240,-90.745610 42.657001,-90.694791 42.637928,-90.664380 42.571391,-90.639219 42.555714,-90.625707 42.528562,-90.638456 42.509363,-90.651899 42.494700,-90.648473 42.475647,-90.605955 42.460564,-90.563711 42.421843,-90.491171 42.388791,-90.441725 42.360083,-90.427809 42.340645,-90.418112 42.263939,-90.407301 42.242661,-90.367858 42.210226,-90.323730 42.197337,-90.231063 42.159741,-90.191702 42.122710,-90.176214 42.120524,-90.166776 42.103767,-90.168226 42.061066,-90.150663 42.033453,-90.142796 41.983989,-90.154645 41.930802,-90.195965 41.806167,-90.255438 41.781769,-90.305016 41.756497,-90.326157 41.722768,-90.341262 41.649122,-90.339476 41.602831,-90.348494 41.586882,-90.423135 41.567305,-90.435098 41.543612,-90.455126 41.527579,-90.540975 41.526003,-90.600838 41.509618,-90.658929 41.462350,-90.708354 41.450093,-90.780042 41.449852,-90.844284 41.444652,-90.949800 41.421263,-91.000842 41.431112,-91.027637 41.423536,-91.055935 41.401407,-91.073429 41.334925,-91.102496 41.267848,-91.101672 41.231552,-91.056466 41.176290,-91.018402 41.165857,-90.990485 41.144404,-90.957930 41.104393,-90.954794 41.070397,-90.960851 40.950541,-90.983419 40.923965,-91.049353 40.879623,-91.089050 40.833767,-91.092895 40.761587,-91.120132 40.705443))'
        return self.get_buffered(wkt.loads(wkt_str))

    @property
    def multipoint(self):
        pt1 = self.point
        pt2 = Point(pt1.x+2, pt1.y-4)
        pt3 = Point(pt1.x+8, pt1.y-2.5)
        return MultiPoint([pt1, pt2, pt3])

    @property
    def multipolygon(self):
        return MultiPolygon([self.nebraska, self.iowa])

    @property
    def possible(self):
        return {
            'point': self.point,
            'polygon': self.nebraska,
            'multipoint': self.multipoint,
            'multipolygon': self.multipolygon,
            'axis_polygon': self.axis_polygon,
            'axis_multipolygon': self.axis_multipolygon,
            'axis_multipoint': self.axis_multipoint
        }

    @property
    def actual_unwrapped(self):
        return {'axis_multipoint': wkt.loads('MULTIPOINT(0.142094 64.038162, 13.393532 37.034936, 344.561927 31.202741)'),
                'axis_multipolygon': wkt.loads('MULTIPOLYGON(((349.212808 72.513808,360.000000 70.577645,360.000000 57.093250,352.418756 57.552720,349.212808 72.513808)),((0.000000 70.577645,10.051467 68.773536,10.051467 56.484070,0.000000 57.093250,0.000000 70.577645)),((7.914168 46.866228,23.409581 34.576763,13.257414 28.699192,6.311195 38.317035,7.914168 46.866228)),((345.472536 40.454333,354.021729 26.027569,339.594966 24.424596,336.923343 34.576763,345.472536 40.454333)))')}

    def get_buffered(self, geom):
        ret = geom.buffer(0)
        assert ret.is_valid
        return ret

    def write_geometries(self, geoms=None):
        geoms = geoms or self.possible
        for geom in geoms:
            _, path = tempfile.mkstemp(prefix=geom.geom_type, suffix='.shp')
            os.remove(path)
            write_geom_dict({1: geom}, path=path)

    def test_init(self):
        Wrapper()

    def test_unwrap(self):
        """Test different geometry types are appropriately unwrapped."""

        wrapper = Wrapper()
        path = tempfile.mkdtemp()
        for desc, geom in self.possible.iteritems():
            unwrapped = wrapper.unwrap(geom)
            if desc in self.actual_unwrapped:
                self.assertTrue(self.actual_unwrapped[desc].almost_equals(unwrapped, decimal=5))
            try:
                self.assertEqual(type(geom), type(unwrapped))
            except AssertionError:
                if desc == 'axis_polygon':
                    # by necessity of being split on the axis, this will come out as a multipolygon
                    self.assertIsInstance(unwrapped, MultiPolygon)
                else:
                    raise

            self.assertFalse(np.any(np.array(unwrapped) < 0.0))
            if isinstance(unwrapped, (MultiPolygon, Polygon)):
                it = get_iter(unwrapped)
                for polygon in it:
                    self.assertFalse(np.any(np.array(polygon.exterior) > 360.0))
            else:
                self.assertFalse(np.any(np.array(unwrapped) > 360.0))

    def test_wrap(self):
        """Test different geometry types are appropriately wrapped."""

        wrapper = Wrapper()
        for desc, geom in self.possible.iteritems():
            unwrapped = wrapper.unwrap(geom)
            wrapped = wrapper.wrap(unwrapped)
            try:
                self.assertTrue(geom.almost_equals(wrapped))
            except AssertionError:
                # the axis multipolygon when wrapped will have an extra polygon as the split portion on the axis will
                # be in two parts
                if desc == 'axis_multipolygon':
                    self.assertNumpyAllClose(np.array(wrapped.bounds), np.array(geom.bounds))
                    self.assertEqual(len(wrapped), 4)
                # polygon will also be split...
                elif desc == 'axis_polygon':
                    self.assertNumpyAllClose(np.array(wrapped.bounds), np.array(geom.bounds))
                    self.assertEqual(len(wrapped), 2)
                else:
                    raise
