from setuptools import setup, Command, find_packages
import sys
from subprocess import check_call
import os
import tempfile


VERSION = '1.1.0'


########################################################################################################################
# commands
########################################################################################################################


class SimpleTestCommand(Command):
    description = 'run a simple test suite to validate installation'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        cwd = os.path.split(os.path.realpath(__file__))[0]
        path = os.path.join(cwd, 'src', 'ocgis', 'test', 'test_simple', 'run_simple.py')
        cmd = ['python', path]
        check_call(cmd)


class UninstallCommand(Command):
    description = "information on how to uninstall OCGIS"
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            import ocgis

            print('To uninstall, manually remove the Python package folder located here: {0}'.format(
                os.path.split(ocgis.__file__)[0]))
        except ImportError:
            raise (ImportError("Either OpenClimateGIS is not installed or not available on the Python path."))


class InstallDependenciesUbuntu(Command):
    description = "install on Ubuntu systems"
    user_options = []

    def run(self):
        cwd = os.getcwd()
        out = 'install_out.log'
        err = 'install_err.log'
        odir = tempfile.mkdtemp()
        stdout = open(out, 'w')
        stderr = open(err, 'w')

        def call(args):
            check_call(args, stdout=stdout, stderr=stderr)

        def install_dependency(odir, url, tarball, edir, config_flags=None, custom_make=None):
            path = tempfile.mkdtemp(dir=odir)
            os.chdir(path)
            print('downloading {0}...'.format(edir))
            call(['wget', url])
            print('extracting {0}...'.format(edir))
            call(['tar', '-xzvf', tarball])
            os.chdir(edir)
            if custom_make is None:
                print('configuring {0}...'.format(edir))
                call(['./configure'] + config_flags)
                print('building {0}...'.format(edir))
                call(['make'])
                print('installing {0}...'.format(edir))
                call(['make', 'install'])
            else:
                print('installing {0}...'.format(edir))
                custom_make()

        print('installing apt packages...')
        call(['apt-get', 'update'])
        call(['apt-get', '-y', 'install', 'g++', 'libz-dev', 'curl', 'wget', 'python-dev', 'python-setuptools',
              'python-gdal'])
        print('installing shapely...')
        call(['easy_install', 'shapely'])
        call(['easy_install', 'fiona'])

        prefix = '/usr/local'

        hdf5 = 'hdf5-1.8.10-patch1'
        hdf5_tarball = '{0}.tar.gz'.format(hdf5)
        hdf5_url = 'http://www.hdfgroup.org/ftp/HDF5/current/src/{0}'.format(hdf5_tarball)
        hdf5_flags = ['--prefix={0}'.format(prefix), '--enable-shared', '--enable-hl']
        install_dependency(odir, hdf5_url, hdf5_tarball, hdf5, hdf5_flags)

        nc4 = 'netcdf-4.2.1'
        nc4_tarball = '{0}.tar.gz'.format(nc4)
        nc4_url = 'ftp://ftp.unidata.ucar.edu/pub/netcdf/{0}'.format(nc4_tarball)
        nc4_flags = ['--prefix={0}'.format(prefix), '--enable-shared', '--enable-dap', '--enable-netcdf-4']
        os.putenv('LDFLAGS', '-L{0}/lib'.format(prefix))
        os.putenv('CPPFLAGS', '-I{0}/include'.format(prefix))
        install_dependency(odir, nc4_url, nc4_tarball, nc4, nc4_flags)
        os.unsetenv('LDFLAGS')
        os.unsetenv('CPPFLAGS')

        nc4p = 'netCDF4-1.0.4'
        nc4p_tarball = '{0}.tar.gz'.format(nc4p)
        nc4p_url = 'http://netcdf4-python.googlecode.com/files/{0}'.format(nc4p_tarball)
        call(['ldconfig'])

        def nc4p_make():
            call(['python', 'setup.py', 'install'])

        install_dependency(odir, nc4p_url, nc4p_tarball, nc4p, custom_make=nc4p_make)

        stdout.close()
        stderr.close()
        # shutil.rmtree(odir)
        os.chdir(cwd)
        print('dependencies installed.')


########################################################################################################################
# check python version
########################################################################################################################


python_version = float(sys.version_info[0]) + float(sys.version_info[1]) / 10
if python_version != 2.7:
    raise (ImportError('This software requires Python version 2.7.x. You have {0}.x'.format(python_version)))


########################################################################################################################
# set up data files for installation
########################################################################################################################

shp_parts = ['state_boundaries.cfg', 'state_boundaries.dbf', 'state_boundaries.prj', 'state_boundaries.shp',
             'state_boundaries.shx']
shp_parts = ['bin/shp/state_boundaries/{0}'.format(element) for element in shp_parts]
bin_files = ['bin/test_csv_calc_conversion_two_calculations.csv']
bin_files += shp_parts
package_data = {'ocgis.test': bin_files}

########################################################################################################################
# setup command
########################################################################################################################

setup(
    name='ocgis',
    version=VERSION,
    author='NESII/CIRES/NOAA-ESRL',
    author_email='ocgis_support@list.woc.noaa.gov',
    url='http://ncpp.github.io/ocgis/install.html#installing-openclimategis',
    license='NCSA License',
    platforms=['all'],
    packages=find_packages(where='./src'),
    package_dir={'': 'src'},
    package_data=package_data,
    cmdclass={'uninstall': UninstallCommand,
              'install_dependencies_ubuntu': InstallDependenciesUbuntu,
              'test': SimpleTestCommand},
    install_requires=['numpy', 'netCDF4', 'fiona', 'shapely'],
)
