import os

from Cython.Build import cythonize
from setuptools import setup, find_packages
from setuptools.extension import Extension

MAJOR = 0
MINOR = 0
MICRO = 4
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)
FULLVERSION = VERSION


def write_version_py(filename=None):
    cnt = """\
version = '%s'
short_version = '%s'
"""
    if not filename:
        filename = os.path.join(
            os.path.dirname(__file__), 'snakebacon', 'version.py')

    a = open(filename, 'w')
    try:
        a.write(cnt % (FULLVERSION, VERSION))
    finally:
        a.close()

write_version_py()

bacon = Extension("snakebacon.mcmcbackends.bacon.baconwrap",
                  sources=["snakebacon/mcmcbackends/bacon/baconwrap.pyx",
                           "snakebacon/mcmcbackends/bacon/input.cpp",
                           "snakebacon/mcmcbackends/bacon/Matrix.cpp",
                           "snakebacon/mcmcbackends/bacon/ranfun.cpp",
                           "snakebacon/mcmcbackends/bacon/vector.cpp",
                           "snakebacon/mcmcbackends/bacon/kernel.cpp"],
                  language="c++",
                  libraries=["gsl", "openblas"],
                  extra_compile_args=["-xc++", "-lstdc++", "-shared-libgcc", "-O2", "-fopenmp"],
                  )

setup_kwargs = dict(name='snakebacon',
                    version=FULLVERSION,
                    description='snakebacon',
                    url='https://github.com/brews/snakebacon',
                    author='S. Brewster Malevich',
                    author_email='malevich@email.arizona.edu',
                    license='GPLv3',
                    classifiers=[
                        'Development Status :: 3 - Alpha',

                        'Operating System :: POSIX',

                        'Intended Audience :: Developers',
                        'Intended Audience :: Science/Research',

                        'Topic :: Scientific/Engineering',
                        'Topic :: Software Development',

                        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

                        'Programming Language :: Python :: 3'],
                    keywords='marine radiocarbon c14',
                    install_requires=['numpy', 'Cython', 'pandas', 'matplotlib', 'scipy'],
                    packages=find_packages(),
                    package_data={'snakebacon': ['tests/*.csv', 'mcmcbackends/bacon/Curves/*.14C']},
                    )

# Deal with bad linking bug with conda in readthedocs builds.
if os.environ.get('READTHEDOCS') != 'True':
    setup_kwargs['ext_modules'] = cythonize([bacon])

setup(**setup_kwargs)
