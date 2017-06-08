from Cython.Build import cythonize
from setuptools import setup, find_packages
from setuptools.extension import Extension

MAJOR = 0
MINOR = 0
MICRO = 1
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


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


if write_version:
    write_version_py()

bacon = Extension("snakebacon.bacon.baconwrap",
                  sources=["snakebacon/bacon/baconwrap.pyx",
                           "snakebacon/bacon/input.cpp",
                           "snakebacon/bacon/Matrix.cpp",
                           "snakebacon/bacon/ranfun.cpp",
                           "snakebacon/bacon/vector.cpp",
                           "snakebacon/bacon/kernel.cpp"],
                  language="c++",
                  libraries=["gsl", "openblas"],
                  extra_compile_args=["-xc++", "-lstdc++", "-shared-libgcc",
                                      "-O2", "-fopenmp"])

setup(name='snakebacon',
      version='0.0.1',
      description='snakebacon',
      url='https://github.com/brews/snakebacon',
      author='S. Brewster Malevich',
      author_email='malevich@email.arizona.edu',
      license='GPLv3',
      classifiers=[
          'Development Status :: 3 - Alpha',

          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering',
          'Topic :: Software Development',

          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

          'Programming Language :: Python :: 3'],
      keywords='marine radiocarbon c14',
      packages=find_packages(exclude=['contrib', 'docs', 'build']),
      install_requires=['setuptools', 'numpy', 'cython', 'pandas', 'matplotlib'],
      extras_require={'dev': ['check-manifest']},
      ext_modules=cythonize([bacon]))
