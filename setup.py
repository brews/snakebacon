from distutils.core import setup, Extension

from Cython.Build import cythonize

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
      ext_modules=cythonize([bacon]))
