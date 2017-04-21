from distutils.core import setup, Extension
from Cython.Build import cythonize


bacon = Extension("baconwrap",
                    sources=["snakebacon/bacon/baconwrap.pyx", "snakebacon/bacon/input.cpp", "snakebacon/bacon/Matrix.cpp", "snakebacon/bacon/ranfun.cpp", "snakebacon/bacon/vector.cpp", "snakebacon/bacon/kernel.cpp"],
                    language="c++",
                    library_dirs=["/usr/local/lib"],  # TODO(brews): Use platform agnostic #includes
                    include_dirs=["snakebacon/bacon"],  # TODO(brews): Use platform agnostic #includes
                    libraries=["gsl", "gslcblas", "m"],
                    extra_compile_args=["-xc++", "-lstdc++", "-shared-libgcc", "-O2", "-fopenmp", "-D_FILE_OFFSET_BITS=64", "-Wno-write-strings"])


setup(name = 'snakebacon',
      # ext_modules=cythonize([bacon], , gdb_debug=True)  # To debug, then do something like `cython --gdb baconwrap.pyx`
      ext_modules=cythonize([bacon]))
