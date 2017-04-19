from distutils.core import setup, Extension
from Cython.Build import cythonize

module1 =  Extension("baconwrap",
                    sources=["snakebacon/bacon/*.pyx"],
                    language="c++",
                    library_dirs=["/usr/local/lib"],
                    include_dirs=["/usr/local/include", "snakebacon/bacon"],
                    libraries=["gsl", "gslcblas", "m"],
                    extra_compile_args=["-O2", "-fopenmp", "-D_FILE_OFFSET_BITS=64"])

setup(name = 'snakebacon',
      ext_modules = cythonize([module1]))
