from distutils.core import setup, Extension
from Cython.Build import cythonize

module1 =  Extension("bacon_wrap",
                    sources=["snakebacon/*.pyx"],#, "snakebacon/*.c"],
                    # sources=["src/snakebacon.pdx", "src/snakebacon.pxy", "src/bacon.c", "src/blt.c", "src/events.c", "src/hist2.c", "src/input.c", "src/Matrix.c", "src/vector.c", "src/kernel.c"],
                    # headers=["src/include/bacon.h", "src/include/blt.h", "src/include/input.h", "src/include/cal.h", "src/include/ranfun.h", "src/include/Matrix.h", "src/include/twalk.h"],
                    # extra_objects=["snakebacon/*.o"],
                    language="c++",
                    library_dirs=["/usr/local/lib"],
                    include_dirs=["/usr/local/include", "snakebacon/"],
                    libraries=["gsl", "gslcblas", "m"],
                    extra_compile_args=["-O2", "-fopenmp", "-D_FILE_OFFSET_BITS=64"])
                    # extra_link_args=[...],       # if needed
                    # cmdclass = {'build_ext': build_ext})

setup(name = 'snakebacon',
      ext_modules = cythonize([module1]))
