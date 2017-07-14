from distutils.core import setup, Extension
from Cython.Build import cythonize


sources = [
    'fem2d/fem_solver.cpp',
    'fem2d/fem2d.pyx',
]

setup(ext_modules = cythonize(Extension(
    "fem2d", sources=sources,
    language="c++", extra_compile_args=['-std=c++11']
)))
