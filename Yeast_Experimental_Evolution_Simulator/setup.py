from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("Yeast_Simulator", ["Yeast_Simulator.pyx"], include_dirs =[ np. get_include ()])]
)

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("Cell_Division", ["Cell_Division.pyx"], include_dirs =[ np. get_include ()])]
)


setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("Sub_Functions", ["Sub_Functions.pyx"], include_dirs =[ np. get_include ()])]
)

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("Mutations", ["Mutations.pyx"], include_dirs =[ np. get_include ()])]
)

