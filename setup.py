from distutils.core import setup

import numpy as np
from Cython.Build import cythonize
from Cython.Distutils import build_ext
from setuptools import Extension

# use glob patterns to find and cythonize all .pyx files that have changed
extensions = [Extension("*", ["classes/*.pyx"]),
              Extension("*", ["classes/lines/*.pyx"]),
              Extension("*", ["classes/InfeasibilitySources/*.pyx"]),
              Extension("*", ["lib/*.pyx"])]

directives = {'linetrace': False, 'language_level': 3}

setup(include_dirs=[np.get_include()],
      ext_modules=cythonize(extensions, compiler_directives=directives),
      cmdclass={'build_ext': build_ext},
      script_args=['build_ext'],
      options={'build_ext': {
          'inplace': True
      }})
