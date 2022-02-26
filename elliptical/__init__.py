# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
elliptical: a package for computing the ellipses on images.
"""

# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import sys

__minimum_python_version__ = "3.6"

class UnsupportedPythonError(Exception):
    pass

if sys.version_info < tuple((int(val) for val in __minimum_python_version__.split('.'))):
    raise UnsupportedPythonError("packagename does not support Python < {}"
                                 .format(__minimum_python_version__))
