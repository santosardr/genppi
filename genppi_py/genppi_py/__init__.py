"""
GenPPI Python Interface - A Python wrapper for GenPPI protein-protein interaction prediction tool.

This package provides an easy-to-use Python interface for GenPPI, automatically handling
executable downloads, model files, and providing a unified command-line interface across
all operating systems.
"""

__version__ = "0.1.7"
__author__ = "GenPPI Team"
__email__ = "santosardr@ufu.br"

from .genppi import run_genppi, get_executable_path, main

__all__ = ["run_genppi", "get_executable_path", "main", "__version__"]