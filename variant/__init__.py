"""
VARIANT: Viral mutAtion trackeR aImed At GeNome and proTein-level

A comprehensive Python framework for analyzing viral mutations, supporting both 
single-segment and multi-segment viruses.
"""

__version__ = "1.0.0"
__author__ = "Rui Wang"
__email__ = "rw3594@nyu.edu"

# Re-export everything from src for backward compatibility
from src.core import *
from src.utils import *

# Import CLI module
from src.cli import commands
