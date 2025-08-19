# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: wang.rui@nyu.edu
FilePath: /VARIANT/./variant/__init__.py
Description: Python module for __init__.
'''


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
