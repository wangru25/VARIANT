# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: wang.rui@nyu.edu
FilePath: /VARIANT/./variant/cli.py
Description: Command line interface for VARIANT tool.
'''


"""
CLI module for VARIANT package.

This module provides a wrapper around the CLI commands to ensure proper imports.
"""

# Import the main function from src.cli.commands
from src.cli.commands import main

# Re-export main for the entry point
__all__ = ['main']
