"""
CLI module for VARIANT package.

This module provides a wrapper around the CLI commands to ensure proper imports.
"""

# Import the main function from src.cli.commands
from src.cli.commands import main

# Re-export main for the entry point
__all__ = ['main']
