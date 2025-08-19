#!/usr/bin/env python3
"""
VARIANT - Virus Mutation Analysis Tool
Author: Rui Wang
Date: 2025-08-19
Email: wang.rui@nyu.edu
Description: Main script for virus mutation parsing with support for multiple viruses
"""

import argparse
import os
import sys
import glob
from datetime import date
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Add src to path
sys.path.append(str(Path(__file__).parent / "src"))

from src.core.mutation_processor import MutationProcessor


def create_argument_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        description="VARIANT - Virus Mutation Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "--virus", type=str, default="SARS-CoV-2",
        help="Virus name (e.g., SARS-CoV-2, HIV, H3N2)"
    )
    parser.add_argument(
        "--genome-id", type=str, default="",
        help="Specific genome ID to process (optional, use --process-all to analyze all genomes)"
    )
    parser.add_argument(
        "--msa-file", type=str,
        help="MSA file name to process (uses default from config if not specified)"
    )
    parser.add_argument(
        "--process-all", action="store_true",
        help="Process all genomes in the MSA file instead of just one"
    )
    parser.add_argument(
        "--list-viruses", action="store_true",
        help="List all available viruses in configuration"
    )
    parser.add_argument(
        "--config", type=str, default="virus_config.yaml",
        help="Path to virus configuration file"
    )
    parser.add_argument(
        "--segment", type=str,
        help="Segment name for multi-segment viruses (e.g., segment_1, segment_2)"
    )
    parser.add_argument(
        "--detect-frameshifts", action="store_true",
        help="Detect potential frameshifting sites"
    )
    
    return parser


def main():
    """Main entry point."""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Initialize processor
    try:
        processor = MutationProcessor(args.virus, args.config)
    except Exception as e:
        print(f"Error initializing mutation processor: {e}")
        return 1
    
    # List available viruses if requested
    if args.list_viruses:
        viruses = processor.list_available_viruses()
        print("Available viruses:")
        for virus in viruses:
            print(f"  - {virus}")
        return 0
    
    # Process mutations
    try:
        processor.process(args)
        return 0
    except Exception as e:
        print(f"Error during processing: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
