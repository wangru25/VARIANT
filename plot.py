#!/usr/bin/env python3
"""
Unified plotting script for VARIANT package visualizations.
"""

import argparse
import sys
import os
from src.visualization.figure1_mutation_analysis import main as mutation_analysis_main
from src.visualization.figure2_row_hot_mutations import main as row_hot_main
from src.visualization.figure3_PRF import main as prf_main

def main():
    parser = argparse.ArgumentParser(
        description='VARIANT - Unified Visualization Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        '--type', 
        type=str, 
        required=False,
        choices=['mutation', 'row-hot', 'prf'],
        help='Type of visualization: mutation (combined analysis), row-hot (row/hot mutations), prf (PRF regions)'
    )
    
    parser.add_argument(
        '--virus', 
        type=str, 
        required=False,
        help='Virus name (e.g., SARS-CoV-2, HIV-1, H3N2, Chikungunya, ZaireEbola)'
    )
    
    parser.add_argument(
        '--genome-id', 
        type=str, 
        help='Specific genome ID to analyze (optional, will auto-detect if not provided)'
    )
    
    parser.add_argument(
        '--output', 
        type=str, 
        help='Output file path (optional, will auto-generate if not provided)'
    )
    
    parser.add_argument(
        '--list-viruses', 
        action='store_true',
        help='List available viruses in the data directory'
    )
    
    args = parser.parse_args()
    
    # List viruses if requested
    if args.list_viruses:
        import os
        data_dir = "data"
        if os.path.exists(data_dir):
            virus_dirs = [d for d in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, d)) and not d.startswith('.')]
            print("Available viruses:")
            for virus in virus_dirs:
                print(f"  - {virus}")
        else:
            print("Data directory not found")
        return
    
    # Check required arguments
    if not args.type:
        print("Error: --type is required (use --list-viruses to see available viruses)")
        sys.exit(1)
    
    if not args.virus:
        print("Error: --virus is required (use --list-viruses to see available viruses)")
        sys.exit(1)
    
    # Route to appropriate visualization
    if args.type == 'mutation':
        print("Generating mutation analysis visualization...")
        # Set up sys.argv for the mutation analysis module
        original_argv = sys.argv.copy()
        sys.argv = ['figure1_mutation_analysis.py']
        if args.virus:
            sys.argv.extend(['--virus', args.virus])
        if args.genome_id:
            sys.argv.extend(['--genome-id', args.genome_id])
        if args.output:
            sys.argv.extend(['--output', args.output])
        
        try:
            mutation_analysis_main()
        finally:
            # Restore original argv
            sys.argv = original_argv
    elif args.type == 'row-hot':
        print("Generating row/hot mutations visualization...")
        # Set up sys.argv for the row-hot module
        original_argv = sys.argv.copy()
        sys.argv = ['figure2_row_hot_mutations.py']
        if args.virus:
            sys.argv.extend(['--virus', args.virus])
        if args.genome_id:
            sys.argv.extend(['--genome-id', args.genome_id])
        if args.output:
            sys.argv.extend(['--output', args.output])
        
        try:
            row_hot_main()
        finally:
            # Restore original argv
            sys.argv = original_argv
    elif args.type == 'prf':
        print("Generating PRF regions visualization...")
        # Set up sys.argv for the PRF module
        original_argv = sys.argv.copy()
        sys.argv = ['figure3_PRF.py']
        if args.virus:
            sys.argv.extend(['--virus', args.virus])
        if args.genome_id:
            sys.argv.extend(['--genome-id', args.genome_id])
        if args.output:
            sys.argv.extend(['--output', args.output])
        
        try:
            prf_main()
        finally:
            # Restore original argv
            sys.argv = original_argv
    else:
        print(f"Unknown visualization type: {args.type}")
        sys.exit(1)

if __name__ == "__main__":
    main()
