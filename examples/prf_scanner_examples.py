#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-01-15 14:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-01-15 14:30:00
Email: rw3594@nyu.edu
FilePath: /VARIANT/examples/prf_scanner_examples.py
Description: Example script demonstrating PRF scanner functionality with different viruses and parameters.
'''

import subprocess
import sys
from pathlib import Path

def run_command(cmd, description):
    '''
    Run a command and print the result.
    
    Args:
        cmd: Command to run as list
        description: Description of what the command does
    '''
    print(f"\n{'='*60}")
    print(f"Example: {description}")
    print(f"{'='*60}")
    print(f"Command: {' '.join(cmd)}")
    print("-" * 60)
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("✅ SUCCESS")
        if result.stdout:
            print("Output:")
            print(result.stdout)
        if result.stderr:
            print("Warnings/Info:")
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print("❌ ERROR")
        print(f"Exit code: {e.returncode}")
        if e.stdout:
            print("Output:")
            print(e.stdout)
        if e.stderr:
            print("Error:")
            print(e.stderr)
    except FileNotFoundError:
        print("❌ ERROR: Command not found. Make sure you're in the VARIANT directory and have activated the conda environment.")

def main():
    '''
    Run PRF scanner examples for different viruses and configurations.
    '''
    print("PRF Scanner Examples for VARIANT")
    print("=" * 60)
    print("This script demonstrates various ways to use the PRF scanner")
    print("to detect programmed ribosomal frameshifting sites in viral genomes.")
    
    # Check if we're in the right directory
    if not Path("src/core/prf_scanner.py").exists():
        print("❌ ERROR: Please run this script from the VARIANT root directory")
        sys.exit(1)
    
    # Example 1: Basic SARS-CoV-2 PRF scanning
    run_command([
        "python", "src/core/prf_scanner.py",
        "--fasta", "data/SARS-CoV-2/refs/NC_045512.fasta",
        "--out", "examples/sars_cov2_basic",
        "--use-rnafold"
    ], "Basic SARS-CoV-2 PRF scanning with RNA structure prediction")
    
    # Example 2: Advanced SARS-CoV-2 PRF scanning with custom parameters
    run_command([
        "python", "src/core/prf_scanner.py",
        "--fasta", "data/SARS-CoV-2/refs/NC_045512.fasta",
        "--out", "examples/sars_cov2_advanced",
        "--spacer-min", "3",
        "--spacer-max", "12",
        "--window", "150",
        "--use-rnafold",
        "--organism", "sars_cov2"
    ], "Advanced SARS-CoV-2 PRF scanning with custom spacer and window parameters")
    
    # Example 3: HIV-1 PRF scanning
    run_command([
        "python", "src/core/prf_scanner.py",
        "--fasta", "data/HIV-1/refs/NC_001802.1.fasta",
        "--out", "examples/hiv1_prf",
        "--use-rnafold",
        "--organism", "human"
    ], "HIV-1 PRF scanning (looks for gag-pol frameshift site)")
    
    # Example 4: Chikungunya PRF scanning
    run_command([
        "python", "src/core/prf_scanner.py",
        "--fasta", "data/Chikungunya/refs/NC_004162.2.fasta",
        "--out", "examples/chikungunya_prf",
        "--use-rnafold",
        "--organism", "human"
    ], "Chikungunya virus PRF scanning")
    
    # Example 5: Zaire Ebola PRF scanning
    run_command([
        "python", "src/core/prf_scanner.py",
        "--fasta", "data/ZaireEbola/refs/NC_002549.1.fasta",
        "--out", "examples/zaire_ebola_prf",
        "--use-rnafold",
        "--organism", "human"
    ], "Zaire Ebola virus PRF scanning")
    
    # Example 6: H3N2 segment 1 PRF scanning
    run_command([
        "python", "src/core/prf_scanner.py",
        "--fasta", "data/H3N2/segment_1/refs/NC_007373.1.fasta",
        "--out", "examples/h3n2_segment1_prf",
        "--use-rnafold",
        "--organism", "human"
    ], "H3N2 influenza segment 1 PRF scanning")
    
    # Example 7: PRF scanning without RNA structure prediction (faster)
    run_command([
        "python", "src/core/prf_scanner.py",
        "--fasta", "data/SARS-CoV-2/refs/NC_045512.fasta",
        "--out", "examples/sars_cov2_no_rnafold"
    ], "SARS-CoV-2 PRF scanning without RNA structure prediction (faster)")
    
    print(f"\n{'='*60}")
    print("PRF Scanner Examples Complete!")
    print(f"{'='*60}")
    print("Generated files:")
    print("- examples/sars_cov2_basic.prf_candidates.csv")
    print("- examples/sars_cov2_basic.prf_candidates.bed")
    print("- examples/sars_cov2_advanced.prf_candidates.csv")
    print("- examples/sars_cov2_advanced.prf_candidates.bed")
    print("- examples/hiv1_prf.prf_candidates.csv")
    print("- examples/hiv1_prf.prf_candidates.bed")
    print("- examples/chikungunya_prf.prf_candidates.csv")
    print("- examples/chikungunya_prf.prf_candidates.bed")
    print("- examples/zaire_ebola_prf.prf_candidates.csv")
    print("- examples/zaire_ebola_prf.prf_candidates.bed")
    print("- examples/h3n2_segment1_prf.prf_candidates.csv")
    print("- examples/h3n2_segment1_prf.prf_candidates.bed")
    print("- examples/sars_cov2_no_rnafold.prf_candidates.csv")
    print("- examples/sars_cov2_no_rnafold.prf_candidates.bed")
    print("\nEach CSV file contains detailed PRF candidate information including:")
    print("- Slippery sequence motifs")
    print("- RNA secondary structure predictions")
    print("- tRNA interaction validation")
    print("- Frame context analysis")
    print("\nEach BED file contains genomic coordinates for visualization.")

if __name__ == "__main__":
    main()
