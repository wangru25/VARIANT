#!/bin/bash
# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-01-15 14:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-01-15 14:30:00
Email: rw3594@nyu.edu
FilePath: /VARIANT/scripts/run_prf_scan.sh
Description: Quick PRF scanning script for all supported viruses.
'''

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to run PRF scanning for a virus
run_prf_scan() {
    local virus=$1
    local fasta_path=$2
    local output_prefix=$3
    local organism=${4:-human}
    
    print_status "Scanning $virus for PRF sites..."
    
    if [ ! -f "$fasta_path" ]; then
        print_error "FASTA file not found: $fasta_path"
        return 1
    fi
    
    # Create output directory if it doesn't exist
    mkdir -p "$(dirname "$output_prefix")"
    
    # Run PRF scanner
    python src/core/prf_scanner.py \
        --fasta "$fasta_path" \
        --out "$output_prefix" \
        --use-rnafold \
        --organism "$organism"
    
    if [ $? -eq 0 ]; then
        print_success "PRF scanning completed for $virus"
        print_status "Output files:"
        print_status "  - ${output_prefix}.prf_candidates.csv"
        print_status "  - ${output_prefix}.prf_candidates.bed"
        
        # Count candidates
        if [ -f "${output_prefix}.prf_candidates.csv" ]; then
            local count=$(tail -n +2 "${output_prefix}.prf_candidates.csv" | wc -l)
            print_status "Found $count PRF candidates"
        fi
    else
        print_error "PRF scanning failed for $virus"
        return 1
    fi
}

# Main script
main() {
    print_status "VARIANT PRF Scanner - Quick Scan for All Viruses"
    print_status "================================================="
    
    # Check if we're in the right directory
    if [ ! -f "src/core/prf_scanner.py" ]; then
        print_error "Please run this script from the VARIANT root directory"
        exit 1
    fi
    
    # Check if conda environment is activated
    if [ -z "$CONDA_DEFAULT_ENV" ] || [ "$CONDA_DEFAULT_ENV" != "variant" ]; then
        print_warning "Conda environment 'variant' not activated. Attempting to activate..."
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate variant
        if [ $? -ne 0 ]; then
            print_error "Failed to activate conda environment 'variant'"
            exit 1
        fi
    fi
    
    # Create results directory
    mkdir -p result/prf_scan_results
    
    # Scan all supported viruses
    print_status "Starting PRF scanning for all supported viruses..."
    
    # SARS-CoV-2
    run_prf_scan "SARS-CoV-2" \
        "data/SARS-CoV-2/refs/NC_045512.fasta" \
        "result/prf_scan_results/sars_cov2_prf" \
        "sars_cov2"
    
    # HIV-1
    run_prf_scan "HIV-1" \
        "data/HIV-1/refs/NC_001802.1.fasta" \
        "result/prf_scan_results/hiv1_prf" \
        "human"
    
    # Chikungunya
    run_prf_scan "Chikungunya" \
        "data/Chikungunya/refs/NC_004162.2.fasta" \
        "result/prf_scan_results/chikungunya_prf" \
        "human"
    
    # Zaire Ebola
    run_prf_scan "ZaireEbola" \
        "data/ZaireEbola/refs/NC_002549.1.fasta" \
        "result/prf_scan_results/zaire_ebola_prf" \
        "human"
    
    # H3N2 segments
    for segment in {1..8}; do
        if [ -f "data/H3N2/segment_${segment}/refs/NC_00737$((segment-1)).1.fasta" ]; then
            run_prf_scan "H3N2-segment-${segment}" \
                "data/H3N2/segment_${segment}/refs/NC_00737$((segment-1)).1.fasta" \
                "result/prf_scan_results/h3n2_segment${segment}_prf" \
                "human"
        fi
    done
    
    print_success "PRF scanning completed for all viruses!"
    print_status "Results are available in: result/prf_scan_results/"
    print_status ""
    print_status "Summary of generated files:"
    ls -la result/prf_scan_results/*.csv 2>/dev/null | while read line; do
        print_status "  $line"
    done
    ls -la result/prf_scan_results/*.bed 2>/dev/null | while read line; do
        print_status "  $line"
    done
}

# Run main function
main "$@"
