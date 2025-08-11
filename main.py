"""
Author: Rui Wang
Date: 2023-10-11 10:24:49
LastModifiedBy: Rui Wang
LastEditTime: 2025-06-30 14:39:25
Email: rw3594@nyu.edu
FilePath: /7_MutParser/main.py
Description: Main script for virus mutation parsing with support for multiple viruses
"""

import argparse
import os
import sys
from datetime import date
from pathlib import Path
from typing import Dict, List, Optional

import yaml

# Add src to path
sys.path.append(str(Path(__file__).parent / "src"))

from src.core.frameshift_detector import FrameshiftDetector
from src.core.genome_processor import GenomeSNPProcessor
from src.core.mutation_detector import (
    GeneMutationDetector,
    MultipleSequenceAlignment,
)
from src.core.protein_analyzer import (
    AlignmentProcessor,
    AminoAcidMutationProcessor,
    ORFProcessor,
    ProMutationDetector,
    Proteome,
)
from src.core.reference_genome import ReferenceGenome
from src.core.snp_processor import SNPProcessor
from src.core.virus_processor import VirusMutationProcessor
from src.utils.mutation_utils import (
    classify_mutation_type,
    convert_key_to_int,
    detect_hot_mutation,
    format_amino_acid_change_for_csv,
    sort_dict_by_consecutive_keys,
    split_multi_protein_mutations,
)

# Constants
DATE_FORMAT = "%Y%m%d"

# Global Variables
# today = date.today().strftime("%Y%m%d")
today = 20250807
def main():
    """
    Main function to process virus mutations.
    Supports both single-segment and multi-segment viruses.
    """
    parser = argparse.ArgumentParser(
        description="Process virus mutations with virus-specific organization"
    )
    parser.add_argument(
        "--virus",
        type=str,
        default="SARS-CoV-2",
        help="Virus name (e.g., SARS-CoV-2, HIV, Pox, H3N2)",
    )
    parser.add_argument(
        "--genome-id",
        type=str,
        default="",
        help="Specific genome ID to process (optional, use --process-all to analyze all genomes)",
    )
    parser.add_argument(
        "--msa-file",
        type=str,
        help="MSA file name to process (uses default from config if not specified)",
    )
    parser.add_argument(
        "--process-all",
        action="store_true",
        help="Process all genomes in the MSA file instead of just one",
    )
    parser.add_argument(
        "--list-viruses",
        action="store_true",
        help="List all available viruses in configuration",
    )
    parser.add_argument(
        "--config",
        type=str,
        default="virus_config.yaml",
        help="Path to virus configuration file",
    )
    parser.add_argument(
        "--segment",
        type=str,
        help="Segment name for multi-segment viruses (e.g., segment_1, segment_2)",
    )
    parser.add_argument(
        "--detect-frameshifts",
        action="store_true",
        help="Detect potential frameshifting sites",
    )



    args = parser.parse_args()

    # Initialize virus processor
    virus_processor = VirusMutationProcessor(args.virus, args.config)

    # List available viruses if requested
    if args.list_viruses:
        available_viruses = virus_processor.list_available_viruses()
        if available_viruses:
            print("Available viruses:")
            for virus in available_viruses:
                print(f"  - {virus}")
        else:
            print("No viruses found in configuration.")
        return

    # Check if this is a multi-segment virus
    if virus_processor.is_multi_segment_virus():
        print(f"Detected multi-segment virus: {args.virus}")
        segments = virus_processor.get_segments()
        print(f"Available segments: {', '.join(segments)}")

        # If segment is specified, process only that segment
        if args.segment:
            if args.segment not in segments:
                print(f"Error: Segment '{args.segment}' not found. Available segments: {', '.join(segments)}")
                return
            segments_to_process = [args.segment]
        else:
            # Process all segments
            segments_to_process = segments

        # Process each segment
        for segment in segments_to_process:
            print(f"\n=== Processing segment: {segment} ===")
            _process_segment(virus_processor, segment, args)

        print(f"\nProcessing complete. Results saved in: {virus_processor.result_path}")

    else:
        # Single-segment virus (traditional processing)
        print(f"Processing single-segment virus: {args.virus}")
        _process_single_segment(virus_processor, args)


def _process_segment(virus_processor: VirusMutationProcessor, segment: str, args) -> None:
    """
    Process a single segment of a multi-segment virus.
    
    Args:
        virus_processor: Virus processor instance
        segment: Segment name to process
        args: Command line arguments
    """
    # Get MSA file for this segment
    msa_file = args.msa_file or virus_processor.get_default_msa_file(segment)

    # Check if MSA file exists
    clustalw_path = virus_processor.get_default_msa_file(segment).rsplit('/', 1)[0]  # Get directory from MSA file path
    msa_file_path = virus_processor.get_default_msa_file(segment)  # Use the full path directly
    if not os.path.exists(msa_file_path):
        print(f"MSA file not found: {msa_file_path}")
        print(f"Available files in {clustalw_path}:")
        if os.path.exists(clustalw_path):
            for file in os.listdir(clustalw_path):
                print(f"  - {file}")
        return

    # Initialize components for this segment
    ref_genome_path = virus_processor.get_reference_genome_path(segment)
    proteome_path = virus_processor.get_proteome_path(segment)

    ref_genome = ReferenceGenome(ref_genome_path)
    proteome = Proteome(proteome_path)
    msa_instance = MultipleSequenceAlignment(msa_file_path)

    # Create processors
    gene_mut_detector = GeneMutationDetector(ref_genome, msa_instance, args.genome_id)
    orf_processor = ORFProcessor(ref_genome)
    aa_mutation_processor = AminoAcidMutationProcessor()
    alignment_processor = AlignmentProcessor(proteome, aa_mutation_processor)
    pro_mutation_detector_instance = ProMutationDetector(
        orf_processor, alignment_processor, aa_mutation_processor, ref_genome
    )

    # Initialize genome SNP processor
    genome_snp_processor_instance = GenomeSNPProcessor(
        clustalw_path,
        msa_instance,
        gene_mut_detector,
        pro_mutation_detector_instance,
    )

    # Initialize SNP processor
    snp_processor_instance = SNPProcessor(
        genome_snp_processor_instance, virus_processor
    )

    # Initialize frameshift detector if requested
    frameshift_detector = None

    if args.detect_frameshifts:
        frameshift_detector = FrameshiftDetector(ref_genome)
        print("Frameshift detection enabled")

    msa_files = [msa_file]

    # Run SNP processing only if not doing frameshift-only analysis
    if not args.detect_frameshifts:
        if args.process_all or not args.genome_id:
            # Process all genomes in the MSA file
            print(f"Processing all genomes in segment {segment}...")
            snp_processor_instance.process_snp_records(msa_files, segment)
        else:
            # Process only one specific genome
            print(f"Processing genome: {args.genome_id} in segment {segment}")
            snp_processor_instance.get_snp_records_for_one_genome(msa_files, args.genome_id, segment)

    # Run frameshift analysis if requested
    if frameshift_detector:
        _run_additional_analyses(
            msa_instance, frameshift_detector,
            virus_processor, segment, args
        )


def _process_single_segment(virus_processor: VirusMutationProcessor, args) -> None:
    """
    Process a single-segment virus (traditional processing).
    
    Args:
        virus_processor: Virus processor instance
        args: Command line arguments
    """
    # Get MSA file
    msa_file = args.msa_file or virus_processor.get_default_msa_file()

    # Check if MSA file exists
    clustalw_path = virus_processor.get_default_msa_file().rsplit('/', 1)[0]  # Get directory from MSA file path
    msa_file_path = virus_processor.get_default_msa_file()  # Use the full path directly
    if not os.path.exists(msa_file_path):
        print(f"MSA file not found: {msa_file_path}")
        print(f"Available files in {clustalw_path}:")
        if os.path.exists(clustalw_path):
            for file in os.listdir(clustalw_path):
                print(f"  - {file}")
        return

    # Initialize components
    ref_genome_path = virus_processor.get_reference_genome_path()
    proteome_path = virus_processor.get_proteome_path()

    ref_genome = ReferenceGenome(ref_genome_path)
    proteome = Proteome(proteome_path)
    msa_instance = MultipleSequenceAlignment(msa_file_path)

    # Create processors
    gene_mut_detector = GeneMutationDetector(ref_genome, msa_instance, args.genome_id)
    orf_processor = ORFProcessor(ref_genome)
    aa_mutation_processor = AminoAcidMutationProcessor()
    alignment_processor = AlignmentProcessor(proteome, aa_mutation_processor)
    pro_mutation_detector_instance = ProMutationDetector(
        orf_processor, alignment_processor, aa_mutation_processor, ref_genome
    )

    # Initialize genome SNP processor
    genome_snp_processor_instance = GenomeSNPProcessor(
        clustalw_path,
        msa_instance,
        gene_mut_detector,
        pro_mutation_detector_instance,
    )

    # Initialize SNP processor
    snp_processor_instance = SNPProcessor(
        genome_snp_processor_instance, virus_processor
    )

    # Initialize frameshift detector if requested
    frameshift_detector = None

    if args.detect_frameshifts:
        frameshift_detector = FrameshiftDetector(ref_genome)
        print("Frameshift detection enabled")

    msa_files = [msa_file]

    # Run SNP processing only if not doing frameshift-only analysis
    if not args.detect_frameshifts:
        if args.process_all or not args.genome_id:
            # Process all genomes in the MSA file
            print("Processing all genomes in the MSA file...")
            snp_processor_instance.process_snp_records(msa_files)
        else:
            # Process only one specific genome
            print(f"Processing genome: {args.genome_id}")
            snp_processor_instance.get_snp_records_for_one_genome(msa_files, args.genome_id)

    # Run frameshift analysis if requested
    if frameshift_detector:
        _run_additional_analyses(
            msa_instance, frameshift_detector,
            virus_processor, None, args
        )

    print(f"Processing complete. Results saved in: {virus_processor.result_path}")


def _run_additional_analyses(
    msa_instance, frameshift_detector,
    virus_processor, segment, args
) -> None:
    """
    Run additional analyses (frameshifting detection).
    
    Args:
        msa_instance: Multiple sequence alignment instance
        frameshift_detector: Frameshift detector instance (or None)
        virus_processor: Virus processor instance
        segment: Segment name (or None for single-segment)
        args: Command line arguments
    """
    print("\n=== Running Additional Analyses ===")

    # Get reference sequence
    ref_seq_id, ref_seq_msa = msa_instance.find_ref_msa("reference")
    if ref_seq_id is None:
        print("Warning: Could not find reference sequence for additional analyses")
        return

    # Remove gaps from reference sequence
    ref_seq_clean = ref_seq_msa.replace('-', '')

    # Run frameshift detection
    if frameshift_detector:
        print("Detecting frameshift sites...")
        frameshift_sites = frameshift_detector.detect_frameshift_sites(ref_seq_clean, start_pos=1)

        result_path = virus_processor.get_result_path(segment)

        # Generate CSV output
        frameshift_output = frameshift_detector.format_frameshift_output(frameshift_sites)
        print(frameshift_output)

        frameshift_file = os.path.join(result_path, f"potential_PRF_{today}.csv")
        with open(frameshift_file, 'w') as f:
            f.write(frameshift_output)
        print(f"Potential PRF analysis saved to: {frameshift_file}")





def extract_mutation_summary_to_csv(virus_name: str, genome_id: str, segment: Optional[str] = None) -> None:
    """
    Extract mutation information from .txt output files and generate a comprehensive CSV file.
    
    Parameters:
    -----
    virus_name : str
        Name of the virus (e.g., 'SARS-CoV-2')
    genome_id : str
        Sample identifier (e.g., 'EPI_ISL_16327572')
    segment : Optional[str]
        Segment name for multi-segment viruses
    """
    from src.utils.mutation_summary import extract_mutation_summary_to_csv as extract_mutations
    
    # Get proteome file path from virus processor
    virus_processor = VirusMutationProcessor(virus_name)
    proteome_file = virus_processor.get_proteome_path(segment)
    
    extract_mutations(virus_name, genome_id, segment, proteome_file)


if __name__ == "__main__":
    main()
