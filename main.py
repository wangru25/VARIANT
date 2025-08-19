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

from src.core.frameshift_detector import FrameshiftDetector
from src.core.genome_processor import GenomeSNPProcessor
from src.core.mutation_detector import GeneMutationDetector, MultipleSequenceAlignment
from src.core.protein_analyzer import (
    AlignmentProcessor, AminoAcidMutationProcessor, ORFProcessor, 
    ProMutationDetector, Proteome
)
from src.core.reference_genome import ReferenceGenome
from src.core.snp_processor import SNPProcessor
from src.core.virus_processor import VirusMutationProcessor
from src.utils.mutation_summary import extract_mutation_summary_to_csv


# Constants
DATE = "20250807"
TODAY = 20250807


class MutationProcessor:
    """
    Main processor class that handles the complete mutation analysis pipeline.
    """
    
    def __init__(self, virus_name: str, config_file: str = "virus_config.yaml"):
        """Initialize the mutation processor."""
        self.virus_processor = VirusMutationProcessor(virus_name, config_file)
        self.virus_name = virus_name
        
    def list_available_viruses(self) -> List[str]:
        """List all available viruses in the configuration."""
        return self.virus_processor.list_available_viruses()
    
    def process(self, args) -> None:
        """
        Main processing method that handles both single and multi-segment viruses.
        
        Args:
            args: Command line arguments
        """
        if self.virus_processor.is_multi_segment_virus():
            self._process_multi_segment_virus(args)
        else:
            self._process_single_segment_virus(args)
    
    def _process_multi_segment_virus(self, args) -> None:
        """Process multi-segment virus."""
        print(f"Detected multi-segment virus: {self.virus_name}")
        segments = self.virus_processor.get_segments()
        print(f"Available segments: {', '.join(segments)}")
        
        # Determine which segments to process
        if args.segment:
            if args.segment not in segments:
                print(f"Error: Segment '{args.segment}' not found. Available segments: {', '.join(segments)}")
                return
            segments_to_process = [args.segment]
        else:
            segments_to_process = segments
        
        # Process each segment
        for segment in segments_to_process:
            print(f"\n=== Processing segment: {segment} ===")
            self._process_segment(segment, args)
        
        print(f"\nProcessing complete. Results saved in: {self.virus_processor.result_path}")
    
    def _process_single_segment_virus(self, args) -> None:
        """Process single-segment virus."""
        print(f"Processing single-segment virus: {self.virus_name}")
        self._process_segment(None, args)
        print(f"Processing complete. Results saved in: {self.virus_processor.result_path}")
    
    def _process_segment(self, segment: Optional[str], args) -> None:
        """
        Process a single segment (or entire virus for single-segment).
        
        Args:
            segment: Segment name (None for single-segment viruses)
            args: Command line arguments
        """
        # Setup paths and validate files
        setup_result = self._setup_processing_environment(segment, args)
        if not setup_result:
            return
            
        msa_file_path, clustalw_path = setup_result
        
        # Initialize all processors
        processors = self._initialize_processors(segment, args, msa_file_path, clustalw_path)
        if not processors:
            return
        
        # Run the processing pipeline
        self._run_processing_pipeline(processors, segment, args, [args.msa_file or msa_file_path])
    
    def _setup_processing_environment(self, segment: Optional[str], args) -> Optional[Tuple[str, str]]:
        """
        Setup and validate the processing environment.
        
        Returns:
            Tuple of (msa_file_path, clustalw_path) if successful, None otherwise
        """
        # Get MSA file path
        if segment:
            msa_file_path = args.msa_file or self.virus_processor.get_default_msa_file(segment)
        else:
            msa_file_path = args.msa_file or self.virus_processor.get_default_msa_file()
        
        # Validate MSA file exists
        if not os.path.exists(msa_file_path):
            clustalw_path = msa_file_path.rsplit('/', 1)[0]
            print(f"MSA file not found: {msa_file_path}")
            print(f"Available files in {clustalw_path}:")
            if os.path.exists(clustalw_path):
                for file in os.listdir(clustalw_path):
                    print(f"  - {file}")
            return None
        
        clustalw_path = msa_file_path.rsplit('/', 1)[0]
        return msa_file_path, clustalw_path
    
    def _initialize_processors(self, segment: Optional[str], args, msa_file_path: str, clustalw_path: str) -> Optional[Dict]:
        """
        Initialize all required processors.
        
        Returns:
            Dictionary of processors if successful, None otherwise
        """
        try:
            # Get file paths
            if segment:
                ref_genome_path = self.virus_processor.get_reference_genome_path(segment)
                proteome_path = self.virus_processor.get_proteome_path(segment)
            else:
                ref_genome_path = self.virus_processor.get_reference_genome_path()
                proteome_path = self.virus_processor.get_proteome_path()
            
            # Initialize core components
            ref_genome = ReferenceGenome(ref_genome_path)
            proteome = Proteome(proteome_path)
            msa_instance = MultipleSequenceAlignment(msa_file_path)
            
            # Initialize mutation detectors
            gene_mut_detector = GeneMutationDetector(ref_genome, msa_instance, args.genome_id)
            orf_processor = ORFProcessor(ref_genome)
            aa_mutation_processor = AminoAcidMutationProcessor()
            alignment_processor = AlignmentProcessor(proteome, aa_mutation_processor)
            pro_mutation_detector = ProMutationDetector(
                orf_processor, alignment_processor, aa_mutation_processor, ref_genome
            )
            
            # Initialize processing pipeline
            genome_snp_processor = GenomeSNPProcessor(
                clustalw_path, msa_instance, gene_mut_detector, pro_mutation_detector
            )
            snp_processor = SNPProcessor(genome_snp_processor, self.virus_processor)
            
            # Initialize frameshift detector if needed
            frameshift_detector = None
            if args.detect_frameshifts:
                frameshift_detector = FrameshiftDetector(ref_genome)
                print("Frameshift detection enabled")
            
            return {
                'ref_genome': ref_genome,
                'proteome': proteome,
                'msa_instance': msa_instance,
                'gene_mut_detector': gene_mut_detector,
                'snp_processor': snp_processor,
                'frameshift_detector': frameshift_detector
            }
            
        except Exception as e:
            print(f"Error initializing processors: {e}")
            return None
    
    def _run_processing_pipeline(self, processors: Dict, segment: Optional[str], args, msa_files: List[str]) -> None:
        """
        Run the complete processing pipeline.
        
        Args:
            processors: Dictionary of initialized processors
            segment: Segment name (None for single-segment)
            args: Command line arguments
            msa_files: List of MSA files to process
        """
        # Step 1: SNP Processing
        if not args.detect_frameshifts:
            self._run_snp_processing(processors['snp_processor'], segment, args, msa_files)
        
        # Step 2: Frameshift Analysis (if requested)
        if processors['frameshift_detector']:
            self._run_frameshift_analysis(processors, segment, args)
        
        # Step 3: Generate Mutation Summaries
        if not args.detect_frameshifts:
            self._generate_mutation_summaries(segment, args)
    
    def _run_snp_processing(self, snp_processor: SNPProcessor, segment: Optional[str], args, msa_files: List[str]) -> None:
        """Run SNP processing."""
        if args.process_all or not args.genome_id:
            # Process all genomes
            if segment:
                print(f"Processing all genomes in segment {segment}...")
                snp_processor.process_snp_records(msa_files, segment)
            else:
                print("Processing all genomes in the MSA file...")
                snp_processor.process_snp_records(msa_files)
        else:
            # Process specific genome
            if segment:
                print(f"Processing genome: {args.genome_id} in segment {segment}")
                snp_processor.get_snp_records_for_one_genome(msa_files, args.genome_id, segment)
            else:
                print(f"Processing genome: {args.genome_id}")
                snp_processor.get_snp_records_for_one_genome(msa_files, args.genome_id)
    
    def _run_frameshift_analysis(self, processors: Dict, segment: Optional[str], args) -> None:
        """Run frameshift analysis."""
        print("\n=== Running Frameshift Analysis ===")
        
        msa_instance = processors['msa_instance']
        frameshift_detector = processors['frameshift_detector']
        
        # Get reference sequence
        ref_seq_id, ref_seq_msa = msa_instance.find_ref_msa("reference")
        if ref_seq_id is None:
            print("Warning: Could not find reference sequence for frameshift analysis")
            return
        
        # Remove gaps and run detection
        ref_seq_clean = ref_seq_msa.replace('-', '')
        frameshift_sites = frameshift_detector.detect_frameshift_sites(ref_seq_clean, start_pos=1)
        
        # Save results
        if segment:
            result_path = self.virus_processor.get_result_path(segment)
        else:
            result_path = self.virus_processor.get_result_path()
            
        frameshift_output = frameshift_detector.format_frameshift_output(frameshift_sites)
        frameshift_file = os.path.join(result_path, f"potential_PRF_{TODAY}.csv")
        
        with open(frameshift_file, 'w') as f:
            f.write(frameshift_output)
        
        print(f"Frameshift analysis saved to: {frameshift_file}")
    
    def _generate_mutation_summaries(self, segment: Optional[str], args) -> None:
        """Generate mutation summaries for processed genomes."""
        print("\n=== Generating Mutation Summaries ===")
        
        if args.process_all or not args.genome_id:
            self._generate_summaries_for_all_genomes(segment, args)
        else:
            self._generate_summary_for_single_genome(segment, args)
    
    def _generate_summaries_for_all_genomes(self, segment: Optional[str], args) -> None:
        """Generate mutation summaries for all processed genomes."""
        print("Generating mutation summaries for all processed genomes...")
        
        # Get result directory
        if segment:
            result_path = self.virus_processor.get_result_path(segment)
        else:
            result_path = self.virus_processor.get_result_path()
        
        # Find all .txt files
        txt_pattern = os.path.join(result_path, f"*_{TODAY}.txt")
        txt_files = glob.glob(txt_pattern)
        
        if not txt_files:
            print(f"No .txt files found in {result_path}")
            return
        
        print(f"Found {len(txt_files)} .txt files to process")
        
        # Process each file
        processed_count = 0
        for txt_file in txt_files:
            genome_id = os.path.basename(txt_file).replace(f"_{TODAY}.txt", "")
            csv_file = os.path.join(result_path, f"{genome_id}_mutation_summary_{TODAY}.csv")
            
            if os.path.exists(csv_file):
                print(f"  ✓ Mutation summary already exists for {genome_id}")
                continue
            
            try:
                print(f"  Generating mutation summary for {genome_id}...")
                self._extract_mutation_summary(self.virus_name, genome_id, segment)
                processed_count += 1
            except Exception as e:
                print(f"  ⚠ Error generating mutation summary for {genome_id}: {e}")
        
        print(f"✓ Generated mutation summaries for {processed_count} genomes")
    
    def _generate_summary_for_single_genome(self, segment: Optional[str], args) -> None:
        """Generate mutation summary for a single genome."""
        print(f"Generating mutation summary for genome: {args.genome_id}")
        
        # Check if required .txt file exists
        if segment:
            result_path = self.virus_processor.get_result_path(segment)
        else:
            result_path = self.virus_processor.get_result_path()
            
        txt_file = os.path.join(result_path, f"{args.genome_id}_{TODAY}.txt")
        
        if not os.path.exists(txt_file):
            print(f"⚠ Skipping mutation summary generation: {txt_file} not found")
            print("  This may indicate that SNP processing did not complete successfully")
            return
        
        try:
            self._extract_mutation_summary(self.virus_name, args.genome_id, segment)
            print("✓ Mutation summary generated successfully")
        except Exception as e:
            print(f"⚠ Error generating mutation summary: {e}")
    
    def _extract_mutation_summary(self, virus_name: str, genome_id: str, segment: Optional[str]) -> None:
        """Extract mutation summary to CSV file."""
        proteome_file = self.virus_processor.get_proteome_path(segment)
        extract_mutation_summary_to_csv(virus_name, genome_id, segment, proteome_file)


def create_argument_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        description="VARIANT - Virus Mutation Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --virus SARS-CoV-2 --genome-id EPI_ISL_16327572
  %(prog)s --virus SARS-CoV-2 --process-all
  %(prog)s --virus H3N2 --segment segment_1 --genome-id sample123
  %(prog)s --list-viruses
        """
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
    
    # Handle list viruses request
    if args.list_viruses:
        available_viruses = processor.list_available_viruses()
        if available_viruses:
            print("Available viruses:")
            for virus in available_viruses:
                print(f"  - {virus}")
        else:
            print("No viruses found in configuration.")
        return 0
    
    # Validate arguments
    if not args.process_all and not args.genome_id:
        print("Error: Either --process-all or --genome-id must be specified")
        return 1
    
    # Run processing
    try:
        processor.process(args)
        return 0
    except Exception as e:
        print(f"Error during processing: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
