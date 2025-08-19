# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: wang.rui@nyu.edu
FilePath: /VARIANT/./src/cli/commands.py
Description: CLI commands and argument parsing.
'''


"""
Command-line interface commands for MutParser.

This module provides the main CLI commands for running virus mutation analysis.
"""

import argparse
import os
import sys

from ..core import (
    AlignmentProcessor,
    AminoAcidMutationProcessor,
    GeneMutationDetector,
    GenomeSNPProcessor,
    MultipleSequenceAlignment,
    ORFProcessor,
    ProMutationDetector,
    Proteome,
    ReferenceGenome,
)


def setup_parser() -> argparse.ArgumentParser:
    """Set up the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="MutParser: Comprehensive Viral Mutation Analysis Framework",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  variant analyze --virus SARS-CoV-2 --msa data/SARS-CoV-2/clustalW/test_msa_1.txt
  variant setup --virus HIV-1 --config virus_config.yaml
  variant prf --genome data/SARS-CoV-2/refs/NC_045512.fasta --output results.csv
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Analyze command
    analyze_parser = subparsers.add_parser("analyze", help="Analyze virus mutations")
    analyze_parser.add_argument(
        "--virus", required=True, help="Virus name (e.g., SARS-CoV-2)"
    )
    analyze_parser.add_argument("--msa", required=True, help="Path to MSA file")
    analyze_parser.add_argument("--output", default="result", help="Output directory")
    analyze_parser.add_argument(
        "--config", default="virus_config.yaml", help="Configuration file"
    )

    # Setup command
    setup_parser = subparsers.add_parser("setup", help="Setup virus dataset")
    setup_parser.add_argument("--virus", required=True, help="Virus name")
    setup_parser.add_argument(
        "--config", default="virus_config.yaml", help="Configuration file"
    )
    # Note: Force parameter not yet implemented in setup method
    # setup_parser.add_argument(
    #     "--force", action="store_true", help="Force overwrite existing data"
    # )

    # PRF command
    prf_parser = subparsers.add_parser("prf", help="Analyze PRF sites")
    prf_parser.add_argument("--genome", required=True, help="Path to genome file")
    prf_parser.add_argument("--mutations", help="Path to mutations file (optional)")
    prf_parser.add_argument("--output", help="Output file for report")

    return parser


def run_analysis(args: argparse.Namespace) -> int:
    """Run the main mutation analysis."""
    try:
        import sys
        import os
        
        # Add the project root to the path to import from main.py
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        if project_root not in sys.path:
            sys.path.insert(0, project_root)
            
        from main import SNPProcessor, VirusMutationProcessor

        # Initialize virus processor
        virus_processor = VirusMutationProcessor(args.virus, args.config)

        # Use the provided MSA file
        if not os.path.exists(args.msa):
            print(f"Error: MSA file not found: {args.msa}")
            return 1

        print(f"Processing MSA file: {args.msa}")

        # Initialize components
        ref_genome_path = virus_processor.get_reference_genome_path()
        proteome_path = virus_processor.get_proteome_path()

        ref_genome = ReferenceGenome(ref_genome_path)
        proteome = Proteome(proteome_path)
        msa = MultipleSequenceAlignment(args.msa)

        # Create processors with MSA (use reference genome name for initialization)
        # The actual processing will create proper detectors for each genome
        gene_mut_detector = GeneMutationDetector(ref_genome, msa, ref_genome.ref_name)
        orf_processor = ORFProcessor(ref_genome)
        aa_mutation_processor = AminoAcidMutationProcessor()
        alignment_processor = AlignmentProcessor(proteome, aa_mutation_processor)
        pro_mut_detector = ProMutationDetector(
            orf_processor, alignment_processor, aa_mutation_processor, ref_genome
        )

        genome_snp_processor = GenomeSNPProcessor(
            virus_processor.data_path, msa, gene_mut_detector, pro_mut_detector
        )

        snp_processor = SNPProcessor(genome_snp_processor, virus_processor)

        # Get SNP records
        snp_records = genome_snp_processor.get_genome_snps(args.msa)

        # Process and save results using the main processor
        from main import MutationProcessor
        
        # Create a mutation processor instance
        processor = MutationProcessor(args.virus, args.config)
        
        # Process the SNP records to generate individual genome files and row hot mutations
        processor._process_root_variants(snp_records, "Combined", virus_processor, segment=None)

        print("Analysis completed successfully!")
        return 0

    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return 1


def run_setup(args: argparse.Namespace) -> int:
    """Run the setup command."""
    try:
        import sys
        import os
        
        # Add the project root to the path to import from scripts
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        if project_root not in sys.path:
            sys.path.insert(0, project_root)
            
        from scripts.setup_virus_dataset import VirusDatasetSetup

        setup = VirusDatasetSetup(args.config)
        # Note: The force parameter is not implemented in the setup method
        setup.setup_virus_dataset(args.virus)

        print(f"Setup completed for virus: {args.virus}")
        return 0

    except Exception as e:
        print(f"Error during setup: {e}")
        return 1


def run_prf_analysis(args: argparse.Namespace) -> int:
    """Run PRF analysis."""
    try:
        import sys
        import os
        from datetime import datetime
        
        # Add the project root to the path
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
        if project_root not in sys.path:
            sys.path.insert(0, project_root)
        
        from ..core.frameshift_detector import FrameshiftDetector
        from ..core.reference_genome import ReferenceGenome
        from Bio import SeqIO
        
        # Read genome file
        if not os.path.exists(args.genome):
            print(f"Error: Genome file not found: {args.genome}")
            return 1
            
        print(f"Loading genome from: {args.genome}")
        
        # Read the genome sequence
        with open(args.genome, 'r') as f:
            record = SeqIO.read(f, 'fasta')
            genome_sequence = str(record.seq)
        
        # Initialize frameshift detector
        ref_genome = ReferenceGenome(args.genome)
        frameshift_detector = FrameshiftDetector(ref_genome)
        
        print("Detecting frameshift sites...")
        frameshift_sites = frameshift_detector.detect_frameshift_sites(genome_sequence)
        
        # Generate output
        output_content = frameshift_detector.format_frameshift_output(frameshift_sites)
        print(output_content)
        
        # Save to file if output specified
        if args.output:
            with open(args.output, 'w') as f:
                f.write(output_content)
            print(f"Results saved to: {args.output}")
        else:
            # Default output file
            today = datetime.now().strftime("%Y%m%d")
            default_output = f"potential_PRF_{today}.csv"
            with open(default_output, 'w') as f:
                f.write(output_content)
            print(f"Results saved to: {default_output}")
        
        return 0
        
    except Exception as e:
        print(f"Error during PRF analysis: {e}")
        import traceback
        traceback.print_exc()
        return 1


def main() -> int:
    """Main CLI entry point."""
    parser = setup_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    if args.command == "analyze":
        return run_analysis(args)
    elif args.command == "setup":
        return run_setup(args)
    elif args.command == "prf":
        return run_prf_analysis(args)
    else:
        print(f"Unknown command: {args.command}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
