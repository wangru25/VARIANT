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
from ..utils.prf_analyzer import PRFAnalyzer


def setup_parser() -> argparse.ArgumentParser:
    """Set up the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="MutParser: Comprehensive Viral Mutation Analysis Framework",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  mutparser analyze --virus SARS-CoV-2 --msa data/SARS-CoV-2/clustalW/msa.txt
  mutparser setup --virus HIV --config virus_config.yaml
  mutparser prf --genome data/refs/NC_045512.fasta
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
    setup_parser.add_argument(
        "--force", action="store_true", help="Force overwrite existing data"
    )

    # PRF command
    prf_parser = subparsers.add_parser("prf", help="Analyze PRF sites")
    prf_parser.add_argument("--genome", required=True, help="Path to genome file")
    prf_parser.add_argument("--mutations", help="Path to mutations file (optional)")
    prf_parser.add_argument("--output", help="Output file for report")

    return parser


def run_analysis(args: argparse.Namespace) -> int:
    """Run the main mutation analysis."""
    try:
        from main import SNPProcessor, VirusMutationProcessor

        # Initialize virus processor
        virus_processor = VirusMutationProcessor(args.virus, args.config)

        # Get MSA files
        msa_files = [
            f
            for f in os.listdir(virus_processor.data_path)
            if f.endswith(".txt") and "msa" in f.lower()
        ]

        if not msa_files:
            print(f"No MSA files found in {virus_processor.data_path}")
            return 1

        # Initialize components
        ref_genome_path = virus_processor.get_reference_genome_path()
        proteome_path = virus_processor.get_proteome_path()

        ref_genome = ReferenceGenome(ref_genome_path)
        proteome = Proteome(proteome_path)

        # Create processors
        gene_mut_detector = GeneMutationDetector(ref_genome, None, None)
        orf_processor = ORFProcessor(ref_genome)
        aa_mutation_processor = AminoAcidMutationProcessor()
        alignment_processor = AlignmentProcessor(proteome, aa_mutation_processor)
        pro_mut_detector = ProMutationDetector(
            orf_processor, alignment_processor, aa_mutation_processor, ref_genome
        )

        genome_snp_processor = GenomeSNPProcessor(
            virus_processor.data_path, None, gene_mut_detector, pro_mut_detector
        )

        snp_processor = SNPProcessor(genome_snp_processor, virus_processor)

        # Process MSA files
        for msa_file in msa_files:
            print(f"Processing {msa_file}...")
            full_msa_path = os.path.join(virus_processor.data_path, msa_file)
            msa = MultipleSequenceAlignment(full_msa_path)

            # Update processors with MSA
            gene_mut_detector.msa = msa
            genome_snp_processor.msa = msa

            # Get SNP records
            snp_records = genome_snp_processor.get_genome_snps(full_msa_path)

            # Process and save results
            snp_processor._process_root_variants(snp_records, msa_file.split(".")[0])

        print("Analysis completed successfully!")
        return 0

    except Exception as e:
        print(f"Error during analysis: {e}")
        return 1


def run_setup(args: argparse.Namespace) -> int:
    """Run the setup command."""
    try:
        from scripts.setup_virus_dataset import VirusDatasetSetup

        setup = VirusDatasetSetup(args.config)
        setup.setup_virus_dataset(args.virus, force=args.force)

        print(f"Setup completed for virus: {args.virus}")
        return 0

    except Exception as e:
        print(f"Error during setup: {e}")
        return 1


def run_prf_analysis(args: argparse.Namespace) -> int:
    """Run PRF analysis."""
    try:
        from ..utils.sequence_utils import read_fasta

        # Read genome
        records = read_fasta(args.genome)
        genome_sequence = records[0]["seq"]

        # Initialize PRF analyzer
        analyzer = PRFAnalyzer(genome_sequence)
        prf_site = analyzer.detect_prf_site()

        # Load mutations if provided
        mutations = []
        if args.mutations:
            # This would need to be implemented based on your mutation file format
            pass

        # Generate report
        report = analyzer.generate_prf_report(prf_site, mutations)

        # Output results
        if args.output:
            import json

            with open(args.output, "w") as f:
                json.dump(report, f, indent=2)
            print(f"Report saved to {args.output}")
        else:
            print("PRF Analysis Results:")
            print(f"  Position: {prf_site.position}")
            print(f"  Sequence: {prf_site.sequence}")
            print(f"  Type: {prf_site.site_type.value}")
            print(f"  Efficiency: {prf_site.efficiency_score:.3f}")
            print(f"  Notes: {prf_site.notes}")

        return 0

    except Exception as e:
        print(f"Error during PRF analysis: {e}")
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
