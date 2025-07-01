#!/usr/bin/env python3
"""
Example script demonstrating PRF analysis in SARS-CoV-2.

This script shows how to use the PRF analyzer to:
1. Detect PRF sites in viral genomes
2. Calculate frameshifting efficiency
3. Analyze mutations affecting PRF
4. Predict protein ratios
"""

import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import json

from src.utils.prf_analyzer import PRFAnalyzer, analyze_prf_in_sequence
from src.utils.sequence_utils import read_fasta


def analyze_reference_genome():
    """Analyze PRF in the reference SARS-CoV-2 genome."""
    print("=== PRF Analysis of Reference SARS-CoV-2 Genome ===\n")

    # Read reference genome
    try:
        ref_fasta = "./data/refs/NC_045512.fasta"
        records = read_fasta(ref_fasta)
        ref_genome = records[0]["seq"]
        print(f"Reference genome loaded: {len(ref_genome)} nucleotides")
        print(f"Genome ID: {records[0]['seqId']}\n")
    except FileNotFoundError:
        print("Reference genome not found. Using simulated data...")
        # Create simulated reference genome with canonical PRF site
        ref_genome = "A" * 13467 + "UUUAAAC" + "GGC" * 20 + "A" * 1000

    # Analyze PRF
    analyzer = PRFAnalyzer(ref_genome)
    prf_site = analyzer.detect_prf_site()

    print("PRF Site Analysis:")
    print(f"  Position: {prf_site.position}")
    print(f"  Sequence: {prf_site.sequence}")
    print(f"  Type: {prf_site.site_type.value}")
    print(f"  Efficiency Score: {prf_site.efficiency_score:.3f}")
    print(f"  Downstream Structure: {prf_site.downstream_structure}\n")

    # Predict protein ratios
    ratios = analyzer.predict_protein_ratios(prf_site.efficiency_score)
    print("Predicted Protein Ratios:")
    print(
        f"  ORF1a:ORF1ab ratio: {ratios['orf1a_ratio']:.3f}:{ratios['orf1ab_ratio']:.3f}"
    )
    print(f"  RdRp availability: {ratios['rdrp_availability']:.3f}")
    print(f"  Replication efficiency: {ratios['replication_efficiency']:.3f}\n")

    return analyzer, prf_site


def analyze_mutations():
    """Analyze mutations affecting PRF."""
    print("=== PRF Mutation Analysis ===\n")

    # Example mutations that might affect PRF
    mutations = [
        {
            "pos": "13468",
            "type": "pointMutation",
            "SNP": "U->A",
            "description": "Mutation in slippery site",
        },
        {
            "pos": "13470:13472",
            "type": "deletion",
            "SNP": "UAA",
            "description": "Deletion in slippery site",
        },
        {
            "pos": "13480:13485",
            "type": "insertion",
            "SNP": "GGGGGG",
            "description": "Insertion in downstream region",
        },
        {
            "pos": "20000",
            "type": "pointMutation",
            "SNP": "A->G",
            "description": "Mutation outside PRF region",
        },
    ]

    # Create analyzer with reference genome
    ref_genome = "A" * 13467 + "UUUAAAC" + "GGC" * 20 + "A" * 1000
    analyzer = PRFAnalyzer(ref_genome)

    # Analyze PRF mutations
    prf_mutations = analyzer.analyze_prf_mutations(mutations)

    print(f"Total mutations: {len(mutations)}")
    print(f"Mutations affecting PRF: {len(prf_mutations)}\n")

    for i, mutation in enumerate(prf_mutations, 1):
        print(f"PRF Mutation {i}:")
        print(f"  Position: {mutation['pos']}")
        print(f"  Type: {mutation['type']}")
        print(f"  SNP: {mutation['SNP']}")
        print(f"  Region: {mutation['prf_region']}")
        print(f"  Impact: {mutation['prf_impact']['severity']}")
        print(f"  Efficiency change: {mutation['prf_impact']['efficiency_change']:.3f}")
        print(f"  Mechanism: {mutation['prf_impact']['mechanism']}")
        print()


def analyze_variant_genomes():
    """Analyze PRF in variant genomes."""
    print("=== PRF Analysis in Variant Genomes ===\n")

    # Example variant genomes with different PRF sites
    variants = {
        "Alpha": "A" * 13467
        + "UUUAAAA"
        + "GGC" * 20
        + "A" * 1000,  # Variant slippery site
        "Beta": "A" * 13467
        + "UUUAAAC"
        + "AAA" * 20
        + "A" * 1000,  # Different downstream
        "Gamma": "A" * 13467 + "AAAAAAA" + "GGC" * 20 + "A" * 1000,  # No PRF site
    }

    for variant_name, genome in variants.items():
        print(f"{variant_name} Variant:")

        try:
            analyzer = PRFAnalyzer(genome)
            prf_site = analyzer.detect_prf_site()

            print(f"  PRF Site: {prf_site.sequence}")
            print(f"  Type: {prf_site.site_type.value}")
            print(f"  Efficiency: {prf_site.efficiency_score:.3f}")

            ratios = analyzer.predict_protein_ratios(prf_site.efficiency_score)
            print(f"  ORF1ab ratio: {ratios['orf1ab_ratio']:.3f}")
            print(f"  Replication efficiency: {ratios['replication_efficiency']:.3f}")

        except Exception as e:
            print(f"  Error: {e}")

        print()


def generate_comprehensive_report():
    """Generate a comprehensive PRF analysis report."""
    print("=== Comprehensive PRF Analysis Report ===\n")

    # Create reference genome
    ref_genome = "A" * 13467 + "UUUAAAC" + "GGC" * 20 + "A" * 1000

    # Example mutations
    mutations = [
        {"pos": "13468", "type": "pointMutation", "SNP": "U->A"},
        {"pos": "13470:13472", "type": "deletion", "SNP": "UAA"},
        {"pos": "20000", "type": "pointMutation", "SNP": "A->G"},
    ]

    # Generate report
    report = analyze_prf_in_sequence(ref_genome, mutations)

    print("Report Summary:")
    print(f"  Total mutations analyzed: {report['summary']['total_mutations']}")
    print(f"  PRF-affecting mutations: {report['summary']['prf_affecting_mutations']}")
    print(f"  Efficiency category: {report['summary']['efficiency_category']}")

    print("\nDetailed Report:")
    print(json.dumps(report, indent=2))


def main():
    """Main function to run all PRF analysis examples."""
    print("SARS-CoV-2 Programmed Ribosomal Frameshifting (PRF) Analysis")
    print("=" * 60)

    try:
        # Analyze reference genome
        analyzer, prf_site = analyze_reference_genome()

        # Analyze mutations
        analyze_mutations()

        # Analyze variant genomes
        analyze_variant_genomes()

        # Generate comprehensive report
        generate_comprehensive_report()

        print("\n=== Analysis Complete ===")
        print("The PRF analyzer successfully identified and analyzed")
        print("programmed ribosomal frameshifting sites in SARS-CoV-2.")

    except Exception as e:
        print(f"Error during analysis: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
