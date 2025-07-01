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
from pathlib import Path
from datetime import date
from typing import Dict, List, Optional

import yaml

# Add src to path
sys.path.append(str(Path(__file__).parent / "src"))

from src.core.genome_processor import GenomeSNPProcessor
from src.core.reference_genome import ReferenceGenome
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
from src.utils.mutation_utils import (
    convert_key_to_int,
    sort_dict_by_consecutive_keys,
    classify_mutation_type,
    detect_hot_mutation,
    split_multi_protein_mutations,
    format_amino_acid_change_for_csv,
)

# Constants
DATE_FORMAT = "%Y%m%d"

# Global Variables
today = date.today().strftime("%Y%m%d")


class VirusMutationProcessor:
    """
    Class for processing virus-specific mutations with organized file structure.
    """

    def __init__(self, virus_name: str, config_file: str = "virus_config.yaml"):
        """
        Initialize virus mutation processor.

        Args:
            virus_name (str): Name of the virus to process.
            config_file (str): Path to the virus configuration file.
        """
        self.virus_name = virus_name
        self.config_file = config_file
        self.config = self._load_config()

        # Create virus-specific directories
        self.base_path = os.path.join("data", virus_name)
        self.data_path = os.path.join(self.base_path, "clustalW")
        self.refs_path = os.path.join(self.base_path, "refs")
        self.result_path = os.path.join("result", virus_name)

        # Ensure directories exist
        os.makedirs(self.data_path, exist_ok=True)
        os.makedirs(self.refs_path, exist_ok=True)
        os.makedirs(self.result_path, exist_ok=True)

    def _load_config(self) -> Dict:
        """
        Load virus configuration from YAML file.

        Returns:
            Dict: Configuration dictionary.
        """
        try:
            with open(self.config_file, "r") as file:
                return yaml.safe_load(file)
        except FileNotFoundError:
            print(f"Configuration file {self.config_file} not found. Using defaults.")
            return {}

    def get_virus_config(self) -> Dict[str, str]:
        """
        Get virus-specific configuration including reference files.

        Returns:
            Dict[str, str]: Configuration dictionary with file paths
        """
        # Default configuration
        default_config = {
            "reference_genome": f"{self.virus_name}_reference.fasta",
            "proteome_file": f"{self.virus_name}_proteome.fasta",
            "codon_table_id": "1",  # Standard genetic code
            "default_msa_file": f"{self.virus_name}_msa.txt",
        }

        # Try to get configuration from YAML file
        if (
            self.config
            and "viruses" in self.config
            and self.virus_name in self.config["viruses"]
        ):
            virus_config = self.config["viruses"][self.virus_name]
            # Merge with defaults for missing values
            for key, value in default_config.items():
                if key not in virus_config:
                    virus_config[key] = value
            return virus_config

        print(f"Warning: Virus '{self.virus_name}' not found in configuration.")
        print(f"Using default configuration. Consider adding it to {self.config_file}")
        return default_config

    def get_reference_genome_path(self) -> str:
        """
        Get the path to the reference genome file.

        Returns:
            str: Path to the reference genome file
        """
        config = self.get_virus_config()
        return os.path.join(self.refs_path, config["reference_genome"])

    def get_proteome_path(self) -> str:
        """
        Get the path to the proteome file.

        Returns:
            str: Path to the proteome file
        """
        config = self.get_virus_config()
        return os.path.join(self.refs_path, config["proteome_file"])

    def get_codon_table_id(self) -> int:
        """
        Get the codon table ID for the virus.

        Returns:
            int: Codon table ID
        """
        config = self.get_virus_config()
        return int(config.get("codon_table_id", 1))

    def get_default_msa_file(self) -> str:
        """
        Get the default MSA file name for the virus.

        Returns:
            str: Default MSA file name
        """
        config = self.get_virus_config()
        return config.get("default_msa_file", f"{self.virus_name}_msa.txt")

    def list_available_viruses(self) -> List[str]:
        """
        List all available viruses in the configuration.

        Returns:
            List[str]: List of available virus names
        """
        if self.config and "viruses" in self.config:
            return list(self.config["viruses"].keys())
        return []


class SNPProcessor:
    """
    Class for processing SNP records from MSA files.
    """

    def __init__(
        self,
        genome_snp_processor: GenomeSNPProcessor,
        virus_processor: VirusMutationProcessor,
    ):
        self.genome_snp_processor = genome_snp_processor
        self.virus_processor = virus_processor

    def get_snp_records_from_one_msa(
        self, msa_files: List[str], index: int
    ) -> List[Dict]:
        """
        Retrieve SNP records from a single MSA file.

        Args:
            msa_files (List[str]): List of MSA file names.
            index (int): Index of the MSA file in the list.

        Returns:
            List[Dict]: List of SNP records.
        """
        msa_file = msa_files[index]
        msa_file_prefix = msa_file.split(".")[0]
        full_msa_path = os.path.join(self.virus_processor.data_path, msa_file)
        snp_records = self.genome_snp_processor.get_genome_snps(full_msa_path)

        output_file = os.path.join(
            self.virus_processor.data_path, f"snpRecords_{msa_file_prefix}.txt"
        )
        self.genome_snp_processor.print_snp_records(snp_records, file_name=output_file)

        return snp_records

    def process_snp_records(self, msa_files: List[str]) -> None:
        """
        Process SNP records from multiple MSA files.

        Args:
            msa_files (List[str]): List of MSA file names.
        """
        all_snp_records = []
        for msa_file in msa_files:
            full_msa_path = os.path.join(self.virus_processor.data_path, msa_file)
            snp_records = self.genome_snp_processor.get_genome_snps(full_msa_path)
            all_snp_records.extend(snp_records)

        msa_file_prefix = f"Combined_{today}"
        self._write_and_print_snp_records(all_snp_records, msa_file_prefix)
        self._process_root_variants(all_snp_records, msa_file_prefix)

    def _write_and_print_snp_records(
        self, snp_records: List[Dict], prefix: str
    ) -> None:
        """
        Write and print SNP records.

        Args:
            snp_records (List[Dict]): List of SNP records.
            prefix (str): Prefix for the output file name.
        """
        # Write SNP records to a text file
        output_file_txt = os.path.join(
            self.virus_processor.data_path, f"snpRecords_{prefix}.txt"
        )
        self.genome_snp_processor.print_snp_records(
            snp_records, file_name=output_file_txt
        )

    def _process_root_variants(self, snp_records: List[Dict], prefix: str) -> None:
        """
        Process and save individual mutation records for each genome.

        Args:
            snp_records (List[Dict]): List of SNP records.
            prefix (str): Prefix for the output file name.
        """
        # Generate individual mutation files for each genome
        for snp_record in snp_records:
            genome_id = snp_record.get("seqId", "unknown")

            # Create individual mutation dictionary for this genome
            individual_mutations = {}
            row_hot_mutations = []  # Track row and hot mutations for CSV
            
            for genome_snpp in snp_record.get("genomeSNPPs", []):
                # Only filter out specific artifacts: deletions at the very beginning (1:1, 1:2, 1:3) 
                # with "Invalid protein sequence" that represent missing sequence data
                pos = genome_snpp.get('pos', '')
                protein_mutation = genome_snpp.get('proteinMutation', [])
                
                # Very specific filter: only remove deletions at positions 1:1, 1:2, 1:3 with invalid protein
                should_filter = False
                if genome_snpp.get('type') == 'deletion' and protein_mutation:
                    # Only filter positions 1:1, 1:2, 1:3 (very beginning artifacts)
                    if pos in ['1:1', '1:2', '1:3']:
                        # Check if any protein mutation has "Invalid protein sequence"
                        for protein_info in protein_mutation:
                            if isinstance(protein_info, dict) and protein_info.get('protein') == 'Invalid protein sequence':
                                should_filter = True
                                break
                
                if should_filter:
                    continue
                
                # Split multi-protein mutations
                separated_mutations = split_multi_protein_mutations(genome_snpp)
                
                for separated_mutation in separated_mutations:
                    # Classify mutation type
                    biological_type = classify_mutation_type(separated_mutation.get('proteinMutation', []))
                    
                    # Create mutation string with biological classification
                    mutation = f"{biological_type} {separated_mutation.get('pos', '')} {separated_mutation.get('SNP', '')} {separated_mutation.get('proteinMutation', '')}"
                    dict_key = convert_key_to_int(str(separated_mutation.get("pos", "")))
                    individual_mutations[dict_key] = mutation
                    
                    # Track row and hot mutations for CSV
                    original_type = separated_mutation.get('type', '')
                    if original_type in ['rowMutation', 'hotMutation']:
                        # Extract protein and amino acid information
                        protein_info = separated_mutation.get('proteinMutation', [])
                        if protein_info and isinstance(protein_info[0], dict):
                            # For row and hot mutations, process all proteins
                            mutation_type = separated_mutation.get('type', '')
                            if mutation_type in ['rowMutation', 'hotMutation']:
                                # Process each protein separately for CSV
                                for protein_data in protein_info:
                                    if isinstance(protein_data, dict):
                                        protein_name = protein_data.get('protein', 'unknown')
                                        amino_acid_change = protein_data.get('mutation', 'NA')
                                        biological_type = classify_mutation_type([protein_data])
                                        
                                        # Always check for hot mutation if SNP is X->Y and length > 1
                                        snp_info = separated_mutation.get('SNP', '')
                                        final_mutation_type = None
                                        if '->' in snp_info:
                                            original_seq, mutated_seq = snp_info.split('->')
                                            if len(original_seq) > 1 and len(mutated_seq) > 1:
                                                # For hot mutation detection, we need to extract the full sequence
                                                # from the reference genome when we have position ranges
                                                pos_str = separated_mutation.get('pos', '')
                                                if ':' in pos_str:
                                                    # Extract position range (e.g., "10447:10449")
                                                    start_pos, end_pos = map(int, pos_str.split(':'))
                                                    # Get the full sequence from reference genome (1-based to 0-based)
                                                    ref_genome = ReferenceGenome("data/refs/NC_045512.fasta")
                                                    if len(original_seq) == 2 and len(mutated_seq) == 2:
                                                        # For 2-nucleotide changes like GC->AA, check if it's a hot mutation
                                                        # by looking at the 3-nucleotide context
                                                        if start_pos + 1 <= len(ref_genome.ref_seq):
                                                            # Get the 3-nucleotide sequence: start_pos, start_pos+1, end_pos
                                                            # For position 10447:10449, we want positions 10447, 10448, 10449
                                                            full_original_seq = ref_genome.ref_seq[start_pos-1:end_pos]
                                                            # Construct mutated sequence: first_nt->first_mut, middle_nt (conserved), last_nt->last_mut
                                                            full_mutated_seq = mutated_seq[0] + full_original_seq[1] + mutated_seq[1]
                                                            is_hot = detect_hot_mutation(full_original_seq, full_mutated_seq)
                                                            final_mutation_type = 'hot' if is_hot else 'row'
                                                        else:
                                                            final_mutation_type = 'row'
                                                    else:
                                                        # For other cases, use the original detection logic
                                                        is_hot = detect_hot_mutation(original_seq, mutated_seq)
                                                        final_mutation_type = 'hot' if is_hot else 'row'
                                                else:
                                                    # For single position mutations, use the original detection logic
                                                    is_hot = detect_hot_mutation(original_seq, mutated_seq)
                                                    final_mutation_type = 'hot' if is_hot else 'row'
                                        
                                        if not final_mutation_type:
                                            final_mutation_type = mutation_type if mutation_type in ['rowMutation', 'hotMutation'] else 'row'
                                        
                                        row_hot_mutations.append({
                                            'genome_id': genome_id,
                                            'mutation_type': final_mutation_type,
                                            'position': separated_mutation.get('pos', ''),
                                            'nucleotide_change': separated_mutation.get('SNP', ''),
                                            'protein_affected': protein_name,
                                            'amino_acid_change': format_amino_acid_change_for_csv(amino_acid_change),
                                            'biological_classification': biological_type
                                        })
                        else:
                            # For other mutation types, use the original logic
                            protein_name = protein_info[0].get('protein', 'unknown')
                            amino_acid_change = protein_info[0].get('mutation', 'NA')
                            biological_type = classify_mutation_type(protein_info)
                            # Always check for hot mutation if SNP is X->Y and length > 1
                            snp_info = separated_mutation.get('SNP', '')
                            mutation_type = None
                            if '->' in snp_info:
                                original_seq, mutated_seq = snp_info.split('->')
                                if len(original_seq) > 1 and len(mutated_seq) > 1:
                                    # For hot mutation detection, we need to extract the full sequence
                                    # from the reference genome when we have position ranges
                                    pos_str = separated_mutation.get('pos', '')
                                    if ':' in pos_str:
                                        # Extract position range (e.g., "10447:10449")
                                        start_pos, end_pos = map(int, pos_str.split(':'))
                                        # Get the full sequence from reference genome (1-based to 0-based)
                                        ref_genome = ReferenceGenome("data/refs/NC_045512.fasta")
                                        full_original_seq = ref_genome.ref_seq[start_pos-1:end_pos]
                                        # Construct the full mutated sequence
                                        if len(original_seq) == 2 and len(mutated_seq) == 2:
                                            # For 2-nucleotide changes like GC->AA, check if it's a hot mutation
                                            # by looking at the 3-nucleotide context
                                            if start_pos + 1 <= len(ref_genome.ref_seq):
                                                # Get the 3-nucleotide sequence: start_pos, start_pos+1, end_pos
                                                # For position 10447:10449, we want positions 10447, 10448, 10449
                                                full_original_seq = ref_genome.ref_seq[start_pos-1:end_pos]
                                                # Construct mutated sequence: first_nt->first_mut, middle_nt (conserved), last_nt->last_mut
                                                full_mutated_seq = mutated_seq[0] + full_original_seq[1] + mutated_seq[1]
                                                is_hot = detect_hot_mutation(full_original_seq, full_mutated_seq)
                                                mutation_type = 'hot' if is_hot else 'row'
                                            else:
                                                is_hot = detect_hot_mutation(original_seq, mutated_seq)
                                                mutation_type = 'hot' if is_hot else 'row'
                                        else:
                                            is_hot = detect_hot_mutation(original_seq, mutated_seq)
                                            mutation_type = 'hot' if is_hot else 'row'
                            if not mutation_type:
                                mutation_type = original_type if original_type in ['rowMutation', 'hotMutation'] else 'row'
                            row_hot_mutations.append({
                                'genome_id': genome_id,
                                'mutation_type': mutation_type,
                                'position': separated_mutation.get('pos', ''),
                                'nucleotide_change': separated_mutation.get('SNP', ''),
                                'protein_affected': protein_name,
                                'amino_acid_change': format_amino_acid_change_for_csv(amino_acid_change),
                                'biological_classification': biological_type
                            })

            # Sort mutations by position
            sorted_mutations = sort_dict_by_consecutive_keys(individual_mutations)

            # Save individual genome mutation file in virus-specific result directory
            safe_genome_id = genome_id.replace("|", "_").replace("/", "_")
            output_file_snp_roots = os.path.join(
                self.virus_processor.result_path, f"{safe_genome_id}_{today}.txt"
            )
            with open(output_file_snp_roots, "w") as file:
                for _key, value in sorted_mutations.items():
                    file.write(f"{value}\n")

            # Save row and hot mutations CSV
            if row_hot_mutations:
                csv_file = os.path.join(
                    self.virus_processor.result_path, f"{safe_genome_id}_row_hot_mutations_{today}.csv"
                )
                with open(csv_file, "w") as file:
                    # Write header
                    file.write("genome_id,mutation_type,position,nucleotide_change,protein_affected,amino_acid_change,biological_classification\n")
                    # Write data
                    for mutation in row_hot_mutations:
                        file.write(f"{mutation['genome_id']},{mutation['mutation_type']},{mutation['position']},{mutation['nucleotide_change']},{mutation['protein_affected']},{mutation['amino_acid_change']},{mutation['biological_classification']}\n")

            print(
                f"Generated mutation file for {genome_id}: {len(sorted_mutations)} mutations"
            )
            if row_hot_mutations:
                print(f"Generated row/hot mutations CSV for {genome_id}: {len(row_hot_mutations)} mutations")

    def extract_unique_snp_records(self, file_path: str) -> None:
        """
        Extract unique SNP records from a text file and save them in a text file.

        Args:
            file_path (str): Path to the text file containing SNP records.
        """
        # Read the text file and parse the SNP records
        snp_records = []

        with open(file_path) as txt_file:
            for line in txt_file:
                line = line.strip()
                if ":" in line and "[" in line and "]" in line:
                    # Parse format: "genome_id: [key1, key2, key3, ...]"
                    parts = line.split(":", 1)
                    if len(parts) == 2:
                        seq_id = parts[0].strip()
                        keys_str = parts[1].strip()

                        # Extract keys from the list format [key1, key2, key3, ...]
                        if keys_str.startswith("[") and keys_str.endswith("]"):
                            keys_content = keys_str[1:-1]  # Remove brackets
                            if keys_content:
                                # Split by comma and convert to integers
                                keys = [
                                    int(k.strip())
                                    for k in keys_content.split(",")
                                    if k.strip().isdigit()
                                ]

                                # Create a record with the expected structure
                                record = {
                                    "seqId": seq_id,
                                    "genomeSNPPs": [{"key": key} for key in keys],
                                }
                                snp_records.append(record)

        unique_keys = set()
        unique_snp_records = []

        for record in snp_records:
            # Extract keys from genomeSNPPs
            keys = [
                int(mut.get("key", 0))
                for mut in record.get("genomeSNPPs", [])
                if "key" in mut
            ]
            key_sum = sum(keys)
            if key_sum not in unique_keys:
                unique_keys.add(key_sum)
                unique_snp_records.append(record)

        print(f"Unique SNP Records: {len(unique_snp_records)}")

        output_file_unique_snps = os.path.join(
            self.virus_processor.data_path, f"Unique_SnpRecords_{today}.txt"
        )
        self.genome_snp_processor.print_snp_records(
            unique_snp_records, file_name=output_file_unique_snps
        )

    def get_snp_records_for_one_genome(
        self, msa_files: List[str], genome_id: str
    ) -> List[Dict]:
        """
        Retrieve SNP records for a single specific genome.

        Args:
            msa_files (List[str]): List of MSA file names.
            genome_id (str): The specific genome ID to process.

        Returns:
            List[Dict]: List of SNP records for the specified genome.
        """
        msa_file = msa_files[0]  # Assuming we're working with the first MSA file

        # Create a new gene mutation detector for the specific genome
        ref_genome_path = self.virus_processor.get_reference_genome_path()
        ref_genome = ReferenceGenome(ref_genome_path)
        proteome_path = self.virus_processor.get_proteome_path()
        proteome = Proteome(proteome_path)
        msa_instance = MultipleSequenceAlignment(os.path.join(self.virus_processor.data_path, msa_file))

        gene_mut_detector = GeneMutationDetector(ref_genome, msa_instance, genome_id)
        orf_processor = ORFProcessor(ref_genome)
        aa_mutation_processor = AminoAcidMutationProcessor()
        alignment_processor = AlignmentProcessor(proteome, aa_mutation_processor)
        pro_mutation_detector_instance = ProMutationDetector(
            orf_processor, alignment_processor, aa_mutation_processor, ref_genome
        )

        # Create a new genome SNP processor for the specific genome
        genome_snp_processor = GenomeSNPProcessor(
            self.virus_processor.data_path,
            msa_instance,
            gene_mut_detector,
            pro_mutation_detector_instance,
        )

        # Get SNP records for the specific genome
        snp_records = genome_snp_processor.get_one_genome_snps(genome_id)

        # Generate snpRoots file (this is useful)
        root_mutation_dict = genome_snp_processor.root_variant([snp_records])

        # Save in virus-specific result directory
        output_file_snp_roots = os.path.join(
            self.virus_processor.result_path, f"{genome_id}_{today}.txt"
        )
        with open(output_file_snp_roots, "w") as file:
            for _key, value in root_mutation_dict.items():
                file.write(f"{value}\n")

        # Also create row/hot mutations CSV for single genome
        row_hot_mutations = []
        for genome_snpp in snp_records.get("genomeSNPPs", []):
            original_type = genome_snpp.get('type', '')
            if original_type in ['rowMutation', 'hotMutation']:
                # Split multi-protein mutations
                separated_mutations = split_multi_protein_mutations(genome_snpp)
                for separated_mutation in separated_mutations:
                    # Extract protein and amino acid information
                    protein_info = separated_mutation.get('proteinMutation', [])
                    if protein_info and isinstance(protein_info[0], dict):
                        # For row and hot mutations, process all proteins
                        mutation_type = separated_mutation.get('type', '')
                        if mutation_type in ['rowMutation', 'hotMutation']:
                            # Process each protein separately for CSV
                            for protein_data in protein_info:
                                if isinstance(protein_data, dict):
                                    protein_name = protein_data.get('protein', 'unknown')
                                    amino_acid_change = protein_data.get('mutation', 'NA')
                                    biological_type = classify_mutation_type([protein_data])
                                    
                                    # Always check for hot mutation if SNP is X->Y and length > 1
                                    snp_info = separated_mutation.get('SNP', '')
                                    final_mutation_type = None
                                    if '->' in snp_info:
                                        original_seq, mutated_seq = snp_info.split('->')
                                        if len(original_seq) > 1 and len(mutated_seq) > 1:
                                            # For hot mutation detection, we need to extract the full sequence
                                            # from the reference genome when we have position ranges
                                            pos_str = separated_mutation.get('pos', '')
                                            if ':' in pos_str:
                                                # Extract position range (e.g., "10447:10449")
                                                start_pos, end_pos = map(int, pos_str.split(':'))
                                                # Get the full sequence from reference genome (1-based to 0-based)
                                                ref_genome = ReferenceGenome("data/refs/NC_045512.fasta")
                                                if len(original_seq) == 2 and len(mutated_seq) == 2:
                                                    # For 2-nucleotide changes like GC->AA, check if it's a hot mutation
                                                    # by looking at the 3-nucleotide context
                                                    if start_pos + 1 <= len(ref_genome.ref_seq):
                                                        # Get the 3-nucleotide sequence: start_pos, start_pos+1, end_pos
                                                        # For position 10447:10449, we want positions 10447, 10448, 10449
                                                        full_original_seq = ref_genome.ref_seq[start_pos-1:end_pos]
                                                        # Construct mutated sequence: first_nt->first_mut, middle_nt (conserved), last_nt->last_mut
                                                        full_mutated_seq = mutated_seq[0] + full_original_seq[1] + mutated_seq[1]
                                                        is_hot = detect_hot_mutation(full_original_seq, full_mutated_seq)
                                                        final_mutation_type = 'hot' if is_hot else 'row'
                                                    else:
                                                        final_mutation_type = 'row'
                                                else:
                                                    # For other cases, use the original detection logic
                                                    is_hot = detect_hot_mutation(original_seq, mutated_seq)
                                                    final_mutation_type = 'hot' if is_hot else 'row'
                                            else:
                                                # For single position mutations, use the original detection logic
                                                is_hot = detect_hot_mutation(original_seq, mutated_seq)
                                                final_mutation_type = 'hot' if is_hot else 'row'
                                        
                                        if not final_mutation_type:
                                            final_mutation_type = mutation_type if mutation_type in ['rowMutation', 'hotMutation'] else 'row'
                                        
                                        row_hot_mutations.append({
                                            'genome_id': genome_id,
                                            'mutation_type': final_mutation_type,
                                            'position': separated_mutation.get('pos', ''),
                                            'nucleotide_change': separated_mutation.get('SNP', ''),
                                            'protein_affected': protein_name,
                                            'amino_acid_change': format_amino_acid_change_for_csv(amino_acid_change),
                                            'biological_classification': biological_type
                                        })
                        else:
                            # For other mutation types, use the original logic
                            protein_name = protein_info[0].get('protein', 'unknown')
                            amino_acid_change = protein_info[0].get('mutation', 'NA')
                            biological_type = classify_mutation_type(protein_info)
                            # Always check for hot mutation if SNP is X->Y and length > 1
                            snp_info = separated_mutation.get('SNP', '')
                            mutation_type = None
                            if '->' in snp_info:
                                original_seq, mutated_seq = snp_info.split('->')
                                if len(original_seq) > 1 and len(mutated_seq) > 1:
                                    # For hot mutation detection, we need to extract the full sequence
                                    # from the reference genome when we have position ranges
                                    pos_str = separated_mutation.get('pos', '')
                                    if ':' in pos_str:
                                        # Extract position range (e.g., "10447:10449")
                                        start_pos, end_pos = map(int, pos_str.split(':'))
                                        # Get the full sequence from reference genome (1-based to 0-based)
                                        ref_genome = ReferenceGenome("data/refs/NC_045512.fasta")
                                        full_original_seq = ref_genome.ref_seq[start_pos-1:end_pos]
                                        # Construct the full mutated sequence
                                        if len(original_seq) == 2 and len(mutated_seq) == 2:
                                            # For 2-nucleotide changes like GC->AA, check if it's a hot mutation
                                            # by looking at the 3-nucleotide context
                                            if start_pos + 1 <= len(ref_genome.ref_seq):
                                                # Get the 3-nucleotide sequence: start_pos, start_pos+1, end_pos
                                                # For position 10447:10449, we want positions 10447, 10448, 10449
                                                full_original_seq = ref_genome.ref_seq[start_pos-1:end_pos]
                                                # Construct mutated sequence: first_nt->first_mut, middle_nt (conserved), last_nt->last_mut
                                                full_mutated_seq = mutated_seq[0] + full_original_seq[1] + mutated_seq[1]
                                                is_hot = detect_hot_mutation(full_original_seq, full_mutated_seq)
                                                mutation_type = 'hot' if is_hot else 'row'
                                            else:
                                                is_hot = detect_hot_mutation(original_seq, mutated_seq)
                                                mutation_type = 'hot' if is_hot else 'row'
                                        else:
                                            is_hot = detect_hot_mutation(original_seq, mutated_seq)
                                            mutation_type = 'hot' if is_hot else 'row'
                            if not mutation_type:
                                mutation_type = original_type if original_type in ['rowMutation', 'hotMutation'] else 'row'
                            row_hot_mutations.append({
                                'genome_id': genome_id,
                                'mutation_type': mutation_type,
                                'position': separated_mutation.get('pos', ''),
                                'nucleotide_change': separated_mutation.get('SNP', ''),
                                'protein_affected': protein_name,
                                'amino_acid_change': format_amino_acid_change_for_csv(amino_acid_change),
                                'biological_classification': biological_type
                            })

        # Save row and hot mutations CSV
        if row_hot_mutations:
            csv_file = os.path.join(
                self.virus_processor.result_path, f"{genome_id}_row_hot_mutations_{today}.csv"
            )
            with open(csv_file, "w") as file:
                # Write header
                file.write("genome_id,mutation_type,position,nucleotide_change,protein_affected,amino_acid_change,biological_classification\n")
                # Write data
                for mutation in row_hot_mutations:
                    file.write(f"{mutation['genome_id']},{mutation['mutation_type']},{mutation['position']},{mutation['nucleotide_change']},{mutation['protein_affected']},{mutation['amino_acid_change']},{mutation['biological_classification']}\n")
            
            print(f"Generated row/hot mutations CSV for {genome_id}: {len(row_hot_mutations)} mutations")

        return [snp_records]


def main():
    """
    Main function to process virus mutations.
    """
    parser = argparse.ArgumentParser(
        description="Process virus mutations with virus-specific organization"
    )
    parser.add_argument(
        "--virus",
        type=str,
        default="SARS-CoV-2",
        help="Virus name (e.g., SARS-CoV-2, HIV, Pox)",
    )
    parser.add_argument(
        "--genome-id",
        type=str,
        default="EPI_ISL_16327572",
        help="Specific genome ID to process (optional)",
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

    # Get MSA file
    msa_file = args.msa_file or virus_processor.get_default_msa_file()

    # Check if MSA file exists
    msa_file_path = os.path.join(virus_processor.data_path, msa_file)
    if not os.path.exists(msa_file_path):
        print(f"MSA file not found: {msa_file_path}")
        print(f"Available files in {virus_processor.data_path}:")
        if os.path.exists(virus_processor.data_path):
            for file in os.listdir(virus_processor.data_path):
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
        virus_processor.data_path,
        msa_instance,
        gene_mut_detector,
        pro_mutation_detector_instance,
    )

    # Initialize SNP processor
    snp_processor_instance = SNPProcessor(
        genome_snp_processor_instance, virus_processor
    )

    msa_files = [msa_file]

    if args.process_all:
        # Process all genomes in the MSA file
        print("Processing all genomes in the MSA file...")
        snp_processor_instance.process_snp_records(msa_files)
    else:
        # Process only one specific genome
        print(f"Processing genome: {args.genome_id}")
        snp_processor_instance.get_snp_records_for_one_genome(msa_files, args.genome_id)

    print(f"Processing complete. Results saved in: {virus_processor.result_path}")


if __name__ == "__main__":
    main()
