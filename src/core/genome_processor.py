"""
Author: Rui Wang
Date: 2023-10-31 12:42:13
LastModifiedBy: Rui Wang
LastEditTime: 2023-12-22 12:43:52
Email: rw3594@nyu.edu
FilePath: /MutParser/src/core/genome_processor.py
Description: 
This code is used to track mutation records from multiple genomes.
"""

"""Genome processing and SNP analysis functionality."""

from typing import Dict, List, Optional

from ..utils.mutation_utils import (
    convert_key_to_int,
    find_common_elements,
    is_next_genome,
    sort_dict_by_consecutive_keys,
    classify_mutation_type,
    split_multi_protein_mutations,
)
from ..utils.sequence_utils import has_only_valid_nts
from .mutation_detector import (
    GeneMutationDetector,
    MultipleSequenceAlignment,
)
from .protein_analyzer import (
    AlignmentProcessor,
    AminoAcidMutationProcessor,
    ORFProcessor,
    ProMutationDetector, 
)


class GenomeSNPProcessor:
    def __init__(
        self,
        path: str,
        msa: MultipleSequenceAlignment,
        gene_mut_detector: GeneMutationDetector,
        pro_mut_detector: ProMutationDetector,
    ):
        self.path = path
        self.msa = msa
        self.gene_mut_detector = gene_mut_detector
        self.pro_mut_detector = pro_mut_detector

    def get_one_genome_snps(self, genome_id: str) -> Dict:
        (
            row_mutations,
            hot_mutations,
            point_mutations,
        ) = self.gene_mut_detector.get_others_mutations()
        del_mutations = self.gene_mut_detector.get_deletions()
        ins_mutations = self.gene_mut_detector.get_insertions()

        # Add protein mutation information to each gene mutation
        # Process point mutations
        point_mutations_with_protein = self.pro_mut_detector.protein_point_row_snps(
            point_mutations
        )
        
        # Process row mutations
        row_mutations_with_protein = self.pro_mut_detector.protein_point_row_snps(
            row_mutations
        )
        
        # Process hot mutations
        hot_mutations_with_protein = self.pro_mut_detector.protein_hot_snps(
            hot_mutations
        )
        
        # Process deletion mutations
        del_mutations_with_protein = self.pro_mut_detector.protein_dels_snps(
            del_mutations
        )
        
        # Process insertion mutations
        ins_mutations_with_protein = self.pro_mut_detector.protein_ins_snps(
            ins_mutations
        )

        # New formatted output with protein mutations
        print(f"Analyzing genome: {genome_id}")
        print("=============ins_mutations=================")
        for mut in ins_mutations_with_protein:
            protein_info = mut.get("proteinMutation", [])
            print(f"{mut['type']} {mut['pos']} {mut['SNP']} {protein_info}")
        print("=============del_mutations=================")
        for mut in del_mutations_with_protein:
            protein_info = mut.get("proteinMutation", [])
            print(f"{mut['type']} {mut['pos']} {mut['SNP']} {protein_info}")
        print("=============point_mutations=================")
        for mut in point_mutations_with_protein:
            protein_info = mut.get("proteinMutation", [])
            print(f"{mut['type']} {mut['pos']} {mut['SNP']} {protein_info}")
        print("=============row_mutations=================")
        for mut in row_mutations_with_protein:
            protein_info = mut.get("proteinMutation", [])
            print(f"{mut['type']} {mut['pos']} {mut['SNP']} {protein_info}")
        print("=============hot_mutations=================")
        for mut in hot_mutations_with_protein:
            protein_info = mut.get("proteinMutation", [])
            print(f"{mut['type']} {mut['pos']} {mut['SNP']} {protein_info}")

        # Combine all mutations with protein information
        all_mutations = (
            point_mutations_with_protein
            + row_mutations_with_protein
            + hot_mutations_with_protein
            + ins_mutations_with_protein
            + del_mutations_with_protein
        )

        return {
            "seqId": genome_id,
            "genomeSNPPs": all_mutations,
        }

    def get_genome_snps(
        self, msa_name: str, genome_id: Optional[str] = None
    ) -> List[Dict]:
        # Use the reference genome name from the gene mutation detector
        ref_genome_name = self.gene_mut_detector.ref_seq.ref_name
        ref_seq_id, ref_seq_msa = self.msa.find_ref_msa(ref_genome_name)

        snp_records = []
        if genome_id:
            seq_id, seq_msa = self.msa.find_seq_msa(genome_id)
            if (
                seq_id is not None
                and has_only_valid_nts(seq_msa)
                and seq_id != ref_seq_id
            ):
                # Create new mutation detectors for this specific genome
                gene_mut_detector = GeneMutationDetector(
                    self.gene_mut_detector.ref_seq, self.msa, genome_id
                )
                
                # Create new components for ProMutationDetector
                orf_processor = ORFProcessor(self.gene_mut_detector.ref_seq)
                aa_mutation_processor = AminoAcidMutationProcessor()
                alignment_processor = AlignmentProcessor(
                    self.pro_mut_detector.alignment_processor.proteome,
                    aa_mutation_processor,
                )
                
                pro_mut_detector = ProMutationDetector(
                    orf_processor,
                    alignment_processor,
                    aa_mutation_processor,
                    self.gene_mut_detector.ref_seq,
                )
                
                # Create a temporary processor for this genome
                temp_processor = GenomeSNPProcessor(
                    self.path, self.msa, gene_mut_detector, pro_mut_detector
                )
                
                snp_record = temp_processor.get_one_genome_snps(genome_id)
                snp_records.append(snp_record)
        else:
            for idx, record in enumerate(self.msa.align):
                if has_only_valid_nts(str(record.seq)) and idx != ref_seq_id:
                    # Create new mutation detectors for this specific genome
                    gene_mut_detector = GeneMutationDetector(
                        self.gene_mut_detector.ref_seq, self.msa, record.id
                    )
                    
                    # Create new components for ProMutationDetector
                    orf_processor = ORFProcessor(self.gene_mut_detector.ref_seq)
                    aa_mutation_processor = AminoAcidMutationProcessor()
                    alignment_processor = AlignmentProcessor(
                        self.pro_mut_detector.alignment_processor.proteome,
                        aa_mutation_processor,
                    )
                    
                    pro_mut_detector = ProMutationDetector(
                        orf_processor,
                        alignment_processor,
                        aa_mutation_processor,
                        self.gene_mut_detector.ref_seq,
                    )
                    
                    # Create a temporary processor for this genome
                    temp_processor = GenomeSNPProcessor(
                        self.path, self.msa, gene_mut_detector, pro_mut_detector
                    )
                    
                    snp_record = temp_processor.get_one_genome_snps(record.id)
                    snp_records.append(snp_record)

        return snp_records

    def next_variants(self, snp_record0: Dict, snp_records: List[Dict]) -> List[str]:
        keys0 = snp_record0["keys"]
        return [
            rec["seqId"] for rec in snp_records if is_next_genome(keys0, rec["keys"])
        ]

    def root_variant(self, snp_records: list) -> dict:
        # Extract all mutation keys from genomeSNPPs for each record
        all_keys_lists = [
            [mut["key"] for mut in rec.get("genomeSNPPs", []) if "key" in mut]
            for rec in snp_records
        ]
        common_keys = find_common_elements(all_keys_lists)
        common_mutation_dict = {}
        for key in common_keys:
            for snp_record in snp_records:
                for genome_snpp in snp_record.get("genomeSNPPs", []):
                    if genome_snpp.get("key") == key:
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
                            common_mutation_dict[dict_key] = mutation
                        
                        break
        return sort_dict_by_consecutive_keys(common_mutation_dict)

    def root_variant_majority(
        self, snp_records: list, threshold_percentage: float = 0.8
    ) -> dict:
        """
        Find mutations present in a majority of genomes (default 80%).
        
        Args:
            snp_records: List of SNP records from all genomes
            threshold_percentage: Minimum percentage of genomes that must have the mutation (default 0.8 = 80%)
            
        Returns:
            Dictionary of mutations present in majority of genomes
        """
        if not snp_records:
            return {}
            
        # Extract all mutation keys from genomeSNPPs for each record
        all_keys_lists = [
            [mut["key"] for mut in rec.get("genomeSNPPs", []) if "key" in mut]
            for rec in snp_records
        ]
        
        # Count how many genomes have each mutation
        mutation_counts = {}
        for keys_list in all_keys_lists:
            for key in keys_list:
                mutation_counts[key] = mutation_counts.get(key, 0) + 1
        
        # Calculate threshold
        total_genomes = len(snp_records)
        threshold_count = int(total_genomes * threshold_percentage)
        
        # Find mutations that meet the threshold
        majority_keys = [
            key for key, count in mutation_counts.items() if count >= threshold_count
        ]
        
        # Build the mutation dictionary
        majority_mutation_dict = {}
        for key in majority_keys:
            for snp_record in snp_records:
                for genome_snpp in snp_record.get("genomeSNPPs", []):
                    if genome_snpp.get("key") == key:
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
                            majority_mutation_dict[dict_key] = mutation
                        
                        break
        
        return sort_dict_by_consecutive_keys(majority_mutation_dict)

    def print_snp_records(self, snp_records: list, file_name: str = None):
        # Print or write all mutation keys for each genome
        all_keys = []
        for rec in snp_records:
            genome_keys = [
                mut["key"] for mut in rec.get("genomeSNPPs", []) if "key" in mut
            ]
            all_keys.append({"seqId": rec.get("seqId"), "keys": genome_keys})
        if file_name:
            with open(file_name, "w") as f:
                for entry in all_keys:
                    f.write(f"{entry['seqId']}: {entry['keys']}\n")
        else:
            for entry in all_keys:
                print(f"{entry['seqId']}: {entry['keys']}")
        return all_keys
