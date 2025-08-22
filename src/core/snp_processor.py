# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-11 15:43:52
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-20 09:45:20
Email: wang.rui@nyu.edu
FilePath: /VARIANT/src/core/snp_processor.py
Description: 
'''

import os
import re
from typing import Dict, List, Optional

from src.core.genome_processor import GenomeSNPProcessor
from src.core.virus_processor import VirusMutationProcessor


class SNPProcessor:
    '''Class for processing SNP records and variants.'''

    def __init__(self, genome_snp_processor: GenomeSNPProcessor, virus_processor: VirusMutationProcessor):
        self.genome_snp_processor = genome_snp_processor
        self.virus_processor = virus_processor

    def get_snp_records_from_one_msa(self, msa_files: List[str], index: int, genome_id: Optional[str] = None) -> List[Dict]:
        '''
        Get SNP records from a single MSA file.
        
        Args:
            msa_files: List of MSA file paths
            index: Index of the MSA file to process
            genome_id: Optional genome ID to process only that specific genome
            
        Returns:
            List of SNP records
        '''
        msa_file = msa_files[index]
        print(f"Processing MSA file: {msa_file}")
        
        snp_records = self.genome_snp_processor.get_genome_snps(msa_file, genome_id)
        
        print(f"Found {len(snp_records)} SNP records")
        return snp_records

    def process_snp_records(self, msa_files: List[str], segment: Optional[str] = None) -> List[Dict]:
        '''
        Process SNP records from MSA files.
        
        Args:
            msa_files: List of MSA file paths
            segment: Segment name for multi-segment viruses
            
        Returns:
            List of SNP records from all MSA files
        '''
        print(f"Processing {len(msa_files)} MSA files...")
        
        all_snp_records = []
        
        for i, msa_file in enumerate(msa_files):
            print(f"\nProcessing MSA file {i+1}/{len(msa_files)}: {msa_file}")
            
            snp_records = self.get_snp_records_from_one_msa(msa_files, i, genome_id=None)
            
            if snp_records:
                all_snp_records.extend(snp_records)
                # Write SNP records to file
                prefix = f"msa_{i+1}"
                self._write_and_print_snp_records(snp_records, prefix, segment)
            else:
                print(f"No SNP records found in {msa_file}")
        
        print(f"Total SNP records from all files: {len(all_snp_records)}")
        return all_snp_records

    def _write_and_print_snp_records(self, snp_records: List[Dict], prefix: str, segment: Optional[str] = None) -> None:
        '''
        Print SNP records summary.
        
        Args:
            snp_records: List of SNP records
            prefix: File prefix
            segment: Segment name for multi-segment viruses
        '''
        print(f"Found {len(snp_records)} SNP records")


    def extract_unique_snp_records(self, file_path: str) -> None:
        '''
        Extract unique SNP records from a file.
        
        Args:
            file_path: Path to the file containing SNP records
        '''
        print(f"Extracting unique SNP records from: {file_path}")
        
        # Read all SNP records
        all_records = []
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    all_records.append(line)
        
        # Extract unique records
        unique_records = list(set(all_records))
        
        # Write unique records back to file
        with open(file_path, "w") as f:
            for record in unique_records:
                f.write(f"{record}\n")
        
        print(f"Original records: {len(all_records)}")
        print(f"Unique records: {len(unique_records)}")
        print(f"Removed {len(all_records) - len(unique_records)} duplicate records")

    def get_snp_records_for_one_genome(self, msa_files: List[str], genome_id: str, segment: Optional[str] = None) -> List[Dict]:
        '''
        Get SNP records for a specific genome.
        
        Args:
            msa_files: List of MSA file paths
            genome_id: Genome identifier
            segment: Segment name for multi-segment viruses
            
        Returns:
            List of SNP records for the specified genome
        '''
        print(f"Getting SNP records for genome: {genome_id}")
        
        all_snp_records = []
        
        for i, msa_file in enumerate(msa_files):
            print(f"Processing MSA file {i+1}/{len(msa_files)}: {msa_file}")
            
            snp_records = self.get_snp_records_from_one_msa(msa_files, i, genome_id=genome_id)
            
            all_snp_records.extend(snp_records)
            print(f"Found {len(snp_records)} records for {genome_id} in {msa_file}")
        
        print(f"Total records for {genome_id}: {len(all_snp_records)}")
        return all_snp_records