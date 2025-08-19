#!/usr/bin/env python3
"""
Mutation Processor Module

This module contains the MutationProcessor class that handles the complete
mutation analysis pipeline for both single and multi-segment viruses.
"""

import os
import sys
import glob
from datetime import date
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .frameshift_detector import FrameshiftDetector
from .genome_processor import GenomeSNPProcessor
from .mutation_detector import GeneMutationDetector, MultipleSequenceAlignment
from .protein_analyzer import (
    AlignmentProcessor, AminoAcidMutationProcessor, ORFProcessor, 
    ProMutationDetector, Proteome
)
from .reference_genome import ReferenceGenome
from .snp_processor import SNPProcessor
from .virus_processor import VirusMutationProcessor
from ..utils.mutation_summary import extract_mutation_summary_to_csv
from ..utils.mutation_utils import split_multi_protein_mutations, classify_mutation_type, convert_key_to_int, sort_dict_by_consecutive_keys, detect_hot_mutation, format_amino_acid_change_for_csv

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
            snp_records = self._run_snp_processing(processors['snp_processor'], segment, args, msa_files)
            
            # Step 1.5: Process root variants (generate individual genome files and row hot mutations)
            if snp_records:
                self._process_root_variants(snp_records, "Combined", self.virus_processor, segment)
        
        # Step 2: Frameshift Analysis (if requested)
        if processors['frameshift_detector']:
            self._run_frameshift_analysis(processors, segment, args)
        
        # Step 3: Generate Mutation Summaries
        if not args.detect_frameshifts:
            self._generate_mutation_summaries(segment, args)
    
    def _run_snp_processing(self, snp_processor: SNPProcessor, segment: Optional[str], args, msa_files: List[str]) -> List[Dict]:
        """Run SNP processing."""
        if args.process_all or not args.genome_id:
            # Process all genomes
            if segment:
                print(f"Processing all genomes in segment {segment}...")
                return snp_processor.process_snp_records(msa_files, segment)
            else:
                print("Processing all genomes in the MSA file...")
                return snp_processor.process_snp_records(msa_files)
        else:
            # Process specific genome
            if segment:
                print(f"Processing genome: {args.genome_id} in segment {segment}")
                return snp_processor.get_snp_records_for_one_genome(msa_files, args.genome_id, segment)
            else:
                print(f"Processing genome: {args.genome_id}")
                return snp_processor.get_snp_records_for_one_genome(msa_files, args.genome_id)
    
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
        all_txt_files = glob.glob(txt_pattern)
        
        # Filter out files that don't contain genome mutation data
        txt_files = []
        for txt_file in all_txt_files:
            try:
                with open(txt_file, 'r') as f:
                    # Read first few lines to determine if this is a genome mutation file
                    lines = [f.readline().strip() for _ in range(5)]
                    
                    # Check if file contains mutation data structure
                    # Genome mutation files have lines with: mutation_type position nucleotide_change [protein_info]
                    has_mutation_structure = False
                    for line in lines:
                        if line and '[' in line and ']' in line and any(char.isdigit() for char in line):
                            # This looks like a mutation entry: "missense 123 A->T [{'protein': 'name', 'mutation': 'desc'}]"
                            has_mutation_structure = True
                            break
                    
                    if has_mutation_structure:
                        txt_files.append(txt_file)
                    else:
                        # Skip files that don't contain mutation data structure
                        continue
            except Exception:
                # If we can't read the file, skip it
                continue
        
        if not txt_files:
            print(f"No .txt files found in {result_path}")
            return
        
        print(f"Found {len(txt_files)} .txt files to process")
        
        # Process each file
        processed_count = 0
        for txt_file in txt_files:
            genome_id = os.path.basename(txt_file).replace(f"_{TODAY}.txt", "")
            
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

    def _process_root_variants(self, snp_records: List[Dict], prefix: str, virus_processor=None, segment: Optional[str] = None) -> None:
        """
        Process and save individual mutation records for each genome.

        Args:
            snp_records: List of SNP records.
            prefix: Prefix for the output file name.
            virus_processor: Virus processor instance.
            segment: Segment name for multi-segment viruses.
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
                    
                    # Create unique key that includes position and protein information to avoid collisions
                    # when multiple proteins are affected at the same position
                    pos_key = convert_key_to_int(str(separated_mutation.get("pos", "")))
                    protein_info = separated_mutation.get('proteinMutation', [])
                    protein_name = "unknown"
                    if protein_info and isinstance(protein_info[0], dict):
                        protein_name = protein_info[0].get('protein', 'unknown')
                    
                    # Create a unique key combining position and protein hash
                    import hashlib
                    protein_hash = hashlib.md5(protein_name.encode()).hexdigest()[:8]
                    unique_key = int(f"{pos_key}{int(protein_hash, 16) % 10000}")
                    
                    individual_mutations[unique_key] = mutation

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
                                                    if virus_processor:
                                                        ref_genome = ReferenceGenome(virus_processor.get_reference_genome_path(segment))
                                                    else:
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

            # Sort mutations by position (extract position from key for proper sorting)
            def sort_by_position(item):
                key, value = item
                # Extract position from the mutation string
                # Format: "missense 26529 G->A [{'protein': 'membrane_glycoprotein', 'mutation': 'D3N'}]"
                try:
                    parts = value.split()
                    if len(parts) >= 2:
                        # Handle position ranges like "10447:10449"
                        pos_str = parts[1]
                        if ':' in pos_str:
                            # For ranges, use the start position
                            pos = int(pos_str.split(':')[0])
                        else:
                            pos = int(pos_str)
                        return pos
                except (ValueError, IndexError):
                    return 0
                return 0
            
            # Sort by position, then by protein name for overlapping positions
            sorted_items = sorted(individual_mutations.items(), key=sort_by_position)
            sorted_mutations = dict(sorted_items)

            # Save individual genome mutation file in virus-specific result directory
            safe_genome_id = genome_id.replace("|", "_").replace("/", "_")
            result_path = self.virus_processor.get_result_path(segment)
            output_file_snp_roots = os.path.join(
                result_path, f"{safe_genome_id}_{TODAY}.txt"
            )
            with open(output_file_snp_roots, "w") as file:
                for _key, value in sorted_mutations.items():
                    file.write(f"{value}\n")

            # Save row and hot mutations CSV
            if row_hot_mutations:
                csv_file = os.path.join(
                    result_path, f"{safe_genome_id}_row_hot_mutations_{TODAY}.csv"
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