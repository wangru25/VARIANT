# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:42:02
Email: wang.rui@nyu.edu
FilePath: /VARIANT/src/core/virus_processor.py
Description: Virus-specific processing module.
'''
import os
import yaml
from typing import Dict, List, Optional


class VirusMutationProcessor:
    """
    Class for processing virus-specific mutations with organized file structure.
    Supports both single-segment and multi-segment viruses.
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

        # Detect virus structure (single vs multi-segment)
        self.is_multi_segment = self._detect_multi_segment_structure()
        self.segments = self._get_segments()

        # Create virus-specific directories
        self.base_path = os.path.join("data", virus_name)
        self.result_path = os.path.join("result", virus_name)

        # Initialize paths based on virus structure
        self._initialize_paths()

        # Ensure directories exist
        self._ensure_directories()

    def _detect_multi_segment_structure(self) -> bool:
        """
        Detect if the virus has a multi-segment structure.
        
        Returns:
            bool: True if multi-segment, False if single-segment
        """
        base_path = os.path.join("data", self.virus_name)

        # Check if segment_1 directory exists
        segment_1_path = os.path.join(base_path, "segment_1")
        if os.path.exists(segment_1_path) and os.path.isdir(segment_1_path):
            return True

        # Check if traditional structure exists (clustalW, refs directories)
        clustalw_path = os.path.join(base_path, "clustalW")
        refs_path = os.path.join(base_path, "refs")
        
        if os.path.exists(clustalw_path) and os.path.exists(refs_path):
            return False

        # Default to single-segment if structure is unclear
        return False

    def _get_segments(self) -> List[str]:
        """
        Get list of segments for the virus.
        
        Returns:
            List[str]: List of segment names
        """
        if not self.is_multi_segment:
            return []

        base_path = os.path.join("data", self.virus_name)
        segments = []
        
        # Look for segment directories
        for item in os.listdir(base_path):
            item_path = os.path.join(base_path, item)
            if os.path.isdir(item_path) and item.startswith("segment_"):
                segments.append(item)

        return sorted(segments)

    def _initialize_paths(self) -> None:
        """Initialize virus-specific paths."""
        if self.is_multi_segment:
            # Multi-segment virus paths
            self.segment_paths = {}
            for segment in self.segments:
                self.segment_paths[segment] = {
                    "base": os.path.join(self.base_path, segment),
                    "result": os.path.join(self.result_path, segment),
                    "clustalw": os.path.join(self.base_path, segment, "clustalW"),
                    "refs": os.path.join(self.base_path, segment, "refs"),
                }
        else:
            # Single-segment virus paths
            self.clustalw_path = os.path.join(self.base_path, "clustalW")
            self.refs_path = os.path.join(self.base_path, "refs")

    def _ensure_directories(self) -> None:
        """Ensure all necessary directories exist."""
        if self.is_multi_segment:
            for segment_paths in self.segment_paths.values():
                os.makedirs(segment_paths["result"], exist_ok=True)
        else:
            os.makedirs(self.result_path, exist_ok=True)

    def _load_config(self) -> Dict:
        """
        Load virus configuration from YAML file.
        
        Returns:
            Dict: Configuration dictionary
        """
        try:
            with open(self.config_file, "r") as f:
                return yaml.safe_load(f)
        except FileNotFoundError:
            print(f"Warning: Config file {self.config_file} not found. Using default configuration.")
            return {}

    def get_virus_config(self, segment: Optional[str] = None) -> Dict[str, str]:
        """
        Get virus configuration for a specific segment.
        
        Args:
            segment (Optional[str]): Segment name for multi-segment viruses
            
        Returns:
            Dict[str, str]: Virus configuration
        """
        if self.is_multi_segment and segment:
            # Multi-segment virus
            if segment not in self.segments:
                raise ValueError(f"Segment {segment} not found for virus {self.virus_name}")
            
            # Get segment-specific config
            segment_config = self.config.get("viruses", {}).get(self.virus_name, {}).get("segments", {}).get(segment, {})
            
            # Merge with virus-level config (excluding segments)
            virus_config = self.config.get("viruses", {}).get(self.virus_name, {})
            virus_level_config = {k: v for k, v in virus_config.items() if k != "segments"}
            
            config = {**virus_level_config, **segment_config}
            
            # Add segment-specific paths
            config.update({
                "data_path": self.segment_paths[segment]["base"],
                "result_path": self.segment_paths[segment]["result"],
                "clustalw_path": self.segment_paths[segment]["clustalw"],
                "refs_path": self.segment_paths[segment]["refs"],
            })
            
            return config
        else:
            # Single-segment virus
            config = self.config.get("viruses", {}).get(self.virus_name, {})
            
            # Add default paths
            config.update({
                "data_path": self.base_path,
                "result_path": self.result_path,
                "clustalw_path": self.clustalw_path,
                "refs_path": self.refs_path,
            })
            
            return config

    def get_reference_genome_path(self, segment: Optional[str] = None) -> str:
        """
        Get reference genome file path.
        
        Args:
            segment (Optional[str]): Segment name for multi-segment viruses
            
        Returns:
            str: Path to reference genome file
        """
        config = self.get_virus_config(segment)
        refs_path = config.get("refs_path")
        if refs_path is None:
            if self.is_multi_segment and segment:
                refs_path = self.segment_paths[segment]["refs"]
            else:
                refs_path = self.refs_path
        ref_genome = config.get("reference_genome", "ref_genome.fasta")
        return os.path.join(refs_path, ref_genome)

    def get_proteome_path(self, segment: Optional[str] = None) -> str:
        """
        Get proteome file path.
        
        Args:
            segment (Optional[str]): Segment name for multi-segment viruses
            
        Returns:
            str: Path to proteome file
        """
        config = self.get_virus_config(segment)
        refs_path = config.get("refs_path")
        if refs_path is None:
            if self.is_multi_segment and segment:
                refs_path = self.segment_paths[segment]["refs"]
            else:
                refs_path = self.refs_path
        proteome = config.get("proteome_file", "proteome.fasta")
        return os.path.join(refs_path, proteome)

    def get_codon_table_id(self, segment: Optional[str] = None) -> int:
        """
        Get codon table ID for translation.
        
        Args:
            segment (Optional[str]): Segment name for multi-segment viruses
            
        Returns:
            int: Codon table ID
        """
        config = self.get_virus_config(segment)
        return config.get("codon_table_id", 1)

    def get_default_msa_file(self, segment: Optional[str] = None) -> str:
        """
        Get default MSA file path.
        
        Args:
            segment (Optional[str]): Segment name for multi-segment viruses
            
        Returns:
            str: Path to default MSA file
        """
        config = self.get_virus_config(segment)
        clustalw_path = config.get("clustalw_path")
        if clustalw_path is None:
            if self.is_multi_segment and segment:
                clustalw_path = self.segment_paths[segment]["clustalw"]
            else:
                clustalw_path = self.clustalw_path
        default_msa = config.get("default_msa_file", "aligned.fasta")
        return os.path.join(clustalw_path, default_msa)

    def get_data_path(self, segment: Optional[str] = None) -> str:
        """
        Get data directory path.
        
        Args:
            segment (Optional[str]): Segment name for multi-segment viruses
            
        Returns:
            str: Data directory path
        """
        config = self.get_virus_config(segment)
        return config.get("data_path", self.base_path)

    def get_result_path(self, segment: Optional[str] = None) -> str:
        """
        Get result directory path.
        
        Args:
            segment (Optional[str]): Segment name for multi-segment viruses
            
        Returns:
            str: Result directory path
        """
        config = self.get_virus_config(segment)
        return config.get("result_path", self.result_path)

    def list_available_viruses(self) -> List[str]:
        """
        List all available viruses in the configuration.
        
        Returns:
            List[str]: List of virus names
        """
        return list(self.config.get("viruses", {}).keys())

    def get_segments(self) -> List[str]:
        """
        Get list of segments for the virus.
        
        Returns:
            List[str]: List of segment names
        """
        return self.segments

    def is_multi_segment_virus(self) -> bool:
        """
        Check if the virus is multi-segment.
        
        Returns:
            bool: True if multi-segment, False otherwise
        """
        return self.is_multi_segment
