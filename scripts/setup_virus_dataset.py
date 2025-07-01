#!/usr/bin/env python3
"""
Author: Rui Wang
Date: 2023-10-11 10:24:49
LastModifiedBy: Rui Wang
LastEditTime: 2023-12-22 12:43:52
Email: rw3594@nyu.edu
FilePath: /7_MutParser/scripts/setup_virus_dataset.py
Description: Script for setting up virus datasets and preparing data for analysis
"""

import argparse
import os
import shutil
import sys
from typing import Dict, List, Optional

import yaml

# Add the project root to the path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class VirusDatasetSetup:
    """Helper class for setting up virus datasets."""
    
    def __init__(self, config_file: str = "virus_config.yaml"):
        """
        Initialize the virus dataset setup.
        
        Args:
            config_file (str): Path to the virus configuration file
        """
        self.config_file = config_file
        self.config = self._load_config()
    
    def _load_config(self) -> Dict:
        """Load the virus configuration file."""
        try:
            with open(self.config_file) as f:
                return yaml.safe_load(f)
        except FileNotFoundError:
            print(f"Error: Configuration file '{self.config_file}' not found.")
            print(
                "Please create a virus_config.yaml file or use the --create-config option."
            )
            return {}
    
    def list_viruses(self) -> None:
        """List all configured viruses."""
        if not self.config or "viruses" not in self.config:
            print("No viruses configured. Use --add-virus to add a new virus.")
            return
        
        print("Configured viruses:")
        print("-" * 50)
        for virus_name, config in self.config["viruses"].items():
            print(f"Name: {virus_name}")
            print(f"Description: {config.get('description', 'No description')}")
            print(
                f"Reference Genome: {config.get('reference_genome', 'Not specified')}"
            )
            print(f"Proteome File: {config.get('proteome_file', 'Not specified')}")
            print(f"Codon Table ID: {config.get('codon_table_id', 1)}")
            print(
                f"Default MSA File: {config.get('default_msa_file', 'Not specified')}"
            )
            print("-" * 50)
    
    def add_virus(
        self,
        virus_name: str,
        reference_genome: str,
        proteome_file: Optional[str] = None,
        codon_table_id: int = 1,
        description: str = "",
        default_msa_file: str = "",
    ) -> None:
        """
        Add a new virus to the configuration.
        
        Args:
            virus_name (str): Name of the virus
            reference_genome (str): Name of the reference genome file
            proteome_file (str, optional): Name of the proteome file
            codon_table_id (int): Codon table ID (default: 1 for standard genetic code)
            description (str): Description of the virus
            default_msa_file (str): Default MSA file name
        """
        if not self.config:
            self.config = {"viruses": {}, "global": {}}
        
        if "viruses" not in self.config:
            self.config["viruses"] = {}
        
        # Add the new virus configuration
        self.config["viruses"][virus_name] = {
            "reference_genome": reference_genome,
            "proteome_file": proteome_file or f"{virus_name}_proteome.fasta",
            "codon_table_id": codon_table_id,
            "description": description or f"{virus_name} virus",
            "default_msa_file": default_msa_file or f"{virus_name}_msa.txt",
        }
        
        # Save the updated configuration
        self._save_config()
        
        # Create directory structure
        self._create_virus_directories(virus_name)
        
        print(f"✓ Added virus '{virus_name}' to configuration")
        print(f"✓ Created directory structure: data/{virus_name}/")
        print(f"  - data/{virus_name}/clustalW/ (for MSA files)")
        print(f"  - data/{virus_name}/refs/ (for reference files)")
        print(f"  - data/{virus_name}/fasta/ (for FASTA files)")
        print(f"  - result/{virus_name}/ (for output files)")
        print("\nNext steps:")
        print(
            f"1. Place your reference genome file in: data/{virus_name}/refs/{reference_genome}"
        )
        if proteome_file:
            print(
                f"2. Place your proteome file in: data/{virus_name}/refs/{proteome_file}"
            )
        print(f"3. Place your MSA files in: data/{virus_name}/clustalW/")
        print(f"4. Run: python main.py --virus {virus_name}")
    
    def setup_virus_dataset(
        self,
        virus_name: str,
                           reference_genome_path: Optional[str] = None,
                           proteome_path: Optional[str] = None,
        msa_files: Optional[List[str]] = None,
    ) -> None:
        """
        Set up a complete virus dataset with file copying.
        
        Args:
            virus_name (str): Name of the virus
            reference_genome_path (str, optional): Path to reference genome file
            proteome_path (str, optional): Path to proteome file
            msa_files (List[str], optional): List of paths to MSA files
        """
        # Check if virus is configured
        if (
            not self.config
            or "viruses" not in self.config
            or virus_name not in self.config["viruses"]
        ):
            print(f"Error: Virus '{virus_name}' is not configured.")
            print("Please add it first using --add-virus")
            return
        
        virus_config = self.config["viruses"][virus_name]
        
        # Create directories
        self._create_virus_directories(virus_name)
        
        # Copy reference genome
        if reference_genome_path:
            target_ref = f"data/{virus_name}/refs/{virus_config['reference_genome']}"
            self._copy_file(reference_genome_path, target_ref, "Reference genome")
        
        # Copy proteome file
        if proteome_path:
            target_prot = f"data/{virus_name}/refs/{virus_config['proteome_file']}"
            self._copy_file(proteome_path, target_prot, "Proteome file")
        
        # Copy MSA files
        if msa_files:
            for msa_file in msa_files:
                target_msa = f"data/{virus_name}/clustalW/{os.path.basename(msa_file)}"
                self._copy_file(
                    msa_file, target_msa, f"MSA file: {os.path.basename(msa_file)}"
                )
        
        print(f"✓ Virus dataset '{virus_name}' setup complete!")
        print(f"✓ Files are ready in data/{virus_name}/")
        print(f"✓ Run analysis with: python main.py --virus {virus_name}")
    
    def _create_virus_directories(self, virus_name: str) -> None:
        """Create the directory structure for a virus."""
        directories = [
            f"data/{virus_name}/clustalW",
            f"data/{virus_name}/refs",
            f"data/{virus_name}/fasta",
            f"result/{virus_name}",
        ]
        
        for directory in directories:
            os.makedirs(directory, exist_ok=True)
    
    def _copy_file(self, source: str, target: str, file_type: str) -> None:
        """Copy a file with error handling."""
        try:
            if os.path.exists(source):
                shutil.copy2(source, target)
                print(f"✓ Copied {file_type}: {source} → {target}")
            else:
                print(f"⚠ Warning: {file_type} not found: {source}")
        except Exception as e:
            print(f"✗ Error copying {file_type}: {e}")
    
    def _save_config(self) -> None:
        """Save the configuration to file."""
        try:
            with open(self.config_file, "w") as f:
                yaml.dump(self.config, f, default_flow_style=False, indent=2)
        except Exception as e:
            print(f"Error saving configuration: {e}")
    
    def validate_virus_dataset(self, virus_name: str) -> bool:
        """
        Validate that a virus dataset is properly set up.
        
        Args:
            virus_name (str): Name of the virus to validate
            
        Returns:
            bool: True if dataset is valid, False otherwise
        """
        if (
            not self.config
            or "viruses" not in self.config
            or virus_name not in self.config["viruses"]
        ):
            print(f"Error: Virus '{virus_name}' is not configured.")
            return False
        
        virus_config = self.config["viruses"][virus_name]
        is_valid = True
        
        print(f"Validating virus dataset: {virus_name}")
        print("-" * 40)
        
        # Check directories
        directories = [
            f"data/{virus_name}/clustalW",
            f"data/{virus_name}/refs",
            f"data/{virus_name}/fasta",
            f"result/{virus_name}",
        ]
        
        for directory in directories:
            if os.path.exists(directory):
                print(f"✓ Directory exists: {directory}")
            else:
                print(f"✗ Missing directory: {directory}")
                is_valid = False
        
        # Check reference genome
        ref_path = f"data/{virus_name}/refs/{virus_config['reference_genome']}"
        if os.path.exists(ref_path):
            print(f"✓ Reference genome: {ref_path}")
        else:
            print(f"✗ Missing reference genome: {ref_path}")
            is_valid = False
        
        # Check proteome file (optional)
        prot_path = f"data/{virus_name}/refs/{virus_config['proteome_file']}"
        if os.path.exists(prot_path):
            print(f"✓ Proteome file: {prot_path}")
        else:
            print(f"⚠ Proteome file not found (optional): {prot_path}")
        
        # Check MSA files
        msa_dir = f"data/{virus_name}/clustalW"
        msa_files = (
            [f for f in os.listdir(msa_dir) if f.endswith((".txt", ".fasta", ".fa"))]
            if os.path.exists(msa_dir)
            else []
        )
        
        if msa_files:
            print(f"✓ Found {len(msa_files)} MSA file(s): {', '.join(msa_files)}")
        else:
            print(f"⚠ No MSA files found in: {msa_dir}")
        
        print("-" * 40)
        if is_valid:
            print(f"✓ Virus dataset '{virus_name}' is ready for analysis!")
        else:
            print(f"✗ Virus dataset '{virus_name}' needs attention.")
        
        return is_valid
    
    def create_config_template(self) -> None:
        """Create a template configuration file."""
        template = {
            "viruses": {
                "ExampleVirus": {
                    "reference_genome": "example_reference.fasta",
                    "proteome_file": "example_proteome.fasta",
                    "codon_table_id": 1,
                    "description": "Example virus for demonstration",
                    "default_msa_file": "example_msa.txt",
                }
            },
            "global": {
                "default_codon_table": 1,
                "supported_codon_tables": {
                    1: "Standard genetic code",
                    2: "Vertebrate mitochondrial code",
                    3: "Yeast mitochondrial code",
                },
            },
        }
        
        with open(self.config_file, "w") as f:
            yaml.dump(template, f, default_flow_style=False, indent=2)
        
        print(f"✓ Created template configuration file: {self.config_file}")
        print("Edit this file to add your virus configurations.")


def main():
    """Main function for the setup script."""
    parser = argparse.ArgumentParser(
        description="Setup virus datasets for MutParser",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List all configured viruses
  python setup_virus_dataset.py --list
  
  # Add a new virus
  python setup_virus_dataset.py --add-virus "MyVirus" --reference "my_ref.fasta"
  
  # Set up a complete dataset with files
  python setup_virus_dataset.py --setup "MyVirus" --reference-path "path/to/ref.fasta" --msa-files "file1.txt,file2.txt"
  
  # Validate a virus dataset
  python setup_virus_dataset.py --validate "MyVirus"
  
  # Create a template configuration file
  python setup_virus_dataset.py --create-config
        """,
    )
    
    parser.add_argument(
        "--list", action="store_true", help="List all configured viruses"
    )
    parser.add_argument(
        "--add-virus", type=str, help="Add a new virus to configuration"
    )
    parser.add_argument(
        "--reference", type=str, help="Reference genome filename (for --add-virus)"
    )
    parser.add_argument(
        "--proteome", type=str, help="Proteome filename (for --add-virus)"
    )
    parser.add_argument(
        "--codon-table", type=int, default=1, help="Codon table ID (default: 1)"
    )
    parser.add_argument(
        "--description",
        type=str,
        default="",
        help="Virus description (for --add-virus)",
    )
    parser.add_argument(
        "--default-msa",
        type=str,
        default="",
        help="Default MSA filename (for --add-virus)",
    )
    
    parser.add_argument("--setup", type=str, help="Set up a complete virus dataset")
    parser.add_argument(
        "--reference-path", type=str, help="Path to reference genome file (for --setup)"
    )
    parser.add_argument(
        "--proteome-path", type=str, help="Path to proteome file (for --setup)"
    )
    parser.add_argument(
        "--msa-files",
        type=str,
        help="Comma-separated list of MSA file paths (for --setup)",
    )
    
    parser.add_argument("--validate", type=str, help="Validate a virus dataset")
    parser.add_argument(
        "--create-config",
        action="store_true",
        help="Create a template configuration file",
    )
    parser.add_argument(
        "--config",
        type=str,
        default="virus_config.yaml",
        help="Configuration file path (default: virus_config.yaml)",
    )
    
    args = parser.parse_args()
    
    setup = VirusDatasetSetup(args.config)
    
    if args.create_config:
        setup.create_config_template()
    elif args.list:
        setup.list_viruses()
    elif args.add_virus:
        if not args.reference:
            print("Error: --reference is required when adding a virus")
            return
        setup.add_virus(
            args.add_virus,
            args.reference,
            args.proteome,
            args.codon_table,
            args.description,
            args.default_msa,
        )
    elif args.setup:
        msa_files = args.msa_files.split(",") if args.msa_files else None
        setup.setup_virus_dataset(
            args.setup, args.reference_path, args.proteome_path, msa_files
        )
    elif args.validate:
        setup.validate_virus_dataset(args.validate)
    else:
        parser.print_help()


if __name__ == "__main__":
    main() 
