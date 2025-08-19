#!/usr/bin/env python3
"""
Master Script: Create All VARIANT Figures
========================================

This script runs all figure generation scripts to create comprehensive
visualizations for the VARIANT paper. Generic version that works with any virus data.
"""

import os
import sys
import subprocess
import glob
from pathlib import Path

def run_figure_script(script_name, description):
    """
    Run a figure generation script and handle any errors.
    
    Args:
        script_name (str): Name of the script to run
        description (str): Description of what the script does
    """
    print(f"\n{'='*60}")
    print(f"Running {script_name}: {description}")
    print(f"{'='*60}")
    
    script_path = Path(script_name)
    if not script_path.exists():
        print(f"❌ Script not found: {script_name}")
        return False
    
    try:
        # Run the script
        result = subprocess.run([sys.executable, script_name], 
                              capture_output=True, text=True, cwd=Path.cwd())
        
        if result.returncode == 0:
            print(f"✅ {script_name} completed successfully")
            if result.stdout:
                print("Output:")
                print(result.stdout)
            return True
        else:
            print(f"❌ {script_name} failed with return code {result.returncode}")
            if result.stderr:
                print("Error output:")
                print(result.stderr)
            return False
            
    except Exception as e:
        print(f"❌ Error running {script_name}: {e}")
        return False

def check_requirements():
    """
    Check if required packages are available.
    
    Returns:
        bool: True if all requirements are met
    """
    required_packages = ['plotly', 'pandas', 'numpy']
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing_packages.append(package)
    
    if missing_packages:
        print(f"❌ Missing required packages: {', '.join(missing_packages)}")
        print("Please install them using: pip install plotly pandas numpy")
        return False
    
    return True

def check_data_availability():
    """
    Check if result data is available.
    
    Returns:
        bool: True if data is available
    """
    result_dir = Path("../result")
    if not result_dir.exists():
        print(f"❌ Result directory not found: {result_dir}")
        return False
    
    virus_dirs = [d for d in result_dir.iterdir() 
                  if d.is_dir() and not d.name.startswith('.')]
    
    if not virus_dirs:
        print("❌ No virus data found in result directory")
        return False
    
    print(f"✅ Found virus data for: {[d.name for d in virus_dirs]}")
    return True

def list_generated_files():
    """
    List all generated figure files.
    """
    print(f"\n{'='*60}")
    print("Generated Figure Files")
    print(f"{'='*60}")
    
    # Find all generated files
    html_files = glob.glob("figure*.html")
    pdf_files = glob.glob("figure*.pdf")
    
    if html_files:
        print("\n📊 HTML Files (Interactive):")
        for file in sorted(html_files):
            print(f"  - {file}")
    
    if pdf_files:
        print("\n📄 PDF Files (Static):")
        for file in sorted(pdf_files):
            print(f"  - {file}")
    
    if not html_files and not pdf_files:
        print("❌ No figure files found")
    
    total_files = len(html_files) + len(pdf_files)
    print(f"\n📈 Total files generated: {total_files}")

def main():
    """Main function to run all figure generation scripts."""
    print("🎨 VARIANT Figure Generation Suite")
    print("=" * 50)
    
    # Check requirements
    if not check_requirements():
        return
    
    # Check data availability
    if not check_data_availability():
        return
    
    # Define all figure scripts
    figure_scripts = [
        ("figure1_workflow_overview.py", "Workflow Overview Diagram"),
        ("figure2_mutation_distribution.py", "Mutation Type Distribution Across Viruses"),
        ("figure3_genome_landscape.py", "Genome-Wide Mutation Landscape"),
        ("figure4_protein_impact.py", "Protein-Level Mutation Impact"),
        ("figure5_hot_mutations.py", "Hot Mutation Analysis"),
        ("figure6_prf_sites.py", "PRF Site Distribution")
    ]
    
    # Track success/failure
    successful_runs = []
    failed_runs = []
    
    # Run each figure script
    for script_name, description in figure_scripts:
        success = run_figure_script(script_name, description)
        if success:
            successful_runs.append(script_name)
        else:
            failed_runs.append(script_name)
    
    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"✅ Successful: {len(successful_runs)}/{len(figure_scripts)}")
    print(f"❌ Failed: {len(failed_runs)}/{len(figure_scripts)}")
    
    if successful_runs:
        print(f"\n✅ Successfully generated figures for:")
        for script in successful_runs:
            print(f"  - {script}")
    
    if failed_runs:
        print(f"\n❌ Failed to generate figures for:")
        for script in failed_runs:
            print(f"  - {script}")
    
    # List generated files
    list_generated_files()
    
    # Final message
    if successful_runs:
        print(f"\n🎉 Figure generation completed!")
        print("📁 Check the generated HTML and PDF files in the current directory.")
        print("🌐 Open HTML files in a web browser for interactive visualizations.")
        print("📄 PDF files are ready for publication.")
    else:
        print(f"\n💥 No figures were generated successfully.")
        print("Please check the error messages above and ensure data is available.")

if __name__ == "__main__":
    main()
