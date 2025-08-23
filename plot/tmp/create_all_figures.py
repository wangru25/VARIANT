# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: rw3594@nyu.edu
FilePath: /VARIANT/./plot/tmp/create_all_figures.py
Description: Script to create all figures for the VARIANT analysis.
'''


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
