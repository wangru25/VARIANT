# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: rw3594@nyu.edu
FilePath: /VARIANT/./plot/tmp/figure4_protein_impact.py
Description: Figure generation script for figure4_protein_impact.py.
'''


    print("Creating Figure 4: Protein-Level Mutation Impact...")
    
    # Get available viruses
    result_dir = "../result"
    if not os.path.exists(result_dir):
        print(f"❌ Result directory not found: {result_dir}")
        return
    
    virus_dirs = [d for d in os.listdir(result_dir) 
                  if os.path.isdir(os.path.join(result_dir, d)) and not d.startswith('.')]
    
    if not virus_dirs:
        print("❌ No virus data found in result directory")
        return
    
    print(f"Available viruses: {virus_dirs}")
    
    # For now, use the first available virus (you can modify this to analyze all or specific ones)
    virus_name = virus_dirs[0]
    print(f"Analyzing {virus_name}...")
    
    # Get protein mutation data
    protein_data, total_files = get_protein_mutation_data(virus_name)
    
    if not protein_data:
        print(f"❌ No protein mutation data found for {virus_name}")
        # Try SARS-CoV-2 if available
        if 'SARS-CoV-2' in virus_dirs and virus_name != 'SARS-CoV-2':
            virus_name = 'SARS-CoV-2'
            print(f"Trying {virus_name}...")
            protein_data, total_files = get_protein_mutation_data(virus_name)
            if not protein_data:
                print(f"❌ No protein mutation data found for {virus_name}")
                return
        else:
            return
    
    print(f"Analyzed {total_files} {virus_name} samples")
    print(f"Found mutations in {len(protein_data)} proteins")
    
    # Create heatmap
    fig1 = create_protein_impact_heatmap(protein_data, total_files, virus_name)
    fig1.write_html(f"figure4a_{virus_name}_protein_impact_heatmap.html")
    fig1.write_image(f"figure4a_{virus_name}_protein_impact_heatmap.pdf", width=800, height=800)
    
    # Create summary chart
    fig2 = create_protein_mutation_summary(protein_data, total_files, virus_name)
    fig2.write_html(f"figure4b_{virus_name}_protein_mutation_summary.html")
    fig2.write_image(f"figure4b_{virus_name}_protein_mutation_summary.pdf", width=800, height=800)
    
    # Create mutation type distribution
    fig3 = create_mutation_type_distribution(protein_data, total_files, virus_name)
    fig3.write_html(f"figure4c_{virus_name}_mutation_type_distribution.html")
    fig3.write_image(f"figure4c_{virus_name}_mutation_type_distribution.pdf", width=800, height=800)
    
    # Create functional impact radar
    fig4 = create_functional_impact_radar(protein_data, total_files, virus_name)
    fig4.write_html(f"figure4d_{virus_name}_functional_impact_radar.html")
    fig4.write_image(f"figure4d_{virus_name}_functional_impact_radar.pdf", width=800, height=800)
    
    print(f"✅ Figure 4 saved as:")
    print(f"   - figure4a_{virus_name}_protein_impact_heatmap.html")
    print(f"   - figure4a_{virus_name}_protein_impact_heatmap.pdf")
    print(f"   - figure4b_{virus_name}_protein_mutation_summary.html")
    print(f"   - figure4b_{virus_name}_protein_mutation_summary.pdf")
    print(f"   - figure4c_{virus_name}_mutation_type_distribution.html")
    print(f"   - figure4c_{virus_name}_mutation_type_distribution.pdf")
    print(f"   - figure4d_{virus_name}_functional_impact_radar.html")
    print(f"   - figure4d_{virus_name}_functional_impact_radar.pdf")
    
    # Print protein summary
    print(f"\n📊 Top 10 Most Mutated Proteins in {virus_name}:")
    protein_totals = []
    for protein, mutations in protein_data.items():
        total = sum(mutations.values())
        if total > 0:
            protein_totals.append((protein, total, mutations))
    
    protein_totals.sort(key=lambda x: x[1], reverse=True)
    for i, (protein, total, mutations) in enumerate(protein_totals[:10]):
        print(f"  {i+1}. {protein}: {total} mutations ({dict(mutations)})")
    
    return fig1, fig2, fig3, fig4

if __name__ == "__main__":
    main()