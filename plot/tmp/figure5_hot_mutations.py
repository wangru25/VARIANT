# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: rw3594@nyu.edu
FilePath: /VARIANT/./plot/tmp/figure5_hot_mutations.py
Description: Figure generation script for figure5_hot_mutations.py.
'''


    print("Creating Figure 5: Hot Mutation Analysis...")
    
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
    
    # Prioritize SARS-CoV-2 if available
    if 'SARS-CoV-2' in virus_dirs:
        virus_name = 'SARS-CoV-2'
        print(f"Prioritizing {virus_name} due to comprehensive data...")
    
    # Get hot mutation data
    hot_mutations, total_files = get_hot_mutation_data(virus_name)
    
    if not hot_mutations:
        print(f"❌ No hot mutation data found for {virus_name}")
        return
    
    print(f"Analyzed {total_files} {virus_name} samples")
    print(f"Found {len(hot_mutations)} hot mutations")
    
    # Create scatter plot
    fig1 = create_hot_mutation_scatter(hot_mutations, total_files, virus_name)
    fig1.write_html(f"figure5a_{virus_name}_hot_mutation_scatter.html")
    fig1.write_image(f"figure5a_{virus_name}_hot_mutation_scatter.pdf", width=800, height=800)
    
    # Create mutation type distribution
    fig2 = create_mutation_type_distribution(hot_mutations, virus_name)
    fig2.write_html(f"figure5b_{virus_name}_mutation_type_distribution.html")
    fig2.write_image(f"figure5b_{virus_name}_mutation_type_distribution.pdf", width=800, height=800)
    
    # Create protein impact analysis
    fig3 = create_protein_impact_analysis(hot_mutations, virus_name)
    fig3.write_html(f"figure5c_{virus_name}_protein_impact_analysis.html")
    fig3.write_image(f"figure5c_{virus_name}_protein_impact_analysis.pdf", width=800, height=800)
    
    # Create biological classification analysis
    fig4 = create_biological_classification_analysis(hot_mutations, virus_name)
    fig4.write_html(f"figure5d_{virus_name}_biological_classification.html")
    fig4.write_image(f"figure5d_{virus_name}_biological_classification.pdf", width=800, height=800)
    
    # Create sample distribution heatmap
    fig5 = create_sample_distribution_heatmap(hot_mutations, total_files, virus_name)
    fig5.write_html(f"figure5e_{virus_name}_sample_distribution_heatmap.html")
    fig5.write_image(f"figure5e_{virus_name}_sample_distribution_heatmap.pdf", width=1000, height=800)
    
    print(f"✅ Figure 5 saved as:")
    print(f"   - figure5a_{virus_name}_hot_mutation_scatter.html")
    print(f"   - figure5a_{virus_name}_hot_mutation_scatter.pdf")
    print(f"   - figure5b_{virus_name}_mutation_type_distribution.html")
    print(f"   - figure5b_{virus_name}_mutation_type_distribution.pdf")
    print(f"   - figure5c_{virus_name}_protein_impact_analysis.html")
    print(f"   - figure5c_{virus_name}_protein_impact_analysis.pdf")
    print(f"   - figure5d_{virus_name}_biological_classification.html")
    print(f"   - figure5d_{virus_name}_biological_classification.pdf")
    print(f"   - figure5e_{virus_name}_sample_distribution_heatmap.html")
    print(f"   - figure5e_{virus_name}_sample_distribution_heatmap.pdf")
    
    # Print summary statistics
    print(f"\n🔥 Hot Mutation Summary for {virus_name}:")
    
    # Count by type
    type_counts = defaultdict(int)
    for mutation in hot_mutations:
        type_counts[mutation['mutation_type']] += 1
    
    for mut_type, count in type_counts.items():
        print(f"  {mut_type.capitalize()}: {count} mutations")
    
    # Count by protein
    protein_counts = defaultdict(int)
    for mutation in hot_mutations:
        protein_counts[mutation['protein_affected']] += 1
    
    print(f"\n📊 Top 5 Proteins with Hot Mutations in {virus_name}:")
    sorted_proteins = sorted(protein_counts.items(), key=lambda x: x[1], reverse=True)
    for i, (protein, count) in enumerate(sorted_proteins[:5]):
        print(f"  {i+1}. {protein}: {count} hot mutations")
    
    return fig1, fig2, fig3, fig4, fig5

if __name__ == "__main__":
    main()
