# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: rw3594@nyu.edu
FilePath: /VARIANT/./plot/tmp/figure3_genome_landscape.py
Description: Figure generation script for figure3_genome_landscape.py.
'''


    print("Creating Figure 3: Genome-Wide Mutation Landscape...")
    
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
    
    # Get mutation data
    mutation_data, total_files = get_mutation_data(virus_name)
    
    if not mutation_data:
        print(f"❌ No mutation data found for {virus_name}")
        return
    
    print(f"Analyzed {total_files} {virus_name} samples")
    
    # Estimate genome length from mutation data
    genome_length = max(mutation_data.keys()) if mutation_data else 30000
    
    # Create main landscape figure
    fig1 = create_genome_landscape_figure(mutation_data, total_files, virus_name, genome_length)
    fig1.write_html(f"figure3a_{virus_name}_genome_landscape.html")
    fig1.write_image(f"figure3a_{virus_name}_genome_landscape.pdf", width=1000, height=1000)
    
    # Create hotspot figure
    fig2 = create_mutation_hotspot_figure(mutation_data, total_files, virus_name)
    fig2.write_html(f"figure3b_{virus_name}_mutation_hotspots.html")
    fig2.write_image(f"figure3b_{virus_name}_mutation_hotspots.pdf", width=800, height=800)
    
    print(f"✅ Figure 3 saved as:")
    print(f"   - figure3a_{virus_name}_genome_landscape.html")
    print(f"   - figure3a_{virus_name}_genome_landscape.pdf")
    print(f"   - figure3b_{virus_name}_mutation_hotspots.html")
    print(f"   - figure3b_{virus_name}_mutation_hotspots.pdf")
    
    # Print hotspot summary
    print(f"\n🔥 Top 10 Mutation Hotspots in {virus_name}:")
    hotspot_threshold = 0.1
    hotspots = []
    for pos, counts in mutation_data.items():
        total_mutations = sum(counts.values())
        frequency = total_mutations / total_files
        if frequency >= hotspot_threshold:
            hotspots.append((pos, frequency, counts))
    
    hotspots.sort(key=lambda x: x[1], reverse=True)
    for i, (pos, freq, counts) in enumerate(hotspots[:10]):
        print(f"  {i+1}. Position {pos}: {freq*100:.1f}% ({dict(counts)})")
    
    return fig1, fig2

if __name__ == "__main__":
    main()
