# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: wang.rui@nyu.edu
FilePath: /VARIANT/./plot/tmp/figure6_prf_sites.py
Description: Figure generation script for figure6_prf_sites.py.
'''


    print("Creating Figure 6: PRF Site Distribution...")
    
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
    
    # Get PRF data
    prf_sites = get_prf_data(virus_name)
    
    if not prf_sites:
        print(f"❌ No PRF data found for {virus_name}")
        return
    
    print(f"Found {len(prf_sites)} PRF sites")
    
    # Create genome map
    fig1 = create_prf_genome_map(prf_sites, virus_name)
    fig1.write_html(f"figure6a_{virus_name}_prf_genome_map.html")
    fig1.write_image(f"figure6a_{virus_name}_prf_genome_map.pdf", width=1000, height=1000)
    
    # Create sequence analysis
    fig2 = create_prf_sequence_analysis(prf_sites, virus_name)
    fig2.write_html(f"figure6b_{virus_name}_prf_sequence_analysis.html")
    fig2.write_image(f"figure6b_{virus_name}_prf_sequence_analysis.pdf", width=1000, height=800)
    
    # Create biological significance analysis
    fig3 = create_prf_biological_significance(prf_sites, virus_name)
    fig3.write_html(f"figure6c_{virus_name}_prf_biological_significance.html")
    fig3.write_image(f"figure6c_{virus_name}_prf_biological_significance.pdf", width=800, height=800)
    
    print(f"✅ Figure 6 saved as:")
    print(f"   - figure6a_{virus_name}_prf_genome_map.html")
    print(f"   - figure6a_{virus_name}_prf_genome_map.pdf")
    print(f"   - figure6b_{virus_name}_prf_sequence_analysis.html")
    print(f"   - figure6b_{virus_name}_prf_sequence_analysis.pdf")
    print(f"   - figure6c_{virus_name}_prf_biological_significance.html")
    print(f"   - figure6c_{virus_name}_prf_biological_significance.pdf")
    
    # Print PRF summary
    print(f"\n🔄 PRF Site Summary for {virus_name}:")
    
    # Count by type
    type_counts = defaultdict(int)
    for site in prf_sites:
        type_counts[site['type']] += 1
    
    for prf_type, count in type_counts.items():
        print(f"  {prf_type}: {count} sites")
    
    # Count by sequence
    sequence_counts = defaultdict(int)
    for site in prf_sites:
        sequence_counts[site['sequence']] += 1
    
    print(f"\n🧬 Top PRF Sequences in {virus_name}:")
    sorted_sequences = sorted(sequence_counts.items(), key=lambda x: x[1], reverse=True)
    for i, (seq, count) in enumerate(sorted_sequences[:5]):
        print(f"  {i+1}. {seq}: {count} occurrences")
    
    # Analyze by protein region
    genome_regions = get_genome_regions(virus_name)
    print(f"\n📊 PRF Sites by Protein Region in {virus_name}:")
    region_counts = defaultdict(int)
    for site in prf_sites:
        pos = site['position']
        for region_name, (start, end) in genome_regions.items():
            if start <= pos <= end:
                region_counts[region_name] += 1
                break
        else:
            region_counts['Intergenic'] += 1
    
    sorted_regions = sorted(region_counts.items(), key=lambda x: x[1], reverse=True)
    for region, count in sorted_regions:
        print(f"  {region}: {count} PRF sites")
    
    return fig1, fig2, fig3

if __name__ == "__main__":
    main()
