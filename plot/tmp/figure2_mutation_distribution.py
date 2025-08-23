# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: rw3594@nyu.edu
FilePath: /VARIANT/./plot/tmp/figure2_mutation_distribution.py
Description: Figure generation script for figure2_mutation_distribution.py.
'''


    print("Creating Figure 2: Mutation Type Distribution Across Viruses...")
    
    # Get virus data
    virus_data = get_virus_data()
    
    if not virus_data:
        print("❌ No virus data found. Please ensure result files exist.")
        return
    
    print(f"Found data for {len(virus_data)} viruses: {list(virus_data.keys())}")
    
    # Create stacked bar chart
    fig1 = create_mutation_distribution_figure(virus_data)
    if fig1:
        fig1.write_html("figure2a_mutation_distribution_stacked.html")
        fig1.write_image("figure2a_mutation_distribution_stacked.pdf", width=800, height=800)
    
    # Create pie charts
    fig2 = create_mutation_type_pie_charts(virus_data)
    if fig2:
        fig2.write_html("figure2b_mutation_distribution_pie.html")
        fig2.write_image("figure2b_mutation_distribution_pie.pdf", width=800, height=800)
    
    print("✅ Figure 2 saved as:")
    print("   - figure2a_mutation_distribution_stacked.html")
    print("   - figure2a_mutation_distribution_stacked.pdf")
    print("   - figure2b_mutation_distribution_pie.html")
    print("   - figure2b_mutation_distribution_pie.pdf")
    
    # Print summary statistics
    print("\n📊 Summary Statistics:")
    for virus, data in virus_data.items():
        total = sum(data['counts'].values())
        print(f"  {virus}: {data['file_count']} samples, {total} total mutations")
        for mt, count in data['counts'].items():
            if count > 0:
                percentage = (count / total) * 100
                print(f"    {mt.capitalize()}: {count} ({percentage:.1f}%)")
    
    return fig1, fig2

if __name__ == "__main__":
    main()
