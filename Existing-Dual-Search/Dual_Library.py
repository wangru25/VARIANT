# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2026-04-15 12:22:57
LastModifiedBy: Rui Wang
LastEditTime: 2026-04-15 15:41:58
Email: wang.rui@nyu.edu
FilePath: /Test-EDS/Dual_Library.py
Description: 
'''
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modified version of Dual_Library.py for processing IDs listed in list_file.txt.
Uses project-relative directories.
"""

import os
import sys
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Add script directory to path so ClassesFunctions / dualGraphs are importable
sys.path.insert(0, _SCRIPT_DIR)

from igraph import *
from ClassesFunctions import *
from dualGraphs import *
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import re

def get_Graphs(fname, DualGraph_lib, NoCT, base_dir="."):
    ct_file = os.path.join(base_dir, 'PDB_DSSR_CT', fname+".ct")
    
    if os.path.isfile(ct_file):
        RNA = getCTInfo(ct_file)
        bpexist = countHelices(RNA)
        
        if bpexist:
            if len(RNA.Helices)>100:
                print('Too many helices!')
                output_file = os.path.join(base_dir, 'PDB_DSSR_Dual', fname+'.txt')
                with open(output_file, 'w') as f:
                    f.write('Too many helices!')
            
            else:                
                print('Extracting dual graphs for '+fname)
                changeHelices(RNA)
                RNA.makeMatrices()
                connectHelices(RNA)
                
                print ("Number of Vertices: " + str(len(RNA.Helices)-1))
                RNA.printAdj()
                RNA.printDeg()
                
                vertexOrder = []
                for i in range(0,len(RNA.adjMatrix)):
                    vertexOrder.append(0)
                
                success, graph=calcEigen(RNA, vertexOrder)
                correctHNumbers(RNA)
                    
                if success == 1:
                    if graph in DualGraph_lib:
                        DualGraph_lib[graph].append(fname)
                    else:
                        DualGraph_lib[graph] = [fname]
                    
                output_file = os.path.join(base_dir, 'PDB_DSSR_Dual', fname+'.txt')
                with open(output_file, 'w') as f:
                    if success == 1:
                        f.write(graph+'\n')
                    elif len(RNA.adjMatrix)==1:
                        f.write('Vertex number < 2\n')
                    elif len(RNA.adjMatrix)>9:
                        f.write('Vertex number > 9\n')
                    elif len(RNA.adjMatrix)==0:
                        f.write('No vertex\n')
                    else:
                        f.write("No matching graph exists (even if the vertex number is between 2 and 9).")
    
        else:
            print(fname+' has no base pair!')
            output_file = os.path.join(BASE_DIR, 'PDB_DSSR_Dual', fname+'.txt')
            with open(output_file, 'w') as f:
                f.write('No base pair')
    
    else:
        NoCT.append(fname)


def run_dual_library(ids, base_dir="."):
    """Run dual graph assignment for a list of structure IDs.

    Args:
        ids: List of structure identifiers (e.g. ['7LYJ'])
        base_dir: Root directory containing PDB_DSSR_2D/, PDB_DSSR_CT/, PDB_DSSR_Dual/

    Returns:
        dict with keys:
          'graph_assignments': {fname: graph_id_string}
          'no_ct': [fnames with missing CT files]
          'errors': [error messages]
    """
    os.makedirs(os.path.join(base_dir, 'PDB_DSSR_Dual'), exist_ok=True)

    DualGraph_lib = {}
    NoCT = []
    graph_assignments = {}

    for id in ids:
        if not id:
            continue

        pdb_2d_file = os.path.join(base_dir, 'PDB_DSSR_2D', id+'.txt')
        if not os.path.isfile(pdb_2d_file):
            print(f"Warning: {pdb_2d_file} not found. Run PDBto2D.py first.")
            continue

        with open(pdb_2d_file, 'r') as f:
            lines = f.readlines()

        for i in range(len(lines)):
            if lines[i] == 'Interacting substructures\n':
                if i+1 < len(lines):
                    c = lines[i+1].split()[1:]
                    fname = id+'_' + ''.join(c)
                    get_Graphs(fname, DualGraph_lib, NoCT, base_dir)
                    i += 5

                    while i < len(lines):
                        if lines[i] != '\n':
                            c = lines[i].split()[1:]
                            fname = id+'_' + ''.join(c)
                            get_Graphs(fname, DualGraph_lib, NoCT, base_dir)
                            i += 4
                        else:
                            i += 1

    # Collect results from written files
    import glob
    for id in ids:
        if id:
            files = glob.glob(os.path.join(base_dir, 'PDB_DSSR_Dual', id+'_*.txt'))
            for f in files:
                with open(f, 'r') as fh:
                    content = fh.read().strip()
                graph_assignments[os.path.basename(f).replace('.txt', '')] = content

    return {
        'graph_assignments': graph_assignments,
        'no_ct': NoCT,
    }


if __name__ == "__main__":
    import glob
    base_dir = _SCRIPT_DIR
    with open(os.path.join(base_dir, 'list_file.txt'), 'r') as f:
        _lines = f.readlines()
    ID = [id.strip() for id in _lines[0].split(',') if id.strip()]

    results = run_dual_library(ID, base_dir)

    print("\n" + "="*60)
    print("Processing complete!")
    print("="*60)
    print(f"\nDual graph assignments saved in: {os.path.join(base_dir, 'PDB_DSSR_Dual')}")
    print(f"\nFiles created:")
    for fname, content in results['graph_assignments'].items():
        print(f"  {fname}.txt: {content}")
    print("No CT:", results['no_ct'])
