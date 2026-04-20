# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2026-04-15 12:22:57
LastModifiedBy: Rui Wang
LastEditTime: 2026-04-15 16:37:16
Email: wang.rui@nyu.edu
FilePath: /Test-EDS/PDBto2D.py
Description: 
'''
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Complete version of PDBto2D.py for processing IDs listed in list_file.txt.
Uses project-relative paths.
"""

import os
import os.path
from difflib import SequenceMatcher
import copy

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def _dot_bracket_to_ct(dbn_file, ct_file):
    """Convert a dot-bracket .out file to CT format in pure Python.

    Input file format (3 lines):
        >name
        SEQUENCE
        DOT-BRACKET   (supports () [] {} for pseudoknots)

    CT format:
        N  name
        i  nt  i-1  i+1  pair_or_0  i
    """
    with open(dbn_file) as f:
        lines = [l.strip() for l in f if l.strip()]

    name = lines[0].lstrip('>') if lines else 'RNA'
    seq  = ''
    db   = ''
    for line in lines[1:]:
        if line[0] in 'ACGUTNacgutn-':
            seq = line
        elif line[0] in '.([{<':
            db = line

    if not seq or not db:
        return False

    n = len(seq)
    pairs = [0] * (n + 1)  # 1-indexed; 0 = unpaired

    # Handle all bracket types: () [] {}
    bracket_pairs = [('(', ')'), ('[', ']'), ('{', '}')]
    for open_b, close_b in bracket_pairs:
        stack = []
        for i, c in enumerate(db, 1):
            if c == open_b:
                stack.append(i)
            elif c == close_b:
                if stack:
                    j = stack.pop()
                    pairs[i] = j
                    pairs[j] = i

    with open(ct_file, 'w') as f:
        f.write(f'{n}  {name}\n')
        for i in range(1, n + 1):
            nt   = seq[i - 1].upper().replace('T', 'U')
            prev = i - 1
            nxt  = i + 1 if i < n else 0
            f.write(f'{i}\t{nt}\t{prev}\t{nxt}\t{pairs[i]}\t{i}\n')

    return True

aaList = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']


def process_id(id, base_dir):
    """Process a single structure ID through the DSSR → 2D → CT pipeline.

    Args:
        id: Structure identifier (e.g. '7LYJ')
        base_dir: Root directory containing PDB_DSSR/, PDB_DSSR_2D/, PDB_DSSR_CT/ subdirs

    Returns:
        List of substructure filenames (without extension) that were processed, or [] on failure.
    """
    Problems = []
    print(f"\nProcessing {id}...")

    name = ''
    seq = ''
    db = ''
    chains = []
    chainType = {}
    subchainPos = {}
    chainBind = {}
    BasePairs = {}
    chainSeq = {}
    chainDB = {}
    metals = {}
    protein = 'NO'
    title = ''
    pseudo = 'NO'
    NT = []
  
    dssr_file = os.path.join(base_dir, 'PDB_DSSR', id+'.out')
    if not os.path.isfile(dssr_file):
        print(f'ERROR: DSSR file not found: {dssr_file}')
        return []
        
    with open(dssr_file, 'r') as f:
        lines = f.readlines()
  
    n = 0    
    while n < len(lines):
        line = lines[n].strip() if n < len(lines) else ''
        
        if 'no. of metals:' in line:
            me = line.split()[-1]
            if me != '0':
                try:
                    me = me.split('[')[1].split(']')[0].split(',')
                    for m in me:
                        metals[m.split('=')[0]] = m.split('=')[1]
                except:
                    pass
            n += 1
              
        elif line == 'Secondary structures in dot-bracket notation (dbn) as a whole and per chain':
            # Skip empty line
            n += 1
            # Get name (next non-empty line starting with >)
            while n < len(lines) and (not lines[n].strip() or not lines[n].strip().startswith('>')):
                n += 1
            if n < len(lines):
                name = lines[n].strip()
                n += 1
            
            # Get sequence (next non-empty line)
            while n < len(lines) and not lines[n].strip():
                n += 1
            if n < len(lines):
                seq = lines[n].strip()
                n += 1
            
            # Get dot-bracket (next non-empty line)
            while n < len(lines) and not lines[n].strip():
                n += 1
            if n < len(lines):
                db = lines[n].strip()
                n += 1
            
            if '[' in db or '{' in db:
                pseudo = 'YES'
            
            # Parse chain information
            while n < len(lines) and lines[n].strip().startswith('>'):
                line_parts = lines[n].strip()[1:].split()
                if len(line_parts) > 0:
                    chain_header = line_parts[0]
                    # Extract chain ID (e.g., "Chain A" -> "A", or "strand_A" -> "A")
                    if 'Chain' in chain_header:
                        c = chain_header.split('Chain')[1].strip()
                    elif 'strand_' in chain_header:
                        c = chain_header.split('strand_')[1]
                    else:
                        c = chain_header
                    
                    n += 1
                    # Get sequence for this chain
                    while n < len(lines) and not lines[n].strip():
                        n += 1
                    if n < len(lines):
                        chain_seq = lines[n].strip()
                        n += 1
                    
                    # Get dot-bracket for this chain
                    while n < len(lines) and not lines[n].strip():
                        n += 1
                    if n < len(lines):
                        chain_db = lines[n].strip()
                        n += 1
                    
                    chains.append(c)
                    chainType[c] = 'RNA'
                    chainSeq[c] = chain_seq
                    chainDB[c] = chain_db
                else:
                    n += 1
              
        elif n < len(lines) and lines[n][0:30] == 'Summary of structural features':
            n += 8
            while n < len(lines) and lines[n].strip() != '':
                NT.append(lines[n])
                n += 1
              
        else:
            n += 1

    # If we don't have chain info from parsing, create it from the whole structure
    if not chains and seq and db:
        # Split by & if multiple chains
        if '&' in seq:
            seq_parts = seq.split('&')
            db_parts = db.split('&')
            for i, (s, d) in enumerate(zip(seq_parts, db_parts)):
                c = chr(65 + i)  # A, B, C, etc.
                chains.append(c)
                chainSeq[c] = s
                chainDB[c] = d
                chainType[c] = 'RNA'
        else:
            # Single chain
            chains = ['A']
            chainSeq['A'] = seq
            chainDB['A'] = db
            chainType['A'] = 'RNA'

    # Write the original whole structure in ct format
    ct_out_file = os.path.join(base_dir, 'PDB_DSSR_CT', id+'.out')
    os.makedirs(os.path.dirname(ct_out_file), exist_ok=True)
    with open(ct_out_file, 'w') as f:
        f.write(f'>{id}\n')
        f.write(seq.replace('&','') + '\n')
        f.write(db.replace('&','') + '\n')
    
    # Convert to CT format (pure Python — no dot2ct binary needed)
    ct_file = os.path.join(base_dir, 'PDB_DSSR_CT', id+'.ct')
    if not _dot_bracket_to_ct(ct_out_file, ct_file):
        print(f"Warning: dot-bracket to CT conversion failed for {id}")

    # Find chain interactions from CT file if it exists
    if os.path.isfile(ct_file):
        relabel = {}
        respos = 1
        current_chain_idx = 0
        current_pos_in_chain = 0
        
        # Create mapping from position to chain
        for c in chains:
            chain_len = len(chainSeq.get(c, ''))
            for i in range(chain_len):
                relabel[respos] = c
                respos += 1
        
        # Read CT file to find base pairs
        with open(ct_file, 'r') as f:
            ctlines = f.readlines()
        
        for lp in range(1, len(ctlines)):
            parts = ctlines[lp].split()
            if len(parts) >= 5:
                try:
                    lpbp = int(parts[4])
                    if lpbp != 0 and lp in relabel and lpbp in relabel:
                        c1 = relabel[lp]
                        c2 = relabel[lpbp]
                        if c1 != c2:  # Inter-chain base pair
                            if c1 not in chainBind:
                                chainBind[c1] = []
                            if c2 not in chainBind:
                                chainBind[c2] = []
                            if c2 not in chainBind[c1]:
                                chainBind[c1].append(c2)
                            if c1 not in chainBind[c2]:
                                chainBind[c2].append(c1)
                except (ValueError, IndexError):
                    pass

    # Get title from PDB file if available
    pdb_file = os.path.join(base_dir, 'PDB_files', id + '.pdb')
    cif_file = os.path.join(base_dir, 'PDB_files', id + '.cif')
    
    if os.path.isfile(pdb_file):
        with open(pdb_file, 'r') as f:
            pdb = f.readlines()
        i = 0
        t = []
        firsttitleline = True
        while i < len(pdb) and protein != 'YES':
            l = pdb[i]
            i += 1
            if l[0:5] == 'TITLE':
                if firsttitleline:
                    t += l.split()[1:]
                    firsttitleline = False
                else:
                    t += l.split()[2:]
            if l[0:4] == 'ATOM':
                for aa in aaList:
                    if aa in l:
                        protein = 'YES'
        for w in t:
            title = title + w + ' '
    elif os.path.isfile(cif_file):
        with open(cif_file, 'r') as f:
            cif = f.readlines()
        for l in cif:
            if 'title' in l.lower() and not title:
                title = ' '.join(l.split()[1:])
            if l[0:4] == 'ATOM':
                for aa in aaList:
                    if aa in l:
                        protein = 'YES'
    
    if not title:
        title = f"{id} RNA structure"

    # Write 2D structure file
    os.makedirs(os.path.join(base_dir, 'PDB_DSSR_2D'), exist_ok=True)
    output_file = os.path.join(base_dir, 'PDB_DSSR_2D', id+'.txt')
    with open(output_file, 'w') as f:
        f.write(f'Title\t\t\t{title}\n')
        
        s = ' '.join(chains)
        f.write(f'Chain order\t\t{s}\n')
        
        f.write(f'Amino acids\t\t{protein}\n')
        
        f.write('Metal')
        if metals:
            for me in metals:
                f.write(f'\t\t\t{me}\t{metals[me]}\n')
        else:
            f.write('\n')
        
        f.write('Chains')
        for c in chains:
            f.write(f'\t\t\t{c}\t{chainType.get(c, "RNA")}\t')
            ss = ' '.join(chainBind.get(c, []))
            f.write(f'{ss}\n')
        
        f.write(f'Pseudoknot\t\t{pseudo}\n')
        
        f.write('\nWhole structure\n')
        f.write(f'>{id}\n')
        f.write(f'{seq}\n')
        f.write(f'{db}\n')
        
        f.write('\nInteracting substructures\n')
        
        # Find interacting substructures
        chainOrder = {c: i for i, c in enumerate(chains)}
        lib = copy.copy(chains)
        subs = []
        
        for c in chains:
            if c in lib:
                ic = list(set(chainBind.get(c, [])))
                if c in ic:
                    ic.remove(c)
                
                substructure = copy.copy(ic)
                new = copy.copy(ic)
                while new:
                    sub_new = copy.copy(substructure)
                    for x in new:
                        sub_new += chainBind.get(x, [])
                    sub_new = list(set(sub_new))
                    if c in sub_new:
                        sub_new.remove(c)
                    new = [s for s in sub_new if s not in substructure]
                    substructure = copy.copy(sub_new)
                
                substructure.append(c)
                substructure = sorted(substructure, key=lambda s: chainOrder.get(s, 999))
                for x in substructure:
                    if x in lib:
                        lib.remove(x)
                subs.append(substructure)
        
        # If no interactions found, treat each chain as separate substructure
        if not subs:
            for c in chains:
                subs.append([c])
        
        # Write substructures
        for substructure in subs:
            ic = ' '.join(substructure)
            subseq = ''
            subdb = ''
            for sc in substructure:
                subseq += chainSeq.get(sc, '')
                subdb += chainDB.get(sc, '')
            f.write(f'Chain {ic}\n')
            f.write(f'{subseq}\n')
            f.write(f'{subdb}\n\n')
    
    # Create CT files for each substructure
    with open(output_file, 'r') as f:
        lines = f.readlines()
    
    for i in range(len(lines)):
        if lines[i].strip() == 'Interacting substructures':
            if i+1 < len(lines):
                c = lines[i+1].split()[1:]
                seq_line = lines[i+2].strip() if i+2 < len(lines) else ''
                db_line = lines[i+3].strip() if i+3 < len(lines) else ''
                fname = id+'_' + ''.join(c)
                
                ct_out = os.path.join(base_dir, 'PDB_DSSR_CT', fname+'.out')
                ct_ct = os.path.join(base_dir, 'PDB_DSSR_CT', fname+'.ct')
                
                with open(ct_out, 'w') as f:
                    f.write(f'>{fname}\n')
                    f.write(f'{seq_line}\n')
                    f.write(f'{db_line}\n')
                
                _dot_bracket_to_ct(ct_out, ct_ct)
                if not os.path.isfile(ct_ct):
                    Problems.append(fname)
                i += 5

                while i < len(lines):
                    if lines[i].strip() != '':
                        c = lines[i].split()[1:]
                        seq_line = lines[i+1].strip() if i+1 < len(lines) else ''
                        db_line = lines[i+2].strip() if i+2 < len(lines) else ''
                        fname = id+'_' + ''.join(c)

                        ct_out = os.path.join(base_dir, 'PDB_DSSR_CT', fname+'.out')
                        ct_ct = os.path.join(base_dir, 'PDB_DSSR_CT', fname+'.ct')

                        with open(ct_out, 'w') as f:
                            f.write(f'>{fname}\n')
                            f.write(f'{seq_line}\n')
                            f.write(f'{db_line}\n')

                        _dot_bracket_to_ct(ct_out, ct_ct)
                        if not os.path.isfile(ct_ct):
                            Problems.append(fname)
                        i += 4
                    else:
                        i += 1
            break
    
    print(f"✓ Processed {id}")
    if Problems:
        print(f"  Problems: {Problems}")
    return Problems


if __name__ == "__main__":
    base_dir = _SCRIPT_DIR
    with open(os.path.join(base_dir, 'list_file.txt'), 'r') as f:
        _lines = f.readlines()
    ids = [x.strip() for x in _lines[0].split(',') if x.strip()]
    for _id in ids:
        process_id(_id, base_dir)
    print("\n" + "="*60)
    print("Processing complete!")
    print("="*60)
    print("\nNext step: Run Dual_Library.py to identify dual graph motifs for processed IDs")
