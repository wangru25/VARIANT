#!/usr/bin/env bash
# Run INSIDE rw3594/dual_rag_if:1.0. Locates NUPACK parameter files, sets
# NUPACKHOME, runs mfe -pseudo, and prints the real .mfe output.
set +e
echo "########## 1. locate the mfe binary and the parameter (.dG) files ##########"
MFE=$(which mfe); echo "mfe = $MFE"
echo "searching for rna1995.dG ..."
DG=$(find /workspace /opt / -name "rna1995.dG" 2>/dev/null | head -1)
echo "rna1995.dG = $DG"
# parameters dir = the directory that CONTAINS the .dG files
PARAMDIR=$(dirname "$DG")
echo "parameters dir = $PARAMDIR"
# NUPACKHOME = parent of the parameters dir
export NUPACKHOME=$(dirname "$PARAMDIR")
echo "NUPACKHOME (set) = $NUPACKHOME"

echo; echo "########## 2. real run: SARS-CoV-2 FSE (89 nt) with -pseudo ##########"
cd /tmp
printf "UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU\n" > fse.in
mfe -material rna -pseudo fse; echo ">>> exit=$?"
echo "=== fse.mfe CONTENTS ==="
cat fse.mfe

echo; echo "########## 3. control: simple hairpin (expect nested, no pseudoknot) ##########"
printf "GGGGGAAAAACCCCC\n" > hp.in
mfe -material rna -pseudo hp >/dev/null 2>&1; echo ">>> exit=$?"
echo "=== hp.mfe CONTENTS ==="
cat hp.mfe

echo; echo "########## 4. designed H-type pseudoknot (expect crossing pairs) ##########"
printf "GGGAAACCCCUUUGGGGAAAUUUCCC\n" > pk.in
mfe -material rna -pseudo pk >/dev/null 2>&1; echo ">>> exit=$?"
echo "=== pk.mfe CONTENTS ==="
cat pk.mfe
