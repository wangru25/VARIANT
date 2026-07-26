#!/usr/bin/env bash
# Run the HIV-1 -1 PRF scan INSIDE the NUPACK container and confirm NUPACK ran.
# Usage (from the repo root on the host):
#   docker run --rm -v "$HOME/Dropbox/github/VARIANT:/v" -w /v \
#       rw3594/dual_rag_if:1.0 bash tools/hiv1_nupack_test.sh
set -e
cd /v
echo "########## environment ##########"
which mfe && python3 -c "import sys; sys.path.insert(0,'src/core/prf'); import nupack_runner as nk; print('nupack_available:', nk.nupack_available())"

echo; echo "########## HIV-1 scan (NUPACK auto-enables) ##########"
python3 src/core/prf/prf_scanner.py \
  --fasta data/HIV-1/refs/NC_001802.1.fasta \
  --out /tmp/HIV-1_nupack --use-rnafold --use-pknots --pknots-top-n 5 \
  --trna data/tRNA_abundance_sample.csv --organism hiv_1

echo; echo "########## did NUPACK populate? ##########"
python3 - << 'PY'
import csv
rows=list(csv.DictReader(open('/tmp/HIV-1_nupack.prf_candidates.csv')))
folded=[r for r in rows if r['pknots_energy'] not in ('','nan')]
print(f"total sites: {len(rows)} | folded (top-N): {len(folded)}")
print(f"{'pos':>6} {'motif':8} {'pknots_E':>9} {'pk?':>5} {'nupack_E':>9} {'nu_pk?':>6} {'call':12} {'conf':8}")
for r in sorted(folded, key=lambda x: float(x['pknots_energy']))[:6]:
    print(f"{r['site_start_1based']:>6} {r['slippery_motif']:8} {r['pknots_energy']:>9} {str(r['pknots_has_pseudoknot']):>5} {str(r['nupack_energy']):>9} {str(r['nupack_has_pseudoknot']):>6} {r['structure_call']:12} {r['structure_confidence']:8}")
n_nu=sum(1 for r in folded if r['nupack_energy'] not in ('','nan'))
print(f"\nNUPACK populated on {n_nu}/{len(folded)} folded sites", "<-- WORKS" if n_nu>0 else "<-- still blank (check mfe/NUPACKHOME)")
PY
