#!/usr/bin/env python3
"""Validate the ACTUAL VARIANT run_nupack() wrapper inside the NUPACK container.
Run: docker run --rm -v "$HOME/Dropbox/github/VARIANT:/v" rw3594/dual_rag_if:1.0 \
        python3 /v/tools/nupack_wrapper_test.py
Prints PASS/FAIL for each case; exercises the real binary + real NUPACKHOME logic."""
import sys, os
sys.path.insert(0, "/v/src/core/prf")
import nupack_runner as nk

print("mfe on PATH:", nk._resolve_mfe("mfe"))
print("resolved NUPACKHOME:", nk._resolve_nupackhome(nk._resolve_mfe("mfe")))

cases = [
    ("SARS-CoV-2 FSE",
     "UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU",
     True,   -29.5),   # expect pseudoknot, ~-29.5
    ("hairpin", "GGGGGAAAAACCCCC", False, -8.3),
]
allpass = True
for name, seq, exp_pk, exp_e in cases:
    db, e, pk = nk.run_nupack(seq)
    ok = (db is not None and pk == exp_pk and e is not None and abs(e - exp_e) < 1.0)
    allpass &= ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: energy={e} has_pseudoknot={pk} (expected pk={exp_pk}, e~={exp_e})")
    print(f"        structure={db}")
print("\nOVERALL:", "PASS -- wrapper works end-to-end" if allpass else "FAIL")
sys.exit(0 if allpass else 1)
