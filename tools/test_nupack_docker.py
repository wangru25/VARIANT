#!/usr/bin/env python3
"""Test the docker-backed NUPACK path FROM THE HOST (not inside the container).
This is the VARIANT web-server scenario: NUPACK is not installed on the host;
the wrapper shells out to the docker image for the fold.

Run on the host (any env with python3; needs docker + the image pulled):
    python3 tools/test_nupack_docker.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "core", "prf"))
import nupack_runner as nk

print("docker client :", nk._docker_bin())
print("image present :", nk._docker_image_present(nk._docker_bin() or "docker", nk.DEFAULT_NUPACK_DOCKER_IMAGE)
      if nk._docker_bin() else False)
print("nupack_available():", nk.nupack_available())
print()

cases = [
    ("SARS-CoV-2 FSE",
     "UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU",
     True, -29.5),
    ("hairpin", "GGGGGAAAAACCCCC", False, -8.3),
]
ok_all = True
for name, seq, exp_pk, exp_e in cases:
    db, e, pk = nk.run_nupack(seq)   # <- goes through docker automatically
    ok = db is not None and pk == exp_pk and e is not None and abs(e - exp_e) < 1.0
    ok_all &= ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: energy={e} has_pseudoknot={pk} (expect pk={exp_pk}, e~{exp_e})")
    print(f"        {db}")
print("\nOVERALL:", "PASS -- docker-backed NUPACK works from the host" if ok_all else "FAIL")
sys.exit(0 if ok_all else 1)
