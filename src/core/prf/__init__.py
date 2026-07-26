# -*- coding: utf-8 -*-
"""
VARIANT programmed ribosomal frameshift (PRF) detection subpackage.

Modules
-------
frameshift_detector : FrameshiftDetector -- the detector wired into the `prf`
                      CLI subcommand. Detects -1 PRF (slippery-heptamer regex)
                      and +1 PRF (codon-in-frame: curated known sites, proline
                      P-site slippage, and leaky-stop candidates).
prf_scanner         : Standalone -1 PRF scanner (scan_seq / CLI). Regex slippery
                      search + RNAfold ranking + PKNOTS/ProbKnot/NUPACK pseudoknot
                      confirmation with a transparent structure consensus.
pknots_runner       : Wrapper around the bundled PKNOTS binary (MFE, pseudoknot-aware).
probknot_runner     : Wrapper around ProbKnot (RNAstructure; probability-based,
                      pseudoknot-aware). Optional -- blank if not installed.
nupack_runner       : Wrapper around NUPACK 3.2.2 (MFE, pseudoknot-aware), reached
                      via a Docker image. Optional -- blank if unavailable.

See README.md in this directory for the full algorithm and output-column reference.
"""
from .frameshift_detector import FrameshiftDetector

__all__ = ["FrameshiftDetector"]
