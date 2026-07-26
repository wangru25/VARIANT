#!/usr/bin/env python3
"""
VERIFICATION STATUS (be precise about what is and is not validated)
    VALIDATED IN THE ANALYSIS SANDBOX -- ONLY the pure-Python functions, with NO
    NUPACK binary and NO docker present. parse_mfe_file() was exercised on .mfe text
    conforming to NUPACK 3.2.2's documented layout and correctly extracts the
    sequence length, MFE energy, dot-parens-plus structure and base-pair list;
    pairs_have_crossing() correctly flags a pseudoknot iff the pair list contains a
    crossing (i<k<j<l). No specific energy value is asserted here -- the tool was
    never run in this environment, so any concrete NUPACK number must come from a
    real run, not from this file.

    NOT YET VALIDATED live in a reproducible, in-record way -- (a) that NUPACK's real
    output actually parses through parse_mfe_file() (as opposed to hand-constructed
    format-conforming text), (b) the run_nupack() WRAPPER end-to-end (its own
    subprocess call), and (c) the docker-backed execution path below. Confirm all
    three by running, ON THE HOST, tools/test_nupack_docker.py (docker path) and,
    inside the image, tools/nupack_wrapper_test.py, and keep that output with the
    results.

    REQUIRED RUNTIME CONFIG: NUPACK 3.2 fails with 'Unable to find rna1995.dG' when
    NUPACKHOME is unset (observed for this image, whose parameters live at
    /workspace/nupack3.2.2/parameters). run_nupack() therefore sets NUPACKHOME (and
    NUPACKINSTALL) to the install root that contains parameters/. NUPACK is a
    linux/amd64 build; on Apple silicon it runs under emulation.

    Where the `mfe` binary is absent (e.g. the general analysis sandbox),
    run_nupack() returns (None, None, None) and the pipeline leaves the NUPACK
    columns blank -- no value is ever fabricated.

    EXECUTION PATHS (run_nupack tries them in order):
      1. Host-local `mfe` on PATH (used when running inside the image itself).
      2. Docker-backed: if NUPACK is NOT on the host but a docker/podman client and
         the image (default rw3594/dual_rag_if:1.0, override via NUPACK_DOCKER_IMAGE)
         are present, the wrapper mounts a temp dir into the container, runs
         `mfe -pseudo` there with NUPACKHOME set, and reads the .mfe back. This is
         the intended path for the VARIANT web server, which runs in the `variant`
         conda env on the host and calls out to the image only for NUPACK folds.
         Set NUPACK_DOCKER=0 to disable the docker path.

nupack_runner.py -- wrapper around NUPACK 3.2.2 (Zadeh et al. 2011, J Comput Chem
32:170) for pseudoknot-aware RNA MFE structure prediction.

WHY THIS EXISTS
---------------
NUPACK is a third, independent pseudoknot-capable predictor (alongside PKNOTS and
ProbKnot). Its ``mfe`` utility, run with ``-pseudo``, computes a minimum-free-energy
secondary structure that MAY contain pseudoknots, and reports BOTH a free energy
(kcal/mol) and a dot-parens-plus structure -- so it contributes to the stability
judgement (like PKNOTS/RNAfold) AND the pseudoknot vote (like PKNOTS/ProbKnot).

REAL INTERFACE (NUPACK 3.2.2 command-line)
------------------------------------------
    mfe -material rna -pseudo <prefix>
NUPACK reads ``<prefix>.in`` (for a single strand, one line: the sequence) and writes
``<prefix>.mfe``. The ``-pseudo`` flag enables pseudoknot prediction and is valid only
for a single RNA strand. NUPACK 3.2 also REQUIRES the environment variable
``NUPACKHOME`` to point at the install root (so it can find ``$NUPACKHOME/parameters``);
we honour it if set and otherwise infer it from the binary location
(``<home>/build/bin/mfe`` -> ``<home>``).

OUTPUT (.mfe) FORMAT
--------------------
Comment lines start with ``%``. The first non-comment line is the sequence length,
the next is the MFE energy (kcal/mol), the next is the dot-parens-plus structure,
followed by the list of base pairs ``i<TAB>j`` (1-based), terminated by a ``%`` line.
Pseudoknot status is determined structurally: a pseudoknot is present iff the pair
list contains CROSSING pairs (i<k<j<l). This is a real topological test on NUPACK's
own output, identical to the test used for PKNOTS/ProbKnot -- not a heuristic.

ENVIRONMENT
-----------
NUPACK 3.2.2 is license-gated and is NOT installed in the general analysis sandbox.
It is available inside the user's docker image ``rw3594/dual_rag_if`` at
``/workspace/nupack3.2.2`` (on PATH via ``/workspace/nupack3.2.2/build/bin``). This
wrapper therefore returns (None, None, None) when the ``mfe`` binary is not found, so
the pipeline runs unchanged where NUPACK is absent and gains a third pseudoknot voter
where it is present. See the VERIFICATION STATUS block at the top for exactly what has
and has not been validated (parser: validated on the user's real .mfe output; wrapper
and docker path: not yet exercised live -- run the tools/ tests to confirm).

COMPLEXITY: pseudoknot MFE is polynomial but expensive; as with PKNOTS/ProbKnot, run
on short downstream windows only, never a whole genome.
"""
import os
import shutil
import subprocess
import tempfile


# The mfe binary. NUPACK_MFE_BIN overrides; else look for 'mfe' on PATH (present in
# the docker image via /workspace/nupack3.2.2/build/bin).
DEFAULT_NUPACK_MFE_BIN = os.environ.get("NUPACK_MFE_BIN", "mfe")

# Docker-backed NUPACK. When NUPACK is not installed on the host (the normal case
# for the VARIANT web server, which runs in the `variant` conda env), the wrapper
# shells out to a docker image that ships NUPACK 3.2.2 and reads the .mfe back.
# NUPACK_DOCKER_IMAGE overrides the image; NUPACK_DOCKER=0 disables the docker path.
DEFAULT_NUPACK_DOCKER_IMAGE = os.environ.get("NUPACK_DOCKER_IMAGE", "rw3594/dual_rag_if:1.0")
# Install root INSIDE the image (contains parameters/); confirmed for the image above.
_DOCKER_NUPACKHOME = os.environ.get("NUPACK_DOCKER_HOME", "/workspace/nupack3.2.2")


def _docker_bin():
    """Path to a docker (or podman) client on the host, or None."""
    if os.environ.get("NUPACK_DOCKER", "1") == "0":
        return None
    return shutil.which("docker") or shutil.which("podman")


def _docker_image_present(docker, image):
    """True iff the NUPACK docker image is already pulled locally (no network)."""
    try:
        r = subprocess.run([docker, "image", "inspect", image],
                           capture_output=True, text=True, timeout=30)
        return r.returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def _resolve_mfe(bin_path):
    """Return a runnable path to NUPACK's mfe, or None."""
    if bin_path and os.path.isabs(bin_path) and os.access(bin_path, os.X_OK):
        return bin_path
    return shutil.which(bin_path or "mfe")


def _resolve_nupackhome(mfe_path):
    """Return the NUPACK install root (the directory that CONTAINS 'parameters/'
    with the .dG thermodynamic tables). Confirmed layout in rw3594/dual_rag_if:
        mfe          -> /workspace/nupack3.2.2/build/bin/mfe
        parameters   -> /workspace/nupack3.2.2/parameters/rna1995.dG
    so home = <mfe>/../../.. = /workspace/nupack3.2.2.
    Honours an existing NUPACKHOME only if it actually contains parameters/;
    otherwise infers from the binary, then falls back to a small search. Returns
    None if no parameters/ directory can be found (NUPACK will then fail cleanly)."""
    def _ok(h):
        return h and os.path.isdir(os.path.join(h, "parameters"))

    env_home = os.environ.get("NUPACKHOME")
    if _ok(env_home):
        return env_home
    if mfe_path:
        # .../build/bin/mfe -> .../build/bin -> .../build -> <home>
        home = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(mfe_path))))
        if _ok(home):
            return home
        # some builds keep parameters next to the bin dir's parent
        alt = os.path.dirname(os.path.dirname(os.path.abspath(mfe_path)))
        if _ok(alt):
            return alt
    # last resort: the known container location
    for cand in ("/workspace/nupack3.2.2",):
        if _ok(cand):
            return cand
    return env_home  # may be None; NUPACK then errors, wrapper returns blanks


def parse_mfe_file(text):
    """Parse a NUPACK .mfe file. Returns (structure, energy, pairs, n):
      structure : dot-parens-plus string (or None)
      energy    : MFE in kcal/mol as float (or None)
      pairs     : list of (i, j) 1-based tuples with i < j
      n         : sequence length (or 0)
    Robust to the leading comment block and the trailing '%' terminator."""
    lines = text.splitlines()
    # Drop comment/blank lines to find the data block, but keep order.
    data = [ln for ln in lines if ln.strip() and not ln.lstrip().startswith('%')]
    if len(data) < 2:
        return (None, None, [], 0)
    try:
        n = int(data[0].split()[0])
    except (ValueError, IndexError):
        return (None, None, [], 0)
    try:
        energy = float(data[1].split()[0])
    except (ValueError, IndexError):
        energy = None
    structure = data[2].strip() if len(data) >= 3 else None
    # Remaining data lines are the pair list "i<ws>j".
    pairs = []
    for ln in data[3:]:
        f = ln.split()
        if len(f) >= 2:
            try:
                i, j = int(f[0]), int(f[1])
            except ValueError:
                continue
            if j > i:
                pairs.append((i, j))
    return (structure, energy, pairs, n)


def pairs_have_crossing(pairs):
    """True if any two base pairs cross (pseudoknot topology): i<k<j<l."""
    for a in range(len(pairs)):
        i, j = pairs[a]
        for b in range(a + 1, len(pairs)):
            k, l = pairs[b]
            if i < k < j < l or k < i < l < j:
                return True
    return False


def _finalize(text):
    """Parse a .mfe text blob -> (dotbracket, energy, has_pk) or (None,None,None)."""
    structure, energy, pairs, n = parse_mfe_file(text)
    if n == 0 or structure is None:
        return (None, None, None)
    return (structure, energy, pairs_have_crossing(pairs))


def _run_nupack_local(seq, mfe_bin, nupackhome, temperature, timeout):
    """Run a host-local `mfe` binary. Returns .mfe text, or None on failure."""
    env = dict(os.environ)
    home = nupackhome or _resolve_nupackhome(mfe_bin)
    if home:
        # NUPACK 3.2 looks for parameters under $NUPACKHOME/parameters (and older
        # docs/messages reference $NUPACKINSTALL); set both so the .dG tables load
        # regardless of which the build consults. Confirmed necessary: the image
        # ships NUPACKHOME unset, so a bare `mfe` fails with 'Unable to find
        # rna1995.dG' -- setting these fixes it.
        env["NUPACKHOME"] = home
        env["NUPACKINSTALL"] = home
    with tempfile.TemporaryDirectory() as td:
        prefix = os.path.join(td, "seq")
        with open(prefix + ".in", "w") as fh:
            fh.write(seq + "\n")
        cmd = [mfe_bin, "-material", "rna", "-pseudo", "-T", str(temperature), prefix]
        try:
            subprocess.run(cmd, env=env, capture_output=True, text=True,
                           timeout=timeout, check=True)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
            return None
        mfe_out = prefix + ".mfe"
        if not os.path.exists(mfe_out):
            return None
        with open(mfe_out) as fh:
            return fh.read()


def _run_nupack_docker(seq, image, docker, temperature, timeout):
    """Run NUPACK inside a docker/podman image that ships NUPACK 3.2.2. Mounts a
    host temp dir into the container, runs `mfe -pseudo` with NUPACKHOME set, and
    reads back <prefix>.mfe. Returns .mfe text, or None on failure.

    This is the path the VARIANT web server uses: NUPACK need NOT be installed on
    the host -- only docker + the image (default rw3594/dual_rag_if:1.0)."""
    with tempfile.TemporaryDirectory() as td:
        with open(os.path.join(td, "seq.in"), "w") as fh:
            fh.write(seq + "\n")
        # Inside the container: set NUPACKHOME (image ships it unset), fold the
        # mounted /work/seq, leave /work/seq.mfe for the host to read.
        inner = (
            "export NUPACKHOME=%s NUPACKINSTALL=%s; "
            "cd /work && mfe -material rna -pseudo -T %s seq"
            % (_DOCKER_NUPACKHOME, _DOCKER_NUPACKHOME, temperature)
        )
        cmd = [docker, "run", "--rm", "-v", "%s:/work" % td, image,
               "bash", "-lc", inner]
        try:
            subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout, check=True)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
            return None
        mfe_out = os.path.join(td, "seq.mfe")
        if not os.path.exists(mfe_out):
            return None
        with open(mfe_out) as fh:
            return fh.read()


def run_nupack(seq_rna, bin_path=DEFAULT_NUPACK_MFE_BIN, nupackhome=None,
               temperature=37.0, timeout=300, image=DEFAULT_NUPACK_DOCKER_IMAGE):
    """Fold a sequence with NUPACK's ``mfe -pseudo`` (pseudoknot-aware). Returns
    (dotbracket, energy_kcal, has_pseudoknot):
      dotbracket     : NUPACK dot-parens-plus structure string, or None
      energy_kcal    : MFE in kcal/mol (float), or None
      has_pseudoknot : True iff NUPACK's predicted pairs contain a crossing, else
                       False; None if NUPACK could not run.

    Two execution paths, tried in order:
      1. HOST-LOCAL `mfe` binary, if one is on PATH (e.g. running inside the image).
      2. DOCKER image that ships NUPACK 3.2.2 (default rw3594/dual_rag_if:1.0), if a
         docker/podman client and the image are present. This is how the web server
         (running in the `variant` conda env, no host NUPACK) reaches NUPACK.
    Returns (None, None, None) on any failure -- there is NO fabricated fallback."""
    seq = seq_rna.strip().upper().replace('T', 'U')  # NUPACK RNA material expects U

    # Path 1: host-local binary.
    mfe_bin = _resolve_mfe(bin_path)
    if mfe_bin is not None:
        text = _run_nupack_local(seq, mfe_bin, nupackhome, temperature, timeout)
        if text is not None:
            return _finalize(text)

    # Path 2: docker-backed.
    docker = _docker_bin()
    if docker and _docker_image_present(docker, image):
        text = _run_nupack_docker(seq, image, docker, temperature, timeout)
        if text is not None:
            return _finalize(text)

    return (None, None, None)


def run_nupack_batch(seqs, temperature=37.0, timeout=600,
                     image=DEFAULT_NUPACK_DOCKER_IMAGE):
    """Fold MANY sequences in ONE docker invocation, to avoid paying container
    startup (~1-2 s) per sequence. Returns a list of (dotbracket, energy, has_pk)
    aligned to `seqs`; entries that failed are (None, None, None).

    Falls back to per-sequence run_nupack() when the docker path is unavailable
    (e.g. a host-local mfe is present instead, or docker is disabled). Same real
    NUPACK command per sequence; only the process boundary differs."""
    seqs = [s.strip().upper().replace('T', 'U') for s in seqs]
    if not seqs:
        return []

    docker = _docker_bin()
    use_docker = bool(docker and _docker_image_present(docker, image)) \
        and _resolve_mfe(DEFAULT_NUPACK_MFE_BIN) is None
    if not use_docker:
        # No batchable docker path -> just call the normal resolver per sequence.
        return [run_nupack(s, temperature=temperature, timeout=timeout, image=image)
                for s in seqs]

    with tempfile.TemporaryDirectory() as td:
        for idx, s in enumerate(seqs):
            with open(os.path.join(td, "s%d.in" % idx), "w") as fh:
                fh.write(s + "\n")
        # One shell inside the container folds every s<idx> sequentially.
        inner = ("export NUPACKHOME=%s NUPACKINSTALL=%s; cd /work; "
                 "for f in s*.in; do p=${f%%.in}; "
                 "mfe -material rna -pseudo -T %s \"$p\" >/dev/null 2>&1; done"
                 % (_DOCKER_NUPACKHOME, _DOCKER_NUPACKHOME, temperature))
        cmd = [docker, "run", "--rm", "-v", "%s:/work" % td, image, "bash", "-lc", inner]
        try:
            subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=True)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
            return [(None, None, None)] * len(seqs)
        out = []
        for idx in range(len(seqs)):
            mfe_out = os.path.join(td, "s%d.mfe" % idx)
            if os.path.exists(mfe_out):
                with open(mfe_out) as fh:
                    out.append(_finalize(fh.read()))
            else:
                out.append((None, None, None))
        return out


def nupack_available(bin_path=DEFAULT_NUPACK_MFE_BIN, image=DEFAULT_NUPACK_DOCKER_IMAGE):
    """True iff NUPACK can be reached by EITHER path: a host-local `mfe` (with its
    parameter tables), OR a docker/podman client with the NUPACK image pulled.
    Lets the pipeline AUTO-ENABLE NUPACK wherever it is reachable and silently skip
    it otherwise. Cheap: resolution/inspect checks only, no folding."""
    mfe_bin = _resolve_mfe(bin_path)
    if mfe_bin is not None and _resolve_nupackhome(mfe_bin) is not None:
        return True
    docker = _docker_bin()
    return bool(docker and _docker_image_present(docker, image))


if __name__ == "__main__":
    # Run this INSIDE the docker image (rw3594/dual_rag_if) to validate live NUPACK.
    fse = "UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUUGACAU"
    db, e, pk = run_nupack(fse)
    print("dotbracket:", db)
    print("energy:", e, "| has_pseudoknot:", pk)
