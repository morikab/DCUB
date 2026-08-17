# Welcome to DCUB – Differential Codon Usage Bias

DCUB is a software tool for designing genetic sequences that are tailored to specific microbial environments.  
With only a few steps, you can obtain a microbiome‑specific version of your gene of interest.

Please follow the steps below to install and run the tool.

---

## Prerequisites

A local Node.js / Python environment if building from source. Check the exact version constraints before you start:

- **Python** — Python version is pinned in [`pyproject.toml`](./pyproject.toml) (`python = "~3.11"`).  
- **Node.js / npm** — see [`ui/DCUB/package.json`](./ui/DCUB/package.json) and [`ui/DCUB/electron/package.json`](./ui/DCUB/electron/package.json) for dependency requirements; a recent Node LTS is recommended.

---

## Installation Guide (Pre‑built binaries)

Pre‑compiled versions of DCUB for macOS, Windows, and Linux are published as `.zip` archives on the [GitHub Releases page](https://github.com/morikab/DCUB/releases/latest):

- https://github.com/morikab/DCUB/releases/latest

They're also available from the project website, which includes per‑OS installation notes and troubleshooting tips:

- https://www.tau.ac.il/~bentulila

Download the `.zip` archive that matches your operating system, extract it, and run the app. See the installation notes linked above for OS‑specific steps (e.g. clearing the macOS quarantine flag, bypassing Windows SmartScreen, or making the Linux AppImage executable).

After installation, follow the usage instructions from the user guide packaged with the distribution.

---

## Building DCUB from Source

If you prefer to build the tool from source (for development or customization), use the following workflow:

0. **Install Python dependencies (Poetry)**

   - Install [Poetry](https://python-poetry.org/) (one-time), then from the project root run:
     ```
     poetry install --with build
     ```
   - This repo is configured so `poetry install` creates an in-project virtualenv at `.venv/`.
   - For notebooks / analysis tooling:
     ```
     poetry install --with build,analysis
     ```
   - To run the test suite:
     ```
     poetry install --with dev
     poetry run pytest app/tests -v
     ```

1. **Build the standalone FastAPI server**

   - From the project root, run:
     ```
     ./build_fastapi.sh
     ```
   - This script produces a standalone FastAPI server executable and copies it to a dedicated backend folder. 
   - You can test the generated executable by starting the sever locally:
        ```
        ./ui/DCUB/backend/fastapi_server
        ```


2. **Build the web UI**

   - From the `ui/DCUB` directory, run:
     ```
     npm run electron-dev
     ```
   - This builds the Next.js application and copies the resulting standalone front‑end build files into the Electron application structure to prepare a working Electron setup.

3. **Build the Electron desktop application**

   - Change directory to the Electron app folder (`cd ui/DCUB/electron`) and run:
     ```
     npm run build-electron
     ```
   - This produces the DCUB desktop application that bundles the UI and the backend server.

---

## Optimization options

### Hotspot avoidance (optional)

With **Advanced Options → Hotspot Avoidance** set to **On**, DCUB runs
[ESO](https://github.com/itamar-menuhin/evolutionary-stability-optimizer) over each
optimized candidate to detect hypermutable sites - replication slippage and
recombination-mediated deletion, by default - and edits them away.

Replacements inside a detected site are chosen using DCUB's own per-codon
preference model, so they still reflect the wanted/unwanted-organism tradeoff.
Every nucleotide outside a detected site is locked, so codon choices elsewhere
cannot drift.

This runs on every ORF-optimization candidate before evaluation, so the
reported scores describe the sequence that actually ships. Z-Score methods
produce several candidates (`1 + ZSCORE_INITIAL_PERMUTATIONS_NUM`, doubled for
`max_CAI_tAI`), so expect a corresponding slowdown with those methods.

While it is enabled, DCUB's own repeat-avoidance ("dedup") heuristic is turned
off for that run - ESO's slippage and recombination detection measures the same
thing directly.

It is **off by default**; runs without it are unaffected.

#### Methylation-motif detection is off by default

Motif detection (`dam`, `dcm`, and ESO's other bundled motifs) is **off by
default**, unlike slippage and recombination detection, which are always on
whenever hotspot avoidance is enabled. It is dominated by false positives:
ESO's PSSM-based motif scanner keeps every position that scores above random
chance against the motif, which is the right behavior for a genuinely
degenerate binding motif but not for a fixed consensus like `dam`'s `GATC` -
Dam methylase only acts on the exact sequence, so a 3-of-4 near-match carries
no real methylation risk. Measured on a 711nt real gene (mCherry against E.
coli/B. subtilis): motif detection reported 83 hits and drove 72 edits
touching **24.9%** of the gene's codons, but of the 37 `dam` hits, only **2**
were genuine `GATC` sites - the rest were near-matches like `GATG`, `GTTC`,
`TATC`, `CATC`. Three of ESO's other bundled motifs
(`shine_dalgarno`, `sigma70_minus35`, `sigma70_minus10`) are regulatory
elements (a ribosome binding site and promoter boxes), not hypermutable
sites, and are excluded even when motif detection is turned on.

The deeper fix (exact-consensus filtering, or tightening the PSSM threshold
in ESO) is deferred to a follow-up; ESO itself is not modified by DCUB.

To opt in, edit the `HOTSPOT_AVOIDANCE` section of
`app/modules/configuration.yaml`:

```yaml
HOTSPOT_AVOIDANCE:
  COMPUTE_MOTIFS: True
  COMMON_MOTIFS: ["dam", "dcm"]   # restricted to genuine methylation motifs
  RECOMBINATION_MODE: "thorough"  # or "fast" - see eso.detection.dispatch
  SLIPPAGE_MODE: "default"        # or "fast"
```

Translation is unaffected either way this setting is configured; this is a
heads-up on edit volume and false-positive rate, not a correctness concern.

#### Sub-codon slippage limitation

ESO discards slippage-avoidance patterns narrower than 3nt whenever locked
regions are present, so DCUB re-expresses sub-codon repeat units (a 15×`A` run
becomes 5×`AAA`) before handing them over. A repeat too short to yield two
codon-width units cannot be disrupted at codon resolution; those are reported
in the run summary's `warnings` rather than silently skipped.

---

## Contact Details

- Email: bentulila@mail.tau.ac.il  
- Website: https://www.tau.ac.il/~bentulila

---

## Credits

- **Source code:**  
  https://github.com/morikab/DCUB

- **Enable Chrome plugin (accessibility):**  
  Free Chrome plugin – https://enable.co.il/
