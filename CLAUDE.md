# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

DCUB (Differential Codon Usage Bias) — a bioinformatics tool that optimizes a gene's coding sequence for expression in specific microbial environments. Given a target gene and a set of "wanted" and "unwanted" organisms (as GenBank `.gb` files), it outputs a codon-optimized sequence that maximizes expression in wanted organisms while minimizing it in unwanted ones.

The tool ships as an **Electron desktop app** bundling a **Next.js frontend** and a **PyInstaller-packaged FastAPI backend**.

---

## Branching strategy

This is a deliberate split — use it to decide which branch a given task belongs on:

- **`main`** — the canonical branch. Always holds the latest version of the model/backend code (`app/`), the wrapping tool (the DCUB Electron + Next.js app, `ui/DCUB/`), and the website together with all of its artifacts and downloadable builds (`docs/`/website content). Product, model, and app changes belong here.
- **`remote-analysis`** — dedicated to running analyses: testing different inputs, model configurations, etc., and generating the summary graphs/analytics that feed the article and thesis work about the DCUB app. Despite the "remote" name, it also contains local scripts and notebooks for running analysis locally — sometimes on inputs or run outputs that were curated remotely (e.g. on TAU's Power9 cluster), not exclusively remote-cluster work.
  - Analysis on `remote-analysis` must always pull and run against the up-to-date model version from `main`, rather than against a stale or diverged copy of the model code on that branch.

Every new feature, model change, or app fix must be implemented on its own dedicated branch checked out from `main` — never committed directly to `main`. When instructed to push to remote, push the branch and open a PR against `main` on the repo's GitHub.

---

## Commands

### Python backend

```bash
# Install all dependencies (creates .venv/ in-project)
poetry install --with build

# Also install Jupyter/analysis tools
poetry install --with build,analysis

# Run the FastAPI server directly (development)
cd app && poetry run python api_server.py
# Server starts at http://127.0.0.1:8000

# Build the standalone FastAPI executable (output → dist/fastapi_server/, then copied to ui/DCUB/backend/)
./build_fastapi.sh

# Test the built executable
./ui/DCUB/backend/fastapi_server/fastapi_server
```

### Next.js frontend (`ui/DCUB/`)

```bash
cd ui/DCUB

npm run dev          # Next.js dev server on :3000
npm run build        # Production build
npm run lint         # ESLint
npm run electron-dev # Full build + copy standalone + launch Electron
```

### Electron desktop app (`ui/DCUB/electron/`)

```bash
cd ui/DCUB/electron
npm run build-electron   # Package the Electron app (uses electron-builder)
```

No automated test suite exists in the repository.

---

## Architecture

### Backend pipeline (`app/`)

All optimization runs through `app/modules/main.py:run_modules()`. The data model that flows between stages is `models.ModuleInput` (defined in `app/modules/models.py`).

Pipeline stages, each implemented as a `*Module` class with a `run_module()` static/class method:

1. **UserInputModule** (`user_IO/user_input.py`) — Parses `UserInput` (Pydantic model from the API) into `ModuleInput`. Reads GenBank genome files, extracts CDS annotations, computes CAI weights via `codon-bias`, looks up tAI weights, and normalizes organism priorities.

2. **InitiationModule** (`initiation/initiation_main.py`) — Optionally optimizes the first N codons (default: 15, set in `configuration.yaml`) for ribosome binding. Supports three modes: `original` (no-op), `external` (hardcoded reference sequence from config), `weak_folding`.

3. **SequenceFamilyModule** (`sequence_family/sequence_family_main.py`) — When `clusters_count > 1`, clusters wanted organisms by codon usage similarity and splits the problem into sub-problems (one `ModuleInput` per cluster). When `clusters_count == 1` (default), passes through unchanged.

4. **ORFModule** (`ORF/orf_main.py`) — Core codon optimization. Runs separately for CAI and/or tAI depending on `orf_optimization_cub_index`. Optimization methods:
   - `single_codon_*` — codon-by-codon greedy substitution
   - `zscore_single_aa_*` / `zscore_bulk_aa_*` — z-score–based iterative optimization; generates multiple random starting permutations and picks the best
   - `single_wanted_organism` — optimize for a single target

5. **EvaluationModule** (`evaluation/evaluation.py`) — Scores each candidate sequence using three metrics: `average_distance`, `weakest_link`, and `ratio`. Scores are z-score normalized against organism-level statistics. The best result across all candidates (CAI vs tAI, multiple random seeds) is selected by the configured `evaluation_score`.

6. **UserOutputModule** (`user_IO/user_output.py`) — Writes the optimized sequence as a FASTA file and bundles it with the run log into a `communique_results.zip`.

### Key shared types (`app/modules/models.py`)

- `UserInput` (Pydantic) — API request body; organism genomes are keyed by name as `Dict[str, OrganismRequest]`
- `ModuleInput` (dataclass) — internal representation passed between modules
- `Organism` — holds computed CAI/tAI profiles, codon frequencies, scores, and priority weight
- `ORFOptimizationMethod`, `ORFOptimizationCubIndex`, `EvaluationScore` — enums controlling optimization behavior

### API (`app/api_server.py`)

Single endpoint: `POST /run-modules` — accepts `RunModulesRequest` (which wraps `UserInput` under the alias `user_input_dict`) and returns a `RunModulesResponse` with the full run summary dict, including `zip_output_file_path` and `final_evaluation`.

### Frontend (`ui/DCUB/`)

- **Next.js 15 + React 19 + TypeScript**, styled with Tailwind CSS v4 and Radix UI primitives
- Global state managed with **Zustand** (`lib/store.ts`), persisted to `localStorage`
- Single page (`app/page.tsx`) — collects DNA sequence, wanted/unwanted organisms, and advanced options, then POSTs to `http://localhost:8000/run-modules`
- File I/O (genome files, expression data, result download) goes through Electron IPC when running in Electron, with a browser fallback (`lib/electron-utils.ts`, `components/electron-file-handler.tsx`)

### Electron (`ui/DCUB/electron/`)

`main.js` starts two child processes on launch:
1. The Next.js standalone server (on port 3000), spawned via Node.js
2. The FastAPI executable (on port 8000), spawned from `../backend/fastapi_server/`

In development (`isDev = !app.isPackaged`), the backend is loaded from `../backend/`; in the packaged app it comes from `process.resourcesPath`.

### Configuration (`app/modules/configuration.yaml`)

Runtime tuning knobs read by `Configuration.get_config()`:
- `INITIATION.NUMBER_OF_CODONS_TO_OPTIMIZE` — how many start codons are handled by the initiation module (default 15)
- `ORF.ZSCORE_MAX_ITERATIONS`, `ORF.ZSCORE_INITIAL_PERMUTATIONS_NUM` — z-score optimization loop controls

### Analysis (`analysis/` directory, primarily on `remote-analysis`)

The `analysis/` directory contains Jupyter notebooks used to run the model across many organism/method combinations and evaluate results. These are typically executed on **TAU's Power9 cluster** (`power9login` host), which requires connecting through Tel Aviv University's VPN first.

Install the analysis extras before running notebooks:
```bash
poetry install --with build,analysis
poetry run jupyter lab
```

Key notebooks in `analysis/orf_model_analysis/`:
- `arabidopsis.ipynb` — sweeps all pairwise organism combinations from `analysis/example_data/arabidopsis_microbiome/` across multiple optimization methods, calling `run_modules()` directly (not via the API)
- `mcherry_variations.ipynb` — analyzes mCherry sequence variants
- `_plot_panel_a_alternatives.py` — standalone script generating panel-A figure alternatives for the arabidopsis dataset; outputs to `analysis/results/arabidopsis/figures/`

Notebooks add the project root to `sys.path` and import `from modules.main import run_modules` directly, so they must be launched from within the repository tree. Results land in `analysis/orf_model_analysis/results/`.

The `analysis/example_data/arabidopsis_microbiome/` directory holds the reference genome set (`.gbff` files) used by the analysis notebooks. Genomes are downloaded from NCBI RefSeq via Entrez and decompressed before use (see `arabidopsis.ipynb`).

---

## Figure Organization

All analysis-generated figures are written to `analysis/results/<dataset>/figures/`. The datasets currently are:

- `analysis/results/arabidopsis/figures/` — figures from arabidopsis microbiome analysis
- `analysis/results/homo_sapiens/figures/` — figures from homo sapiens analysis

**Rules:**
- Scripts and notebooks must always write figure output to `analysis/results/<dataset>/figures/`, never to the repo root or any other location.
- Figure files must never be placed at the repo root. If you find figure files (`.svg`, `.pdf`, `.png`) at the root, move them to the appropriate `analysis/results/<dataset>/figures/` directory.
- Draft/alternative figures (files prefixed with `_sample_`, `_supplementary_`, etc.) belong in the same `analysis/results/<dataset>/figures/` folder as the finalized versions — not at the root.
- Analysis scripts that generate figures (e.g., `_plot_panel_a_alternatives.py`) belong in `analysis/orf_model_analysis/`, not at the repo root.
- When a figure is finalized for the article, copy it to `figures/article/` in the thesis folder at `/Users/shimka/Documents/Moran's Thesis/figures/article/`. Do not keep separate "article copy" versions inside the repo.
- `docs/assets/` is for website/documentation assets only (logos, UI screenshots), not for article figures.

---

### Build notes

- `build_fastapi.sh` calls PyInstaller and bundles the `codonbias` genetic code CSV as a data file; the path is resolved dynamically from the installed package.
- `ui/DCUB/scripts/copy-standalone-next.js` copies the Next.js `.next/standalone` output into the Electron app folder before packaging.
- Python requirement: `~3.9` (pinned in `pyproject.toml`).
