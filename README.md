# Welcome to DCUB – Differential Codon Usage Bias

DCUB is a software tool for designing genetic sequences that are tailored to specific microbial environments.  
With only a few steps, you can obtain a microbiome‑specific version of your gene of interest.

Please follow the steps below to install and run the tool.

---

## Prerequisites

A local Node.js / Python environment if building from source. Check the exact version constraints before you start:

- **Python** — pinned in [`pyproject.toml`](./pyproject.toml) (`python = "~3.9"`); install via [Poetry](https://python-poetry.org/).
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

## Contact Details

- Email: bentulila@mail.tau.ac.il  
- Website: https://www.tau.ac.il/~bentulila

---

## Credits

- **Source code:**  
  https://github.com/morikab/DCUB

- **Enable Chrome plugin (accessibility):**  
  Free Chrome plugin – https://enable.co.il/
