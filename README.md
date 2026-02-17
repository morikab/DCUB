# Welcome to DCUB – Differential Codon Usage Bias

DCUB is a software tool for designing genetic sequences that are tailored to specific microbial environments.  
With only a few steps, you can obtain a microbiome‑specific version of your gene of interest.

Please follow the steps below to install and run the tool.

---

## Prerequisites

A local Node.js / Python environment if building from source.  

---

## Installation Guide (Pre‑built binaries)

Pre‑compiled versions of DCUB for the major operating systems are available from the project website:

- https://www.tau.ac.il/~bentulila

Download the installer or archive that matches your operating system and follow the platform‑specific instructions provided on the website.

After installation, follow the usage instructions from the user guide packaged with the distribution.

---

## Building DCUB from Source

If you prefer to build the tool from source (for development or customization), use the following workflow:

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
  Free Chrome plugin – https://www.enable.co.il/tos/
