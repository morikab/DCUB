#!/usr/bin/env bash

# From project root: /Users/.../Igem_TAU_2021
set -euo pipefail

csv_path=$(poetry run python - << 'EOF'
import codonbias, os
print(os.path.join(os.path.dirname(codonbias.__file__), "genetic_code_ncbi.csv"))
EOF
)

eso_data_path=$(poetry run python - << 'EOF'
import os
import eso
print(os.path.join(os.path.dirname(eso.__file__), "data"))
EOF
)

echo "Using genetic_code_ncbi.csv at: $csv_path"
echo "Using eso data directory at: $eso_data_path"

# ESO pulls in DNAChisel, and both DNAChisel and python_codon_tables read data
# files out of their own package directories AT IMPORT TIME - dnachisel's
# biotools/data/complements.csv, and python_codon_tables' os.listdir of its
# tables directory. PyInstaller does not collect package data on its own, so
# without these the frozen server dies during import with a FileNotFoundError
# before it ever binds port 8000, and Electron reports it as "Unable to connect
# to the optimization server". Nothing exercises this in the dev workflow,
# where both packages are read from site-packages.
poetry run pyinstaller \
  --noconfirm \
  --onedir \
  --name fastapi_server \
  --add-data="${csv_path}:codonbias" \
  --add-data="${eso_data_path}:eso/data" \
  --add-data="app/modules/configuration.yaml:modules" \
  --collect-data dnachisel \
  --collect-data python_codon_tables \
  app/api_server.py
  

backend_path="ui/DCUB/backend"
echo "Copying backend executable to: $backend_path"
# Create backend directory if it does not exist and copy the executable there
mkdir -p $backend_path
echo "Removing previous executable if exists from: $backend_path"
rm -rf $backend_path/fastapi_server
cp -r dist/fastapi_server $backend_path/fastapi_server
