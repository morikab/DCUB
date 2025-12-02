#!/usr/bin/env zsh

# From project root: /Users/.../Igem_TAU_2021

csv_path=$(python - << 'EOF'
import codonbias, os
print(os.path.join(os.path.dirname(codonbias.__file__), "genetic_code_ncbi.csv"))
EOF
)

echo "Using genetic_code_ncbi.csv at: $csv_path"

pyinstaller \
  --onefile \
  --name fastapi_server \
  --add-data="${csv_path}:codonbias" \
  --add-data="app/modules/configuration.yaml:modules" \
  app/api_server.py
  # --add-data "app/modules/configuration.yaml:modules" \
  

