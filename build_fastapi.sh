#!/usr/bin/env zsh

# From project root: /Users/.../Igem_TAU_2021

csv_path=$(python - << 'EOF'
import codonbias, os
print(os.path.join(os.path.dirname(codonbias.__file__), "genetic_code_ncbi.csv"))
EOF
)

echo "Using genetic_code_ncbi.csv at: $csv_path"

pyinstaller \
  --onedir \
  --name fastapi_server \
  --add-data="${csv_path}:codonbias" \
  --add-data="app/modules/configuration.yaml:modules" \
  app/api_server.py
  

backend_path="ui/DCUB/backend"
echo "Copying backend executable to: $backend_path"
# Create backend directory if it does not exist and copy the executable there
mkdir -p $backend_path
cp -r dist/fastapi_server $backend_path/fastapi_server
