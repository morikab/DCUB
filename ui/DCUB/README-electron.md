# DNA Sequence Optimizer - Electron Desktop App

This application can be packaged as a desktop application using Electron.

## Development Setup

1. Install dependencies for both Next.js and Electron:
\`\`\`bash
npm install
\`\`\`

2. Build the Next.js application:
\`\`\`bash
npm run build
npm run export
\`\`\`

3. Run the Electron app in development mode:
\`\`\`bash
npm run electron-dev
\`\`\`

## Building for Production

1. Build the complete Electron application:
\`\`\`bash
npm run build-electron
\`\`\`

This will create distributable packages in the `electron/dist` directory.

## Features in Electron Version

- **Native File Dialogs**: Use system file dialogs for selecting FASTA files, genome files, and CSV files
- **Offline Operation**: Works completely offline once installed
- **Native Network Requests**: Handles CORS issues when connecting to localhost:8000
- **Auto-updater Ready**: Can be extended with auto-update functionality
- **Cross-platform**: Builds for Windows, macOS, and Linux

## API Endpoint

The application sends POST requests to:
\`\`\`
http://localhost:8000/run_modules
\`\`\`

### Request Format

\`\`\`json
{
  "dna_sequence": {
    "content": "ATGCGATCG...", // or null if file uploaded
    "file_name": "sequence.fasta", // or null
    "file_content": "base64_encoded_content" // or null
  },
  "wanted_organisms": [
    {
      "id": "unique_id",
      "genome_path": "/path/to/genome.fasta",
      "priority": 75,
      "expression_data_path": "/path/to/expression.csv"
    }
  ],
  "unwanted_organisms": [
    {
      "id": "unique_id", 
      "genome_path": "/path/to/genome.fasta",
      "priority": 25,
      "expression_data_path": null
    }
  ],
  "advanced_options": {
    "tuning_parameter": 50,
    "optimization_method": "single_codon_diff",
    "cub_index": "CAI"
  },
  "timestamp": "2024-01-20T10:30:00.000Z"
}
\`\`\`

## File Structure

\`\`\`
project/
├── app/                    # Next.js app
├── components/            # React components
├── lib/                   # Utilities and stores
├── electron/              # Electron-specific files
│   ├── main.js           # Main Electron process
│   ├── preload.js        # Preload script
│   ├── package.json      # Electron dependencies
│   └── assets/           # App icons and resources
└── out/                  # Built Next.js app (for Electron)
