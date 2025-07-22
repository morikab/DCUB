# Commuique

A DNA sequence optimization tool that uses GenBank genome files to configure optimization parameters for target organisms.

## Features

- **DNA Sequence Input**: Manual FASTA entry or file upload
- **Organism Configuration**: Add wanted and unwanted organisms with GenBank genome files (.gb/.gbf)
- **Priority System**: Set optimization priorities (1-100) for each organism
- **Expression Data**: Optional CSV files with expression data
- **Advanced Options**: Tuning parameters, optimization methods, and CUB index selection
- **Backend Integration**: POST requests to localhost:8000/run_modules endpoint

## Getting Started

1. Install dependencies:
\`\`\`bash
npm install
\`\`\`

2. Run the development server:
\`\`\`bash
npm run dev
\`\`\`

3. Open [http://localhost:3000](http://localhost:3000) in your browser.

## GenBank File Requirements

- Genome files must be in GenBank format (.gb or .gbf extensions)
- Files should contain complete genome sequences with annotations
- Both wanted and unwanted organisms require GenBank genome files

## API Endpoint

The application sends POST requests to:
\`\`\`
http://localhost:8000/run_modules
\`\`\`

### Request Format

\`\`\`json
{
  "dna_sequence": {
    "content": "ATGCGATCG...",
    "file_name": "sequence.fasta",
    "file_content": "base64_encoded_content"
  },
  "wanted_organisms": [
    {
      "id": "unique_id",
      "genome_path": "/path/to/genome.gb",
      "priority": 75,
      "expression_data_path": "/path/to/expression.csv"
    }
  ],
  "unwanted_organisms": [
    {
      "id": "unique_id", 
      "genome_path": "/path/to/genome.gb",
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

## Development

Make sure your backend server has CORS enabled for localhost:3000 during development.

## Build

\`\`\`bash
npm run build
\`\`\`

## License

MIT
