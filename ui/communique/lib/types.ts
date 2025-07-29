export interface Organism {
  id: string
  name: string // Add organism name field
  genomePath: string
  priority: number
  expressionDataPath?: string
}

export interface ValidationResult {
  isValid: boolean
  error?: string
  errors?: string[]
}

export interface SubmissionData {
  dnaSequence: string
  sequenceFile: File | null
  wantedOrganisms: Organism[]
  unwantedOrganisms: Organism[]
}

export interface OptimizationResult {
  optimized_sequence: string
  evaluation_scores: {
    cai_score: number
    gc_content: number
    codon_usage_bias: number
  }
  original_sequence: string
  optimization_parameters: {
    tuning_parameter: number
    optimization_method: string
    cub_index: string
  }
  processing_time: number
  timestamp: string
}
