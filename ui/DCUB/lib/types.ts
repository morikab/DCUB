export interface Organism {
  id: string
  name: string // Add organism name field
  genomePath: string
  priority: number
  expressionDataPath?: string
  expressionDataFormat?: "csv" | "json"
  expressionDataType?: "protein_abundance" | "mrna_levels"
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
    average_distance_score: number
    ratio_score: number
    weakest_link_score: number
  }
  original_sequence: string
  optimization_parameters: {
    tuning_parameter: number
    optimization_method: string
    cub_index: string
  }
  processing_time: number
  timestamp: string
  organisms_dist_scores: Array<{
    name: string
    is_wanted: boolean
    dist_score: number
  }>
  hotspot_avoidance?: HotspotAvoidanceResult
}

export interface HotspotAvoidanceResult {
  enabled: boolean
  sequence_before: string
  sequence_after: string
  num_edits: number
  detected_sites: { recombination: number; slippage: number; motifs: number }
  warnings: string[]
}
