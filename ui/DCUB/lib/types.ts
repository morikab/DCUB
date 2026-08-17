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

/** One organism as the backend's `user_input_dict.organisms` map expects it. */
export interface OrganismRequestPayload {
  genome_path: string
  optimized: boolean
  expression_data_type?: string
  expression_file_format?: string
  expression_file_path: string | null
  optimization_priority: number
}

/**
 * The server's raw JSON, before parseOptimizationResponse normalizes it.
 * Only the fields that parser reads are declared.
 *
 * `final_evaluation` is required rather than optional on purpose: the parser
 * dereferences it directly and relies on the surrounding try/catch to turn a
 * missing one into "Invalid response format from server". Making it optional
 * here would force optional chaining and silently produce a result object
 * full of defaults instead of surfacing the bad response.
 */
export interface RawOptimizationResponse {
  final_evaluation: {
    final_sequence?: string
    average_distance_score?: number
    ratio_score?: number
    weakest_link_score?: number
    organisms_dist_scores?: Array<{
      name: string
      is_wanted: boolean
      dist_score: number
    }>
  }
  original_sequence?: string
  processing_time?: number
  timestamp?: string
  hotspot_avoidance?: {
    enabled?: boolean
    sequence_before?: string
    sequence_after?: string
    num_edits?: number
    detected_sites?: {
      recombination?: number
      slippage?: number
      motifs?: number
    }
    warnings?: string[]
  }
}
