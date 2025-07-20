export interface Organism {
  id: string
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
