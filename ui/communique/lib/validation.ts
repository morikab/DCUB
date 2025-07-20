import type { ValidationResult, SubmissionData } from "./types"

export function validateFastaSequence(sequence: string): ValidationResult {
  if (!sequence.trim()) {
    return { isValid: false, error: "DNA sequence is required" }
  }

  const lines = sequence.trim().split("\n")

  // Check if it starts with a header
  if (!lines[0].startsWith(">")) {
    return { isValid: false, error: "FASTA sequence must start with a header line (>)" }
  }

  // Validate sequence content
  const sequenceLines = lines.slice(1).join("").replace(/\s/g, "").toUpperCase()
  const validBases = /^[ATGCNRYSWKMBDHV]*$/

  if (!validBases.test(sequenceLines)) {
    return {
      isValid: false,
      error: "Sequence contains invalid characters. Only DNA bases (A, T, G, C) and IUPAC codes are allowed",
    }
  }

  if (sequenceLines.length === 0) {
    return { isValid: false, error: "Sequence cannot be empty" }
  }

  return { isValid: true }
}

export function validateOrganism(genomePath: string, priority: number): ValidationResult {
  const errors: string[] = []

  if (!genomePath.trim()) {
    errors.push("Genome file path is required")
  }

  if (priority < 1 || priority > 100) {
    errors.push("Priority must be between 1 and 100")
  }

  return {
    isValid: errors.length === 0,
    errors,
  }
}

export function validateSubmission(data: SubmissionData): ValidationResult {
  const errors: string[] = []

  // Validate DNA sequence
  if (!data.dnaSequence.trim() && !data.sequenceFile) {
    errors.push("DNA sequence is required (either manual entry or file upload)")
  }

  if (data.dnaSequence.trim()) {
    const sequenceValidation = validateFastaSequence(data.dnaSequence)
    if (!sequenceValidation.isValid) {
      errors.push(`DNA sequence: ${sequenceValidation.error}`)
    }
  }

  // Validate organisms
  if (data.wantedOrganisms.length === 0 && data.unwantedOrganisms.length === 0) {
    errors.push("At least one wanted or unwanted organism must be specified")
  }
  // Validate each organism
  ;[...data.wantedOrganisms, ...data.unwantedOrganisms].forEach((organism, index) => {
    const validation = validateOrganism(organism.genomePath, organism.priority)
    if (!validation.isValid) {
      validation.errors?.forEach((error) => {
        errors.push(`Organism ${index + 1}: ${error}`)
      })
    }
  })

  return {
    isValid: errors.length === 0,
    errors,
  }
}
