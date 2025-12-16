import type { ValidationResult, SubmissionData } from "./types"

export function validateFastaSequence(sequence: string): ValidationResult {
  if (!sequence.trim()) {
    return { isValid: false, error: "DNA sequence is required" }
  }

  const lines = sequence.trim().split("\n")

  let sequenceContent: string

  // Check if it starts with a FASTA header
  if (lines[0].startsWith(">")) {
    // FASTA format: extract sequence from lines after the header
    sequenceContent = lines.slice(1).join("").replace(/\s/g, "").toUpperCase()
  } else {
    // Plain sequence format: use all lines as sequence
    sequenceContent = lines.join("").replace(/\s/g, "").toUpperCase()
  }

  // Validate sequence content
  const validBases = /^[ATGCNRYSWKMBDHV]*$/

  if (!validBases.test(sequenceContent)) {
    return {
      isValid: false,
      error: "Sequence contains invalid characters. Only DNA bases (A, T, G, C) and IUPAC codes are allowed",
    }
  }

  if (sequenceContent.length === 0) {
    return { isValid: false, error: "Sequence cannot be empty" }
  }

  return { isValid: true }
}

export function validateGenBankPath(genomePath: string): ValidationResult {
  if (!genomePath.trim()) {
    return { isValid: false, error: "GenBank genome file path is required" }
  }

  const path = genomePath.trim().toLowerCase()
  if (!path.endsWith(".gb") && !path.endsWith(".gbf")) {
    return {
      isValid: false,
      error: "Genome file must be in GenBank format (.gb or .gbf extension)",
    }
  }

  return { isValid: true }
}

export function validateOrganism(genomePath: string, priority: number): ValidationResult {
  const errors: string[] = []

  // Validate GenBank file path
  const pathValidation = validateGenBankPath(genomePath)
  if (!pathValidation.isValid) {
    errors.push(pathValidation.error!)
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
