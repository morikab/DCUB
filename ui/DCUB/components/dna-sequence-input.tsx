"use client"

import type React from "react"

import { useState, useRef } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Textarea } from "@/components/ui/textarea"
import { Button } from "@/components/ui/button"
import { Label } from "@/components/ui/label"
import { Progress } from "@/components/ui/progress"
import { Upload, FileText, X, Download, Info } from "lucide-react"
import { useOptimizationStore } from "@/lib/store"
import { validateFastaSequence } from "@/lib/validation"

export function DnaSequenceInput() {
  const [uploadProgress, setUploadProgress] = useState(0)
  const [isUploading, setIsUploading] = useState(false)
  const [validationError, setValidationError] = useState("")
  const [sequenceInfo, setSequenceInfo] = useState<{ title: string; length: number } | null>(null)
  const [displayText, setDisplayText] = useState("") // Store the original FASTA text for display
  const fileInputRef = useRef<HTMLInputElement>(null)

  const { dnaSequence, sequenceFile, setDnaSequence, setSequenceFile } = useOptimizationStore()

  // Function to parse FASTA format and extract sequence only
  const parseFastaSequence = (fastaText: string): { sequence: string; title: string } => {
    const lines = fastaText.trim().split("\n")
    let title = ""
    let sequence = ""

    for (const line of lines) {
      const trimmedLine = line.trim()
      if (trimmedLine.startsWith(">")) {
        // Extract title (remove the '>' character)
        title = trimmedLine.substring(1).trim()
      } else if (trimmedLine) {
        // Accumulate sequence lines, removing any whitespace
        sequence += trimmedLine.replace(/\s/g, "").toUpperCase()
      }
    }

    return { sequence, title }
  }

  const handleTextChange = (value: string) => {
    setSequenceFile(null)
    setSequenceInfo(null)
    setDisplayText(value) // Always store the original input for display

    if (value.trim()) {
      const validation = validateFastaSequence(value)
      if (validation.isValid) {
        // Parse FASTA and extract sequence
        const { sequence, title } = parseFastaSequence(value)
        setDnaSequence(sequence) // Store only the sequence for API
        setSequenceInfo({
          title: title || "Untitled sequence",
          length: sequence.length,
        })
        setValidationError("")
      } else {
        // For invalid FASTA, check if it's just a raw sequence
        const cleanSequence = value.replace(/\s/g, "").toUpperCase()
        const validBases = /^[ATGCNRYSWKMBDHV]*$/

        if (validBases.test(cleanSequence) && cleanSequence.length > 0) {
          // It's a raw sequence without FASTA header
          setDnaSequence(cleanSequence)
          setSequenceInfo({
            title: "Raw sequence (no header)",
            length: cleanSequence.length,
          })
          setValidationError("")
        } else {
          setDnaSequence("") // Don't store invalid sequences
          setValidationError(validation.error || "")
        }
      }
    } else {
      setDnaSequence("")
      setValidationError("")
    }
  }

  const handleFileUpload = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0]
    if (!file) return

    // Validate file type
    if (!file.name.toLowerCase().endsWith(".fasta") && !file.name.toLowerCase().endsWith(".fa")) {
      setValidationError("Please upload a FASTA file (.fasta or .fa)")
      return
    }

    // Check file size (10MB limit)
    if (file.size > 10 * 1024 * 1024) {
      setValidationError("File size must be less than 10MB")
      return
    }

    setIsUploading(true)
    setUploadProgress(0)
    setValidationError("")
    setSequenceInfo(null)

    try {
      // Simulate upload progress
      const progressInterval = setInterval(() => {
        setUploadProgress((prev) => {
          if (prev >= 90) {
            clearInterval(progressInterval)
            return 90
          }
          return prev + 10
        })
      }, 100)

      // Read file content
      const text = await file.text()

      // Validate FASTA content
      const validation = validateFastaSequence(text)
      if (!validation.isValid) {
        setValidationError(validation.error || "Invalid FASTA format")
        setIsUploading(false)
        setUploadProgress(0)
        return
      }

      // Parse FASTA and extract sequence
      const { sequence, title } = parseFastaSequence(text)

      // Complete upload
      setUploadProgress(100)
      setSequenceFile(file)
      setDnaSequence(sequence) // Store only the sequence for API
      setDisplayText("") // Clear manual entry when file is uploaded
      setSequenceInfo({
        title: title || file.name,
        length: sequence.length,
      })

      setTimeout(() => {
        setIsUploading(false)
        setUploadProgress(0)
      }, 500)
    } catch (error) {
      setValidationError("Error reading file. Please try again.")
      setIsUploading(false)
      setUploadProgress(0)
    }
  }

  const removeFile = () => {
    setSequenceFile(null)
    setSequenceInfo(null)
    setValidationError("")
    setDnaSequence("")
    setDisplayText("")
    if (fileInputRef.current) {
      fileInputRef.current.value = ""
    }
  }

  const downloadSampleFasta = () => {
    const sampleFasta = `>Sample DNA Sequence for DCUB Optimization
ATGAAAGTTCTGTTCCAGGGCCCGCCCGCGCCGCTGCTGCTGCTGCTGCTGCTGCTGCTG
CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG
CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG
CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGTAG`

    const blob = new Blob([sampleFasta], { type: "text/plain" })
    const url = URL.createObjectURL(blob)
    const a = document.createElement("a")
    a.href = url
    a.download = "sample_sequence.fasta"
    document.body.appendChild(a)
    a.click()
    document.body.removeChild(a)
    URL.revokeObjectURL(url)
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <FileText className="w-5 h-5" />
          DNA Sequence Input
        </CardTitle>
        <CardDescription>Enter your DNA sequence manually or upload a FASTA file for optimization</CardDescription>
      </CardHeader>
      <CardContent className="space-y-4">
        {/* Manual Input */}
        <div className="space-y-2">
          <div className="flex items-center justify-between">
            <Label htmlFor="dna-sequence">Manual Entry (FASTA Format or Raw Sequence)</Label>
            <Button
              variant="ghost"
              size="sm"
              onClick={downloadSampleFasta}
              className="text-blue-600 hover:text-blue-800"
            >
              <Download className="w-4 h-4 mr-1" />
              Download Sample
            </Button>
          </div>
          <Textarea
            id="dna-sequence"
            placeholder={`>sequence_name
ATGCGATCGATCGATCGATCG...

Or just enter raw sequence:
ATGCGATCGATCGATCGATCG...`}
            value={displayText} // Show the original input
            onChange={(e) => handleTextChange(e.target.value)}
            disabled={!!sequenceFile}
            className="min-h-32 font-mono text-sm"
          />
        </div>

        {/* File Upload */}
        <div className="space-y-2">
          <Label>Or Upload FASTA File</Label>
          <div className="border-2 border-dashed border-gray-300 rounded-lg p-6 text-center hover:border-gray-400 transition-colors">
            {sequenceFile ? (
              <div className="space-y-2">
                <div className="flex items-center justify-center gap-2 text-green-600">
                  <FileText className="w-5 h-5" />
                  <span className="font-medium">{sequenceFile.name}</span>
                  <Button variant="ghost" size="sm" onClick={removeFile} className="text-red-500 hover:text-red-700">
                    <X className="w-4 h-4" />
                  </Button>
                </div>
                <p className="text-sm text-gray-500">File size: {(sequenceFile.size / 1024).toFixed(1)} KB</p>
              </div>
            ) : (
              <>
                <Upload className="w-8 h-8 mx-auto mb-2 text-gray-400" />
                <p className="text-gray-600 mb-2">Click to upload or drag and drop your FASTA file</p>
                <Button
                  variant="outline"
                  onClick={() => fileInputRef.current?.click()}
                  disabled={isUploading || !!displayText.trim()}
                >
                  <Upload className="w-4 h-4 mr-2" />
                  Choose File
                </Button>
              </>
            )}
          </div>

          <input ref={fileInputRef} type="file" accept=".fasta,.fa" onChange={handleFileUpload} className="hidden" />

          {isUploading && (
            <div className="space-y-2">
              <div className="flex justify-between text-sm">
                <span>Uploading...</span>
                <span>{uploadProgress}%</span>
              </div>
              <Progress value={uploadProgress} className="w-full" />
            </div>
          )}
        </div>

        {/* Sequence Information Display */}
        {sequenceInfo && (
          <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
            <div className="flex items-start gap-2">
              <Info className="w-5 h-5 text-blue-600 mt-0.5 flex-shrink-0" />
              <div className="space-y-2">
                <h4 className="font-medium text-blue-900">Sequence Information</h4>
                <div className="text-sm text-blue-800 space-y-1">
                  <p>
                    <strong>Title:</strong> {sequenceInfo.title}
                  </p>
                  <p>
                    <strong>Length:</strong> {sequenceInfo.length.toLocaleString()} nucleotides
                  </p>
                  <p>
                    <strong>Status:</strong> Ready for optimization (sequence content extracted)
                  </p>
                </div>
              </div>
            </div>
          </div>
        )}

        {/* Extracted Sequence Preview (only show if different from display) */}
        {dnaSequence && displayText && displayText !== dnaSequence && (
          <div className="space-y-2">
            <Label>Extracted Sequence (will be sent to API)</Label>
            <div className="bg-gray-50 p-3 rounded-md border max-h-32 overflow-y-auto">
              <p className="text-xs font-mono text-gray-700 break-all">
                {dnaSequence.length > 200
                  ? `${dnaSequence.substring(0, 200)}... (${dnaSequence.length} total nucleotides)`
                  : dnaSequence}
              </p>
            </div>
          </div>
        )}

        {/* Validation Error */}
        {validationError && (
          <div className="text-red-600 text-sm bg-red-50 p-3 rounded-md border border-red-200">{validationError}</div>
        )}

        {/* Info */}
        <div className="text-sm text-gray-500 bg-gray-50 p-3 rounded-md">
          <p>
            <strong>Input Processing:</strong>
          </p>
          <ul className="list-disc list-inside mt-1 space-y-1">
            <li>Your original input (including headers) is preserved in the text box</li>
            <li>Only the DNA sequence content is extracted and sent for optimization</li>
            <li>Supports both FASTA format and raw sequence input</li>
            <li>Headers and formatting are automatically parsed but not included in processing</li>
            <li>Sequences are cleaned and validated before optimization</li>
            <li>File size limit: 10MB</li>
          </ul>
        </div>
      </CardContent>
    </Card>
  )
}
