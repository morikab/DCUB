"use client"

import type React from "react"

import { useState, useRef } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Textarea } from "@/components/ui/textarea"
import { Button } from "@/components/ui/button"
import { Label } from "@/components/ui/label"
import { Progress } from "@/components/ui/progress"
import { Upload, FileText, X, Download } from "lucide-react"
import { useOptimizationStore } from "@/lib/store"
import { validateFastaSequence } from "@/lib/validation"

export function DnaSequenceInput() {
  const [uploadProgress, setUploadProgress] = useState(0)
  const [isUploading, setIsUploading] = useState(false)
  const [validationError, setValidationError] = useState("")
  const fileInputRef = useRef<HTMLInputElement>(null)

  const { dnaSequence, sequenceFile, setDnaSequence, setSequenceFile } = useOptimizationStore()

  const handleTextChange = (value: string) => {
    setDnaSequence(value)
    setSequenceFile(null)

    if (value.trim()) {
      const validation = validateFastaSequence(value)
      setValidationError(validation.isValid ? "" : validation.error || "")
    } else {
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
        return
      }

      // Complete upload
      setUploadProgress(100)
      setSequenceFile(file)
      setDnaSequence("")

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
    setValidationError("")
    if (fileInputRef.current) {
      fileInputRef.current.value = ""
    }
  }

  const downloadSampleFasta = () => {
    const sampleFasta = `>Sample DNA Sequence
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
            <Label htmlFor="dna-sequence">Manual Entry (FASTA Format)</Label>
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
ATGCGATCGATCGATCGATCG...`}
            value={dnaSequence}
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
                  disabled={isUploading || !!dnaSequence.trim()}
                >
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

        {/* Validation Error */}
        {validationError && (
          <div className="text-red-600 text-sm bg-red-50 p-3 rounded-md border border-red-200">{validationError}</div>
        )}

        {/* Info */}
        <div className="text-sm text-gray-500 bg-gray-50 p-3 rounded-md">
          <p>
            <strong>FASTA Format Requirements:</strong>
          </p>
          <ul className="list-disc list-inside mt-1 space-y-1">
            <li>Must start with a header line beginning with {">"}</li>
            <li>Sequence should contain only valid DNA bases (A, T, G, C, N)</li>
            <li>IUPAC nucleotide codes are supported</li>
            <li>Multiple sequences are supported</li>
            <li>File size limit: 10MB</li>
          </ul>
        </div>
      </CardContent>
    </Card>
  )
}
