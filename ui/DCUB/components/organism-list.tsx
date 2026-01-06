"use client"

import type React from "react"

import { useState, useRef } from "react"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { Badge } from "@/components/ui/badge"
import { Progress } from "@/components/ui/progress"
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"
import { Plus, Trash2, Download, Upload, FileText, X } from "lucide-react"
import { useOptimizationStore } from "@/lib/store"
import type { Organism } from "@/lib/types"

const downloadExpressionSample = (
  format: "csv" | "json" = "csv", 
  dataType: "protein_abundance" | "mrna_levels" = "mrna_levels"
) => {
  let content = ""
  const fileExtension = format
  let mimeType = `text/${format}`

  if (format === 'json') {
    mimeType = 'application/json'
    if (dataType === 'protein_abundance') {
      content = JSON.stringify([
        {"gene": "gene_001", "protein_abundance": 1250.5},
        {"gene": "gene_002", "protein_abundance": 890.3},
        {"gene": "gene_003", "protein_abundance": 2100.8},
      ], null, 2);
    } else { // mrna_levels
      content = JSON.stringify([
        {"gene": "gene_001", "mRNA_level": 45.2},
        {"gene": "gene_002", "mRNA_level": 28.7},
        {"gene": "gene_003", "mRNA_level": 67.9},
      ], null, 2);
    }
  } else { // csv
    if (dataType === 'protein_abundance') {
      content = `gene,protein_abundance
gene_001,1250.5
gene_002,890.3
gene_003,2100.8`
    } else { // mrna_levels
      content = `gene,mRNA_level
gene_001,1250.5
gene_002,890.3
gene_003,2100.8`
    }
  }

  const blob = new Blob([content], { type: mimeType })
  const url = URL.createObjectURL(blob)
  const a = document.createElement("a")
  a.href = url
  a.download = `sample_expression_${dataType}.${fileExtension}`
  document.body.appendChild(a)
  a.click()
  document.body.removeChild(a)
  URL.revokeObjectURL(url)
}

const downloadSampleGenBank = () => {
  const sampleGenBank = `LOCUS       SAMPLE_GENOME           1000 bp    DNA     linear   BCT 01-JAN-2024
DEFINITION  Sample GenBank file for Commuique optimization.
ACCESSION   SAMPLE001
VERSION     SAMPLE001.1
KEYWORDS    sample, genome, optimization.
SOURCE      Sample organism
  ORGANISM  Sample organism
            Bacteria; Sample phylum; Sample class; Sample order;
            Sample family; Sample genus.
FEATURES             Location/Qualifiers
     source          1..1000
                     /organism="Sample organism"
                     /mol_type="genomic DNA"
     gene            1..300
                     /gene="sampleGene1"
     CDS             1..300
                     /gene="sampleGene1"
                     /product="sample protein 1"
                     /translation="MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"
ORIGIN      
        1 atgaaacagc ataaagcaat gattgtcgct ttgattgtga tttgtattac tgctgttgtt
       61 gctgctttgg ttactcgtaa agatctttgt gaagtgcata ttcgtactgg tcagactgaa
      121 gtggctgtat tttagatcga tgcacgtatt ggtcagattt tgcgtgaagt tgctggtgaa
      181 cgtggtgaag ttgctggtga acgtggtgaag ttgctggtg aacgtggtga agttgctggt
      241 gaacgtggtg aagttgctgg tgaacgtggtgaa gttgctggtg aacgtgg tgaagttgct
      301 ggtgaacgtg gtgaagttgc tggtgaacgt ggtgaagttg ctggtgaacg tggtgaagtt
      361 gctggtgaac gtggtgaagt tgctggtgaa cgtggtgaag ttgctggtga acgtggtgaa
      421 gttgctggtg aacgtggtga agttgctggt gaacgtggtg aagttgctgg tgaacgtggt
      481 gaagttgctg gtgaacgtgg tgaagttgct ggtgaacgtg gtgaagttgc tggtgaacgt
      541 ggtgaagttg ctggtgaacg tggtgaagtt gctggtgaac gtggtgaagt tgctggtgaa
      601 cgtggtgaag ttgctggtga acgtggtgaa gttgctggtg aacgtggtga agttgctggt
      661 gaacgtggtg aagttgctgg tgaacgtggt gaagttgctg gtgaacgtgg tgaagttgct
      721 ggtgaacgtg gtgaagttgc tggtgaacgt ggtgaagttg ctggtgaacg tggtgaagtt
      781 gctggtgaac gtggtgaagt tgctggtgaa cgtggtgaag ttgctggtga acgtggtgaa
      841 gttgctggtg aacgtggtga agttgctggt gaacgtggtg aagttgctgg tgaacgtggt
      901 gaagttgctg gtgaacgtgg tgaagttgct ggtgaacgtg gtgaagttgc tggtgaacgt
      961 ggtgaagttg ctggtgaacg tggtgaagtt gctggtgaac
//`

  const blob = new Blob([sampleGenBank], { type: "text/plain" })
  const url = URL.createObjectURL(blob)
  const a = document.createElement("a")
  a.href = url
  a.download = "sample_genome.gb"
  document.body.appendChild(a)
  a.click()
  document.body.removeChild(a)
  URL.revokeObjectURL(url)
}

interface OrganismListProps {
  type: "wanted" | "unwanted"
}

export function OrganismList({ type }: OrganismListProps) {
  const [newOrganism, setNewOrganism] = useState<Partial<Organism>>({
    name: "", // Add name field
    genomePath: "",
    priority: undefined,
    expressionDataPath: "",
    expressionDataFormat: "csv",
    expressionDataType: "protein_abundance",
  })

  const {
    wantedOrganisms,
    unwantedOrganisms,
    addWantedOrganism,
    addUnwantedOrganism,
    removeWantedOrganism,
    removeUnwantedOrganism,
    updateWantedOrganism,
    updateUnwantedOrganism,
  } = useOptimizationStore()

  const organisms = type === "wanted" ? wantedOrganisms : unwantedOrganisms
  const addOrganism = type === "wanted" ? addWantedOrganism : addUnwantedOrganism
  const removeOrganism = type === "wanted" ? removeWantedOrganism : removeUnwantedOrganism
  const updateOrganism = type === "wanted" ? updateWantedOrganism : updateUnwantedOrganism

  const calculateDefaultPriority = () => {
    if (organisms.length === 0) return 50
    const validPriorities = organisms.filter((org) => org.priority !== undefined).map((org) => org.priority!)
    if (validPriorities.length === 0) return 50
    return Math.round(validPriorities.reduce((sum, p) => sum + p, 0) / validPriorities.length)
  }

  const handleAddOrganism = () => {
    if (!newOrganism.genomePath?.trim() || !newOrganism.name?.trim()) return

    const organism: Organism = {
      id: Date.now().toString(),
      name: newOrganism.name.trim(),
      genomePath: newOrganism.genomePath.trim(),
      priority: newOrganism.priority || calculateDefaultPriority(),
      expressionDataPath: newOrganism.expressionDataPath?.trim() || undefined,
      expressionDataFormat: newOrganism.expressionDataFormat,
      expressionDataType: newOrganism.expressionDataType,
    }

    addOrganism(organism)
    setNewOrganism({
      name: "",
      genomePath: "",
      priority: undefined,
      expressionDataPath: "",
      expressionDataFormat: "csv",
      expressionDataType: "protein_abundance",
    })
  }

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between">
        <div>
          <h3 className="text-lg font-semibold">{type === "wanted" ? "Wanted" : "Unwanted"} Organisms</h3>
          <p className="text-sm text-gray-600">
            Configure organisms to {type === "wanted" ? "optimize for" : "avoid during optimization"}
          </p>
        </div>
        <Badge variant={type === "wanted" ? "default" : "destructive"}>
          {organisms.length} organism{organisms.length !== 1 ? "s" : ""}
        </Badge>
      </div>

      {/* Existing Organisms */}
      <div className="space-y-3">
        {organisms.map((organism) => (
          <OrganismCard
            key={organism.id}
            organism={organism}
            onUpdate={(updated) => updateOrganism(organism.id, updated)}
            onRemove={() => removeOrganism(organism.id)}
          />
        ))}
      </div>

      {/* Add New Organism */}
      <Card className="border-dashed border-2">
        <CardHeader className="pb-4">
          <CardTitle className="text-base flex items-center gap-2">
            <Plus className="w-4 h-4" />
            Add New Organism
          </CardTitle>
        </CardHeader>
        <CardContent className="space-y-4">
          <div className="space-y-2 mb-4">
            <Label htmlFor="organism-name">Organism Name *</Label>
            <Input
              id="organism-name"
              placeholder="e.g., Escherichia-coli, Bacillus-subtilis"
              value={newOrganism.name || ""}
              onChange={(e) => setNewOrganism((prev) => ({ ...prev, name: e.target.value }))}
            />
          </div>
          <div className="space-y-2">
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
              <div className="space-y-2">
                <div className="flex items-center justify-between">
                  <Label htmlFor="genome-path">GenBank Genome File (.gb/.gbff) *</Label>
                  <Button
                    variant="ghost"
                    size="sm"
                    onClick={downloadSampleGenBank}
                    className="text-blue-600 hover:text-blue-800"
                  >
                    <Download className="w-4 h-4 mr-1" />
                    Sample .gb
                  </Button>
                </div>
                <div className="flex-1 min-h-[120px]">
                  <FileInput
                    id="genome-path"
                    placeholder="/path/to/genome.gb"
                    value={newOrganism.genomePath || ""}
                    onChange={(value) => setNewOrganism((prev) => ({ ...prev, genomePath: value }))}
                    accept=".gb,.gbf,.gbff,.gbk"
                    fileType="GenBank"
                    onOrganismNameSuggestion={(suggestedName) => {
                      // Only update if the name field is empty
                      if (!newOrganism.name?.trim()) {
                        setNewOrganism((prev) => ({ ...prev, name: suggestedName }))
                      }
                    }}
                  />
                </div>
                
              </div>
              <div className="space-y-2">
                <div/>
                <div>
                   <Label htmlFor="priority">Priority Score (1-100)</Label>
                </div>
                <div className="min-h-[47px] flex items-end">
                  <Input
                    id="priority"
                    type="number"
                    min="1"
                    max="100"
                    accept="1-100"
                    placeholder={`Default: ${calculateDefaultPriority()}`}
                    value={newOrganism.priority || ""}
                    onChange={(e) =>
                      setNewOrganism((prev) => ({
                        ...prev,
                        priority: e.target.value ? Number.parseInt(e.target.value) : undefined,
                      }))
                    }
                  />
                </div>
              </div>
            </div>
            <div className="border rounded-lg p-4 space-y-4">
              <div>
                <h3 className="text-base font-semibold">Expression Data (Optional)</h3>
                <p className="text-sm text-gray-600">
                  {newOrganism.expressionDataFormat === 'json' ?
                    (newOrganism.expressionDataType === 'protein_abundance' ?
                    'For protein abundance, we expect each json record to contain "gene" and "protein_abundance" fields.' :
                    'For mRNA levels, we expect each json record to contain "gene" and "mRNA_level" fields.') : 
                    (newOrganism.expressionDataType === 'protein_abundance' ?
                      'For protein abundance, we expect the csv to contain "gene" and "protein_abundance" columns.' :
                      'For mRNA levels, we expect the csv to contain "gene" "mRNA_level" columns.')
                  }
                </p>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-2 items-start">
                {/* Left side */}
                <div className="space-y-4">
                  <div className="space-y-2">
                    <Label htmlFor="expression-format">File Format</Label>
                    <Select
                      value={newOrganism.expressionDataFormat}
                      onValueChange={(value: "csv" | "json") =>
                        setNewOrganism((prev) => ({ ...prev, expressionDataFormat: value }))
                      }
                    >
                      <SelectTrigger>
                        <SelectValue placeholder="Select format" />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="csv">CSV</SelectItem>
                        <SelectItem value="json">JSON</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>
                  <div className="space-y-2">
                    <Label htmlFor="expression-type">Data Type</Label>
                    <Select
                      value={newOrganism.expressionDataType}
                      onValueChange={(value: "protein_abundance" | "mrna_levels") =>
                        setNewOrganism((prev) => ({ ...prev, expressionDataType: value }))
                      }
                    >
                      <SelectTrigger>
                        <SelectValue placeholder="Select type" />
                      </SelectTrigger>
                      <SelectContent>
                        <SelectItem value="protein_abundance">Protein Abundance</SelectItem>
                        <SelectItem value="mrna_levels">mRNA Levels</SelectItem>
                      </SelectContent>
                    </Select>
                  </div>
                </div>

                {/* Right side */}
                <div className="space-y-2">
                  <div className="flex items-center justify-between">
                    <Label>Expression File</Label>
                    <Button
                      variant="ghost"
                      size="sm"
                      onClick={downloadExpressionSample}
                      className="text-blue-600 hover:text-blue-800"
                    >
                      <Download className="w-4 h-4 mr-1" />
                      Sample
                    </Button>
                  </div>
                  <FileInput
                    id="expression-path"
                    placeholder={`/path/to/expression_data.${newOrganism.expressionDataFormat}`}
                    value={newOrganism.expressionDataPath || ""}
                    onChange={(value) => setNewOrganism((prev) => ({ ...prev, expressionDataPath: value }))}
                    accept={`.${newOrganism.expressionDataFormat}`}
                    fileType={newOrganism.expressionDataFormat?.toUpperCase() || "DATA"}
                  />
                </div>
              </div>
            </div>
          </div>
          <Button
            onClick={handleAddOrganism}
            disabled={!newOrganism.genomePath?.trim() || !newOrganism.name?.trim()}
            className="w-full"
          >
            <Plus className="w-4 h-4 mr-2" />
            Add Organism
          </Button>
        </CardContent>
      </Card>
    </div>
  )
}

interface OrganismCardProps {
  organism: Organism
  onUpdate: (organism: Organism) => void
  onRemove: () => void
}

function OrganismCard({ organism, onUpdate, onRemove }: OrganismCardProps) {
  const getDisplayPath = (fullPath: string) => {
    if (!fullPath) return ""
    return fullPath.split("/").pop() || fullPath
  }

  return (
    <Card>
      <CardContent className="pt-6">
        <div className="flex items-start justify-between mb-4">
          <div className="flex-1">
            <div className="flex items-center gap-2 mb-2">
              <Badge variant="secondary" className="font-medium">
                {organism.name}
              </Badge>
              <Badge variant="outline">Priority: {organism.priority}</Badge>
              <Button variant="ghost" size="sm" onClick={onRemove} className="text-red-500 hover:text-red-700">
                <Trash2 className="w-4 h-4" />
              </Button>
            </div>
          </div>
        </div>

        <div className="space-y-4">
          {/* Organism Name */}
          <div className="space-y-2">
            <Label>Organism Name</Label>
            <Input
              placeholder="e.g., Escherichia-coli"
              value={organism.name}
              onChange={(e) => onUpdate({ ...organism, name: e.target.value })}
            />
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="space-y-2">
              <Label>GenBank Genome File (.gb/.gbf/.gbff)</Label>
              <FileInput
                placeholder="/path/to/genome.gb"
                value={organism.genomePath}
                onChange={(value) => onUpdate({ ...organism, genomePath: value })}
                accept=".gb,.gbf,.gbk,.gbff"
                fileType="GenBank"
                displayValue={getDisplayPath(organism.genomePath)}
                onOrganismNameSuggestion={(suggestedName) => {
                  // Update organism name with suggested name
                  onUpdate({ ...organism, name: suggestedName })
                }}
              />
            </div>

            <div className="space-y-2">
              <Label>Priority Score (1-100)</Label>
              <Input
                type="number"
                min="1"
                max="100"
                value={organism.priority}
                onChange={(e) =>
                  onUpdate({
                    ...organism,
                    priority: Number.parseInt(e.target.value) || 50,
                  })
                }
              />
            </div>
          </div>

          <div className="border rounded-lg p-4 space-y-4">
            <div>
              <h3 className="text-base font-semibold">Expression Data (Optional)</h3>
              <p className="text-sm text-gray-600">
                {(organism.expressionDataFormat || 'csv') === 'json' ?
                  ((organism.expressionDataType || 'protein_abundance') === 'protein_abundance' ?
                 'For protein abundance, we expect each json record to contain "gene" and "protein_abundance" fields.' :
                    'For mRNA levels, we expect each json record to contain "gene" and "mRNA_level" fields.') : 
                    ((organism.expressionDataType || 'protein_abundance') === 'protein_abundance' ?
                      'For protein abundance, we expect the csv to contain "gene" and "protein_abundance" columns.' :
                      'For mRNA levels, we expect the csv to contain "gene" "mRNA_level" columns.')
                  }
              </p>
            </div>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-2 items-start">
              {/* Left side */}
              <div className="space-y-4">
                <div className="space-y-2">
                  <Label htmlFor="expression-format">File Format</Label>
                  <Select
                    value={organism.expressionDataFormat || 'csv'}
                    onValueChange={(value: "csv" | "json") =>
                      onUpdate({ ...organism, expressionDataFormat: value })
                    }
                  >
                    <SelectTrigger>
                      <SelectValue placeholder="Select format" />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="csv">CSV</SelectItem>
                      <SelectItem value="json">JSON</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
                <div className="space-y-2">
                  <Label htmlFor="expression-type">Data Type</Label>
                  <Select
                    value={organism.expressionDataType || 'protein_abundance'}
                    onValueChange={(value: "protein_abundance" | "mrna_levels") =>
                      onUpdate({ ...organism, expressionDataType: value })
                    }
                  >
                    <SelectTrigger>
                      <SelectValue placeholder="Select type" />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="protein_abundance">Protein Abundance</SelectItem>
                      <SelectItem value="mrna_levels">mRNA Levels</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
              </div>

              {/* Right side */}
              <div className="space-y-2">
                <div className="flex items-center justify-between">
                  <Label>Expression File</Label>
                  <Button
                    variant="ghost"
                    size="sm"
                    onClick={() => downloadExpressionSample(newOrganism.expressionDataFormat, newOrganism.expressionDataType)}
                    className="text-blue-600 hover:text-blue-800"
                  >
                    <Download className="w-4 h-4 mr-1" />
                    Sample
                  </Button>
                </div>
                <FileInput
                  placeholder={`/path/to/expression_data.${organism.expressionDataFormat || 'csv'}`}
                  value={organism.expressionDataPath || ""}
                  onChange={(value) =>
                    onUpdate({
                      ...organism,
                      expressionDataPath: value || undefined,
                    })
                  }
                  accept={`.${organism.expressionDataFormat || 'csv'}`}
                  fileType={(organism.expressionDataFormat || "DATA").toUpperCase()}
                  displayValue={organism.expressionDataPath ? getDisplayPath(organism.expressionDataPath) : ""}
                />
              </div>
            </div>
          </div>
        </div>
      </CardContent>
    </Card>
  )
}

interface FileInputProps {
  id?: string
  placeholder: string
  value: string
  onChange: (value: string) => void
  accept: string
  fileType: string
  displayValue?: string
  onOrganismNameSuggestion?: (name: string) => void 
}

function FileInput({ 
  id, 
  placeholder, 
  value, 
  onChange, 
  accept, 
  fileType, 
  displayValue,
  onOrganismNameSuggestion
 }: FileInputProps) {
  const [isUploading, setIsUploading] = useState(false)
  const [uploadProgress, setUploadProgress] = useState(0)
  const [uploadedFile, setUploadedFile] = useState<File | null>(null)
  const [validationError, setValidationError] = useState("")
  const fileInputRef = useRef<HTMLInputElement>(null)

  const handleFileUpload = async (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0]
    if (!file) return

    // Validate file type
    const fileExtension = file.name.toLowerCase().split(".").pop()
    const acceptedExtensions = accept.split(",").map((ext) => ext.trim().replace(".", ""))

    if (!acceptedExtensions.includes(fileExtension || "")) {
      setValidationError(`Please upload a ${fileType} file (${accept})`)
      return
    }

    // Check file size (50MB limit for GenBank files, 10MB for CSV)
    const maxSize = fileType === "GenBank" ? 50 * 1024 * 1024 : 10 * 1024 * 1024
    if (file.size > maxSize) {
      setValidationError(`File size must be less than ${fileType === "GenBank" ? "50MB" : "10MB"}`)
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

      // For GenBank files, do basic validation
      if (fileType === "GenBank") {
        const text = await file.text()
        if (!text.includes("LOCUS") || !text.includes("ORIGIN")) {
          setValidationError("Invalid GenBank format. File must contain LOCUS and ORIGIN sections.")
          setIsUploading(false)
          setUploadProgress(0)
          return
        }
      }

      // For CSV files, validate basic structure
      if (fileType === "CSV") {
        const text = await file.text()
        const lines = text.split("\n")
        if (lines.length < 2) {
          setValidationError("CSV file must contain at least a header row and one data row.")
          setIsUploading(false)
          setUploadProgress(0)
          return
        }
      }

      // Complete upload
      setUploadProgress(100)
      setUploadedFile(file)
      // For web browsers, we can't get the full path due to security restrictions
      // But we can simulate it or use the webkitRelativePath if available
      const fullPath = (file as any).path || file.webkitRelativePath || `/path/to/${file.name}`
      onChange(fullPath)

      if (fileType === "GenBank" && id === "genome-path") {
        const filename = file.name
        const nameWithoutExtension = filename.replace(/\.(gb|gbff|gbk)$/i, "")
        // Trigger organism name update through a callback
        if (onOrganismNameSuggestion) {
          onOrganismNameSuggestion(nameWithoutExtension)
        }
      } 
      // onChange(file.name) // TODO - Store just the filename, in real app this would be the uploaded file path

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
    setUploadedFile(null)
    setValidationError("")
    onChange("")
    if (fileInputRef.current) {
      fileInputRef.current.value = ""
    }
  }

  const handleBrowseClick = () => {
    fileInputRef.current?.click()
  }

  return (
    <div className="space-y-2">
      <div className="flex gap-2">
        <Input
          id={id}
          placeholder={placeholder}
          value={displayValue || value} // Use displayValue if provided
          onChange={(e) => onChange(e.target.value)}
          className="flex-1"
        />
        <Button
          type="button"
          onClick={handleBrowseClick}
          variant="outline"
          size="sm"
          className="px-3 bg-transparent"
          disabled={isUploading}
        >
          <Upload className="w-4 h-4" />
        </Button>
      </div>

      <input ref={fileInputRef} type="file" accept={accept} onChange={handleFileUpload} className="hidden" />

      {/* Upload Progress */}
      {isUploading && (
        <div className="space-y-2">
          <div className="flex justify-between text-sm">
            <span>Uploading {fileType} file...</span>
            <span>{uploadProgress}%</span>
          </div>
          <Progress value={uploadProgress} className="w-full h-2" />
        </div>
      )}

      {/* Uploaded File Display */}
      {uploadedFile && !isUploading && (
        <div className="flex items-center gap-2 p-2 bg-green-50 border border-green-200 rounded-md">
          <FileText className="w-4 h-4 text-green-600" />
          <span className="text-sm text-green-800 flex-1">{uploadedFile.name}</span>
          <span className="text-xs text-green-600">{(uploadedFile.size / 1024).toFixed(1)} KB</span>
          <Button
            variant="ghost"
            size="sm"
            onClick={removeFile}
            className="text-red-500 hover:text-red-700 h-6 w-6 p-0"
          >
            <X className="w-3 h-3" />
          </Button>
        </div>
      )}

      {/* Validation Error */}
      {validationError && (
        <div className="text-red-600 text-sm bg-red-50 p-2 rounded-md border border-red-200">{validationError}</div>
      )}

      {/* File Type Info */}
      <div className="text-xs text-gray-500">
        {fileType === "GenBank" ? (
          <span>Accepted formats: .gb, .gbf, .gbff, .gbk (Max: 50MB)</span>
        ) : (
          <span>Accepted format: .csv (Max: 10MB)</span>
        )}
      </div>
    </div>
  )
}
