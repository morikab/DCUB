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
import { Plus, Trash2, Download, Upload, FileText, X, ChevronDown, ChevronUp } from "lucide-react"
import { useOptimizationStore } from "@/lib/store"
import type { Organism } from "@/lib/types"
import { ErrorDialog } from "@/components/error-dialog"
import type { ElectronFile, ElectronFileOperationResult } from "@/lib/electron-utils"

const downloadExpressionSample = (
  format: "csv" | "json" = "json",
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
    expressionDataFormat: "json",
    expressionDataType: "protein_abundance",
  })
  const [formResetKey, setFormResetKey] = useState(0)
  const [isExpressionDataOpen, setIsExpressionDataOpen] = useState(false)
  const [errorDialogOpen, setErrorDialogOpen] = useState(false)
  const [errorDialogText, setErrorDialogText] = useState("")
  const nameInputRef = useRef<HTMLInputElement>(null)

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
  const opposingOrganisms = type === "wanted" ? unwantedOrganisms : wantedOrganisms
  const addOrganism = type === "wanted" ? addWantedOrganism : addUnwantedOrganism
  const removeOrganism = type === "wanted" ? removeWantedOrganism : removeUnwantedOrganism
  const updateOrganism = type === "wanted" ? updateWantedOrganism : updateUnwantedOrganism

  const calculateDefaultPriority = () => {
    if (organisms.length === 0) return 50
    const validPriorities = organisms.filter((org) => org.priority !== undefined).map((org) => org.priority!)
    if (validPriorities.length === 0) return 50
    return Math.round(validPriorities.reduce((sum, p) => sum + p, 0) / validPriorities.length)
  }

  const trimmedNewName = newOrganism.name?.trim().toLowerCase()
  const isDuplicateNewName = !!trimmedNewName && organisms.some((org) => org.name.toLowerCase() === trimmedNewName)

  const trimmedNewPath = newOrganism.genomePath?.trim() || ""
  const isDuplicateInSameList = !!trimmedNewPath && organisms.some((org) => org.genomePath.trim() === trimmedNewPath)
  const isDuplicateInOpposingList = !!trimmedNewPath && opposingOrganisms.some((org) => org.genomePath.trim() === trimmedNewPath)
  const isDuplicateGenomePath = isDuplicateInSameList || isDuplicateInOpposingList

  const handleGenomePathChange = (value: string) => {
    setNewOrganism((prev) => ({ ...prev, genomePath: value }))
    const trimmed = value.trim()
    if (trimmed && opposingOrganisms.some((org) => org.genomePath.trim() === trimmed)) {
      setErrorDialogText(
        `This genome file is already used by an organism in the ${type === "wanted" ? "unwanted" : "wanted"} organisms list. Each organism must have a unique genome file.`
      )
      setErrorDialogOpen(true)
    }
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
      expressionDataFormat: "json",
      expressionDataType: "protein_abundance",
    })
    setFormResetKey((key) => key + 1)
    setIsExpressionDataOpen(false)
    nameInputRef.current?.focus()
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

      {/* Weight distribution summary */}
      {organisms.length >= 2 && (() => {
        const totalPriority = organisms.reduce((sum, org) => sum + org.priority, 0)
        const segmentColors = ["bg-blue-400", "bg-emerald-400", "bg-violet-400", "bg-orange-400", "bg-pink-400", "bg-teal-400"]
        return (
          <div className="space-y-1">
            <div className="flex h-2 rounded-full overflow-hidden gap-px">
              {organisms.map((org, i) => {
                const w = totalPriority > 0 ? org.priority / totalPriority : 0
                return (
                  <div
                    key={org.id}
                    className={segmentColors[i % segmentColors.length]}
                    style={{ width: `${w * 100}%` }}
                    title={`${org.name}: ${w.toFixed(2)}`}
                  />
                )
              })}
            </div>
            <p className="text-xs text-gray-500">
              {organisms.map((org) => {
                const w = totalPriority > 0 ? org.priority / totalPriority : 0
                return `${org.name} ${w.toFixed(2)}`
              }).join(" · ")}
            </p>
          </div>
        )
      })()}

      {/* Existing Organisms */}
      <div className="space-y-3">
        {(() => {
          const totalPriority = organisms.reduce((sum, org) => sum + org.priority, 0)
          return organisms.map((organism) => (
            <OrganismCard
              key={organism.id}
              organism={organism}
              effectiveWeight={totalPriority > 0 ? organism.priority / totalPriority : 0}
              otherOrganismNames={organisms
                .filter((other) => other.id !== organism.id)
                .map((other) => other.name.toLowerCase())}
              otherGenomePaths={[
                ...organisms.filter((other) => other.id !== organism.id).map((other) => other.genomePath.trim()),
                ...opposingOrganisms.map((other) => other.genomePath.trim()),
              ].filter(Boolean)}
              onUpdate={(updated) => updateOrganism(organism.id, updated)}
              onRemove={() => removeOrganism(organism.id)}
            />
          ))
        })()}
      </div>

      <ErrorDialog
        open={errorDialogOpen}
        title="Duplicate genome file"
        errorText={errorDialogText}
        onClose={() => setErrorDialogOpen(false)}
      />

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
              ref={nameInputRef}
              placeholder="e.g., Escherichia-coli, Bacillus-subtilis"
              value={newOrganism.name || ""}
              onChange={(e) => setNewOrganism((prev) => ({ ...prev, name: e.target.value }))}
            />
            {isDuplicateNewName && (
              <p className="text-sm text-amber-600">
                An organism named &quot;{newOrganism.name?.trim()}&quot; already exists in this list.
              </p>
            )}
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
                    key={`genome-${formResetKey}`}
                    id="genome-path"
                    placeholder="/path/to/genome.db"
                    value={newOrganism.genomePath || ""}
                    onChange={handleGenomePathChange}
                    accept=".gb,.gbf,.gbff,.gbk"
                    fileType="GenBank"
                    onOrganismNameSuggestion={(suggestedName) => {
                      if (!newOrganism.name?.trim()) {
                        setNewOrganism((prev) => ({ ...prev, name: suggestedName }))
                      }
                    }}
                  />
                </div>
                {isDuplicateInSameList && (
                  <p className="text-sm text-amber-600">
                    This genome file is already used by another organism in this list.
                  </p>
                )}
                {isDuplicateInOpposingList && !isDuplicateInSameList && (
                  <p className="text-sm text-amber-600">
                    This genome file is already used by an organism in the {type === "wanted" ? "unwanted" : "wanted"} organisms list.
                  </p>
                )}
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
              <button
                type="button"
                onClick={() => setIsExpressionDataOpen((open) => !open)}
                className="flex items-center justify-between w-full text-left"
              >
                <h3 className="text-base font-semibold">Expression Data (Optional)</h3>
                {isExpressionDataOpen ? (
                  <ChevronUp className="w-4 h-4 text-gray-500" />
                ) : (
                  <ChevronDown className="w-4 h-4 text-gray-500" />
                )}
              </button>

              {isExpressionDataOpen && (
                <>
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
                          onClick={() => downloadExpressionSample(newOrganism.expressionDataFormat, newOrganism.expressionDataType)}
                          className="text-blue-600 hover:text-blue-800"
                        >
                          <Download className="w-4 h-4 mr-1" />
                          Sample
                        </Button>
                      </div>
                      <FileInput
                        key={`expression-${formResetKey}`}
                        id="expression-path"
                        placeholder={`/path/to/expression_data.${newOrganism.expressionDataFormat}`}
                        value={newOrganism.expressionDataPath || ""}
                        onChange={(value) => setNewOrganism((prev) => ({ ...prev, expressionDataPath: value }))}
                        accept={`.${newOrganism.expressionDataFormat}`}
                        fileType={newOrganism.expressionDataFormat?.toUpperCase() || "DATA"}
                      />
                    </div>
                  </div>
                </>
              )}
            </div>
          </div>
          <Button
            onClick={handleAddOrganism}
            disabled={!newOrganism.genomePath?.trim() || !newOrganism.name?.trim() || isDuplicateNewName || isDuplicateGenomePath}
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
  effectiveWeight: number
  otherOrganismNames: string[]
  otherGenomePaths: string[]
  onUpdate: (organism: Organism) => void
  onRemove: () => void
}

function OrganismCard({ organism, effectiveWeight, otherOrganismNames, otherGenomePaths, onUpdate, onRemove }: OrganismCardProps) {
  const [isExpanded, setIsExpanded] = useState(false)
  const [isExpressionDataOpen, setIsExpressionDataOpen] = useState(!!organism.expressionDataPath)

  const getDisplayPath = (fullPath: string) => {
    if (!fullPath) return ""
    return fullPath.split("/").pop() || fullPath
  }

  const isDuplicateName = otherOrganismNames.includes(organism.name.trim().toLowerCase())
  const isDuplicateGenomePath = !!organism.genomePath.trim() && otherGenomePaths.includes(organism.genomePath.trim())

  if (!isExpanded) {
    return (
      <Card>
        <CardContent className="py-3">
          <div className="flex items-center justify-between gap-2">
            <button
              type="button"
              onClick={() => setIsExpanded(true)}
              className="flex flex-1 items-center gap-2 min-w-0 text-left"
            >
              <ChevronDown className="w-4 h-4 text-gray-400 flex-shrink-0" />
              <Badge variant="secondary" className="font-medium">
                {organism.name}
              </Badge>
              <Badge variant="outline">Priority: {organism.priority}</Badge>
              <Badge variant="outline" className="text-gray-500">Weight: {effectiveWeight.toFixed(2)}</Badge>
              <span className="text-sm text-gray-500 truncate">{getDisplayPath(organism.genomePath)}</span>
              {organism.expressionDataPath && (
                <Badge variant="outline" className="text-green-700 border-green-300">
                  Expression data
                </Badge>
              )}
            </button>
            <Button variant="ghost" size="sm" onClick={onRemove} className="text-red-500 hover:text-red-700">
              <Trash2 className="w-4 h-4" />
            </Button>
          </div>
        </CardContent>
      </Card>
    )
  }

  return (
    <Card>
      <CardContent className="pt-6">
        <div className="flex items-start justify-between mb-4">
          <button
            type="button"
            onClick={() => setIsExpanded(false)}
            className="flex flex-1 items-center gap-2 text-left"
          >
            <ChevronUp className="w-4 h-4 text-gray-400" />
            <Badge variant="secondary" className="font-medium">
              {organism.name}
            </Badge>
            <Badge variant="outline">Priority: {organism.priority}</Badge>
            <Badge variant="outline" className="text-gray-500">Weight: {effectiveWeight.toFixed(2)}</Badge>
          </button>
          <Button variant="ghost" size="sm" onClick={onRemove} className="text-red-500 hover:text-red-700">
            <Trash2 className="w-4 h-4" />
          </Button>
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
            {isDuplicateName && (
              <p className="text-sm text-amber-600">
                An organism named &quot;{organism.name.trim()}&quot; already exists in this list.
              </p>
            )}
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="space-y-2">
              <Label>GenBank Genome File (.gb/.gbf/.gbff)</Label>
              <FileInput
                placeholder="/path/to/genome.db"
                value={organism.genomePath}
                onChange={(value) => onUpdate({ ...organism, genomePath: value })}
                accept=".gb,.gbf,.gbk,.gbff"
                fileType="GenBank"
                displayValue={getDisplayPath(organism.genomePath)}
                onOrganismNameSuggestion={(suggestedName) => {
                  onUpdate({ ...organism, name: suggestedName })
                }}
              />
              {isDuplicateGenomePath && (
                <p className="text-sm text-amber-600">
                  This genome file is already used by another organism.
                </p>
              )}
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
              <p className="text-xs text-gray-500">Effective weight: {effectiveWeight.toFixed(2)}</p>
            </div>
          </div>

          <div className="border rounded-lg p-4 space-y-4">
            <button
              type="button"
              onClick={() => setIsExpressionDataOpen((open) => !open)}
              className="flex items-center justify-between w-full text-left"
            >
              <h3 className="text-base font-semibold">Expression Data (Optional)</h3>
              {isExpressionDataOpen ? (
                <ChevronUp className="w-4 h-4 text-gray-500" />
              ) : (
                <ChevronDown className="w-4 h-4 text-gray-500" />
              )}
            </button>

            {isExpressionDataOpen && (
              <>
                <p className="text-sm text-gray-600">
                  {(organism.expressionDataFormat || 'json') === 'json' ?
                    ((organism.expressionDataType || 'protein_abundance') === 'protein_abundance' ?
                   'For protein abundance, we expect each json record to contain "gene" and "protein_abundance" fields.' :
                      'For mRNA levels, we expect each json record to contain "gene" and "mRNA_level" fields.') :
                      ((organism.expressionDataType || 'protein_abundance') === 'protein_abundance' ?
                        'For protein abundance, we expect the csv to contain "gene" and "protein_abundance" columns.' :
                        'For mRNA levels, we expect the csv to contain "gene" "mRNA_level" columns.')
                    }
                </p>

                <div className="grid grid-cols-1 md:grid-cols-2 gap-x-4 gap-y-2 items-start">
                  {/* Left side */}
                  <div className="space-y-4">
                    <div className="space-y-2">
                      <Label htmlFor="expression-format">File Format</Label>
                      <Select
                        value={organism.expressionDataFormat || 'json'}
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
                        onClick={() => downloadExpressionSample(organism.expressionDataFormat || 'json', organism.expressionDataType || 'protein_abundance')}
                        className="text-blue-600 hover:text-blue-800"
                      >
                        <Download className="w-4 h-4 mr-1" />
                        Sample
                      </Button>
                    </div>
                    <FileInput
                      placeholder={`/path/to/expression_data.${organism.expressionDataFormat || 'json'}`}
                      value={organism.expressionDataPath || ""}
                      onChange={(value) =>
                        onUpdate({
                          ...organism,
                          expressionDataPath: value || undefined,
                        })
                      }
                      accept={`.${organism.expressionDataFormat || 'json'}`}
                      fileType={(organism.expressionDataFormat || "DATA").toUpperCase()}
                      displayValue={organism.expressionDataPath ? getDisplayPath(organism.expressionDataPath) : ""}
                    />
                  </div>
                </div>
              </>
            )}
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
  const [isDragOver, setIsDragOver] = useState(false)
  const fileInputRef = useRef<HTMLInputElement>(null)

  const processFile = async (file: File) => {
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

      // Electron exposes file.path; browsers don't. In dev (browser) mode, upload the
      // file to the backend so it gets a real server-side path for /run-modules.
      let fullPath: string
      const electronPath = (file as ElectronFile).path
      if (electronPath) {
        fullPath = electronPath
      } else {
        const formData = new FormData()
        formData.append("file", file)
        const uploadRes = await fetch("http://localhost:8000/upload-file", { method: "POST", body: formData })
        if (!uploadRes.ok) throw new Error("File upload to dev server failed")
        const { file_path } = await uploadRes.json()
        fullPath = file_path
      }
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
    } catch {
      setValidationError("Error reading file. Please try again.")
      setIsUploading(false)
      setUploadProgress(0)
    }
  }

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0]
    if (file) processFile(file)
  }

  const removeFile = () => {
    setUploadedFile(null)
    setValidationError("")
    onChange("")
    if (fileInputRef.current) {
      fileInputRef.current.value = ""
    }
  }

  const handleBrowseClick = async () => {
    const electronAPI = typeof window !== "undefined" ? window.electronAPI : null
    if (electronAPI?.isElectron) {
      const extensions = accept.split(",").map((e: string) => e.trim().replace(".", ""))
      let result: ElectronFileOperationResult | undefined
      try {
        result = await electronAPI.fileOperations?.("select-file", {
          filters: [{ name: `${fileType} Files`, extensions }],
        })
      } catch {
        setValidationError("Error opening file dialog. Please try again.")
        return
      }
      // Cancelled, unknown-operation and thrown-error replies all come back
      // as success:false; only the "picked a file" reply carries these three.
      if (!result?.success || !result.fileName || !result.filePath || result.content === undefined) return

      const { filePath, fileName, content } = result

      const fileExtension = fileName.toLowerCase().split(".").pop() || ""
      const acceptedExtensions = accept.split(",").map((e: string) => e.trim().replace(".", ""))
      if (!acceptedExtensions.includes(fileExtension)) {
        setValidationError(`Please upload a ${fileType} file (${accept})`)
        return
      }

      setIsUploading(true)
      setUploadProgress(50)
      setValidationError("")

      if (fileType === "GenBank") {
        if (!content.includes("LOCUS") || !content.includes("ORIGIN")) {
          setValidationError("Invalid GenBank format. File must contain LOCUS and ORIGIN sections.")
          setIsUploading(false)
          setUploadProgress(0)
          return
        }
      }
      if (fileType === "CSV") {
        const lines = (content as string).split("\n")
        if (lines.length < 2) {
          setValidationError("CSV file must contain at least a header row and one data row.")
          setIsUploading(false)
          setUploadProgress(0)
          return
        }
      }

      setUploadProgress(100)
      onChange(filePath)

      if (fileType === "GenBank" && id === "genome-path" && onOrganismNameSuggestion) {
        onOrganismNameSuggestion((fileName as string).replace(/\.(gb|gbff|gbk)$/i, ""))
      }

      setTimeout(() => { setIsUploading(false); setUploadProgress(0) }, 500)
      return
    }

    fileInputRef.current?.click()
  }

  const handleDrop = (event: React.DragEvent<HTMLDivElement>) => {
    event.preventDefault()
    setIsDragOver(false)
    if (isUploading) return
    const file = event.dataTransfer.files?.[0]
    if (file) processFile(file)
  }

  const handleDragOver = (event: React.DragEvent<HTMLDivElement>) => {
    event.preventDefault()
    setIsDragOver(true)
  }

  const handleDragLeave = (event: React.DragEvent<HTMLDivElement>) => {
    event.preventDefault()
    setIsDragOver(false)
  }

  const baseName = (path: string) => path.split("/").pop() || path
  const hasValue = !!value
  const selectedLabel = displayValue || (uploadedFile ? uploadedFile.name : value ? baseName(value) : "")
  const acceptedFormatsLabel =
    fileType === "GenBank"
      ? "Accepted formats: .gb, .gbf, .gbff, .gbk (Max: 50MB)"
      : `Accepted format: ${accept} (Max: 10MB)`

  return (
    <div className="space-y-2">
      {isUploading ? (
        <div className="space-y-2 rounded-md border-2 border-dashed border-gray-300 p-4">
          <div className="flex justify-between text-sm">
            <span>Uploading {fileType} file...</span>
            <span>{uploadProgress}%</span>
          </div>
          <Progress value={uploadProgress} className="w-full h-2" />
        </div>
      ) : (
        <div
          onClick={handleBrowseClick}
          onDrop={handleDrop}
          onDragOver={handleDragOver}
          onDragLeave={handleDragLeave}
          role="button"
          tabIndex={0}
          className={`flex flex-col items-center justify-center gap-1 rounded-md border-2 border-dashed p-4 text-center cursor-pointer transition-colors ${
            isDragOver
              ? "border-blue-500 bg-blue-50"
              : hasValue
                ? "border-green-300 bg-green-50"
                : "border-gray-300 hover:border-gray-400 hover:bg-gray-50"
          }`}
        >
          {hasValue ? (
            <div className="flex w-full items-center gap-2">
              <FileText className="w-5 h-5 text-green-600 flex-shrink-0" />
              <div className="flex-1 min-w-0 text-left">
                <p className="text-sm text-green-800 truncate">{selectedLabel}</p>
                <p className="text-xs text-green-600">
                  {uploadedFile ? `${(uploadedFile.size / 1024).toFixed(1)} KB · ` : ""}
                  Click or drop to replace
                </p>
              </div>
              <Button
                type="button"
                variant="ghost"
                size="sm"
                onClick={(e) => {
                  e.stopPropagation()
                  removeFile()
                }}
                className="h-6 w-6 p-0 text-red-500 hover:text-red-700 flex-shrink-0"
              >
                <X className="w-3 h-3" />
              </Button>
            </div>
          ) : (
            <>
              <Upload className="w-6 h-6 text-gray-400" />
              <p className="text-sm font-medium text-gray-700">Click to upload or drag and drop</p>
              <p className="text-xs text-gray-500">{acceptedFormatsLabel}</p>
            </>
          )}
        </div>
      )}

      <input ref={fileInputRef} type="file" accept={accept} onChange={handleFileUpload} className="hidden" />

      {/* Manual path entry, as a fallback to drag-and-drop / browse */}
      <Input
        id={id}
        placeholder={`Or paste a file path, e.g. ${placeholder}`}
        value={value}
        onChange={(e) => onChange(e.target.value)}
        disabled={isUploading}
        className="text-sm"
      />

      {/* Validation Error */}
      {validationError && (
        <div className="text-red-600 text-sm bg-red-50 p-2 rounded-md border border-red-200">{validationError}</div>
      )}
    </div>
  )
}
