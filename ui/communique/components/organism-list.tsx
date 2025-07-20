"use client"

import type React from "react"

import { useState } from "react"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { Badge } from "@/components/ui/badge"
import { Plus, Trash2, Upload } from "lucide-react"
import { useOptimizationStore } from "@/lib/store"
import type { Organism } from "@/lib/types"

interface OrganismListProps {
  type: "wanted" | "unwanted"
}

export function OrganismList({ type }: OrganismListProps) {
  const [newOrganism, setNewOrganism] = useState<Partial<Organism>>({
    genomePath: "",
    priority: undefined,
    expressionDataPath: "",
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
    if (!newOrganism.genomePath?.trim()) return

    const organism: Organism = {
      id: Date.now().toString(),
      genomePath: newOrganism.genomePath.trim(),
      priority: newOrganism.priority || calculateDefaultPriority(),
      expressionDataPath: newOrganism.expressionDataPath?.trim() || undefined,
    }

    addOrganism(organism)
    setNewOrganism({
      genomePath: "",
      priority: undefined,
      expressionDataPath: "",
    })
  }

  const handleFileUpload = async (organismId: string, fileType: "genome" | "expression", file: File) => {
    // In a real application, you would upload the file to a server
    // For now, we'll just store the file path/name
    const filePath = `/uploads/${file.name}`

    const organism = organisms.find((org) => org.id === organismId)
    if (!organism) return

    const updatedOrganism = {
      ...organism,
      [fileType === "genome" ? "genomePath" : "expressionDataPath"]: filePath,
    }

    updateOrganism(organismId, updatedOrganism)
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
            onFileUpload={(fileType, file) => handleFileUpload(organism.id, fileType, file)}
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
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="space-y-2">
              <Label htmlFor="genome-path">Genome File Path *</Label>
              <Input
                id="genome-path"
                placeholder="/path/to/genome.fasta"
                value={newOrganism.genomePath || ""}
                onChange={(e) => setNewOrganism((prev) => ({ ...prev, genomePath: e.target.value }))}
              />
            </div>
            <div className="space-y-2">
              <Label htmlFor="priority">Priority Score (1-100)</Label>
              <Input
                id="priority"
                type="number"
                min="1"
                max="100"
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
          <div className="space-y-2">
            <Label htmlFor="expression-path">Expression Data CSV Path (Optional)</Label>
            <Input
              id="expression-path"
              placeholder="/path/to/expression_data.csv"
              value={newOrganism.expressionDataPath || ""}
              onChange={(e) => setNewOrganism((prev) => ({ ...prev, expressionDataPath: e.target.value }))}
            />
          </div>
          <Button onClick={handleAddOrganism} disabled={!newOrganism.genomePath?.trim()} className="w-full">
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
  onFileUpload: (fileType: "genome" | "expression", file: File) => void
}

function OrganismCard({ organism, onUpdate, onRemove, onFileUpload }: OrganismCardProps) {
  const handleFileChange = (fileType: "genome" | "expression") => (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0]
    if (file) {
      onFileUpload(fileType, file)
    }
  }

  return (
    <Card>
      <CardContent className="pt-6">
        <div className="flex items-start justify-between mb-4">
          <div className="flex-1">
            <div className="flex items-center gap-2 mb-2">
              <Badge variant="outline">Priority: {organism.priority}</Badge>
              <Button variant="ghost" size="sm" onClick={onRemove} className="text-red-500 hover:text-red-700">
                <Trash2 className="w-4 h-4" />
              </Button>
            </div>
          </div>
        </div>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <div className="space-y-2">
            <Label>Genome File Path</Label>
            <div className="flex gap-2">
              <Input
                value={organism.genomePath}
                onChange={(e) => onUpdate({ ...organism, genomePath: e.target.value })}
                placeholder="/path/to/genome.fasta"
              />
              <Button
                variant="outline"
                size="sm"
                onClick={() => document.getElementById(`genome-${organism.id}`)?.click()}
              >
                <Upload className="w-4 h-4" />
              </Button>
              <input
                id={`genome-${organism.id}`}
                type="file"
                accept=".fasta,.fa,.gbk,.gb"
                onChange={handleFileChange("genome")}
                className="hidden"
              />
            </div>
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

        <div className="mt-4 space-y-2">
          <Label>Expression Data CSV (Optional)</Label>
          <div className="flex gap-2">
            <Input
              value={organism.expressionDataPath || ""}
              onChange={(e) =>
                onUpdate({
                  ...organism,
                  expressionDataPath: e.target.value || undefined,
                })
              }
              placeholder="/path/to/expression_data.csv"
            />
            <Button
              variant="outline"
              size="sm"
              onClick={() => document.getElementById(`expression-${organism.id}`)?.click()}
            >
              <Upload className="w-4 h-4" />
            </Button>
            <input
              id={`expression-${organism.id}`}
              type="file"
              accept=".csv"
              onChange={handleFileChange("expression")}
              className="hidden"
            />
          </div>
        </div>
      </CardContent>
    </Card>
  )
}
