"use client"

import { useState } from "react"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Input } from "@/components/ui/input"
import { Label } from "@/components/ui/label"
import { Badge } from "@/components/ui/badge"
import { Plus, Trash2, Download } from "lucide-react"
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

  const downloadSampleCSV = () => {
    const sampleCSV = `gene_id,expression_level,tpm,fpkm
gene_001,1250.5,45.2,32.1
gene_002,890.3,28.7,21.4
gene_003,2100.8,67.9,48.3
gene_004,456.2,15.8,11.2
gene_005,1780.4,58.1,41.7`

    const blob = new Blob([sampleCSV], { type: "text/csv" })
    const url = URL.createObjectURL(blob)
    const a = document.createElement("a")
    a.href = url
    a.download = "sample_expression_data.csv"
    document.body.appendChild(a)
    a.click()
    document.body.removeChild(a)
    URL.revokeObjectURL(url)
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
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="space-y-2">
              <Label htmlFor="genome-path">GenBank Genome File Path (.gb/.gbf) *</Label>
              <Input
                id="genome-path"
                placeholder="/path/to/genome.gb"
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
            <div className="flex items-center justify-between">
              <Label htmlFor="expression-path">Expression Data CSV Path (Optional)</Label>
              <Button
                variant="ghost"
                size="sm"
                onClick={downloadSampleCSV}
                className="text-blue-600 hover:text-blue-800"
              >
                <Download className="w-4 h-4 mr-1" />
                Sample CSV
              </Button>
            </div>
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
}

function OrganismCard({ organism, onUpdate, onRemove }: OrganismCardProps) {
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
            <Label>GenBank Genome File Path (.gb/.gbf)</Label>
            <Input
              value={organism.genomePath}
              onChange={(e) => onUpdate({ ...organism, genomePath: e.target.value })}
              placeholder="/path/to/genome.gb"
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

        <div className="mt-4 space-y-2">
          <Label>Expression Data CSV (Optional)</Label>
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
        </div>
      </CardContent>
    </Card>
  )
}
