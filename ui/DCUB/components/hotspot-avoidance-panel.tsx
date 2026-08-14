"use client"

import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Alert, AlertDescription } from "@/components/ui/alert"
import { CheckCircle, ShieldAlert, AlertTriangle } from "lucide-react"
import type { HotspotAvoidanceResult } from "@/lib/types"

interface HotspotAvoidancePanelProps {
  result: HotspotAvoidanceResult
}

const LINE_LENGTH = 60

/**
 * Split both sequences into aligned fixed-width lines, marking which indexes
 * differ. The two sequences are always the same length - the backend preserves
 * translation and locks everything outside a hotspot - so an index-wise diff is
 * exact and no alignment algorithm is needed.
 */
function buildDiffLines(before: string, after: string) {
  const lines: { start: number; before: string; after: string; changed: boolean[] }[] = []
  for (let start = 0; start < before.length; start += LINE_LENGTH) {
    const beforeSlice = before.slice(start, start + LINE_LENGTH)
    const afterSlice = after.slice(start, start + LINE_LENGTH)
    lines.push({
      start,
      before: beforeSlice,
      after: afterSlice,
      changed: Array.from(beforeSlice, (character, index) => character !== afterSlice[index]),
    })
  }
  return lines
}

function DiffRow({
  label,
  sequence,
  changed,
  changedClassName,
}: {
  label: string
  sequence: string
  changed: boolean[]
  changedClassName: string
}) {
  return (
    <div className="flex gap-3">
      <span className="w-14 shrink-0 select-none text-xs text-gray-400">{label}</span>
      <span>
        {Array.from(sequence, (character, index) => (
          <span key={index} className={changed[index] ? changedClassName : undefined}>
            {character}
          </span>
        ))}
      </span>
    </div>
  )
}

export function HotspotAvoidancePanel({ result }: HotspotAvoidancePanelProps) {
  const { recombination, slippage, motifs } = result.detected_sites
  const diffLines = buildDiffLines(result.sequence_before, result.sequence_after)

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <ShieldAlert className="w-5 h-5 text-blue-600" />
          Hypermutable Site Avoidance
        </CardTitle>
        <CardDescription>
          Detected hypermutable sites were edited away. Codon choices outside those sites were locked.
        </CardDescription>
      </CardHeader>

      <CardContent className="space-y-4">
        <div className="flex flex-wrap gap-2">
          <Badge variant="secondary">Recombination: {recombination}</Badge>
          <Badge variant="secondary">Slippage: {slippage}</Badge>
          <Badge variant="secondary">Methylation motifs: {motifs}</Badge>
          <Badge variant={result.num_edits > 0 ? "default" : "secondary"}>
            {result.num_edits} nucleotide {result.num_edits === 1 ? "edit" : "edits"}
          </Badge>
        </div>

        {result.num_edits === 0 ? (
          <div className="flex items-center gap-2 rounded-lg border border-green-200 bg-green-50 p-4 text-sm text-green-800">
            <CheckCircle className="w-4 h-4 shrink-0" />
            No hypermutable sites needed to be edited in the optimized sequence.
          </div>
        ) : (
          <div className="overflow-x-auto rounded-lg border bg-gray-50 p-4">
            <pre className="text-sm font-mono text-gray-800">
              {diffLines.map((line) => (
                <div key={line.start} className="mb-3">
                  <DiffRow
                    label="before"
                    sequence={line.before}
                    changed={line.changed}
                    changedClassName="bg-red-100 text-red-800 line-through"
                  />
                  <DiffRow
                    label="after"
                    sequence={line.after}
                    changed={line.changed}
                    changedClassName="bg-green-100 text-green-800 font-semibold"
                  />
                </div>
              ))}
            </pre>
          </div>
        )}

        {result.warnings.length > 0 && (
          <Alert variant="destructive">
            <AlertTriangle className="w-4 h-4" />
            <AlertDescription>
              <p className="mb-2 font-medium">
                {result.warnings.length} site{result.warnings.length === 1 ? "" : "s"} could not be cleared:
              </p>
              <ul className="list-disc space-y-1 pl-4 text-xs">
                {result.warnings.map((warning, index) => (
                  <li key={index}>{warning}</li>
                ))}
              </ul>
            </AlertDescription>
          </Alert>
        )}
      </CardContent>
    </Card>
  )
}
