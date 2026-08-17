"use client"

import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Badge } from "@/components/ui/badge"
import { Alert, AlertDescription } from "@/components/ui/alert"
import { CheckCircle, ShieldAlert, AlertTriangle } from "lucide-react"
import type { HotspotAvoidanceResult, HotspotKind, HotspotRegion } from "@/lib/types"

interface HotspotAvoidancePanelProps {
  result: HotspotAvoidanceResult
}

const LINE_LENGTH = 60

// Precedence for colouring a nucleotide covered by more than one window, so the
// rendering is stable no matter what order the backend listed the regions in.
const KIND_ORDER: HotspotKind[] = ["recombination", "slippage", "motifs"]

const KIND_LABELS: Record<HotspotKind, string> = {
  recombination: "Recombination",
  slippage: "Slippage",
  motifs: "Methylation motif",
}

const KIND_CLASSES: Record<HotspotKind, string> = {
  recombination: "bg-purple-100 text-purple-900",
  slippage: "bg-amber-100 text-amber-900",
  motifs: "bg-sky-100 text-sky-900",
}

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

/**
 * Map each nucleotide of the pre-repair sequence to the hotspot kind covering
 * it, or undefined. Region ends are EXCLUSIVE, matching the backend, so the
 * loop stops one short of `end`.
 */
function buildHotspotIndex(length: number, regions: HotspotRegion[]): (HotspotKind | undefined)[] {
  const marks = new Array<HotspotKind | undefined>(length).fill(undefined)
  for (const region of regions) {
    const start = Math.max(0, region.start)
    const end = Math.min(length, region.end)
    for (let index = start; index < end; index += 1) {
      const current = marks[index]
      if (current === undefined || KIND_ORDER.indexOf(region.kind) < KIND_ORDER.indexOf(current)) {
        marks[index] = region.kind
      }
    }
  }
  return marks
}

function SequenceRow({
  label,
  sequence,
  className,
}: {
  label: string
  sequence: React.ReactNode
  className?: string
}) {
  return (
    <div className={`flex gap-3 ${className ?? ""}`}>
      <span className="w-14 shrink-0 select-none text-xs text-gray-400">{label}</span>
      <span>{sequence}</span>
    </div>
  )
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
    <SequenceRow
      label={label}
      sequence={Array.from(sequence, (character, index) => (
        <span key={index} className={changed[index] ? changedClassName : undefined}>
          {character}
        </span>
      ))}
    />
  )
}

export function HotspotAvoidancePanel({ result }: HotspotAvoidancePanelProps) {
  const { recombination, slippage, motifs } = result.detected_sites
  const regions = result.detected_regions ?? []
  const diffLines = buildDiffLines(result.sequence_before, result.sequence_after)

  const hotspotMarks = buildHotspotIndex(result.sequence_before.length, regions)
  const presentKinds = KIND_ORDER.filter((kind) => regions.some((region) => region.kind === kind))

  const detectedLines: { start: number; sequence: string }[] = []
  for (let start = 0; start < result.sequence_before.length; start += LINE_LENGTH) {
    detectedLines.push({ start, sequence: result.sequence_before.slice(start, start + LINE_LENGTH) })
  }

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

        {regions.length > 0 && (
          <div className="space-y-2">
            <div className="flex flex-wrap items-center gap-x-4 gap-y-1">
              <p className="text-sm font-medium text-gray-800">
                Detected sites, shown on the sequence before repair
              </p>
              {presentKinds.map((kind) => (
                <span key={kind} className="flex items-center gap-1.5 text-xs text-gray-600">
                  <span className={`inline-block h-3 w-3 rounded-sm ${KIND_CLASSES[kind]}`} />
                  {KIND_LABELS[kind]}
                </span>
              ))}
            </div>
            <div className="overflow-x-auto rounded-lg border bg-gray-50 p-4">
              <pre className="text-sm font-mono text-gray-800">
                {detectedLines.map((line) => (
                  <SequenceRow
                    key={line.start}
                    label={String(line.start + 1)}
                    className="mb-1"
                    sequence={Array.from(line.sequence, (character, offset) => {
                      const kind = hotspotMarks[line.start + offset]
                      return (
                        <span key={offset} className={kind ? KIND_CLASSES[kind] : undefined}>
                          {character}
                        </span>
                      )
                    })}
                  />
                ))}
              </pre>
            </div>
            <ul className="space-y-0.5 text-xs text-gray-600">
              {regions.map((region, index) => (
                <li key={`${region.kind}-${region.start}-${region.end}-${index}`}>
                  <span className={`rounded-sm px-1 ${KIND_CLASSES[region.kind]}`}>
                    {KIND_LABELS[region.kind]}
                  </span>{" "}
                  {/* Displayed 1-indexed and inclusive - how a biologist reads a
                      position - from the backend's 0-indexed, exclusive-end pair. */}
                  {region.start + 1}&ndash;{region.end} ({region.end - region.start} nt)
                </li>
              ))}
            </ul>
          </div>
        )}

        {result.num_edits === 0 ? (
          <div
            className={`flex items-center gap-2 rounded-lg border p-4 text-sm ${
              regions.length > 0
                ? "border-amber-200 bg-amber-50 text-amber-900"
                : "border-green-200 bg-green-50 text-green-800"
            }`}
          >
            {regions.length > 0 ? (
              <>
                <AlertTriangle className="w-4 h-4 shrink-0" />
                Sites were detected but none could be edited at codon resolution - the sequence is
                unchanged. See the notes below.
              </>
            ) : (
              <>
                <CheckCircle className="w-4 h-4 shrink-0" />
                No hypermutable sites were detected in the optimized sequence.
              </>
            )}
          </div>
        ) : (
          <div className="space-y-2">
            <p className="text-sm font-medium text-gray-800">Edits applied</p>
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
