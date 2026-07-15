import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { FlaskConical } from "lucide-react"

interface OrganismDistScore {
  name: string
  is_wanted: boolean
  dist_score: number
}

interface OrganismDistChartProps {
  organismsDistScores: OrganismDistScore[]
}

function getBarColor(isWanted: boolean, distScore: number): string {
  if (isWanted) {
    return distScore >= 0 ? "bg-green-500" : "bg-green-200"
  }
  return distScore < 0 ? "bg-red-500" : "bg-red-200"
}

function OrganismRow({
  organism,
  maxAbs,
}: {
  organism: OrganismDistScore
  maxAbs: number
}) {
  const barWidthPct = maxAbs > 0 ? (Math.abs(organism.dist_score) / maxAbs) * 50 : 0
  const isPositive = organism.dist_score >= 0
  const barColor = getBarColor(organism.is_wanted, organism.dist_score)

  return (
    <div className="flex items-center gap-3 py-1.5">
      {/* Name column */}
      <div className="w-1/4 min-w-0 flex items-center gap-2">
        <span
          className="text-xs font-medium px-1.5 py-0.5 rounded flex-shrink-0"
          style={{
            backgroundColor: organism.is_wanted ? "#dcfce7" : "#fee2e2",
            color: organism.is_wanted ? "#166534" : "#991b1b",
          }}
        >
          {organism.is_wanted ? "W" : "U"}
        </span>
        <span className="text-sm text-gray-800 truncate" title={organism.name}>
          {organism.name}
        </span>
      </div>

      {/* Bar area */}
      <div className="flex-1 relative h-6">
        {/* Center zero line */}
        <div className="absolute left-1/2 top-0 bottom-0 w-px bg-gray-300" />

        {/* Bar */}
        {barWidthPct > 0 && (
          <div
            className={`absolute top-1 bottom-1 rounded-sm ${barColor}`}
            style={
              isPositive
                ? { left: "50%", width: `${barWidthPct}%` }
                : { right: "50%", width: `${barWidthPct}%` }
            }
          />
        )}
      </div>

      {/* Score value */}
      <div className="w-16 text-right">
        <span className="text-sm font-mono text-gray-700">
          {organism.dist_score >= 0 ? "+" : ""}
          {organism.dist_score.toFixed(2)}
        </span>
      </div>
    </div>
  )
}

export function OrganismDistChart({ organismsDistScores }: OrganismDistChartProps) {
  if (!organismsDistScores || organismsDistScores.length === 0) return null

  const wanted = organismsDistScores.filter((o) => o.is_wanted)
  const unwanted = organismsDistScores.filter((o) => !o.is_wanted)
  const maxAbs = Math.max(...organismsDistScores.map((o) => Math.abs(o.dist_score)), 0)

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <FlaskConical className="w-5 h-5" />
          Per-Organism CUB Difference
        </CardTitle>
        <CardDescription>
          z-score of CUB change relative to each organism&apos;s proteome variability — positive = improved
          expression, negative = reduced.
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-4">
        {wanted.length > 0 && (
          <div>
            <p className="text-xs font-semibold text-gray-500 uppercase tracking-wide mb-2">
              Wanted Organisms
            </p>
            <div className="space-y-0.5">
              {wanted.map((org) => (
                <OrganismRow key={org.name} organism={org} maxAbs={maxAbs} />
              ))}
            </div>
          </div>
        )}

        {wanted.length > 0 && unwanted.length > 0 && (
          <div className="border-t border-gray-100" />
        )}

        {unwanted.length > 0 && (
          <div>
            <p className="text-xs font-semibold text-gray-500 uppercase tracking-wide mb-2">
              Unwanted Organisms
            </p>
            <div className="space-y-0.5">
              {unwanted.map((org) => (
                <OrganismRow key={org.name} organism={org} maxAbs={maxAbs} />
              ))}
            </div>
          </div>
        )}
      </CardContent>
    </Card>
  )
}
