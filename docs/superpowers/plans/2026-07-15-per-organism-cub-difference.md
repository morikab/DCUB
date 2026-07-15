# Per-Organism CUB Difference Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Surface per-organism CUB dist_scores (already computed in the backend) through the API and display them as a diverging horizontal bar chart on the results page.

**Architecture:** The backend change adds the organism-level breakdown to `EvaluationModuleResult.summary`, which flows automatically into `final_evaluation` in the run summary and therefore into the API response. The frontend parses this new field, stores it in `OptimizationResult`, and renders it in a new `OrganismDistChart` component — a full-width card inserted between the top result grid and the Optimization Details card.

**Tech Stack:** Python / FastAPI (backend), Next.js 15 + React 19 + TypeScript + Tailwind CSS v4 (frontend). No new dependencies.

## Global Constraints

- Python `~3.9` (no walrus operator, no `match` statements)
- No new npm packages — chart is pure CSS/Tailwind
- Tailwind CSS v4 (utility classes, no `@apply` in new code)
- Color scheme: wanted = green, unwanted = red; muted variant when score direction is wrong
- Organisms rendered in two groups: Wanted first, then Unwanted, with a section label above each group
- `maxAbs` bar scaling is computed across ALL organisms (both groups share the same scale)
- Guard against `maxAbs === 0` (all scores are zero) to avoid division by zero
- No automated test suite — each task ends with a manual verification step

---

### Task 1: Backend — expose per-organism dist_scores in `final_evaluation`

**Files:**
- Modify: `app/modules/evaluation/models.py`
- Modify: `app/modules/evaluation/evaluation.py`

**Interfaces:**
- Produces: `EvaluationModuleResult.organisms_dist_scores: List[Dict]` where each dict has `name: str`, `is_wanted: bool`, `dist_score: float`
- Produces: `final_evaluation.organisms` in the API response (list of the same dicts)

---

- [ ] **Step 1: Add `organisms_dist_scores` field to `EvaluationModuleResult`**

Open `app/modules/evaluation/models.py`. The current file looks like:

```python
import typing
from dataclasses import dataclass
from modules import models


@dataclass
class EvaluationModuleResult:
    sequence: str
    average_distance_score: float
    weakest_link_score: float
    ratio_score: float

    @property
    def summary(self) -> typing.Dict[str, typing.Any]:
        return {
            "final_sequence": self.sequence,
            "average_distance_score": self.average_distance_score,
            "weakest_link_score": self.weakest_link_score,
            "ratio_score": self.ratio_score,
        }

    def get_score(self, score_type: models.EvaluationScore) -> float:
        score_value = score_type.value
        if score_value == models.EvaluationScore.average_distance.value:
            return self.average_distance_score
        if score_value == models.EvaluationScore.weakest_link.value:
            return self.weakest_link_score
        if score_value == models.EvaluationScore.ratio.value:
            return self.ratio_score
        raise ValueError(F"score type {score_type} (value: {score_value}) is not supported")
```

Replace the entire file with:

```python
import typing
from dataclasses import dataclass, field

from modules import models


@dataclass
class EvaluationModuleResult:
    sequence: str
    average_distance_score: float
    weakest_link_score: float
    ratio_score: float
    organisms_dist_scores: typing.List[typing.Dict[str, typing.Any]] = field(default_factory=list)

    @property
    def summary(self) -> typing.Dict[str, typing.Any]:
        return {
            "final_sequence": self.sequence,
            "average_distance_score": self.average_distance_score,
            "weakest_link_score": self.weakest_link_score,
            "ratio_score": self.ratio_score,
            "organisms": self.organisms_dist_scores,
        }

    def get_score(self, score_type: models.EvaluationScore) -> float:
        score_value = score_type.value
        if score_value == models.EvaluationScore.average_distance.value:
            return self.average_distance_score
        if score_value == models.EvaluationScore.weakest_link.value:
            return self.weakest_link_score
        if score_value == models.EvaluationScore.ratio.value:
            return self.ratio_score
        raise ValueError(F"score type {score_type} (value: {score_value}) is not supported")
```

- [ ] **Step 2: Pass `organisms_dist_scores` when constructing `EvaluationModuleResult`**

Open `app/modules/evaluation/evaluation.py`. Find the block that constructs `evaluation_result` (around line 101). Currently:

```python
        evaluation_result = models.EvaluationModuleResult(
            sequence=final_sequence,
            average_distance_score=average_distance_score,
            weakest_link_score=weakest_link_score,
            ratio_score=ratio_score,
        )
```

Replace with:

```python
        evaluation_result = models.EvaluationModuleResult(
            sequence=final_sequence,
            average_distance_score=average_distance_score,
            weakest_link_score=weakest_link_score,
            ratio_score=ratio_score,
            organisms_dist_scores=[
                {
                    "name": o["name"],
                    "is_wanted": o["is_wanted"],
                    "dist_score": o["dist_score"],
                }
                for o in organisms_evaluation_summary
            ],
        )
```

`organisms_evaluation_summary` is built just above this block (the loop over `module_input.organisms`) and already contains the `name`, `is_wanted`, and `dist_score` keys.

- [ ] **Step 3: Verify the backend change manually**

Start the FastAPI server:
```bash
cd app && poetry run python api_server.py
```

In a separate terminal, send a minimal test request (adjust paths and sequence to match your local test data):
```bash
curl -s -X POST http://localhost:8000/run-modules \
  -H "Content-Type: application/json" \
  -d '{
    "user_input_dict": {
      "sequence": "ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTT",
      "organisms": {
        "E. coli": {
          "genome_path": "/path/to/ecoli.gb",
          "optimized": true,
          "optimization_priority": 1
        }
      },
      "tuning_param": 0.5,
      "clusters_count": 1,
      "orf_optimization_method": "single_codon_zscore_ratio",
      "orf_optimization_cub_index": "both",
      "initiation_optimization_method": "original",
      "output_path": "results/test",
      "evaluation_score": "average_distance"
    }
  }' | python3 -m json.tool | grep -A 20 '"final_evaluation"'
```

Expected: `final_evaluation` in the response contains an `"organisms"` array with objects like:
```json
{ "name": "E. coli", "is_wanted": true, "dist_score": 1.23 }
```

- [ ] **Step 4: Commit**

```bash
git add app/modules/evaluation/models.py app/modules/evaluation/evaluation.py
git commit -m "feat: expose per-organism dist_scores in final_evaluation API response"
```

---

### Task 2: Frontend — extend types and parse per-organism data

**Files:**
- Modify: `ui/DCUB/lib/types.ts`
- Modify: `ui/DCUB/app/page.tsx`

**Interfaces:**
- Consumes: `result.final_evaluation.organisms` — array of `{ name, is_wanted, dist_score }` from Task 1
- Produces: `OptimizationResult.organisms_dist_scores: Array<{ name: string; is_wanted: boolean; dist_score: number }>` — consumed by Task 3

---

- [ ] **Step 1: Extend `OptimizationResult` type**

Open `ui/DCUB/lib/types.ts`. Add the new field to `OptimizationResult`:

```typescript
export interface OptimizationResult {
  optimized_sequence: string
  evaluation_scores: {
    average_distance_score: number
    ratio_score: number
    weakest_link_score: number
  }
  original_sequence: string
  optimization_parameters: {
    tuning_parameter: number
    optimization_method: string
    cub_index: string
  }
  processing_time: number
  timestamp: string
  organisms_dist_scores: Array<{
    name: string
    is_wanted: boolean
    dist_score: number
  }>
}
```

- [ ] **Step 2: Extract `organisms_dist_scores` in the response parser**

Open `ui/DCUB/app/page.tsx`. Find `parseOptimizationResponse`. It currently builds and returns an object with fields like `optimized_sequence`, `evaluation_scores`, etc. Add the new field to the returned object:

```typescript
organisms_dist_scores: (optimization_result.final_evaluation?.organisms ?? []).map((o: any) => ({
  name: o.name,
  is_wanted: o.is_wanted,
  dist_score: o.dist_score,
})),
```

Place this after `timestamp: ...` in the return statement. The full return object should now include this field.

- [ ] **Step 3: Verify TypeScript compiles without errors**

```bash
cd ui/DCUB && npm run build 2>&1 | head -40
```

Expected: build completes with no type errors. (There will be a "missing prop" error once Task 3's component is imported and used, but until then this should be clean.)

- [ ] **Step 4: Commit**

```bash
git add ui/DCUB/lib/types.ts ui/DCUB/app/page.tsx
git commit -m "feat: parse per-organism dist_scores from API response"
```

---

### Task 3: Create `OrganismDistChart` component

**Files:**
- Create: `ui/DCUB/components/organism-dist-chart.tsx`

**Interfaces:**
- Consumes: `organisms_dist_scores: Array<{ name: string; is_wanted: boolean; dist_score: number }>` from `OptimizationResult` (Task 2)
- Produces: `<OrganismDistChart organismsDistScores={...} />` — consumed by Task 4

---

- [ ] **Step 1: Create the component file**

Create `ui/DCUB/components/organism-dist-chart.tsx` with the following content:

```typescript
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
```

- [ ] **Step 2: Verify TypeScript compiles**

```bash
cd ui/DCUB && npm run build 2>&1 | head -40
```

Expected: no errors from `organism-dist-chart.tsx`. (There may be a warning that the component is not yet imported anywhere — that's fine.)

- [ ] **Step 3: Commit**

```bash
git add ui/DCUB/components/organism-dist-chart.tsx
git commit -m "feat: add OrganismDistChart component for per-organism CUB difference"
```

---

### Task 4: Mount `OrganismDistChart` in the results screen

**Files:**
- Modify: `ui/DCUB/components/results-screen.tsx`

**Interfaces:**
- Consumes: `result.organisms_dist_scores` from `OptimizationResult` (Task 2)
- Consumes: `<OrganismDistChart />` from Task 3

---

- [ ] **Step 1: Import the new component**

Open `ui/DCUB/components/results-screen.tsx`. Add the import at the top with the other component imports:

```typescript
import { OrganismDistChart } from "@/components/organism-dist-chart"
```

- [ ] **Step 2: Mount the chart between the top grid and Optimization Details**

In `results-screen.tsx`, find the JSX structure. It currently looks like:

```tsx
{/* Main Results */}
<div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
  {/* Evaluation Scores */}
  ...
  {/* Optimized Sequence */}
  ...
</div>

{/* Optimization Parameters */}
<Card>
  ...
</Card>
```

Insert `<OrganismDistChart>` between these two blocks:

```tsx
{/* Main Results */}
<div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
  {/* Evaluation Scores */}
  ...
  {/* Optimized Sequence */}
  ...
</div>

{/* Per-Organism CUB Difference */}
<OrganismDistChart organismsDistScores={result.organisms_dist_scores} />

{/* Optimization Parameters */}
<Card>
  ...
</Card>
```

- [ ] **Step 3: Verify TypeScript compiles cleanly**

```bash
cd ui/DCUB && npm run build 2>&1 | head -40
```

Expected: build completes with zero type errors.

- [ ] **Step 4: Run the dev server and visually verify**

```bash
cd ui/DCUB && npm run dev
```

Open `http://localhost:3000`. Run an optimization with at least one wanted and one unwanted organism. On the results page, verify:

1. The "Per-Organism CUB Difference" card appears between the Evaluation Scores/Sequence row and the Optimization Details card.
2. Wanted organisms appear in a "Wanted Organisms" section (green bars/badges), unwanted in an "Unwanted Organisms" section (red bars/badges).
3. A wanted organism with a positive dist_score shows a solid `green-500` bar extending right from center.
4. An unwanted organism with a negative dist_score shows a solid `red-500` bar extending left from center.
5. A muted color (`green-200` or `red-200`) appears when the score went the wrong direction.
6. The center zero line is visible as a thin vertical gray rule.
7. Score values are shown to 2 decimal places with a `+` prefix for positive values.
8. Organism names that are too long to fit truncate with an ellipsis; hovering shows the full name.

- [ ] **Step 5: Commit**

```bash
git add ui/DCUB/components/results-screen.tsx
git commit -m "feat: mount per-organism CUB difference chart in results screen"
```
