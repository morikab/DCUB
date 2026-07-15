# Per-Organism CUB Difference — Design Spec

**Date:** 2026-07-15
**Branch:** fix/lab-meeting-feedback

## Summary

Add a per-organism CUB difference breakdown to the results page. For the winning optimized sequence, show each organism's individual dist_score (the z-score of the CUB change relative to that organism's proteome variability) in a diverging horizontal bar chart, grouped by role (wanted / unwanted).

---

## Background

The evaluation module already computes `dist_score = (final_cub_score - initial_cub_score) / sigma` for every organism during `EvaluationModule.run_module()`. This data is stored in `organisms_evaluation_summary` and appended to the run_summary under `"evaluation"`, but the `"final_evaluation"` key (what the API surfaces to the frontend) only exposes the three aggregate scores. The per-organism breakdown is computed but discarded before reaching the UI.

---

## Backend Changes

### `app/modules/evaluation/models.py`

Add a new field to `EvaluationModuleResult`:

```python
organisms_dist_scores: typing.List[typing.Dict[str, typing.Any]]
```

Extend the `summary` property to include it:

```python
"organisms": self.organisms_dist_scores,
```

### `app/modules/evaluation/evaluation.py`

Pass `organisms_evaluation_summary` into `EvaluationModuleResult`, stripping the raw CUB scores — only `name`, `is_wanted`, and `dist_score` are needed:

```python
organisms_dist_scores=[
    {"name": o["name"], "is_wanted": o["is_wanted"], "dist_score": o["dist_score"]}
    for o in organisms_evaluation_summary
]
```

### Resulting API shape

`result.final_evaluation` in the API response will include:

```json
{
  "average_distance_score": 1.2,
  "weakest_link_score": 0.8,
  "ratio_score": 1.5,
  "final_sequence": "ATG...",
  "organisms": [
    { "name": "E. coli", "is_wanted": true, "dist_score": 1.4 },
    { "name": "B. subtilis", "is_wanted": false, "dist_score": -2.1 }
  ]
}
```

No other backend changes needed.

---

## Frontend Changes

### `ui/DCUB/lib/types.ts`

Extend `OptimizationResult` with:

```typescript
organisms_dist_scores: Array<{
  name: string
  is_wanted: boolean
  dist_score: number
}>
```

### `ui/DCUB/app/page.tsx`

In `parseOptimizationResponse`, extract the new field:

```typescript
organisms_dist_scores: (optimization_result.final_evaluation.organisms || []).map((o: any) => ({
  name: o.name,
  is_wanted: o.is_wanted,
  dist_score: o.dist_score,
})),
```

### New component: `ui/DCUB/components/organism-dist-chart.tsx`

A self-contained component that receives `organisms_dist_scores` as a prop and renders a full-width card.

**Card:**
- Title: **"Per-Organism CUB Difference"**
- Subtitle: *"z-score of CUB change relative to each organism's proteome variability — positive = improved expression, negative = reduced."*

**Organism grouping:**
Render organisms in two sections separated by a subtle divider:
1. **Wanted Organisms** (section label)
2. **Unwanted Organisms** (section label)

Within each section, organisms are listed in the order returned by the backend.

**Row layout (per organism):**
```
[name + badge (~25%)] | [diverging bar area (~60%)] | [score value (~15%)]
```

- Left: organism name + small role badge ("Wanted" in green, "Unwanted" in red)
- Center: diverging bar with hairline center zero line
- Right: dist_score formatted to 2 decimal places

**Bar rendering:**
- `maxAbs = Math.max(...allOrganisms.map(o => Math.abs(o.dist_score)))` — single scale across both groups
- Bar width = `|dist_score| / maxAbs * 50%` of the center column
- Positive bars: `position: absolute; left: 50%` (extends right)
- Negative bars: `position: absolute; right: 50%` (extends left)
- Center line: `position: absolute; left: 50%; width: 1px; height: 100%; bg-gray-300`

**Color scheme (role + outcome):**

| Role | Score direction | Meaning | Color |
|------|----------------|---------|-------|
| Wanted | positive | CUB improved ✓ | `bg-green-500` |
| Wanted | negative | CUB worsened ✗ | `bg-green-200` (muted) |
| Unwanted | negative | CUB reduced ✓ | `bg-red-500` |
| Unwanted | positive | CUB increased ✗ | `bg-red-200` (muted) |

### `ui/DCUB/components/results-screen.tsx`

Add the new `OrganismDistChart` card between the top grid row (Evaluation Scores / Optimized Sequence) and the Optimization Details card.

Pass `result.organisms_dist_scores` as the prop.

---

## File Change Summary

| File | Change type |
|------|-------------|
| `app/modules/evaluation/models.py` | Add field + extend summary |
| `app/modules/evaluation/evaluation.py` | Pass organisms_dist_scores to result |
| `ui/DCUB/lib/types.ts` | Extend OptimizationResult type |
| `ui/DCUB/app/page.tsx` | Extract organisms_dist_scores in parser |
| `ui/DCUB/components/organism-dist-chart.tsx` | New component |
| `ui/DCUB/components/results-screen.tsx` | Mount new card |

---

## Out of Scope

- Showing per-organism data for runner-up candidates (only winning sequence)
- Axis tick labels / gridlines (CSS approach; can upgrade to SVG later if needed)
- Sorting organisms by score magnitude within a group
