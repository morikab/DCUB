import { create } from "zustand"
import { persist } from "zustand/middleware"
import type { Organism } from "./types"

interface OptimizationState {
  // DNA Sequence
  dnaSequence: string
  sequenceFile: File | null
  sequenceFilePath: string

  // Organisms
  wantedOrganisms: Organism[]
  unwantedOrganisms: Organism[]

  // Advanced Options
  tuningParameter: number
  optimizationMethod: string
  cubIndex: string
  enableHotspotAvoidance: boolean
  enableMotifDetection: boolean

  // Actions
  setDnaSequence: (sequence: string) => void
  setSequenceFile: (file: File | null) => void
  setSequenceFilePath: (path: string) => void

  addWantedOrganism: (organism: Organism) => void
  addUnwantedOrganism: (organism: Organism) => void
  removeWantedOrganism: (id: string) => void
  removeUnwantedOrganism: (id: string) => void
  updateWantedOrganism: (id: string, organism: Organism) => void
  updateUnwantedOrganism: (id: string, organism: Organism) => void

  // Advanced Options Actions
  setTuningParameter: (value: number) => void
  setOptimizationMethod: (method: string) => void
  setCubIndex: (index: string) => void
  setEnableHotspotAvoidance: (enabled: boolean) => void
  setEnableMotifDetection: (enabled: boolean) => void

  reset: () => void
}

const initialState = {
  dnaSequence: "",
  sequenceFile: null,
  sequenceFilePath: "",
  wantedOrganisms: [],
  unwantedOrganisms: [],
  tuningParameter: 50,
  optimizationMethod: "single_codon_diff",
  cubIndex: "CAI",
  enableHotspotAvoidance: false,
  enableMotifDetection: false,
}

/** Dropdown values that predate the fix for the 422 on every Z-Score method.
 *  A browser that ran the old build has one of these persisted, and would keep
 *  sending it after the fix, so rehydration maps them to the enum values the
 *  backend accepts. */
const RENAMED_OPTIMIZATION_METHODS: Record<string, string> = {
  zscore_bulk_diff: "zscore_bulk_aa_diff",
  zscore_bulk_ratio: "zscore_bulk_aa_ratio",
  zscore_single_diff: "zscore_single_aa_diff",
  zscore_single_ratio: "zscore_single_aa_ratio",
}

export const useOptimizationStore = create<OptimizationState>()(
  persist(
    (set) => ({
      ...initialState,

      setDnaSequence: (sequence) => set({ dnaSequence: sequence }),
      setSequenceFile: (file) => set({ sequenceFile: file }),
      setSequenceFilePath: (path) => set({ sequenceFilePath: path }),

      addWantedOrganism: (organism) =>
        set((state) => ({
          wantedOrganisms: [...state.wantedOrganisms, organism],
        })),

      addUnwantedOrganism: (organism) =>
        set((state) => ({
          unwantedOrganisms: [...state.unwantedOrganisms, organism],
        })),

      removeWantedOrganism: (id) =>
        set((state) => ({
          wantedOrganisms: state.wantedOrganisms.filter((org) => org.id !== id),
        })),

      removeUnwantedOrganism: (id) =>
        set((state) => ({
          unwantedOrganisms: state.unwantedOrganisms.filter((org) => org.id !== id),
        })),

      updateWantedOrganism: (id, organism) =>
        set((state) => ({
          wantedOrganisms: state.wantedOrganisms.map((org) => (org.id === id ? organism : org)),
        })),

      updateUnwantedOrganism: (id, organism) =>
        set((state) => ({
          unwantedOrganisms: state.unwantedOrganisms.map((org) => (org.id === id ? organism : org)),
        })),

      setTuningParameter: (value) => set({ tuningParameter: value }),
      setOptimizationMethod: (method) => set({ optimizationMethod: method }),
      setCubIndex: (index) => set({ cubIndex: index }),
      setEnableHotspotAvoidance: (enabled) => set({ enableHotspotAvoidance: enabled }),
      setEnableMotifDetection: (enabled) => set({ enableMotifDetection: enabled }),

      reset: () => set(initialState),
    }),
    {
      name: "dna-optimization-storage",
      version: 1,
      migrate: (persistedState, version) => {
        if (version >= 1 || typeof persistedState !== "object" || persistedState === null) {
          return persistedState as OptimizationState
        }
        const state = persistedState as Partial<OptimizationState>
        const renamed = state.optimizationMethod
          ? RENAMED_OPTIMIZATION_METHODS[state.optimizationMethod]
          : undefined
        return (renamed ? { ...state, optimizationMethod: renamed } : state) as OptimizationState
      },
      partialize: (state) => ({
        dnaSequence: state.dnaSequence,
        wantedOrganisms: state.wantedOrganisms,
        unwantedOrganisms: state.unwantedOrganisms,
        tuningParameter: state.tuningParameter,
        optimizationMethod: state.optimizationMethod,
        cubIndex: state.cubIndex,
        enableHotspotAvoidance: state.enableHotspotAvoidance,
        enableMotifDetection: state.enableMotifDetection,
      }),
    },
  ),
)
