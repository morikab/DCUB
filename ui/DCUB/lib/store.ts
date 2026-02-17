import { create } from "zustand"
import { persist } from "zustand/middleware"
import type { Organism } from "./types"

interface OptimizationState {
  // DNA Sequence
  dnaSequence: string
  sequenceFile: File | null

  // Organisms
  wantedOrganisms: Organism[]
  unwantedOrganisms: Organism[]

  // Advanced Options
  tuningParameter: number
  optimizationMethod: string
  cubIndex: string

  // Actions
  setDnaSequence: (sequence: string) => void
  setSequenceFile: (file: File | null) => void

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

  reset: () => void
}

const initialState = {
  dnaSequence: "",
  sequenceFile: null,
  wantedOrganisms: [],
  unwantedOrganisms: [],
  tuningParameter: 50,
  optimizationMethod: "single_codon_diff",
  cubIndex: "CAI",
}

export const useOptimizationStore = create<OptimizationState>()(
  persist(
    (set) => ({
      ...initialState,

      setDnaSequence: (sequence) => set({ dnaSequence: sequence }),
      setSequenceFile: (file) => set({ sequenceFile: file }),

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

      reset: () => set(initialState),
    }),
    {
      name: "dna-optimization-storage",
      partialize: (state) => ({
        dnaSequence: state.dnaSequence,
        wantedOrganisms: state.wantedOrganisms,
        unwantedOrganisms: state.unwantedOrganisms,
        tuningParameter: state.tuningParameter,
        optimizationMethod: state.optimizationMethod,
        cubIndex: state.cubIndex,
      }),
    },
  ),
)
