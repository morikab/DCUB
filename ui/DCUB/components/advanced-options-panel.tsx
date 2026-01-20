"use client"

import { useState, useEffect } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Label } from "@/components/ui/label"
import { Slider } from "@/components/ui/slider"
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"
import { RadioGroup, RadioGroupItem } from "@/components/ui/radio-group"
import { X, Settings, Info } from "lucide-react"
import { useOptimizationStore } from "@/lib/store"
import { cn } from "@/lib/utils"

interface AdvancedOptionsPanelProps {
  isOpen: boolean
  onClose: () => void
}

export function AdvancedOptionsPanel({ isOpen, onClose }: AdvancedOptionsPanelProps) {
  const { tuningParameter, optimizationMethod, cubIndex, setTuningParameter, setOptimizationMethod, setCubIndex } =
    useOptimizationStore()

  const [localTuningParameter, setLocalTuningParameter] = useState([tuningParameter])

  useEffect(() => {
    setLocalTuningParameter([tuningParameter])
  }, [tuningParameter])

  const handleTuningParameterChange = (value: number[]) => {
    setLocalTuningParameter(value)
    setTuningParameter(value[0])
  }

  const optimizationMethods = [
    { value: "single_codon_diff", label: "Single Codon Diff" },
    { value: "single_codon_ratio", label: "Single Codon Ratio" },
    { value: "zscore_bulk_diff", label: "Z-Score Bulk Diff" },
    { value: "zscore_bulk_ratio", label: "Z-Score Bulk Ratio" },
    { value: "zscore_single_diff", label: "Z-Score Single Diff" },
    { value: "zscore_single_ratio", label: "Z-Score Single Ratio" },
  ]

  const resetToDefaults = () => {
    setTuningParameter(50)
    setOptimizationMethod("single_codon_diff")
    setCubIndex("CAI")
    setLocalTuningParameter([50])
  }

  return (
    <>
      {/* Backdrop */}
      <div
        className={cn(
          "fixed inset-0 bg-black/50 z-40 transition-opacity duration-300",
          isOpen ? "opacity-100" : "opacity-0 pointer-events-none",
        )}
        onClick={onClose}
      />

      {/* Side Panel */}
      <div
        className={cn(
          "fixed right-0 top-0 h-full w-96 bg-white shadow-2xl z-50 transform transition-transform duration-300 ease-in-out overflow-y-auto",
          isOpen ? "translate-x-0" : "translate-x-full",
        )}
      >
        <div className="p-6 space-y-6">
          {/* Header */}
          <div className="flex items-center justify-between border-b pb-4">
            <div className="flex items-center gap-2">
              <Settings className="w-5 h-5 text-blue-600" />
              <h2 className="text-xl font-semibold">Advanced Options</h2>
            </div>
            <Button variant="ghost" size="sm" onClick={onClose}>
              <X className="w-4 h-4" />
            </Button>
          </div>

          {/* Tuning Parameter */}
          <Card>
            <CardHeader className="pb-4">
              <CardTitle className="text-base flex items-center gap-2">
                Tuning Parameter
                <div className="group relative">
                  <Info className="w-4 h-4 text-gray-400 cursor-help" />
                  <div className="absolute left-6 top-0 w-64 p-2 bg-gray-900 text-white text-xs rounded opacity-0 group-hover:opacity-100 transition-opacity z-10">
                    Controls the tradeoff between optimization and deoptimization of wanted and unwanted hosts. Higher
                    values favor wanted organisms more strongly.
                  </div>
                </div>
              </CardTitle>
              <CardDescription>Range: 1-100 (Default: 50)</CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="space-y-3">
                <div className="flex justify-between text-sm text-gray-600">
                  <span>Deoptimization Focus</span>
                  <span>Optimization Focus</span>
                </div>
                <Slider
                  value={localTuningParameter}
                  onValueChange={handleTuningParameterChange}
                  max={100}
                  min={1}
                  step={1}
                  className="w-full"
                />
                <div className="flex justify-between text-xs text-gray-500">
                  <span>1</span>
                  <span className="font-medium text-blue-600">Current: {localTuningParameter[0]}</span>
                  <span>100</span>
                </div>
              </div>
            </CardContent>
          </Card>

          {/* Optimization Method */}
          <Card>
            <CardHeader className="pb-4">
              <CardTitle className="text-base flex items-center gap-2">
                Optimization Method
                <div className="group relative">
                  <Info className="w-4 h-4 text-gray-400 cursor-help" />
                  <div className="absolute left-6 top-0 w-64 p-2 bg-gray-900 text-white text-xs rounded opacity-0 group-hover:opacity-100 transition-opacity z-10">
                    Different algorithms for calculating codon optimization scores and making sequence modifications.
                  </div>
                </div>
              </CardTitle>
              <CardDescription>Select the algorithm for sequence optimization</CardDescription>
            </CardHeader>
            <CardContent>
              <Select value={optimizationMethod} onValueChange={setOptimizationMethod}>
                <SelectTrigger>
                  <SelectValue placeholder="Select optimization method" />
                </SelectTrigger>
                <SelectContent>
                  {optimizationMethods.map((method) => (
                    <SelectItem key={method.value} value={method.value}>
                      {method.label}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>
            </CardContent>
          </Card>

          {/* CUB Index */}
          <Card>
            <CardHeader className="pb-4">
              <CardTitle className="text-base flex items-center gap-2">
                Codon Usage Bias (CUB) Index
                <div className="group relative">
                  <Info className="w-4 h-4 text-gray-400 cursor-help" />
                  <div className="absolute left-6 top-0 w-64 p-2 bg-gray-900 text-white text-xs rounded opacity-0 group-hover:opacity-100 transition-opacity z-10">
                    CAI: Codon Adaptation Index - measures codon usage bias relative to highly expressed genes. tAI:
                    tRNA Adaptation Index - considers tRNA availability and efficiency.
                  </div>
                </div>
              </CardTitle>
              <CardDescription>Choose the index for measuring codon usage bias</CardDescription>
            </CardHeader>
            <CardContent>
              <RadioGroup value={cubIndex} onValueChange={setCubIndex} className="space-y-3">
                <div className="flex items-center space-x-2">
                  <RadioGroupItem value="CAI" id="cai" />
                  <Label htmlFor="cai" className="cursor-pointer">
                    <div>
                      <div className="font-medium">Codon Adaptation Index (CAI)</div>
                      <div className="text-sm text-gray-500">Based on codon usage in highly expressed genes</div>
                    </div>
                  </Label>
                </div>
                <div className="flex items-center space-x-2">
                  <RadioGroupItem value="tAI" id="tai" />
                  <Label htmlFor="tai" className="cursor-pointer">
                    <div>
                      <div className="font-medium">tRNA Adaptation Index (tAI)</div>
                      <div className="text-sm text-gray-500">Considers tRNA availability and decoding efficiency</div>
                    </div>
                  </Label>
                </div>
              </RadioGroup>
            </CardContent>
          </Card>

          {/* Current Settings Summary */}
          <Card className="bg-blue-50 border-blue-200">
            <CardHeader className="pb-3">
              <CardTitle className="text-base text-blue-900">Current Settings</CardTitle>
            </CardHeader>
            <CardContent className="text-sm space-y-2">
              <div className="flex justify-between">
                <span className="text-blue-700">Tuning Parameter:</span>
                <span className="font-medium text-blue-900">{tuningParameter}</span>
              </div>
              <div className="flex justify-between">
                <span className="text-blue-700">Method:</span>
                <span className="font-medium text-blue-900">
                  {optimizationMethods.find((m) => m.value === optimizationMethod)?.label}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-blue-700">CUB Index:</span>
                <span className="font-medium text-blue-900">{cubIndex}</span>
              </div>
            </CardContent>
          </Card>

          {/* Action Buttons */}
          <div className="space-y-3 pt-4 border-t">
            <Button onClick={resetToDefaults} variant="outline" className="w-full bg-transparent">
              Reset to Defaults
            </Button>
            <Button onClick={onClose} className="w-full">
              Apply Settings
            </Button>
          </div>
        </div>
      </div>
    </>
  )
}
