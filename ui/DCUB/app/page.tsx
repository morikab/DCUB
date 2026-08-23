"use client"

import { useState } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"
import { Alert, AlertDescription } from "@/components/ui/alert"
import { Settings, HelpCircle } from "lucide-react"
import { useOptimizationStore } from "@/lib/store"
import { DnaSequenceInput } from "@/components/dna-sequence-input"
import { OrganismList } from "@/components/organism-list"
import { validateSubmission } from "@/lib/validation"
import { AdvancedOptionsPanel } from "@/components/advanced-options-panel"
import { LoadingScreen } from "@/components/loading-screen"
import { ResultsScreen } from "@/components/results-screen"
import { ErrorDialog } from "@/components/error-dialog"
import type { OptimizationResult, OrganismRequestPayload, RawOptimizationResponse } from "@/lib/types"

export default function DNAOptimizerPage() {
  const [isSubmitting, setIsSubmitting] = useState(false)
  const [submitStatus, setSubmitStatus] = useState<"idle" | "success">("idle")
  const [submitMessage, setSubmitMessage] = useState("")
  const [errorDialogOpen, setErrorDialogOpen] = useState(false)
  const [errorText, setErrorText] = useState("")
  const [showAdvancedOptions, setShowAdvancedOptions] = useState(false)
  const [showHelp, setShowHelp] = useState(false)
  const [showHelpTooltip, setShowHelpTooltip] = useState(false)
  const [showAdvTooltip, setShowAdvTooltip] = useState(false)
  const [optimizationResult, setOptimizationResult] = useState<OptimizationResult | null>(null)
  const [showResults, setShowResults] = useState(false)

  const { dnaSequence, sequenceFile, sequenceFilePath, wantedOrganisms, unwantedOrganisms, reset } =
    useOptimizationStore()

  const handleSubmit = async () => {
    setIsSubmitting(true)
    setSubmitStatus("idle")
    setOptimizationResult(null)
    setShowResults(false)

    try {
      // Validate all inputs
      const validation = validateSubmission({
        dnaSequence,
        sequenceFile,
        wantedOrganisms,
        unwantedOrganisms,
      })

      if (!validation.isValid) {
        setErrorText(validation.errors?.join(", ") || "Validation failed")
        setErrorDialogOpen(true)
        return
      }

      // Get current advanced options
      const currentState = useOptimizationStore.getState()

      // Prepare organisms object in the new format
      const organismsObject: Record<string, OrganismRequestPayload> = {}

      // Add wanted organisms (optimized: true)
      wantedOrganisms.forEach((org) => {
        // Use organism name as key, fallback to generated name if empty
        const organismKey = org.name || `wanted-organism-${org.id}`
        organismsObject[organismKey] = {
          genome_path: org.genomePath, // Use full path
          optimized: true,
          expression_data_type: org.expressionDataType,
          expression_file_format: org.expressionDataFormat,
          expression_file_path: org.expressionDataPath || null, // Use full path
          optimization_priority: org.priority,
        }
      })

      // Add unwanted organisms (optimized: false)
      unwantedOrganisms.forEach((org) => {
        // Use organism name as key, fallback to generated name if empty
        const organismKey = org.name || `unwanted-organism-${org.id}`
        organismsObject[organismKey] = {
          genome_path: org.genomePath, // Use full path
          optimized: false,
          expression_data_type: org.expressionDataType,
          expression_file_format: org.expressionDataFormat,
          expression_file_path: org.expressionDataPath || null, // Use full path
          optimization_priority: org.priority,
        }
      })

      // Prepare JSON payload in the new structure
      const optimizationPayload = {
        user_input_dict: {
          sequence_file_path: sequenceFile ? sequenceFilePath : null,
          sequence: dnaSequence,
          tuning_param: currentState.tuningParameter / 100, // Convert to 0-1 range
          organisms: organismsObject,
          clusters_count: 1,
          orf_optimization_method: currentState.optimizationMethod,
          orf_optimization_cub_index: currentState.cubIndex,
          initiation_optimization_method: "original",
          output_path: `results/DCUB/${Date.now()}`,
          evaluation_score: "average_distance",
          enable_hotspot_avoidance: currentState.enableHotspotAvoidance,
          enable_motif_detection: currentState.enableMotifDetection,
        },
      }

      console.log("Sending optimization request:", optimizationPayload)

      // Send POST request to backend with extended timeout
      const controller = new AbortController()
      const timeoutId = setTimeout(() => controller.abort(), 300000) // 5 minute timeout

      const response = await fetch("http://localhost:8000/run-modules", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Accept: "application/json",
        },
        body: JSON.stringify(optimizationPayload),
        signal: controller.signal,
      })

      clearTimeout(timeoutId)

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({ detail: "Unknown error" }))
        throw new Error(errorData.detail || errorData.error || `HTTP ${response.status}: ${response.statusText}`)
      }

      const raw_response = await response.json()
      const result = raw_response?.result
      console.log("Optimization result:", result)

      // Parse and validate the response
      const parsedResult = parseOptimizationResponse(result)
      setOptimizationResult(parsedResult)
      setShowResults(true)
      setSubmitStatus("success")
      setSubmitMessage("DNA optimization completed successfully!")
      
      // Clear all form fields after successful optimization
      reset()
    } catch (error) {
      console.error("Optimization error:", error)

      if (error instanceof Error && error.name === "AbortError") {
        setErrorText("Request timed out. The optimization process took too long to complete.")
      } else if (error instanceof TypeError && error.message.includes("fetch")) {
        setErrorText(
          "Unable to connect to the optimization server. Please ensure the backend is running on localhost:8000 and CORS is properly configured.",
        )
      } else {
        setErrorText(error instanceof Error ? error.message : "Unknown error")
      }
      setErrorDialogOpen(true)
    } finally {
      setIsSubmitting(false)
    }
  }

  // Helper function to parse optimization response
  const parseOptimizationResponse = (optimization_result: RawOptimizationResponse): OptimizationResult => {
    console.info(optimization_result.final_evaluation)
    try {
      return {
        optimized_sequence: optimization_result.final_evaluation.final_sequence || "",
        // optimized_sequence: response.optimized_sequence || response.sequence || "",
        evaluation_scores: {
          average_distance_score: optimization_result.final_evaluation.average_distance_score || 0,
          ratio_score: optimization_result.final_evaluation?.ratio_score || 0,
          weakest_link_score: optimization_result.final_evaluation?.weakest_link_score || 0,
        },
        // evaluation_scores: {
        //  cai_score: response.evaluation_scores?.cai_score || response.cai_score || 0,
        //  gc_content: response.evaluation_scores?.gc_content || response.gc_content || 0,
        //  codon_usage_bias: response.evaluation_scores?.codon_usage_bias || response.codon_usage_bias || 0,
        // },
        original_sequence: optimization_result.original_sequence || dnaSequence || sequenceFile?.name || "",
        optimization_parameters: {
          tuning_parameter: useOptimizationStore.getState().tuningParameter,
          optimization_method: useOptimizationStore.getState().optimizationMethod,
          cub_index: useOptimizationStore.getState().cubIndex,
        },
        processing_time: optimization_result.processing_time || 0,
        timestamp: optimization_result.timestamp || new Date().toISOString(),
        organisms_dist_scores: (optimization_result.final_evaluation?.organisms_dist_scores ?? []).map((o) => ({
          name: o.name,
          is_wanted: o.is_wanted,
          dist_score: o.dist_score,
        })),
        hotspot_avoidance: optimization_result.hotspot_avoidance
          ? {
              enabled: optimization_result.hotspot_avoidance.enabled ?? false,
              sequence_before: optimization_result.hotspot_avoidance.sequence_before || "",
              sequence_after: optimization_result.hotspot_avoidance.sequence_after || "",
              num_edits: optimization_result.hotspot_avoidance.num_edits || 0,
              detected_sites: {
                recombination: optimization_result.hotspot_avoidance.detected_sites?.recombination || 0,
                slippage: optimization_result.hotspot_avoidance.detected_sites?.slippage || 0,
                motifs: optimization_result.hotspot_avoidance.detected_sites?.motifs || 0,
              },
              detected_regions: optimization_result.hotspot_avoidance.detected_regions ?? [],
              rounds: optimization_result.hotspot_avoidance.rounds ?? 0,
              residual_regions: optimization_result.hotspot_avoidance.residual_regions ?? [],
              warnings: optimization_result.hotspot_avoidance.warnings ?? [],
            }
          : undefined,
      }
    } catch (error) {
      console.error("Error parsing response:", error)
      console.error("Result:", optimization_result)
      throw new Error("Invalid response format from server")
    }
  }

  const handleReset = () => {
    reset()
    setSubmitStatus("idle")
    setSubmitMessage("")
    setOptimizationResult(null)
    setShowResults(false)
  }

  const handleBackToForm = () => {
    setShowResults(false)
    setOptimizationResult(null)
    setSubmitStatus("idle")
    setSubmitMessage("")
    // Clear all form fields when going back to form
    reset()
  }

  // Show loading screen while processing
  if (isSubmitting) {
    return <LoadingScreen />
  }

  // Show results screen after successful optimization
  if (showResults && optimizationResult) {
    return <ResultsScreen result={optimizationResult} onBackToForm={handleBackToForm} onReset={handleReset} />
  }

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100 p-4">
      {/* Advanced Options Side Button */}
      <div className="fixed right-4 top-1/2 transform -translate-y-1/2 z-50">
        <div className="relative">
          <Button
            onClick={() => setShowAdvancedOptions(true)}
            onMouseEnter={() => setShowAdvTooltip(true)}
            onMouseLeave={() => setShowAdvTooltip(false)}
            className="rounded-full w-12 h-12 shadow-lg bg-blue-600 hover:bg-blue-700"
            size="sm"
          >
            <Settings className="w-5 h-5" />
          </Button>
          {showAdvTooltip && (
            <div style={{ position: "absolute", right: "100%", marginRight: 8, top: "50%", transform: "translateY(-50%)", backgroundColor: "#1f2937", color: "white", fontSize: "0.75rem", borderRadius: 4, padding: "2px 8px", whiteSpace: "nowrap", pointerEvents: "none" }}>
              Advanced options
            </div>
          )}
        </div>
      </div>

      {/* Advanced Options Panel */}
      <AdvancedOptionsPanel isOpen={showAdvancedOptions} onClose={() => setShowAdvancedOptions(false)} />

      <div className="max-w-6xl mx-auto space-y-6">
        {/* Header */}
        <div className="text-center space-y-2">
          <h1 className="text-4xl font-bold text-gray-900">DCUB</h1>
          <h1 className="text-xl font-bold text-gray-600">Your Differential Codon Usage Bias tool!</h1>
          <p className="text-lg text-gray-600 max-w-2xl mx-auto">
            Configure your DNA sequence optimization parameters by specifying target organisms, priorities, and
            expression data for comprehensive genetic optimization using GenBank genome files.
          </p>
        </div>

        {/* Main Content */}
        <Card className="shadow-lg">
          <CardHeader>
            <div className="flex items-center justify-between">
              <CardTitle className="flex items-center gap-2">
                <div className="w-2 h-2 bg-blue-500 rounded-full"></div>
                Optimization Configuration
              </CardTitle>
              <div className="relative">
                <Button
                  variant="ghost"
                  size="sm"
                  onClick={() => setShowHelp((v) => !v)}
                  onMouseEnter={() => setShowHelpTooltip(true)}
                  onMouseLeave={() => setShowHelpTooltip(false)}
                  className="text-gray-400 hover:text-gray-600"
                  aria-label="Toggle help"
                >
                  <HelpCircle className="w-5 h-5" />
                </Button>
                {showHelpTooltip && (
                  <div style={{ position: "absolute", right: 0, top: "100%", marginTop: 4, backgroundColor: "#1f2937", color: "white", fontSize: "0.75rem", borderRadius: 4, padding: "2px 8px", whiteSpace: "nowrap", zIndex: 10, pointerEvents: "none" }}>
                    Instructions
                  </div>
                )}
              </div>
            </div>
            <CardDescription>Configure all parameters for your DNA sequence optimization process</CardDescription>
          </CardHeader>

          {showHelp && (
            <div className="px-6 pb-2">
              <div className="rounded-lg border border-blue-200 bg-blue-50 p-4 text-sm space-y-3">
                <div className="flex gap-3">
                  <span className="inline-flex items-center justify-center w-5 h-5 rounded-full text-white text-xs font-bold flex-shrink-0 mt-0.5" style={{ backgroundColor: "#9ca3af" }}>1</span>
                  <div className="text-blue-800">
                    <span className="font-semibold">DNA Sequence — </span>
                    Enter a raw sequence or FASTA format directly, or upload a .fasta / .fa file.
                  </div>
                </div>
                <div className="flex gap-3">
                  <span className="inline-flex items-center justify-center w-5 h-5 rounded-full text-white text-xs font-bold flex-shrink-0 mt-0.5" style={{ backgroundColor: "#16a34a" }}>2</span>
                  <div className="text-blue-800">
                    <span className="font-semibold">Wanted Organisms — </span>
                    Organisms you want your gene expressed in. Each needs a GenBank genome file (.gb), a priority score (1–100; higher = more weight), and optionally a JSON or CSV file with gene-level expression data.
                  </div>
                </div>
                <div className="flex gap-3">
                  <span className="inline-flex items-center justify-center w-5 h-5 rounded-full text-white text-xs font-bold flex-shrink-0 mt-0.5" style={{ backgroundColor: "#dc2626" }}>3</span>
                  <div className="text-blue-800">
                    <span className="font-semibold">Unwanted Organisms — </span>
                    Organisms where expression should be minimized. Same fields as step 2.
                  </div>
                </div>
                <div className="border-t border-blue-200 pt-3 text-blue-700 text-xs">
                  <strong>Advanced Options</strong> — use the ⚙ button (floating, right side) to tune the optimization algorithm and parameters.
                </div>
              </div>
            </div>
          )}

          <CardContent className="space-y-6">
            {/* DNA Sequence Input */}
            <DnaSequenceInput />

            {/* Organism Configuration */}
            <Tabs defaultValue="wanted" className="w-full">
              <TabsList className="grid w-full grid-cols-2">
                <TabsTrigger value="wanted" className="flex items-center gap-2">
                  <span className="inline-flex items-center justify-center w-5 h-5 rounded-full text-white text-xs font-bold flex-shrink-0" style={{ backgroundColor: "#16a34a" }}>2</span>
                  Wanted Organisms
                </TabsTrigger>
                <TabsTrigger value="unwanted" className="flex items-center gap-2">
                  <span className="inline-flex items-center justify-center w-5 h-5 rounded-full text-white text-xs font-bold flex-shrink-0" style={{ backgroundColor: "#dc2626" }}>3</span>
                  Unwanted Organisms
                </TabsTrigger>
              </TabsList>

              <TabsContent value="wanted" className="mt-6">
                <OrganismList type="wanted" />
              </TabsContent>

              <TabsContent value="unwanted" className="mt-6">
                <OrganismList type="unwanted" />
              </TabsContent>
            </Tabs>

            {/* Submit Status */}
            {submitStatus === "success" && (
              <Alert className="border-green-200 bg-green-50">
                <AlertDescription className="text-green-800">{submitMessage}</AlertDescription>
              </Alert>
            )}

            {/* Action Buttons */}
            <div className="flex flex-col sm:flex-row gap-4 pt-6 border-t">
              <Button
                onClick={handleSubmit}
                disabled={isSubmitting}
                className="flex-1 h-12 text-lg font-medium"
                size="lg"
              >
                {isSubmitting ? "Processing..." : "Optimize DNA Sequence"}
              </Button>
              <Button
                onClick={handleReset}
                variant="outline"
                className="flex-1 h-12 text-lg font-medium bg-transparent"
                size="lg"
              >
                Reset Configuration
              </Button>
            </div>
          </CardContent>
        </Card>

      </div>

      <ErrorDialog open={errorDialogOpen} errorText={errorText} onClose={() => setErrorDialogOpen(false)} />
    </div>
  )
}
