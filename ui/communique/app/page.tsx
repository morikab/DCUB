"use client"

import { useState } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"
import { Alert, AlertDescription } from "@/components/ui/alert"
import { CheckCircle, AlertCircle, Settings } from "lucide-react"
import { useOptimizationStore } from "@/lib/store"
import { DnaSequenceInput } from "@/components/dna-sequence-input"
import { OrganismList } from "@/components/organism-list"
import { validateSubmission } from "@/lib/validation"
import { AdvancedOptionsPanel } from "@/components/advanced-options-panel"
import { LoadingScreen } from "@/components/loading-screen"
import { ResultsScreen } from "@/components/results-screen"
import type { OptimizationResult } from "@/lib/types"

export default function DNAOptimizerPage() {
  const [isSubmitting, setIsSubmitting] = useState(false)
  const [submitStatus, setSubmitStatus] = useState<"idle" | "success" | "error">("idle")
  const [submitMessage, setSubmitMessage] = useState("")
  const [showAdvancedOptions, setShowAdvancedOptions] = useState(false)
  const [optimizationResult, setOptimizationResult] = useState<OptimizationResult | null>(null)
  const [showResults, setShowResults] = useState(false)

  const { dnaSequence, sequenceFile, wantedOrganisms, unwantedOrganisms, reset } = useOptimizationStore()

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
        setSubmitStatus("error")
        setSubmitMessage(validation.errors.join(", "))
        return
      }

      // Get current advanced options
      const currentState = useOptimizationStore.getState()

      // Prepare organisms object in the new format
      const organismsObject: Record<string, any> = {}

      // Add wanted organisms (optimized: true)
      wantedOrganisms.forEach((org) => {
        // Use organism name as key, fallback to generated name if empty
        const organismKey = org.name || `wanted-organism-${org.id}`
        organismsObject[organismKey] = {
          genome_path: org.genomePath, // Use full path
          optimized: true,
          expression_csv_type: "protein_abundance",
          expression_csv: org.expressionDataPath || null, // Use full path
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
          expression_csv_type: "protein_abundance",
          expression_csv: org.expressionDataPath || null, // Use full path
          optimization_priority: org.priority,
        }
      })

      // Prepare JSON payload in the new structure
      const optimizationPayload = {
        user_input_dict: {
          sequence_file_path: sequenceFile ? sequenceFile.name : null,
          sequence: sequenceFile ? null : dnaSequence,
          tuning_param: currentState.tuningParameter / 100, // Convert to 0-1 range
          organisms: organismsObject,
          clusters_count: 1,
          orf_optimization_method: currentState.optimizationMethod,
          orf_optimization_cub_index: currentState.cubIndex,
          initiation_optimization_method: "original",
          output_path: `results/commuique/${Date.now()}`,
          evaluation_score: "average_distance",
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
        const errorData = await response.json().catch(() => ({ error: "Unknown error" }))
        throw new Error(errorData.error || `HTTP ${response.status}: ${response.statusText}`)
      }

      const result = await response.json()
      console.log("Optimization result:", result)

      // Parse and validate the response
      const parsedResult = parseOptimizationResponse(result)
      setOptimizationResult(parsedResult)
      setShowResults(true)
      setSubmitStatus("success")
      setSubmitMessage("DNA optimization completed successfully!")
    } catch (error) {
      console.error("Optimization error:", error)
      setSubmitStatus("error")

      if (error instanceof Error && error.name === "AbortError") {
        setSubmitMessage("Request timed out. The optimization process took too long to complete.")
      } else if (error instanceof TypeError && error.message.includes("fetch")) {
        setSubmitMessage(
          "Unable to connect to the optimization server. Please ensure the backend is running on localhost:8000 and CORS is properly configured.",
        )
      } else {
        setSubmitMessage(`Optimization failed: ${error instanceof Error ? error.message : "Unknown error"}`)
      }
    } finally {
      setIsSubmitting(false)
    }
  }

  // Helper function to parse optimization response
  const parseOptimizationResponse = (response: any): OptimizationResult => {
    try {
      return {
        optimized_sequence: response.final_evaluation?.final_sequence || "",
        // optimized_sequence: response.optimized_sequence || response.sequence || "",
        evaluation_scores: {
          average_distance_score: response.final_evaluation?.average_distance_score || 0,
          ratio_score: response.final_evaluation?.ratio_score || 0,
          weakest_link_score: response.final_evaluation?.weakest_link_score || 0,
        },
        // evaluation_scores: {
        //  cai_score: response.evaluation_scores?.cai_score || response.cai_score || 0,
        //  gc_content: response.evaluation_scores?.gc_content || response.gc_content || 0,
        //  codon_usage_bias: response.evaluation_scores?.codon_usage_bias || response.codon_usage_bias || 0,
        // },
        original_sequence: response.original_sequence || dnaSequence || sequenceFile?.name || "",
        optimization_parameters: {
          tuning_parameter: useOptimizationStore.getState().tuningParameter,
          optimization_method: useOptimizationStore.getState().optimizationMethod,
          cub_index: useOptimizationStore.getState().cubIndex,
        },
        processing_time: response.processing_time || 0,
        timestamp: response.timestamp || new Date().toISOString(),
      }
    } catch (error) {
      console.error("Error parsing response:", error)
      console.error("Response:", response)
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
        <Button
          onClick={() => setShowAdvancedOptions(true)}
          className="rounded-full w-12 h-12 shadow-lg bg-blue-600 hover:bg-blue-700"
          size="sm"
        >
          <Settings className="w-5 h-5" />
        </Button>
      </div>

      {/* Advanced Options Panel */}
      <AdvancedOptionsPanel isOpen={showAdvancedOptions} onClose={() => setShowAdvancedOptions(false)} />

      <div className="max-w-6xl mx-auto space-y-6">
        {/* Header */}
        <div className="text-center space-y-2">
          <h1 className="text-4xl font-bold text-gray-900">Commuique</h1>
          <p className="text-lg text-gray-600 max-w-2xl mx-auto">
            Configure your DNA sequence optimization parameters by specifying target organisms, priorities, and
            expression data for comprehensive genetic optimization using GenBank genome files.
          </p>
        </div>

        {/* Main Content */}
        <Card className="shadow-lg">
          <CardHeader>
            <CardTitle className="flex items-center gap-2">
              <div className="w-2 h-2 bg-blue-500 rounded-full"></div>
              Optimization Configuration
            </CardTitle>
            <CardDescription>Configure all parameters for your DNA sequence optimization process</CardDescription>
          </CardHeader>
          <CardContent className="space-y-6">
            {/* DNA Sequence Input */}
            <DnaSequenceInput />

            {/* Organism Configuration */}
            <Tabs defaultValue="wanted" className="w-full">
              <TabsList className="grid w-full grid-cols-2">
                <TabsTrigger value="wanted" className="flex items-center gap-2">
                  <CheckCircle className="w-4 h-4" />
                  Wanted Organisms
                </TabsTrigger>
                <TabsTrigger value="unwanted" className="flex items-center gap-2">
                  <AlertCircle className="w-4 h-4" />
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
            {submitStatus !== "idle" && (
              <Alert
                className={submitStatus === "success" ? "border-green-200 bg-green-50" : "border-red-200 bg-red-50"}
              >
                <AlertDescription className={submitStatus === "success" ? "text-green-800" : "text-red-800"}>
                  {submitMessage}
                </AlertDescription>
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

        {/* Instructions */}
        <Card className="bg-blue-50 border-blue-200">
          <CardHeader>
            <CardTitle className="text-blue-900">Instructions</CardTitle>
          </CardHeader>
          <CardContent className="text-blue-800 space-y-2">
            <p>
              <strong>1. DNA Sequence:</strong> Enter your sequence manually or upload a FASTA file
            </p>
            <p>
              <strong>2. Wanted Organisms:</strong> Add organisms you want to optimize for with their GenBank genome
              files (.gb/.gbf) and priorities
            </p>
            <p>
              <strong>3. Unwanted Organisms:</strong> Add organisms you want to avoid with their GenBank genome files
              (.gb/.gbf) and priorities
            </p>
            <p>
              <strong>4. Priority Scores:</strong> Use values between 1-100 (higher = more important)
            </p>
            <p>
              <strong>5. Expression Data:</strong> Optionally provide CSV files with expression data for each organism
            </p>
            <p>
              <strong>6. Advanced Options:</strong> Click the settings button to configure tuning parameters and
              optimization methods
            </p>
          </CardContent>
        </Card>

        {/* CORS Notice for Development */}
        <Card className="bg-yellow-50 border-yellow-200">
          <CardHeader>
            <CardTitle className="text-yellow-900 text-base">Development Notice</CardTitle>
          </CardHeader>
          <CardContent className="text-yellow-800 text-sm">
            <p>
              <strong>For local development:</strong> Ensure your backend server at localhost:8000 has CORS enabled to
              allow requests from this web application.
            </p>
            <p className="mt-2">
              <strong>GenBank Files:</strong> Genome files must be in GenBank format (.gb or .gbf extensions).
            </p>
            <p className="mt-2">
              <strong>Processing Time:</strong> Optimization can take several minutes to complete. Please be patient
              while the server processes your request.
            </p>
            <p className="mt-2">
              <strong>API Endpoint:</strong> Requests are sent to http://localhost:8000/run-modules
            </p>
          </CardContent>
        </Card>
      </div>
    </div>
  )
}
