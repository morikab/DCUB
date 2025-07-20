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

export default function DNAOptimizerPage() {
  const [isSubmitting, setIsSubmitting] = useState(false)
  const [submitStatus, setSubmitStatus] = useState<"idle" | "success" | "error">("idle")
  const [submitMessage, setSubmitMessage] = useState("")
  const [showAdvancedOptions, setShowAdvancedOptions] = useState(false)

  const { dnaSequence, sequenceFile, wantedOrganisms, unwantedOrganisms, reset } = useOptimizationStore()

  const handleSubmit = async () => {
    setIsSubmitting(true)
    setSubmitStatus("idle")

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
        if (validation.errors) {
          setSubmitMessage(validation.errors.join(", "))
        } else {
          setSubmitMessage("Invalid submission data")
        }
        return
      }

      // Prepare JSON payload for backend
      const optimizationPayload = {
        dna_sequence: {
          content: sequenceFile ? null : dnaSequence,
          file_name: sequenceFile?.name || null,
          file_content: sequenceFile ? await fileToBase64(sequenceFile) : null,
        },
        wanted_organisms: wantedOrganisms.map((org) => ({
          id: org.id,
          genome_path: org.genomePath,
          priority: org.priority,
          expression_data_path: org.expressionDataPath || null,
        })),
        unwanted_organisms: unwantedOrganisms.map((org) => ({
          id: org.id,
          genome_path: org.genomePath,
          priority: org.priority,
          expression_data_path: org.expressionDataPath || null,
        })),
        advanced_options: {
          tuning_parameter: useOptimizationStore.getState().tuningParameter,
          optimization_method: useOptimizationStore.getState().optimizationMethod,
          cub_index: useOptimizationStore.getState().cubIndex,
        },
        timestamp: new Date().toISOString(),
      }

      console.log("Sending optimization request:", optimizationPayload)

      // Send POST request to backend
      const response = await fetch("http://localhost:8000/run_modules", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Accept: "application/json",
        },
        body: JSON.stringify(optimizationPayload),
      })

      if (!response.ok) {
        const errorData = await response.json().catch(() => ({ error: "Unknown error" }))
        throw new Error(errorData.error || `HTTP ${response.status}: ${response.statusText}`)
      }

      const result = await response.json()
      console.log("Optimization result:", result)

      setSubmitStatus("success")
      setSubmitMessage("DNA optimization completed successfully! Check the console for detailed results.")
    } catch (error) {
      console.error("Optimization error:", error)
      setSubmitStatus("error")

      if (error instanceof TypeError && error.message.includes("fetch")) {
        setSubmitMessage(
          "Unable to connect to the optimization server. Please ensure the backend is running on localhost:8000.",
        )
      } else {
        setSubmitMessage(`Optimization failed: ${error instanceof Error ? error.message : "Unknown error"}`)
      }
    } finally {
      setIsSubmitting(false)
    }
  }

  // Helper function to convert file to base64
  const fileToBase64 = (file: File): Promise<string> => {
    return new Promise((resolve, reject) => {
      const reader = new FileReader()
      reader.readAsDataURL(file)
      reader.onload = () => {
        const result = reader.result as string
        // Remove the data URL prefix (e.g., "data:text/plain;base64,")
        const base64 = result.split(",")[1]
        resolve(base64)
      }
      reader.onerror = (error) => reject(error)
    })
  }

  const handleReset = () => {
    reset()
    setSubmitStatus("idle")
    setSubmitMessage("")
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
          <h1 className="text-4xl font-bold text-gray-900">DNA Sequence Optimizer</h1>
          <p className="text-lg text-gray-600 max-w-2xl mx-auto">
            Configure your DNA sequence optimization parameters by specifying target organisms, priorities, and
            expression data for comprehensive genetic optimization.
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
              <strong>2. Wanted Organisms:</strong> Add organisms you want to optimize for with their genome paths and
              priorities
            </p>
            <p>
              <strong>3. Unwanted Organisms:</strong> Add organisms you want to avoid with their genome paths and
              priorities
            </p>
            <p>
              <strong>4. Priority Scores:</strong> Use values between 1-100 (higher = more important)
            </p>
            <p>
              <strong>5. Expression Data:</strong> Optionally provide CSV files with expression data for each organism
            </p>
          </CardContent>
        </Card>
      </div>
    </div>
  )
}
