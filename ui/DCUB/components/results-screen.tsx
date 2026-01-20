"use client"

import { useState } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { Button } from "@/components/ui/button"
import { Badge } from "@/components/ui/badge"
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs"
import { ArrowLeft, Download, Copy, CheckCircle, BarChart3, Dna, Clock, RefreshCw } from "lucide-react"
import type { OptimizationResult } from "@/lib/types"

interface ResultsScreenProps {
  result: OptimizationResult
  onBackToForm: () => void
  onReset: () => void
}

export function ResultsScreen({ result, onBackToForm, onReset }: ResultsScreenProps) {
  const [copiedSequence, setCopiedSequence] = useState(false)
  const [copiedScores, setCopiedScores] = useState(false)

  const copyToClipboard = async (text: string, type: "sequence" | "scores") => {
    try {
      await navigator.clipboard.writeText(text)
      if (type === "sequence") {
        setCopiedSequence(true)
        setTimeout(() => setCopiedSequence(false), 2000)
      } else {
        setCopiedScores(true)
        setTimeout(() => setCopiedScores(false), 2000)
      }
    } catch (error) {
      console.error("Failed to copy to clipboard:", error)
    }
  }

  const downloadResults = () => {
    const resultsData = {
      optimized_sequence: result.optimized_sequence,
      evaluation_scores: result.evaluation_scores,
      optimization_parameters: result.optimization_parameters,
      processing_time: result.processing_time,
      timestamp: result.timestamp,
    }

    const blob = new Blob([JSON.stringify(resultsData, null, 2)], { type: "application/json" })
    const url = URL.createObjectURL(blob)
    const a = document.createElement("a")
    a.href = url
    a.download = `commuique_results_${new Date().toISOString().split("T")[0]}.json`
    document.body.appendChild(a)
    a.click()
    document.body.removeChild(a)
    URL.revokeObjectURL(url)
  }

  const downloadFasta = () => {
    const fastaContent = `>Optimized_Sequence_${new Date().toISOString().split("T")[0]}
${result.optimized_sequence}`

    const blob = new Blob([fastaContent], { type: "text/plain" })
    const url = URL.createObjectURL(blob)
    const a = document.createElement("a")
    a.href = url
    a.download = `optimized_sequence_${new Date().toISOString().split("T")[0]}.fasta`
    document.body.appendChild(a)
    a.click()
    document.body.removeChild(a)
    URL.revokeObjectURL(url)
  }

  const formatSequence = (sequence: string, lineLength = 60) => {
    return sequence.match(new RegExp(`.{1,${lineLength}}`, "g"))?.join("\n") || sequence
  }

  const getScoreColor = (score: number, type: "average_distance" | "ratio" | "weakest_link") => {
    if (type === "average_distance") {
      if (score >= 0.8) return "text-green-600"
      if (score >= 0.6) return "text-yellow-600"
      return "text-red-600"
    }
    if (type === "ratio") {
      if (score >= 40 && score <= 60) return "text-green-600"
      if (score >= 30 && score <= 70) return "text-yellow-600"
      return "text-red-600"
    }
    if (type === "weakest_link") {
      if (score >= 0.7) return "text-green-600"
      if (score >= 0.5) return "text-yellow-600"
      return "text-red-600"
    }
    return "text-gray-600"
  }

  const getScoreBadgeVariant = (score: number, type: "average_distance" | "ratio" | "weakest_link") => {
    const color = getScoreColor(score, type)
    if (color.includes("green")) return "default"
    if (color.includes("yellow")) return "secondary"
    return "destructive"
  }

  return (
    <div className="min-h-screen bg-gradient-to-br from-green-50 to-blue-100 p-4">
      <div className="max-w-6xl mx-auto space-y-6">
        {/* Header */}
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-4">
            <Button onClick={onBackToForm} variant="outline" className="flex items-center gap-2 bg-transparent">
              <ArrowLeft className="w-4 h-4" />
              Back to Form
            </Button>
            <div>
              <h1 className="text-3xl font-bold text-gray-900">Optimization Results</h1>
              <p className="text-gray-600">Your DNA sequence has been successfully optimized</p>
            </div>
          </div>
          <div className="flex gap-2">
            <Button onClick={downloadResults} variant="outline" className="flex items-center gap-2 bg-transparent">
              <Download className="w-4 h-4" />
              Download JSON
            </Button>
            <Button onClick={onReset} variant="outline" className="flex items-center gap-2 bg-transparent">
              <RefreshCw className="w-4 h-4" />
              New Optimization
            </Button>
          </div>
        </div>

        {/* Success Banner */}
        <Card className="bg-green-50 border-green-200">
          <CardContent className="pt-6">
            <div className="flex items-center gap-3">
              <CheckCircle className="w-8 h-8 text-green-600" />
              <div>
                <h3 className="text-lg font-semibold text-green-900">Optimization Complete!</h3>
                <p className="text-green-700">
                  Processing completed in {result.processing_time ? `${result.processing_time}s` : "unknown time"}
                </p>
              </div>
            </div>
          </CardContent>
        </Card>

        {/* Main Results */}
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Evaluation Scores */}
          <Card className="lg:col-span-1">
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <BarChart3 className="w-5 h-5" />
                Evaluation Scores
              </CardTitle>
              <CardDescription>Quality metrics for the optimized sequence</CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <div className="space-y-3">
                <div className="flex items-center justify-between p-3 bg-gray-50 rounded-lg">
                  <div>
                    <p className="font-medium text-gray-900">Average Distance Score</p>
                  </div>
                  <Badge variant={getScoreBadgeVariant(result.evaluation_scores.average_distance_score, "average_distance")}>
                    {result.evaluation_scores.average_distance_score.toFixed(3)}
                  </Badge>
                </div>

                <div className="flex items-center justify-between p-3 bg-gray-50 rounded-lg">
                  <div>
                    <p className="font-medium text-gray-900">Ratio Score</p>
                    {/*<p className="text-xs text-gray-600">Ratio of gemetric means between organism groups</p>*/}
                  </div>
                  <Badge variant={getScoreBadgeVariant(result.evaluation_scores.ratio_score, "ratio")}>
                    {result.evaluation_scores.ratio_score.toFixed(3)}
                  </Badge>
                </div>

                <div className="flex items-center justify-between p-3 bg-gray-50 rounded-lg">
                  <div>
                    <p className="font-medium text-gray-900">Weakest Link Score</p>
                    {/*<p className="text-xs text-gray-600">Minimal differetial distance between organism groups</p>*/}
                  </div>
                  <Badge variant={getScoreBadgeVariant(result.evaluation_scores.weakest_link_score, "weakest_link")}>
                    {result.evaluation_scores.weakest_link_score.toFixed(3)}
                  </Badge>
                </div>
              </div>

              <Button
                onClick={() => copyToClipboard(JSON.stringify(result.evaluation_scores, null, 2), "scores")}
                variant="outline"
                className="w-full"
              >
                {copiedScores ? (
                  <>
                    <CheckCircle className="w-4 h-4 mr-2" />
                    Copied!
                  </>
                ) : (
                  <>
                    <Copy className="w-4 h-4 mr-2" />
                    Copy Scores
                  </>
                )}
              </Button>
            </CardContent>
          </Card>

          {/* Optimized Sequence */}
          <Card className="lg:col-span-2">
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Dna className="w-5 h-5" />
                Optimized Sequence
              </CardTitle>
              <CardDescription>Length: {result.optimized_sequence.length} nucleotides</CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              <Tabs defaultValue="formatted" className="w-full">
                <TabsList className="grid w-full grid-cols-2">
                  <TabsTrigger value="formatted">Formatted</TabsTrigger>
                  <TabsTrigger value="raw">Raw Sequence</TabsTrigger>
                </TabsList>

                <TabsContent value="formatted" className="mt-4">
                  <div className="bg-gray-50 p-4 rounded-lg border">
                    <pre className="text-sm font-mono text-gray-800 whitespace-pre-wrap break-all">
                      {formatSequence(result.optimized_sequence)}
                    </pre>
                  </div>
                </TabsContent>

                <TabsContent value="raw" className="mt-4">
                  <div className="bg-gray-50 p-4 rounded-lg border">
                    <p className="text-sm font-mono text-gray-800 break-all">{result.optimized_sequence}</p>
                  </div>
                </TabsContent>
              </Tabs>

              <div className="flex gap-2">
                <Button
                  onClick={() => copyToClipboard(result.optimized_sequence, "sequence")}
                  variant="outline"
                  className="flex-1"
                >
                  {copiedSequence ? (
                    <>
                      <CheckCircle className="w-4 h-4 mr-2" />
                      Copied!
                    </>
                  ) : (
                    <>
                      <Copy className="w-4 h-4 mr-2" />
                      Copy Sequence
                    </>
                  )}
                </Button>
                <Button onClick={downloadFasta} variant="outline" className="flex-1 bg-transparent">
                  <Download className="w-4 h-4 mr-2" />
                  Download FASTA
                </Button>
              </div>
            </CardContent>
          </Card>
        </div>

        {/* Optimization Parameters */}
        <Card>
          <CardHeader>
            <CardTitle className="flex items-center gap-2">
              <Clock className="w-5 h-5" />
              Optimization Details
            </CardTitle>
            <CardDescription>Parameters used for this optimization</CardDescription>
          </CardHeader>
          <CardContent>
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
              <div className="bg-blue-50 p-3 rounded-lg">
                <p className="text-sm font-medium text-blue-900">Tuning Parameter</p>
                <p className="text-lg font-bold text-blue-700">{result.optimization_parameters.tuning_parameter}</p>
              </div>
              <div className="bg-purple-50 p-3 rounded-lg">
                <p className="text-sm font-medium text-purple-900">Method</p>
                <p className="text-sm font-bold text-purple-700">
                  {result.optimization_parameters.optimization_method.replace(/_/g, " ").toUpperCase()}
                </p>
              </div>
              <div className="bg-green-50 p-3 rounded-lg">
                <p className="text-sm font-medium text-green-900">CUB Index</p>
                <p className="text-lg font-bold text-green-700">{result.optimization_parameters.cub_index}</p>
              </div>
              <div className="bg-orange-50 p-3 rounded-lg">
                <p className="text-sm font-medium text-orange-900">Completed</p>
                <p className="text-sm font-bold text-orange-700">{new Date(result.timestamp).toLocaleString()}</p>
              </div>
            </div>
          </CardContent>
        </Card>
      </div>
    </div>
  )
}
