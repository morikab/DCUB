"use client"

import { useEffect, useState } from "react"
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card"
import { Progress } from "@/components/ui/progress"
import { Loader2, Dna, Cpu, BarChart3 } from "lucide-react"

export function LoadingScreen() {
  const [progress, setProgress] = useState(0)
  const [currentStep, setCurrentStep] = useState(0)
  const [elapsedTime, setElapsedTime] = useState(0)

  const steps = [
    { icon: Dna, label: "Analyzing DNA sequence", description: "Parsing and validating input sequence" },
    { icon: Cpu, label: "Processing organisms", description: "Loading GenBank files and expression data" },
    { icon: BarChart3, label: "Optimizing codons", description: "Calculating optimal codon usage patterns" },
    { icon: Loader2, label: "Finalizing results", description: "Generating evaluation scores and output" },
  ]

  useEffect(() => {
    const startTime = Date.now()

    // Update elapsed time every second
    const timeInterval = setInterval(() => {
      setElapsedTime(Math.floor((Date.now() - startTime) / 1000))
    }, 1000)

    // Simulate progress through steps
    const progressInterval = setInterval(() => {
      setProgress((prev) => {
        const newProgress = prev + Math.random() * 2

        // Update current step based on progress
        if (newProgress < 25) setCurrentStep(0)
        else if (newProgress < 50) setCurrentStep(1)
        else if (newProgress < 75) setCurrentStep(2)
        else setCurrentStep(3)

        return Math.min(newProgress, 95) // Don't reach 100% until actual completion
      })
    }, 1000)

    return () => {
      clearInterval(progressInterval)
      clearInterval(timeInterval)
    }
  }, [])

  const formatTime = (seconds: number) => {
    const mins = Math.floor(seconds / 60)
    const secs = seconds % 60
    return `${mins}:${secs.toString().padStart(2, "0")}`
  }

  const CurrentIcon = steps[currentStep].icon

  return (
    <div className="min-h-screen bg-gradient-to-br from-blue-50 to-indigo-100 flex items-center justify-center p-4">
      <Card className="w-full max-w-2xl shadow-2xl">
        <CardHeader className="text-center pb-6">
          <CardTitle className="text-3xl font-bold text-gray-900 mb-2">Optimizing DNA Sequence</CardTitle>
          <p className="text-gray-600">
            Please wait while we process your optimization request. This may take several minutes.
          </p>
        </CardHeader>
        <CardContent className="space-y-8">
          {/* Main Loading Animation */}
          <div className="flex flex-col items-center space-y-4">
            <div className="relative">
              <div className="w-24 h-24 border-4 border-blue-200 rounded-full animate-spin border-t-blue-600"></div>
              <div className="absolute inset-0 flex items-center justify-center">
                <CurrentIcon className="w-8 h-8 text-blue-600 animate-pulse" />
              </div>
            </div>

            <div className="text-center">
              <h3 className="text-lg font-semibold text-gray-800 mb-1">{steps[currentStep].label}</h3>
              <p className="text-sm text-gray-600">{steps[currentStep].description}</p>
            </div>
          </div>

          {/* Progress Bar */}
          <div className="space-y-2">
            <div className="flex justify-between text-sm text-gray-600">
              <span>Progress</span>
              <span>{Math.round(progress)}%</span>
            </div>
            <Progress value={progress} className="w-full h-2" />
          </div>

          {/* Processing Steps */}
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
            {steps.map((step, index) => {
              const StepIcon = step.icon
              const isActive = index === currentStep
              const isCompleted = index < currentStep

              return (
                <div
                  key={index}
                  className={`flex flex-col items-center p-3 rounded-lg transition-all ${
                    isActive
                      ? "bg-blue-100 border-2 border-blue-300"
                      : isCompleted
                        ? "bg-green-100 border-2 border-green-300"
                        : "bg-gray-100 border-2 border-gray-200"
                  }`}
                >
                  <StepIcon
                    className={`w-6 h-6 mb-2 ${
                      isActive ? "text-blue-600 animate-pulse" : isCompleted ? "text-green-600" : "text-gray-400"
                    }`}
                  />
                  <span
                    className={`text-xs text-center font-medium ${
                      isActive ? "text-blue-800" : isCompleted ? "text-green-800" : "text-gray-600"
                    }`}
                  >
                    {step.label}
                  </span>
                </div>
              )
            })}
          </div>

          {/* Status Information */}
          <div className="bg-gray-50 rounded-lg p-4 space-y-2">
            <div className="flex justify-between items-center">
              <span className="text-sm font-medium text-gray-700">Elapsed Time:</span>
              <span className="text-sm font-mono text-gray-900">{formatTime(elapsedTime)}</span>
            </div>
            <div className="flex justify-between items-center">
              <span className="text-sm font-medium text-gray-700">Status:</span>
              <span className="text-sm text-blue-600 font-medium">Processing...</span>
            </div>
            <div className="flex justify-between items-center">
              <span className="text-sm font-medium text-gray-700">Server:</span>
              <span className="text-sm text-green-600 font-medium">Connected</span>
            </div>
          </div>

          {/* Tips */}
          <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
            <h4 className="text-sm font-semibold text-blue-900 mb-2">💡 Processing Tips</h4>
            <ul className="text-xs text-blue-800 space-y-1">
              <li>• Complex sequences may take longer to optimize</li>
              <li>• Multiple organisms increase processing time</li>
              <li>• Large GenBank files require additional parsing time</li>
              <li>• Please keep this tab open during processing</li>
            </ul>
          </div>
        </CardContent>
      </Card>
    </div>
  )
}
