"use client"

import type React from "react"

import { useState } from "react"
import { Button } from "@/components/ui/button"
import { Upload, Download } from "lucide-react"
import { isElectron, handleElectronFileOperation } from "@/lib/electron-utils"

interface ElectronFileHandlerProps {
  onFileSelect?: (file: { name: string; content: string; path: string }) => void
  onFileSave?: (content: string, defaultName?: string) => void
  acceptedTypes?: string[]
  children?: React.ReactNode
}

export function ElectronFileHandler({
  onFileSelect,
  onFileSave,
  acceptedTypes = ["*"],
  children,
}: ElectronFileHandlerProps) {
  const [isLoading, setIsLoading] = useState(false)

  const handleFileSelect = async () => {
    if (!isElectron()) {
      console.warn("File selection only available in Electron app")
      return
    }

    setIsLoading(true)
    try {
      const filters = acceptedTypes.includes("*")
        ? [{ name: "All Files", extensions: ["*"] }]
        : [
            { name: "FASTA Files", extensions: ["fasta", "fa"] },
            { name: "CSV Files", extensions: ["csv"] },
            { name: "GenBank Files", extensions: ["gb", "gbk"] },
            { name: "All Files", extensions: ["*"] },
          ]

      const result = await handleElectronFileOperation("select-file", { filters })

      if (result?.success && !result.canceled) {
        onFileSelect?.({
          name: result.fileName,
          content: result.content,
          path: result.filePath,
        })
      }
    } catch (error) {
      console.error("File selection failed:", error)
    } finally {
      setIsLoading(false)
    }
  }

  const handleFileSave = async (content: string, defaultName = "file.txt") => {
    if (!isElectron()) {
      console.warn("File saving only available in Electron app")
      return
    }

    setIsLoading(true)
    try {
      const result = await handleElectronFileOperation("save-file", {
        content,
        defaultName,
        filters: [
          { name: "JSON Files", extensions: ["json"] },
          { name: "Text Files", extensions: ["txt"] },
          { name: "All Files", extensions: ["*"] },
        ],
      })

      if (result?.success) {
        console.log("File saved successfully:", result.filePath)
      }
    } catch (error) {
      console.error("File saving failed:", error)
    } finally {
      setIsLoading(false)
    }
  }

  if (!isElectron()) {
    return children || null
  }

  return (
    <div className="flex gap-2">
      {onFileSelect && (
        <Button onClick={handleFileSelect} disabled={isLoading} variant="outline" size="sm">
          <Upload className="w-4 h-4 mr-2" />
          {isLoading ? "Loading..." : "Select File"}
        </Button>
      )}

      {onFileSave && (
        <Button
          onClick={() => handleFileSave("", "optimization-config.json")}
          disabled={isLoading}
          variant="outline"
          size="sm"
        >
          <Download className="w-4 h-4 mr-2" />
          Save Config
        </Button>
      )}

      {children}
    </div>
  )
}
