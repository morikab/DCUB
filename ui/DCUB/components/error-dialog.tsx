"use client"

import { useState } from "react"
import { AlertCircle, ChevronDown, ChevronUp } from "lucide-react"
import { Dialog, DialogContent, DialogFooter, DialogHeader, DialogTitle } from "@/components/ui/dialog"
import { Button } from "@/components/ui/button"

function extractHeadline(errorText: string): string {
  const lines = errorText.trim().split("\n").filter((line) => line.trim())
  return lines.length > 0 ? lines[lines.length - 1].trim() : errorText
}

interface ErrorDialogProps {
  open: boolean
  errorText: string
  onClose: () => void
  title?: string
}

export function ErrorDialog({ open, errorText, onClose, title = "DNA optimization failed" }: ErrorDialogProps) {
  const [showDetails, setShowDetails] = useState(false)

  return (
    <Dialog
      open={open}
      onOpenChange={(nextOpen) => {
        if (!nextOpen) {
          setShowDetails(false)
          onClose()
        }
      }}
    >
      <DialogContent>
        <DialogHeader>
          <DialogTitle className="flex items-center gap-2 text-destructive">
            <AlertCircle className="w-5 h-5" />
            {title}
          </DialogTitle>
        </DialogHeader>

        <p className="text-sm text-foreground">{extractHeadline(errorText)}</p>

        <Button
          variant="ghost"
          size="sm"
          className="w-fit -ml-2 text-muted-foreground"
          onClick={() => setShowDetails((prev) => !prev)}
        >
          {showDetails ? <ChevronUp className="w-4 h-4" /> : <ChevronDown className="w-4 h-4" />}
          {showDetails ? "Hide technical details" : "Show technical details"}
        </Button>

        {showDetails && (
          <pre className="max-h-64 overflow-auto rounded-md border bg-muted p-3 text-xs whitespace-pre-wrap break-all">
            {errorText}
          </pre>
        )}

        <DialogFooter>
          <Button onClick={onClose}>Close</Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
