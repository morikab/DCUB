// The preload script exposes these on `window`; both are absent in the web build.
/**
 * What the main process's "file-operations" IPC handler returns. It answers
 * several operations with different shapes, so the fields other than
 * `success` are optional (see electron/main.js):
 *   select-file, picked    -> { success: true, filePath, fileName, content }
 *   save-file, saved       -> { success: true, filePath }
 *   either one, cancelled  -> { success: false, canceled: true }
 *   unknown op or a throw  -> { success: false, error }
 * Callers must check `success` and then the field they need.
 */
export interface ElectronFileOperationResult {
  success: boolean
  canceled?: boolean
  error?: string
  filePath?: string
  fileName?: string
  content?: string
}

export interface ElectronAPI {
  isElectron?: boolean
  fileOperations?: (operation: string, data?: unknown) => Promise<ElectronFileOperationResult>
  makeRequest?: (url: string, options: RequestInit) => Promise<Response>
}

/** Electron's renderer adds `path` to File objects; browsers do not. */
export type ElectronFile = File & { path?: string }

declare global {
  interface Window {
    electronAPI?: ElectronAPI
  }
}

export const isElectron = (): boolean => {
  if (typeof window === "undefined") return false
  // Electron's renderer sets process.type. `window.process` resolves to
  // @types/node's `Process`, which does not declare `type` - reading it
  // directly is a TS2339 error - so narrow to just the field being tested
  // rather than widening the whole object to `any`.
  const rendererProcess = (window as { process?: { type?: string } }).process
  return rendererProcess?.type === "renderer"
}

export const getElectronAPI = (): ElectronAPI | null => {
  if (typeof window !== "undefined" && window.electronAPI) {
    return window.electronAPI
  }
  return null
}

// Helper to handle file operations in Electron
export const handleElectronFileOperation = async (operation: string, data?: unknown) => {
  const electronAPI = getElectronAPI()
  if (electronAPI && electronAPI.fileOperations) {
    try {
      return await electronAPI.fileOperations(operation, data)
    } catch (error) {
      console.error("Electron file operation failed:", error)
      return null
    }
  }
  return null
}

// Network request helper for Electron
export const makeElectronRequest = async (url: string, options: RequestInit) => {
  const electronAPI = getElectronAPI()

  // If running in Electron, use the main process for network requests
  if (electronAPI && electronAPI.makeRequest) {
    try {
      return await electronAPI.makeRequest(url, options)
    } catch (error) {
      console.error("Electron request failed:", error)
      throw error
    }
  }

  // Fallback to regular fetch for web version
  return fetch(url, options)
}
