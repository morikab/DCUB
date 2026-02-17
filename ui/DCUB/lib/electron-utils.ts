export const isElectron = (): boolean => {
  return typeof window !== "undefined" && window.process && window.process.type === "renderer"
}

export const getElectronAPI = () => {
  if (typeof window !== "undefined" && (window as any).electronAPI) {
    return (window as any).electronAPI
  }
  return null
}

// Helper to handle file operations in Electron
export const handleElectronFileOperation = async (operation: string, data?: any) => {
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
