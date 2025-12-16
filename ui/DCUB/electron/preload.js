const { contextBridge, ipcRenderer } = require("electron")

// Expose protected methods that allow the renderer process to use
// the ipcRenderer without exposing the entire object
contextBridge.exposeInMainWorld("electronAPI", {
  // File operations
  fileOperations: (operation, data) => ipcRenderer.invoke("file-operations", operation, data),

  // Network requests
  makeRequest: (url, options) => ipcRenderer.invoke("make-request", url, options),

  // App information
  getAppInfo: () => ipcRenderer.invoke("app-info"),

  // Platform detection
  platform: process.platform,
  isElectron: true,
})

// Add window controls for custom title bar if needed
contextBridge.exposeInMainWorld("windowControls", {
  minimize: () => ipcRenderer.invoke("window-minimize"),
  maximize: () => ipcRenderer.invoke("window-maximize"),
  close: () => ipcRenderer.invoke("window-close"),
})
