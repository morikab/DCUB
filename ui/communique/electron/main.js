const { app, BrowserWindow, ipcMain, dialog } = require("electron")
const path = require("path")
const fs = require("fs").promises
const fetch = require("node-fetch")

let mainWindow

function createWindow() {
  mainWindow = new BrowserWindow({
    width: 1400,
    height: 900,
    minWidth: 800,
    minHeight: 600,
    webPreferences: {
      nodeIntegration: false,
      contextIsolation: true,
      enableRemoteModule: false,
      preload: path.join(__dirname, "preload.js"),
    },
    icon: path.join(__dirname, "assets", "icon.png"), // Add your app icon
    titleBarStyle: "default",
    show: false,
  })

  // Load the Next.js app
  const isDev = process.env.NODE_ENV === "development"

  if (isDev) {
    mainWindow.loadURL("http://localhost:3000")
    mainWindow.webContents.openDevTools()
  } else {
    mainWindow.loadFile(path.join(__dirname, "../out/index.html"))
  }

  mainWindow.once("ready-to-show", () => {
    mainWindow.show()
  })

  mainWindow.on("closed", () => {
    mainWindow = null
  })
}

app.whenReady().then(createWindow)

app.on("window-all-closed", () => {
  if (process.platform !== "darwin") {
    app.quit()
  }
})

app.on("activate", () => {
  if (BrowserWindow.getAllWindows().length === 0) {
    createWindow()
  }
})

// IPC handlers for file operations
ipcMain.handle("file-operations", async (event, operation, data) => {
  try {
    switch (operation) {
      case "select-file":
        const result = await dialog.showOpenDialog(mainWindow, {
          properties: ["openFile"],
          filters: data.filters || [{ name: "All Files", extensions: ["*"] }],
        })

        if (!result.canceled && result.filePaths.length > 0) {
          const filePath = result.filePaths[0]
          const fileContent = await fs.readFile(filePath, "utf-8")
          return {
            success: true,
            filePath,
            fileName: path.basename(filePath),
            content: fileContent,
          }
        }
        return { success: false, canceled: true }

      case "save-file":
        const saveResult = await dialog.showSaveDialog(mainWindow, {
          defaultPath: data.defaultName || "optimization-results.json",
          filters: data.filters || [
            { name: "JSON Files", extensions: ["json"] },
            { name: "All Files", extensions: ["*"] },
          ],
        })

        if (!saveResult.canceled) {
          await fs.writeFile(saveResult.filePath, data.content, "utf-8")
          return { success: true, filePath: saveResult.filePath }
        }
        return { success: false, canceled: true }

      default:
        return { success: false, error: "Unknown operation" }
    }
  } catch (error) {
    console.error("File operation error:", error)
    return { success: false, error: error.message }
  }
})

// IPC handler for network requests
ipcMain.handle("make-request", async (event, url, options) => {
  try {
    const response = await fetch(url, {
      ...options,
      headers: {
        "Content-Type": "application/json",
        Accept: "application/json",
        ...options.headers,
      },
    })

    const data = await response.json()

    return {
      ok: response.ok,
      status: response.status,
      statusText: response.statusText,
      data: data,
    }
  } catch (error) {
    console.error("Network request error:", error)
    throw error
  }
})

// Handle app updates and other Electron-specific features
ipcMain.handle("app-info", async () => {
  return {
    version: app.getVersion(),
    name: app.getName(),
    platform: process.platform,
    arch: process.arch,
  }
})
