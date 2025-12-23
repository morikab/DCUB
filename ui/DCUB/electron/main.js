const { app, BrowserWindow, ipcMain, dialog } = require("electron")
const path = require("path")
const { spawn, exec } = require('child_process');
const http = require('http');
const fs = require("fs").promises
const fsSync = require('fs');
const fetch = require("node-fetch")

let mainWindow;
let nextProcess;
let backendProcess;

function trace(msg) {
  fsSync.appendFileSync("/tmp/electron-trace.txt", `[${Date.now()}] ${msg}\n`);
}

// Function to kill processes on port 3000
function killProcessOnPort(port) {
  return new Promise((resolve) => {
    exec(`lsof -ti:${port}`, (error, stdout) => {
      if (stdout.trim()) {
        const pids = stdout.trim().split('\n');
        console.log(`Killing processes on port ${port}: ${pids.join(', ')}`);
        
        pids.forEach(pid => {
          exec(`kill -9 ${pid}`, (killError) => {
            if (killError) {
              console.log(`Could not kill process ${pid}:`, killError.message);
            } else {
              console.log(`Killed process ${pid}`);
            }
          });
        });
        
        // Wait a bit for processes to be killed
        setTimeout(resolve, 1000);
      } else {
        resolve();
      }
    });
  });
}

function waitForServer(url, timeout = 15000) {
  return new Promise((resolve, reject) => {
    // Force IPv4 to avoid localhost resolving to ::1 on some systems
    if (url.includes('localhost')) {
      url = url.replace('localhost', '127.0.0.1')
    }

    const start = Date.now();
    let attempts = 0;
    let currentRequest = null;
    let timeoutId = null;
    let isResolved = false;

    const cleanup = () => {
      if (currentRequest) {
        currentRequest.removeAllListeners();
        currentRequest.destroy();
        currentRequest = null;
      }
      if (timeoutId) {
        clearTimeout(timeoutId);
        timeoutId = null;
      }
    };

    const check = () => {
      if (isResolved) return;
      
      attempts++;
      
      // Cleanup any previous request
      cleanup();
      
      try {
        currentRequest = http.get(url, (res) => {
          if (isResolved) return;
          isResolved = true;
          cleanup();
          resolve();
        });
        
        currentRequest.on('error', (err) => {
          if (isResolved) return;
          
          cleanup();
          
          const elapsed = Date.now() - start;
          
          if (elapsed > timeout) {
            isResolved = true;
            reject(new Error(`Server did not start within ${timeout}ms after ${attempts} attempts: ${err.message}`));
            return;
          }
          
          // Retry after 1 second
          timeoutId = setTimeout(() => {
            if (!isResolved) {
              check();
            }
          }, 1000);
        });
        
        // Set request timeout
        currentRequest.setTimeout(2000, () => {
          if (isResolved || !currentRequest) return;
          
          cleanup();
          
          const elapsed = Date.now() - start;
          
          if (elapsed > timeout) {
            isResolved = true;
            reject(new Error(`Server did not start within ${timeout}ms after ${attempts} attempts: request timeout`));
            return;
          }
          
          // Retry after 1 second
          timeoutId = setTimeout(() => {
            if (!isResolved) {
              check();
            }
          }, 1000);
        });
        
        // Handle abort/close
        currentRequest.on('close', () => {
          if (isResolved || !currentRequest) return;
          cleanup();
        });
        
      } catch (err) {
        cleanup();
        if (!isResolved) {
          isResolved = true;
          reject(new Error(`Failed to check server: ${err.message}`));
        }
      }
    };

    check();
  });
}

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
      webSecurity: false,  // TODO - Allow local files
    },
    icon: path.join(__dirname, "assets", "icon.png"), // Add your app icon
    titleBarStyle: "default",
    show: false,
  })

  mainWindow.loadFile(
    path.join(__dirname, "loading.html")
  );

  mainWindow.once("ready-to-show", () => {
    mainWindow.show();
    mainWindow.focus();
  })

  mainWindow.on("closed", () => {
    shutdown();
    mainWindow = null
  })
}

// Kill child processes when Electron quits
const shutdown = () => {
  if (nextProcess) {
    nextProcess.kill();
    nextProcess = null;
  }
  if (backendProcess) {
    backendProcess.kill();
    backendProcess = null;
  }
};

// Start FastAPI backend server as a child process
async function startBackendServer() {
  try {
    // Path to packaged FastAPI backend executable, relative to this file
    const isDev = !app.isPackaged;

    const basePath = isDev
      ? path.join(__dirname, "../backend")
      : process.resourcesPath;

    const backendExecutable = path.join(basePath, "fastapi_server", "fastapi_server");
    
    // FIXMW
    // switch (process.platform) {
    //   case "darwin": return path.join(basePath, "fastapi_server_macos");
    //   case "win32":  return path.join(basePath, "fastapi_server_win.exe");
    //   case "linux":  return path.join(basePath, "fastapi_server_linux");
    // }
    

    console.log("Starting FastAPI backend executable...", { backendExecutable });

    backendProcess = spawn(backendExecutable, [], {
      cwd: path.dirname(backendExecutable),
      shell: false,
      stdio: ["ignore", "pipe", "pipe"],
      env: {
        ...process.env,
      },
    });

    backendProcess.stdout.on("data", (data) => {
      console.log(`FastAPI stdout: ${data.toString().trim()}`);
    });

    backendProcess.stderr.on("data", (data) => {
      console.error(`FastAPI stderr: ${data.toString().trim()}`);
    });

    backendProcess.on("error", (err) => {
      console.error("Failed to start FastAPI backend executable:", err);
    });

    backendProcess.on("exit", (code) => {
      console.log(`FastAPI backend executable exited with code ${code}`);
    });
  } catch (err) {
    console.error("Error while starting FastAPI backend executable:", err);
  }
}

// app.whenReady().then(createWindow)
app.on('ready', async () => {
  trace("app ready event fired");
  try {
    const projectRoot = path.join(__dirname, '..');
    const distDir = path.join(projectRoot, '.next');
    const standaloneDir = path.join(distDir, 'standalone');
    const embeddedNextDir = path.join(standaloneDir, '.next');
    const standaloneServer = path.join(standaloneDir, 'server.js');

    if (!fsSync.existsSync(standaloneServer)) {
      console.error('Standalone server not found. Please run "npm run build" first.');
      app.quit();
      return;
    }

    if (!fsSync.existsSync(embeddedNextDir)) {
      console.error('Missing ".next/standalone/.next" assets. Run "npm run copy-standalone-next" first.');
      app.quit();
      return;
    }

    createWindow();
    trace("window created");

    console.log('Starting Next.js standalone server...');
    // Kill any existing processes on port 3000
    await killProcessOnPort(3000);
    
    // Start the standalone server using Electron's Node.js runtime
    // Use absolute paths to avoid PATH issues when launched from Finder
    const absoluteStandaloneServer = path.resolve(standaloneServer);
    const absoluteStandaloneDir = path.resolve(standaloneDir);
    
    // Build environment variables - explicitly remove NODE_OPTIONS in packaged apps
    const env = { ...process.env };
    
    // Remove NODE_OPTIONS in packaged apps (Electron's Node.js doesn't support it)
    if (app.isPackaged) {
      delete env.NODE_OPTIONS;
    } else {
      // Only set NODE_OPTIONS in dev mode
      env.NODE_OPTIONS = '--dns-result-order=ipv4first';
    }
    
    env.PORT = '3000';
    env.HOSTNAME = '127.0.0.1';  // Force IPv4
    
    // Ensure PATH includes common directories for Finder launches
    // if (app.isPackaged && process.platform === 'darwin') {
    //   const commonPaths = ['/usr/local/bin', '/usr/bin', '/bin', '/opt/homebrew/bin'];
    //   const currentPath = env.PATH || '';
    //   env.PATH = [...commonPaths, currentPath].join(':');
    // }
    
    console.log('Spawning Next.js server:', { absoluteStandaloneServer, absoluteStandaloneDir });
    
    nextProcess = spawn(process.execPath, [absoluteStandaloneServer], {
      cwd: absoluteStandaloneDir,
      stdio: ['ignore', 'pipe', 'pipe'],
      env: {
        ...env,
        ELECTRON_RUN_AS_NODE: '1',
      }
    });
    trace("next spawn called");

    // Handle process output
    nextProcess.stdout.on('data', (data) => {
      trace("next stdout first chunk");
      const output = data.toString().trim();
      console.log(`Next.js stdout: ${output}`);
      
      // Check for server ready message
      if (output.includes('Ready') || output.includes('started server') || output.includes('Local:')) {
        console.log('Server appears to be ready based on output');
      }
    });

    nextProcess.stderr.on('data', (data) => {
      const output = data.toString().trim();
      console.error(`Next.js stderr: ${output}`);
      
      // Check for specific errors
      if (output.includes('EADDRINUSE')) {
        console.error('Port 3000 is already in use!');
      } else if (output.includes('ECONNREFUSED')) {
        console.error('Connection refused - server may not be starting properly');
      }
    });

    nextProcess.on('error', (err) => {
      console.error('Failed to start Next.js server:', err);
      app.quit();
    });

    nextProcess.on('exit', (code) => {
      if (code !== 0) {
        console.error(`Next.js server exited with code ${code}`);
        fsSync.writeFileSync('/tmp/electron-window-started.txt', `Next.js server exited with code ${code}`);
        app.quit();
      }
    });
    trace("waitForServer started");
    // Wait for server to be ready
    await waitForServer('http://127.0.0.1:3000');
    trace("waitForServer resolved");
    mainWindow.loadURL("http://127.0.0.1:3000");
    trace("loadURL called");
    // Start Python FastAPI backend in parallel
    await startBackendServer();
    
  } catch (err) {
    console.error('Failed to start application:', err);
    app.quit();
  }
});

app.on("window-all-closed", () => {
  shutdown();
  if (process.platform !== "darwin") {
    app.quit()
  }
})

app.on("activate", () => {
  if (BrowserWindow.getAllWindows().length === 0) {
    createWindow()
  } else {
    mainWindow.show();
    mainWindow.focus();
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
