const fs = require('fs')
const path = require('path')

const fsp = fs.promises

const projectRoot = process.cwd()
const nextDir = path.join(projectRoot, '.next')
const standaloneDir = path.join(nextDir, 'standalone')
const destNextDir = path.join(standaloneDir, '.next')

const dirsToCopy = ['server', 'static']
const filesToCopy = [
  'BUILD_ID',
  'app-build-manifest.json',
  'app-path-routes-manifest.json',
  'build-manifest.json',
  'images-manifest.json',
  'prerender-manifest.json',
  'react-loadable-manifest.json',
  'required-server-files.json',
  'routes-manifest.json',
]

async function pathExists(targetPath) {
  try {
    await fsp.access(targetPath)
    return true
  } catch {
    return false
  }
}

async function copyFile(src, dest) {
  await fsp.mkdir(path.dirname(dest), { recursive: true })
  await fsp.copyFile(src, dest)
  console.log(`Copied file: ${path.relative(projectRoot, dest)}`)
}

async function copyDir(src, dest) {
  await fsp.rm(dest, { recursive: true, force: true })
  await fsp.mkdir(path.dirname(dest), { recursive: true })
  await fsp.cp(src, dest, { recursive: true })
  console.log(`Copied directory: ${path.relative(projectRoot, dest)}`)
}

async function main() {
  const standaloneServer = path.join(standaloneDir, 'server.js')

  if (!(await pathExists(standaloneServer))) {
    throw new Error(
      'Standalone server not found. Run "npm run build" before copying assets.',
    )
  }

  await fsp.mkdir(destNextDir, { recursive: true })

  for (const dir of dirsToCopy) {
    const srcDir = path.join(nextDir, dir)
    const destDir = path.join(destNextDir, dir)

    if (await pathExists(srcDir)) {
      await copyDir(srcDir, destDir)
    } else {
      console.warn(`Skipping missing directory: ${path.relative(projectRoot, srcDir)}`)
    }
  }

  for (const file of filesToCopy) {
    const srcFile = path.join(nextDir, file)
    const destFile = path.join(destNextDir, file)

    if (await pathExists(srcFile)) {
      await copyFile(srcFile, destFile)
    } else {
      console.warn(`Skipping missing file: ${path.relative(projectRoot, srcFile)}`)
    }
  }

  console.log('Standalone .next assets copied successfully.')
}

main().catch((err) => {
  console.error(err.message)
  process.exit(1)
})

