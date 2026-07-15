# DCUB vX.Y.Z

<!-- Replace vX.Y.Z with the release version and fill in the changelog below -->

## What's new

- 

---

## Download

| Platform | File |
|----------|------|
| macOS (Apple Silicon) | `DCUB-X.Y.Z-arm64.dmg` |
| Windows | `DCUB.Setup.X.Y.Z.exe` |
| Linux | `DCUB-X.Y.Z.AppImage` |

---

## Installation notes

### macOS — "App is damaged and can't be opened"

macOS quarantine flags unsigned apps downloaded from the internet. After moving the app to `/Applications`, run once in Terminal:

```bash
xattr -d com.apple.quarantine /Applications/DCUB.app
```

Then open the app normally.

### Windows — SmartScreen warning

If Windows Defender SmartScreen blocks the installer, click **More info** → **Run anyway**.

### Linux — AppImage won't launch

Make the file executable before running:

```bash
chmod +x DCUB-*.AppImage && ./DCUB-*.AppImage
```
