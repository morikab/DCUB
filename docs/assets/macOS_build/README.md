## macOS Installation Note

If macOS reports that the app is damaged and cannot be opened, run the following command in Terminal:

```bash
xattr -d com.apple.quarantine /Applications/DCUB.app
```

After running this command, try opening the app again.
