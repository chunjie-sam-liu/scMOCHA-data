---
name: pptx-to-pdf-powerpoint
description: Convert PowerPoint .pptx/.ppt files to PDFs on macOS by driving Microsoft PowerPoint via AppleScript. Use when the user asks for batch conversion, directory-based conversion, or a CLI workflow for PPTX/PPT to PDF using Microsoft PowerPoint.
---

# Pptx To Pdf Powerpoint

## Overview

Convert .pptx/.ppt files to PDFs on macOS using Microsoft PowerPoint with the bundled AppleScript at `scripts/pptx2pdf.applescript`.

## Workflow

1. Confirm prerequisites: macOS, Microsoft PowerPoint installed, files are local.
2. Choose input mode: file list (`-p`) or directory (`-d`).
3. Run the script with `osascript` or execute it directly.
4. Report the success count and any errors.

## Quick Start

```bash
osascript scripts/pptx2pdf.applescript -p "/path/file1.pptx" "/path/file2.pptx"
osascript scripts/pptx2pdf.applescript -d "/path/to/dir"
```

## Notes

- Output PDFs are written alongside inputs with a `.pdf` extension.
- Ignores temporary files starting with `~$`.
- Activates PowerPoint and processes files sequentially.
- Missing directories are logged to stderr; conversion continues for other inputs.
- If no valid files are found, the script returns a short message.

## Troubleshooting

- If macOS prompts for automation permissions, allow `osascript` to control Microsoft PowerPoint.
- If large files time out (open wait is ~5 seconds), increase the wait loop in `scripts/pptx2pdf.applescript`.

## Resources

- `scripts/pptx2pdf.applescript`
