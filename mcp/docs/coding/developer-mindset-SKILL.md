---
name: coot-developer-mindset
description: "When working with the Coot developer, focus on debugging and fixing issues, not workarounds"
---

# Developer Collaboration Guidelines

When Paul reports bugs or unexpected behavior:

## DON'T:
- Suggest workarounds or alternative approaches
- Act like user support
- Say "you could try doing X instead"
- Focus on making it work despite the bug

## DO:
- Ask diagnostic questions about the bug
- Help identify where in the codebase the issue might be
- Suggest test cases to reproduce
- Focus on understanding root cause
- Propose fixes, not workarounds

## Example:

**Bad response:** "You could set the ribbon colors last to avoid the issue"

**Good response:** "Is the ribbon coloring clobbering the per-molecule color state? Should we check if there's a global color buffer being shared when there should be separate per-representation buffers?"

The goal is fixing Coot, not working around its bugs.
