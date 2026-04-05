---
description: "Use when editing ZDIpy_WebUI Python code, Gradio WebUI flow, CLI entrypoints, JSON config loading, or ZDI pipeline integration. Enforces project structure, naming, safe edit boundaries, and run/verification habits."
name: "ZDIpy WebUI Project Rules"
applyTo: ["**/*.py", "config/**/*.json", "README.md", "start_webui.command"]
---
# ZDIpy WebUI Instructions

- Keep entrypoint logic in `app.py`, CLI execution in `zdi_runner.py`, and UI behavior in `webui/app.py`.
- Treat `core/` as scientific engine code. Do not refactor it broadly unless the task is specifically about numerical correctness or API compatibility.
- Route user-facing parameters through `config/config.json` and `config_loader.py` (`ZDIConfig`) instead of adding ad hoc config sources.
- Preserve the existing naming style used in this project (camelCase-heavy function and class names in core modules). Match local style before introducing new names.
- Keep imports stable for scientific stack (`numpy as np`, SciPy and matplotlib usage) and avoid changing numerical defaults without a clear reason.
- For path handling, prefer project-root-relative paths and existing root-resolution patterns already used by entry scripts.
- For WebUI changes, keep Gradio interactions deterministic and thread-safe; preserve run locking and clear user-facing error messages.
- For CLI changes, preserve parity with WebUI config behavior so both modes consume equivalent settings.
- Do not silently change observational data assumptions. Treat `LSDprof/` file usage and observation path mapping as explicit config behavior.
- Keep edits minimal and localized. Avoid broad style-only rewrites in legacy scientific modules.
- When behavior changes, update `README.md` run instructions and flags if user-facing commands or defaults are affected.
- Prefer validating with the documented flow: `pip install -r requirements.txt`, then `python app.py` (WebUI) and `python app.py --cli --config config/config.json` (CLI).
