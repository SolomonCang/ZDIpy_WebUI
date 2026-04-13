#!/usr/bin/env python3
"""Python launcher that mirrors start_webui.command behavior."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import threading
import webbrowser
from pathlib import Path


ROOT = Path(__file__).resolve().parent
VENV_DIR = ROOT / ".venv"
VENV_PYTHON = VENV_DIR / "bin" / "python"
RELAUNCH_ENV = "ZDI_WEBUI_LAUNCHER_REEXEC"


def _is_running_in_project_venv() -> bool:
    return Path(sys.executable).resolve() == VENV_PYTHON.resolve()


def _prompt_before_exit() -> None:
    if not sys.stdin.isatty():
        return
    try:
        input("Press Enter to close...")
    except EOFError:
        pass


def _ensure_supported_python() -> None:
    if sys.version_info < (3, 8):
        print("Error: Python 3.8 or newer is required.", file=sys.stderr)
        _prompt_before_exit()
        raise SystemExit(1)


def _ensure_virtualenv() -> None:
    if VENV_DIR.is_dir():
        return
    print("Creating virtual environment...")
    subprocess.check_call([sys.executable, "-m", "venv", str(VENV_DIR)], cwd=ROOT)


def _relaunch_inside_venv() -> None:
    if _is_running_in_project_venv():
        return
    env = os.environ.copy()
    env[RELAUNCH_ENV] = "1"
    raise SystemExit(
        subprocess.call([str(VENV_PYTHON), str(Path(__file__).resolve()), *sys.argv[1:]], cwd=ROOT, env=env)
    )


def _install_dependencies() -> None:
    print("Installing dependencies (if needed)...")
    with open(os.devnull, "wb") as devnull:
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "--upgrade", "pip"],
            cwd=ROOT,
            stdout=devnull,
            stderr=subprocess.STDOUT,
        )
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "-r", "requirements.txt"],
            cwd=ROOT,
            stdout=devnull,
            stderr=subprocess.STDOUT,
        )


def _schedule_browser_open(port: int) -> None:
    url = f"http://127.0.0.1:{port}"
    timer = threading.Timer(2.0, lambda: webbrowser.open(url))
    timer.daemon = True
    timer.start()


def _run_app(port: int, share: bool) -> int:
    cmd = [sys.executable, "app.py", "--port", str(port)]
    if share:
        cmd.append("--share")
    return subprocess.call(cmd, cwd=ROOT)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Bootstrap the ZDIpy WebUI environment and launch the frontend.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--port",
        type=int,
        default=int(os.environ.get("PORT", "7860")),
        help="Port to run the web server on",
    )
    parser.add_argument(
        "--share",
        action="store_true",
        help="Bind the server to 0.0.0.0 via app.py --share",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not open the browser automatically",
    )
    args = parser.parse_args()

    os.chdir(ROOT)
    _ensure_supported_python()
    _ensure_virtualenv()
    if os.environ.get(RELAUNCH_ENV) != "1":
        _relaunch_inside_venv()

    _install_dependencies()
    print(f"Starting ZDIpy WebUI on http://127.0.0.1:{args.port}")
    if not args.no_browser:
        _schedule_browser_open(args.port)

    try:
        return _run_app(args.port, args.share)
    finally:
        print()


if __name__ == "__main__":
    exit_code = 1
    try:
        exit_code = main()
    except subprocess.CalledProcessError as exc:
        print(f"Launcher failed with exit code {exc.returncode}.", file=sys.stderr)
        exit_code = exc.returncode
    except KeyboardInterrupt:
        print("\nLauncher interrupted.")
        exit_code = 130

    print(f"WebUI exited with code {exit_code}.")
    _prompt_before_exit()
    raise SystemExit(exit_code)
