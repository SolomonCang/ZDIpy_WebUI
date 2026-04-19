#!/usr/bin/env python3
"""
lsd_interactive.py — 启动 LSD Interactive Web 应用（Voigt / UR 模型交互式拟合）。

用法（在项目根目录 ZDIpy_WebUI/ 下）：
    python scripts/lsd_interactive.py [--port PORT] [--no-browser]

功能：
  - 加载 config.json 与 LSD 观测数据
  - 启动 FastAPI 后端 + Web 前端
  - 自动在浏览器中打开交互界面
  - 支持 Voigt / Unno-Rachkovsky 模型切换
  - 支持 Auto-estimate line strength（最小化 RMS）
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import threading
import time
import webbrowser
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
_SERVER = _SCRIPT_DIR / "lsd_interactive_web" / "server.py"


def main():
    parser = argparse.ArgumentParser(
        description="LSD Interactive Fit — 启动 Web 服务")
    parser.add_argument("--host", default="127.0.0.1",
                        help="绑定地址 (default: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=8051,
                        help="端口号 (default: 8051)")
    parser.add_argument("--no-browser", action="store_true",
                        help="不自动打开浏览器")
    args = parser.parse_args()

    url = f"http://{args.host}:{args.port}"

    if not args.no_browser:
        def _open():
            time.sleep(1.5)
            webbrowser.open(url)
        threading.Thread(target=_open, daemon=True).start()

    print(f"Starting LSD Interactive Web at {url}")
    print("Press Ctrl+C to stop.\n")

    subprocess.run(
        [sys.executable, str(_SERVER),
         "--host", args.host, "--port", str(args.port)],
        cwd=str(_SCRIPT_DIR.parent),   # ZDIpy_WebUI/
    )


if __name__ == "__main__":
    main()
