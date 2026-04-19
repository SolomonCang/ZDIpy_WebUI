#!/usr/bin/env python3
"""启动 LSD Interactive Web 应用的便捷入口。"""
import subprocess
import sys
from pathlib import Path

if __name__ == "__main__":
    server = Path(__file__).parent / "server.py"
    subprocess.run([sys.executable, str(server)] + sys.argv[1:])
