#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

PORT="${PORT:-7860}"

if ! command -v python3 >/dev/null 2>&1; then
  echo "Error: python3 is not installed. Please install Python 3 first."
  read -n 1 -r -p "Press any key to close..." _
  echo
  exit 1
fi

if [ ! -d ".venv" ]; then
  echo "Creating virtual environment..."
  python3 -m venv .venv
fi

source .venv/bin/activate

echo "Installing dependencies (if needed)..."
python -m pip install --upgrade pip >/dev/null
python -m pip install -r requirements.txt >/dev/null

echo "Starting ZDIpy WebUI on http://127.0.0.1:${PORT}"

# Open browser shortly after server starts.
( sleep 2 && open "http://127.0.0.1:${PORT}" ) >/dev/null 2>&1 &

python app.py --port "$PORT"

EXIT_CODE=$?
echo
echo "WebUI exited with code ${EXIT_CODE}."
read -n 1 -r -p "Press any key to close..." _
echo
exit "$EXIT_CODE"