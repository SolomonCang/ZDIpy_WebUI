#!/usr/bin/env python3
"""
app.py - Main entry point for ZDIpy WebUI

Run with:
    python app.py              # start FastAPI/uvicorn server on localhost:7860
    python app.py --share      # bind to 0.0.0.0 (network-accessible)
    python app.py --port 8080  # use custom port
    python app.py --cli --config config/config.json  # CLI mode (no WebUI)
"""

import argparse
import sys
import os

# Ensure the project root is in sys.path
_ROOT = os.path.dirname(os.path.abspath(__file__))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)


def main():
    parser = argparse.ArgumentParser(
        description="ZDIpy WebUI – Zeeman Doppler Imaging web application",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--share",
        action="store_true",
        help="Create a public shareable Gradio link",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=7860,
        help="Port to run the web server on",
    )
    parser.add_argument(
        "--cli",
        action="store_true",
        help="Run in CLI mode (no WebUI, just run the ZDI pipeline)",
    )
    parser.add_argument(
        "--config",
        default=os.path.join(_ROOT, "config", "config.json"),
        help="Path to config.json (used in --cli mode)",
    )
    parser.add_argument(
        "--forward-only",
        action="store_true",
        help="Run only forward model (no inversion, used in --cli mode)",
    )
    parser.add_argument(
        "--verbose",
        type=int,
        default=1,
        choices=[0, 1, 2],
        help="Verbosity level",
    )
    args = parser.parse_args()

    if args.cli:
        # CLI mode: run ZDI pipeline directly
        from zdi_runner import run_zdi
        result = run_zdi(
            config_path=args.config,
            forward_only=args.forward_only,
            verbose=args.verbose,
        )
        print("\nResult:")
        for k, v in result.items():
            print(f"  {k}: {v}")
    else:
        # WebUI mode: launch FastAPI server via uvicorn
        import uvicorn
        host = "0.0.0.0" if args.share else "127.0.0.1"
        url = f"http://127.0.0.1:{args.port}"
        print(f"Starting ZDIpy WebUI on {url}")
        print(f"  API docs: {url}/docs")
        if args.share:
            print(f"  Network access enabled (binding to 0.0.0.0:{args.port})")
        uvicorn.run(
            "api.server:app",
            host=host,
            port=args.port,
            reload=False,
            log_level="warning",
        )


if __name__ == "__main__":
    main()
