#!/usr/bin/env python3
"""
zdi_runner.py - Main ZDI runner that accepts a config.json file.

Usage:
    python zdi_runner.py                          # uses config.json
    python zdi_runner.py --config my_config.json  # custom config path
    python zdi_runner.py --forward-only           # run forward model only
    python zdi_runner.py --verbose 2              # verbose output

This script replaces the legacy zdipy.py + inzdi.dat workflow with a
cleaner config.json-based approach while preserving full compatibility
with the ZDIpy physical models.
"""

__version__ = "1.0.0"

import argparse
import logging
import os
import sys

_log = logging.getLogger("zdipy")
if not _log.handlers:
    _log.addHandler(logging.NullHandler())  # library-style: no default output


def run_zdi(config_path: str,
            forward_only: bool = False,
            verbose: int = 1,
            progress_callback=None,
            stop_event=None):
    """
    Run the ZDI forward model and/or inversion from a config.json file.

    Parameters
    ----------
    config_path : str
        Path to the config.json file.
    forward_only : bool
        If True, only compute the forward model (no MEM inversion).
    verbose : int
        Verbosity level (0=silent, 1=normal, 2=detailed).
    progress_callback : callable, optional
        Called with a single ``str`` argument for every progress message.
        If None, messages go to stdout only.

    Returns
    -------
    ZDIResult
        Result object.  Scalar summary fields are available directly
        (``iterations``, ``entropy``, ``chi2``, ``test``, ``converged``).
        Extended metadata (``mean_bright``, ``mean_mag``, …) lives in
        ``result.metadata``.  Full maps and profiles are in
        ``bright_map``, ``mag_coeffs``, ``synthetic_profiles``, etc.
    """
    import config_loader as cl
    from pipeline.pipeline import ZDIPipeline

    def _emit(msg: str) -> None:
        """Print to stdout and forward to callback when set."""
        if verbose >= 1:
            print(msg)
        if progress_callback is not None:
            progress_callback(msg)

    _emit(f"ZDIpy_WebUI v{__version__}")
    _emit(f"Loading configuration from: {config_path}")

    # --- Load configuration (all paths resolved to absolute inside ZDIConfig) ---
    par = cl.ZDIConfig(config_path)

    pipeline = ZDIPipeline(
        par,
        forward_only=forward_only,
        verbose=verbose,
        callback=progress_callback,
        stop_event=stop_event,
    )
    result = pipeline.run()
    result.metadata["config_path"] = config_path
    return result


def main():
    parser = argparse.ArgumentParser(
        description="ZDIpy_WebUI – Zeeman Doppler Imaging with config.json",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config",
        default="config.json",
        help="Path to the JSON configuration file",
    )
    parser.add_argument(
        "--forward-only",
        action="store_true",
        help="Run only the forward model (no MEM inversion)",
    )
    parser.add_argument(
        "--verbose",
        type=int,
        default=1,
        choices=[0, 1, 2],
        help="Verbosity level",
    )
    args = parser.parse_args()

    config_path = args.config
    if not os.path.isfile(config_path):
        print(f"Error: config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)

    result = run_zdi(config_path,
                     forward_only=args.forward_only,
                     verbose=args.verbose)
    print("\nResult summary:")
    summary = result.to_serializable()
    for k in ("iterations", "entropy", "chi2", "test", "converged"):
        print(f"  {k}: {summary[k]}")


if __name__ == "__main__":
    main()
