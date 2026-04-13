"""Regression tests for ha_num model EW-fitting behavior.

Run with:
    python -m pytest tests/test_ha_num_regression.py -q
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import core.line_models as line_models
from config_loader import ZDIConfig
from pipeline.pipeline import ZDIPipeline


def test_ha_num_skips_voigt_ew_fitting(monkeypatch):
    """ha_num must not call Voigt EW fitter even when estimate_strength=1."""
    config_path = os.path.join(os.path.dirname(__file__), "..", "config.json")
    par = ZDIConfig(config_path)

    par.line_model_type = "ha_num"
    par.estimateStrenght = 1

    # Keep runtime short while still traversing the full model-selection path.
    par.numIterations = 0

    def _fail_if_called(*_args, **_kwargs):
        raise AssertionError(
            "fitLineStrength should not be called for ha_num model")

    monkeypatch.setattr(line_models, "fitLineStrength", _fail_if_called)

    logs = []
    result = ZDIPipeline(par,
                         forward_only=True,
                         verbose=0,
                         callback=logs.append).run()

    assert result.iterations == 0
    assert any("estimate_strength skipped" in msg for msg in logs)
