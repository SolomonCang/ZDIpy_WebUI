"""core/plotting/base.py — PlotBackend abstract base class."""

from abc import ABC, abstractmethod
from typing import Any

from core.plotting.data import ProfilePlotData, SurfaceMapData, LightCurvePlotData


class PlotBackend(ABC):
    """Plotting backend interface.

    MatplotlibBackend returns matplotlib.figure.Figure objects (for CLI).
    PlotlyBackend returns JSON-serialisable dicts (for the web API).
    """
    @abstractmethod
    def plot_profiles(self, data: ProfilePlotData) -> Any:
        """Render Stokes I/V line profile comparison plot."""

    @abstractmethod
    def plot_surface_map(self, data: SurfaceMapData) -> Any:
        """Render stellar surface map (equal-area scatter projection)."""

    @abstractmethod
    def plot_light_curve(self, data: LightCurvePlotData) -> Any:
        """Render light curve comparison plot."""
