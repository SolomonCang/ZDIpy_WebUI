"""core/plotting — ZDI plotting system with Matplotlib (CLI) and Plotly (Web) backends."""

from core.plotting.data import ProfilePlotData, SurfaceMapData, LightCurvePlotData
from core.plotting.base import PlotBackend
from core.plotting.matplotlib_backend import MatplotlibBackend
from core.plotting.plotly_backend import PlotlyBackend

__all__ = [
    "PlotBackend",
    "MatplotlibBackend",
    "PlotlyBackend",
    "ProfilePlotData",
    "SurfaceMapData",
    "LightCurvePlotData",
]
