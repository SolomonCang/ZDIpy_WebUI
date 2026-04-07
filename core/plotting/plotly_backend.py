"""core/plotting/plotly_backend.py — Plotly backend (web API use).

All methods return a dict with ``"data"`` (list of traces) and ``"layout"`` keys
suitable for direct use with ``Plotly.newPlot(div, resp.data, resp.layout)``.
"""

import numpy as np
from core.plotting.base import PlotBackend
from core.plotting.data import ProfilePlotData, SurfaceMapData, LightCurvePlotData


class PlotlyBackend(PlotBackend):
    """Produce JSON-serialisable Plotly dicts for the web API."""

    def plot_profiles(self, data: ProfilePlotData) -> dict:
        """Return a single Plotly figure containing all phases.

        Each phase contributes Obs and Model traces for both Stokes I (top
        panel) and Stokes V (bottom panel).  All but the first phase are set
        to ``"legendonly"`` so the chart is readable by default; click the
        legend to show/hide phases.  The two panels share the same x-axis.
        """
        traces = []
        vel = data.vel_grid.tolist()
        for i, phase in enumerate(data.phases):
            visible = True if i == 0 else "legendonly"

            # ── Stokes I (top panel — yaxis) ─────────────────────────────
            obs_I_trace: dict = {
                "x": vel,
                "y": data.obs_I[i].tolist(),
                "mode": "markers",
                "marker": {
                    "size": 3,
                    "color": "#f0f6fc"
                },
                "name": f"I obs  φ={phase:.3f}",
                "legendgroup": f"phase{i}",
                "visible": visible,
            }
            if data.obs_I_sigma and len(data.obs_I_sigma[i]) > 0:
                obs_I_trace["error_y"] = {
                    "type": "data",
                    "array": data.obs_I_sigma[i].tolist(),
                    "visible": True,
                    "color": "#484f58",
                    "thickness": 1,
                }
            traces.append(obs_I_trace)
            traces.append({
                "x": vel,
                "y": data.mod_I[i].tolist(),
                "mode": "lines",
                "line": {
                    "color": "#f85149",
                    "width": 1.5
                },
                "name": f"model  φ={phase:.3f}",
                "legendgroup": f"phase{i}",
                "visible": visible,
            })

            # ── Stokes V (bottom panel — yaxis2) ─────────────────────────
            obs_V_trace: dict = {
                "x": vel,
                "y": data.obs_V[i].tolist(),
                "mode": "markers",
                "marker": {
                    "size": 3,
                    "color": "#5d6be5"
                },
                "name": f"V obs  φ={phase:.3f}",
                "legendgroup": f"phase{i}",
                "visible": visible,
                "yaxis": "y2",
                "showlegend": False,
            }
            if data.obs_V_sigma and len(data.obs_V_sigma[i]) > 0:
                obs_V_trace["error_y"] = {
                    "type": "data",
                    "array": data.obs_V_sigma[i].tolist(),
                    "visible": True,
                    "color": "#484f58",
                    "thickness": 1,
                }
            traces.append(obs_V_trace)
            traces.append({
                "x": vel,
                "y": data.mod_V[i].tolist(),
                "mode": "lines",
                "line": {
                    "color": "#f85149",
                    "width": 1.5
                },
                "name": f"V mod  φ={phase:.3f}",
                "legendgroup": f"phase{i}",
                "visible": visible,
                "yaxis": "y2",
                "showlegend": False,
            })

        layout = {
            "xaxis": {
                "title": {
                    "text": "v (km/s)"
                }
            },
            "yaxis": {
                "title": {
                    "text": "I/Ic"
                },
                "domain": [0.45, 1.0],
            },
            "yaxis2": {
                "title": {
                    "text": "V/Ic"
                },
                "domain": [0.0, 0.40],
                "anchor": "x",
                "zeroline": True,
                "zerolinecolor": "#30363d",
            },
            "template": "plotly_dark",
            "margin": {
                "t": 30,
                "r": 10,
                "b": 40,
                "l": 50
            },
            "legend": {
                "orientation": "v"
            },
        }
        return {"data": traces, "layout": layout}

    def plot_surface_map(self, data: SurfaceMapData) -> dict:
        """Return a scatter-based equal-area surface map in Plotly JSON."""
        lon_deg = np.degrees(data.lon).tolist()
        lat_deg = (90.0 - np.degrees(data.clat)).tolist()
        colorscale = "RdBu_r" if "B" in data.map_type else "Hot"
        trace = {
            "type": "scatter",
            "x": lon_deg,
            "y": lat_deg,
            "mode": "markers",
            "marker": {
                "color": data.values.tolist(),
                "colorscale": colorscale,
                "reversescale": data.map_type == "brightness",
                "size": 4,
                "colorbar": {
                    "title": {
                        "text": data.map_type
                    }
                },
                "cmin": data.vmin,
                "cmax": data.vmax,
            },
            "name": data.map_type,
        }
        layout = {
            "xaxis": {
                "title": {
                    "text": "Longitude (deg)"
                },
                "range": [0, 360]
            },
            "yaxis": {
                "title": {
                    "text": "Latitude (deg)"
                },
                "range": [-90, 90]
            },
            "template": "plotly_dark",
            "margin": {
                "t": 30,
                "r": 10,
                "b": 40,
                "l": 50
            },
        }
        return {"data": [trace], "layout": layout}

    def plot_light_curve(self, data: LightCurvePlotData) -> dict:
        jd = data.jdates.tolist()
        traces = [
            {
                "x": jd,
                "y": data.obs_flux.tolist(),
                "mode": "markers",
                "error_y": {
                    "type": "data",
                    "array": data.sigma.tolist(),
                    "visible": True,
                    "color": "#484f58",
                    "thickness": 1,
                },
                "marker": {
                    "color": "#f0f6fc",
                    "size": 5
                },
                "name": "Observed",
            },
            {
                "x": jd,
                "y": data.mod_flux.tolist(),
                "mode": "lines",
                "line": {
                    "color": "#f85149",
                    "width": 2
                },
                "name": "Model",
            },
        ]
        layout = {
            "xaxis": {
                "title": {
                    "text": "JD"
                }
            },
            "yaxis": {
                "title": {
                    "text": "Flux"
                }
            },
            "template": "plotly_dark",
            "margin": {
                "t": 30,
                "r": 10,
                "b": 40,
                "l": 50
            },
        }
        return {"data": traces, "layout": layout}
