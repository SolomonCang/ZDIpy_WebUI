"""core/plotting/matplotlib_backend.py — Matplotlib backend (CLI use)."""

import numpy as np
from core.plotting.base import PlotBackend
from core.plotting.data import ProfilePlotData, SurfaceMapData, LightCurvePlotData


class MatplotlibBackend(PlotBackend):
    """Produce matplotlib Figure objects for CLI/offline use."""
    def plot_profiles(self, data: ProfilePlotData):
        import matplotlib.pyplot as plt  # deferred import — only needed for CLI

        n = len(data.phases)
        fig, axes = plt.subplots(2, n, figsize=(3 * n, 5), sharex=True)
        if n == 1:
            axes = axes.reshape(2, 1)
        for i, phase in enumerate(data.phases):
            ax_I, ax_V = axes[0, i], axes[1, i]
            ax_I.plot(data.vel_grid, data.obs_I[i], 'k.', ms=2, label='Obs')
            ax_I.plot(data.vel_grid, data.mod_I[i], 'r-', lw=1, label='Model')
            if data.obs_I_sigma:
                ax_I.fill_between(data.vel_grid,
                                  data.obs_I[i] - data.obs_I_sigma[i],
                                  data.obs_I[i] + data.obs_I_sigma[i],
                                  alpha=0.3,
                                  color='grey')
            ax_V.plot(data.vel_grid, data.obs_V[i], 'k.', ms=2)
            ax_V.plot(data.vel_grid, data.mod_V[i], 'r-', lw=1)
            ax_I.set_title(f'φ={phase:.3f}', fontsize=8)
        axes[0, 0].set_ylabel('I/Ic')
        axes[1, 0].set_ylabel('V/Ic')
        fig.tight_layout()
        return fig

    def plot_surface_map(self, data: SurfaceMapData):
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        lon_deg = np.degrees(data.lon)
        lat_deg = 90.0 - np.degrees(data.clat)
        colormap = 'RdBu_r' if 'B' in data.map_type else 'hot_r'
        sc = ax.scatter(lon_deg,
                        lat_deg,
                        c=data.values,
                        vmin=data.vmin,
                        vmax=data.vmax,
                        cmap=colormap,
                        s=4,
                        rasterized=True)
        plt.colorbar(sc, ax=ax, label=data.map_type)
        ax.set_xlabel('Longitude (deg)')
        ax.set_ylabel('Latitude (deg)')
        ax.set_title(data.map_type)
        fig.tight_layout()
        return fig

    def plot_light_curve(self, data: LightCurvePlotData):
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7, 3))
        ax.errorbar(data.jdates,
                    data.obs_flux,
                    yerr=data.sigma,
                    fmt='k.',
                    ms=4,
                    capsize=2,
                    label='Obs')
        ax.plot(data.jdates, data.mod_flux, 'r-', lw=1.5, label='Model')
        ax.set_xlabel('JD')
        ax.set_ylabel('Flux')
        ax.legend()
        fig.tight_layout()
        return fig
