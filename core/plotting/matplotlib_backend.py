"""core/plotting/matplotlib_backend.py — Matplotlib backend (CLI use)."""

import numpy as np
from core.plotting.base import PlotBackend
from core.plotting.data import ProfilePlotData, SurfaceMapData, LightCurvePlotData, MagneticPolarData, BrightnessPolarData, BrightnessPolarData


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

    def plot_magnetic_polar(self, data: MagneticPolarData):
        """Three-panel polar projection magnetic field map (Radial/Azimuthal/Meridional).

        Returns:
            matplotlib Figure — save with fig.savefig(...) before plt.close().
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors

        npLon = len(data.lon_grid)

        # ---- coordinate conversion ----------------------------------------
        lat_1d = np.pi / 2.0 - data.clat_grid  # latitude in radians
        lat_deg_1d = lat_1d * 180.0 / np.pi  # latitude in degrees
        r_1d = 90.0 - lat_deg_1d  # polar r: 0=pole, 90=equator
        fullLonPlt_1d = np.pi - data.lon_grid  # east-west flip for display

        # 2-D meshes (npClat × npLon)
        lon_2d, r_2d = np.meshgrid(fullLonPlt_1d, r_1d)
        lat_deg_2d = np.tile(lat_deg_1d.reshape(-1, 1), (1, npLon))

        # ---- southern hemisphere mask  (lat < -30°) ------------------------
        mask = lat_deg_2d < -30.0
        Br_m = np.where(mask, np.nan, data.Br)
        Blon_m = np.where(mask, np.nan, data.Blon)
        Blat_m = np.where(mask, np.nan, -data.Blat)  # sign convention flip

        # ---- unified colour scale -----------------------------------------
        BMax = max(
            np.nanmax(np.abs(Br_m)),
            np.nanmax(np.abs(Blon_m)),
            np.nanmax(np.abs(Blat_m)),
        )
        if BMax == 0:
            BMax = 1.0

        n = data.discrete_levels
        levels = np.linspace(-BMax, BMax, 2 * n + 1)
        cmap_base = cm.get_cmap('RdBu_r', 2 * n)
        colors_arr = cmap_base(np.linspace(0, 1, 2 * n))
        # two central colors → white (near-zero values)
        colors_arr[n - 1] = [1.0, 1.0, 1.0, 1.0]
        colors_arr[n] = [1.0, 1.0, 1.0, 1.0]
        new_cmap = mcolors.ListedColormap(colors_arr)
        norm = mcolors.BoundaryNorm(levels, new_cmap.N)

        # ---- figure -------------------------------------------------------
        fig = plt.figure(figsize=(13, 5))
        titles = ['Radial', 'Azimuthal', 'Meridional']
        fields = [Br_m, Blon_m, Blat_m]
        RLIM = (0, 97)
        TICK_INNER = 91.5
        theta_circle = np.linspace(0, 2.0 * np.pi, 360)

        for i, (title, field) in enumerate(zip(titles, fields)):
            ax = fig.add_subplot(1, 3, i + 1, projection='polar')
            ax.set_theta_zero_location('S')
            ax.set_theta_direction(-1)

            ax.pcolormesh(lon_2d,
                          r_2d,
                          field,
                          cmap=new_cmap,
                          norm=norm,
                          shading='auto')

            # equator ring
            ax.plot(theta_circle, np.full(360, 90.0), 'k-', linewidth=0.7)

            # observation phase tick marks
            for phase in data.obs_phases:
                theta_t = np.pi - 2.0 * np.pi * float(phase)
                ax.plot([theta_t, theta_t], [TICK_INNER, RLIM[1]],
                        'k-',
                        linewidth=1.5)

            ax.set_rlim(*RLIM)
            ax.set_rticks([30, 60, 90])
            ax.set_yticklabels(['60°', '30°', '0°'], fontsize=7)
            ax.tick_params(axis='x', labelsize=8)
            ax.set_title(title, pad=10, fontsize=11)

        # ---- shared horizontal colorbar -----------------------------------
        cbar_ax = fig.add_axes([0.15, 0.08, 0.70, 0.035])
        sm = plt.cm.ScalarMappable(cmap=new_cmap, norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm,
                          cax=cbar_ax,
                          orientation='horizontal',
                          ticks=levels[::2])
        cb.set_label('Magnetic field (G)', fontsize=9)
        cb.ax.tick_params(labelsize=8)

        fig.subplots_adjust(left=0.02,
                            right=0.98,
                            top=0.90,
                            bottom=0.18,
                            wspace=0.35)
        return fig

    def plot_brightness_polar(self, data: BrightnessPolarData):
        """Single-panel polar projection brightness map.

        Returns:
            matplotlib Figure — save with fig.savefig(...) before plt.close().
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors

        npLon = len(data.lon_grid)

        # ---- coordinate conversion ----------------------------------------
        lat_1d = np.pi / 2.0 - data.clat_grid
        lat_deg_1d = lat_1d * 180.0 / np.pi
        r_1d = 90.0 - lat_deg_1d
        fullLonPlt_1d = np.pi - data.lon_grid

        lon_2d, r_2d = np.meshgrid(fullLonPlt_1d, r_1d)
        lat_deg_2d = np.tile(lat_deg_1d.reshape(-1, 1), (1, npLon))

        # ---- southern hemisphere mask (lat < -30°) -------------------------
        mask = lat_deg_2d < -30.0
        bri_m = np.where(mask, np.nan, data.brightness)

        # ---- discrete colormap (reversed Hot) --------------------------------
        n = data.discrete_levels
        bri_min = max(0.0, float(np.nanmin(bri_m)))
        bri_max = max(float(np.nanmax(bri_m)), bri_min + 1e-6)
        levels = np.linspace(bri_min, bri_max, 2 * n + 1)
        cmap_base = cm.get_cmap('hot_r', 2 * n)
        new_cmap = mcolors.ListedColormap(cmap_base(np.linspace(0, 1, 2 * n)))
        norm = mcolors.BoundaryNorm(levels, new_cmap.N)

        # ---- figure -------------------------------------------------------
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(1, 1, 1, projection='polar')
        ax.set_theta_zero_location('S')
        ax.set_theta_direction(-1)

        ax.pcolormesh(lon_2d,
                      r_2d,
                      bri_m,
                      cmap=new_cmap,
                      norm=norm,
                      shading='auto')

        RLIM = (0, 97)
        TICK_INNER = 91.5
        theta_circle = np.linspace(0, 2.0 * np.pi, 360)

        # equator ring
        ax.plot(theta_circle, np.full(360, 90.0), 'k-', linewidth=0.7)

        # observation phase tick marks
        for phase in data.obs_phases:
            theta_t = np.pi - 2.0 * np.pi * float(phase)
            ax.plot([theta_t, theta_t], [TICK_INNER, RLIM[1]],
                    'k-',
                    linewidth=1.5)

        ax.set_rlim(*RLIM)
        ax.set_rticks([30, 60, 90])
        ax.set_yticklabels(['60°', '30°', '0°'], fontsize=8)
        ax.tick_params(axis='x', labelsize=9)
        ax.set_title('Brightness Map', pad=12, fontsize=12)

        # ---- colorbar at bottom -------------------------------------------
        cbar_ax = fig.add_axes([0.15, 0.06, 0.70, 0.035])
        sm = plt.cm.ScalarMappable(cmap=new_cmap, norm=norm)
        sm.set_array([])
        cb = fig.colorbar(sm,
                          cax=cbar_ax,
                          orientation='horizontal',
                          ticks=levels[::2])
        cb.set_label('Brightness (normalised)', fontsize=9)
        cb.ax.tick_params(labelsize=8)

        fig.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.14)
        return fig
