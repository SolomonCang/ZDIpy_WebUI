import numpy as np
import json
import os
import sys

# Add project root to sys.path to resolve 'core' module when running script directly
sys.path.insert(0,
                os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from core.line_models import halpha_preproc
from core import readObs


def test_auto_estimate_halpha_params():
    # Read config.json for observation files and velocities
    config_path = os.path.join(os.path.dirname(__file__), '../config.json')
    with open(config_path, 'r') as f:
        config = json.load(f)
    obs_files = config["observations"]["files"]
    file_list = [
        os.path.join(os.path.dirname(__file__), "../" + f["filename"])
        for f in obs_files
    ]
    vel_rs = np.array([f["vel_center_kms"] for f in obs_files])

    # Load observation files
    obsSet = readObs.obsProfSet(file_list)

    # Run auto parameter estimation
    result = halpha_preproc.auto_estimate_halpha_params(obsSet, vel_rs, None)
    print("Auto-estimated parameters:", result["params"])
    print(
        "Fit status:", "converged"
        if result["fit_ok"] else "morphological fallback (fit_ok=False)")
    if not result["fit_ok"]:
        print(
            "  Note: curve_fit did not converge; morphological estimates are used."
        )

    print("Parameter estimation:")
    for k, v in result["params"].items():
        print(f"  {k}: {v}")


# Visualization: plot the median Stokes I profile and the fitted components
    import matplotlib.pyplot as plt
    plot_data = result.get("plot_data")
    if plot_data is not None:
        for trace in plot_data["data"]:
            label = trace.get("name", "")
            # Skip individual observation traces to avoid clutter
            if label.lower().startswith("obs "):
                continue

            x = trace.get("x")
            y = trace.get("y")
            if x is not None and y is not None:
                # Translate Chinese labels to English
                if "中值" in label:
                    label = "Median Stokes I"
                elif "拟合模型" in label:
                    label = "Fitted Model"
                elif "发射" in label:
                    label = "Emission Component"
                elif "吸收" in label:
                    label = "Absorption Component"

                # Use dashed line for components to make them distinct
                linestyle = '--' if 'Component' in label else '-'
                plt.plot(x, y, label=label, linestyle=linestyle)

        plt.xlabel("Velocity (km/s)")
        plt.ylabel("Stokes I")
        plt.title("H-alpha Auto Init Fit Result")
        plt.legend()
        plt.tight_layout()

        # Save the figure to a file just in case GUI doesn't pop up
        save_path = os.path.join(os.path.dirname(__file__),
                                 "test_halpha_auto_init.png")
        plt.savefig(save_path, dpi=300)
        print(f"\nPlot saved to: {save_path}")

        plt.show()

if __name__ == "__main__":
    test_auto_estimate_halpha_params()
