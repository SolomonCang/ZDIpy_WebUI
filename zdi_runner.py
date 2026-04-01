#!/usr/bin/env python3
"""
zdi_runner.py - Main ZDI runner that accepts a config.json file.

Usage:
    python zdi_runner.py                          # uses config/config.json
    python zdi_runner.py --config my_config.json  # custom config path
    python zdi_runner.py --forward-only           # run forward model only
    python zdi_runner.py --verbose 2              # verbose output

This script replaces the legacy zdipy.py + inzdi.dat workflow with a
cleaner config.json-based approach while preserving full compatibility
with the ZDIpy physical models.
"""

__version__ = "1.0.0"

import argparse
import os
import sys
from pathlib import Path


def run_zdi(config_path: str, forward_only: bool = False, verbose: int = 1) -> dict:
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

    Returns
    -------
    dict
        Result dictionary with keys:
        - 'iterations': number of fitting iterations
        - 'entropy': final entropy value
        - 'chi2': final chi-squared per data point
        - 'test': final Skilling-Bryan test statistic
        - 'mean_bright': mean brightness
        - 'mean_bright_diff': mean brightness deviation from default
        - 'mean_mag': mean magnetic field strength (Gauss)
        - 'converged': True if convergence criterion was met
        - 'config_path': path to the config file used
    """
    import config_loader as cl
    import core.mainFuncs as mf
    import core.readObs as readObs
    import core.geometryStellar as geometryStellar
    import core.magneticGeom as magneticGeom
    import core.brightnessGeom as brightnessGeom
    import core.lineprofileVoigt as lineprofile
    import core.memSimple3 as memSimple

    if verbose >= 1:
        print(f"ZDIpy_WebUI v{__version__}")
        print(f"Loading configuration from: {config_path}")

    # --- Load configuration -----------------------------------------------
    par = cl.ZDIConfig(config_path)

    # Change working directory to config file directory so relative paths work
    config_dir = str(Path(config_path).parent.resolve())
    original_dir = os.getcwd()
    if config_dir != original_dir:
        os.chdir(config_dir)

    try:
        result = _run_pipeline(par, forward_only, verbose, mf, readObs,
                               geometryStellar, magneticGeom, brightnessGeom,
                               lineprofile, memSimple)
    finally:
        os.chdir(original_dir)

    result["config_path"] = config_path
    return result


def _run_pipeline(par, forward_only, verbose, mf, readObs, geometryStellar,
                  magneticGeom, brightnessGeom, lineprofile, memSimple):
    """Internal pipeline execution."""

    # --- Prepare parameters -----------------------------------------------
    par.setTarget()
    par.setCalcdIdV(verbose)
    par.calcCycles(verbose)

    if forward_only:
        par.numIterations = 0
        if verbose >= 1:
            print("Forward-only mode: numIterations set to 0")

    # --- Load line model data ---------------------------------------------
    model_file = par.model_file if hasattr(par, "model_file") else "model-voigt-line.dat"
    lineData = lineprofile.lineData(model_file, par.instrumentRes)

    # --- Load observed spectra --------------------------------------------
    obsSet = readObs.obsProfSetInRange(
        par.fnames, par.velStart, par.velEnd, par.velRs
    )

    # --- Build wavelength grid --------------------------------------------
    wlSynSet, nDataTot = mf.getWavelengthGrid(par.velRs, obsSet, lineData, verbose)

    # Scale Stokes I error bars
    for obs in obsSet:
        obs.scaleIsig(par.chiScaleI)

    # --- Initialize stellar grid ------------------------------------------
    sGrid = geometryStellar.starGrid(
        par.nRingsStellarGrid, par.period, par.mass, par.radius, verbose
    )

    # --- Initialize magnetic geometry -------------------------------------
    magGeom = magneticGeom.SetupMagSphHarmoics(
        sGrid, par.initMagFromFile, par.initMagGeomFile,
        par.lMax, par.magGeomType, verbose
    )

    # --- Initialize brightness map ----------------------------------------
    briMap = brightnessGeom.SetupBrightMap(
        sGrid, par.initBrightFromFile, par.initBrightFile,
        par.defaultBright, verbose
    )

    # --- Pre-calculate geometry for each phase ----------------------------
    listGridView = geometryStellar.getListGridView(par, sGrid)

    # --- Magnetic vectors and derivatives ---------------------------------
    vecMagCart = magGeom.getAllMagVectorsCart()
    if par.fitMag == 1:
        dMagCart0 = magGeom.getAllMagDerivsCart()
    else:
        dMagCart0 = 0.0

    # --- Optionally estimate line strength from EW -----------------------
    if par.estimateStrenght == 1:
        meanEW = readObs.getObservedEW(obsSet, lineData, verbose)
        mf.fitLineStrength(meanEW, par, listGridView, vecMagCart,
                           dMagCart0, briMap, lineData, wlSynSet, verbose)

    # --- Initialize synthetic spectra objects -----------------------------
    setSynSpec = lineprofile.getAllProfDiriv(
        par, listGridView, vecMagCart, dMagCart0, briMap, lineData, wlSynSet
    )

    # --- MEM constants and data vectors -----------------------------------
    constMem = memSimple.constantsMEM(par, briMap, magGeom, nDataTot)
    Data, sig2 = memSimple.packDataVector(obsSet, par.fitBri, par.fitMag)

    if (par.calcDV == 1) and (par.calcDI != 1):
        allModeldIdV = memSimple.packResponseMatrix(
            setSynSpec, nDataTot, constMem.npBriMap,
            magGeom, par.magGeomType, par.fitBri, par.fitMag
        )
        par.calcDV = 0
    else:
        allModeldIdV = 0

    weightEntropy = memSimple.setEntropyWeights(par, magGeom, sGrid)

    # --- Run main fitting loop -------------------------------------------
    iIter, entropy, chi2, test, meanBright, meanBrightDiff, meanMag = \
        mf.mainFittingLoop(
            par, lineData, wlSynSet, sGrid, briMap, magGeom,
            listGridView, dMagCart0, setSynSpec, constMem,
            nDataTot, Data, sig2, allModeldIdV, weightEntropy, verbose
        )

    # --- Save outputs -----------------------------------------------------
    out_mag = getattr(par, "outMagCoeffFile", "outMagCoeff.dat")
    out_bri = getattr(par, "outBrightMapFile", "outBrightMap.dat")
    out_bri_gd = getattr(par, "outBrightMapGDarkFile", "outBrightMapGDark.dat")
    out_spec = getattr(par, "outLineModelsFile", "outLineModels.dat")
    out_obs = getattr(par, "outObservedFile", "outObserved.dat")

    magGeom.saveToFile(out_mag, compatibility=True)
    brightnessGeom.saveMap(briMap, out_bri)

    briMapGDark = brightnessGeom.brightMap(sGrid.clat, sGrid.long)
    briMapGDark.bright = briMap.bright * sGrid.gravityDarkening(lineData.gravDark)
    brightnessGeom.saveMap(briMapGDark, out_bri_gd)

    mf.saveModelProfs(par, setSynSpec, lineData, out_spec)
    mf.saveObsUsed(obsSet, out_obs)

    if verbose >= 1:
        print(f"\nFitting complete after {iIter} iterations")
        print(f"  Entropy:     {entropy:.5f}")
        print(f"  chi2/dof:    {chi2:.6f}")
        print(f"  Test:        {test:.6f}")
        print(f"  Mean bright: {meanBright:.7f}")
        print(f"  Mean |B|:    {meanMag:.4f} G")

    chi_aim = par.chiTarget * float(constMem.nDataTotIV)
    converged = (
        (chi2 * constMem.nDataTotIV <= chi_aim * 1.001) and (test < par.test_aim)
        if par.fixedEntropy == 0
        else (entropy >= par.ent_aim * 1.001) and (test < par.test_aim)
    )

    return {
        "iterations": iIter,
        "entropy": float(entropy),
        "chi2": float(chi2),
        "test": float(test),
        "mean_bright": float(meanBright),
        "mean_bright_diff": float(meanBrightDiff),
        "mean_mag": float(meanMag),
        "converged": bool(converged),
    }


def main():
    parser = argparse.ArgumentParser(
        description="ZDIpy_WebUI – Zeeman Doppler Imaging with config.json",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config",
        default="config/config.json",
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

    result = run_zdi(config_path, forward_only=args.forward_only, verbose=args.verbose)
    print("\nResult summary:")
    for k, v in result.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
