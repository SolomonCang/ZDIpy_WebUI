"""pipeline/pipeline.py - ZDIPipeline: full ZDI execution encapsulated in a class.

Replaces the bare ``_run_pipeline`` function in ``zdi_runner.py``.  The class
owns the module imports (deferred until ``run()`` is called) and a
``_log()`` helper that writes to both stdout and an optional callback,
giving the WebUI a clean push interface without touching sys.stdout.
"""

import logging
from typing import Callable, Optional
import threading

import numpy as np

_pipeline_log = logging.getLogger("zdipy.pipeline")


class ZDIPipeline:
    """
    Encapsulates a single ZDI inversion pipeline run.

    Parameters
    ----------
    par : ZDIConfig
        Loaded and fully path-resolved configuration object.
    forward_only : bool
        If True, skip MEM inversion (run forward model only).
    verbose : int
        0 = silent, 1 = normal, 2 = detailed.
    callback : callable, optional
        Called with a single ``str`` for every progress message.
        Runs on the caller's thread (usually a background worker thread).
    stop_event : threading.Event, optional
        When set, the pipeline will stop gracefully at the next iteration
        boundary and return early with ``cancelled=True``.
    """

    def __init__(
        self,
        par,
        forward_only: bool = False,
        verbose: int = 1,
        callback: Optional[Callable[[str], None]] = None,
        stop_event: Optional[threading.Event] = None,
    ):
        self.par = par
        self.forward_only = forward_only
        self.verbose = verbose
        self.callback = callback
        self.stop_event = stop_event

    # ------------------------------------------------------------------
    # Logging helper
    # ------------------------------------------------------------------

    def _log(self, msg: str) -> None:
        """Forward to Python logging and optional callback."""
        if self.verbose >= 1:
            _pipeline_log.info(msg)
        if self.callback is not None:
            self.callback(msg)

    # ------------------------------------------------------------------
    # Main entry point
    # ------------------------------------------------------------------

    def run(self) -> "ZDIResult":
        """Execute the full ZDI pipeline and return a ZDIResult object."""
        import core.fitting as mf
        import core.readObs as readObs
        import core.geometry as geometryStellar
        import core.magneticGeom as magneticGeom
        import core.brightnessGeom as brightnessGeom
        import core.line_models as lineprofile
        import core.mem.zdi_adapter as memSimple

        par = self.par

        # --- Prepare parameters -------------------------------------------
        par.setTarget()
        par.setCalcdIdV(self.verbose)
        par.calcCycles(self.verbose)

        if self.forward_only:
            par.numIterations = 0
            self._log("Forward-only mode: numIterations set to 0")

        # --- Build line model data directly from config.json ---------------
        _model_type = getattr(par, 'line_model_type', 'voigt')
        if _model_type == 'unno':
            lineData = lineprofile.lineDataUnno.from_parameters(
                wavelength_nm=par.line_wavelength_nm,
                line_strength=par.line_strength,
                gauss_width_kms=par.line_gauss_width_kms,
                lorentz_width_fraction=par.line_lorentz_width_fraction,
                lande_g=par.line_lande_g,
                limb_darkening=par.line_limb_darkening,
                gravity_darkening=par.line_gravity_darkening,
                instRes=par.instrumentRes,
                beta=getattr(par, 'unno_beta', -1.0),
                filling_factor_I=getattr(par, 'unno_filling_factor_I', 1.0),
                filling_factor_V=getattr(par, 'unno_filling_factor_V', 1.0),
            )
            self._log(
                "Line model: Unno-Rachkovsky (Milne-Eddington, full polarised RT)"
            )
        elif _model_type == 'halpha_compound':
            lineData = lineprofile.lineDataHalpha.from_parameters(
                wavelength_nm=par.line_wavelength_nm,
                lande_g=getattr(par, 'halpha_lande_g', 1.048),
                limb_darkening=par.line_limb_darkening,
                gravity_darkening=par.line_gravity_darkening,
                emission_strength=getattr(par, 'halpha_emission_strength',
                                          2.5),
                emission_gauss_kms=getattr(par, 'halpha_emission_gauss_kms',
                                           80.0),
                emission_lorentz_ratio=getattr(
                    par, 'halpha_emission_lorentz_ratio', 0.15),
                absorption_strength=getattr(par, 'halpha_absorption_strength',
                                            0.0),
                absorption_gauss_kms=getattr(par,
                                             'halpha_absorption_gauss_kms',
                                             25.0),
                absorption_lorentz_ratio=getattr(
                    par, 'halpha_absorption_lorentz_ratio', 0.10),
                filling_factor_V=getattr(par, 'halpha_filling_factor_V', 1.0),
                instRes=par.instrumentRes,
            )
            self._log(
                "Line model: H-alpha compound double-Voigt (weak-field approximation)"
            )
        else:
            lineData = lineprofile.lineData.from_parameters(
                wavelength_nm=par.line_wavelength_nm,
                line_strength=par.line_strength,
                gauss_width_kms=par.line_gauss_width_kms,
                lorentz_width_fraction=par.line_lorentz_width_fraction,
                lande_g=par.line_lande_g,
                limb_darkening=par.line_limb_darkening,
                gravity_darkening=par.line_gravity_darkening,
                instRes=par.instrumentRes,
            )
            self._log("Line model: Voigt (weak-field approximation)")

        # --- Load observed spectra ----------------------------------------
        obsSet = readObs.obsProfSetInRange(par.fnames, par.velStart,
                                           par.velEnd, par.velRs)

        # --- Build wavelength grid ----------------------------------------
        wlSynSet, nDataTot = readObs.getWavelengthGrid(par.velRs, obsSet,
                                                       lineData, self.verbose)

        # Scale Stokes I error bars
        for obs in obsSet:
            obs.scaleIsig(par.chiScaleI)
            obs.scaleVsig(getattr(par, 'chiScaleV', 1.0))

        # --- Initialize stellar grid --------------------------------------
        sGrid = geometryStellar.starGrid(par.nRingsStellarGrid, par.period,
                                         par.mass, par.radius, self.verbose)

        # --- Initialize magnetic geometry ---------------------------------
        magGeom = magneticGeom.SetupMagSphHarmoics(sGrid, par.initMagFromFile,
                                                   par.initMagGeomFile,
                                                   par.lMax, par.magGeomType,
                                                   self.verbose)

        # --- Initialize brightness map ------------------------------------
        briMap = brightnessGeom.SetupBrightMap(sGrid, par.initBrightFromFile,
                                               par.initBrightFile,
                                               par.defaultBright, self.verbose)

        # --- Pre-calculate geometry for each phase -----------------------
        listGridView = geometryStellar.build_batch_visible_grid(par, sGrid)

        # --- Magnetic vectors and derivatives ----------------------------
        vecMagCart = magGeom.getAllMagVectorsCart()
        if par.fitMag == 1:
            dMagCart0 = magGeom.getAllMagDerivsCart()
        else:
            dMagCart0 = 0.0

        # --- Optionally estimate line strength from EW -------------------
        if par.estimateStrenght == 1:
            meanEW = readObs.getObservedEW(obsSet, lineData, self.verbose)
            lineprofile.fitLineStrength(meanEW, par, listGridView, vecMagCart,
                                        dMagCart0, briMap, lineData, wlSynSet,
                                        self.verbose)

        # --- Initialize synthetic spectra objects ------------------------
        _model_type_spec = getattr(par, 'line_model_type', 'voigt')
        if _model_type_spec == 'unno':
            setSynSpec = lineprofile.getAllProfDirivUnno(
                par, listGridView, vecMagCart, dMagCart0, briMap, lineData,
                wlSynSet)
        elif _model_type_spec == 'halpha_compound':
            setSynSpec = lineprofile.getAllProfDirivHalpha(
                par, listGridView, vecMagCart, dMagCart0, briMap, lineData,
                wlSynSet)
        else:
            setSynSpec = lineprofile.getAllProfDiriv(par, listGridView,
                                                     vecMagCart, dMagCart0,
                                                     briMap, lineData,
                                                     wlSynSet)

        # --- MEM constants and data vectors ------------------------------
        constMem = memSimple.constantsMEM(par, briMap, magGeom, nDataTot)
        Data, sig2 = memSimple.packDataVector(obsSet, par.fitBri, par.fitMag)

        # For linear (Voigt) models R = dV/dCoeff is independent of B; computing
        # it once is exact and efficient. For non-linear (UR) models R(B) changes
        # with B at every iteration, so calcDV must remain 1 so the fitting loop
        # recomputes R after each MEM step (Gauss-Newton linearization).
        _model_is_linear = getattr(par, 'line_model_type',
                                   'voigt') in ('voigt', 'halpha_compound')
        if (par.calcDV == 1) and (par.calcDI != 1) and _model_is_linear:
            allModeldIdV = memSimple.packResponseMatrix(
                setSynSpec, nDataTot, constMem.npBriMap, magGeom,
                par.magGeomType, par.fitBri, par.fitMag)
            par.calcDV = 0
        else:
            allModeldIdV = 0

        weightEntropy = memSimple.setEntropyWeights(par, magGeom, sGrid)

        # --- Optionally append light-curve constraint ---------------------
        lc_obs_flux = None
        lc_jdates = None
        if getattr(par, 'fitLightCurve', 0) == 1 and getattr(
                par, 'lcFiles', []):
            import core.light_curve as lc_module
            from core.geometry import getCyclesClat, BatchVisibleGrid

            lc_entries = par.lcFiles
            lc_jdates = np.array([e['jdate'] for e in lc_entries])
            lc_obs_flux = np.array([e['flux'] for e in lc_entries])
            lc_sigma = np.array([e['flux_sigma'] for e in lc_entries])

            # Build BatchVisibleGrid for LC phases
            lc_cycles = getCyclesClat(par.period, par.dOmega, lc_jdates,
                                      par.jDateRef, sGrid.clat)
            lc_batch = BatchVisibleGrid(
                star_grid=sGrid,
                inclination=par.incRad,
                cycles_at_clat=lc_cycles,
                period=par.period,
                dOmega=par.dOmega,
            )
            lc_jac = lc_module.compute_lc_jacobian(
                sGrid,
                lc_batch,
                briMap,
                float(lineData.limbDark[0]),
                float(lineData.gravDark[0]),
            )  # (N_lc, N_cells)
            lc_syn = lc_module.compute_synthetic_photometry(
                sGrid,
                lc_batch,
                briMap,
                float(lineData.limbDark[0]),
                float(lineData.gravDark[0]),
            )  # (N_lc,)

            # Append to data vector and variance
            Data = np.concatenate([Data, lc_obs_flux])
            sig2 = np.concatenate([sig2, lc_sigma**2 * par.chi2ScaleLc])

            # Append Jacobian rows (brightness columns only; mag columns = 0)
            if allModeldIdV != 0:
                n_mag_cols = allModeldIdV.shape[1] - constMem.npBriMap
                lc_jac_full = np.concatenate(
                    [lc_jac, np.zeros((len(lc_entries), n_mag_cols))], axis=1)
                allModeldIdV = np.concatenate([allModeldIdV, lc_jac_full],
                                              axis=0)
        else:
            lc_syn = None

        # --- Run main fitting loop ---------------------------------------
        self._log("Starting MEM inversion...")
        iIter, entropy, chi2, test, meanBright, meanBrightDiff, meanMag = \
            mf.mainFittingLoop(
                par, lineData, wlSynSet, sGrid, briMap, magGeom,
                listGridView, dMagCart0, setSynSpec, constMem,
                nDataTot, Data, sig2, allModeldIdV, weightEntropy, self.verbose,
                stop_event=self.stop_event,
                callback=self.callback,
            )

        # --- Save outputs ------------------------------------------------
        import os as _os
        for _f in [
                par.outMagCoeffFile, par.outBrightMapFile,
                par.outBrightMapGDarkFile, par.outLineModelsFile,
                par.outObservedFile
        ]:
            _os.makedirs(_os.path.dirname(_f), exist_ok=True)
        magGeom.saveToFile(par.outMagCoeffFile, compatibility=True)
        brightnessGeom.saveMap(briMap, par.outBrightMapFile)

        briMapGDark = brightnessGeom.brightMap(sGrid.clat, sGrid.long)
        briMapGDark.bright = (briMap.bright *
                              sGrid.gravityDarkening(lineData.gravDark))
        brightnessGeom.saveMap(briMapGDark, par.outBrightMapGDarkFile)

        mf.saveModelProfs(par, setSynSpec, lineData, par.outLineModelsFile)
        mf.saveObsUsed(obsSet, par.outObservedFile)

        self._log(f"\nFitting complete after {iIter} iterations")
        self._log(f"  Entropy:     {entropy:.5f}")
        self._log(f"  chi2/dof:    {chi2:.6f}")
        self._log(f"  Test:        {test:.6f}")
        self._log(f"  Mean bright: {meanBright:.7f}")
        self._log(f"  Mean |B|:    {meanMag:.4f} G")

        # --- Magnetic field energy analysis (after Donati/Lehmann-Petit) ---
        mag_energy = {}
        if par.fitMag == 1 and magGeom.nTot > 0:
            _l = magGeom.l
            _m = magGeom.m
            _lTerm = _l / (_l + 1)
            _m0 = (_m == 0).astype(float)

            def _E(coeff):
                e = 0.5 * coeff * np.conj(coeff)
                e += _m0 * 0.25 * (coeff**2 + np.conj(coeff)**2)
                return np.real(e)

            alphaEs = _E(magGeom.alpha)
            betaEs = _E(magGeom.beta) * _lTerm
            gammaEs = _E(magGeom.gamma) * _lTerm

            totE = float(np.sum(alphaEs) + np.sum(betaEs) + np.sum(gammaEs))
            totEpol = float(np.sum(alphaEs) + np.sum(betaEs))
            totEtor = float(np.sum(gammaEs))

            if totE > 0:
                # per-l energy buckets
                Epol_l = {}
                Etor_l = {}
                for i in range(magGeom.nTot):
                    li = int(_l[i])
                    Epol_l[li] = Epol_l.get(
                        li, 0.0) + float(alphaEs[i] + betaEs[i])
                    Etor_l[li] = Etor_l.get(li, 0.0) + float(gammaEs[i])

                # axisymmetric
                totEaxi = float(np.sum((alphaEs + betaEs + gammaEs) * _m0))
                polEaxi = float(np.sum((alphaEs + betaEs) * _m0))
                torEaxi = float(np.sum(gammaEs * _m0))

                mag_energy = {
                    "totE": totE,
                    "pct_pol": totEpol / totE,
                    "pct_tor": totEtor / totE,
                    "pct_axi": totEaxi / totE,
                    "pct_pol_axi": polEaxi / totEpol if totEpol > 0 else 0.0,
                    "pct_tor_axi": torEaxi / totEtor if totEtor > 0 else 0.0,
                    "pol_by_l": {
                        str(li):
                        (Epol_l.get(li, 0.0) / totEpol if totEpol > 0 else 0.0)
                        for li in range(1, magGeom.nl + 1)
                    },
                    "tor_by_l": {
                        str(li):
                        (Etor_l.get(li, 0.0) / totEtor if totEtor > 0 else 0.0)
                        for li in range(1, magGeom.nl + 1)
                    },
                }

                self._log("")
                self._log("Magnetic Field Energy Analysis")
                self._log("-" * 44)
                self._log(
                    f"  Poloidal:        {mag_energy['pct_pol']:7.3%}  (% tot)"
                )
                self._log(
                    f"  Toroidal:        {mag_energy['pct_tor']:7.3%}  (% tot)"
                )
                self._log(
                    f"  Axisymmetric:    {mag_energy['pct_axi']:7.3%}  (% tot)"
                )
                if totEpol > 0:
                    self._log(
                        f"  Pol axisym:      {mag_energy['pct_pol_axi']:7.3%}  (% pol)"
                    )
                if totEtor > 0:
                    self._log(
                        f"  Tor axisym:      {mag_energy['pct_tor_axi']:7.3%}  (% tor)"
                    )
                self._log("  Poloidal by order:")
                for li in range(1, min(magGeom.nl + 1, 6)):
                    pct = mag_energy['pol_by_l'].get(str(li), 0.0)
                    if totEpol > 0 and pct > 0:
                        _lnames = {
                            1: 'dipole',
                            2: 'quadrupole',
                            3: 'octopole',
                            4: 'l=4',
                            5: 'l=5'
                        }
                        self._log(
                            f"    {_lnames.get(li, f'l={li}'):15s} {pct:7.3%}  (% pol)"
                        )
                self._log("  Toroidal by order:")
                for li in range(1, min(magGeom.nl + 1, 6)):
                    pct = mag_energy['tor_by_l'].get(str(li), 0.0)
                    if totEtor > 0 and pct > 0:
                        self._log(f"    l={li:<13d} {pct:7.3%}  (% tor)")

        chi_aim = par.chiTarget * float(constMem.nDataTotIV)
        converged = ((chi2 * constMem.nDataTotIV <= chi_aim * 1.001) and
                     (test < par.test_aim) if par.fixedEntropy == 0 else
                     (entropy >= par.ent_aim * 1.001) and
                     (test < par.test_aim))

        # --- Build ZDIResult ----------------------------------------------
        from pipeline.result import ZDIResult

        _c_kms = 2.99792458e5

        syn_profiles: list[dict] = []
        for nPhase, spec in enumerate(setSynSpec):
            vel = ((spec.wl - lineData.wl0[0]) / lineData.wl0[0] * _c_kms +
                   par.velRs[nPhase]).tolist()
            syn_profiles.append({
                "phase": float(par.cycleList[nPhase]),
                "vel": vel,
                "I_mod": spec.IIc.tolist(),
                "V_mod": spec.VIc.tolist(),
            })

        obs_profiles: list[dict] = []
        for nPhase, obs in enumerate(obsSet):
            obs_profiles.append({
                "phase": float(par.cycleList[nPhase]),
                "vel": obs.wl.tolist(),
                "I_obs": obs.specI.tolist(),
                "I_sig": obs.specIsig.tolist(),
                "V_obs": obs.specV.tolist(),
                "V_sig": obs.specVsig.tolist(),
            })

        def _coeff_list(arr) -> list:
            return (arr.tolist() if hasattr(arr, "tolist") else list(arr))

        result = ZDIResult(
            iterations=iIter,
            entropy=float(entropy),
            chi2=float(chi2),
            test=float(test),
            converged=bool(converged),
            bright_map=briMap.bright.copy(),
            mag_coeffs={
                "alpha": _coeff_list(magGeom.alpha),
                "beta": _coeff_list(magGeom.beta),
                "gamma": _coeff_list(magGeom.gamma),
            },
            synthetic_profiles=syn_profiles,
            observed_profiles=obs_profiles,
            light_curve_synthetic=(lc_syn if lc_syn is not None else None),
            metadata={
                "config_path":
                str(par.config_path),
                "mean_bright":
                float(meanBright),
                "mean_bright_diff":
                float(meanBrightDiff),
                "mean_mag":
                float(meanMag),
                "kL_fitted":
                float(lineData.str[0]),
                "clat":
                sGrid.clat.tolist(),
                "lon":
                sGrid.long.tolist(),
                "lc_jdates":
                (lc_jdates.tolist() if lc_jdates is not None else []),
                "lc_obs_flux":
                (lc_obs_flux.tolist() if lc_obs_flux is not None else []),
                "mag_energy":
                mag_energy,
            },
        )
        return result
