"""core.fitting — ZDI 主拟合循环和结果输出函数。

本模块从原 ``core/mainFuncs.py`` 中提取了核心运算职责：

- ``mainFittingLoop``  : MEM 迭代拟合循环（entropy/chi^2 双目标，含取消机制）
- ``saveModelProfs``   : 将合成轮廓写入输出文件（批量或单独 LSD 格式）
- ``saveObsUsed``      : 将观测轮廓写入汇总文件

``core/mainFuncs.py`` 保留向后兼容 shim，直接从本模块 re-export 以上符号。
"""
from __future__ import annotations

import threading
from typing import Optional

import numpy as np

_C_KMS = 2.99792458e5  # speed of light in km/s


def _chi2_components(par, nDataTot: int, fmodel: np.ndarray, data: np.ndarray,
                     sig2: np.ndarray) -> tuple[float, float, float]:
    """Return weighted chi^2 contributions for Stokes I, Stokes V, and extra terms."""
    n_use = min(fmodel.shape[0], data.shape[0], sig2.shape[0])
    if n_use <= 0:
        return 0.0, 0.0, 0.0

    chi2_i = 0.0
    chi2_v = 0.0
    idx = 0

    if par.fitBri == 1:
        end_i = min(idx + nDataTot, n_use)
        if end_i > idx:
            chi2_i = float(
                np.sum((fmodel[idx:end_i] - data[idx:end_i])**2 /
                       sig2[idx:end_i]))
        idx = end_i

    if par.fitMag == 1:
        end_v = min(idx + nDataTot, n_use)
        if end_v > idx:
            chi2_v = float(
                np.sum((fmodel[idx:end_v] - data[idx:end_v])**2 /
                       sig2[idx:end_v]))
        idx = end_v

    chi2_extra = 0.0
    if n_use > idx:
        chi2_extra = float(
            np.sum((fmodel[idx:n_use] - data[idx:n_use])**2 / sig2[idx:n_use]))

    return chi2_i, chi2_v, chi2_extra


# ---------------------------------------------------------------------------
# 主拟合循环
# ---------------------------------------------------------------------------


def mainFittingLoop(par,
                    lineData,
                    wlSynSet: list,
                    sGrid,
                    briMap,
                    magGeom,
                    listGridView,
                    dMagCart0,
                    setSynSpec: list,
                    coMem,
                    nDataTot: int,
                    Data: np.ndarray,
                    sig2: np.ndarray,
                    allModeldIdV,
                    weightEntropy: np.ndarray,
                    verbose: int = 1,
                    stop_event: Optional[threading.Event] = None,
                    callback=None):
    """MEM 迭代拟合循环。

    该函数在 Stokes I 亮度图和 Stokes V 磁场球谐系数上交替执行 MEM 更新，
    直到达到目标 chi^2（或目标熵），或超过最大迭代次数。

    Parameters
    ----------
    par
        参数对象（ZDIConfig 或 readParamsZDI），需含 chiTarget、fixedEntropy、
        ent_aim、numIterations、test_aim 等属性。
    lineData
        谱线参数对象（``lineData``）。
    wlSynSet : list
        各相位合成波长网格。
    sGrid
        ``starGrid`` 对象。
    briMap
        ``brightMap`` 亮度图（原地更新）。
    magGeom
        ``magSphHarmonics`` 磁场球谐对象（原地更新）。
    listGridView
        ``BatchVisibleGrid`` 或列表，按索引返回单相位几何。
    dMagCart0
        初始磁场导数（``fitMag=0`` 时传 0）。
    setSynSpec : list
        初始化完毕的 ``diskIntProfAndDeriv`` 列表。
    coMem
        ``constantsMEM`` 对象（由 memSimple.constantsMEM 构造）。
    nDataTot : int
        数据总点数。
    Data : (N_data,) ndarray
        观测数据向量（Stokes I + V）。
    sig2 : (N_data,) ndarray
        方差向量。
    allModeldIdV
        响应矩阵（或 0 表示未预计算）。
    weightEntropy : ndarray
        熵权重向量。
    verbose : int, optional
        输出级别：0=静默，1=每迭代打印（默认 1）。
    stop_event : threading.Event, optional
        设置后在下一迭代边界处提前返回（WebUI 取消机制）。

    Returns
    -------
    tuple
        (iIter, entropy, chi2nu, test, meanBright, meanBrightDiff, meanMag)
    """
    import core.mem.zdi_adapter as memSimple  # noqa: PLC0415

    chi_aim = par.chiTarget * float(coMem.nDataTotIV)
    target_aim = par.ent_aim if par.fixedEntropy == 1 else chi_aim

    # UR 模型在 B=0 处 dV/dCoeff=0，响应矩阵退化，MEM 无法建立搜索方向。
    # 若磁场系数全为零且使用 UR 模型，用随机小扰动初始化：
    #   - 量级取 defaultBent * 1e-2（足以使响应矩阵非零，同时远小于 defaultBent，
    #     避免 Hobson-Lasenby 熵梯度在大振幅处趋近于 ±1 形成局部极小值）；
    #   - 随机复相位打破所有系数同相位的对称性，避免 B∥LOS 辐条形奇异初始场。
    if (par.fitMag == 1 and getattr(par, 'line_model_type', 'voigt') == 'unno'
            and np.all(magGeom.alpha == 0) and np.all(magGeom.beta == 0)
            and np.all(magGeom.gamma == 0)):
        rng = np.random.default_rng(seed=42)
        _binit = par.defaultBent * 1e-2
        n = magGeom.nTot
        magGeom.alpha[:] = _binit * (rng.standard_normal(n) +
                                     1j * rng.standard_normal(n))
        magGeom.beta[:] = _binit * (rng.standard_normal(n) +
                                    1j * rng.standard_normal(n))
        magGeom.gamma[:] = _binit * (rng.standard_normal(n) +
                                     1j * rng.standard_normal(n))

    with open('outFitSummary.txt', 'w') as fOutFitSummary:

        # 初始化收敛参数
        Chi2 = chi2nu = 0.0
        chi2I = chi2V = chi2Extra = 0.0
        entropy = 0.0
        test = 1.0
        meanBright = meanBrightDiff = meanMag = 0.0
        iIter = 0
        bConverged = False

        # 主迭代循环
        while iIter < par.numIterations and not bConverged:

            if stop_event is not None and stop_event.is_set():
                if verbose >= 1:
                    _msg = f'Run cancelled after {iIter} iterations.'
                    print(_msg)
                    if callback is not None:
                        callback(_msg)
                return iIter, entropy, chi2nu, test, meanBright, meanBrightDiff, meanMag

            iIter += 1
            if iIter > 1:
                memSimple.unpackImageVector(
                    Image,  # noqa: F821 — first assigned below
                    briMap,
                    magGeom,
                    par.magGeomType,
                    par.fitBri,
                    par.fitMag)

            # 当前相位磁场向量
            vecMagCart = magGeom.getAllMagVectorsCart()

            # ── 逐历元亮度变化修正（可选，基于CTTSzdi2算法）────────────
            from core.brightnessGeom import epoch_brightness_scale  # noqa: PLC0415

            for nObs, _phase in enumerate(par.cycleList):
                sGridView = listGridView[nObs]
                spec = setSynSpec[nObs]

                # 若启用逐历元亮度变化（dT修正），应用相位依赖调制
                if getattr(par, 'epochVarEnabled', False):
                    obs_dT = getattr(par, 'obsDT', None)
                    if obs_dT is not None and len(obs_dT) > nObs:
                        dt_val = float(obs_dT[nObs])
                    else:
                        dt_val = float(getattr(par, 'epochVarDT', 1.0))
                    ccq, _ = epoch_brightness_scale(
                        briMap.bright, _phase, getattr(par, 'epochVarMode', 1),
                        getattr(par, 'epochVarT', 100.0),
                        getattr(par, 'epochVarKK', 10.0),
                        getattr(par, 'epochVarDphi', 0.0), dt_val)
                    # 创建临时亮度图副本，避免污染原始 briMap
                    import copy
                    briMapEpoch = copy.copy(briMap)
                    briMapEpoch.bright = ccq
                else:
                    briMapEpoch = briMap

                spec.updateIntProfDeriv(sGridView, vecMagCart, dMagCart0,
                                        briMapEpoch, lineData, par.calcDI,
                                        par.calcDV)
                spec.convolveIGnumpy(par.instrumentRes)

            # 打包 MEM 输入向量
            allModelIV = memSimple.packModelVector(setSynSpec, par.fitBri,
                                                   par.fitMag)
            Image = memSimple.packImageVector(briMap, magGeom, par.magGeomType,
                                              par.fitBri, par.fitMag)

            if par.calcDI == 1 or par.calcDV == 1:
                allModeldIdV = memSimple.packResponseMatrix(
                    setSynSpec, nDataTot, coMem.npBriMap, magGeom,
                    par.magGeomType, par.calcDI, par.calcDV)

            # 执行一步 MEM 迭代
            entropy, Chi2, test, Image, entStand, entFF, entMag = \
                memSimple.mem_iter(
                    coMem.n1Model, coMem.n2Model, coMem.nModelTot,
                    Image, Data, allModelIV, sig2, allModeldIdV,
                    weightEntropy, par.defaultBright, par.defaultBent,
                    par.maximumBright, target_aim, par.fixedEntropy)

            chi2I, chi2V, chi2Extra = _chi2_components(par, nDataTot,
                                                       allModelIV, Data, sig2)

            # 统计量
            meanBright = (np.sum(briMap.bright * sGrid.area) /
                          np.sum(sGrid.area))
            meanBrightDiff = np.sum(
                np.abs(briMap.bright - par.defaultBright) *
                sGrid.area) / np.sum(sGrid.area)
            absMagCart = np.sqrt(vecMagCart[0, :]**2 + vecMagCart[1, :]**2 +
                                 vecMagCart[2, :]**2)
            meanMag = np.sum(absMagCart * sGrid.area) / np.sum(sGrid.area)

            # 收敛判断
            if par.fixedEntropy == 1:
                bConverged = (entropy >= par.ent_aim * 1.001
                              and test < par.test_aim and iIter > 2)
            else:
                bConverged = (Chi2 <= chi_aim * 1.001 and test < par.test_aim)

            chi2nu = Chi2 / max(float(coMem.nDataTotIV), 1.0)

            _chi2_parts = []
            if par.fitBri == 1:
                _chi2_parts.append(f'chi2I {chi2I:10.4f}')
            if par.fitMag == 1:
                _chi2_parts.append(f'chi2V {chi2V:10.4f}')
            if chi2Extra > 0.0:
                _chi2_parts.append(f'chi2X {chi2Extra:10.4f}')
            _chi2_txt = ('  ' + '  '.join(_chi2_parts)) if _chi2_parts else ''

            _summary = ('it {:3n}  entropy {:13.5f}  chi2 {:10.6f}{}  '
                        'Test {:10.6f} meanBright {:10.7f} '
                        'meanSpot {:10.7f} meanMag {:10.4f}'.format(
                            iIter, entropy, chi2nu, _chi2_txt, test,
                            meanBright, meanBrightDiff, meanMag))
            if verbose == 1:
                print(_summary)
            if callback is not None:
                callback(_summary)
            fOutFitSummary.write(_summary + '\n')

            if verbose == 1 and bConverged:
                _conv_msg = 'Success: sufficiently small value of Test achieved'
                print(_conv_msg)
                if callback is not None:
                    callback(_conv_msg)

        # 零迭代时仍计算模型诊断
        if iIter == 0:
            allModelIV = memSimple.packModelVector(setSynSpec, par.fitBri,
                                                   par.fitMag)
            chi2I, chi2V, chi2Extra = _chi2_components(par, nDataTot,
                                                       allModelIV, Data, sig2)
            chi2nu = (np.sum((allModelIV - Data)**2 / sig2) /
                      max(float(coMem.nDataTotIV), 1.0))
            Image = memSimple.packImageVector(briMap, magGeom, par.magGeomType,
                                              par.fitBri, par.fitMag)

            if par.fitBri == 1 or par.fitMag == 1:
                entropy, _, _, _, _, _, _, _ = memSimple.get_s_grads(
                    coMem.n1Model, coMem.n2Model, coMem.nModelTot, Image,
                    weightEntropy, par.defaultBright, par.defaultBent,
                    par.maximumBright)

            meanBright = (np.sum(briMap.bright * sGrid.area) /
                          np.sum(sGrid.area))
            meanBrightDiff = np.sum(
                np.abs(briMap.bright - par.defaultBright) *
                sGrid.area) / np.sum(sGrid.area)
            vecMagCart = magGeom.getAllMagVectorsCart()
            absMagCart = np.sqrt(vecMagCart[0, :]**2 + vecMagCart[1, :]**2 +
                                 vecMagCart[2, :]**2)
            meanMag = np.sum(absMagCart * sGrid.area) / np.sum(sGrid.area)

        # 末尾打印（非逐迭代模式或零迭代）
        if (verbose != 1 and par.numIterations > 0) or (verbose == 1
                                                        and iIter == 0):
            _chi2_parts = []
            if par.fitBri == 1:
                _chi2_parts.append(f'chi2I {chi2I:10.4f}')
            if par.fitMag == 1:
                _chi2_parts.append(f'chi2V {chi2V:10.4f}')
            if chi2Extra > 0.0:
                _chi2_parts.append(f'chi2X {chi2Extra:10.4f}')
            _chi2_txt = ('  ' + '  '.join(_chi2_parts)) if _chi2_parts else ''

            _summary = ('it {:3n}  entropy {:13.5f}  chi2 {:10.6f}{}  '
                        'Test {:10.6f} meanBright {:10.7f} '
                        'meanSpot {:10.7f} meanMag {:10.4f}'.format(
                            iIter, entropy, chi2nu, _chi2_txt, test,
                            meanBright, meanBrightDiff, meanMag))
            print(_summary)
            if callback is not None:
                callback(_summary)
            fOutFitSummary.write(_summary + '\n')

    return iIter, entropy, chi2nu, test, meanBright, meanBrightDiff, meanMag


# ---------------------------------------------------------------------------
# 结果输出
# ---------------------------------------------------------------------------


def saveModelProfs(par, setSynSpec: list, lineData, saveName: str) -> None:
    """将合成轮廓写入汇总文件并各自保存为 LSD 格式 .model 文件。

    Parameters
    ----------
    par
        参数对象（需有 ``cycleList``, ``velRs``, ``fnames``）。
    setSynSpec : list
        ``diskIntProfAndDeriv`` 列表。
    lineData
        谱线参数对象（需有 ``wl0``）。
    saveName : str
        汇总输出文件路径。
    """
    with open(saveName, 'w') as fOutputSpec:
        for nPhase, spec in enumerate(setSynSpec):
            wl = ((spec.wl - lineData.wl0) / lineData.wl0 * _C_KMS +
                  par.velRs[nPhase])  # 观测者静止系速度
            fOutputSpec.write('#cycle %f\n' % par.cycleList[nPhase])
            for i in range(spec.wl.shape[0]):
                fOutputSpec.write('%e %e %e\n' %
                                  (wl[i], spec.IIc[i], spec.VIc[i]))
            fOutputSpec.write('\n')

            # 同时写入各相位 LSD 格式 .model 文件
            sigmaOut = 1e-8
            with open(par.fnames[nPhase].strip() + '.model', 'w') as fLSD:
                fLSD.write('#synthetic prof cycle %f\n' %
                           par.cycleList[nPhase])
                fLSD.write('%i %i \n' % (spec.wl.shape[0], 6))
                for i in range(spec.wl.shape[0]):
                    fLSD.write('%f %e %e %e %e %e %e\n' %
                               (wl[i], spec.IIc[i], sigmaOut, spec.VIc[i],
                                sigmaOut, 0.0, sigmaOut))


def saveObsUsed(obsSet, fName: str) -> None:
    """将使用的观测轮廓写入文件。

    Parameters
    ----------
    obsSet : iterable of obsProf
        观测轮廓集合。
    fName : str
        输出文件路径。
    """
    with open(fName, 'w') as outObserved:
        for tmpObs in obsSet:
            for i in range(tmpObs.wl.shape[0]):
                outObserved.write(
                    '%e %e %e %e %e\n' %
                    (tmpObs.wl[i], tmpObs.specI[i], tmpObs.specIsig[i],
                     tmpObs.specV[i], tmpObs.specVsig[i]))
            outObserved.write('\n')
