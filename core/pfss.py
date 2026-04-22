"""core/pfss.py — Potential Field Source Surface (PFSS) extrapolation.

将 ZDI 反演所得磁场球谐系数作为输入，利用 pfsspy 计算 PFSS 外推，
并追踪三维磁力线，最终返回 Plotly 3D 交互图的 JSON 字典。

主要公开接口
-----------
compute_pfss_plotly(mag_coeffs, *, nlon, nclat, nrho, rss,
                    n_lat_seeds, n_lon_seeds) -> dict
    接收 mag_coeffs（与 pipeline 输出格式一致的 alpha/beta/gamma 列表），
    返回 {"data": [...Plotly traces], "layout": {...}} 供前端直接传给 Plotly.newPlot。

依赖
----
    pfsspy>=1.0   sunpy>=4.0   astropy>=5.0
以上依赖为可选；若未安装，调用时会抛出含安装提示的 ImportError。
core/ 内部不得 import api/ 或 frontend/ 的任何符号。
"""

from __future__ import annotations

import math
import numpy as np


def _check_imports():
    """延迟检查可选依赖，未安装时给出清晰的错误信息。"""
    missing = []
    try:
        import pfsspy  # noqa: F401
    except ImportError:
        missing.append("pfsspy")
    try:
        import sunpy  # noqa: F401
    except ImportError:
        missing.append("sunpy")
    try:
        import astropy  # noqa: F401
    except ImportError:
        missing.append("astropy")
    if missing:
        raise ImportError(f"PFSS 计算需要以下依赖包，请先安装：{', '.join(missing)}\n"
                          "  pip install " + " ".join(missing))


def _to_complex(lst: list) -> np.ndarray:
    """将 [[real, imag], ...] 格式还原为复数数组（与 pipeline 输出格式一致）。"""
    return np.array([r + 1j * i for r, i in lst], dtype=complex)


def _build_br_map(mag_coeffs: dict, nlon: int, nclat: int) -> np.ndarray:
    """从球谐系数重建恒星表面径向磁场 Br 图。

    Parameters
    ----------
    mag_coeffs : dict
        含 "alpha", "beta", "gamma" 键，每个值为 [[real, imag], ...] 列表，
        与 ZDIPipeline.run() 输出的 result["mag_coeffs"] 格式相同。
    nlon, nclat : int
        经度/余纬度格点数（建议 360 × 180）。

    Returns
    -------
    Br : ndarray, shape (nclat, nlon)
        单位 Gauss，经度方向从 -π 到 π，余纬度方向从 0 到 π。
        注意：pfsspy 要求纬度从南极到北极排列，因此返回的 Br 已按
        latitude 方向翻转（shape 行 = nclat 从 -90° 到 +90°）。
    """
    import core.magneticGeom as mG  # noqa: PLC0415

    lon_list = np.linspace(-np.pi, np.pi, nlon, endpoint=False)
    lat_list = np.linspace(-np.pi / 2, np.pi / 2, nclat)
    lon_2d, lat_2d = np.meshgrid(lon_list, lat_list)
    # ZDI 使用余纬度 (colatitude)：clat = pi/2 - lat
    clat_2d = np.pi / 2 - lat_2d

    _nTot = len(mag_coeffs["alpha"])
    _nl = int(round((-3.0 + math.sqrt(9.0 + 8.0 * _nTot)) / 2.0))

    mag_geom = mG.magSphHarmonics(_nl)
    mag_geom.alpha = _to_complex(mag_coeffs["alpha"])
    mag_geom.beta = _to_complex(mag_coeffs["beta"])
    mag_geom.gamma = _to_complex(mag_coeffs["gamma"])
    mag_geom.initMagGeom(clat_2d.flatten(), lon_2d.flatten())

    vec_b = mag_geom.getAllMagVectors()  # (3, nclat*nlon)
    Br_flat = vec_b[0]

    # reshape → (nclat, nlon)，经度方向翻转使其从 0 → 2π（pfsspy 需要）
    Br = np.reshape(Br_flat, (nclat, nlon))[:, ::-1]
    return Br


def _run_pfss(Br: np.ndarray, nrho: int, rss: float):
    """用 pfsspy 运行 PFSS 外推。

    Parameters
    ----------
    Br : ndarray, shape (nclat, nlon)
        恒星表面径向磁场，单位 Gauss。行方向对应纬度从 -90° 到 90°。
    nrho, rss : int, float
        PFSS 网格径向格点数与源面半径（单位：恒星半径）。

    Returns
    -------
    pfss_out : pfsspy.Output
    input_map : sunpy.map.Map  （用于绘制表面图）
    """
    import pfsspy
    import sunpy.map
    from astropy.time import Time

    nclat, nlon = Br.shape
    header = pfsspy.utils.carr_cea_wcs_header(Time("2000-1-1"), (nlon, nclat))
    input_map = sunpy.map.Map((Br, header))
    pfss_in = pfsspy.Input(input_map, nrho, rss)
    pfss_out = pfsspy.pfss(pfss_in)
    return pfss_out, input_map


def _trace_field_lines(pfss_out, n_lat_seeds: int, n_lon_seeds: int):
    """在均匀种子点网格上追踪磁力线。

    Returns
    -------
    field_lines : list of pfsspy field line objects
    """
    import astropy.constants as const
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from pfsspy import tracing

    r = 1.01 * const.R_sun
    lat = np.linspace(-np.pi / 2, np.pi / 2, n_lat_seeds, endpoint=False)
    lon = np.linspace(0, 2 * np.pi, n_lon_seeds, endpoint=False)
    lat_2d, lon_2d = np.meshgrid(lat, lon, indexing="ij")
    lat_seeds = lat_2d.ravel() * u.rad
    lon_seeds = lon_2d.ravel() * u.rad

    seeds = SkyCoord(lon_seeds, lat_seeds, r, frame=pfss_out.coordinate_frame)
    # FortranTracer requires the compiled `streamtracer` extension.
    # Fall back to PythonTracer when streamtracer is unavailable.
    try:
        tracer = tracing.FortranTracer()
    except Exception:
        tracer = tracing.PythonTracer()
    return tracer.trace(seeds, pfss_out)


def _build_sphere_surface(Br: np.ndarray) -> dict:
    """构建带 Br 着色的球面 Plotly Surface trace（单位半径）。"""
    nclat, nlon = Br.shape
    lats = np.linspace(-90, 90, nclat)
    lons = np.linspace(-180, 180, nlon)

    x = np.outer(np.cos(np.radians(lons)), np.cos(np.radians(lats)))
    y = np.outer(np.sin(np.radians(lons)), np.cos(np.radians(lats)))
    z = np.outer(np.ones(nlon), np.sin(np.radians(lats)))

    abs_max = float(np.max(np.abs(Br))) or 1.0

    return {
        "type": "surface",
        "x": x.tolist(),
        "y": y.tolist(),
        "z": z.tolist(),
        "surfacecolor": Br.T.tolist(),
        "colorscale": "RdBu",
        "cmin": -abs_max,
        "cmax": abs_max,
        "showscale": True,
        "colorbar": {
            "title": "Br (G)",
            "titlefont": {
                "size": 11
            },
            "tickfont": {
                "size": 9
            },
            "len": 0.6,
            "thickness": 12,
        },
        "opacity": 1.0,
        "name": "Surface Br",
        "showlegend": False,
    }


def _build_field_line_traces(field_lines) -> list:
    """将追踪到的磁力线转换为 Plotly Scatter3d trace 列表。"""
    import astropy.constants as const

    # 极性颜色映射：开放正极 → 红，开放负极 → 蓝，封闭 → 灰
    _color_map = {0: "rgba(160,160,160,0.7)", -1: "#3b82f6", 1: "#ef4444"}

    traces = []
    for fl in field_lines:
        color = _color_map.get(fl.polarity, "rgba(160,160,160,0.5)")
        coords = fl.coords
        coords.representation_type = "cartesian"
        x_line = (coords.x / const.R_sun).value.tolist()
        y_line = (coords.y / const.R_sun).value.tolist()
        z_line = (coords.z / const.R_sun).value.tolist()
        traces.append({
            "type": "scatter3d",
            "x": x_line,
            "y": y_line,
            "z": z_line,
            "mode": "lines",
            "line": {
                "color": color,
                "width": 2
            },
            "showlegend": False,
            "hoverinfo": "skip",
        })
    return traces


def compute_pfss_plotly(
    mag_coeffs: dict,
    *,
    nlon: int = 360,
    nclat: int = 180,
    nrho: int = 30,
    rss: float = 2.5,
    n_lat_seeds: int = 18,
    n_lon_seeds: int = 8,
) -> dict:
    """从 ZDI 磁场系数计算 PFSS 外推并返回 Plotly 3D 图 JSON。

    Parameters
    ----------
    mag_coeffs : dict
        ZDI 反演结果中的 ``result["mag_coeffs"]``，含 "alpha", "beta", "gamma"，
        每项为 [[real, imag], ...] 列表。
    nlon, nclat : int
        表面磁场重建网格经度/余纬度格点数，默认 360 × 180。
    nrho : int
        PFSS 径向网格层数（越大越精确，越慢），默认 30。
    rss : float
        源面半径，单位恒星半径，默认 2.5。
    n_lat_seeds : int
        磁力线追踪纬度方向种子数，默认 18。
    n_lon_seeds : int
        磁力线追踪经度方向种子数，默认 8。

    Returns
    -------
    dict
        Plotly JSON 字典 ``{"data": [...traces], "layout": {...}}``，
        可直接传给前端的 ``Plotly.newPlot(div, resp.data, resp.layout)``。

    Raises
    ------
    ImportError
        若 pfsspy / sunpy / astropy 未安装。
    ValueError
        若 mag_coeffs 格式不正确或为空。
    """
    _check_imports()

    if not mag_coeffs or not mag_coeffs.get("alpha"):
        raise ValueError("mag_coeffs 为空或缺少 alpha 键，请先完成 ZDI 反演。")

    # 1. 重建表面 Br 图
    Br = _build_br_map(mag_coeffs, nlon, nclat)

    # 2. 运行 PFSS
    pfss_out, _ = _run_pfss(Br, nrho, rss)

    # 3. 追踪磁力线
    field_lines = _trace_field_lines(pfss_out, n_lat_seeds, n_lon_seeds)

    # 4. 构建 Plotly traces
    sphere_trace = _build_sphere_surface(Br)
    fl_traces = _build_field_line_traces(field_lines)

    # 5. 组装 layout（暗色主题，隐藏坐标轴）
    layout = {
        "template": "plotly_dark",
        "scene": {
            "xaxis": {
                "visible": False,
                "range": [-rss * 1.1, rss * 1.1]
            },
            "yaxis": {
                "visible": False,
                "range": [-rss * 1.1, rss * 1.1]
            },
            "zaxis": {
                "visible": False,
                "range": [-rss * 1.1, rss * 1.1]
            },
            "aspectmode": "cube",
            "camera": {
                "eye": {
                    "x": 1.4,
                    "y": 0.0,
                    "z": 0.7
                }
            },
            "bgcolor": "#0d1117",
        },
        "margin": {
            "r": 10,
            "b": 10,
            "l": 10,
            "t": 40
        },
        "showlegend": False,
        "title": {
            "text": f"PFSS 磁力线外推  (rss={rss} R★, nrho={nrho})",
            "font": {
                "size": 14
            },
            "x": 0.5,
        },
        "paper_bgcolor": "#0d1117",
    }

    # 添加源面半径辅助球（半透明线框，仅 wireframe 效果）
    _theta = np.linspace(0, np.pi, 20)
    _phi = np.linspace(0, 2 * np.pi, 40)
    _t, _p = np.meshgrid(_theta, _phi)
    _xs = rss * np.sin(_t) * np.cos(_p)
    _ys = rss * np.sin(_t) * np.sin(_p)
    _zs = rss * np.cos(_t)
    source_surface = {
        "type":
        "surface",
        "x":
        _xs.tolist(),
        "y":
        _ys.tolist(),
        "z":
        _zs.tolist(),
        "colorscale": [[0, "rgba(100,100,200,0.08)"],
                       [1, "rgba(100,100,200,0.08)"]],
        "showscale":
        False,
        "opacity":
        0.12,
        "name":
        f"Source surface r={rss}",
        "showlegend":
        False,
        "hoverinfo":
        "skip",
        "surfacecolor":
        np.zeros_like(_xs).tolist(),
    }

    data = [sphere_trace, source_surface] + fl_traces
    return {"data": data, "layout": layout}
