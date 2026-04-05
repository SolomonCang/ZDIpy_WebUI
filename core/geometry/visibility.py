"""core.geometry.visibility — 批量多相位可见性几何。

提供：
- _SinglePhaseView  : 单相位视图（与旧版 visibleGrid 接口兼容）
- BatchVisibleGrid  : 多相位批量可见性，所有数组形状 (N_phases, N_cells)
- build_batch_visible_grid : 从参数对象快速构造 BatchVisibleGrid
"""
from __future__ import annotations

import numpy as np

from core.geometry.differential_rotation import getCyclesClat


class _SinglePhaseView:
    """BatchVisibleGrid 的单相位切片，提供与旧版 visibleGrid 兼容的接口。

    下游代码（diskIntProfAndDeriv、fitLineStrength 等）通过属性名访问
    几何量，无需修改即可兼容新旧两套接口。
    """
    __slots__ = (
        'vViewCart',
        'velRotProj',
        'projArea',
        'visible',
        'viewAngle',
        'gravityDarkening',
    )

    def __init__(
        self,
        v_view_cart: np.ndarray,  # (3, N_cells)
        vel_rot_proj: np.ndarray,  # (N_cells,)
        proj_area: np.ndarray,  # (N_cells,)
        visible: np.ndarray,  # (N_cells,), bool
        view_angle: np.ndarray,  # (N_cells,)
        gravity_darkening_fn,
    ) -> None:
        self.vViewCart = v_view_cart
        self.velRotProj = vel_rot_proj
        self.projArea = proj_area
        self.visible = visible.astype(int)
        self.viewAngle = view_angle
        self.gravityDarkening = gravity_darkening_fn


class BatchVisibleGrid:
    """多相位批量可见性几何对象。

    所有内部数组均为形状 (N_phases, N_cells)，无 Python 相位循环。
    单相位兼容视图通过 ``batch[i]`` 返回 :class:`_SinglePhaseView`。

    Parameters
    ----------
    star_grid : starGrid
        恒星表面格点对象（含预计算法向量和旋转速度）。
    inclination : float
        恒星自转轴倾角（弧度）。
    cycles_at_clat : (N_phases, N_cells) ndarray
        各相位各余纬处的旋转圈数（由 getCyclesClat 生成）。
    period : float
        自转周期（天），用于差分自转因子计算。
    dOmega : float
        赤道-极角速度差 (rad/day)。
    """
    def __init__(
        self,
        star_grid,  # starGrid 实例
        inclination: float,
        cycles_at_clat: np.ndarray,  # (N_phases, N_cells)
        period: float,
        dOmega: float,
    ) -> None:
        self._star_grid = star_grid
        n_phases, n_cells = cycles_at_clat.shape

        # 各格点的视向方向（考虑差分自转引起的经度偏移）
        # view_long shape: (N_phases, N_cells)
        view_long = -cycles_at_clat * 2.0 * np.pi

        si = np.sin(inclination)
        ci = np.cos(inclination)

        view_x = si * np.cos(view_long)  # (N_phases, N_cells)
        view_y = si * np.sin(view_long)  # (N_phases, N_cells)
        view_z = np.full((n_phases, n_cells), ci)  # (N_phases, N_cells)

        # v_view_cart shape: (N_phases, 3, N_cells)
        v_view_cart = np.stack([view_x, view_y, view_z], axis=1)
        # 归一化（理论上已是单位向量，数值归一化更稳健）
        len_view = np.sqrt(np.sum(v_view_cart**2, axis=1, keepdims=True))
        v_view_cart = v_view_cart / len_view

        # 可见性和投影面积
        # normals: (3, N_cells) 广播为 (N_phases, 3, N_cells)
        cos_view = np.sum(v_view_cart * star_grid._normals[np.newaxis, :, :],
                          axis=1)  # (N_phases, N_cells)

        visible = cos_view > 0.0  # bool
        proj_area = star_grid.area[
            np.newaxis, :] * cos_view  # (N_phases, N_cells)
        view_angle = np.arccos(np.clip(cos_view, -1.0,
                                       1.0))  # (N_phases, N_cells)

        # 投影旋转速度
        diffrot = star_grid.get_diffrot_factor(period, dOmega)  # (N_cells,)
        rot_vel = star_grid._rot_vel * diffrot[np.newaxis, :]  # (3, N_cells)
        vel_rot_proj = -np.sum(v_view_cart * rot_vel[np.newaxis, :, :],
                               axis=1)  # (N_phases, N_cells)

        # 公开属性
        self.v_view_cart: np.ndarray = v_view_cart  # (N_phases, 3, N_cells)
        self.visible: np.ndarray = visible  # (N_phases, N_cells) bool
        self.proj_area: np.ndarray = proj_area  # (N_phases, N_cells)
        self.view_angle: np.ndarray = view_angle  # (N_phases, N_cells)
        self.vel_rot_proj: np.ndarray = vel_rot_proj  # (N_phases, N_cells)

    def __getitem__(self, phase_idx: int) -> _SinglePhaseView:
        """返回第 phase_idx 相位的 :class:`_SinglePhaseView`（向后兼容接口）。"""
        return _SinglePhaseView(
            v_view_cart=self.v_view_cart[phase_idx],  # (3, N_cells)
            vel_rot_proj=self.vel_rot_proj[phase_idx],  # (N_cells,)
            proj_area=self.proj_area[phase_idx],  # (N_cells,)
            visible=self.visible[phase_idx],  # (N_cells,) bool
            view_angle=self.view_angle[phase_idx],  # (N_cells,)
            gravity_darkening_fn=self._star_grid.gravityDarkening,
        )

    def __len__(self) -> int:
        return self.visible.shape[0]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


def build_batch_visible_grid(par, s_grid) -> BatchVisibleGrid:
    """从参数对象和格点一步构造 BatchVisibleGrid。

    Parameters
    ----------
    par
        包含 ``period``、``dOmega``、``jDates``、``jDateRef``、``incRad``
        属性的参数对象（ZDIConfig 或 readParamsZDI 均可）。
    s_grid : starGrid
        恒星表面格点对象。

    Returns
    -------
    BatchVisibleGrid
        所有相位的批量可见性几何，形状 (N_phases, N_cells)。
    """
    cycles_at_clat = getCyclesClat(par.period, par.dOmega, par.jDates,
                                   par.jDateRef,
                                   s_grid.clat)  # (N_phases, N_cells)
    return BatchVisibleGrid(
        star_grid=s_grid,
        inclination=par.incRad,
        cycles_at_clat=cycles_at_clat,
        period=par.period,
        dOmega=par.dOmega,
    )
