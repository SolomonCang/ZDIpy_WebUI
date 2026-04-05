"""
验证 BatchVisibleGrid 与旧逐相位 visibleGrid 逻辑的数值一致性。

旧 visibleGrid 代码已迁移到本文件内作为内联参考实现；不需要保留原始类。
运行方式：
    python scripts/check_geom_consistency.py
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from math import pi
from core.geometry import (
    starGrid,
    BatchVisibleGrid,
    getCyclesClat,
    calcVelDiffrotFactor,
)


# ---------------------------------------------------------------------------
# Reference single-phase implementation (reproduces old visibleGrid logic)
# ---------------------------------------------------------------------------
def _ref_single_phase(sg: starGrid, inclination: float,
                      cycle_at_clat: np.ndarray, period: float,
                      dOmega: float) -> dict:
    """Reproduce the former visibleGrid.__init__ as a pure function."""
    view_long = -cycle_at_clat * 2.0 * pi
    view_x = np.sin(inclination) * np.cos(view_long)
    view_y = np.sin(inclination) * np.sin(view_long)
    view_z = np.cos(inclination) * np.ones(view_long.shape)
    v_view = np.array([view_x, view_y, view_z])
    len_view = np.linalg.norm(v_view, axis=0)
    v_view_cart = v_view / len_view  # (3, N_cells)

    v_normals = sg.getCartesianNormals()
    dot_prod = np.einsum('ij,ij->j', v_normals, v_view)
    cos_view = dot_prod / len_view
    view_angle = np.arccos(cos_view)

    visible = np.zeros(sg.numPoints, dtype=int)
    for i in range(sg.numPoints):
        if 0.0 < view_angle[i] < pi * 0.5:
            visible[i] = 1

    proj_area = sg.area * cos_view

    v_vel = sg.GetCartesianRotVel()
    vel_diffrot = calcVelDiffrotFactor(period, dOmega, sg.clat)
    v_vel = v_vel * vel_diffrot
    vel_rot_proj = np.einsum('ij,ij->j', -v_view_cart, v_vel)

    return dict(
        vViewCart=v_view_cart,
        visible=visible,
        projArea=proj_area,
        velRotProj=vel_rot_proj,
        viewAngle=view_angle,
    )


# ---------------------------------------------------------------------------
# Test parameters
# ---------------------------------------------------------------------------
N_RINGS = 20
PERIOD = 5.0
D_OMEGA = 0.05
J_DATE_REF = 2456000.0
J_DATES = np.array([2456000.0, 2456001.25, 2456002.5, 2456003.75])
INC_RAD = np.deg2rad(60.0)

print(f"Building starGrid with {N_RINGS} rings …")
sg = starGrid(N_RINGS, period=PERIOD, verbose=0)

cycles = getCyclesClat(PERIOD, D_OMEGA, J_DATES, J_DATE_REF, sg.clat)
# (N_phases, N_cells)

print("Building BatchVisibleGrid …")
batch = BatchVisibleGrid(sg, INC_RAD, cycles, PERIOD, D_OMEGA)

# ---------------------------------------------------------------------------
# Compare phase-by-phase
# ---------------------------------------------------------------------------
ATOL_BOOL = 0
ATOL_FLOAT = 1e-12

print(f"\nComparing {len(J_DATES)} phases …")
for i, jd in enumerate(J_DATES):
    ref = _ref_single_phase(sg, INC_RAD, cycles[i], PERIOD, D_OMEGA)
    sv = batch[i]

    assert np.array_equal(sv.visible, ref['visible']), \
        f"visible mismatch at phase {i}"
    assert np.allclose(sv.projArea, ref['projArea'], atol=ATOL_FLOAT), \
        f"projArea mismatch at phase {i}: max_err={np.max(np.abs(sv.projArea - ref['projArea']))}"
    assert np.allclose(sv.velRotProj, ref['velRotProj'], atol=ATOL_FLOAT), \
        f"velRotProj mismatch at phase {i}: max_err={np.max(np.abs(sv.velRotProj - ref['velRotProj']))}"  # noqa: E501
    assert np.allclose(sv.vViewCart, ref['vViewCart'], atol=ATOL_FLOAT), \
        f"vViewCart mismatch at phase {i}"
    print(
        f"  Phase {i} (JD={jd:.2f}) … OK  "
        f"(projArea ΔMax={np.max(np.abs(sv.projArea - ref['projArea'])):.2e}, "
        f"velRot ΔMax={np.max(np.abs(sv.velRotProj - ref['velRotProj'])):.2e})"
    )

# ---------------------------------------------------------------------------
# Shape checks
# ---------------------------------------------------------------------------
N_P, N_C = len(J_DATES), sg.numPoints
assert batch.visible.shape == (N_P, N_C), "visible shape wrong"
assert batch.proj_area.shape == (N_P, N_C), "proj_area shape wrong"
assert batch.vel_rot_proj.shape == (N_P, N_C), "vel_rot_proj shape wrong"
assert batch.v_view_cart.shape == (N_P, 3, N_C), "v_view_cart shape wrong"
print("\nShape checks … OK")

# ---------------------------------------------------------------------------
# Spherical-star flux conservation  (solid-body rotation, fine grid)
# ---------------------------------------------------------------------------
# Note: the finite-grid discretization at limb cells means exact conservation
# only holds in the limit of infinite resolution.  We use a 50-ring grid
# (rel_var ~1e-5) to verify the batch code produces physically reasonable flux.
sg50 = starGrid(50, verbose=0)
cycles_sb = getCyclesClat(PERIOD, 0.0, J_DATES, J_DATE_REF, sg50.clat)
batch_sb = BatchVisibleGrid(sg50, INC_RAD, cycles_sb, PERIOD, 0.0)
proj_sum_sb = np.sum(batch_sb.proj_area * batch_sb.visible, axis=1)
rel_var_sb = np.std(proj_sum_sb) / np.mean(proj_sum_sb)
assert rel_var_sb < 5e-5, f"Flux not conserved (dOmega=0, 50 rings): rel_var={rel_var_sb:.2e}"
print(
    f"Flux conservation (dOmega=0, 50 rings) … OK  (rel_var={rel_var_sb:.2e})")

# ---------------------------------------------------------------------------
# Diffrot cache
# ---------------------------------------------------------------------------
d1 = sg.get_diffrot_factor(PERIOD, D_OMEGA)
d2 = sg.get_diffrot_factor(PERIOD, D_OMEGA)
assert d1 is d2, "get_diffrot_factor should return the cached object"
print("Diffrot cache … OK")

print("\n✓ All geometry consistency checks passed.")
