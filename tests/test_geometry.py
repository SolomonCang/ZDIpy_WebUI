"""
Unit tests for the refactored core/geometryStellar.py (BatchVisibleGrid).

Run with:
    python -m pytest tests/test_geometry.py -v
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pytest
from math import pi

from core.geometry import (
    starGrid,
    BatchVisibleGrid,
    _SinglePhaseView,
    getCyclesClat,
    calcVelDiffrotFactor,
)


# ---------------------------------------------------------------------------
# Reference single-phase implementation (same as in check_geom_consistency.py)
# ---------------------------------------------------------------------------
def _ref_single_phase(sg, inclination, cycle_at_clat, period, dOmega):
    view_long = -cycle_at_clat * 2.0 * pi
    view_x = np.sin(inclination) * np.cos(view_long)
    view_y = np.sin(inclination) * np.sin(view_long)
    view_z = np.cos(inclination) * np.ones(view_long.shape)
    v_view = np.array([view_x, view_y, view_z])
    len_view = np.linalg.norm(v_view, axis=0)
    v_view_cart = v_view / len_view

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
    )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def sg20():
    return starGrid(20, verbose=0)


@pytest.fixture(scope="module")
def j_dates():
    return np.linspace(2456000.0, 2456015.0, 16)


@pytest.fixture(scope="module")
def batch_params():
    return dict(period=5.0,
                dOmega=0.05,
                j_date_ref=2456000.0,
                inclination=np.deg2rad(60.0))


@pytest.fixture(scope="module")
def batch16(sg20, j_dates, batch_params):
    cycles = getCyclesClat(batch_params['period'], batch_params['dOmega'],
                           j_dates, batch_params['j_date_ref'], sg20.clat)
    return BatchVisibleGrid(sg20, batch_params['inclination'], cycles,
                            batch_params['period'], batch_params['dOmega'])


# ---------------------------------------------------------------------------
# 1. Shape correctness
# ---------------------------------------------------------------------------
class TestShapes:
    def test_visible_shape(self, batch16, sg20, j_dates):
        assert batch16.visible.shape == (len(j_dates), sg20.numPoints)

    def test_proj_area_shape(self, batch16, sg20, j_dates):
        assert batch16.proj_area.shape == (len(j_dates), sg20.numPoints)

    def test_vel_rot_proj_shape(self, batch16, sg20, j_dates):
        assert batch16.vel_rot_proj.shape == (len(j_dates), sg20.numPoints)

    def test_v_view_cart_shape(self, batch16, sg20, j_dates):
        assert batch16.v_view_cart.shape == (len(j_dates), 3, sg20.numPoints)


# ---------------------------------------------------------------------------
# 2. Numerical consistency with reference implementation
# ---------------------------------------------------------------------------
class TestNumericalConsistency:
    @pytest.mark.parametrize("phase_idx", range(4))
    def test_visible(self, sg20, j_dates, batch_params, batch16, phase_idx):
        cycles = getCyclesClat(batch_params['period'], batch_params['dOmega'],
                               j_dates, batch_params['j_date_ref'], sg20.clat)
        ref = _ref_single_phase(sg20, batch_params['inclination'],
                                cycles[phase_idx], batch_params['period'],
                                batch_params['dOmega'])
        sv = batch16[phase_idx]
        assert np.array_equal(sv.visible, ref['visible'])

    @pytest.mark.parametrize("phase_idx", range(4))
    def test_proj_area(self, sg20, j_dates, batch_params, batch16, phase_idx):
        cycles = getCyclesClat(batch_params['period'], batch_params['dOmega'],
                               j_dates, batch_params['j_date_ref'], sg20.clat)
        ref = _ref_single_phase(sg20, batch_params['inclination'],
                                cycles[phase_idx], batch_params['period'],
                                batch_params['dOmega'])
        sv = batch16[phase_idx]
        assert np.allclose(sv.projArea, ref['projArea'], atol=1e-12)

    @pytest.mark.parametrize("phase_idx", range(4))
    def test_vel_rot_proj(self, sg20, j_dates, batch_params, batch16,
                          phase_idx):
        cycles = getCyclesClat(batch_params['period'], batch_params['dOmega'],
                               j_dates, batch_params['j_date_ref'], sg20.clat)
        ref = _ref_single_phase(sg20, batch_params['inclination'],
                                cycles[phase_idx], batch_params['period'],
                                batch_params['dOmega'])
        sv = batch16[phase_idx]
        assert np.allclose(sv.velRotProj, ref['velRotProj'], atol=1e-12)

    @pytest.mark.parametrize("phase_idx", range(4))
    def test_v_view_cart(self, sg20, j_dates, batch_params, batch16,
                         phase_idx):
        cycles = getCyclesClat(batch_params['period'], batch_params['dOmega'],
                               j_dates, batch_params['j_date_ref'], sg20.clat)
        ref = _ref_single_phase(sg20, batch_params['inclination'],
                                cycles[phase_idx], batch_params['period'],
                                batch_params['dOmega'])
        sv = batch16[phase_idx]
        assert np.allclose(sv.vViewCart, ref['vViewCart'], atol=1e-12)


# ---------------------------------------------------------------------------
# 3. Flux conservation for a spherical star with uniform brightness
# ---------------------------------------------------------------------------
class TestFluxConservation:
    def test_projected_area_constant(self, sg20, j_dates, batch_params):
        # Finite-grid discretization at limb cells limits conservation precision.
        # Use a 50-ring grid (rel_var ~1e-5) with dOmega=0 (solid-body).
        sg50 = starGrid(50, verbose=0)
        cycles_sb = getCyclesClat(batch_params['period'], 0.0, j_dates,
                                  batch_params['j_date_ref'], sg50.clat)
        batch_sb = BatchVisibleGrid(sg50, batch_params['inclination'],
                                    cycles_sb, batch_params['period'], 0.0)
        proj_sum = np.sum(batch_sb.proj_area * batch_sb.visible, axis=1)
        rel_var = np.std(proj_sum) / np.mean(proj_sum)
        assert rel_var < 5e-5, f"Flux varies across phases: rel_var={rel_var:.2e}"

    def test_oblate_grid_area_nonzero_per_cell(self):
        sg_oblate = starGrid(60,
                             period=0.4232,
                             mass=0.66,
                             radiusEq=0.72,
                             verbose=0)
        assert np.all(sg_oblate.area > 0.0)
        assert np.isclose(sg_oblate.area.sum(), 12.945489703083938, atol=1e-10)


# ---------------------------------------------------------------------------
# 4. N_phases=1 vs. N_phases=16 first row
# ---------------------------------------------------------------------------
class TestSinglePhaseConsistency:
    def test_single_vs_batch_first_phase(self, sg20, j_dates, batch_params,
                                         batch16):
        j_single = j_dates[:1]
        cycles_single = getCyclesClat(batch_params['period'],
                                      batch_params['dOmega'], j_single,
                                      batch_params['j_date_ref'], sg20.clat)
        batch1 = BatchVisibleGrid(sg20, batch_params['inclination'],
                                  cycles_single, batch_params['period'],
                                  batch_params['dOmega'])
        assert np.allclose(batch1.visible[0], batch16.visible[0], atol=1e-14)
        assert np.allclose(batch1.proj_area[0],
                           batch16.proj_area[0],
                           atol=1e-14)
        assert np.allclose(batch1.vel_rot_proj[0],
                           batch16.vel_rot_proj[0],
                           atol=1e-14)


# ---------------------------------------------------------------------------
# 5. Diffrot cache returns the same object
# ---------------------------------------------------------------------------
class TestDiffrotCache:
    def test_cache_identity(self, sg20):
        d1 = sg20.get_diffrot_factor(5.0, 0.05)
        d2 = sg20.get_diffrot_factor(5.0, 0.05)
        assert d1 is d2

    def test_cache_invalidation(self, sg20):
        d1 = sg20.get_diffrot_factor(5.0, 0.05)
        d2 = sg20.get_diffrot_factor(5.0, 0.10)
        assert d1 is not d2


# ---------------------------------------------------------------------------
# 6. _SinglePhaseView attribute interface
# ---------------------------------------------------------------------------
class TestSinglePhaseView:
    def test_attributes_exist(self, batch16):
        sv = batch16[0]
        assert isinstance(sv, _SinglePhaseView)
        assert hasattr(sv, 'vViewCart')
        assert hasattr(sv, 'velRotProj')
        assert hasattr(sv, 'projArea')
        assert hasattr(sv, 'visible')
        assert hasattr(sv, 'gravityDarkening')

    def test_gravity_darkening_callable(self, batch16):
        sv = batch16[0]
        gd = sv.gravityDarkening(0.0)
        assert gd is not None
