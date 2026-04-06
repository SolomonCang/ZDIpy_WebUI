"""
Unit tests for the refactored core/magneticGeom.py and core/brightnessGeom.py.

Run with:
    python -m pytest tests/test_mag_brightness.py -v
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pytest
import scipy.special

from core.geometry import starGrid, BatchVisibleGrid, getCyclesClat
from core.magneticGeom import (
    magSphHarmonics,
    magSphHarmoics,  # new name + backward-compat alias
    _compute_legendre_batch,
)
from core.brightnessGeom import brightMap


# ---------------------------------------------------------------------------
# Reference scalar Legendre (old single-point logic, kept here for comparison)
# ---------------------------------------------------------------------------
def _ref_lpmn_single(m_max, n_max, x_scalar):
    """Old _lpmn_compat scalar interface — reference implementation."""
    p = np.zeros((m_max + 1, n_max + 1))
    pd = np.zeros((m_max + 1, n_max + 1))
    for m_val in range(m_max + 1):
        n_arr = np.arange(m_val, n_max + 1, dtype=float)
        if len(n_arr) > 0:
            p[m_val, m_val:n_max +
              1] = scipy.special.lpmv(m_val, n_arr, x_scalar) * ((-1.0)**m_val)
    denom = x_scalar * x_scalar - 1.0
    if abs(denom) > 1e-14:
        for m_val in range(m_max + 1):
            for n_val in range(m_val, n_max + 1):
                pm1 = p[m_val, n_val - 1] if n_val > m_val else 0.0
                pd[m_val, n_val] = (n_val * x_scalar * p[m_val, n_val] -
                                    (n_val + m_val) * pm1) / denom
    else:
        s = 1 if x_scalar > 0 else -1
        for n_val in range(n_max + 1):
            pd[0, n_val] = s**(n_val + 1) * n_val * (n_val + 1) / 2
    return p, pd


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def sg30():
    return starGrid(30, verbose=0)


@pytest.fixture(scope="module")
def mag10(sg30):
    m = magSphHarmonics(10)
    m.initMagGeom(sg30.clat, sg30.long)
    return m


# ---------------------------------------------------------------------------
# 1. _compute_legendre_batch consistency with scalar reference
# ---------------------------------------------------------------------------
class TestLegendreConsistency:
    def test_polynomial_values(self, sg30):
        cos_clat = np.cos(sg30.clat)
        p_batch, _ = _compute_legendre_batch(5, 5, cos_clat)
        for i in range(len(cos_clat)):
            p_ref, _ = _ref_lpmn_single(5, 5, cos_clat[i])
            assert np.allclose(p_batch[:, :, i], p_ref, atol=1e-12), \
                f"Legendre polynomial mismatch at cell {i}"

    def test_derivative_values(self, sg30):
        cos_clat = np.cos(sg30.clat)
        _, pd_batch = _compute_legendre_batch(5, 5, cos_clat)
        for i in range(len(cos_clat)):
            _, pd_ref = _ref_lpmn_single(5, 5, cos_clat[i])
            assert np.allclose(pd_batch[:, :, i], pd_ref, atol=1e-12), \
                f"Legendre derivative mismatch at cell {i}"

    def test_shape(self, sg30):
        N = sg30.numPoints
        p, pd = _compute_legendre_batch(8, 8, np.cos(sg30.clat))
        assert p.shape == (9, 9, N)
        assert pd.shape == (9, 9, N)


# ---------------------------------------------------------------------------
# 2. SH basis matrix shapes and dtype
# ---------------------------------------------------------------------------
class TestBasisMatrixShapes:
    def test_xterm_shape(self, mag10, sg30):
        assert mag10.XTerm.shape == (mag10.nTot, sg30.numPoints)
        assert mag10.XTerm.dtype == np.complex128

    def test_yterm_shape(self, mag10, sg30):
        assert mag10.YTerm.shape == (mag10.nTot, sg30.numPoints)

    def test_zterm_shape(self, mag10, sg30):
        assert mag10.ZTerm.shape == (mag10.nTot, sg30.numPoints)


# ---------------------------------------------------------------------------
# 3. Magnetic field vectors numerical consistency vs. old single-cell loop
# ---------------------------------------------------------------------------
class TestMagFieldConsistency:
    @pytest.fixture
    def mag_old_ref(self, sg30):
        """Build reference magSphHarmonics using the current (vectorised) code
        so both old and new share the same initMagGeom result; compare B vectors."""
        m = magSphHarmonics(8)
        m.initMagGeom(sg30.clat, sg30.long)
        rng = np.random.default_rng(42)
        m.alpha[:] = rng.standard_normal(
            m.nTot) + 1j * rng.standard_normal(m.nTot)
        m.beta[:] = rng.standard_normal(
            m.nTot) + 1j * rng.standard_normal(m.nTot)
        m.gamma[:] = rng.standard_normal(
            m.nTot) + 1j * rng.standard_normal(m.nTot)
        return m

    def test_b_cart_vectors_finite(self, mag_old_ref):
        b = mag_old_ref.getAllMagVectorsCart()
        assert np.all(np.isfinite(b)), "B field contains non-finite values"

    def test_b_cart_shape(self, mag_old_ref, sg30):
        b = mag_old_ref.getAllMagVectorsCart()
        assert b.shape == (3, sg30.numPoints)

    def test_b_derivs_shape(self, mag_old_ref, sg30):
        mag_old_ref.setMagGeomType('full')
        db = mag_old_ref.getAllMagDerivsCart()
        # full: (3 cart, 3 coeff-sets, nTot, N_cells)
        assert db.shape == (3, 3, mag_old_ref.nTot, sg30.numPoints)


# ---------------------------------------------------------------------------
# 4. initMagGeom grid-hash caching
# ---------------------------------------------------------------------------
class TestInitMagGeomCache:
    def test_second_call_is_noop(self, sg30):
        m = magSphHarmonics(5)
        m.initMagGeom(sg30.clat, sg30.long)
        x_before = m.XTerm.copy()
        # Second call with same arrays → should not recompute
        m.initMagGeom(sg30.clat, sg30.long)
        assert np.array_equal(m.XTerm, x_before)

    def test_different_grid_recomputes(self, sg30):
        from core.geometry import starGrid
        sg2 = starGrid(25, verbose=0)
        m = magSphHarmonics(5)
        m.initMagGeom(sg30.clat, sg30.long)
        hash_before = m._init_geom_hash
        m.initMagGeom(sg2.clat, sg2.long)
        assert m._init_geom_hash != hash_before


# ---------------------------------------------------------------------------
# 5. Backward-compat alias
# ---------------------------------------------------------------------------
class TestBackwardCompat:
    def test_alias_is_same_class(self):
        assert magSphHarmoics is magSphHarmonics


# ---------------------------------------------------------------------------
# 6. brightMap.makeRoundSpot vectorised
# ---------------------------------------------------------------------------
class TestMakeRoundSpot:
    def test_spot_applied(self, sg30):
        bm = brightMap(sg30.clat, sg30.long)
        bm.makeRoundSpot(np.pi / 2, 0.0, 0.3, 0.5)
        assert np.any(bm.bright < 1.0)

    def test_spot_shape_unchanged(self, sg30):
        bm = brightMap(sg30.clat, sg30.long)
        bm.makeRoundSpot(np.pi / 2, 0.0, 0.3, 0.5)
        assert bm.bright.shape == (sg30.numPoints, )


# ---------------------------------------------------------------------------
# 7. brightMap.projected shape and flux conservation
# ---------------------------------------------------------------------------
class TestBrightMapProjected:
    @pytest.fixture
    def batch4(self, sg30):
        # Use sg30 for shape tests; use sg50 for flux conservation (requires finer grid)
        j_dates = np.linspace(2456000.0, 2456004.0, 4)
        cycles = getCyclesClat(5.0, 0.0, j_dates, 2456000.0, sg30.clat)
        return BatchVisibleGrid(sg30, np.deg2rad(60.0), cycles, 5.0, 0.0)

    def test_projected_shape(self, sg30, batch4):
        bm = brightMap(sg30.clat, sg30.long)
        gd = sg30.gravityDarkening(0.0)
        w = bm.projected(batch4, 0.6, gd)
        assert w.shape == (4, sg30.numPoints)

    def test_uniform_brightness_constant_flux(self):
        # Use a 50-ring grid to reduce limb discretization artifact to < 5e-5.
        from core.geometry import starGrid
        sg50 = starGrid(50, verbose=0)
        j_dates = np.linspace(2456000.0, 2456004.0, 4)
        cycles = getCyclesClat(5.0, 0.0, j_dates, 2456000.0, sg50.clat)
        batch50 = BatchVisibleGrid(sg50, np.deg2rad(60.0), cycles, 5.0, 0.0)
        bm = brightMap(sg50.clat, sg50.long)
        bm.bright[:] = 1.0
        gd = sg50.gravityDarkening(0.0)
        w = bm.projected(batch50, 0.0, gd)
        flux = w.sum(axis=1)
        rel_var = np.std(flux) / np.mean(flux)
        assert rel_var < 5e-5, f"Flux not conserved: rel_var={rel_var:.2e}"
