"""core.geometry — 恒星表面几何计算子包。

公开接口
--------
StarGrid（或 starGrid）
    球形/扁球形恒星表面网格，含法向量、面积、重力昏暗等。
BatchVisibleGrid
    多相位批量可见性几何（向量化，无 Python 循环）。
_SinglePhaseView
    BatchVisibleGrid 单相位切片，与旧版 visibleGrid 接口兼容。
build_batch_visible_grid
    从参数对象一步构造 BatchVisibleGrid。
calcOmegaClat, getCyclesClat, calcVelDiffrotFactor
    差分自转辅助函数。

向后兼容
--------
``from core.geometryStellar import starGrid, ...`` 仍可使用；
``core/geometryStellar.py`` 已变为薄 shim，转发到本子包。
"""
from core.geometry.differential_rotation import (  # noqa: F401
    calcOmegaClat, getCyclesClat, calcVelDiffrotFactor,
)
from core.geometry.stellar_grid import starGrid  # noqa: F401
from core.geometry.visibility import (  # noqa: F401
    _SinglePhaseView, BatchVisibleGrid, build_batch_visible_grid,
)

__all__ = [
    "starGrid",
    "BatchVisibleGrid",
    "_SinglePhaseView",
    "build_batch_visible_grid",
    "calcOmegaClat",
    "getCyclesClat",
    "calcVelDiffrotFactor",
]
