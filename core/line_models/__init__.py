"""core.line_models — 谱线轮廓模块（新旧架构统一出口）。

新架构（推荐）
--------------
- ``LineModel``       : 抽象基类 (core.line_models.base)
- ``VoigtLineModel``  : Humlicek Voigt 弱场模型 (core.line_models.voigt)
- ``disk_integrate``  : 独立盘积分函数 (core.line_models.disk_integration)

向后兼容接口（pipeline 级别）
------------------------------
- ``lineData``, ``localProfileAndDeriv``, ``diskIntProfAndDeriv``,
  ``getAllProfDiriv`` — 来自 core.line_models.profile
- ``limbDarkening``, ``calcSynEW``, ``fitLineStrength``,
  ``equivWidComp2``   — 来自 core.line_models.line_utils

``core/lineprofileVoigt.py`` 现为薄 shim，从本子包 re-export 以上符号。
"""
from core.line_models.base import LineModel  # noqa: F401
from core.line_models.voigt import VoigtLineModel  # noqa: F401
from core.line_models.disk_integration import (  # noqa: F401
    disk_integrate, normalize_by_continuum,
)

# ---------------------------------------------------------------------------
# 向后兼容接口：从新位置 re-export 旧名称
# ---------------------------------------------------------------------------
from core.line_models.profile import (  # noqa: F401
    lineData, localProfileAndDeriv, diskIntProfAndDeriv, getAllProfDiriv,
    explicitConvolution,
)
from core.line_models.line_utils import (  # noqa: F401
    limbDarkening, calcSynEW, fitLineStrength, equivWidComp2,
)

__all__ = [
    # 新架构
    "LineModel",
    "VoigtLineModel",
    "disk_integrate",
    "normalize_by_continuum",
    # 向后兼容
    "lineData",
    "localProfileAndDeriv",
    "diskIntProfAndDeriv",
    "getAllProfDiriv",
    "explicitConvolution",
    "limbDarkening",
    "calcSynEW",
    "fitLineStrength",
    "equivWidComp2",
]
