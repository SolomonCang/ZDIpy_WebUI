# core.geometry — 恒星表面几何 模块说明

> 代码位置：`core/geometry/`
> 算法基础：Domiciano de Souza et al. (2002)；Donati & Brown (1997)

---

## 1. 基本原理

### 1.1 问题背景

ZDI 反演需要将恒星表面离散化为有限个网格元（surface elements），并在每个观测相位计算每个网格元对视线方向的可见性、投影面积和旋转速度投影。
这些几何量是后续谱线盘积分（disk integration）的核心权重因子，直接决定合成谱的精度。

恒星快速自转时表面呈扁球形（Roche 外形），临边昏暗和重力昏暗均会调制局部辐射强度；
差分自转（differential rotation）导致不同纬度的旋转相位积累速率不同，须在相位计算中逐纬度修正。

### 1.2 数学模型

**球面坐标约定**：以余纬 $\theta \in [0,\pi]$（colatitude）和经度 $\phi \in [0,2\pi)$ 定义表面位置，自转轴为 $z$ 轴，视线方向由倾角 $i$ 确定。

**等面积网格**：赤道处经度格点数 $N_{\rm eq} = 2N_{\rm rings}$；余纬圈 $k$ 的经度格点数为：
$$
N_\phi^{(k)} = \mathrm{round}\!\left(N_{\rm eq}\,\sin\theta_k\right), \quad \theta_k = \frac{\pi(k+0.5)}{N_{\rm rings}}
$$
总格点数约为 $N_{\rm cells} \approx 2N_{\rm rings}^2$，各格点面积接近相等。

**球形面积**：格点 $k$ 的面积：
$$
dA_k = d\phi_k\bigl[\cos(\theta_k - \tfrac{d\theta}{2}) - \cos(\theta_k + \tfrac{d\theta}{2})\bigr]
$$

**临边昏暗（线性律）**：
$$
I(\mu) = 1 - \varepsilon(1 - \mu), \quad \mu = \cos\theta_{\rm view}
$$

**差分自转律**（太阳型）：
$$
\Omega(\theta) = \Omega_{\rm eq} - d\Omega\,\cos^2\theta
$$
各余纬处不同的旋转周期导致旋转相位 $\phi(t)$ 随时间积累速率不同。

**投影旋转速度**：格点的旋转速度向量在视线方向的投影（用于 Doppler 移位）：
$$
v_{\rm rot,proj} = -\hat{v}_{\rm rot} \cdot \hat{v}_{\rm los}
$$

### 1.3 算法流程

1. 由 `nRings` 生成等面积余纬-经度网格，计算格点坐标 $(\theta_k, \phi_k)$。
2. 若 `mass > 0, radiusEq > 0`，计算 Roche 扁球形半径 $r(\theta)$ 及重力昏暗因子。
3. 预计算格点的笛卡尔法向量 `_normals` 和单位旋转速度向量 `_rot_vel`。
4. 对每个观测时刻 $t_j$，调用 `getCyclesClat` 计算各格点旋转圈数（考虑差分自转）。
5. 在 `BatchVisibleGrid` 中向量化计算所有相位、所有格点的可见性、投影面积、视角和旋转速度投影（形状均为 $(N_{\rm phases}, N_{\rm cells})$）。

---

## 2. 模块概述

`core/geometry/` 是 ZDIpy 管线中处理**恒星表面离散化与多相位可见性几何**的独立子包，上游不依赖任何物理模型，下游被 `brightnessGeom`、`magneticGeom`、`line_models` 和 `pipeline` 调用。

```
ZDIConfig (观测参数)
       │
       ▼
 starGrid (stellar_grid.py)
 ─ 等面积网格、法向量、面积、旋转速度
       │
       ▼
 getCyclesClat (differential_rotation.py)
 ─ (N_phases, N_cells) 旋转圈数矩阵
       │
       ▼
 BatchVisibleGrid (visibility.py)
 ─ visible / proj_area / view_angle / vel_rot_proj
 ─ 形状 (N_phases, N_cells)，全向量化
       │
       ▼
 diskIntProfAndDeriv / brightMap.projected (下游消费)
```

---

## 3. 子模块一览

| 文件 | 主要内容 | 关键类/函数 |
|------|---------|------------|
| `stellar_grid.py` | 球形/扁球形恒星表面网格，含预计算法向量、面积、重力昏暗 | `starGrid` |
| `visibility.py` | 多相位批量可见性（全向量化），单相位兼容视图 | `BatchVisibleGrid`, `_SinglePhaseView`, `build_batch_visible_grid` |
| `differential_rotation.py` | 太阳型差分自转辅助函数 | `calcOmegaClat`, `getCyclesClat`, `calcVelDiffrotFactor` |
| `__init__.py` | 包级公共 API 统一入口 | 见第 4 节 |

---

## 4. 公共接口（`core.geometry`）

```python
from core.geometry import (
    # 网格
    starGrid,

    # 可见性
    BatchVisibleGrid,
    _SinglePhaseView,
    build_batch_visible_grid,

    # 差分自转
    calcOmegaClat,
    getCyclesClat,
    calcVelDiffrotFactor,
)
```

---

## 5. 核心类说明

### 5.1 `starGrid` (`stellar_grid.py`)

恒星表面格点对象：离散化恒星表面，预计算与相位无关的几何量。

#### 构造参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `numPoints` | `int` | 余纬方向格点（环）数，总格点数约 $2N^2$ |
| `period` | `float` | 自转周期（天），用于扁球形 Roche 半径 |
| `mass` | `float` | 恒星质量（太阳质量），`> 0` 时启用扁球形几何 |
| `radiusEq` | `float` | 赤道半径（太阳半径），`> 0` 时启用扁球形几何 |
| `verbose` | `int` | `0` 静默，`1` 打印网格信息（默认 `1`） |

#### 主要属性

| 属性 | 形状 | 含义 |
|------|------|------|
| `clat` | `(N_cells,)` | 各格点余纬坐标（弧度） |
| `long` | `(N_cells,)` | 各格点经度坐标（弧度） |
| `area` | `(N_cells,)` | 各格点表面积 |
| `_normals` | `(3, N_cells)` | 笛卡尔单位外法向量（预计算） |
| `_rot_vel` | `(3, N_cells)` | 以赤道线速度为单位的旋转速度向量（预计算） |
| `numPoints` | `int` | 总格点数 |

#### 主要方法

```python
area = sGrid.GetSurfaceArea()           # -> (N_cells,) 表面积
normals = sGrid.getCartesianNormals()   # -> (3, N_cells) 法向量
rot_vel = sGrid.GetCartesianRotVel()    # -> (3, N_cells) 旋转速度向量
factor = sGrid.get_diffrot_factor(period, dOmega)  # -> (N_cells,) 差分自转因子
gdark = sGrid.gravityDarkening(grav_coeff)         # -> (N_cells,) 重力昏暗
```

---

### 5.2 `BatchVisibleGrid` (`visibility.py`)

多相位批量可见性几何对象，所有数组形状为 $(N_{\rm phases}, N_{\rm cells})$，无 Python 相位循环。

#### 构造参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `star_grid` | `starGrid` | 恒星表面格点对象 |
| `inclination` | `float` | 自转轴倾角（弧度） |
| `cycles_at_clat` | `(N_phases, N_cells) ndarray` | 各相位各格点旋转圈数 |
| `period` | `float` | 自转周期（天） |
| `dOmega` | `float` | 赤道-极角速度差（rad/day） |

#### 公开属性

| 属性 | 形状 | 含义 |
|------|------|------|
| `visible` | `(N_phases, N_cells)` bool | 格点是否可见 |
| `proj_area` | `(N_phases, N_cells)` | 投影面积（含 $\cos\theta_{\rm view}$） |
| `view_angle` | `(N_phases, N_cells)` | 视角（弧度） |
| `vel_rot_proj` | `(N_phases, N_cells)` | 投影旋转速度（单位：$v\sin i$） |
| `v_view_cart` | `(N_phases, 3, N_cells)` | 笛卡尔视线方向单位向量 |

```python
# 单相位兼容视图（与旧版 visibleGrid 接口兼容）
single = batch[phase_idx]   # -> _SinglePhaseView
single.visible              # (N_cells,) int
single.projArea             # (N_cells,)
single.viewAngle            # (N_cells,)
single.velRotProj           # (N_cells,)
single.vViewCart            # (3, N_cells)
```

---

### 5.3 `build_batch_visible_grid` (`visibility.py`)

从参数对象一步构造 `BatchVisibleGrid` 的便捷工厂函数。

```python
batch = build_batch_visible_grid(par, s_grid)
# par 需含属性：period, dOmega, jDates, jDateRef, incRad
# -> BatchVisibleGrid，形状 (N_phases, N_cells)
```

---

### 5.4 差分自转函数 (`differential_rotation.py`)

```python
# 各余纬角速度 (rad/day)
omega = calcOmegaClat(period, dOmega, clat)   # (N_cells,)

# 各相位各余纬旋转圈数
cycles = getCyclesClat(period, dOmega, jDates, jDateRef, clat)
# -> (N_phases, N_cells)

# 相对赤道的归一化速度因子
factor = calcVelDiffrotFactor(period, dOmega, clat)   # (N_cells,)
```

---

## 6. 文件结构

```
core/
└── geometry/                     ← 恒星表面几何子包（本文档描述范围）
    ├── __init__.py               ← 公共 API 统一入口
    ├── stellar_grid.py           ← starGrid：网格、法向量、面积、重力昏暗
    ├── visibility.py             ← BatchVisibleGrid：批量可见性
    └── differential_rotation.py  ← 差分自转辅助函数
```

---

## 7. 参考文献

- Donati, J.-F. & Brown, S. F. (1997). Zeeman-Doppler imaging of active stars. V. *A&A*, 326, 1135–1142.
- Domiciano de Souza, A. et al. (2002). The spinning-top Be star Achernar. *A&A*, 407, L47–L50.
- Gray, D. F. (2005). *The Observation and Analysis of Stellar Photospheres*, 3rd ed. Cambridge University Press. (临边昏暗律 Eq. 17.11)
- Folsom, C. P. et al. (2018). The large-scale magnetic field and differential rotation of ξ Bootis A. *MNRAS*, 481, 5286–5295.
