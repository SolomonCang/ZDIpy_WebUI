# core.brightnessGeom — 亮度图 模块说明

> 代码位置：`core/brightnessGeom.py`
> 算法基础：Collier Cameron (1992)；Gray (2005)

---

## 1. 基本原理

### 1.1 问题背景

恒星表面并非均匀明亮——活动区（如黑子）的温度低于周围光球，导致局部辐射强度减弱。
ZDI 亮度图反演从 Stokes I（总强度）轮廓中提取这种表面亮度不均匀性：每个网格元拥有一个归一化亮度值 $b_k$（均匀恒星为 1.0，黑子区域 $b_k < 1$）。

同时，临边昏暗（limb darkening）和重力昏暗（gravity darkening）使不同视角位置的格点呈现不同权重，盘积分时须逐格点加权。

### 1.2 数学模型

**线性临边昏暗**（Gray 2005, Eq. 17.11）：
$$
I(\mu) = 1 - \varepsilon(1 - \mu), \quad \mu = \cos\theta_{\rm view}
$$
其中 $\varepsilon$ 为临边昏暗系数（`limb_darkening` 配置项）。

**重力昏暗**（von Zeipel 定律，适用于快速旋转的扁球形恒星）：
$$
I_{\rm grav} \propto g_{\rm eff}^{\,\beta_{\rm gd}}
$$
其中有效重力 $g_{\rm eff}$ 由 `starGrid.gravityDarkening()` 预计算并存储，对球形星该因子恒为 1。

**加权投影亮度**（单相位）：
$$
B_k^{(j)} = b_k \cdot \bigl[1 - \varepsilon(1-\mu_k^{(j)})\bigr] \cdot g_k \cdot dA_k^{(j)} \cdot \mathbb{1}[\text{visible}]
$$
其中 $dA_k^{(j)} = dA_k\,\mu_k^{(j)}$ 为投影面积，$g_k$ 为重力昏暗因子。

整合为矩阵形式（批量相位）：
$$
\mathbf{B}^{(j,k)} = \mathbf{b}_k \cdot w^{(j,k)}_{\rm limb} \cdot g_k \cdot A^{(j,k)}_{\rm proj} \cdot V^{(j,k)}
$$
形状 $(N_{\rm phases}, N_{\rm cells})$，中间无 Python 循环。

### 1.3 算法流程

1. 初始化 `brightMap`：每个格点 $b_k = b_{\rm default}$（均匀亮度，通常 1.0）；或从文件读入。
2. MEM 迭代过程中，`mem_iter` 更新 `brightMap.bright` 数组。
3. 每次更新合成谱时，调用 `brightMap.projected(visible_batch, limb_dark, grav_dark)` 获取 $(N_{\rm phases}, N_{\rm cells})$ 加权亮度矩阵。
4. 该矩阵作为 `diskIntProfAndDeriv.updateProfDeriv()` 的亮度权重，参与 Stokes I 积分。

---

## 2. 模块概述

`core/brightnessGeom.py` 是 ZDIpy 管线中**恒星表面亮度图的存储与加权投影**模块，位于几何层（`core/geometry/`）之上，为谱线盘积分（`core/line_models/`）提供每格点亮度权重。

```
starGrid (几何)
  ├─ sGrid.clat, sGrid.long    ─→  brightMap.__init__
  └─ sGrid.gravityDarkening()  ─→  brightMap.projected(grav_dark)

BatchVisibleGrid (可见性)
  ├─ visible_batch.view_angle  ─→  limb_weight 计算
  └─ visible_batch.proj_area   ─→  投影面积

brightMap
  └─ projected() → (N_phases, N_cells) 加权亮度
           │
           ▼
  diskIntProfAndDeriv (盘积分)
```

---

## 3. 子模块一览

`core/brightnessGeom.py` 为单文件模块，包含以下公共符号：

| 符号 | 类型 | 主要内容 |
|------|------|---------|
| `brightMap` | 类 | 亮度图存储与加权投影 |
| `SetupBrightMap` | 函数 | 初始化工厂（从文件或默认值） |
| `saveMap` | 函数 | 将亮度图写入文本文件 |
| `readMap` | 函数 | 从文本文件读入亮度图（含坐标校验） |

---

## 4. 公共接口（`core.brightnessGeom`）

```python
from core.brightnessGeom import (
    brightMap,
    SetupBrightMap,
    saveMap,
    readMap,
)
```

---

## 5. 核心类说明

### 5.1 `brightMap` (`brightnessGeom.py`)

像素级亮度图：存储每个网格元的归一化亮度值，并计算包含临边/重力昏暗的加权投影。

#### 构造参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `clat` | `(N_cells,) ndarray` | 格点余纬坐标（来自 `starGrid.clat`） |
| `lon` | `(N_cells,) ndarray` | 格点经度坐标（来自 `starGrid.long`） |

初始化后 `bright[k] = 1.0`（均匀亮度）。

#### 主要属性

| 属性 | 形状 | 含义 |
|------|------|------|
| `bright` | `(N_cells,)` | 各格点归一化亮度（MEM 反演的目标变量之一） |
| `clat` | `(N_cells,)` | 余纬坐标（弧度） |
| `lon` | `(N_cells,)` | 经度坐标（弧度） |
| `npt` | `int` | 总格点数 |

#### 主要方法

```python
# 加权投影亮度（核心调用）
weighted = bri_map.projected(
    visible_batch,       # BatchVisibleGrid
    limb_dark_coeff,     # float，临边昏暗系数 ε
    grav_dark_factor,    # (N_cells,)，重力昏暗因子
)
# -> (N_phases, N_cells)，不可见格点为 0
```

```python
# 放置圆形黑子（测试用途）
bri_map.makeRoundSpot(
    clat,       # float，黑子中心余纬
    lon,        # float，黑子中心经度
    radius,     # float，黑子角半径（弧度）
    brightness, # float，黑子亮度值
)
```

---

### 5.2 `SetupBrightMap` (`brightnessGeom.py`)

从文件或默认常数初始化亮度图的工厂函数。

```python
bri_map = SetupBrightMap(
    sGrid,               # starGrid 实例
    initBrightFromFile,  # int，1 = 从文件读取，0 = 使用默认值
    initBrightFile,      # str，初始化文件路径（initBrightFromFile=1 时使用）
    defaultBright,       # float，默认亮度值（通常 1.0）
    verbose=1,           # int，是否打印加载信息
)
# -> brightMap 实例
```

---

### 5.3 `saveMap` / `readMap`

```python
# 保存：每行写入 "clat  lon  bright"
saveMap(bri_map, fileName)

# 读取：校验坐标系与输入网格一致后返回 brightMap，不一致时返回 None
bri_map = readMap(fileName, userClat, userLon, verbose=1)
# -> brightMap | None
```

---

## 6. 文件结构

```
core/
└── brightnessGeom.py     ← 亮度图模块（本文档描述范围，单文件）
```

---

## 7. 参考文献

- Collier Cameron, A. (1992). Stellar surface imaging using total light observations. In *ASP Conf. Ser.*, 26, 265.
- Gray, D. F. (2005). *The Observation and Analysis of Stellar Photospheres*, 3rd ed. Cambridge University Press. (临边昏暗 Eq. 17.11)
- von Zeipel, H. (1924). The radiative equilibrium of a rotating system of gaseous masses. *MNRAS*, 84, 665–683. (重力昏暗定律)
- Donati, J.-F. et al. (2006). Spectropolarimetric observations of the Herbig Ae/Be star HD 190073. *MNRAS*, 370, 629–644.
