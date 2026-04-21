# core.line_models — 谱线轮廓模型 模块说明

> 代码位置：`core/line_models/`
> 算法基础：Humlicek (1982)；Unno (1956)；Rachkovsky (1967)；Landi Degl'Innocenti & Landolfi (2004)

---

## 1. 基本原理

### 1.1 问题背景

ZDI 合成 LSD 线轮廓的核心步骤是**盘积分（disk integration）**：将恒星表面所有可见格点的局部线轮廓加权求和，权重包含临边昏暗、重力昏暗、投影面积和局部亮度。
每个格点的局部轮廓形状由谱线物理模型决定，其位置由格点的旋转速度 Doppler 移位决定，其偏振分量（Stokes V）由格点磁场的视线分量决定。

本子包提供两种谱线物理模型，通过配置项 `line_model.model_type` 选择：

- **Voigt 弱场模型（默认）**：假设磁场分裂远小于谱线宽度（弱场近似），Stokes V 与视线磁场分量成正比，正反比于 Stokes I 的导数。
- **Unno-Rachkovsky 模型（完整偏振辐射转移）**：基于 Milne-Eddington 大气的 Unno (1956) 解，Stokes I 本身也依赖磁场，适用于强磁场或需要精确偏振转移的场景。

### 1.2 数学模型

#### Voigt 弱场近似

局部 Stokes I（Voigt 轮廓乘以亮度）：
$$
I_{\rm loc}(\lambda) = b_k \cdot \bigl[1 - d \cdot V_H(\lambda - \lambda_{\rm Doppler})\bigr]
$$
其中 $V_H$ 为 Voigt 轮廓（Humlicek 算法），$d$ 为谱线强度，$\lambda_{\rm Doppler} = \lambda_0(1 + v_{\rm rot,proj}/c)$ 为 Doppler 移位后的线心。

局部 Stokes V（弱场近似）：
$$
V_{\rm loc}(\lambda) = -\Delta\lambda_Z \cdot B_{\rm los} \cdot \frac{\partial I_{\rm loc}}{\partial\lambda}
$$
其中 Zeeman 移位 $\Delta\lambda_Z = 4.6686\times10^{-12}\,g_{\rm eff}\,\lambda_0^2$（nm，场强单位 G），$B_{\rm los}$ 为视线方向磁场分量。

**Humlicek (1982) Voigt 函数**（复误差函数 $w_4$，精度约 $10^{-4}$）：
$$
V_H(x) = \mathrm{Re}\bigl[w_4(y - ix)\bigr], \quad x = \frac{\lambda_0 - \lambda}{\Delta\lambda_D},\quad y = \frac{\Delta\lambda_L}{\Delta\lambda_D}
$$
其中 $\Delta\lambda_D$ 为 Doppler（Gaussian）宽度，$\Delta\lambda_L$ 为 Lorentzian 宽度。

#### Unno-Rachkovsky 模型（Milne-Eddington）

在 Milne-Eddington 大气中，偏振辐射转移方程（Landi Degl'Innocenti & Landolfi 2004, Eq. 9.32）的解析解为：

$$
\frac{I}{I_c} = f_I \cdot \frac{1 + \eta_I}{(1+\eta_I)^2 - \eta_U^2 - \eta_Q^2 - \eta_V^2 + 2(\eta_Q\eta_V - \eta_U(1+\eta_I))\Delta} + (1-f_I) \cdot \text{(安静区)}
$$

$$
\frac{V}{I_c} = f_V \cdot \frac{\eta_V + (...)}{(\text{分母})} + \text{(安静区)}
$$

其中吸收系数 $\eta_{I,Q,U,V}$ 由 Voigt 和 Faraday-Voigt 函数的磁场分量组合给出；$f_I$、$f_V$ 为 Stokes I/V 各自的填充因子，$\beta$ 为 Planck 函数坡度。

### 1.3 算法流程

**Voigt 模型（主路径）**

```
lineData.from_parameters(...)      ← 构造谱线参数对象
        │
        ▼
getAllProfDiriv(par, sGrid, ...)    ← 批量为每个相位构造 diskIntProfAndDeriv
        │  （内部：每格点 Doppler 移位 → 每格点 localProfileAndDeriv）
        ▼
diskIntProfAndDeriv.updateProfDeriv()
  ├─ for each phase: 可见格点加权求和 → IIc(vel), VVIc(vel)
  └─ 计算导数：dIIc/d(bright_k), dVVIc/d(mag_coeff)
        │
        ▼
packDataVector / packResponseMatrix (core/mem/zdi_adapter.py)
```

**Unno-Rachkovsky 模型（可选路径）**

相同调用链，将 `lineData` 换为 `lineDataUnno`，`getAllProfDiriv` 换为 `getAllProfDirivUnno`；仪器展宽改用显式卷积（积分后 Gaussian 卷积）。

---

## 2. 模块概述

`core/line_models/` 是 ZDIpy 管线中**谱线局部轮廓计算与盘积分**的独立子包，上游依赖 `core/geometry/`（可见性权重），下游为 `core/mem/zdi_adapter.py` 提供数据向量和响应矩阵。

```
BatchVisibleGrid (几何)
 ├─ vel_rot_proj → Doppler 移位
 ├─ view_angle   → 临边昏暗 / B_los 投影
 └─ proj_area    → 积分权重

brightMap (亮度)         magSphHarmonics (磁场)
     │                         │
     └─────────┬───────────────┘
               ▼
         diskIntProfAndDeriv / diskIntProfAndDerivUnno
           (每相位：局部轮廓 → 加权盘积分)
               │
               ├─ IIc(vel)           ← 合成 Stokes I/Ic
               ├─ VVIc(vel)          ← 合成 Stokes V/Ic
               ├─ dIIc/d(bright_k)   ← 亮度 Jacobian（含可见性、临边昏暗）
               └─ dVVIc/d(mag_coeff) ← 磁场 Jacobian
```

---

## 3. 子模块一览

| 文件 | 主要内容 | 关键类/函数 |
|------|---------|------------|
| `base.py` | 谱线模型抽象基类（新架构） | `LineModel` |
| `voigt.py` | Humlicek Voigt + 弱场 Stokes V 实现 | `VoigtLineModel`, `_humlicek_voigt` |
| `profile.py` | 弱场模型主类（管线使用） | `lineData`, `localProfileAndDeriv`, `diskIntProfAndDeriv`, `getAllProfDiriv` |
| `unno.py` | Unno-Rachkovsky 完整偏振辐射转移 | `lineDataUnno`, `localProfileAndDerivUnno`, `diskIntProfAndDerivUnno`, `getAllProfDirivUnno` |
| `disk_integration.py` | 独立盘积分函数（与模型解耦） | `disk_integrate`, `normalize_by_continuum` |
| `line_utils.py` | 辅助函数：临边昏暗、等效宽度、谱线强度拟合 | `limbDarkening`, `calcSynEW`, `fitLineStrength`, `equivWidComp2` |
| `__init__.py` | 包级公共 API 统一入口 | 见第 4 节 |

---

## 4. 公共接口（`core.line_models`）

```python
from core.line_models import (
    # 新架构
    LineModel,
    VoigtLineModel,
    disk_integrate,
    normalize_by_continuum,

    # Voigt 弱场（管线级，主路径）
    lineData,
    localProfileAndDeriv,
    diskIntProfAndDeriv,
    getAllProfDiriv,
    explicitConvolution,

    # Unno-Rachkovsky（可选路径）
    lineDataUnno,
    localProfileAndDerivUnno,
    diskIntProfAndDerivUnno,
    getAllProfDirivUnno,

    # 辅助函数
    limbDarkening,
    calcSynEW,
    fitLineStrength,
    equivWidComp2,
)
```

---

## 5. 核心类说明

### 5.1 `LineModel` (`base.py`)

谱线模型抽象基类，规定新架构的标准接口。

#### 抽象方法

```python
# 局部 Stokes I/V 轮廓（不做盘积分）
I, V = model.compute_profile(
    vel_grid,   # (N_vel,)
    B_los,      # (N_cells,) kG
    brightness, # (N_cells,)
    view_angle, # (N_cells,) rad
)
# I, V: (N_cells, N_vel)

# 导数字典
derivs = model.compute_derivatives(
    vel_grid, B_los, dB_los_d_coeff, brightness, view_angle,
)
# keys: 'dI_dBright', 'dV_dBlos', 'dV_dCoeff'
```

---

### 5.2 `lineData` (`profile.py`)

Voigt 弱场模型的谱线参数容器，支持文件读取或直接参数构造。

#### 构造参数（直接参数）

| 参数 | 类型 | 说明 |
|------|------|------|
| `wl0` | `float` | 中心波长（nm） |
| `lineStr` | `float` | 谱线强度（0–1，即谱线深度） |
| `widthGauss` | `float` | Gaussian Doppler 宽度（km/s） |
| `widthLorentz` | `float` | Lorentzian 宽度 / Gaussian 宽度（比值） |
| `landeG` | `float` | 有效 Landé $g$ 因子 |
| `limbDark` | `float` | 临边昏暗系数 $\varepsilon$ |
| `gravDark` | `float` | 重力昏暗系数 |
| `instRes` | `float` | 仪器分辨率 $R$（$< 0$ 表示不施加仪器展宽） |

```python
ld = lineData(
    wl0=630.15, lineStr=0.55, widthGauss=5.0, widthLorentz=0.5,
    landeG=1.2, limbDark=0.6, gravDark=0.0, instRes=65000.0,
)
# 或从 config 对象构造：
ld = lineData.from_parameters(par)
```

---

### 5.3 `diskIntProfAndDeriv` (`profile.py`)

单相位盘积分对象：给定视图几何和亮度/磁场，计算合成 Stokes I/V 轮廓及全部 Jacobian。

```python
spec = diskIntProfAndDeriv(
    visibleGrid,  # _SinglePhaseView（单相位几何）
    lineData,     # lineData 对象
    vGrid,        # (N_vel,) km/s 速度网格
    briMap,       # brightMap
    magGeom,      # magSphHarmonics（已 initMagGeom）
    magGeomType,  # str，磁几何类型
    fitBri=1,     # 是否计算亮度导数
    fitMag=1,     # 是否计算磁场导数
)
```

**主要属性（调用 `updateProfDeriv()` 后）**

| 属性 | 形状 | 含义 |
|------|------|------|
| `IIc` | `(N_vel,)` | 归一化合成 Stokes I/Ic |
| `VVIc` | `(N_vel,)` | 归一化合成 Stokes V/Ic |
| `dIIc` | `(N_cells, N_vel)` | $\partial(I/I_c)/\partial b_k$（亮度 Jacobian） |
| `dVVIc` | `(n_mag_coeff, N_vel)` | $\partial(V/I_c)/\partial \alpha_{\ell m}$ 等（磁场 Jacobian） |

---

### 5.4 `lineDataUnno` (`unno.py`)

Unno-Rachkovsky 模型参数容器，在 `lineData` 参数基础上额外携带：

| 额外参数 | 类型 | 说明 |
|----------|------|------|
| `beta` | `float` | Planck 函数梯度 $\beta$（$\leq 0$ 时由 `limbDark` 自动推导） |
| `fillingFactorI` | `float` | Stokes I 填充因子 $f_I$（磁区面积占比） |
| `fillingFactorV` | `float` | Stokes V 填充因子 $f_V$ |

```python
ldu = lineDataUnno.from_parameters(par)
# 对应 pipeline 中 par.line_model_type == 'unno' 的分支
```

---

### 5.5 `disk_integrate` (`disk_integration.py`)

独立于具体线模型的通用盘积分函数。

```python
integrated = disk_integrate(
    profiles,       # (N_phases, N_cells, N_vel) 或 (N_cells, N_vel)
    proj_areas,     # (N_phases, N_cells) 或 (N_cells,)
    limb_weights,   # 同 proj_areas 形状
    inst_fwhm_pix=0.0,  # 仪器 FWHM（像素），0 = 不卷积
)
# -> (N_phases, N_vel) 或 (N_vel,)
```

---

### 5.6 辅助函数 (`line_utils.py`)

```python
# 线性临边昏暗权重
w = limbDarkening(coeffEta, angle)    # (N_cells,) 或标量

# 合成轮廓等效宽度（梯形积分）
ew = calcSynEW(spec)                  # float，单位同 spec.wl

# 通过等效宽度自动估计谱线强度（scipy.optimize 驱动）
fitLineStrength(par, setSynSpec, lineData, obsSet)
```

---

## 6. 文件结构

```
core/
└── line_models/                ← 谱线轮廓模型子包（本文档描述范围）
    ├── __init__.py             ← 公共 API 统一入口
    ├── base.py                 ← LineModel 抽象基类（新架构）
    ├── voigt.py                ← VoigtLineModel（Humlicek Voigt + 弱场 V）
    ├── profile.py              ← lineData / diskIntProfAndDeriv / getAllProfDiriv
    ├── unno.py                 ← Unno-Rachkovsky 完整偏振辐射转移
    ├── disk_integration.py     ← disk_integrate（与模型解耦的独立盘积分）
    └── line_utils.py           ← limbDarkening / calcSynEW / fitLineStrength
```

---

## 7. 参考文献

- Humlicek, J. (1982). Optimized computation of the Voigt and complex probability functions. *JQSRT*, 27, 437–444.
- Unno, W. (1956). Line formation of a normal Zeeman triplet. *PASJ*, 8, 108–125.
- Rachkovsky, D. N. (1967). Magnetic rotation effects in spectral lines. *Izv. Krymsk. Astrofiz. Obs.*, 36, 3–11.
- Landi Degl'Innocenti, E. & Landolfi, M. (2004). *Polarization in Spectral Lines*. Springer. (Eqs. 9.32, 9.110–9.112)
- Gray, D. F. (2005). *The Observation and Analysis of Stellar Photospheres*, 3rd ed. Cambridge University Press.
- Donati, J.-F. & Brown, S. F. (1997). Zeeman-Doppler imaging of active stars. V. *A&A*, 326, 1135–1142.
