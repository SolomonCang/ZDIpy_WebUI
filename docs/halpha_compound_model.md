# core.line_models.halpha — H-alpha 双峰复合 Voigt 模型 模块说明

> 代码位置：`core/line_models/halpha.py`
> 算法基础：Humlicek (1982)；Landi Degl'Innocenti & Landolfi (2004, §5.3 弱场近似)

---

## 1. 基本原理

### 1.1 问题背景

H-alpha（λ₀ = 656.28 nm，$g_{\rm eff} \approx 1.048$）在活跃恒星（T Tauri、dMe 等）的轮廓
通常表现为**宽色球发射 + 窄中心自吸收**的双峰形态：

- **宽发射成分**：色球/过渡区热等离子体（$T \sim 10^4$ K），发射峰半宽约 50–200 km/s；
- **窄自吸收成分**：冷色球/光球层对发射光子的再吸收，抑制线心约 20–80 km/s。

标准单峰 Voigt 模型（`voigt`）和 Milne-Eddington Unno-Rachkovsky 模型（`unno`）均假设单峰
轮廓，无法直接拟合该双峰形态。`halpha_compound` 模型通过**双 Voigt 叠加**参数化地描述
两个成分，同时保留弱场近似对 Stokes V 的解析表达，使响应矩阵可以预计算。

### 1.2 数学模型

#### 局部 Stokes I（双峰复合 Voigt）

$$
I_{\rm loc}(\lambda) = 1 + A_{\rm em} \cdot V_H\!\left(\lambda;\,\lambda_{\rm Dop},\,\Delta\lambda_{D,{\rm em}},\,a_{\rm em}\right)
  - A_{\rm abs} \cdot V_H\!\left(\lambda;\,\lambda_{\rm Dop},\,\Delta\lambda_{D,{\rm abs}},\,a_{\rm abs}\right)
$$

其中 $V_H(u, a) = \operatorname{Re}[w_4(a - iu)]$（Humlicek w4 算法），各符号含义如下：

| 符号 | 含义 |
|------|------|
| $A_{\rm em} > 0$ | 发射成分峰值强度（连续谱归一单位） |
| $A_{\rm abs} \geq 0$ | 自吸收成分峰值深度（= 0 则退化为单峰发射） |
| $\Delta\lambda_{D,c} = w_{G,c}/c \cdot \lambda_0$ | 各分量 Gaussian Doppler 宽度（nm） |
| $a_c$ | 各分量 Lorentz/Gauss 宽度比 |
| $\lambda_{\rm Dop} = \lambda_0(1 + v_{\rm rot}/c)$ | Doppler 移位后的线心 |
| $u = (\lambda_0 - \lambda)/\Delta\lambda_D$ | 归一化坐标（正值对应蓝翼） |

#### 局部 Stokes V（弱场近似，双分量叠加）

H-alpha 的 Zeeman 分裂通常远小于谱线宽度（弱场极限），对两个 Voigt 分量独立应用弱场
近似后线性叠加：

$$
V_{\rm loc}(\lambda) = -\Delta\lambda_Z \cdot f_V \cdot B_{\rm LOS}
  \left[A_{\rm em} \frac{\partial V_{H,{\rm em}}}{\partial\lambda}
  - A_{\rm abs} \frac{\partial V_{H,{\rm abs}}}{\partial\lambda}\right]
$$

其中 $\Delta\lambda_Z = 4.66864 \times 10^{-12}\, g_{\rm eff}\, \lambda_0^2$（nm/G）为 Zeeman
分裂系数，$f_V$ 为 Stokes V 填充因子。

> **注意**：两分量导数符号相反（发射 +、自吸收 −），两者在线心附近部分抵消，
> V 振幅低于纯单成分 Voigt 模型，这是物理上正确的结果。

#### Voigt 函数解析导数

对 Humlicek $w_4(z)$，其中 $z = a - iu$：

$$
\frac{\partial V_H}{\partial \lambda}
= \frac{\partial \operatorname{Re}[w_4]}{\partial \lambda}
= \frac{1}{\Delta\lambda_D} \cdot \frac{\partial}{\partial u}
\operatorname{Re}[w_4] \cdot \left(-1\right)
= \frac{2}{\Delta\lambda_D}\left(u \cdot H(a,u) - a \cdot F(a,u)\right)
$$

其中 $H = \operatorname{Re}[w_4]$（Voigt 函数），$F = \operatorname{Im}[w_4]$（Faraday-Voigt 函数）。
导数完全解析，无需数值差分，与 `core/line_models/profile.py`（Voigt 模型）符号约定一致。

#### 与 Milne-Eddington 框架的对比

Unno-Rachkovsky 模型的线性源函数假设 $S = B_0(1 + \beta\tau)$ 不能描述双层大气（发射 + 自吸收）。
`halpha_compound` 有意退回弱场 Voigt 叠加，用直接参数化换取对双峰形状的灵活控制；
若未来需要完整偏振转移，需引入双层 ME 大气（Landi Degl'Innocenti & Landolfi 2004, §11）。

### 1.3 算法流程

```
1. 从 lineDataHalpha 读取参数
       (wl0, g, A_em, widthGauss_em, widthLorentz_em,
        A_abs, widthGauss_abs, widthLorentz_abs, fV, instRes)

2. 对每个观测相位，由旋转速度构造各格点 Doppler 移位波长网格
       wlCells[N_vel, N_cells] = wl_grid ⊗ (1 + v_rot/c)⁻¹

3. _voigt_component() × 2 — 分别计算发射/自吸收分量
       输入: wl0, width_gauss_kms, width_lorentz, amplitude, wlCells
       输出: comp_I[N_vel, N_cells],  dcomp_dlambda[N_vel, N_cells]

4. 组合为轮廓和导数核（localProfileAndDerivHalpha.__init__）
       Iunscaled  = 1 + I_em − I_abs
       VsigKernel = fV · gCoeff · (dI_em/dλ − dI_abs/dλ)
       gCoeff     = −Δλ_Z · λ₀² · g_eff  [nm/G, 负数]

5. updateProfDeriv(): 盘积分（加权求和）
       I = Σ(contin × Iunscaled) / ΣIcont
       V = Σ(Blos × contin × VsigKernel) / ΣIcont
       dV/dCoeff = einsum(dBlos/dCoeff, VsigKernel × contin) / ΣIcont
       dI/dBright = (surfaceScale / ΣIcont) × (Iunscaled − I)

6. convolveIGnumpy(): 对 I/V 及全部一阶导数施加显式高斯仪器卷积

7. 返回 diskIntProfAndDerivHalpha 对象，进入 pipeline MEM 迭代
```

---

## 2. 模块概述

`halpha.py` 是 `core/line_models/` 子包中的 H-alpha 专用谱线模型，在 ZDI 管线中作为
`voigt` / `unno` 模型的可替换替代，通过 `config.json` 中 `line_model.model_type = "halpha_compound"` 激活。

```
config.json (model_type="halpha_compound")
        │
        ▼
config_loader.ZDIConfig._parse()
        │ 读取 halpha_* 字段并赋给 par
        ▼
pipeline.ZDIPipeline.run()
        │ lineprofile.lineDataHalpha.from_parameters()
        │ lineprofile.getAllProfDirivHalpha()
        ▼
core/line_models/halpha.py
  ├── lineDataHalpha          — 参数容器
  ├── localProfileAndDerivHalpha — 格点级轮廓与导数核
  ├── diskIntProfAndDerivHalpha  — 盘积分轮廓与全部一阶导数
  └── getAllProfDirivHalpha    — 批量初始化（各相位）
        │
        ▼
core/mem/zdi_adapter.py  (packResponseMatrix / packDataVector)
        │
        ▼
core/fitting.py  (MEM 迭代)
```

**接口边界**：
- 与 Python 其他模块的唯一耦合是 `limbDarkening`（`line_utils.py`）和 `_voigt_faraday_humlicek`（`unno.py`，复用 Humlicek w4）。
- 对外输出的 `diskIntProfAndDerivHalpha` 与 `diskIntProfAndDeriv` / `diskIntProfAndDerivUnno` 接口完全对等，`core/fitting.py` 无需区分模型类型。

---

## 3. 子模块一览

本模型全部代码位于单一文件 `core/line_models/halpha.py`，内部按功能拆分为以下层次：

| 层次 | 主要内容 | 关键符号 |
|------|---------|----------|
| 辅助函数 | 单分量 Voigt 轮廓及解析导数 | `_voigt_component()` |
| 参数容器 | H-alpha 线参数的结构化存储与构造 | `lineDataHalpha` |
| 格点级 | 预计算双峰轮廓、Stokes V 导数核、盘积分 | `localProfileAndDerivHalpha` |
| 观测级 | 完整波长网格、仪器卷积、一阶导数 | `diskIntProfAndDerivHalpha` |
| 批量初始化 | 各相位 `diskIntProfAndDerivHalpha` 列表 | `getAllProfDirivHalpha()` |

`__init__.py` 通过以下语句将以上四个公共符号重导出至包级：
```python
from core.line_models.halpha import (
    lineDataHalpha, localProfileAndDerivHalpha,
    diskIntProfAndDerivHalpha, getAllProfDirivHalpha,
)
```

---

## 4. 公共接口（`core.line_models`）

```python
from core.line_models import (
    # H-alpha 双峰复合 Voigt 模型
    lineDataHalpha,               # 参数容器
    localProfileAndDerivHalpha,   # 格点级轮廓与导数核
    diskIntProfAndDerivHalpha,    # 盘积分轮廓与全部一阶导数
    getAllProfDirivHalpha,         # 批量初始化辅助
)
```

也可直接从子模块导入：

```python
from core.line_models.halpha import lineDataHalpha, getAllProfDirivHalpha
```

---

## 5. 核心类/函数说明

### 5.1 `_voigt_component()` (`halpha.py`)

计算单个 Voigt 分量的轮廓值及对波长的解析导数。被 `localProfileAndDerivHalpha.__init__`
调用两次（发射/自吸收各一次）。

**签名**：
```python
def _voigt_component(
    wl0: float,               # 线心波长（nm）
    width_gauss_kms: float,   # Gaussian 半宽（km/s）
    width_lorentz: float,     # Lorentz/Gauss 比
    amplitude: float,         # 分量强度（A_em 或 A_abs）
    wls: np.ndarray,          # (N_vel, N_cells) Doppler 移位后波长网格
) -> tuple[np.ndarray, np.ndarray]:
    """返回 (comp_I, dcomp_dlambda)，形状均为 (N_vel, N_cells)。"""
```

**导数公式**（与 `profile.py` 符号一致，$u = (\lambda_0 - \lambda)/\Delta\lambda_D$）：

$$
\frac{d(\text{comp\_I})}{d\lambda}
= \frac{2 A}{\Delta\lambda_D}\bigl(u \cdot H(a,u) - a \cdot F(a,u)\bigr)
$$

---

### 5.2 `lineDataHalpha` (`halpha.py`)

H-alpha 双峰复合 Voigt 谱线参数容器。所有数值字段均存储为形状 `(1,)` 的 ndarray，
与 `lineData` / `lineDataUnno` 接口兼容。

#### 构造方法

```python
ldata = lineDataHalpha.from_parameters(
    wavelength_nm, lande_g, limb_darkening, gravity_darkening,
    emission_strength, emission_gauss_kms, emission_lorentz_ratio,
    absorption_strength, absorption_gauss_kms, absorption_lorentz_ratio,
    filling_factor_V=1.0, instRes=-1.0,
)
```

#### 关键字段

| 字段 | 类型 | 含义 |
|------|------|------|
| `wl0` | `ndarray (1,)` | 线心波长（nm），固定 656.28 |
| `g` | `ndarray (1,)` | 有效 Landé 因子，H-alpha 典型值 1.048 |
| `limbDark` | `ndarray (1,)` | 临边昏暗系数 ε（色球层建议 ≤ 0.3） |
| `gravDark` | `ndarray (1,)` | 重力昏暗系数（通常 0） |
| `A_em` | `ndarray (1,)` | 发射成分幅值（> 0） |
| `widthGauss_em` | `ndarray (1,)` | 发射成分 Gaussian 宽度（km/s） |
| `widthLorentz_em` | `ndarray (1,)` | 发射成分 Lorentz/Gauss 比 |
| `A_abs` | `ndarray (1,)` | 自吸收成分深度（≥ 0；构造时强制 ≥ 0） |
| `widthGauss_abs` | `ndarray (1,)` | 自吸收成分 Gaussian 宽度（km/s） |
| `widthLorentz_abs` | `ndarray (1,)` | 自吸收成分 Lorentz/Gauss 比 |
| `fV` | `ndarray (1,)` | Stokes V 填充因子（构造时强制 ≥ 1e-6） |
| `instRes` | `float` | 仪器分辨率 R = λ/Δλ；< 0 = 不卷积 |
| `str`, `widthGauss`, `widthLorentz` | `ndarray (1,)` | 兼容别名（指向发射成分参数） |

---

### 5.3 `localProfileAndDerivHalpha` (`halpha.py`)

格点级双峰轮廓核及 Stokes V 导数核的预计算容器。在 `diskIntProfAndDerivHalpha.__init__`
中构造一次，后由 `updateProfDeriv()` 在每次 MEM 迭代更新亮度/磁场时复用。

#### 构造参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `ldata` | `lineDataHalpha` | 线参数容器 |
| `numPts` | `int` | 波长格点数（保留接口兼容性） |
| `wl_grid` | `ndarray (N_vel, N_cells)` | 各格点 Doppler 移位后波长网格 |

#### 预计算属性

| 属性 | 形状 | 含义 |
|------|------|------|
| `Iunscaled` | `(N_vel, N_cells)` | $1 + A_{\rm em}H_{\rm em} - A_{\rm abs}H_{\rm abs}$，与亮度无关 |
| `VsigKernel` | `(N_vel, N_cells)` | $f_V \cdot g_{\rm Coeff} \cdot (dI_{\rm em}/d\lambda - dI_{\rm abs}/d\lambda)$，与 B 无关 |

#### 主要方法

```python
prof.updateProfDeriv(
    ldata, Blos, dBlos_d, brightMap, numPts, wl_grid, surfaceScale,
    calcDI: int, calcDV: int,
) -> None
```

执行盘积分并填充以下输出属性：

| 属性 | 形状 | 含义 |
|------|------|------|
| `I` | `(N_vel,)` | 归一化盘积分 Stokes I |
| `V` | `(N_vel,)` | 归一化盘积分 Stokes V |
| `dVsum` | `(n_types, nTot, N_vel)` 或 `0` | dV/d(mag coeff)（`calcDV=1`） |
| `dIcsum` | `(N_cells, N_vel)` 或 `0` | dI/d(brightness)（`calcDI=1`） |
| `dVdBright` | `(N_cells, N_vel)` 或 `0` | dV/d(brightness)（两者均为 1） |

---

### 5.4 `diskIntProfAndDerivHalpha` (`halpha.py`)

单观测相位的完整盘积分对象，持有波长网格、所有 Stokes 量及全部一阶导数。

#### 构造参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `visible_grid` | `BatchVisibleGrid` 项 | 该相位可见性几何（面积、视角、投影速度等） |
| `v_mag_cart` | `ndarray (3, N_cells)` | 磁场笛卡尔分量 |
| `d_mag_cart` | `ndarray` 或 `0` | 磁场球谐系数导数 |
| `bright_map` | `BrightMap` | 亮度图对象 |
| `ldata` | `lineDataHalpha` | 线参数容器 |
| `vel_eq` | `float` | 赤道旋转速度（km/s） |
| `wl_grid` | `ndarray (N_vel,)` | 合成波长网格 |
| `calcDI`, `calcDV` | `int` | 亮度/磁场导数计算开关 |

#### 关键输出属性（与 `diskIntProfAndDeriv` / `diskIntProfAndDerivUnno` 对等）

| 属性 | 形状 | 含义 |
|------|------|------|
| `IIc` | `(N_vel,)` | 归一化 Stokes I |
| `VIc` | `(N_vel,)` | 归一化 Stokes V |
| `dVIc` | `(n_types, nTot, N_vel)` | dV/d(mag coeff) |
| `dIIc` | `(N_cells, N_vel)` | dI/d(brightness) |
| `dVdBri` | `(N_cells, N_vel)` | dV/d(brightness) |
| `dImag` | `0` | dI/d(mag coeff)（弱场近似下恒为 0） |

#### 主要方法

```python
spec.convolveIGnumpy(fwhm: float) -> None
```

对 `IIc`、`VIc` 及全部导数数组施加显式高斯仪器卷积（宽度 FWHM = λ/fwhm）。
H-alpha 模型**始终使用显式卷积**（仪器展宽不折入 Gaussian 参数），与 `unno` 策略相同。
`fwhm ≤ 0` 时跳过。

---

### 5.5 `getAllProfDirivHalpha()` (`halpha.py`)

为所有观测相位批量构造 `diskIntProfAndDerivHalpha` 列表，并对每个相位调用
`convolveIGnumpy`（仪器卷积）。

**签名**：
```python
def getAllProfDirivHalpha(
    par,              # ZDIConfig（需有 cycleList, velEq, instrumentRes, calcDV）
    list_grid_view,   # 各相位 BatchVisibleGrid 列表
    vec_mag_cart,     # (3, N_cells) 初始磁场笛卡尔分量
    d_mag_cart0,      # 磁场导数（fitMag=0 时传 0）
    bri_map,          # 亮度图对象
    ldata,            # lineDataHalpha
    wl_syn_set,       # 各相位合成波长网格列表
) -> list[diskIntProfAndDerivHalpha]
```

---

## 6. 配置接口

在 `config.json` 的 `line_model` 节中设置：

```json
{
  "line_model": {
    "model_type": "halpha_compound",
    "wavelength_nm": 656.28,
    "lande_g": 1.048,
    "limb_darkening": 0.0,
    "gravity_darkening": 0.0,
    "emission_strength": 2.5,
    "emission_gauss_kms": 80.0,
    "emission_lorentz_ratio": 0.15,
    "absorption_strength": 1.2,
    "absorption_gauss_kms": 25.0,
    "absorption_lorentz_ratio": 0.10,
    "filling_factor_V": 1.0
  }
}
```

`config_loader.ZDIConfig._parse()` 将上述字段映射为 `par.halpha_*` 属性，
由 `pipeline.ZDIPipeline.run()` 传入 `lineDataHalpha.from_parameters()` 构造线参数对象。

`absorption_strength = 0` 可退化为单峰（纯发射），适用于弱活跃度恒星；
`filling_factor_V < 1` 用于模拟未分辨的混合极性磁场对 Stokes V 的稀释效应。

---

## 7. 线性模型特性与响应矩阵预计算

由于弱场近似下 $\partial I/\partial B_{\rm LOS} = 0$（`dImag = 0`），
且 `VsigKernel` 与磁场无关，**响应矩阵 $R = \partial V/\partial\text{coeff}$
在每次 MEM 迭代内不变化**。`pipeline.py` 通过以下判断将 `calcDV` 在首次迭代后置 0，
避免重复计算：

```python
_model_is_linear = getattr(par, 'line_model_type', 'voigt') in ('voigt', 'halpha_compound')
if (par.calcDV == 1) and (par.calcDI != 1) and _model_is_linear:
    allModeldIdV = memSimple.packResponseMatrix(...)
    par.calcDV = 0
```

这与 `voigt` 模型的优化策略相同，使 `halpha_compound` 在纯磁场反演（`fitMag=1, fitBri=0`）
时具有与 Voigt 模型相当的计算效率。

---

## 8. 与其他谱线模型的对比

| 特性 | `voigt` | `unno` | `halpha_compound` |
|------|---------|--------|--------------------|
| 轮廓形态 | 单峰吸收 | 单峰吸收（磁场依赖） | 双峰（发射 + 自吸收） |
| Stokes V 方法 | 弱场近似 | 完整偏振辐射转移 | 弱场近似（双分量） |
| 响应矩阵线性 | ✓（可预计算） | ✗（每轮迭代更新） | ✓（可预计算） |
| 仪器卷积 | 折入 Gaussian 宽度 | 显式卷积 | 显式卷积 |
| `dImag`（dI/dB） | 0 | 非零 | 0 |
| 适用场景 | LSD 光球线，通常情况 | 强磁场/大填充因子 | H-alpha 活跃恒星 |

---

## 9. 参考文献

- Humlicek, J. (1982). *J. Quant. Spectrosc. Radiat. Transfer*, 27, 437 — Voigt/Faraday-Voigt 函数 w4 算法
- Unno, W. (1956). *Publ. Astron. Soc. Japan*, 8, 108 — Milne-Eddington 偏振辐射转移解（Unno-Rachkovsky 模型基础）
- Landi Degl'Innocenti, E. & Landolfi, M. (2004). *Polarization in Spectral Lines*. Kluwer. §5.3（弱场近似）、§11（双层大气讨论）
- Donati, J.-F. et al. (1997). *MNRAS*, 291, 658 — ZDI 方法综述与 LSD 技术
