# core.magneticGeom — 磁场图 模块说明

> 代码位置：`core/magneticGeom.py`
> 算法基础：Donati et al. (2006)；Petit et al. (2002)；Hobson & Lasenby (1998)

---

## 1. 基本原理

### 1.1 问题背景

ZDI 将恒星表面矢量磁场展开为实数球谐函数（SSH）的有限级数，以有限个展开系数（$\alpha_{\ell m}$、$\beta_{\ell m}$、$\gamma_{\ell m}$）来描述整个磁场拓扑。
这种参数化方式天然地压缩了自由度（从 $N_{\rm cells}$ 个格点分量降为 $\sim 3\ell_{\rm max}^2$ 个系数），并允许对磁场拓扑施加物理约束（如纯极向场、无散场等）。

MEM 反演的目标就是通过最大化磁系数上的熵来寻找最接近"零高阶"的光滑磁场解，同时拟合观测 LSD Stokes V 轮廓。

### 1.2 数学模型

**球谐展开**（遵循 Donati & Brown 1997 惯例）：

$$
B_r = \sum_{\ell=1}^{\ell_{\rm max}}\sum_{m=0}^{\ell} \mathrm{Re}\!\left[\alpha_{\ell m}\,Y_\ell^m(\theta,\phi)\right]
$$
$$
B_\theta = -\sum_{\ell,m} \mathrm{Re}\!\left[\frac{1}{\ell+1}\!\left(\beta_{\ell m}\frac{\partial Y_\ell^m}{\partial\theta} + \gamma_{\ell m}\frac{im}{\sin\theta} Y_\ell^m\right)\right]
$$
$$
B_\phi = -\sum_{\ell,m} \mathrm{Re}\!\left[\frac{1}{\ell+1}\!\left(\beta_{\ell m}\frac{im}{\sin\theta} Y_\ell^m - \gamma_{\ell m}\frac{\partial Y_\ell^m}{\partial\theta}\right)\right]
$$

其中归一化球谐函数：
$$
Y_\ell^m(\theta,\phi) = c_{\ell m}\, P_\ell^m(\cos\theta)\,e^{im\phi}
$$
$$
c_{\ell m} = \sqrt{\frac{2\ell+1}{4\pi}\frac{(\ell-m)!}{(\ell+m)!}}
$$

所有计算（$P_\ell^m$、导数、方位角相位因子）均在 `initMagGeom()` 中**批量预计算**并缓存为矩阵 `YTerm`、`XTerm`、`ZTerm`（形状 $(n_{\rm tot}, N_{\rm cells})$），后续的磁场计算化简为矩阵-向量乘法。

**Legendre 多项式导数递推**（向量化）：
$$
\frac{dP_\ell^m}{dx} = \frac{\ell x\,P_\ell^m - (\ell+m)P_{\ell-1}^m}{x^2 - 1}, \quad x = \cos\theta
$$

### 1.3 算法流程

1. 构造 `magSphHarmonics(nHarmics)`：分配系数数组 $\alpha, \beta, \gamma \in \mathbb{C}^{n_{\rm tot}}$，其中 $n_{\rm tot} = \ell_{\rm max}(\ell_{\rm max}+1)/2 + \ell_{\rm max}$。
2. 调用 `setMagGeomType(type)` 设置几何约束（Full / Poloidal / PotTor / Potential）。
3. 调用 `initMagGeom(clat, lon)`：批量计算 Legendre 多项式及导数（`_compute_legendre_batch`），构建并缓存 `YTerm`、`XTerm`、`ZTerm` 矩阵；结果按网格哈希缓存，重复调用零开销。
4. `getAllMagVectorsCart()` 通过矩阵点积计算球坐标磁场分量，再转换为笛卡尔分量。
5. `getAllMagDerivsCart()` 计算磁场对球谐系数的 Jacobian 矩阵，供 MEM 响应矩阵组装。

---

## 2. 模块概述

`core/magneticGeom.py` 是 ZDIpy 管线中**磁场球谐展开计算**的核心模块，封装了系数存储、Legendre 多项式批量预计算以及球坐标→笛卡尔坐标转换，为 `core/mem/zdi_adapter.py` 提供磁场向量和 Jacobian 矩阵。

```
ZDIConfig (lmax, magGeomType, 系数文件)
       │
       ▼
 magSphHarmonics (磁形状系数)
 ─ alpha[nTot], beta[nTot], gamma[nTot] (complex)
       │
       │ initMagGeom(clat, lon)
       ▼
 预计算缓存
 ─ YTerm (nTot, N_cells)  ← B_r 基函数
 ─ XTerm (nTot, N_cells)  ← B_φ/B_θ 方位角项
 ─ ZTerm (nTot, N_cells)  ← B_θ/B_φ 经向项
       │
       ├─ getAllMagVectorsCart()  → B_xyz (3, N_cells)
       └─ getAllMagDerivsCart()   → dB/dCoeff (n_deriv, nTot, N_cells)
              │
              ▼
       packResponseMatrix (core/mem/zdi_adapter.py)
```

---

## 3. 子模块一览

`core/magneticGeom.py` 为单文件模块，包含以下公共符号：

| 符号 | 类型 | 主要内容 |
|------|------|---------|
| `_compute_legendre_batch` | 函数 | 批量向量化 Legendre 多项式及导数 |
| `magSphHarmonics` | 类 | 系数存储、基函数预计算、磁场/Jacobian 计算 |
| `SetupMagSphHarmonics` | 函数 | 初始化工厂（从文件或零系数） |
| `magSphHarmoicsFromFile` | 函数 | 从 JF Donati 格式实系数文件读入 |
| `magSphHarmoics` | 别名 | `magSphHarmonics` 的向后兼容别名 |
| `dLegendre_dTheta` | 函数 | 单点标量 Legendre 导数（备用） |

---

## 4. 公共接口（`core.magneticGeom`）

```python
from core.magneticGeom import (
    magSphHarmonics,
    SetupMagSphHarmonics,
    magSphHarmoicsFromFile,
)
```

---

## 5. 核心类说明

### 5.1 `magSphHarmonics` (`magneticGeom.py`)

磁场球谐系数存储对象；预计算表面基函数矩阵并支持高效磁场/Jacobian 查询。

#### 构造参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `nHarmics` | `int` | 球谐展开最大阶数 $\ell_{\rm max}$；总系数数 $n_{\rm tot} = \ell_{\rm max}(\ell_{\rm max}+1)/2 + \ell_{\rm max}$ |

#### 主要属性

| 属性 | 形状/类型 | 含义 |
|------|----------|------|
| `nl` | `int` | $\ell_{\rm max}$ |
| `nTot` | `int` | 总 $(\ell,m)$ 对数 |
| `l`, `m` | `(nTot,)` | 各系数对应的 $\ell$, $m$ 指标 |
| `alpha` | `(nTot,) complex` | 径向场系数 |
| `beta` | `(nTot,) complex` | 极向切向场系数（Full/Poloidal 模式独立） |
| `gamma` | `(nTot,) complex` | 环向场系数（Full/PotTor 模式独立） |
| `magGeomType` | `str` | `'full'`/`'poloidal'`/`'pottor'`/`'potential'` |
| `YTerm`, `XTerm`, `ZTerm` | `(nTot, N_cells) complex` | 预计算基函数矩阵（调用 `initMagGeom` 后有效） |

#### 主要方法

```python
# 设置磁几何约束（在 initMagGeom 之前调用）
magGeom.setMagGeomType(magGeomType)  # 'Full'|'Poloidal'|'PotTor'|'Potential'

# 预计算基函数矩阵（输入格点坐标）
magGeom.initMagGeom(clat, lon)       # clat, lon: (N_cells,)

# 球坐标磁场向量 (3, N_cells): [Br, Bclat, Blon]
vMag_sph = magGeom.getAllMagVectors()

# 笛卡尔磁场向量 (3, N_cells): [Bx, By, Bz]
vMag_cart = magGeom.getAllMagVectorsCart()

# Jacobian：球谐系数的磁场导数（球坐标）
#   'full'/'pottor': (5, nTot, N_cells)
#   'poloidal'/'potential': (3, nTot, N_cells)
dB = magGeom.getAllMagDerivs()

# Jacobian（笛卡尔坐标）
dB_cart = magGeom.getAllMagDerivsCart()
```

**磁几何类型约束**

| `magGeomType` | 约束 | 独立系数 |
|---------------|------|---------|
| `full` | 无 | $\alpha, \beta, \gamma$ |
| `poloidal` | $\gamma = 0$ | $\alpha, \beta$ |
| `pottor` | $\beta = \alpha$ | $\alpha, \gamma$ |
| `potential` | $\beta = \alpha,\; \gamma = 0$ | $\alpha$ |

---

### 5.2 `SetupMagSphHarmonics` (`magneticGeom.py`)

从文件或零系数初始化磁场对象的工厂函数。

```python
magGeom = SetupMagSphHarmonics(
    sGrid,              # starGrid 实例（提供 clat, long）
    nHarmics,           # int，ℓ_max
    magGeomType,        # str，磁几何类型
    initMagFromFile,    # int，1 = 从文件读取，0 = 零初始化
    initMagFile,        # str，系数文件路径
    verbose=1,
)
# -> magSphHarmonics 实例（已调用 initMagGeom）
```

---

### 5.3 `_compute_legendre_batch` (`magneticGeom.py`)

批量向量化关联 Legendre 多项式及其导数，通过 `scipy.special.lpmv` 向量化驱动，避免逐格点 Python 循环。

```python
p, pd = _compute_legendre_batch(m_max, n_max, x)
# x: (N_cells,)  cos(colatitude)
# p:  (m_max+1, n_max+1, N_cells)  P_n^m(x)
# pd: (m_max+1, n_max+1, N_cells)  dP_n^m/dx
```

> **注意**：内部乘以 $(-1)^m$ 消除 `scipy.special.lpmv` 的 Condon-Shortley 相位，与参考代码 `lpmn` 保持一致。

---

## 6. 文件结构

```
core/
└── magneticGeom.py     ← 磁场球谐展开模块（本文档描述范围，单文件）
```

---

## 7. 参考文献

- Donati, J.-F. & Brown, S. F. (1997). Zeeman-Doppler imaging of active stars. V. *A&A*, 326, 1135–1142.
- Donati, J.-F. et al. (2006). Spectropolarimetric observations of the Herbig Ae/Be star HD 190073. *MNRAS*, 370, 629–644.
- Petit, P. et al. (2002). Differential rotation of κ Ceti. *MNRAS*, 334, 374–382.
- Hobson, M. P. & Lasenby, A. N. (1998). The entropic prior for distributions with positive and negative values. *MNRAS*, 298, 905–908.
- Humlicek, J. (1982). Optimized computation of the Voigt and complex probability functions. *JQSRT*, 27, 437–444.
