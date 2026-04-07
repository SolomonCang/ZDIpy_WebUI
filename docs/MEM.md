# MEM (Maximum Entropy Method) 模块说明

> 代码位置：`core/mem/`
> 算法基础：Skilling & Bryan (1984, MNRAS 211, 111)

---

## 1. 基本原理

### 1.1 问题背景

Zeeman Doppler Imaging（ZDI）是一种间接成像反演技术：从恒星在不同自转相位下采集的 Stokes I/V LSD 线轮廓中，重建二维恒星表面的**亮度图**和**矢量磁场图**。
由于观测数量远小于未知图像像素数，该反演问题严重欠定（ill-posed），存在无穷多解能拟合观测数据，因此需要正则化约束来选取"最平滑"（最优不偏）的解。

最大熵法（Maximum Entropy Method，MEM）通过**最大化图像熵**来实现正则化：
在所有满足数据拟合约束的解中，选取熵值最高（信息量最低、最接近"均匀背景"）的那一个。

### 1.2 数学模型

设观测数据向量 $\mathbf{d} \in \mathbb{R}^{N_d}$，正向模型 $\mathbf{F}(\mathbf{x})$ 给出合成谱，$\sigma_i$ 为各点误差，图像向量 $\mathbf{x}$ 包含亮度和磁场系数。

**数据约束**（归一化 $\chi^2$）：
$$
\mathcal{C}(\mathbf{x}) = \frac{1}{N_d}\sum_{i=1}^{N_d} \frac{(d_i - F_i(\mathbf{x}))^2}{\sigma_i^2}
$$
目标为使 $\mathcal{C} = \mathcal{C}_{\rm aim}$（通常 $\mathcal{C}_{\rm aim} = 1$，即每自由度拟合误差为 $1\sigma$）。

**熵泛函**（以标准正熵为例）：
$$
S(\mathbf{x}) = -\sum_i w_i\bigl[x_i\bigl(\ln\tfrac{x_i}{m_i}-1\bigr)+m_i\bigr]
$$
其中 $m_i$ 为默认（背景）图像，$w_i$ 为各像素权重。

**优化目标**：
$$
\max_{\mathbf{x}}\; S(\mathbf{x}) \quad \text{subject to} \quad \mathcal{C}(\mathbf{x}) = \mathcal{C}_{\rm aim}
$$
等价于最大化 Lagrangian $\mathcal{L} = S - \lambda \mathcal{C}$，最优点满足 $\nabla S = \lambda \nabla \mathcal{C}$。

### 1.3 算法流程（Skilling-Bryan 1984）

Skilling-Bryan 算法在每次迭代中构建一个低维次空间（最多 `max_search_dirs` 个方向），在该次空间内精确求解约束优化问题：

1. 计算当前熵梯度 $\nabla S$ 和约束梯度 $\nabla \mathcal{C}$（及二阶导数 $\nabla^2 S$）。
2. 以 $\nabla S$、$\nabla \mathcal{C}$、$\nabla^2 S \cdot \nabla S$、$\nabla^2 S \cdot \nabla \mathcal{C}$（及更高阶乘积）构造搜索方向基，Gram-Schmidt 正交化。
3. 在该基张成的子空间内，通过二分法求解 $\mathcal{C}(\mathbf{x} + \delta\mathbf{x}) = \mathcal{C}_{\rm aim}$，同时约束步长 $\|\delta\mathbf{x}\|^2 \leq l_0^2 \cdot I_{\rm tot}$（步长稳定性限制）。
4. 更新图像 $\mathbf{x} \leftarrow \mathbf{x} + \delta\mathbf{x}$，施加边界约束（充填因子 $f \in (0, \hat{f})$、亮度 $x > 0$ 等）。
5. 计算收敛检验统计量：
$$
T = \frac{|\nabla S - \lambda^*\nabla \mathcal{C}|}{|\nabla S| + |\lambda^*\nabla \mathcal{C}|} \to 0 \quad (\text{收敛})
$$
当 $T < T_{\rm aim}$（配置项 `inversion.test_aim`）时停止迭代。

---

## 2. 概述

`core/mem/` 是 ZDIpy 管线中负责**最大熵法反演（MEM inversion）**的独立子包。
它将纯数学优化算法与 ZDI 物理层解耦，通过清晰的回调接口和打包工具链接两者。

```
输入观测数据 (Data, sig2)
       │
       ▼
  ZDI 物理模型
  (BrightnessGeom / MagneticGeom → SyntheticSpectra)
       │  packDataVector / packImageVector / packResponseMatrix
       ▼
 ┌─────────────────────────────┐
 │  core/mem/  MEM 子包        │
 │  ┌─────────────────────┐   │
 │  │  generic.py         │   │  ← 纯算法核心（Skilling-Bryan）
 │  │  MEMOptimizer       │   │
 │  └──────────┬──────────┘   │
 │             │ 回调接口       │
 │  ┌──────────▼──────────┐   │
 │  │  zdi_adapter.py     │   │  ← ZDI 专用打包/解包 + mem_iter
 │  │  mem_iter           │   │
 │  └──────────┬──────────┘   │
 │             │               │
 │  ┌──────────▼──────────┐   │
 │  │  iteration_manager  │   │  ← 迭代循环控制器
 │  │  monitoring         │   │  ← 历史记录 / ETA
 │  │  optimization       │   │  ← 响应矩阵缓存 / 稳定性检查
 │  └─────────────────────┘   │
 └─────────────────────────────┘
       │  unpackImageVector
       ▼
  更新后的物理参数 (briMap, magGeom)
```

---

## 3. 子模块一览

| 文件 | 主要内容 | 关键类/函数 |
|------|---------|------------|
| `generic.py` | 纯 MEM 算法（无 ZDI 依赖），Skilling-Bryan 全流程 | `MEMOptimizer` |
| `monitoring.py` | 迭代历史记录 + 实时进度/ETA | `IterationHistory`, `ProgressMonitor` |
| `optimization.py` | 响应矩阵 LRU 缓存，数值稳定性检查 | `ResponseMatrixCache`, `StabilityMonitor`, `CacheStats` |
| `iteration_manager.py` | 迭代循环高级控制器，收敛判断 | `IterationManager`, `ConvergenceChecker`, `create_iteration_manager_from_config` |
| `zdi_adapter.py` | ZDI 专用：数据/图像/响应矩阵打包，`mem_iter` 入口 | `mem_iter`, `packDataVector`, `packImageVector`, `packResponseMatrix`, `unpackImageVector`, `constantsMEM`, `setEntropyWeights` |
| `saim_adapter.py` | SAIM 熵目标控制程序封装 | `getSaimQuad`, `control` |
| `__init__.py` | 包级公共 API，统一重导出 | 见下节 |

---

## 4. 公共接口（`core.mem`）

```python
from core.mem import (
    # 算法核心
    MEMOptimizer,

    # 监控
    IterationHistory,
    ProgressMonitor,

    # 优化工具
    ResponseMatrixCache,
    StabilityMonitor,
    CacheStats,

    # 迭代控制
    IterationManager,
    ConvergenceChecker,
    create_iteration_manager_from_config,

    # ZDI 管线函数
    mem_iter,
    get_s_grads,
    packDataVector,
    packModelVector,
    packImageVector,
    unpackImageVector,
    packResponseMatrix,
    constantsMEM,
    setEntropyWeights,
)
```

---

## 5. 核心类说明

### 5.1 `MEMOptimizer` (`generic.py`)

通用 MEM 优化器，通过**回调函数**接入任意熵形式和数据拟合方式。

#### 构造参数

| 参数 | 类型 | 说明 |
|------|------|------|
| `compute_entropy_callback` | `Callable` | `(Image, weights, n1, n2, ntot) → (S0, gradS, gradgradS, Itot)` |
| `compute_constraint_callback` | `Callable` | `(Data, Fmodel, sig2, Resp) → (C0, gradC)` |
| `boundary_constraint_callback` | `Callable` *(可选)* | `(Image, n1, n2, ntot) → Image_corrected` |
| `max_search_dirs` | `int` | 最大搜索方向数（推荐 3–10，默认 10） |
| `step_length_factor` | `float` | 步长因子 `l₀² = L_fac × Itot`（默认 0.3） |
| `convergence_tol` | `float` | 二分收敛容差（默认 1e-5） |

#### 方法

```python
optimizer.iterate(
    Image, Fmodel, Data, sig2, Resp,
    weights, n1, n2, ntot,
    fixEntropy=0,    # 0=拟合χ²，1=拟合熵目标
    targetAim=None,
) -> (entropy, chi2, test, Image_new)
```

**返回值：**

| 名称 | 含义 |
|------|------|
| `entropy` | 当前总熵值 S₀ |
| `chi2` | 当前 χ²（C₀） |
| `test` | Skilling-Bryan 收敛检验统计量（接近 0 表示最优） |
| `Image_new` | 更新后的图像向量 |

---

### 5.2 ZDI 适配层 (`zdi_adapter.py`)

#### 图像向量分段约定

ZDIpy 使用三段式图像向量（长度 `ntot`）：

```
Image[0 : n1]     — 亮度图（标准熵，无上界）
Image[n1 : n2]    — 充填因子（有界熵，0 < f < ffIMax）
Image[n2 : ntot]  — 磁场球谐系数（正负熵，可正可负）
```

由 `constantsMEM` 自动计算 `n1`, `n2`, `ntot`。

#### 关键函数签名

```python
# 组装数据向量
Data, sig2 = packDataVector(obsSet, fitBri, fitMag)

# 组装图像向量
Image = packImageVector(briMap, magGeom, magGeomType, fitBri, fitMag)

# 组装响应矩阵 dModel/dImage
Resp = packResponseMatrix(setSynSpec, nDataTot, npBriMap, magGeom,
                          magGeomType, fitBri, fitMag)

# 一步 MEM 迭代
entropy, chi2, test, Image, entStand, entFF, entMag = mem_iter(
    n1, n2, ntot, Image, Data, Fmodel, sig2, Resp,
    weights, defImg, defIpm, ffIMax, targetAim, fixEntropy)

# 写回物理参数
unpackImageVector(Image, briMap, magGeom, magGeomType, fitBri, fitMag)
```

---

### 5.3 `IterationManager` (`iteration_manager.py`)

高级迭代循环控制器，整合收敛检查、进度监控和历史记录。

```python
manager = create_iteration_manager_from_config(config, verbose=1)

for iIter in range(par.numIterations):
    manager.start_iteration()
    # … 计算 entropy, chi2 …
    manager.record_iteration(chi2=chi2, entropy=entropy)
    should_stop, reason = manager.should_stop(chi2)
    if should_stop:
        break

summary = manager.get_summary()  # -> Dict
history = manager.get_iteration_history()  # -> IterationHistory | None
```

**停止条件（按优先级）：**
1. 达到 `max_iterations`
2. 连续 `stall_threshold` 次迭代的相对 Δχ²/χ² < `convergence_threshold`

---

### 5.4 `ResponseMatrixCache` (`optimization.py`)

对 ZDI 响应矩阵做 LRU 缓存，避免在相同磁几何条件下重复计算：

```python
cache = ResponseMatrixCache(max_size=5, verbose=1)
Resp = cache.get_or_compute(magGeom, obsSet, compute_fn)
stats = cache.get_stats()  # -> CacheStats
```

缓存键基于磁场系数（`alpha`, `beta`, `gamma` 数组的 MD5 哈希）和 `obsSet` 对象 id。

---

## 6. 熵形式

ZDIpy 对图像向量不同分段使用三种熵形式：

| 分段 | 名称 | 表达式 | 参考 |
|------|------|--------|------|
| `[0:n1]` | 标准正熵 | $S = -\sum w_i \bigl[I_i(\ln(I_i/m_i)-1)+m_i\bigr]$ | Skilling & Bryan (1984) Eq. 6 |
| `[n1:n2]` | 充填因子熵 | $S = -\sum w_j \bigl[f_j\ln(f_j/m)+(\hat f-f_j)\ln\tfrac{\hat f-f_j}{\hat f-m}\bigr]$ | Unruh & Collier Cameron (1995) Eq. 4 |
| `[n2:ntot]` | 磁场正负熵 | $S = \sum w_k\bigl[\psi_k - 2q - a_k\ln\tfrac{\psi_k+a_k}{2q}\bigr]$，$\psi_k=\sqrt{a_k^2+4q^2}$ | Hobson & Lasenby (1998) Eq. 8 |

其中 $m$ = `defImg`（默认图像强度），$q$ = `defIpm`（磁系数默认值），$\hat f$ = `ffIMax`（充填因子上限）。

---


## 7. 文件结构

```
core/
└── mem/                        ← MEM 子包（本文档描述范围）
    ├── __init__.py             ← 公共 API 统一入口
    ├── generic.py              ← 纯 MEM 算法核心
    ├── monitoring.py           ← 迭代历史 / 进度监控
    ├── optimization.py         ← 缓存 / 数值稳定性工具
    ├── iteration_manager.py    ← 迭代循环控制器
    ├── zdi_adapter.py          ← ZDI 打包/解包 + mem_iter
    └── saim_adapter.py         ← SAIM 熵目标包装
```

---

## 8. 参考文献

- Skilling, J. & Bryan, R. K. (1984). Maximum entropy image reconstruction: general algorithm. *MNRAS*, 211, 111–124.
- Unruh, Y. C. & Collier Cameron, A. (1995). The sensitivity of Zeeman Doppler imaging to field orientation. *MNRAS*, 273, 1–16.
- Hobson, M. P. & Lasenby, A. N. (1998). The entropic prior for distributions with positive and negative values. *MNRAS*, 298, 905–908.
- Folsom, C. P. et al. (2018). The large-scale magnetic field and differential rotation of ξ Bootis A. *MNRAS*, 481, 5286–5295.
