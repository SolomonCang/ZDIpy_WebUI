# 纯色球 ZDI 实现方案

**状态**：设计草案，尚未实现  
**关联文档**：`docs/future/ctts_accretion_model.md`（CTTS 双分量框架，本文不涉及）

---

## 1. 应用场景与定义

"纯色球 ZDI" 指以**色球发射线**（如 He I D3 587.6 nm、Ca II K 393.4 nm
发射核、Hα 发射核）为输入，反演恒星表面磁场分布的流程。
与标准光球吸收线 ZDI 的核心区别在于：

| 特征 | 光球吸收线 | 色球发射线 |
|------|-----------|-----------|
| 归一化轮廓 $I/I_c$ | < 1（吸收） | > 1（发射） |
| 线强符号 $\kappa_L$ | 正 | **负**（方案 B1）或 $\beta<0$（方案 B2）|
| 临边消光 | 临边昏暗（$\varepsilon = 0.3$–$0.8$） | 临边增亮或平坦（$\varepsilon \lesssim 0$ 或 $\beta<0$）|
| Stokes V 方向 | 与磁场 LOS 分量一致 | **反号**（B1 方案中 $\partial I/\partial v$ 反号） |

---

## 2. 方案 B1：弱场近似 + 负线强度（Voigt 模型）

### 2.1 物理原理

弱场近似下，Stokes V 由谱线轮廓的导数给出：

$$V(v) = -\frac{e \lambda_0^2 g_\text{eff}}{4\pi m_e c^2} \cdot B_\text{LOS} \cdot \frac{\partial I}{\partial v}$$

轮廓由 Voigt 函数决定：

$$I(v) = 1 - \kappa_L \cdot H(a, v)$$

**当 $\kappa_L < 0$ 时**（发射线）：

- $I(v) > 1$（轮廓为发射峰）
- $\partial I/\partial v$ 方向与吸收线相反  
- 因此对于给定的 $B_\text{LOS} > 0$，Stokes V 的**正负号翻转**

这与物理上正确的发射线 Zeeman 效应 Stokes V 签名一致：
若用 Queloz et al. (1995) 定义（$V \propto B_\text{LOS} \cdot \partial I/\partial v$），
则发射线的 V 轮廓正好是吸收线的镜像。

### 2.2 需要的代码修改

#### 2.2.1 `core/line_models/profile.py` — 允许负线强

**当前约束**（[profile.py#L77](../../core/line_models/profile.py)）：

```python
# __init__ 中仅有类型断言，无正负约束——代入负值数学上可行
assert lineStr is not None
self.str = np.array([float(lineStr)])   # 允许负值，无需修改
```

`localLineProfile.__init__` 中的核心计算（[profile.py#L235](../../core/line_models/profile.py)）：

```python
self.Iunscaled = 1.0 - lineStr * w4.real
```

当 `lineStr < 0` 时，`1.0 - (负数) * H(v)` 自然给出 $>1$ 的发射峰，
**此行无需修改**。

**`fitLineStrength` 的搜索区间问题**（[line_utils.py#L146](../../core/line_models/line_utils.py)）：

```python
# 当前代码
fitResult = minimize_scalar(
    equivWidComp2,
    bounds=(kL0 * 1e-4, kL0 * 1e4),   # ← kL0 < 0 时区间反转，导致 ValueError
    method='bounded',
    ...
)
```

**需修改**：当 `kL0 < 0` 时调整搜索区间方向：

```python
if kL0 < 0.0:
    lo, hi = kL0 * 1e4, kL0 * 1e-4    # 负值：放大绝对值才能扩大区间
else:
    lo, hi = kL0 * 1e-4, kL0 * 1e4
fitResult = minimize_scalar(
    equivWidComp2,
    bounds=(lo, hi),
    method='bounded',
    ...
)
```

**`equivWidComp2` 的等效宽度符号问题**（[line_utils.py#L79](../../core/line_models/line_utils.py)）：

```python
def equivWidComp2(lineStr, meanEquivWidObs, setSynSpec, lineData):
    ...
    scaleI = 1.0 - (1.0 - spec.IIc) * (lineStr / lineStr0)
    ew = float(np.sum((1.0 - scaleI[:-1]) * dw))   # 吸收线时为正，发射线时为负
```

对发射线，`(1 - I/Ic)` 为负，等效宽度应取绝对值做比较（`meanEquivWidObs` 也需约定为绝对值）。**需修改**为对比绝对值：

```python
return float(np.abs(meanEquivWidObs - np.abs(meanEW)))
```

若不使用 `estimate_strength`（推荐对发射线禁用自动估计），则以上 `fitLineStrength` 问题可绕过。

#### 2.2.2 `config.json` 参数设置

```json
"line_model": {
    "model_type": "voigt",
    "estimate_strength": 0,
    "wavelength_nm": 587.6,
    "line_strength": -0.5,       // 负值 = 发射线（绝对深度由观测校准）
    "gauss_width_kms": 12.0,
    "lorentz_width_fraction": 0.2,
    "lande_g": 1.5,
    "limb_darkening": 0.3
}
```

### 2.3 需修改的文件汇总

| 文件 | 改动 | 估计行数 |
|------|------|---------|
| `core/line_models/line_utils.py` | `fitLineStrength` 搜索区间 + `equivWidComp2` 符号 | ~6 行 |
| `config.json` | `line_strength` 设为负值 | 1 行（用户配置） |

**若禁用 `estimate_strength`（推荐）**，则 `line_utils.py` 无需改动。

### 2.4 验证方法

1. 用合成数据（已知 $B$ 分布）生成发射线轮廓，反演磁场，与真实值比较
2. 与 C 版 `zdipot` 对同一发射线输入的结果对比（若有）
3. 检查 Stokes V 签名方向：发射线磁场 LOS 分量为正时，V 蓝翼应为**负**（与吸收线反号）

---

## 3. 方案 B2：Milne-Eddington 发射（UR 模型 $\beta < 0$）

### 3.1 物理原理

Milne-Eddington 大气中，连续谱归一化的 Stokes I 由下式给出（Landi Degl'Innocenti & Landolfi 2004, Eq. 9.110）：

$$\frac{I}{I_c} = \frac{1 + \beta\mu \cdot \eta_I (\eta_I^2 + \rho_Q^2 + \rho_V^2) / \Delta}{1 + \beta\mu}$$

其中 $\beta = B_1 / B_0$ 为 Planck 函数坡度（$B_\nu(\tau) = B_0 + B_1 \tau$）。

- **$\beta > 0$**（标准光球）：温度随光深增加，$I/I_c < 1$（吸收线）
- **$\beta < 0$**（温度逆转 / 色球）：上层比连续谱形成层更热，$I/I_c > 1$（发射线）

$\beta < 0$ 时，整个 UR 计算（包括 Stokes V）在数学上完全自洽，
$\Delta$ 分母不变号，无需额外处理。

**与方案 B1 的根本区别**：B2 通过修改源函数结构来实现发射，轮廓形状由完整 ME 大气决定（不对称性、临边效应均自然包含）；B1 则是简单翻转轮廓深度。

### 3.2 临边增亮与 $\beta$、$\varepsilon$ 的关系

ME 大气的临边昏暗系数 $\varepsilon$ 与 $\beta$ 的关系：

$$\varepsilon = \frac{\beta}{1 + \beta}  \iff  \beta = \frac{\varepsilon}{1 - \varepsilon}$$

- $\varepsilon \in (0, 1)$：临边昏暗，$\beta > 0$（正常光球）
- $\varepsilon < 0$：临边增亮，$\beta < 0$（色球温度倒转）
- $\varepsilon > 1$：强临边增亮，$\beta > 1$（理论上也可能）

色球发射线一般对应 $\varepsilon \approx -0.5$ 至 $-2$（$\beta \approx -0.33$ 至 $-2$）。

### 3.3 需要的代码修改

#### 3.3.1 `core/line_models/unno.py` — 允许显式负 $\beta$

**当前 `from_parameters` 中的自动推导逻辑**（[unno.py#L289](../../core/line_models/unno.py)）：

```python
beta_val = float(beta)
if beta_val <= 0.0:                          # ← 问题：β<0 被当作"未指定"
    limbd = float(limb_darkening)
    beta_val = limbd / (1.0 - limbd)         # 用 ε 推导，忽略了显式负值
```

**需修改**：使用哨兵值 `None`（或特殊浮点数）区分"未指定"和"显式负值"：

```python
@classmethod
def from_parameters(
    cls,
    ...,
    beta: float | None = None,    # None = 从 limb_darkening 自动推导
    ...
) -> "lineDataUnno":
    ...
    if beta is None:
        limbd = float(limb_darkening)
        if abs(limbd - 1.0) < 1e-10:
            limbd = 1.0 - 1e-6
        beta_val = limbd / (1.0 - limbd)
        print(f'lineDataUnno: beta derived from limb_darkening: {beta_val:.4f}')
    else:
        beta_val = float(beta)     # 允许负值，直接使用
    obj.beta = np.array([beta_val])
```

同时需更新 `config_loader.py` 和 `pipeline.py`：`unno_beta` 字段应支持 `null`（JSON null → Python `None`）以触发自动推导路径。

#### 3.3.2 `from_file` 中的 $\beta$ 解析

**当前代码**（[unno.py#L340](../../core/line_models/unno.py)）：

```python
beta_raw = float(parts[2])
if beta_raw <= 0.0:
    limbd = float(parts[6])
    beta_raw = limbd / (1.0 - limbd)    # ← 同样将负值误认为"未设定"
```

**需修改**：用文件中的哨兵值（如 `NaN` 或 `-999`）表示"从 ε 推导"，
其他负值直接使用：

```python
beta_raw = float(parts[2])
if np.isnan(beta_raw):               # NaN = 从 limbDark 推导
    limbd = float(parts[6])
    if abs(limbd - 1.0) < 1e-10:
        limbd = 1.0 - 1e-6
    beta_raw = limbd / (1.0 - limbd)
# 否则（包括负值）直接使用 beta_raw
obj.beta = np.append(obj.beta, beta_raw)
```

#### 3.3.3 `config.json` 参数设置

```json
"line_model": {
    "model_type": "unno",
    "estimate_strength": 0,
    "wavelength_nm": 587.6,
    "line_strength": 0.5,
    "gauss_width_kms": 12.0,
    "lorentz_width_fraction": 0.2,
    "lande_g": 1.5,
    "limb_darkening": -0.5,        // 临边增亮（color = chrom.）
    "unno_beta": -0.33,            // 显式负 β（或设 null 由 limb_darkening 推导）
    "unno_filling_factor_I": 1.0,
    "unno_filling_factor_V": 1.0
}
```

若 `limb_darkening < 0`，则自动推导给出 `β = ε/(1-ε) < 0`，
可直接沿用当前推导式（无需修改），**仅需去掉 `β ≤ 0 则推导` 的条件**：
对 `limb_darkening = -0.5` 推导出 `β = -0.5/(1-(-0.5)) = -1/3`，符合预期。

#### 3.3.4 `_unno_profile` 的数值稳定性检查

当 `β < 0` 且 `μ < |1/β|` 时，`1 + β·μ` 可能接近零（临边格点）。
现有代码已有防护（[unno.py#L213](../../core/line_models/unno.py)）：

```python
safe_denom = np.where(np.abs(one_plus_bmu) > 1e-14, one_plus_bmu, 1e-14)
```

但对 `β ≈ -1` 时大量格点的 `μ ≈ 1` 区域，`1 + β·μ ≈ 0` 会导致轮廓发散。
建议增加更严格的阈值检查，或在 `β < 0` 路径下对 `|1+βμ| < 0.01` 的格点输出警告。

### 3.4 需修改的文件汇总

| 文件 | 改动 | 估计行数 |
|------|------|---------|
| `core/line_models/unno.py` | `from_parameters` β 哨兵逻辑 + `from_file` β 解析 | ~15 行 |
| `config_loader.py` | `unno_beta` 支持 `None` | ~3 行 |
| `pipeline/pipeline.py` | `from_parameters` 传 `beta=None` 而非 `-1.0` | ~2 行 |

### 3.5 验证方法

1. 对 `β = -0.33`（$\varepsilon = -0.5$）合成轮廓，确认 $I/I_c > 1$
2. 临边格点稳定性：确认 `μ → 0` 时轮廓收敛而非发散
3. Stokes V 方向：发射轮廓的 $V$ 签名应与吸收线反号（因 $\partial I/\partial \lambda$ 方向不变，但 $\beta \cdot \mu/(1+\beta\mu)$ 变号）
4. 与方案 B1 的交叉验证：对同一简单情况（均匀 $B$，小 $\kappa_L$），两方案的 $V$ 幅度应相近

---

## 4. 方案比较

| 维度 | B1：Voigt 负线强 | B2：UR $\beta<0$ |
|------|----------------|----------------|
| 轮廓模型 | 弱场 Voigt | 完整 ME / UR |
| 数值精度 | 弱场近似（适用 $B \lesssim 1$ kG）| 任意场强 |
| 轮廓形状 | 对称 Voigt | ME 大气的反对称效应、临边效应 |
| Stokes V 方向 | 自动（$\partial I/\partial v$ 反号）| 需检查（见 § 3.5） |
| 改动量 | 极小（0–6 行）| 中等（~20 行）|
| 适用谱线 | 弱场、LSD 叠加线 | 单条强线、有明确 ME 物理意义 |
| `estimate_strength` | 需禁用或修改符号处理 | 需禁用 |

### 使用建议

- **快速验证 / LSD 叠加发射线**：优先用 B1。`line_strength` 直接设负值，
  禁用 `estimate_strength`，数分钟内可运行完整反演。

- **精确建模单条发射线（He I D3、Ca II K 核）**：选 B2，
  结合已知 ME 参数（$\varepsilon$、$\beta$）。
  改动集中在 `unno.py` 的参数解析逻辑，前向模型本身无需改动。

- **联合光球 + 色球反演**：参见 `docs/future/ctts_accretion_model.md` 的 B3/CTTS 框架，
  将两层 $\beta$ 参数化后联合使用。

---

## 5. 变更文件一览（两方案合并）

| 文件 | 方案 | 改动类型 | 估计行数 |
|------|------|---------|---------|
| `core/line_models/line_utils.py` | B1 | 搜索区间 + EW 符号 | ~6 行 |
| `core/line_models/unno.py` | B2 | β 哨兵值逻辑 | ~15 行 |
| `config_loader.py` | B2 | `unno_beta` 允许 `None` | ~3 行 |
| `pipeline/pipeline.py` | B2 | 参数透传 | ~2 行 |
| `config.json` | B1/B2 | 用户配置调整 | —（仅配置）|

---

## 6. 参考文献

| 主题 | 参考 |
|------|------|
| Milne-Eddington UR 方程 | Landi Degl'Innocenti & Landolfi 2004, §9.4 |
| 弱场近似 Stokes V | Queloz et al. 1995, A&A 298, 164 |
| He I D3 色球 ZDI | Donati et al. 2010, MNRAS 409, 1347 |
| Ca II K 发射核应用 | Hussain et al. 2000, MNRAS 318, 961 |
| 临边增亮与色球 ME 参数 | Auer & Heasley 1978, A&A 64, 67 |
