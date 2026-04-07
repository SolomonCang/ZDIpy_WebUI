# CTTS 吸积模型移植设计方案

**参考实现**：`/Users/tianqi/Documents/Codes_collection/ZDI_and/CTTSzdi2/unno_ctts.c`  
**目标模块**：`core/line_models/unno.py`  
**状态**：设计草案，尚未实现

---

## 1. 物理背景与动机

### 1.1 经典 T Tauri 星（CTTS）的特殊性

经典 T Tauri 星（CTTS）处于磁控吸积阶段：来自吸积盘的物质沿磁力线
下落至星面，在吸积热斑（accretion hot spot）处激发一个额外的**色球/吸积发射分量**，
叠加在光球吸收轮廓之上。与弱 T Tauri 星（WTTS）和普通 MS 星不同，
CTTS 的谱线轮廓由至少两个独立大气层贡献：

| 分量 | 来源 | 代表亮度 | 轮廓性质 |
|------|------|----------|----------|
| 光球（磁区）| 恒星光球 + 表面磁场 | `Cq` | 磁 Milne-Eddington 吸收 |
| 吸积/色球 | 热斑、色球发射 | `Cm` | 仅 π 分量（无有效 Zeeman 分裂）|

`Cm` 在空间上逐格点变化——吸积热斑处 `Cm` 极大，无吸积区域 `Cm ≈ 0`。
正是通过联合反演 `Cm` 图和磁场图，才能把两种贡献解耦。

### 1.2 WTTS 模式对照

弱 T Tauri 星（WTTS）同样可能存在色球贡献，但其物理图像更接近 MS 星，
轮廓振幅仅依赖光球亮度对比度：

```
# WTTS 模式（C 版被注释的代码）
pmul = imul = 1 + 2*(Cm - defC)/defC   # defC = 0.999
```

WTTS 模式的色球贡献比 CTTS 弱得多，此文档仅关注 CTTS 实现；
WTTS 分支可在后续作为 `fac > 0` 时的扩展。

---

## 2. 数学与物理原理

### 2.1 CTTS 双分量轮廓公式

设 `fac < 0` 为 CTTS 模式开关（`fac` 称为"轮廓乘子符号位"，
其绝对值通常为 1，负号仅作标志）。

**连续谱归一化局部轮廓**：

$$I(\lambda) = I_\text{mag}(\lambda) \cdot f_\text{fff} + I_\pi(\lambda) \cdot (1 - f_\text{fff})$$

其中：
- $f_\text{fff} = $ `fI`：磁活动区面积填充因子（通常 0.8–0.9）
- $I_\text{mag}$：完整 Unno-Rachkovsky 磁轮廓（三 Zeeman 分量）
- $I_\pi$：仅 $\pi$ 分量的非磁 Milne-Eddington 轮廓（色球/吸积模拟）

**混合后的重缩放**（CTTS 核心操作）：

$$I_\text{CTTS}(\lambda) = p_\text{mul} \cdot |fac| \cdot \bigl[I(\lambda) - I_c\bigr] + I_c$$

$$V_\text{CTTS}(\lambda) = fac \cdot p_\text{mul} \cdot V_\text{UR}(\lambda)$$

其中 $I_c = f_\text{fff} \cdot cont + (1 - f_\text{fff}) \cdot ccc$（双分量连续谱基准），
$p_\text{mul}$ 为振幅乘子：

$$p_\text{mul} = \frac{10 \cdot (C_\text{def} - C_m)}{C_\text{def}} + 1, \qquad C_\text{def} = 0.999$$

$$C_m = \begin{cases} \min(C_m^i,\; C_\text{def}) & C_m^i \geq 0 \\ C_q^i & C_m^i < 0 \end{cases}$$

**物理解读**：当 `Cm` 接近 `defC`（最大亮度，强吸积）时，`pmul → 1`，
轮廓受到强烈重缩放，等效于发射填充使吸收深度变浅；
当 `Cm → 0`（无吸积）时，`pmul → 11`，轮廓深度正常。

### 2.2 π-only 非磁 Milne-Eddington 轮廓

第二分量 $I_\pi$ 仅使用 $\pi$ 线（Zeeman 分裂为零的线心分量）：

$$I_\pi(\lambda) = \frac{1 + \beta / (1 + \eta_0 H_\pi) }{1 + \beta}$$

其中 $H_\pi = H(a, \nu)$ 为不做 Zeeman 移位的 Voigt 函数，
$\eta_0$ = `ldata.str`，$\beta$ = `ldata.beta`。
这模拟了一个无磁场、但有相同热力学结构（ME 大气）的色球层。

### 2.3 Stokes I 和 V 对 Cm 的导数

由于 `pmul` 依赖 $C_m$，对 $C_m$ 的导数为：

$$\frac{\partial p_\text{mul}}{\partial C_m} = -\frac{10}{C_\text{def}}$$

因此：

$$\frac{\partial I_\text{CTTS}}{\partial C_m^i} = |fac| \cdot \frac{\partial p_\text{mul}}{\partial C_m} \cdot (I(\lambda) - I_c) \cdot \delta_{ii'}$$

$$\frac{\partial V_\text{CTTS}}{\partial C_m^i} = fac \cdot \frac{\partial p_\text{mul}}{\partial C_m} \cdot V_\text{UR}(\lambda) \cdot \delta_{ii'}$$

这些导数（逐格点，形状 `(N_cells, N_vel)`）需打包进 MEM 响应矩阵，
与现有 $\partial I / \partial C_q^i$（即 `dIcsum`）并列为图像向量中的两列块。

---

## 3. 代码架构方案

### 3.1 现有路由机制（不改动）

Pipeline 通过 `lineData` 对象类型选择盘积分类：

```python
# pipeline.py — 构造时选模型
if par.line_model_type == 'unno':
    lineData = lineDataUnno.from_parameters(...)   # 触发 UR 路径
else:
    lineData = lineData.from_parameters(...)        # 触发 Voigt 路径

# line_utils.py — 运行时路由（已存在）
if isinstance(lineData, lineDataUnno):
    DiskIntCls = diskIntProfAndDerivUnno
else:
    DiskIntCls = diskIntProfAndDeriv
```

添加 CTTS 后，**此路由代码无需任何修改**。
CTTS 模式通过 `lineDataUnno.fac < 0` 在 `_unno_profile` 内部分支，
对 pipeline 完全透明。

### 3.2 图像向量扩展

C 版图像向量结构（5 个分量块）：

```
image = [ Cq_0..N,  Cm_0..N,  Br_0..N,  Bt_0..N,  Bp_0..N ]
          imsw[0]   imsw[1]   imsw[2]   imsw[3]   imsw[4]
```

Python 版现有图像向量只含 `Cq`（亮度）+ 球谐系数。
引入 CTTS 后需增加 `Cm` 块，MEM 响应矩阵 $R$ 增加 $N_\text{cells}$ 列：

```
R_extended = [ R_Cq | R_Cm | R_Bsph ]
               (m×N)  (m×N)  (m×K)
```

其中 $m$ = 总数据点数，$N$ = 格点数，$K$ = 球谐系数数。

### 3.3 熵函数选择

`Cm` 在物理上有界于 $[0, C_\text{def}]$，
应使用 MEM 的**有界熵**（对应 C 版 `mem_nl.c` 中 `n1 ≤ j < n2` 段）：

$$S_{C_m} = -\sum_i w_i \left[ C_m^i \ln\frac{C_m^i}{C_C} + (A - C_m^i)\ln\frac{A - C_m^i}{A - C_C} \right]$$

其中 $C_C = C_\text{def} \cdot A$。
这与现有 `mem.py` 的 `n1`/`n2` 分区对应，增加一个区段即可。

---

## 4. 预期实现目标

1. **前向模型正确性**：对给定 `(Cq, Cm, B, γ, χ)` 的合成轮廓与 C 版 `zdipot` 输出一致（数值偏差 < 0.1%）。

2. **`Cm` 参与反演**：`fitcm > 0` 时 `Cm` 图经 MEM 迭代优化，与磁场图联合反演。

3. **`fitcm = 0` 向后兼容**：当 `fac >= 0` 或 `fitcm = 0` 时，代码退化为现有标准 UR 路径，不影响非 CTTS 数据集。

4. **无 pipeline 改动**：`pipeline.py`、`fitting.py`、`mem/` 等模块的对外接口保持不变。

5. **填充因子 `fI` 可配置**：现有 C 版 `fff` 在代码中按星名硬编码；Python 版通过 `lineDataUnno.fI`（config.json 字段）动态配置。

---

## 5. 具体改动方案

### 5.1 `lineDataUnno`（`core/line_models/unno.py`）

**新增字段**：

| 字段 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `fac` | `float` | `1.0` | `< 0` 启用 CTTS 模式；`|fac|` 通常为 1 |
| `Cm_default` | `float` | `0.0` | `fitcm=0` 时全星面使用的固定 `Cm` 值 |

**`from_parameters` 新增参数**：

```python
@classmethod
def from_parameters(
    cls,
    ...,                            # 现有参数不变
    fac: float = 1.0,               # 新增
    Cm_default: float = 0.0,        # 新增
) -> "lineDataUnno":
```

**`from_file` 新增列解析**：文件格式第 11 列读取 `fac`，第 12 列读取 `Cm_default`。

---

### 5.2 `_unno_profile`（内部函数）

**新增参数**：

```python
def _unno_profile(
    width_gauss: float,
    width_lorentz: float,
    ldata: "lineDataUnno",
    wls: np.ndarray,          # (N_vel, N_cells)
    view_angle: np.ndarray,   # (N_cells,)
    Bmod: np.ndarray,         # (N_cells,)
    Btheta: np.ndarray,       # (N_cells,)
    Cm: np.ndarray | None = None,   # (N_cells,) — 新增，CTTS 吸积亮度
) -> tuple[np.ndarray, np.ndarray]:
```

**新增 CTTS 分支**（在返回前插入）：

```python
_DEFC = 0.999
fac_val = float(ldata.fac)

if fac_val < 0.0 and Cm is not None:
    # 1. 计算 π-only 非磁轮廓（第二分量）
    etaI_pi = 1.0 + eta_half * eta_P          # (N_vel, N_cells)，仅 π 分量
    beta_mu = float(ldata.beta[0]) * mu       # (N_cells,)
    safe_denom = np.where(np.abs(1.0 + beta_mu) > 1e-14, 1.0 + beta_mu, 1e-14)
    sI_pi = (1.0 + beta_mu / etaI_pi) / safe_denom   # (N_vel, N_cells)

    # 2. 双分量混合
    fI_val = float(ldata.fI[0])
    I_mix = sI * fI_val + sI_pi * (1.0 - fI_val)     # (N_vel, N_cells)

    # 3. 连续谱基准（两分量等效于 1.0，因为 sI 和 sI_pi 均已归一化）
    I_cont = fI_val * 1.0 + (1.0 - fI_val) * 1.0     # = 1.0，保留以便扩展

    # 4. 按亮度对比计算逐格点 pmul
    br = np.where(Cm >= 0.0, np.minimum(Cm, _DEFC), np.clip(-Cm, 0.0, _DEFC))
    # 注：Cm < 0 退回 Cq 的逻辑在调用侧处理，此处只确保在界内
    pmul = 10.0 * (_DEFC - br) / _DEFC + 1.0          # (N_cells,)

    # 5. 重缩放（|fac| 通常为 1，fac 为负只作标志）
    abs_fac = abs(fac_val)
    sI_out = pmul * abs_fac * (I_mix - I_cont) + I_cont
    sV_out = fac_val * pmul * sV                       # fac < 0 → V 翻号

    return sI_out, sV_out

# 标准路径（原有代码，不变）
return sI, sV
```

**返回 pmul 供导数使用**：若后续需要解析导数，
可将返回签名改为 `(sI, sV, pmul)`；若仅用数值差分则不需要。

---

### 5.3 `localProfileAndDerivUnno.updateProfDeriv`

**新增参数** `Cm_map: np.ndarray | None = None`（形状 `(N_cells,)`）。

**修改 Stokes I/V 混合段**：

```python
# 原有
I_unscaled = sI * fI_val + self.sI0 * (1.0 - fI_val)
V_unscaled = sV * fV_val

# 替换为：
if getattr(ldata, 'fac', 1.0) < 0.0 and Cm_map is not None:
    sI, sV = _unno_profile(..., Cm=Cm_map)  # 含 CTTS 重缩放
    I_unscaled = sI   # 已在 _unno_profile 内完成混合
    V_unscaled = sV
else:
    # 原有逻辑
    I_unscaled = sI * fI_val + self.sI0 * (1.0 - fI_val)
    V_unscaled = sV * fV_val
```

**新增 `dI/dCm` 和 `dV/dCm` 导数输出**（数值差分）：

```python
if ctts_mode and calcDI:
    eps_cm = 0.01
    # 对整体 Cm_map + eps 计算差分（逐格点可用稀疏扰动优化）
    sI_p, sV_p = _unno_profile(..., Cm=Cm_map + eps_cm)
    # dI/dCm[i] 为 (N_vel,) 向量，i 为格点索引
    self.dIdCm = contin[np.newaxis,:] * (sI_p - sI) / eps_cm * inv_ic  # (N_vel, N_cells)
    self.dVdCm = contin[np.newaxis,:] * (sV_p - sV) / eps_cm * inv_ic
else:
    self.dIdCm = 0
    self.dVdCm = 0
```

**注意**：上述全局扰动实际上计算的是所有格点同时偏移，
而 MEM 需要的是 $\partial F_k / \partial C_m^i$（单格点偏移）。
正确做法是逐格点差分，但 $N_\text{cells}$ 次 `_unno_profile` 调用代价高昂。
由于 `pmul` 对各格点独立，解析导数更高效：

$$\frac{\partial I_k}{\partial C_m^i} = |fac| \cdot \frac{\partial p_\text{mul}^i}{\partial C_m^i} \cdot (I_\text{mix}^i(\lambda_k) - I_c) \cdot \text{scale}^i / I_c^\text{tot}$$

$$= |fac| \cdot \left(-\frac{10}{C_\text{def}}\right) \cdot (I_\text{mix}^i - I_c) \cdot \text{scale}^i / I_c^\text{tot}$$

形状为 `(N_cells, N_vel)`，与 `dIcsum` 结构完全一致，可直接追加到 `dIcsum` 后。

---

### 5.4 `diskIntProfAndDerivUnno`

**新增属性**：

```python
self.dIIcCm = self.prof.dIdCm   # (N_cells, N_vel) — dI/dCm，供 packResponseMatrix
self.dVdCm  = self.prof.dVdCm   # (N_cells, N_vel)
```

`updateIntProfDeriv` 中透传 `Cm_map`：

```python
def updateIntProfDeriv(self, ..., Cm_map=None):
    ...
    self.prof.updateProfDeriv(..., Cm_map=Cm_map)
    ...
    self.dIIcCm = self.prof.dIdCm
    self.dVdCm  = self.prof.dVdCm
```

---

### 5.5 `pipeline.py`

仅在 `unno` 路径的 `lineDataUnno.from_parameters(...)` 调用中增加两个参数：

```python
lineData = lineprofile.lineDataUnno.from_parameters(
    ...,                                               # 现有参数
    fac=getattr(par, 'unno_fac', 1.0),                # 新增
    Cm_default=getattr(par, 'unno_Cm_default', 0.0),  # 新增
)
```

当 `fac >= 0`（默认值 1.0）时，代码完全走标准 UR 路径，无任何行为变化。

---

### 5.6 `config.json` 新增字段

在 `"spectral_line"` 块中追加（所有字段均有缺省值，向后兼容）：

```json
"spectral_line": {
    "model_type": "unno",
    ...
    "unno_fac": -1.0,
    "unno_Cm_default": 0.0,
    "unno_fit_Cm": true
}
```

---

## 6. 实现顺序建议

1. 添加 `lineDataUnno.fac` 和 `Cm_default` 字段（低风险，不影响现有路径）
2. 在 `_unno_profile` 中实现 CTTS 分支（仅前向模型，无导数）
3. 编写单元测试，与 C 版 `zdipot23.out` / `zdipot25.out` 的合成轮廓对比
4. 在 `localProfileAndDerivUnno` 中实现 `dI/dCm`、`dV/dCm` 解析导数
5. 在 `diskIntProfAndDerivUnno` 和 `fighting.py` / `mem/` 中接通响应矩阵扩展
6. 端到端测试：以 TWA12 数据运行完整反演，与 C 版结果比较

---

## 7. 变更文件一览

| 文件 | 改动类型 | 估计行数 |
|------|----------|---------|
| `core/line_models/unno.py` | 扩展（`fac`、CTTS 分支、$\partial/\partial C_m$）| +120 行 |
| `pipeline/pipeline.py` | 参数透传（`fac`, `Cm_default`）| +4 行 |
| `config.json` | 新增可选字段 | +3 行 |
| `core/fitting.py` | 响应矩阵扩展（`Cm` 列块）| +40 行 |
| `core/mem/*.py` | 图像向量分区扩展（`Cm` 有界熵段）| +30 行 |
| **合计** | | **约 200 行** |

`pipeline.py` 以下的所有路由代码（`isinstance` 判断、`getAllProfDiriv*` 等）**均无需修改**。

---

## 8. 参考文献与代码对应关系

| 符号 | C 变量 | Python 字段/变量 | 参考文献 |
|------|--------|-----------------|---------|
| $\eta_0$ | `eta` | `ldata.str` | ME 模型线连比 |
| $\beta$ | `bb` | `ldata.beta` | Planck 函数坡度 |
| $\varepsilon$ | `eps` | `ldata.limbDark`（推导）| 热化参数 |
| $f_\text{fff}$ | `fff`（硬编码）| `ldata.fI` | 填充因子 |
| $C_m$ | `Cm[i]`（逐格点）| `Cm_map` | 吸积亮度图 |
| $p_\text{mul}$ | `pmul` | 局部变量 | — |
| $fac$ | `fac`（`< 0` = CTTS）| `ldata.fac` | — |
| $C_\text{def}$ | `defC = 0.999` | `_DEFC = 0.999` | 亮度上限 |

**C 版关键文件**（`CTTSzdi2/`）：

- `unno_ctts.c`：`unno()` + `unno_deriv()`（谱线轮廓 + 导数）
- `potential.c`：主程序，`fitcm` 控制 `Cm` 是否进入图像向量
- `resp_pot.c`：`doppl_deriv()` + `cp_image()` + `mk_image()`（响应矩阵构建）
- `mem_nl.c`：`calc_grad()` 中的 `n1..n2` 有界熵段（`Cm` 使用此段）
