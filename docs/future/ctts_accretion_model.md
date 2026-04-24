# CTTS 双 Unno-Rachkovsky 吸积模型设计方案

**参考实现**：`/Users/tianqi/Documents/Codes_collection/ZDI_and/CTTSzdi2/unno_ctts.c`（仅作兼容基线）  
**主实现模块**：`core/line_models/unno.py`  
**联动模块**：`pipeline/pipeline.py`、`config_loader.py`、`config.json`、`core/fitting.py`、`core/mem/zdi_adapter.py`、`frontend/js/config.js`（必要时还有 `api/routes/config.py` 与 `frontend/index.html`）  
**状态**：设计草案，尚未实现

---

## 1. 物理背景与动机

### 1.1 经典 T Tauri 星（CTTS）的双层线形成

经典 T Tauri 星（CTTS）处于磁控吸积阶段：来自吸积盘的物质沿磁力线下落至星面，
在吸积热斑处激发额外的色球/吸积发射分量，并与光球吸收轮廓叠加。
若色球层本身可观测到磁场，则局部轮廓不再适合写成“磁光球 + 非磁填充项”，
而应视为两个**均可磁化**的大气层共同贡献：

| 分量 | 来源 | 代表图像 | 轮廓性质 |
|------|------|----------|----------|
| 光球分量 | 恒星光球 + 表面磁场 | `Cq` | 磁 Milne-Eddington 吸收（UR） |
| 色球/吸积分量 | 热斑、后激波、磁化色球 | `Cm` | 磁 Milne-Eddington 发射/吸收（UR） |

`Cq` 主要追踪光球亮度或吸收深度，`Cm` 追踪色球/吸积层的局部发射强度或覆盖度。
两者都可逐格点变化，并与同一套表面磁几何共同决定最终 Stokes I/V 轮廓。

### 1.2 为什么放弃“π-only 非磁色球”作为主模型

旧方案把第二分量写成“仅 π 分量、无有效 Zeeman 分裂”的 Milne-Eddington 轮廓，
这在色球偏振不可测时可作为低阶近似；一旦色球磁场本身进入观测，
该近似就不再自洽，因为它默认：

1. 色球只改变 Stokes I，不直接生成自己的 Stokes V；
2. 所有偏振都来自光球磁区；
3. 色球与光球共享相同热力学参数，却不共享磁场响应。

因此，本文档将主方案改写为：

- 光球分量：完整 Unno-Rachkovsky（UR）
- 色球分量：完整 Unno-Rachkovsky（UR）
- 旧的 π-only 分支：仅保留为 `legacy_pi` 兼容模式，用于历史对比

### 1.3 建模原则

为兼顾物理自洽与实现复杂度，第一阶段采用以下原则：

1. **双 UR、单磁拓扑**：光球和色球都使用 UR 局部轮廓，但默认共用同一套表面磁场几何；
   色球磁场通过独立的响应缩放控制，而不是立刻引入第二套独立磁图。
2. **I 与 V 的比例允许不同**：Stokes I 的光球/色球权重与 Stokes V 的有效偏振权重分离。
3. **色球参数独立**：色球分量允许拥有独立的 `beta`、Landé g、线强、展宽、连续谱比例。
4. **逐步扩展**：若后续数据质量支持，再从“单磁拓扑”推广到“独立 chromospheric magnetic map”。

### 1.4 工程边界说明

这次改动的最小闭环并不是只改 `core/line_models/unno.py`。
从当前代码结构看，双 UR 模式至少会贯穿以下几层：

1. **前向模型层**：`core/line_models/unno.py`，负责参数容器、局部轮廓、盘积分和导数。
2. **流水线装配层**：`pipeline/pipeline.py`，负责根据 `line_model.model_type='unno'` 构造 `lineDataUnno` 并调用 `getAllProfDirivUnno(...)`。
3. **配置解析层**：`config_loader.py` 与 `config.json`，负责让 CLI / WebUI 都能读到新增 CTTS 参数。
4. **反演接口层**：`core/fitting.py` 与 `core/mem/zdi_adapter.py`，负责把新增的 `C_m` 导数块和双 UR 磁响应装入响应矩阵。
5. **前端 UI 层**：`frontend/js/config.js`，负责把新增参数暴露到 WebUI 表单；若字段组织方式变化，还需同步 `api/routes/config.py` 的保存/回填逻辑。

因此，下文中的“目标模块”应理解为**以 `unno.py` 为核心的一组联动模块**，而不是单文件实现任务。

---

## 2. 数学与物理原理

### 2.1 双 UR 局部轮廓公式

对每个可见表面单元，分别计算光球和色球两个 UR 局部解：

$$
I_{\rm ph}(\lambda),\ V_{\rm ph}(\lambda),\ I_{\rm ch}(\lambda),\ V_{\rm ch}(\lambda)
$$

局部总 Stokes I 写为连续谱加权平均：

$$
I_{\rm loc}(\lambda) = \frac{c_{\rm ph} I_{\rm ph}(\lambda) + c_{\rm ch} I_{\rm ch}(\lambda)}{c_{\rm ph} + c_{\rm ch}}
$$

局部总 Stokes V 写为有效偏振通量之和，再用总连续谱归一化：

$$
V_{\rm loc}(\lambda) = \frac{p_{\rm ph} V_{\rm ph}(\lambda) + p_{\rm ch} V_{\rm ch}(\lambda)}{c_{\rm ph} + c_{\rm ch}}
$$

其中：

- $c_{\rm ph}$：光球连续谱权重，建议写成 $c_{\rm ph} = C_q \cdot f_{\rm ph}^I$
- $c_{\rm ch}$：色球连续谱/发射权重，建议写成 $c_{\rm ch} = C_m \cdot f_{\rm ch}^I$
- $p_{\rm ph}$：光球有效偏振权重，建议写成 $p_{\rm ph} = C_q \cdot f_{\rm ph}^V$
- $p_{\rm ch}$：色球有效偏振权重，建议写成 $p_{\rm ch} = C_m \cdot f_{\rm ch}^V$

这里最关键的是：

$$
\frac{p_{\rm ph}}{p_{\rm ch}} \neq \frac{c_{\rm ph}}{c_{\rm ch}}
$$

这正是“Stokes I 和 V 的光球/色球比例不必相同”的数学实现。
物理上，`p` 与 `c` 分离意味着未分辨磁区抵消、不同形成高度、不同偏振效率都可以被吸收进模型。

### 2.2 光球与色球各自的 UR 分量

对任一分量 $j \in \{{\rm ph}, {\rm ch}\}$，其局部轮廓由标准 UR 方程给出：

$$
\frac{I_j}{I_{c,j}} = \frac{1 + \beta_j \mu \cdot \eta_{I,j}(\eta_{I,j}^2 + \rho_{Q,j}^2 + \rho_{V,j}^2)/\Delta_j}{1 + \beta_j \mu}
$$

$$
\frac{V_j}{I_{c,j}} = -\frac{\beta_j \mu}{1 + \beta_j \mu} \cdot \frac{\eta_{V,j}\eta_{I,j}^2 + \rho_{V,j}(\eta_{Q,j}\rho_{Q,j} + \eta_{V,j}\rho_{V,j})}{\Delta_j}
$$

其中每个分量都拥有独立参数：

- `beta_ph`, `beta_ch`
- `line_strength_ph`, `line_strength_ch`
- `gauss_width_ph_kms`, `gauss_width_ch_kms`
- `lorentz_width_ph_fraction`, `lorentz_width_ch_fraction`
- `lande_g_ph`, `lande_g_ch`
- `continuum_scale_ph`, `continuum_scale_ch`

色球分量通常满足：

- `beta_ch < 0`：温度逆转，对应发射线形成
- `limb_darkening_ch <= 0`：可能出现临边增亮
- `gauss_width_ch_kms > gauss_width_ph_kms`：色球/吸积区速度场更宽

### 2.3 单磁拓扑近似与后续扩展

第一阶段不立刻引入独立色球磁图，而是令两分量共享同一套表面磁几何：

$$
\mathbf{B}_{\rm ch} = \kappa_{\rm ch} \cdot \mathbf{B}_{\rm surf}
$$

其中 $\kappa_{\rm ch}$ 为色球磁响应缩放因子，可取常数，也可后续推广为图像量。

这样做的优点：

1. 不需要把球谐系数数量立刻翻倍；
2. 可以先验证“双 UR + 分离 I/V 权重”是否足以解释数据；
3. 保持与现有 MEM/响应矩阵接口更接近。

若未来需要独立色球磁图，则再改为：

$$
\mathbf{B}_{\rm ph} \neq \mathbf{B}_{\rm ch}
$$

此时图像向量和响应矩阵才需要显式拆成两套磁场系数块。

### 2.4 对 `Cq`、`Cm` 和磁参数的导数

双 UR 结构下，`Cm` 不再只是旧 `pmul` 的缩放器，而直接进入权重：

$$
I_{\rm loc} = \frac{c_{\rm ph} I_{\rm ph} + c_{\rm ch} I_{\rm ch}}{c_{\rm ph} + c_{\rm ch}}
$$

因此对 `Cq`、`Cm` 的导数需要同时包含：

1. 权重本身变化的贡献；
2. 归一化分母变化的贡献；
3. 若 `Cq`、`Cm` 还控制局部连续谱，则额外的连续谱项；
4. 若 `\kappa_{\rm ch}` 与 `Cm` 耦合，则还有磁响应链式导数。

对共享磁拓扑的第一阶段，响应矩阵仍可写成：

$$
R_{\rm phase1} = [R_{C_q} \mid R_{C_m} \mid R_{B_{\rm sph}}]
$$

其中：

- $R_{C_q}$：光球亮度/吸收图导数块
- $R_{C_m}$：色球发射/覆盖图导数块
- $R_{B_{\rm sph}}$：共享磁球谐系数导数块，但其内容已同时包含光球和色球两分量的 UR 响应

若未来进入独立色球磁图阶段，则升级为：

$$
R_{\rm phase2} = [R_{C_q} \mid R_{C_m} \mid R_{B_{\rm ph}} \mid R_{B_{\rm ch}}]
$$

---

## 3. 代码架构方案

### 3.1 现有路由机制保持不变

Pipeline 仍然通过 `lineData` 对象类型选择盘积分实现：

```python
# pipeline.py — 构造时选模型
if par.line_model_type == 'unno':
    lineData = lineDataUnno.from_parameters(...)
else:
    lineData = lineData.from_parameters(...)

# line_utils.py — 运行时路由（已存在）
if isinstance(lineData, lineDataUnno):
    DiskIntCls = diskIntProfAndDerivUnno
else:
    DiskIntCls = diskIntProfAndDeriv
```

CTTS 的双 UR 模式不需要改变这层路由；只需在 `lineDataUnno` 内部区分：

- `standard_unno`
- `dual_ur_ctts`
- `legacy_pi`（兼容近似）

### 3.2 `lineDataUnno` 的参数重构

现有 `lineDataUnno` 只有一组 UR 参数，不足以表达双层大气。
建议引入“分量参数”概念，例如：

```python
@dataclass
class URComponentParams:
    beta: float
    line_strength: float
    gauss_width_kms: float
    lorentz_width_fraction: float
    lande_g: float
    limb_darkening: float
    continuum_scale: float
```

`lineDataUnno` 则持有：

```python
self.photosphere: URComponentParams
self.chromosphere: URComponentParams | None
self.filling_factor_I_ph: float
self.filling_factor_I_ch: float
self.filling_factor_V_ph: float
self.filling_factor_V_ch: float
self.chromosphere_B_scale: float
self.ctts_mode: str   # 'none' | 'dual_ur' | 'legacy_pi'
```

这样做比继续堆叠 `fac`、`Cm_default` 等历史字段更清晰，也更接近真实物理。

### 3.3 `_unno_profile` 的拆分方式

建议把当前 `_unno_profile` 拆成两层：

```python
def _unno_profile_component(
    params: URComponentParams,
    wls: np.ndarray,
    view_angle: np.ndarray,
    Bmod: np.ndarray,
    Btheta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    ...

def _unno_profile_ctts(
    ldata: lineDataUnno,
    wls: np.ndarray,
    view_angle: np.ndarray,
    Bmod: np.ndarray,
    Btheta: np.ndarray,
    Cq: np.ndarray,
    Cm: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    ...
```

其中：

1. `_unno_profile_component` 只负责单一 UR 分量求解；
2. `_unno_profile_ctts` 调两次 `_unno_profile_component`，分别得到光球和色球轮廓；
3. 最后按照 `c_ph`, `c_ch`, `p_ph`, `p_ch` 组合成总 `I` 和 `V`。

`legacy_pi` 则作为 `_unno_profile_ctts(..., mode='legacy_pi')` 的一个退化分支保留。

### 3.4 `localProfileAndDerivUnno.updateProfDeriv`

该函数需从“单分量 + 混合填充因子”改成“显式双分量求和”：

```python
Iph, Vph = _unno_profile_component(... photosphere ...)
Ich, Vch = _unno_profile_component(... chromosphere ...)

c_ph = bright_q * fI_ph
c_ch = bright_m * fI_ch
p_ph = bright_q * fV_ph
p_ch = bright_m * fV_ch

I_unscaled = (c_ph * Iph + c_ch * Ich) / (c_ph + c_ch)
V_unscaled = (p_ph * Vph + p_ch * Vch) / (c_ph + c_ch)
```

随后导数需要同步更新：

- `dI/dCq`
- `dI/dCm`
- `dV/dCq`
- `dV/dCm`
- `dI/dB`
- `dV/dB`

与旧方案不同，`dV/dCm` 不再只是标量振幅导数，而包含色球磁分量本身的偏振响应。

### 3.5 图像向量与 MEM 的阶段划分

第一阶段建议**不修改磁场图像向量维度**，仅保留：

```text
image = [ Cq_0..N, Cm_0..N, Br/Bt/Bp spherical coefficients ]
```

这样能把开发风险限制在前向模型和导数组合上。

第二阶段若需要独立色球磁图，再升级为：

```text
image = [ Cq_0..N, Cm_0..N, Bph_sph_coeff, Bch_sph_coeff ]
```

因此本文档把“双 UR”与“双磁图”明确区分：

- **当前要实现的是双 UR**
- **独立双磁图是后续扩展，不是本次最小闭环**

### 3.6 `pipeline.py` 与配置系统

`pipeline.py` 本身不需要改动模型选择框架，但仍属于**必须同步修改**的链路节点：

1. `pipeline/pipeline.py` 需要把新的色球参数透传给 `lineDataUnno.from_parameters(...)`；
2. `config_loader.py` 需要把新增 `line_model.*` 字段解析为 `par` 属性；
3. `config.json` 需要提供默认值；
4. `frontend/js/config.js` 需要把这些字段渲染到 WebUI；
5. 若前端采用新的嵌套结构或别名，还要同步 `api/routes/config.py` 的读写映射。

建议在 `config.json` / `config_loader.py` / `frontend/js/config.js` 中新增：

```json
"line_model": {
  "model_type": "unno",
  "ctts_mode": "dual_ur",
  "chromosphere_beta": -0.5,
  "chromosphere_line_strength": 0.6,
  "chromosphere_gauss_width_kms": 20.0,
  "chromosphere_lorentz_width_fraction": 0.02,
  "chromosphere_lande_g": 1.0,
  "chromosphere_limb_darkening": -0.3,
  "chromosphere_continuum_scale": 1.0,
  "filling_factor_I_ph": 1.0,
  "filling_factor_I_ch": 1.0,
  "filling_factor_V_ph": 1.0,
  "filling_factor_V_ch": 0.5,
  "chromosphere_B_scale": 1.0,
  "legacy_pi_mode": false
}
```

其中：

- `filling_factor_I_*` 控制 Stokes I 的层权重
- `filling_factor_V_*` 控制 Stokes V 的有效偏振权重
- `chromosphere_B_scale` 控制色球分量对共享磁图的响应强度

---

## 4. 预期实现目标

1. **主模型升级为双 UR**：光球和色球都用 UR 局部解，不再默认色球非磁。
2. **I/V 权重可分离**：允许色球在 Stokes I 中强、在 Stokes V 中弱，或反之。
3. **`Cm` 同时进入 I 和 V**：`Cm` 不再只是吸收深度的经验缩放器，而直接参与色球分量的连续谱和偏振权重。
4. **保持向后兼容**：当 `ctts_mode='none'` 或所有色球权重为 0 时，退化为现有标准 UR；当 `legacy_pi_mode=true` 时，保留历史近似用于对照。
5. **第一阶段不翻倍磁自由度**：先共用同一套磁场几何，验证物理收益；必要时再扩展为独立色球磁图。

---

## 5. 具体改动方案

### 5.1 `core/line_models/unno.py`

#### 5.1.1 参数层

- 为 `lineDataUnno` 增加 `photosphere` 与 `chromosphere` 两套参数；
- 新增 `ctts_mode`、`filling_factor_I_ph/ch`、`filling_factor_V_ph/ch`、`chromosphere_B_scale`；
- 将旧 `fac`、`Cm_default` 归入 `legacy_pi` 路径，不再作为主路径核心参数。

#### 5.1.2 轮廓层

- 抽取 `_unno_profile_component()` 作为单分量求解器；
- 新增 `_unno_profile_ctts()` 负责双分量组合；
- `legacy_pi` 分支仅在显式启用时调用旧逻辑。

#### 5.1.3 导数层

- `localProfileAndDerivUnno` 输出 `dI/dCq`、`dI/dCm`、`dV/dCq`、`dV/dCm`；
- 共享磁图阶段，`dI/dB`、`dV/dB` 已自然包含光球和色球两层贡献；
- 若 `chromosphere_B_scale` 可调，也可为其增加全局导数，用于后续参数拟合。

### 5.2 `diskIntProfAndDerivUnno`

需要新增或重定义以下输出：

```python
self.dIIcQ    # dI/dCq
self.dIIcM    # dI/dCm
self.dVdQ     # dV/dCq
self.dVdM     # dV/dCm
self.dVIc     # dV/d(mag coeff)
self.dImag    # dI/d(mag coeff)
```

其中 `dIIcM` 与 `dVdM` 将直接接入 MEM 响应矩阵的 `Cm` 列块。

### 5.3 `core/fitting.py` 与 `core/mem/`

第一阶段的响应矩阵扩展为：

```text
R = [ R_Cq | R_Cm | R_Bsph ]
```

与旧方案相比，变化在于：

- `R_Cm` 现在同时影响 I 和 V；
- `R_Bsph` 现在包含双 UR 的联合磁响应；
- 熵函数部分仍可保持 `Cq` 与 `Cm` 两段图像块设计。

`Cm` 仍建议使用有界熵，因为它表示色球发射强度或覆盖度，应满足正定和上限约束。

这里建议把 `core/mem/` 的影响进一步具体化为 `core/mem/zdi_adapter.py`：

- `packImageVector()` / `unpackImageVector()` 需要确认 `Cq`、`Cm` 图像块与现有亮度图像布局兼容；
- `packResponseMatrix()` 需要明确 `R_Cm` 进入 I 与 V 两个观测块，而不是只影响 I；
- 若 `Cm` 继续使用填充因子熵，还要确认 `n1/n2` 分段与权重设置仍然正确。

### 5.4 `pipeline.py`

只需透传新增参数，不改变模型选择逻辑：

```python
lineData = lineprofile.lineDataUnno.from_parameters(
    ...,  # 现有光球参数
    ctts_mode=getattr(par, 'ctts_mode', 'none'),
    chromosphere_beta=getattr(par, 'chromosphere_beta', -0.5),
    chromosphere_line_strength=getattr(par, 'chromosphere_line_strength', 0.0),
    chromosphere_gauss_width_kms=getattr(par, 'chromosphere_gauss_width_kms', 20.0),
    chromosphere_lorentz_width_fraction=getattr(par, 'chromosphere_lorentz_width_fraction', 0.02),
    chromosphere_lande_g=getattr(par, 'chromosphere_lande_g', 1.0),
    filling_factor_I_ph=getattr(par, 'filling_factor_I_ph', 1.0),
    filling_factor_I_ch=getattr(par, 'filling_factor_I_ch', 1.0),
    filling_factor_V_ph=getattr(par, 'filling_factor_V_ph', 1.0),
    filling_factor_V_ch=getattr(par, 'filling_factor_V_ch', 0.5),
    chromosphere_B_scale=getattr(par, 'chromosphere_B_scale', 1.0),
)
```

### 5.5 配置与前端

这部分不是外围收尾，而是第一阶段的**必做项**。需要同步更新：

- `config.json`
- `config_loader.py`
- `frontend/js/config.js`
- 如涉及 API 保存/读取字段映射，再同步 `api/routes/config.py`
- 若前端静态资源版本号需要刷新，再同步 `frontend/index.html`

其中前端 UI 需要明确解决两件事：

1. 让 `ctts_mode='none' | 'dual_ur' | 'legacy_pi'` 驱动条件显隐；
2. 将光球与色球参数分组展示，避免把新参数直接堆入现有单层 UR 字段区，造成语义混淆。

对外暴露的最低必要字段是：

1. `ctts_mode`
2. `chromosphere_beta`
3. `chromosphere_line_strength`
4. `chromosphere_gauss_width_kms`
5. `chromosphere_lande_g`
6. `filling_factor_I_ph`
7. `filling_factor_I_ch`
8. `filling_factor_V_ph`
9. `filling_factor_V_ch`
10. `chromosphere_B_scale`

### 5.6 兼容模式

为了与历史结果对照，建议保留两个兼容开关：

```text
ctts_mode = 'none'       -> 现有标准 UR
ctts_mode = 'dual_ur'    -> 新主模型
ctts_mode = 'legacy_pi'  -> 旧 π-only 非磁色球近似
```

这样既能做 A/B 测试，也能量化“色球磁化”对拟合的提升幅度。

### 5.7 模式对应关系与退化条件

这三种模式之间应明确区分“严格等价”和“参数极限下近似等价”：

| 目标行为 | 推荐模式 | 是否严格等价 | 说明 |
|----------|----------|--------------|------|
| 现有标准单层 UR | `ctts_mode='none'` | 是 | 直接关闭色球分量，完全退化回当前 `lineDataUnno` 路径 |
| 历史 `π-only` 非磁色球近似 | `ctts_mode='legacy_pi'` | 是 | 保留旧公式与旧组合方式，用于回归测试与历史结果复现 |
| 双 UR 主模型 | `ctts_mode='dual_ur'` | 是 | 光球与色球都使用 UR 分量，是新的主物理模型 |
| 用双 UR 模式逼近旧 `π-only` 行为 | `ctts_mode='dual_ur'` + 参数约束 | 否，仅近似 | 可让色球几乎不产生 Zeeman 分裂和 Stokes V，但代数上不等同于“只保留 π 分量” |

原因在于：`dual_ur` 的色球分量仍然通过完整 UR 方程计算，
即使令 Zeeman 分裂趋近于零，它也会退化成“零分裂极限下的完整非磁 UR”，
而不是字面意义上从不透明度表达式中删去两个 `\sigma` 分量。
因此，若目标是逐公式复现旧模型，必须保留 `legacy_pi`；
若目标只是让新架构在数值行为上逼近旧结果，则可通过参数极限实现。

#### 5.7.1 `dual_ur` 到 `none` 的严格退化

当满足以下条件时，`dual_ur` 应与 `none` 模式数值一致：

```text
filling_factor_I_ch = 0
filling_factor_V_ch = 0
chromosphere_line_strength = 0
```

更稳妥的实现方式是：一旦检测到色球分量所有权重为零，就直接短路到现有单层 UR 路径。
这样可以避免数值噪声，也让回归测试更直接。

#### 5.7.2 `dual_ur` 逼近 `legacy_pi` 的推荐参数

若希望在不切换到 `legacy_pi` 分支的前提下，让双 UR 行为尽可能接近旧的“仅 `\pi` 分量、无有效 Zeeman 分裂”色球近似，可使用：

```text
ctts_mode = 'dual_ur'
chromosphere_B_scale = 0
chromosphere_lande_g = 0
filling_factor_V_ch = 0
beta_ch < 0
chromosphere_line_strength = 重新标定
```

这组参数的含义是：

1. `chromosphere_B_scale = 0`：色球不响应共享磁场；
2. `chromosphere_lande_g = 0`：即使实现层未显式短路，也让 Zeeman 分裂量趋于零；
3. `filling_factor_V_ch = 0`：强制色球不直接贡献 Stokes V；
4. `beta_ch < 0`：保留色球发射层的热力学特征；
5. `chromosphere_line_strength` 需要重新拟合，因为“零分裂极限的完整 UR”与“字面 `π-only` 近似”在线强归一上 generally 不完全一致。

#### 5.7.3 为什么这只能是近似，不是严格相同

旧 `legacy_pi` 路径的本质是“人为指定第二分量只含 `\pi` 线”；
而 `dual_ur` 的零分裂极限本质是“`\pi` 与两个 `\sigma` 分量在同一线心完全并合”。

因此二者会在以下方面保留差异：

1. 局部吸收系数的组合权重；
2. 线强参数与连续谱归一化的等效定义；
3. 极端参数下的线心深度/峰强标定；
4. 若实现保留 `rho` 项或数值差分导数，导数结构也不完全相同。

结论是：

- **要严格复现历史 `π-only` 行为，用 `legacy_pi`**
- **要让新架构平滑逼近旧行为，用 `dual_ur` 参数极限**
- **要做物理分析和正式拟合，优先使用 `dual_ur`**

---

## 6. 实现顺序建议

1. 抽取 `_unno_profile_component()`，保证现有单 UR 路径完全不变。
2. 为 `lineDataUnno` 加入色球参数和 `ctts_mode='dual_ur'`。
3. 在前向模型中实现双 UR 组合，但先不接导数。
4. 用合成数据验证：当 `filling_factor_I_ch=filling_factor_V_ch=0` 时，必须退化回现有 UR。
5. 在 `localProfileAndDerivUnno` 中接入 `dI/dCq`、`dI/dCm`、`dV/dCq`、`dV/dCm`。
6. 扩展 `core/fitting.py` / `core/mem/zdi_adapter.py`，把 `Cm` 列块接入响应矩阵。
7. 同步更新 `config_loader.py`、`config.json`、`frontend/js/config.js`，打通 CLI 与 WebUI 的参数链路。
8. 若 UI 字段组织发生变化，再同步 `api/routes/config.py` 与 `frontend/index.html` 的缓存版本号。
9. 加入 `legacy_pi` 对照测试，量化双 UR 相对旧近似的改进。
10. 如数据确有需求，再进入第二阶段：独立色球磁图。

---

## 7. 变更文件一览（最小实现闭环）

| 文件 | 改动类型 | 估计行数 |
|------|----------|---------|
| `core/line_models/unno.py` | 双 UR 参数容器、单分量/双分量求解器、导数扩展 | +180 行 |
| `pipeline/pipeline.py` | 新参数透传 | +15 行 |
| `config_loader.py` | 新字段解析 | +20 行 |
| `config.json` | 新增色球 UR 配置字段 | +12 行 |
| `frontend/js/config.js` | 新配置项渲染、条件显隐与字段分组 | +30 行 |
| `api/routes/config.py` | 如前后端字段结构不一致时，补充保存/回填映射 | +0~20 行 |
| `frontend/index.html` | 如修改前端资源，更新缓存版本号 | +1 行 |
| `core/fitting.py` | 迭代流程中接入双 UR 响应块 | +20 行 |
| `core/mem/zdi_adapter.py` | `C_m` 图像块、I/V 响应矩阵拼装与熵分段适配 | +30~60 行 |
| `tests/*` | 退化条件、兼容模式和响应矩阵回归测试 | +40~80 行 |
| `docs/line_models.md` | 对外说明 Unno 模型已扩展到 CTTS 双 UR | +15 行 |
| **合计** | | **约 350~440 行** |

与此相比，`line_utils.py` 的运行时类型路由预计无需修改；但这不改变本方案是一个**跨核心层与前端层的联动改动**这一事实。

---

## 8. 参考文献与变量对应关系

| 物理量 | 旧近似中的角色 | 新双 UR 结构中的角色 |
|--------|----------------|----------------------|
| `Cq` | 光球亮度图 | 光球连续谱/吸收权重图 |
| `Cm` | 经验发射缩放或填充图 | 色球连续谱/发射权重图 |
| `pmul` | 旧 C 版经验振幅缩放 | 仅保留在 `legacy_pi`，不再是主模型核心 |
| `fI` | 单一 Stokes I 填充因子 | 拆分为 `filling_factor_I_ph/ch` |
| `fV` | 单一 Stokes V 填充因子/Zeeman 修正 | 拆分为 `filling_factor_V_ph/ch` |
| `Bsurf` | 唯一磁场图 | 第一阶段共享于光球和色球 |
| `kappa_ch` | 无 | 色球磁响应缩放 |
| `beta` | 单层 ME 坡度 | 拆分为 `beta_ph` 与 `beta_ch` |

**参考代码与文档**：

- `CTTSzdi2/unno_ctts.c`：旧 CTTS 近似实现，仅用于回归对照
- `docs/future/chromospheric_zdi.md`：色球 UR 发射线的单分量物理基础
- `docs/line_models.md`：当前项目中的 Unno-Rachkovsky 框架说明

---

## 9. 结论

这一步的核心不是“在旧 π-only 框架上打补丁”，而是把 CTTS 主模型正式改写为：

**光球磁分量 + 色球磁分量** 的双 Unno-Rachkovsky 结构。

在这个结构下：

1. 色球磁场可以直接贡献自己的 Stokes V；
2. Stokes I 和 V 的光球/色球比例可以不同；
3. `Cm` 既影响发射强度，也影响偏振权重；
4. 当前代码架构仍可在不改变 `lineData` 路由方式的前提下逐步接入；
5. 真正复杂的“独立色球磁图”被明确推迟到第二阶段，而不是阻塞当前实现。

这使得模型既比旧方案物理上更自洽，也仍然保留一个可控、可实现、可回归验证的开发路径。
