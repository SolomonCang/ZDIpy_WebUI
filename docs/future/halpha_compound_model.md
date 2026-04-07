# H-alpha 双峰发射线 ZDI — 复合轮廓模型方案（路径 3）

**状态**：设计草案，尚未实现  
**关联文档**：`docs/future/chromospheric_zdi.md`（纯色球 ZDI 通用框架）

---

## 1. 问题背景

H-alpha（656.28 nm，$g_\text{eff} = 1.048$）在活跃恒星（T Tauri、dMe 等）上
通常表现为**色球发射 + 中心自吸收双峰**轮廓。
标准 Voigt 单峰模型和 Milne-Eddington Unno-Rachkovsky 模型均假设单峰线轮廓，
不能直接拟合该形态。

物理图像：

- **宽发射成分**：色球/过渡区热等离子体（$T \sim 10^4$ K），发射峰半宽约 50–200 km/s
- **窄自吸收成分**：冷色球/光球层对发射光子的再吸收，抑制线心约 20–80 km/s

两个成分形成于不同高度，但均受同一地表磁场调制，因此原则上均携带 Zeeman 信息。

---

## 2. 复合轮廓模型（双 Voigt）

### 2.1 局部 Stokes I 公式

$$
I_\text{local}(\lambda) = 1
  + A_\text{em} \cdot H(a_\text{em},\, x_\text{em})
  - A_\text{abs} \cdot H(a_\text{abs},\, x_\text{abs})
$$

其中

| 符号 | 含义 |
|------|------|
| $A_\text{em} > 0$ | 宽发射成分强度（连续谱单位）|
| $A_\text{abs} > 0$ | 窄自吸收成分强度 |
| $a_\text{em},\, a_\text{abs}$ | 各分量 Lorentz-to-Gauss 比 |
| $H(a, x)$ | Humlicek Voigt 函数（实部）|
| $x_c = (\lambda - \lambda_\text{Doppler}) / \Delta\lambda_{D,c}$ | 各分量归一化坐标 |

轮廓在无磁场时呈双峰（$A_\text{em}$ 宽，$A_\text{abs}$ 窄，中心挖坑）；
两个分量的 Doppler 中心均跟踪格点旋转速度。

### 2.2 局部 Stokes V（弱场近似，两分量叠加）

H-alpha 的磁场分裂通常远小于谱线宽度（$g_\text{eff} \approx 1.05$，弱场近似成立），
因此对每个分量独立应用弱场近似再叠加：

$$
V_\text{local}(\lambda) = -\Delta\lambda_Z \cdot B_\text{LOS}
  \left[
    A_\text{em} \cdot \frac{\partial H_\text{em}}{\partial\lambda}
    - A_\text{abs} \cdot \frac{\partial H_\text{abs}}{\partial\lambda}
  \right]
$$

其中 $\Delta\lambda_Z = 4.6686\times10^{-12}\, g_\text{eff}\, \lambda_0^2$（nm/G）。

> **注意**：两个 Voigt 导数的符号相反（一正一负），会在靠近线心处部分抵消，
> 导致 V 振幅低于纯发射/纯吸收单分量模型。这是物理上正确的结果。

### 2.3 关于 Milne-Eddington 框架的适用性

Unno-Rachkovsky 模型的 ME 假设（线性源函数 $S = B_0(1 + \beta\tau)$）不能描述
双峰双层大气。本设计故意退回到弱场 Voigt 叠加，舍弃 ME 约束，换取
对双峰形状的直接参数化控制。

若未来需要完整偏振转移，需引入双层 ME 大气（Landi Degl'Innocenti &
Landolfi 2004, §11）——设计复杂度远高于本方案，留作未来扩展。

---

## 3. 参数设计

### 3.1 `lineDataHalpha` 容器类

新建 `core/line_models/halpha.py`，参数容器继承或平行于 `lineDataUnno`：

```python
class lineDataHalpha:
    """H-alpha 双峰复合 Voigt 线参数容器。

    字段（均为长度-1 的 ndarray，保持与 lineData/lineDataUnno 接口一致）：
    ---------------------------------------------------------------
    wl0             : 线心波长（nm），固定 656.28
    g               : 有效 Landé 因子，H-alpha 取 1.048
    limbDark        : 临边昏暗系数（H-alpha 色球层通常取 0 或负值）
    gravDark        : 重力昏暗系数（通常 0）
    instRes         : 仪器分辨率 R = λ/Δλ（< 0 = 不卷积）

    发射成分（em）：
    A_em            : 发射峰强度（连续谱单位，> 0）
    widthGauss_em   : Gaussian 半宽（km/s）
    widthLorentz_em : Lorentz/Gauss 比（无量纲）

    自吸收成分（abs）：
    A_abs           : 自吸收深度（0 = 无自吸收，> 0 = 双峰）
    widthGauss_abs  : Gaussian 半宽（km/s）
    widthLorentz_abs: Lorentz/Gauss 比

    fV              : 整体 Stokes V 填充因子（默认 1.0）
    """
```

**参数初始化建议（从平均轮廓预拟合）**：

用 `scipy.optimize.curve_fit` 对所有观测相位平均轮廓拟合上述双 Voigt 形式，
将拟合结果固定为反演的线参数（固定轮廓形状，只优化 B 和亮度）。
这保证了 ZDI 反演中轮廓形状不会自由浮动，避免过拟合。

### 3.2 从 `config.json` 构造

在 `config.json` 的 `line_model` 节增加：

```json
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
  "filling_factor_V": 1.0,
  "inst_resolution": 65000.0
}
```

---

## 4. 代码实现规划

### 4.1 新建文件：`core/line_models/halpha.py`

**关键函数**：

```python
_C_KMS: float = 2.99792458e5
_ZEEMAN_CONST: float = 4.66864e-12  # nm/G


def _halpha_compound_profile(
    ldata: "lineDataHalpha",
    wls: np.ndarray,    # (N_vel, N_cells) 已 Doppler 移位
    Blos: np.ndarray,   # (N_cells,)  视线磁场分量（G）
) -> tuple[np.ndarray, np.ndarray]:
    """返回 (sI, sV)，形状均为 (N_vel, N_cells)。"""

    wl0   = float(ldata.wl0[0])
    g     = float(ldata.g[0])
    A_em  = float(ldata.A_em[0])
    A_abs = float(ldata.A_abs[0])
    wG_em  = float(ldata.widthGauss_em[0])
    wG_abs = float(ldata.widthGauss_abs[0])
    aL_em  = float(ldata.widthLorentz_em[0])
    aL_abs = float(ldata.widthLorentz_abs[0])
    fV     = float(ldata.fV[0])

    dZ = _ZEEMAN_CONST * g * wl0**2    # nm/G，Zeeman 分裂系数

    # --- 发射成分 ---
    dwl_em  = wG_em  / _C_KMS * wl0   # Gaussian 宽度（nm）
    x_em    = (wls - wl0) / dwl_em    # (N_vel, N_cells)
    W_em    = _voigt_faraday_humlicek(x_em, aL_em)
    norm_em = 1.0 / (np.sqrt(np.pi) * dwl_em)
    H_em    = W_em.real * norm_em
    # 解析导数：dH/dλ = (1/Δλ_D) * 2*(x*H - F)  （Humlicek w4 虚部即 F）
    dH_em_dlambda = (2.0 * (x_em * W_em.real - W_em.imag) / dwl_em) * norm_em

    # --- 自吸收成分 ---
    dwl_abs  = wG_abs / _C_KMS * wl0
    x_abs    = (wls - wl0) / dwl_abs
    W_abs    = _voigt_faraday_humlicek(x_abs, aL_abs)
    norm_abs = 1.0 / (np.sqrt(np.pi) * dwl_abs)
    H_abs    = W_abs.real * norm_abs
    dH_abs_dlambda = (2.0 * (x_abs * W_abs.real - W_abs.imag) / dwl_abs) * norm_abs

    sI = 1.0 + A_em * H_em - A_abs * H_abs

    # 弱场近似 Stokes V：V = -dZ * Blos * dI/dλ
    dI_dlambda = A_em * dH_em_dlambda - A_abs * dH_abs_dlambda
    sV = -dZ * fV * Blos[np.newaxis, :] * dI_dlambda

    return sI, sV
```

#### Voigt 函数解析导数推导

对 Humlicek `w4(z)` 其中 $z = a - ix$（$a$ = Lorentz/Gauss 比，$x$ = 归一化坐标）：

$$
\frac{\partial H}{\partial x} = \frac{\partial \operatorname{Re}[w_4]}{\partial x}
= 2\bigl(x \cdot H(a,x) - F(a,x)\bigr)
$$

其中 $F = \operatorname{Im}[w_4]$ 为 Faraday-Voigt 函数。
该关系由 Faddeeva 函数微分恒等式 $w'(z) = -2zw(z) + 2i/\sqrt{\pi}$ 的实部导出，
与 del Toro Iniesta (2003) §5.4 一致。

因此 $\partial H / \partial\lambda = (\partial H / \partial x) / \Delta\lambda_D$，
代入上式即得代码中的 `dH_dlambda`，**复用 `_voigt_faraday_humlicek` 的输出，
无需额外函数调用**。

### 4.2 `localProfileAndDerivHalpha` 类

接口与 `localProfileAndDerivUnno` 一致：

- `__init__`：预计算 B=0 时的参考轮廓 `sI0`（`A_abs > 0` 时即为双峰形状）
- `updateProfDeriv`：
  - 调用 `_halpha_compound_profile` 计算 `(sI, sV)`
  - `dV/dBlos` 解析可得（见 §4.1），无需数值差分
  - `dV/dCoeff = (dV/dBlos) * (dBlos/dCoeff)`，形状 `(n_types, nTot, N_vel)`
    → `np.einsum('jl,kl->jk', dBlos_dCoeff, dV_dBlos_per_cell)`
  - `dI/dBright` 与 Voigt/Unno 模型完全相同

与 `localProfileAndDerivUnno` 相比，由于弱场近似下 `dV/dBlos` 解析（而非数值），
每次迭代节省两次 `_halpha_compound_profile` 调用（约 30% 加速）。

### 4.3 `diskIntProfAndDerivHalpha` 类

与 `diskIntProfAndDerivUnno` 结构相同，区别：

- `_BlosProjected` 替代 `_BdBprojected`：只需视线分量
  ```python
  def _BlosProjected(self, v_view, v_mag_cart, d_mag_cart):
      # Blos = dot(v_view, v_mag_cart), 逐格点
      Blos = np.einsum('ij,ij->j', v_view, v_mag_cart)
      if self.calcDV == 1:
          dBlos_d = np.einsum('il,ijkl->jkl', v_view, d_mag_cart, optimize=True)
      else:
          dBlos_d = 0
      return Blos, dBlos_d
  ```
- 仪器展宽：直接复用 `convolveIGnumpy`（与 Unno 版本完全相同）

### 4.4 注册到 pipeline

在 `config_loader.py` 中，当 `model_type == "halpha_compound"` 时构造
`lineDataHalpha`；`getAllProfDirivHalpha` 函数接口与 `getAllProfDirivUnno` 对齐。

`pipeline/pipeline.py` 对 `diskIntProfAndDeriv*` 对象仅通过以下属性访问，
**无需修改**：

| 属性 | 形状 | 含义 |
|------|------|------|
| `IIc` | `(N_vel,)` | 归一化 Stokes I |
| `VIc` | `(N_vel,)` | 归一化 Stokes V |
| `dVIc` | `(n_types, nTot, N_vel)` | dV/d(mag coeff) |
| `dIIc` | `(N_cells, N_vel)` | dI/d(brightness) |
| `dVdBri` | `(N_cells, N_vel)` | dV/d(brightness) |

---

## 5. 参数确定工作流

```
1. 读取全部相位观测文件（LSD 格式，.prof.norm）
   ↓
2. 计算相位平均轮廓 <I(λ)>
   ↓
3. scipy.optimize.curve_fit 拟合双 Voigt：
      I_fit = 1 + A_em * H(a_em, x_em) - A_abs * H(a_abs, x_abs)
   → 得到 6 个参数：A_em, wG_em, aL_em, A_abs, wG_abs, aL_abs
   ↓
4. 写入 config.json 对应字段（或单独保存为 halpha_params.json）
   ↓
5. 运行 ZDI 反演（固定轮廓形状，优化 B 系数 + 亮度图）
```

辅助脚本建议放置于 `scripts/fit_halpha_template.py`，
输出结果既可打印也可直接更新 config.json 对应字段。

---

## 6. 已知局限与风险

| 局限 | 说明 |
|------|------|
| 弱场假设 | $g_\text{eff} \approx 1.05$，H-alpha 在 5 kG 以内弱场近似误差 < 5%；M 矮星典型场强可达 1–3 kG，可接受 |
| 双峰分离依赖 $v\sin i$ | 若 $v\sin i \lesssim$ 自吸收宽度（< 20 km/s），双峰不可分辨，模型退化为单宽峰（仍可工作，令 `A_abs ≈ 0`）|
| Stokes V 信噪比极低 | H-alpha V 振幅约为光球 LSD 的 1/10–1/5，实际约束力有限；建议与 LSD 联合反演（见 chromospheric_zdi.md §4）|
| 固定轮廓形状 | 双峰深度随活动水平变化；若各相位轮廓形状差异大，固定模板引入系统误差 |
| 色球不均匀性 | 真实色球 H-alpha 轮廓依赖局部温度/密度结构，与"全球均匀双 Voigt"假设矛盾；可用亮度图部分弥补 |
| ME 框架缺失 | 放弃了 Milne-Eddington 源函数约束，磁场对 Stokes I 的贡献（$dI/dB$）未建模，因此联合亮度+磁场精度低于 Unno 模型 |

---

## 7. 测试策略

1. **合成测试**：用已知磁场（纯偶极，$B_\text{pol} = 500$ G）生成合成
   H-alpha Stokes I/V，验证模型能以 < 5% 误差恢复
2. **单元测试**（建议放于 `tests/test_halpha_model.py`）：
   - `_halpha_compound_profile` 在 B=0 时返回关于线心对称的双峰
   - Stokes V 在 B→0 时趋于 0
   - `dV/dBlos` 与有限差分数值导数一致（相对误差 < 1e-4）
   - `convolveIGnumpy` 对 I 和 V 的宽度均正确展宽
3. **与 LSD 联合反演对比**：用同一颗活跃恒星的光球 LSD 反演结果作为基准，
   对比 H-alpha 单独约束的磁场图，评估定性一致性

---

## 8. 参考文献

- Donati et al. (1997) MNRAS 291, 658 — ZDI 原始公式，Stokes V 弱场近似
- Johns-Krull et al. (2004) ApJ 616, L139 — T Tauri 星 H-alpha 偏振观测
- Zaire et al. (2022) MNRAS 517, 3 — 联合光球/色球 ZDI 实践
- Landi Degl'Innocenti & Landolfi (2004) §5.9 — 复合线轮廓 Stokes 参数
- del Toro Iniesta (2003) *Introduction to Spectropolarimetry* §5.4 — Faddeeva 函数导数恒等式
- Humlicek (1982) JQSRT 27, 437 — Voigt 函数高效近似算法
