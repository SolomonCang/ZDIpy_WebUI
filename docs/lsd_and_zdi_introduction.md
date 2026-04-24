# LSD 与 ZDI 方法介绍

本文档介绍了恒星偏振光谱观测中两个核心技术方法：最小二乘反卷积（LSD）和塞曼多普勒成像（ZDI），及其在 ZDIpy 中的应用。

## 目录

1. [物理背景](#物理背景)
2. [最小二乘反卷积（LSD）](#最小二乘反卷积lsd)
3. [塞曼多普勒成像（ZDI）](#塞曼多普勒成像zdi)
4. [完整工作流](#完整工作流)
5. [软件工具](#软件工具)

---

## 物理背景

### 塞曼效应与偏振信号

恒星磁场通过**塞曼效应**作用于谱线，使能级分裂为若干子能级。在外磁场 $B$ 作用下，能级分裂量为：

$$
\Delta E = g_J \mu_B m_J B
$$

其中 $g_J$ 为朗德因子，$\mu_B$ 为玻尔磁子。这导致谱线分裂为 $\sigma^+$、$\pi$、$\sigma^-$ 三个分量。

#### 塞曼分裂量

分裂产生的波长差为：

$$
\Delta\lambda_B = 4.67 \times 10^{-13}\, g_{\text{eff}}\, \lambda^2\, B
$$

其中：
- $\lambda$ 为波长（nm）
- $g_{\text{eff}}$ 为有效朗德因子
- $B$ 为磁场强度（Gauss）
- $\Delta\lambda_B$ 为波长位移（nm）

**典型量级**：对于 G/K 型星（$B \sim 100$ G，$\lambda \sim 600$ nm，$g_{\text{eff}} \sim 1.2$），分裂量 $\Delta\lambda_B \sim 0.02$ nm，低于高分辨率光谱仪的像素分辨率（~0.01 nm）。因此分裂体现为**偏振**信号，而非强度信号。

### 圆偏振信号的起源

**纵向塞曼效应**（磁场沿视线方向）：
- $\sigma^+$ 和 $\sigma^-$ 成分分别为右旋和左旋圆偏振
- 在弱场近似下，圆偏振参数 Stokes V 的相对大小为：

$$
\frac{V}{I} \approx -g_{\text{eff}} \frac{\Delta\lambda_B}{\lambda} \cos\theta_B \cdot \frac{d\ln I}{d\ln\lambda}
$$

其中 $\theta_B$ 为磁场与视线的夹角。

- Stokes V 轮廓呈**反对称的"S"形**，对纵向磁场分量 $B_\parallel = B\cos\theta_B$ 敏感

**横向塞曼效应**（磁场在天空平面）：
- 产生线偏振 Stokes Q 和 U
- 信号幅度 $\propto (B\sin\theta_B)^2$，比 Stokes V 弱 1-2 个量级
- 在恒星盘积分观测中**几乎完全抵消**

### 盘积分与大尺度场的可见性

恒星不同区域的磁场方向各异。盘积分观测中：

- **小尺度随机场**：正向和负向磁通完全相消，不可见
- **大尺度有序场**（如偶极场）：南北半球不完全对称，可测到净圆偏振信号

这是为什么 ZDI 只能重建 $\ell \lesssim 10$-15 的大尺度磁场结构。

---

## 最小二乘反卷积（LSD）

### 问题提出

单条谱线的 Stokes V 信噪比低：
- Stokes V 信号幅度：$\sim 10^{-3}$
- 单条谱线噪声：$\sigma \sim 10^{-3}$-$10^{-2}$
- 信噪比远不足以可靠探测

**解决方案**：假设所有谱线的偏振轮廓形态相似，通过**加权叠加数千条谱线**来提升 SNR。

### 数学框架

#### 自相似假设

在弱场近似下，假设所有谱线的圆偏振轮廓可用一个统一的"平均轮廓"缩放表示：

$$
V(\lambda) \approx \sum_k w_k \cdot Z(\lambda - \lambda_k)
$$

其中：
- $V(\lambda)$ 为观测到的 Stokes V 谱
- $Z(v)$ 为待求的"平均偏振轮廓"（LSD 轮廓，速度空间）
- $w_k = g_{\text{eff},k} \cdot \lambda_k \cdot d_k$ 为谱线权重
  - $g_{\text{eff},k}$ 为第 $k$ 条谱线的有效朗德因子
  - $\lambda_k$ 为谱线中心波长
  - $d_k$ 为谱线深度（归一化强度）

#### 矩阵反演

转换为矩阵形式进行最小二乘反演：

$$
\mathbf{V} = \mathbf{M} \cdot \mathbf{Z}
$$

解为：

$$
\mathbf{Z} = (\mathbf{M}^T \mathbf{W} \mathbf{M})^{-1} \mathbf{M}^T \mathbf{W} \mathbf{V}
$$

其中 $\mathbf{M}$ 为线表矩阵，$\mathbf{W}$ 为噪声权重矩阵。

### LSD 效果

**SNR 提升**：

$$
\text{SNR}_{\text{LSD}} \approx \sqrt{N_{\text{lines}}} \times \text{SNR}_{\text{single}}
$$

- 典型 G/K 型星：$N_{\text{lines}} \sim 1000$-5000 条可用谱线
- SNR 提升 **10-70 倍**
- 单个速度 bin 的噪声：$\sigma \sim 10^{-4}$-$10^{-5}$

**输出产品**：速度空间的平均 Stokes 轮廓（$I$、$V$、噪声）

### 纵向磁场：重心法（COG）

获得高信噪比的 LSD 轮廓后，可直接通过**重心法**计算视线方向的盘积分纵向磁场 $\langle B_\ell \rangle$：

$$
\langle B_\ell \rangle = -2.14 \times 10^{-11} \frac{\displaystyle\int v\, V(v)\, \mathrm{d}v}{\lambda_0\, g_\text{eff}\, c \displaystyle\int \left[I_c - I(v)\right] \mathrm{d}v}
$$

其中：
- $v$ 为相对于谱线中心的速度（cm s$^{-1}$）
- $V(v)$、$I(v)$ 为 LSD 得到的平均轮廓
- $I_c$ 为连续谱强度（归一化后 $I_c = 1$）
- $g_\text{eff}$ 为等效朗德因子
- $\lambda_0$ 为参考波长，$c$ 为光速
- 常数 $2.14 \times 10^{-11}$ 为单位换算因子，结果以 Gauss 为单位

**物理意义**：$\langle B_\ell \rangle$ 为圆偏振轮廓的**一阶矩**，代表可见半球视线方向的净磁通量。

#### 典型参数取值

| 谱线 | $g_{\text{eff}}$ | $\lambda_0$ (nm) |
|------|-----------------|-----------------|
| LSD 光球轮廓 | 1.24 | 650 |
| H$\alpha$（色球） | 1.0 | 656.28 |
| Ca IRT（中色球） | 0.968 | 858.56 |

### LSD 的局限性

| 假设/限制 | 影响 |
|----------|------|
| 弱场近似 | 强磁场（> 数百 G）或强谱线分裂时误差增大 |
| 自相似假设 | 实际不同谱线轮廓形态略有差异，引入系统误差 |
| 线表完备性 | 缺失谱线或错误的 $g_{\text{eff}}$ 引入误差 |
| 速度分辨率 | LSD 轮廓分辨率受限于最大可用谱线数 |

---

## 塞曼多普勒成像（ZDI）

### 问题提出

LSD 给出每个旋转相位下的平均 Stokes V 轮廓（盘积分信号），但如何恢复恒星表面的**二维磁场分布**？

**关键物理**：恒星自转使不同经度的磁场区域对应不同的多普勒速度。因此 Stokes V 轮廓的速度结构携带空间信息。

### 磁场参数化：球谐展开

将表面磁场分解为**极向（poloidal）** 和 **环向（toroidal）** 分量，用球谐函数展开：

$$
B_r(\theta,\phi) = \sum_{\ell=1}^{L} \sum_{m=0}^{\ell} \text{Re}\!\left[\alpha_{\ell m}\, Y_{\ell m}(\theta,\phi)\right]
$$

$$
B_\theta(\theta,\phi) = -\sum_{\ell=1}^{L} \sum_{m=0}^{\ell} \text{Re}\!\left[\beta_{\ell m}\, Z_{\ell m}(\theta,\phi) + \gamma_{\ell m}\, X_{\ell m}(\theta,\phi)\right]
$$

$$
B_\phi(\theta,\phi) = -\sum_{\ell=1}^{L} \sum_{m=0}^{\ell} \text{Re}\!\left[\beta_{\ell m}\, X_{\ell m}(\theta,\phi) - \gamma_{\ell m}\, Z_{\ell m}(\theta,\phi)\right]
$$

其中辅助函数为：

$$
Y_{\ell m}(\theta,\phi) = c_{\ell m}\, P_{\ell m}(\cos\theta)\, e^{im\phi}$$

$$
X_{\ell m}(\theta,\phi) = \frac{c_{\ell m}}{\ell+1}\frac{im}{\sin\theta}\, P_{\ell m}(\cos\theta)\, e^{im\phi}$$

$$
Z_{\ell m}(\theta,\phi) = \frac{c_{\ell m}}{\ell+1}\frac{\partial P_{\ell m}(\cos\theta)}{\partial\theta}\, e^{im\phi}$$

**展开系数意义**：
- $\alpha_{\ell m}$：极向场**径向**分量系数
- $\beta_{\ell m}$：极向场**子午**分量系数
- $\gamma_{\ell m}$：环向场分量系数
- $L$（最大阶次 $\ell_{\max}$）通常取 10-15，决定重建的最小空间尺度

### 反演问题

#### 正向模型

给定一组球谐系数，通过辐射转移理论在弱场近似下计算各旋转相位的 Stokes V 轮廓。

#### 优化目标

通过最小化目标函数进行反演：

$$
\min_{\alpha,\beta,\gamma} \left[ \chi^2(\text{obs vs model}) - \lambda_\text{reg} \cdot S(\alpha,\beta,\gamma) \right]
$$

其中：
- $\chi^2$ 为拟合优度
- $S$ 为**信息熵**（表示解的平滑度）
- $\lambda_\text{reg} > 0$ 为正则化参数

### 最大熵正则化（MEM）

ZDI 中最常用的正则化方案是**最大熵法**（Maximum Entropy Method）。

#### 核心思想

在所有能与观测数据相容（$\chi^2_r \approx 1$）的磁场分布中，选取**信息熵最大**的解——即最"平滑"、最接近先验默认值的解，避免在数据不足以约束的区域引入虚假结构。

#### 熵函数定义

对于离散化的表面像素集合：

$$
S = -\sum_{k} \left( f_k \ln\frac{f_k}{q_k} - f_k + q_k \right)
$$

其中：
- $f_k$ 为第 $k$ 个像素的（正定）磁通量密度
- $q_k$ 为对应的先验默认值（通常取全局均匀背景）

#### MEM 优化等价问题

$$
\max_{\alpha,\beta,\gamma}\; S(\alpha,\beta,\gamma) \quad \text{s.t.} \quad \chi^2(\text{obs vs model}) = 1
$$

在 Lagrange 乘子意义下等价于最小化目标函数。

#### MEM 的关键性质

| 性质 | 说明 |
|------|------|
| 偏向平滑解 | 倾向于将磁能集中于大尺度（低 $\ell$）球谐分量 |
| 正则化自适应 | 正则化参数 $\lambda_\text{reg}$ 随迭代自动调整 |
| 低估磁场强度 | 平滑偏置导致强磁极被"抹平"，重建的 $\langle\|B\|\rangle$ 通常偏低 |
| 解的唯一性 | 在给定先验下，满足 $\chi^2_r=1$ 的最大熵解唯一 |

### ZDI 迭代流程

```
初始猜测（零磁场）
    │
    ▼
计算正向模型 Stokes V（所有旋转相位）
    │
    ▼
与观测比较，计算 χ²
    │
    ▼
梯度下降更新球谐系数（最大熵 / L-BFGS）
    │
    ▼
收敛判据：χ²_r ≈ 1 或 δS/δα → 0
    │
    ▼
输出：表面磁场分布图（Br、Bθ、Bφ）
```

### ZDI 的输出产品

- **表面磁场分布图**（经纬度坐标）：$B_r(\theta, \phi)$、$B_\theta(\theta, \phi)$、$B_\phi(\theta, \phi)$
- **磁场能量谱**：各球谐阶次 $\ell$ 的磁能占比（偶极/四极/更高阶成分分解）
- **平均磁场强度**：$\langle|B|\rangle$，偶极分量、轴对称分量等统计

### ZDI 的局限性

| 局限性 | 说明 |
|--------|------|
| 空间分辨率 | 只能重建 $\ell \lesssim v\sin i / (P_{\text{rot}} \cdot \delta v)$ 的结构 |
| 南北简并 | $i < 90°$ 时南北半球磁场的某些分量无法区分 |
| 正则化偏差 | 最大熵倾向于压制小尺度结构，低估磁场总强度 |
| 仅大尺度场 | 小尺度随机场在 Stokes V 中完全抵消，不可见 |
| 参数依赖性 | $i$、$v\sin i$、$P_{\text{rot}}$、$T_{\text{eff}}$、$\log g$ 的误差直接影响重建 |

---

## 完整工作流

从原始偏振光谱到最终磁场分布的完整流程：

```
原始偏振光谱（多历元 Stokes I、V）
        │
        ▼ 数据预处理（Pipeline：Libre-ESpRIT / REDUCE）
Stokes I/V/N 归一化谱（每个旋转相位）
        │
        ▼ LSD 反卷积（VALD 线表 + 反卷积算法）
多相位 LSD 平均轮廓（速度空间 Stokes I/V）
        │
        ▼ 确定恒星参数
恒星参数：v sin i、磁极倾角 i、自转周期 Prot、
有效温度 Teff、重力加速度 log g
        │
        ▼ ZDI 球谐反演（最大熵正则化）
表面磁场分布图（Br, Bθ, Bφ）
        │
        ▼ 磁场分析与科学应用
磁场拓扑 / 磁能谱 / 长期演化 / 磁活动周期
```

### 关键输入参数

ZDI 反演的精度依赖于恒星参数的准确性：

| 参数 | 获取方法 | 影响 |
|------|---------|------|
| $v\sin i$ | Stokes I 谱线宽度拟合 | 直接影响空间分辨率 |
| 磁极倾角 $i$ | Stokes V 轮廓时间变化拟合 | 决定南北半球的分解能力 |
| 自转周期 $P_{\text{rot}}$ | Stokes V 轮廓周期性 | 影响相位采样和多普勒映射 |
| $T_{\text{eff}}$、$\log g$ | 光谱分类或其他观测 | 影响线表选择和 Stokes 轮廓理论计算 |

---

## 软件工具

### 专用工具

| 软件 | 功能 | 特点 |
|------|------|------|
| **LSDpy** / **iLSD** | LSD 反卷积 | 部分开源，广泛使用 |
| **ZDIpy** | ZDI 球谐反演 | 本项目的核心工具，Python 实现 |
| **INVERS** | ZDI 参考实现 | Piskunov/Kochukhov 组，学术受限 |
| **starpolator** | 偏振谱正向模型 | 学术共享 |

### 辅助资源

| 资源 | 用途 |
|------|------|
| **VALD3** | 原子线表数据库，提供 $\lambda$、$g_{\text{eff}}$、$d$ 等参数 |
| **Libre-ESpRIT** / **REDUCE** | 数据预处理和 Stokes 轮廓提取 |

---

## 参考文献与进一步阅读

- Donati, J.-F., et al. (2006). Stellar Magnetism. *Annual Review of Astronomy and Astrophysics*, 44, 269–302.
- Rees, D. E., & Semel, M. D. (1979). Line formation in magnetic fields. *Astronomy & Astrophysics*, 74, 1–5.
- Skilling, J., & Bryan, R. K. (1984). Maximum entropy image reconstruction. *Monthly Notices of the Royal Astronomical Society*, 211, 111–124.
- Kochukhov, O. (2021). Magnetic field topology of stars. *Annual Review of Astronomy and Astrophysics*, 59, 47–90.

---

## ZDIpy 中的实现

本项目的 ZDIpy 工具基于上述理论框架，提供以下功能：

- **LSD 轮廓管理**：读取和处理多相位 LSD 数据
- **球谐反演**：使用最大熵法进行 ZDI 反演
- **磁场可视化**：生成磁场分布图、能量谱、拓扑图
- **参数优化**：恒星参数和 ZDI 超参数的调整
- **WebUI 界面**：交互式反演和结果查看

更多技术细节，请参阅项目的其他文档：
- [stellar_geometry.md](stellar_geometry.md) - 恒星几何与坐标变换
- [magnetic_map.md](magnetic_map.md) - 磁场分布图的生成与解释
- [brightness_map.md](brightness_map.md) - 亮度图与磁场的联合反演
- [MEM.md](MEM.md) - 最大熵正则化的实现细节
