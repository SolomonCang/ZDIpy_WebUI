# 逐历元温度变化算法说明（Epoch Temperature Variation, dT）

## 1. 设计目标

本文档说明 ZDIpy_WebUI 中新增的逐历元温度变化算法。该算法用于在不同观测历元（相位）上，对恒星表面亮度图进行相位依赖修正，以近似描述温度活动区（例如冷斑、热点）随时间演化的效应。

核心思想：

- 温度变化不直接拟合为温度场，而是通过亮度参数变化进行映射。
- 采用关系：dT 与 dCq 反相关。即温度降低通常对应亮度降低。
- 在每个观测历元中，对基础亮度图 Cq 应用相位包络函数，得到历元亮度 Cq_epoch。

当前实现为前向模型扩展：

- 已支持逐历元亮度调制并进入合成谱计算。
- 未将 dT 本身作为独立 MEM 反演参数。

---

## 2. 数学模型

设：

- Cq 为基础亮度图（每个表面网格点一个值）。
- phase 为当前观测历元相位。
- dphi 为相位偏移。
- t 为特征寿命（单位：rotation cycles）。
- kk 为时间尺度调制参数。

定义归一化相位偏移：

$$
 s = \frac{t}{2\sqrt{\ln 2}}, \quad dd = \frac{phase - dphi}{s}
$$

算法返回：

- Cq_epoch：该历元下修正后的亮度图。
- dCq_dCq0：局部 Jacobian（对基础亮度的偏导）。

### 2.1 模式 1：Gaussian

$$
 f(dd) = e^{-dd^2}
$$
$$
 Cq_{epoch} = 1 - (1 - Cq) f(dd)
$$
$$
 \frac{\partial Cq_{epoch}}{\partial Cq} = f(dd)
$$

特点：平滑、对称，适合描述缓慢出现和衰减的活动区。

### 2.2 模式 2：Cosine squared

$$
 f(dd) = \cos^2\left(\frac{\pi dd}{2}\right), \quad |dd| \le 1
$$
$$
 f(dd)=0, \quad |dd|>1
$$
$$
 Cq_{epoch} = 1 - (1 - Cq) f(dd)
$$

特点：有界支持区间，窗口外无调制。

### 2.3 模式 3：Tanh

实现使用平滑双曲正切形式，避免数值发散并处理 Cq 约等于 1 时的除零风险。

该模式在代码中加入了数值保护：

- 当中间变量接近 0 时使用安全下限替代，避免 NaN。
- Jacobian 采用保守近似（常数 1）以保证稳定性。

特点：可用于更平滑的非线性过渡，同时保持数值鲁棒。

### 2.4 模式 4：Exponential

$$
 ph0 = kk\ln(kk)
$$
$$
 Cq_{epoch} = Cq \cdot e^{dd\,ph0}
$$
$$
 \frac{\partial Cq_{epoch}}{\partial Cq} = e^{dd\,ph0}
$$

特点：变化快，适合快速增强或衰减的活动区。

---

## 3. 代码集成位置

### 3.1 亮度调制函数

新增函数：epoch_brightness_scale

作用：输入基础亮度图和当前相位，输出历元亮度与 Jacobian。

### 3.2 拟合循环注入点

在主拟合循环每个观测历元中：

1. 读取当前 phase。
2. 若启用 epoch variation，则计算 Cq_epoch。
3. 用临时亮度图副本参与该历元的谱线积分与导数计算。
4. 不修改全局基础亮度图，避免跨历元污染。

这保证了每个历元使用自己的温度修正，同时不破坏原有 MEM 状态向量流程。

---

## 4. 配置结构

在 brightness 节点下新增：

- epoch_variation.enabled
- epoch_variation.var_mode
- epoch_variation.lifetime_cycles
- epoch_variation.kk
- epoch_variation.dphi

含义：

- enabled：是否启用逐历元修正。
- var_mode：1/2/3/4，对应四种调制模型。
- lifetime_cycles：特征寿命（旋转周数）。
- kk：时间尺度参数。
- dphi：相位偏移。

默认行为：

- enabled 为 false 时，算法对现有流程完全透明。

---

## 5. Web 配置面板兼容

前端配置使用平面字段，后端配置文件使用嵌套字段。

已实现双向映射：

- GET 配置时：后端将嵌套结构展开为前端可用的平面字段。
- PUT 配置时：后端将平面字段回写为嵌套结构并持久化。

这样既保证了前端表单系统兼容性，也保证了配置文件语义清晰。

---

## 6. 数值稳定性与边界行为

### 6.1 稳定性处理

- 模式 3 中已处理除零风险，避免 NaN。
- 输入亮度图不被原地改写，避免副作用。
- 所有模式在常见参数范围内均可返回有限值。

### 6.2 边界说明

- 模式 4 可能在较大参数组合下导致亮度放大较快，应结合 max_bright 约束与物理先验使用。
- 模式 2 在窗口外直接归零调制，可能产生较强分段行为。

---

## 7. 对反演流程的影响

### 7.1 当前阶段（已实现）

- 作为前向模型扩展，改变每个历元的合成谱生成。
- 不改变 MEM 状态向量维度。
- 不新增 dT 参数的反演自由度。

### 7.2 下一阶段（可扩展）

若要将 dT 参数并入反演，可在响应矩阵中引入链式法则：

$$
 \frac{\partial I}{\partial dT}
 = \frac{\partial I}{\partial Cq}
 \cdot \frac{\partial Cq}{\partial dT}
$$

其中第一项已有亮度导数基础，第二项可由选定的温度-亮度映射关系给出。

---

## 8. 典型使用建议

1. 初次启用建议使用模式 1（Gaussian），参数取默认值。
2. 若希望时窗更明确可用模式 2。
3. 若出现非线性过强或数值敏感问题，优先尝试模式 3。
4. 模式 4 建议用于验证快速演化假设，避免直接用于最终物理解读。

建议实验顺序：

1. enabled=false 跑一组基线结果。
2. enabled=true 且 var_mode=1 跑对照。
3. 比较残差、chi2 和重建图稳定性。
4. 再逐步扫描 var_mode 与 lifetime_cycles。

---

## 9. 验证结论摘要

本次扩展验证覆盖了：

- 配置解析正确性。
- 四种调制模型均可计算。
- 拟合循环中逐历元调用逻辑正确。
- 原始亮度图不会被历元修正污染。
- 模式 3 的 NaN 问题已修复。

结论：逐历元温度变化算法已可用于生产流程中的前向建模与对比实验。

---

## 10. 相关模块

- core/brightnessGeom.py
- core/fitting.py
- config_loader.py
- api/routes/config.py
- frontend/js/config.js
- config.json
