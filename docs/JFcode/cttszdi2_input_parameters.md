# CTTSzdi2 输入参数总表

本文档整理 `/Users/tianqi/Documents/Codes_collection/ZDI_and/CTTSzdi2` 中
`zdipot` 程序真正会读取或依赖的输入信息，重点回答两个问题：

1. 用户在运行时需要提供哪些参数？
2. 这些参数在源码中的入口、提示语、物理含义和样例值分别是什么？

说明：

- 这里的“样例值”优先取自 `zdipot23.in`、`zdipot23.out` 和 `TWA12-23.ss1`。
- 本文把参数分成四类：命令行开关、主交互/批处理输入、观测文件字段、源码中硬编码但实际影响反演的谱线模型常量。
- 最后一类严格说不是“外部输入参数”，但它们会实质改变程序行为，因此单独列出。

---

## 1. 命令行开关

| 源码位置 | 提示语 | 物理含义 | 样例值 |
|----------|--------|----------|--------|
| `utilities.c:parse_args()` | 无显式提示，命令行直接传入 | `-b`：批处理模式，从标准输入连续读取预设答案 | `./zdipot -b < zdipot23.in` |
| `utilities.c:parse_args()` | 无显式提示，命令行直接传入 | `-m`：把速度单位视为 `m/s` | 未使用 |
| `utilities.c:parse_args()` | 无显式提示，命令行直接传入 | `-k`：把速度单位视为 `km/s`，内部乘以 1000 转成 `m/s` | 默认行为 |

---

## 2. 主交互/批处理输入参数

下表按 `potential.c` 主流程中的实际问答顺序排列。

| 源码位置 | 提示语 | 物理含义 | 样例值 |
|----------|--------|----------|--------|
| `utilities.c:init_image()` | `Do you want to make a magnetic image (y/n) :` | 是否进入磁 ZDI 路径；若为 `n`，程序会认为当前工具不适用 | `y` |
| `utilities.c:init_image()` | `Do you want to read an old image file (y/n) :` | 是否从旧图像文件读取初始表面图与几何参数 | `n` |
| `utilities.c:init_image()` | `Name of file to read ...`（由 `read_bright()` 触发） | 旧图像文件路径；若为 `y`，会跳过手动输入 `ngrid/incl/vsin(i)/vrad` | 样例未使用 |
| `utilities.c:init_image()` | `Input ngrid, incl, vsin(i) and vrad :` | 网格单元数、恒星倾角、投影自转速度、径向速度 | `500 80 17.9 0.0` |
| `utilities.c:get_imsw2()` | `Select type of image reconstruction (any combination of QMB) :` | 选择拟合的图像块：`Q`=`Cq` 光球亮度，`M`=`Cm` 第二图像量，`B`=磁场；`V/H/L/E` 为 `Cm` 的时间变化模式；`C` 为特殊亮度模式 | `QB` |
| `utilities.c:get_imsw2()` | `Relative weights of B vs Cq/Cm :` | 当同时拟合亮度块和磁场块时，设置磁场块相对亮度块的 MEM 权重 | `1.0` |
| `utilities.c:get_imsw2()` | `Relative weights of Cm vs Cq :` | 当同时拟合 `Cq` 与 `Cm` 且不拟合磁场时，设置 `Cm` 相对 `Cq` 的权重 | 样例未使用 |
| `utilities.c:get_imsw2()` | `Relative weights of Cm vs Cq and B vs Cq :` | 当 `Cq`、`Cm`、`B` 都参与拟合时，同时设置两类相对权重 | 样例未使用 |
| `potential.c:71-92` | `Type of reconstruction:` 后接数值输入 | 选择磁场球谐约束类型：0=无约束，1=极向+环向，2=广义势场，3=广义力自由，4=势场，5=线性力自由，负值为几种历史变体 | `-4` |
| `potential.c:94-96` | `Order of spherical harmonics expansion :` | 球谐展开最高阶 `lmax` | `10` |
| `potential.c:98-103` | `LFF alpha (default = 1.0) :` | 仅对 `potential=3/5/-2` 生效的线性力自由参数 `alpha` | 样例未使用 |
| `medium_ttssp.c:image_defaults()` | `Input default spot brightness :` | 亮度图 `Cq` 的默认值；也是部分模式下 `Cm` 的初始参考值 | `1.0` |
| `medium_ttssp.c:image_defaults()` | `What is the default ... :` | 第二图像量默认值；只有 `Cm_used` 为真时才出现 | 样例未使用 |
| `potential.c:126-128` | `What is a typical value for the magnetic field strength :` | 典型磁场尺度 `Bscale`，用于初始化和 MEM 常数标度 | `1000` |
| `grid.c:init_diffrot()` | `Beta and gamma for differential rotation (0 for default) :` | 差分自转两个系数 `beta`、`gamma` | `0 0` |
| `potential.c:177-188` | `Name of file to read spectral data from :` | 观测输入文件路径，通常是 `.ss1` 风格的多观测文本文件 | `TWA12-23.ss1` |
| `potential.c:189-206` | `Are profiles affected by moon pollution (y/n) :` | 是否启用月光污染误差修正 | `n` |
| `potential.c:200-214` | `File containing info on moon pollution -` | 月光污染信息文件；用于重设误差条 | 样例未使用 |
| `potential.c:get_lines_to_fit()` | `Do you want to use the Mean observations (y/n) :` | 是否拟合名为 `Mean` 的观测组 | `y` |
| `potential.c:get_lines_to_fit()` | `Do you want to use the Accr observations (y/n) :` | 是否拟合名为 `Accr` 的观测组 | 样例未使用 |
| `potential.c:get_lines_to_fit()` | `Do you want to use the phot observations (y/n) :` | 是否拟合光度观测组 | `y` |
| `utilities.c:get_datsw()` | `Select Stokes parameters to use (any combination of ...) :` | 从可用 Stokes 分量中选择实际参与拟合的子集 | `IV` |
| `potential.c:224-226` | 隐式 `scanf`，随后打印 `Reference angle for QU (rad):` | 当同时拟合 `Q` 和 `U` 时输入参考偏振角 | 样例未使用 |
| `potential.c:get_photometry_info()` | `Unspotted V magnitude :` | 无斑点参考星等，用于把 photometric 观测转换到拟合量 | `1.0` |
| `potential.c:get_photometry_info()` | `Weight on photometry :` | 光度数据相对光谱数据的权重 | `1.0` |
| `potential.c:231-233` | `What is the largest phase smearing allowed :` | 允许的最大相位展宽，用于确定每条观测的相位采样数 `nsample` | `1.0` |
| `potential.c:235-237` | `What is the oversampling factor :` | 每条谱线的波长过采样因子，用于正演和卷积 | `1` |
| `potential.c:259-263` | `Input the number of iterations per response matrix calculation :` | 仅 quasi-linear 模式下使用；控制多久重算一次响应矩阵 | 样例未使用 |
| `potential.c:274-281` | `Do you want to weight the spherical harmonics (y/n) :` | 是否按球谐阶次对磁系数施加额外权重 | `y` |
| `potential.c:287-294` | `Enter value for L_fac (0 for default) :` | MEM 步长/正则控制常数 `L_fac` | `0.1` |
| `potential.c:299-303` | `Input values for Caim and Maxit :` | 目标拟合值 `Caim` 与最大迭代次数 `Maxit` | `1.0 61` |
| `potential.c:603-610` | `Name of file to save brightness data to :` | 输出亮度/图像数据文件名 | `TWA12-23.m1` |
| `potential.c:611-620` | `Do you want to save the synthetic spectra (y/n) :` | 是否输出合成光谱 | `y` |
| `potential.c:614-620` | `Name of file to save spectral data to :` | 合成光谱输出文件名 | `TWA12-23.s1` |
| `potential.c:621-631` | `Do you want to save the SHaDe coefficients (y/n) :` | 是否输出球谐系数 | `y` |
| `potential.c:624-631` | `Name of file to save mode coefficients to :` | 球谐系数输出文件名 | `TWA12-23.c1` |

### 2.1 `QMB` 模式字母的内部含义

`get_imsw2()` 不只是简单的 `Q/M/B` 三选一，它还把若干字母映射到 `fitcm` 和 `var` 的特殊模式：

| 源码位置 | 提示语 | 物理含义 | 样例值 |
|----------|--------|----------|--------|
| `utilities.c:get_imsw2()` | 同 `QMB` 入口 | `Q`：拟合 `Cq`，即光球亮度图 | `Q` |
| `utilities.c:get_imsw2()` | 同 `QMB` 入口 | `M`：拟合 `Cm`，即第二图像量/吸积或填充图 | 样例未使用 |
| `utilities.c:get_imsw2()` | 同 `QMB` 入口 | `B`：拟合磁场球谐系数 | `B` |
| `utilities.c:get_imsw2()` + `medium_ttssp.c:calcvar()` | 同 `QMB` 入口 | `V`：启用 `Cm` 的高斯型时间变化，内部 `var=1` | 样例未使用 |
| `utilities.c:get_imsw2()` + `medium_ttssp.c:calcvar()` | 同 `QMB` 入口 | `H`：启用余弦平方时间变化，内部 `var=2` | 样例未使用 |
| `utilities.c:get_imsw2()` + `medium_ttssp.c:calcvar()` | 同 `QMB` 入口 | `L`：启用线性时间变化，内部 `var=3` | 样例未使用 |
| `utilities.c:get_imsw2()` + `medium_ttssp.c:calcvar()` | 同 `QMB` 入口 | `E`：启用指数时间变化，内部 `var=4` | 样例未使用 |
| `utilities.c:get_imsw2()` | 同 `QMB` 入口 | `C`：特殊模式，开启 `Cq`，并让 `Cm` 走固定负值路径 `fitcm=-1` | 样例未使用 |

---

## 3. 观测输入文件字段

`zdipot` 通过 `file.c:read_spec()` 读取一个文本观测文件。
`TWA12-23.ss1` 的前两条记录显示，这个文件不仅包含谱线数组，还包含每条观测的元数据。

### 3.1 文件级字段

| 源码位置 | 提示语 | 物理含义 | 样例值 |
|----------|--------|----------|--------|
| `file.c:read_spec()` | 无提示，文件第 1 行 | 注释/标题行 | `LSD TWA12 11.7 - 08jun23/...` |
| `file.c:read_spec()` | 无提示，文件第 2 行 | 观测条目数 `nobs` | `81` |

### 3.2 每条观测的元数据字段

| 源码位置 | 提示语 | 物理含义 | 样例值 |
|----------|--------|----------|--------|
| `file.c:read_spec()` | 无提示，观测首行第 1 列 | 谱线名称 `line_name`，用于匹配内置线模型 | `Mean` |
| `file.c:read_spec()` | 无提示，观测首行第 2 列 | 名义线心波长 `line_centre` | `1.750000000e+03` |
| `file.c:read_spec()` | 无提示，观测首行第 3 列 | 是否是光度观测 `photometric` | `0` |
| `file.c:read_spec()` | 无提示，观测首行第 4 列 | `datsw[0]`，是否有 Stokes I | `1` |
| `file.c:read_spec()` | 无提示，观测首行第 5 列 | `datsw[1]`，是否有 Stokes Q | `0` |
| `file.c:read_spec()` | 无提示，观测首行第 6 列 | `datsw[2]`，是否有 Stokes U | `0` |
| `file.c:read_spec()` | 无提示，观测首行第 7 列 | `datsw[3]`，是否有 Stokes V | `0` |
| `file.c:read_spec()` | 无提示，观测次行第 1 列 | 光谱点数 `nspec` | `40` |
| `file.c:read_spec()` | 无提示，观测次行第 2 列 | 曝光开始对应的旋转相位 `phasec` | `0.010003` |
| `file.c:read_spec()` | 无提示，观测次行第 3 列 | 曝光结束对应的旋转相位 `phasef` | `0.010003` |
| `file.c:read_spec()` | 无提示，观测次行第 4 列 | 波长网格起点 `laml` | `1.749768256e+03` |
| `file.c:read_spec()` | 无提示，观测次行第 5 列 | 波长网格终点 `lamu` | `1.750223571e+03` |
| `file.c:read_spec()` | 无提示，观测次行第 6 列 | 信噪比 `SN` | `1.118708e+03` |
| `file.c:read_spec()` | 无提示，观测次行第 7 列 | 连续谱基线 `continuum` | `1.003000e+00` |
| `file.c:read_spec()` | 无提示，观测次行第 8 列 | 仪器展宽 `FWHM` | `2.692308e-02` |
| `file.c:read_spec()` | 无提示，观测次行第 9 列 | 偏振器角度 `polariser`，读入后转为弧度 | `0.000000` |

### 3.3 每条观测的数组数据

| 源码位置 | 提示语 | 物理含义 | 样例值 |
|----------|--------|----------|--------|
| `file.c:read_spec()` | 无提示，后续浮点数组 | 若 `datsw[0]=1`，读取该观测的 Stokes I 数组 | `1.00212e+00 ... 1.00187e+00` |
| `file.c:read_spec()` | 无提示，后续浮点数组 | 若 `datsw[1]=1`，读取该观测的 Stokes Q 数组 | 首条记录无 |
| `file.c:read_spec()` | 无提示，后续浮点数组 | 若 `datsw[2]=1`，读取该观测的 Stokes U 数组 | 首条记录无 |
| `file.c:read_spec()` | 无提示，后续浮点数组 | 若 `datsw[3]=1`，读取该观测的 Stokes V 数组 | 首条记录无 |

---

## 4. 程序内置、但会实质影响拟合的硬编码谱线参数

这一组参数并不是运行时输入，而是直接写死在 `medium_ttssp.c` 里。
如果只是想复现实验结果，必须同时记录这些默认值；如果想把 CTTSzdi2 迁移到新框架，
这些量应该显式转成配置字段。

| 源码位置 | 提示语 | 物理含义 | 样例值 |
|----------|--------|----------|--------|
| `medium_ttssp.c:line_names[]` | 无外部提示 | 可用谱线/观测组名称 | `Mean`, `Accr`, `phot` |
| `medium_ttssp.c:l0[]` | 无外部提示 | 线心波长 | `1750.0`, `850.0`, `1750.0` |
| `medium_ttssp.c:dl[]` | 无外部提示 | Doppler 宽度 | `120.0E-4`, `200.0E-4` |
| `medium_ttssp.c:eps[]` | 无外部提示 | 临边参数 `epsilon` | `0.30`, `0.0`, `0.30` |
| `medium_ttssp.c:bb[]` | 无外部提示 | Milne-Eddington 参数 `bb` | `3.0`, `3.0` |
| `medium_ttssp.c:eta[]` | 无外部提示 | 线心不透明度比 `eta` | `2.0`, `2.0` |
| `medium_ttssp.c:aa[]` | 无外部提示 | 阻尼常数 `a` | `0.95`, `0.02` |
| `medium_ttssp.c:g[]` | 无外部提示 | Landé g 因子 | `1.2`, `1.0` |
| `medium_ttssp.c:zeta_rt[]` | 无外部提示 | 径向-切向宏湍流参数 | `0.0`, `0.0` |
| `medium_ttssp.c:fff[]` | 无外部提示 | filling factor | `0.20`, `0.20` |
| `medium_ttssp.c:pmf[]` | 无外部提示 | 轮廓乘子 `pmf`，决定是否走 CTTS 相关特殊路径 | `0.50`, `-1.1` |

---

## 5. 运行流程摘要

把这些输入放回程序控制流里，可以把 `zdipot` 理解为下面四层：

1. **启动层**：命令行只决定批处理与速度单位；
2. **反演设置层**：通过标准输入依次选择图像块、磁场几何、迭代控制、输出文件；
3. **观测层**：通过 `.ss1` 文件提供每条观测的线名、相位、波长范围、S/N、Stokes 数组；
4. **模型层**：通过 `medium_ttssp.c` 中的硬编码常量定义每种谱线的物理参数。

因此，CTTSzdi2 真正“考虑”的参数并不都来自 `.in` 文件：

- `.in` 文件负责回答主程序的交互问题；
- `.ss1` 文件负责提供观测元数据和数据数组；
- 源码里的静态数组负责给出谱线物理模型常数。

如果后续要把这套旧程序完整迁移到 `ZDIpy_WebUI`，最值得优先显式化的就是：

1. `QMB` 模式与 `var` 的对应关系；
2. `potential` 编号与磁场几何类型的映射；
3. `medium_ttssp.c` 中全部硬编码谱线参数；
4. `.ss1` 观测文件的元数据结构。

---

## 6. CTTSzdi2 参数到 ZDIpy_WebUI `config.json` 的映射表

说明：

- 这里的“当前字段”使用 `config.json` 的点路径写法，例如 `star.inclination_deg`。
- “直接”表示语义基本一致；“换算/折叠”表示需要改单位、改结构或合并多个旧参数；“无直接字段”表示当前公开 schema 还不能严格复现旧行为。
- 对旧程序的历史模式，本文优先给出“当前最接近的 config 表达”，并明确指出缺失的自由度。

### 6.1 恒星、磁场与反演控制

| CTTSzdi2 参数 | 当前 `config.json` 字段 | 映射状态 | 说明 |
|---------------|-------------------------|----------|------|
| `incl` | `star.inclination_deg` | 直接 | 倾角单位都可直接用度。 |
| `vsin(i)` | `star.vsini_kms` | 直接 | 都是 `km/s`。 |
| `vrad` | `observations.files[].vel_center_kms` | 换算/折叠 | 旧程序是一个全局径向速度；新 schema 是逐观测的速度中心。若所有历元共用同一系统速度，可把同一个值写入每条观测。 |
| `ngrid` | `grid.nRings` | 无直接字段 | 旧程序直接输入表面单元总数；当前 schema 用纬向环数生成网格，不存在严格的一一对应，只能选一个接近的角分辨率。 |
| `potential = 0 / -3 / -4` | `magnetic.geometry_type = "Full"` | 换算/折叠 | 这几类旧模式都保留三组磁系数自由度；当前最接近的是 `Full`。其中 `-3/-4` 的历史特化关系不会被保留。 |
| `potential = 1` | `magnetic.geometry_type = "PotTor"` | 直接 | 旧程序的 “poloidal field plus toroidal field” 与当前 `PotTor` 的拓扑最接近，都是把环向分量保留、并对 poloidal 部分施加势场关系。 |
| `potential = 4` | `magnetic.geometry_type = "Potential"` | 直接 | 都对应纯势场约束。 |
| `potential = 2 / -1` | `magnetic.geometry_type = "Poloidal"` | 换算/折叠 | 当前 schema 没有“generalised potential”专用类型；若只保留“无环向场”的拓扑，可近似落到 `Poloidal`，但旧代码中的广义关系不会完整保留。 |
| `potential = 3 / 5 / -2` 与 `alpha` | 无直接字段 | 无直接字段 | 这些模式依赖 force-free/LFF 的额外 `alpha` 参数。当前 `config.json` 只有 `geometry_type`，没有对应的 `alpha` 配置，因此不能精确等价。 |
| `lmax` | `magnetic.l_max` | 直接 | 球谐最高阶一一对应。 |
| `Bscale` | `magnetic.default_bent` | 直接 | 都是磁场初始/默认尺度，单位为高斯量级。 |
| 读旧磁图文件开关 | `magnetic.init_from_file` + `magnetic.init_file` | 直接 | 对磁场系数初始化是直接映射。 |
| 读旧亮度图文件开关 | `brightness.init_from_file` + `brightness.init_file` | 直接 | 对亮度图初始化是直接映射。 |
| `Q` | `brightness.fit_brightness = 1` | 直接 | 拟合亮度图。 |
| `B` | `magnetic.fit_magnetic = 1` | 直接 | 拟合磁图。 |
| `QB` | `brightness.fit_brightness = 1` 且 `magnetic.fit_magnetic = 1` | 直接 | 同时拟合亮度和磁场。 |
| `M` / `C` | 无直接字段 | 无直接字段 | 当前公开 schema 只有一个亮度图块，没有旧代码的第二图像量 `Cm` 或 `fitcm = -1` 特殊路径。 |
| `V/H/L/E` 时间变化模式 | `brightness.epoch_variation.enabled` + `brightness.epoch_variation.var_mode` | 换算/折叠 | `var_mode = 1/2/3/4` 可分别对应高斯、余弦平方、线性、指数变化；但当前实现作用在单亮度图的历元变化上，不等价于旧 `Cm` 通道。 |
| `default spot brightness` | `brightness.default_bright` | 直接 | 默认亮度一一对应。 |
| `What is the default ...` (`Cm` 默认值) | 无直接字段 | 无直接字段 | 当前 schema 没有单独的 `Cm` 默认值。 |
| `Relative weights of B vs Cq/Cm`、`Cm vs Cq` | 无直接字段 | 无直接字段 | 当前公开 schema 没有按图像块分别设置 MEM 相对权重的通用字段。`brightness.entropy_scale`、`chi2_scale_I/V` 只能做近似调权，语义并不等价。 |
| `Beta and gamma for differential rotation` | `star.differential_rotation_rad_per_day` | 换算/折叠 | 旧程序用两个系数描述差分自转；当前 schema 只有单个 `dOmega` 风格参数，无法完整保留两参数形式。 |
| `Caim` | `inversion.target_form = "C"` + `inversion.target_value` | 直接 | 当目标形式选 `C` 时，`target_value` 就是旧的目标拟合值。 |
| `Maxit` | `inversion.num_iterations` | 直接 | 最大迭代数一一对应。 |
| `L_fac` | 无直接字段 | 无直接字段 | 当前公开 schema 没有暴露旧 MEM 步长常数。 |
| `weight the spherical harmonics` | 无直接字段 | 无直接字段 | 当前 schema 没有“按阶次额外加权球谐系数”的布尔开关。 |
| `iterations per response matrix calculation` | 无直接字段 | 无直接字段 | 当前 schema 没有 quasi-linear 重算响应矩阵频率这一控制量。 |

### 6.2 观测文件、相位与光度数据

| CTTSzdi2 参数 | 当前 `config.json` 字段 | 映射状态 | 说明 |
|---------------|-------------------------|----------|------|
| `.ss1` 单文件多观测输入 | `observations.files[]` | 换算/折叠 | 旧程序把所有历元装在一个文本文件里；当前 schema 倾向于“一条观测一个文件”，需要先拆分。 |
| `phasec` / `phasef` | `observations.jdate_ref` + `observations.files[].jdate` + `star.period_days` | 换算/折叠 | 当前 schema 不直接存旋转相位，而是存 Julian date。若要从旧相位恢复，可用 `jdate = jdate_ref + phase * period_days`；曝光起止相位宽度本身没有独立字段。 |
| `line_name` (`Mean` / `Accr` / `phot`) | 无直接字段 | 无直接字段 | 当前 schema 不把观测组名作为统一配置字段；谱线类型主要由输入文件和 `line_model.model_type` 决定。 |
| `line_centre` | `line_model.wavelength_nm` | 直接 | 当一组观测共享同一条谱线时可以直接映射。 |
| 每条观测的 `FWHM` | `instrument.spectral_resolution` | 换算/折叠 | 旧程序是逐观测输入展宽；当前 schema 只有全局分辨率。若 `FWHM` 基本恒定，可近似用 `R \approx \lambda / \mathrm{FWHM}` 转成一个统一值。 |
| `datsw[0..3]` Stokes 开关 | 无直接字段 | 无直接字段 | 当前配置里没有单独的 Stokes 选择掩码；可用分量取决于输入观测文件实际包含的列以及当前管线支持。 |
| `Select Stokes parameters to use` | 无直接字段 | 无直接字段 | 同上，当前 schema 没有再提供一次“从已有分量里二次选择”的显式开关。 |
| `Reference angle for QU` | 无直接字段 | 无直接字段 | 当前公开 schema 没有 QU 参考偏振角字段。 |
| `Do you want to use the phot observations` | `light_curve.fit_light_curve` + `light_curve.files[]` | 换算/折叠 | 当前光度数据走独立的 `light_curve` 分支，而不是夹在 `.ss1` 里。 |
| `Weight on photometry` | `light_curve.chi2_scale_lc` | 近似 | 都是给光度数据一个相对权重；新字段是对光变曲线残差的整体缩放。 |
| `Unspotted V magnitude` | 无直接字段 | 无直接字段 | 当前光变模块默认使用归一化通量，没有单独保存“无斑参考星等”。 |
| `largest phase smearing allowed` | 无直接字段 | 无直接字段 | 当前 schema 不暴露旧代码里由相位展宽决定 `nsample` 的控制量。 |
| `oversampling factor` | 无直接字段 | 无直接字段 | 当前 schema 不暴露旧的谱线过采样控制。 |
| `moon pollution` 开关与文件 | 无直接字段 | 无直接字段 | 当前公开配置没有月光污染误差修正输入。 |

### 6.3 谱线模型常量到 `line_model.*` 的映射

| CTTSzdi2 常量 | 当前 `config.json` 字段 | 映射状态 | 说明 |
|---------------|-------------------------|----------|------|
| `l0[]` | `line_model.wavelength_nm` | 直接 | 线心波长直接落到当前公共谱线字段。 |
| `dl[]` | `line_model.gauss_width_kms` | 换算 | 旧值是波长宽度，当前字段是速度宽度。换算关系为 `gauss_width_kms = c_kms * dl / l0`。 |
| `eps[]` | `line_model.limb_darkening` | 直接 | 临边昏暗系数一一对应。 |
| `bb[]` | `line_model.unno_beta` | 直接 | 都是 Unno/Milne-Eddington 源函数斜率参数。 |
| `eta[]` | `line_model.line_strength` | 近似 | 两者都控制线心强度/不透明度量级；当前字段更偏“统一强度参数”，对旧 `eta` 是最接近的承接位置。 |
| `aa[]` | `line_model.lorentz_width_fraction` | 直接 | 都是 Lorentz 相对宽度参数。 |
| `g[]` | `line_model.lande_g` | 直接 | Lande g 因子一一对应。 |
| `zeta_rt[]` | `line_model.unno_macro_turb_kms` | 换算 | 都表示额外展宽；若旧值本身就是速度尺度，可直接填入，否则需要先统一单位。 |
| `fff[]` | `line_model.unno_filling_factor_I` | 近似 | 旧填充因子最接近当前 Unno 的 Stokes I 填充因子。若还要单独控制 V，则需额外设置 `line_model.unno_filling_factor_V`。 |
| `pmf[]` | 无直接字段 | 无直接字段 | 旧代码用它切换 CTTS 特化轮廓路径；当前公开 schema 没有对应的 CTTS/双组分模式开关。 |
| `line_names[]` | 无直接字段 | 无直接字段 | 当前 schema 不维护“内置线族名称表”，而是由输入数据和所选 line model 决定。 |

### 6.4 输出文件名的映射

| CTTSzdi2 输出项 | 当前 `config.json` 字段 | 映射状态 | 说明 |
|-----------------|-------------------------|----------|------|
| 亮度图输出文件 | `output.bright_map_file` | 直接 | 亮度图文件名直接对应。 |
| SHaDe / 球谐系数输出文件 | `output.mag_coeff_file` | 直接 | 都保存磁场球谐系数。 |
| 合成谱输出文件 | `output.line_models_file` | 近似 | 当前最接近的是合成线轮廓输出文件；如果还要单独保存“用于拟合的观测子集”，则另有 `output.observed_used_file`。 |

### 6.5 一个最小迁移思路

如果你的目标是把一份典型 `zdipot23.in + .ss1 + medium_ttssp.c` 组合转成当前 `config.json`，最稳妥的顺序是：

1. 先把 `incl`、`vsin(i)`、`lmax`、`Bscale`、`Caim`、`Maxit` 直接填到 `star`、`magnetic`、`inversion`；
2. 再把 `potential` 先映射到最接近的 `magnetic.geometry_type`，并把所有“旧历史模式的额外自由度”记为暂时丢失；
3. 把 `medium_ttssp.c` 里的 `l0/dl/eps/bb/eta/aa/g/fff` 显式抄到 `line_model.*`；
4. 最后把 `.ss1` 拆成当前 `observations.files[]` 或 `light_curve.files[]` 认识的格式。

如果要做到“严格物理等价”，当前 schema 还至少缺三类字段：

1. `Cm` / 双图像块相关配置；
2. force-free / LFF 的 `alpha`；
3. CTTS 特化轮廓开关（旧 `pmf` 一类参数）。