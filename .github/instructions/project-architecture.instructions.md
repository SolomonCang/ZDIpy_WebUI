---
description: "ZDIpy_WebUI 完整项目架构总览，供所有贡献者快速理解模块职责、调用链与依赖关系。"
name: "Project Architecture Overview"
applyTo: ["**"]
---

# ZDIpy_WebUI 项目架构总览

> 本文档由 `zdi-document-architecture` prompt 自动生成，如有变更请重新运行该 prompt 以更新。

---

## 1. 项目简介

ZDIpy_WebUI 是一个 Zeeman Doppler Imaging（ZDI）反演工具，用于从恒星 Stokes I/V LSD 线轮廓观测数据中重建恒星表面亮度图和磁场图。  
项目同时提供 **WebUI**（FastAPI + 原生 JS，默认端口 7860）和 **CLI** 两种运行方式，底层共享同一套科学计算引擎和配置系统。  
输入格式为 Donati LSD 格式的 `.prof`/`.norm` 文件；输出为磁场球谐系数、亮度图数组以及拟合结果摘要文件。

---

## 2. 顶层目录结构

```
ZDIpy_WebUI/
├── app.py                  # 主入口：解析 CLI 参数，选择 WebUI / CLI 模式
├── zdi_runner.py           # run_zdi() — CLI 兼容的 ZDI 运行器
├── config_loader.py        # ZDIConfig — 将 config.json 解析为参数对象
├── requirements.txt        # 依赖声明
├── start_webui.command     # macOS 双击启动脚本
├── config/
│   └── config.json         # 默认配置文件（所有可调参数）
├── api/                    # FastAPI 层（REST + SSE）
│   ├── server.py           # FastAPI 应用工厂，挂载所有路由和静态文件
│   ├── models.py           # Pydantic 请求/响应模型
│   └── routes/             # 按资源拆分的路由模块
│       ├── config.py       # GET/PUT /api/config
│       ├── run.py          # POST /api/run，GET /api/run/status|stream (SSE)
│       ├── observations.py # GET/POST/DELETE /api/observations
│       ├── results.py      # GET /api/results/{profiles|magnetic|brightness}
│       └── plots.py        # GET /api/plots/{profiles|surface_map|light_curve}
├── core/                   # 科学计算引擎（不得 import api/ 或 frontend/）
│   ├── geometryStellar.py  # 恒星表面网格与多相位可见性批量计算
│   ├── brightnessGeom.py   # 亮度图（像素级，limb + 重力昏暗加权）
│   ├── magneticGeom.py     # 磁场球谐展开与 Legendre 多项式批量计算
│   ├── lineprofileVoigt.py # Voigt 轮廓谱线模型与 Stokes I/V 合成
│   ├── mainFuncs.py        # 遗留参数解析器 + 主拟合循环 + 文件 I/O
│   ├── memSimple3.py       # 向后兼容 shim，转发至 core/mem/zdi_adapter
│   ├── observations.py     # SpectralObservation / LightCurveObservation 数据类
│   ├── readObs.py          # obsProf — LSD 文件读取器（Donati 格式）
│   ├── mem/                # MEM 优化引擎子包
│   │   ├── generic.py      # MEMOptimizer — 纯数学 MEM 算法（Skilling & Bryan 1984）
│   │   ├── zdi_adapter.py  # ZDI 特定接口：向量打包/解包、响应矩阵组装
│   │   ├── saim_adapter.py # SAIM（熵目标）控制程序向后兼容封装
│   │   ├── optimization.py # ResponseMatrixCache（LRU）、StabilityMonitor
│   │   ├── iteration_manager.py  # ConvergenceChecker、IterationManager
│   │   └── monitoring.py   # IterationHistory、ProgressMonitor
│   └── plotting/           # 绘图后端子包
│       ├── base.py         # PlotBackend 抽象基类
│       ├── data.py         # ProfilePlotData / SurfaceMapData / LightCurvePlotData
│       ├── matplotlib_backend.py  # CLI 绘图 → matplotlib Figure
│       └── plotly_backend.py      # WebAPI 绘图 → Plotly JSON dict
├── pipeline/
│   ├── pipeline.py         # ZDIPipeline — 完整反演流程的封装类
│   └── result.py           # ZDIResult — 反演结果数据类
├── frontend/               # 原生 HTML/CSS/JS 单页应用
│   ├── index.html          # SPA 入口，含 4 个功能 Tab 的骨架
│   ├── css/app.css         # 深色主题样式
│   └── js/
│       ├── app.js          # Tab 导航、apiFetch() 封装、Plotly 公共布局
│       ├── config.js       # 配置表单（根据 CONFIG_SCHEMA 动态渲染）
│       ├── observations.js # LSD 文件上传/列表/删除 + 观测条目编辑器
│       ├── run.js          # 运行控制 + SSE 实时日志流
│       └── results.js      # Plotly 结果渲染（线轮廓 / 表面图 / 光变曲线）
├── LSDprof/                # 示例 LSD 轮廓文件（.prof 格式）
├── webui_legacy/           # 已废弃的 Gradio WebUI（保留供参考，不再使用）
├── scripts/                # 独立实用脚本
│   └── check_geom_consistency.py
├── tests/                  # 单元测试
│   ├── test_geometry.py
│   └── test_mag_brightness.py
└── docs/
    └── MEM.md              # MEM 算法说明
```

---

## 3. 启动路径

### 3.1 WebUI 模式

```
python app.py [--port 7860] [--share]
       │
       └─→ uvicorn.run("api.server:app", host="127.0.0.1", port=7860)
                │
                ├─ GET /          → 返回 frontend/index.html（静态文件）
                ├─ /css, /js      → 静态资源挂载
                └─ /api/*         → FastAPI 路由（见第 7 节）
```

浏览器打开 `http://127.0.0.1:7860`，通过 Tab 页面与 REST API 交互；`/docs` 可查看 Swagger 文档。

### 3.2 CLI 模式

```
python app.py --cli [--config config/config.json] [--forward-only] [--verbose 1]
       │
       └─→ zdi_runner.run_zdi(config_path, forward_only, verbose)
                │
                ├─→ config_loader.ZDIConfig(config_path)     # 加载并验证配置
                └─→ pipeline.ZDIPipeline(par, ...).run()      # 执行反演流程
                        │
                        └─→ 打印最终指标（iterations, chi2, entropy…）
```

也可以作为独立脚本直接运行：`python zdi_runner.py --config my_config.json`

---

## 4. 配置系统

`config_loader.ZDIConfig` 读取 `config/config.json` 并将所有路径解析为绝对路径。  
`ZDIConfig` 接口与遗留的 `core/mainFuncs.readParamsZDI` 完全兼容，下游代码无需修改。

| 顶层键 | 类型 | 含义 |
|--------|------|------|
| `star.inclination_deg` | float | 恒星自转轴倾角（度） |
| `star.vsini_kms` | float | 投影自转速度（km/s） |
| `star.period_days` | float | 自转周期（天） |
| `star.differential_rotation_rad_per_day` | float | 较差自转系数（rad/day） |
| `star.mass_msun` | float | 质量（太阳质量） |
| `star.radius_rsun` | float | 赤道半径（太阳半径） |
| `grid.nRings` | int | 表面网格纬度环数（越大精度越高、越慢） |
| `inversion.target_form` | `"C"` / `"E"` | 拟合目标：C = 目标 χ²，E = 目标熵 |
| `inversion.target_value` | float | 目标值（对应 target_form） |
| `inversion.num_iterations` | int | 最大迭代次数 |
| `inversion.test_aim` | float | Skilling-Bryan 收敛判据阈值 |
| `magnetic.fit_magnetic` | 0/1 | 是否拟合磁场图 |
| `magnetic.l_max` | int | 球谐展开截断阶数 ℓ_max |
| `magnetic.default_bent` | float | 磁场熵的默认斜率（G） |
| `magnetic.geometry_type` | `Full`/`Poloidal`/`PotTor`/`Potential` | 磁场几何类型约束 |
| `magnetic.init_from_file` | 0/1 | 从文件初始化磁场系数 |
| `magnetic.init_file` | str | 初始磁场系数文件路径 |
| `brightness.fit_brightness` | 0/1 | 是否拟合亮度图 |
| `brightness.chi2_scale_I` | float | Stokes I 的 χ² 权重缩放 |
| `brightness.entropy_scale` | float | 亮度熵的缩放系数 |
| `brightness.entropy_form` | 1/2 | 熵公式：1 = 传统图像熵，2 = 填充因子熵 |
| `brightness.default_bright` | float | 亮度默认值（均匀分布初始值） |
| `brightness.max_bright` | float | 亮度上限 |
| `line_model.estimate_strength` | 0/1 | 从等效宽度自动估计谱线强度 |
| `line_model.model_file` | str | 谱线模型参数文件 |
| `line_model.wavelength_nm` | float | 中心波长（nm） |
| `line_model.line_strength` | float | 谱线强度（深度） |
| `line_model.gauss_width_kms` | float | 高斯 Doppler 宽度（km/s） |
| `line_model.lorentz_width_fraction` | float | Lorentz 宽度因子 γ/σ |
| `line_model.lande_g` | float | 有效 Landé g 因子 |
| `line_model.limb_darkening` | float | 临边昏暗系数 ε |
| `line_model.gravity_darkening` | float | 重力昏暗系数 |
| `instrument.spectral_resolution` | float | 仪器分辨率 R |
| `velocity_grid.vel_start_kms` | float | 速度网格左端（km/s） |
| `velocity_grid.vel_end_kms` | float | 速度网格右端（km/s） |
| `observations.jdate_ref` | float | 参考儒略日（用于计算旋转相位） |
| `observations.files[]` | list | 观测文件列表：`filename`, `jdate`, `vel_center_kms` |
| `output.*` | str | 各输出文件路径（mag_coeff_file, bright_map_file 等） |
| `light_curve` (可选) | dict | 光变曲线拟合控制：`fit_light_curve`, `chi2_scale_lc`, `files[]` |

---

## 5. 科学核心模块（`core/`）

> 规则：`core/` 内模块**不得** import `api/` 或 `frontend/` 中的任何符号。

### 5.1 几何与网格（`geometryStellar.py`）

| 类/函数 | 职责 |
|---------|------|
| `starGrid` | 构建等面积球面网格（余纬度/经度/半径），支持较差自转速度场和重力昏暗预计算 |
| `visibleGrid` | 单相位可见性（视角、投影面积、视线速度） |
| `BatchVisibleGrid` | 多相位批量可见性数组 `(N_phases, N_cells)`：`visible`, `proj_area`, `vel_rot_proj`, `view_angle` |
| `build_batch_visible_grid(par, sGrid)` | 便捷工厂函数，一次构建所有观测相位的几何 |

### 5.2 亮度图（`brightnessGeom.py`）

| 类/函数 | 职责 |
|---------|------|
| `brightMap` | 保存每个网格单元的亮度值（初始化为均匀亮度 1.0） |
| `brightMap.makeRoundSpot()` | 在亮度图上放置圆形斑点 |
| `brightMap.projected(visible_batch, limb_dark, grav_dark)` | 返回 `(N_phases, N_cells)` 加权亮度（limb + 重力昏暗 × 投影面积） |
| `SetupBrightMap(sGrid, ...)` | 从文件或默认值初始化亮度图 |
| `saveMap(briMap, filename)` | 将亮度图写入文件 |

### 5.3 磁场几何（`magneticGeom.py`）

| 类/函数 | 职责 |
|---------|------|
| `_compute_legendre_batch(m_max, n_max, x)` | 批量向量化计算关联 Legendre 多项式及其导数 `(m+1, n+1, N_cells)` |
| `magSphHarmonics` | 球谐系数存储；`initMagGeom()` 预计算所有基函数并缓存 |
| `getAllMagVectorsCart()` | 返回全部网格单元的 B 场笛卡尔分量 |
| `getAllMagDerivsCart()` | 返回 Jacobian 矩阵（用于 MEM 反演） |
| `SetupMagSphHarmonics(sGrid, ...)` | 初始化工厂（从文件或默认值） |
| `saveToFile(filename)` | 保存磁场球谐系数 |

支持的磁场几何类型：`Full`（完整）、`Poloidal`（纯极向）、`PotTor`（势场+环形）、`Potential`（势场）。

### 5.4 谱线轮廓（`lineprofileVoigt.py`）

| 类/函数 | 职责 |
|---------|------|
| `lineData` | 从文件读取谱线参数（波长、强度、高斯/Lorentz 宽度、Landé g、临边/重力昏暗系数） |
| `syntheticProf` | 单相位合成 Stokes I/V 轮廓 |
| `getAllProfDiriv(par, listGridView, vecMagCart, ...)` | 计算所有相位的合成谱和导数，返回 `setSynSpec` 集合 |

使用弱场近似下的 Voigt 轮廓；仪器展宽可内嵌于局部轮廓（`explicitConvolution=False`）。

### 5.5 主控函数（`mainFuncs.py`）

| 类/函数 | 职责 |
|---------|------|
| `readParamsZDI` | 遗留 `inzdi.dat` 格式参数解析器（接口参考，已被 `ZDIConfig` 替代） |
| `getWavelengthGrid()` | 根据观测 Vr 构建速度网格并统计数据点总数 |
| `fitLineStrength()` | 从观测等效宽度自动估计谱线强度 |
| `mainFittingLoop()` | ZDI 主拟合循环：依次调用 MEM 迭代，返回最终 `(iIter, entropy, chi2, test, meanBright, meanBrightDiff, meanMag)` |
| `saveModelProfs()` | 保存合成轮廓到文件 |
| `saveObsUsed()` | 保存实际使用的观测数据到文件 |

### 5.6 观测数据读取（`readObs.py` + `observations.py`）

| 模块/类 | 职责 |
|---------|------|
| `readObs.obsProf` | 读取单个 LSD 轮廓文件（Donati 格式：2 列仅 I，或 7 列 I/V/N） |
| `readObs.obsProfSetInRange()` | 批量读取并截取至速度范围 `[velStart, velEnd]`，应用 Vr 中心偏移 |
| `readObs.getObservedEW()` | 计算观测平均等效宽度 |
| `observations.SpectralObservation` | 数据类：单次 LSD 观测（速度网格 + Stokes I/V + 误差） |
| `observations.LightCurveObservation` | 数据类：单次光变曲线观测点 |

### 5.7 MEM 优化引擎（`core/mem/`）

| 模块/类 | 职责 |
|---------|------|
| `generic.MEMOptimizer` | 纯数学 MEM 算法（Skilling & Bryan 1984）；通过回调函数传入熵、约束、边界 |
| `zdi_adapter` | ZDI 专用接口：`packDataVector`、`packImageVector`、`packResponseMatrix`、`constantsMEM`、`setEntropyWeights`、`mem_iter`、`get_s_grads`、`updateImg` |
| `saim_adapter` | SAIM（目标熵）控制程序的向后兼容封装 |
| `optimization.ResponseMatrixCache` | 响应矩阵 LRU 缓存（按磁场几何 hash 键控） |
| `optimization.StabilityMonitor` | 梯度和步长的数值稳定性监控 |
| `iteration_manager.ConvergenceChecker` | 基于相对 χ² 变化的收敛检测 |
| `iteration_manager.IterationManager` | 迭代循环总控制器（集成进度监控和历史记录） |
| `monitoring.IterationHistory` | 逐迭代状态记录（熵、χ²、梯度、时间戳） |
| `monitoring.ProgressMonitor` | ETA 估计、停滞检测 |

`memSimple3.py` 是 `core.mem.zdi_adapter` 的向后兼容 shim，所有历史调用代码均无需修改。

### 5.8 绘图系统（`core/plotting/`）

| 模块/类 | 职责 |
|---------|------|
| `base.PlotBackend` | 抽象基类：`plot_profiles`、`plot_surface_map`、`plot_light_curve` |
| `data.ProfilePlotData` | Stokes I/V 线轮廓观测 vs 模型数据容器 |
| `data.SurfaceMapData` | 恒星表面等面积投影数据容器（亮度或 B 场分量） |
| `data.LightCurvePlotData` | 光变曲线数据容器 |
| `MatplotlibBackend` | CLI 后端 → 返回 `matplotlib.figure.Figure` |
| `PlotlyBackend` | WebAPI 后端 → 返回 Plotly JSON dict（可直接传入 `Plotly.newPlot()`） |

---

## 6. 反演流水线（`pipeline/`）

`pipeline.ZDIPipeline.run()` 按以下顺序执行：

```
ZDIPipeline.run()
  │
  ├─ 1. par.setTarget()          → 确定拟合目标（χ² / 熵）
  ├─ 2. par.setCalcdIdV()        → 确定计算 I / V 分量
  ├─ 3. par.calcCycles()         → 从儒略日计算旋转相位
  │
  ├─ 4. lineprofile.lineData(model_file)          → 加载谱线参数
  ├─ 5. readObs.obsProfSetInRange(fnames, ...)    → 加载并截取观测数据
  ├─ 6. mf.getWavelengthGrid()                   → 构建速度网格
  │
  ├─ 7. geometryStellar.starGrid(nRings, ...)     → 构建恒星表面网格
  ├─ 8. magneticGeom.SetupMagSphHarmonics(...)    → 初始化磁场球谐展开
  ├─ 9. brightnessGeom.SetupBrightMap(...)        → 初始化亮度图
  ├─ 10. geometryStellar.build_batch_visible_grid() → 批量计算所有相位可见性
  │
  ├─ 11. magGeom.getAllMagVectorsCart()           → 计算 B 场向量
  ├─ 12. magGeom.getAllMagDerivsCart()            → 计算 Jacobian（仅 fitMag=1）
  │
  ├─ 13. [可选] mf.fitLineStrength()              → 从等效宽度自动估计谱线强度
  ├─ 14. lineprofile.getAllProfDiriv()            → 生成所有相位合成谱和导数
  │
  ├─ 15. memSimple.constantsMEM()                → MEM 维度常数
  ├─ 16. memSimple.packDataVector()              → 数据向量 (Data, sig2)
  ├─ 17. memSimple.setEntropyWeights()           → 熵权重
  │
  ├─ 18. mf.mainFittingLoop()                    → MEM 主拟合循环
  │         └─ 每次迭代：mem_iter() → updateImg() → 更新合成谱
  │
  ├─ 19. 保存输出文件（磁场系数、亮度图、模型轮廓、观测摘要）
  │
  └─ 20. 构建并返回 pipeline.result.ZDIResult
```

`ZDIResult` 字段：`iterations`, `entropy`, `chi2`, `test`, `converged`, `bright_map` `(N_cells,)`, `mag_coeffs` `{alpha/beta/gamma: [...]}`, `synthetic_profiles`, `observed_profiles`, `light_curve_synthetic`, `metadata`。  
方法：`to_serializable()` → JSON 兼容 dict；`save(path)` → NPZ 压缩存档。

---

## 7. REST API 层（`api/`）

应用由 `api/server.py` 创建，入口为 `api.server:app`（供 uvicorn 使用）。  
根路由 `GET /` 返回 `frontend/index.html`；`/css`, `/js` 挂载静态文件。

| 方法 | 路径 | 模块 | 说明 |
|------|------|------|------|
| GET | `/api/config` | `routes/config.py` | 返回当前 config.json 内容 |
| PUT | `/api/config` | `routes/config.py` | 验证并持久化修改后的配置 |
| POST | `/api/run` | `routes/run.py` | 在后台线程启动 ZDI 流水线（防并发锁） |
| GET | `/api/run/status` | `routes/run.py` | 轮询运行状态与最新 50 条日志 |
| GET | `/api/run/stream` | `routes/run.py` | SSE 实时日志流（EventSource 接口） |
| GET | `/api/observations` | `routes/observations.py` | 列出 `LSDprof/` 目录中的文件 |
| POST | `/api/observations/upload` | `routes/observations.py` | 上传 LSD 轮廓文件 |
| DELETE | `/api/observations/{fname}` | `routes/observations.py` | 删除指定 LSD 文件 |
| POST | `/api/observations/validate` | `routes/observations.py` | 验证文件是否存在 |
| GET | `/api/results/profiles` | `routes/results.py` | 原始相位数据（前端降级回退） |
| GET | `/api/results/magnetic` | `routes/results.py` | 恒星表面磁场数据（前端降级回退） |
| GET | `/api/results/brightness` | `routes/results.py` | 亮度图数据（前端降级回退） |
| GET | `/api/plots/profiles` | `routes/plots.py` | Plotly JSON：谱线轮廓对比图 |
| GET | `/api/plots/surface_map` | `routes/plots.py` | Plotly JSON：恒星表面图（`?map_type=brightness\|radial_B\|meridional_B\|azimuthal_B`） |
| GET | `/api/plots/light_curve` | `routes/plots.py` | Plotly JSON：光变曲线对比图 |
| GET | `/docs` | FastAPI 内建 | Swagger 交互式 API 文档 |

**Pydantic 模型**（`api/models.py`）：`RunRequest`、`RunStatus`、`ObservationFile`、`ObservationsList`、`ConfigSaveResponse`、`ObsFileInfo`

**运行状态机**：`idle` → `running` → `done` / `error`，由 `_state_lock` 保护的共享 dict 维护，`_run_lock` 防止并发运行。

---

## 8. 前端（`frontend/`）

单页应用（SPA），纯原生 HTML + CSS + JS，无任何前端框架依赖（仅依赖 CDN 加载的 `Plotly.js`）。

| 文件 | 职责 |
|------|------|
| `index.html` | 骨架：4 个 Tab（Config / Observations / Run / Results）+ 顶部状态栏 |
| `js/app.js` | Tab 切换、`apiFetch()` 统一 fetch 封装（自动 JSON 序列化/错误提取）、`setStatus()` 状态栏工具函数、Plotly 公共暗色主题布局 |
| `js/config.js` | 根据 `CONFIG_SCHEMA` 数组动态渲染所有参数的手风琴表单，通过 `GET /api/config` 加载，`PUT /api/config` 保存 |
| `js/observations.js` | 拖放/浏览上传（`POST /api/observations/upload`）、文件列表/删除、观测条目结构化编辑器（read-modify-write via config API）、文件存在性验证 |
| `js/run.js` | 运行按钮（`POST /api/run`）、`EventSource` SSE 实时日志流（`GET /api/run/stream`）、运行状态徽章 |
| `js/results.js` | 调用 `/api/plots/*` 主端点获取 Plotly JSON 并渲染；若主端点 404 则降级至 `/api/results/*` 手动构造图表 |
| `css/app.css` | GitHub 深色主题样式 |

**前后端通信方式**：REST JSON（配置/观测/结果读取）+ SSE（运行日志实时推送，无 WebSocket）。  
文件上传使用 `multipart/form-data`（`FormData`），其余接口使用 `application/json`。

---

## 9. 遗留 WebUI（`webui_legacy/`）

`webui_legacy/app.py` 是早期基于 **Gradio** 的 WebUI 实现，已被 `api/` + `frontend/` 取代，**不再参与正常启动流程**。主要差异：
- 使用 Gradio 组件渲染 UI，不暴露独立 REST API
- 使用 matplotlib（Agg 后端）生成图像并以 Base64 返回
- 与现有 `frontend/` 独立并列，互不影响
- 保留供历史对比或紧急回退，正常开发不应修改

---

## 10. 关键依赖

| 包 | 版本约束 | 用途 |
|----|----------|------|
| `numpy` | ≥ 1.24 | 所有数值数组运算，向量化几何和光谱计算 |
| `scipy` | ≥ 1.10 | 关联 Legendre 多项式（`scipy.special.lpmv`）、线性代数（`scipy.linalg`） |
| `matplotlib` | ≥ 3.7 | CLI 模式绘图（`MatplotlibBackend`），遗留 WebUI 图像生成 |
| `fastapi` | ≥ 0.110 | REST API 框架（路由、Pydantic 验证、SSE 响应） |
| `uvicorn[standard]` | ≥ 0.29 | ASGI 服务器（运行 FastAPI 应用） |
| `python-multipart` | ≥ 0.0.9 | FastAPI 文件上传支持（`UploadFile`） |

> `gradio` 仅在使用 `webui_legacy/` 时需要，不在 `requirements.txt` 中声明；`plotly` 通过 CDN 加载，不在后端依赖中。

---

## 11. 模块依赖方向

单向依赖，严禁反向引用：

```
┌─────────────────────────────────────────────────────────────┐
│  frontend (HTML/CSS/JS)                                     │
│        │  REST + SSE (HTTP)                                 │
│        ↓                                                    │
│  api/routes/  ──────────────────────────→  zdi_runner.py   │
│        │                                        │           │
│        │  (读取 _state，注入 ZDIResult)           │           │
│        ↓                                        ↓           │
│  api/models.py (Pydantic)           pipeline/pipeline.py   │
│                                             │               │
│                              ┌──────────────┘               │
│                              ↓                              │
│                     core/ (科学引擎)                         │
│                    ┌──────────────────────┐                 │
│                    │ geometryStellar      │                 │
│                    │ brightnessGeom       │                 │
│                    │ magneticGeom         │                 │
│                    │ lineprofileVoigt     │                 │
│                    │ mainFuncs            │                 │
│                    │ readObs              │                 │
│                    │ mem/  (MEM 引擎)     │                 │
│                    │ plotting/            │                 │
│                    └──────────────────────┘                 │
│                              │                              │
│                              ↓                              │
│                    config_loader.ZDIConfig                  │
│                              │                              │
│                              ↓                              │
│                    config/config.json                       │
└─────────────────────────────────────────────────────────────┘
```

**关键约束**：
- `core/` 内任何模块**不得** import `api/`、`frontend/` 的任何符号
- `pipeline/` 可 import `core/`，但**不得** import `api/`
- `api/routes/` 通过共享状态 dict（`_state` + `_state_lock`）向前端传递结果，而非直接调用 `pipeline/`
- `config_loader` 是连接配置文件与科学引擎的唯一桥梁
