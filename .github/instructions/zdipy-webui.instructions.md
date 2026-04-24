---
description: "Use when editing ZDIpy_WebUI Python code, FastAPI/REST backend, native-JS frontend, CLI entrypoints, JSON config loading, or ZDI pipeline integration. Enforces project structure, naming, safe edit boundaries, and run/verification habits."
name: "ZDIpy WebUI Project Rules"
applyTo: ["**/*.py", "**/*.js", "**/*.css", "**/*.html", "config.json", "README.md", "start_webui.command"]
---
# ZDIpy WebUI Instructions

## 项目入口与启动

- 唯一入口是根目录的 `app.py`；WebUI 启动路径为 `python app.py [--port 7860]`，CLI 路径为 `python app.py --cli [--config config.json]`。
- `start_webui.command` 调用 `python app.py`，通过 `api.server:app`（FastAPI + uvicorn）提供服务，默认端口 7860。
- `webui/app.py` 是已废弃的 Gradio UI，**不参与任何启动流程，不应被引用或修改**。Gradio 不在 `requirements.txt` 中。

## 配置系统

- 所有用户可调参数均通过根目录 `config.json` + `config_loader.py`（`ZDIConfig`）传递；禁止添加临时配置来源。
- 新增配置字段时，必须同步更新：`config.json`（新键）、`config_loader.py`（`_parse` 方法解析）、`pipeline/pipeline.py`（消费侧）、`frontend/js/config.js`（`CONFIG_SCHEMA` 中的字段定义）。
- `config.json` 位于根目录（非 `config/` 子目录）。

## 科学核心 `core/`

- 将 `core/` 视为科学引擎：不得 import `api/` 或 `frontend/` 的任何符号。
- 谱线模型通过 `line_model.model_type` 选择：`"voigt"`（默认，弱场近似）或 `"unno"`（Unno-Rachkovsky，Milne-Eddington 完整偏振辐射转移）。`pipeline.py` 负责根据 `par.line_model_type` 分支选择 `lineData` / `lineDataUnno` 和对应的 `getAllProfDiriv*`。
- `core/geometry/`、`core/line_models/`、`core/fitting.py` 是重构后的主代码；`core/geometryStellar.py`、`core/lineprofileVoigt.py`、`core/mainFuncs.py` 是向后兼容 shim，不要直接修改这些 shim 的实质逻辑。
- 保留科学栈既有命名风格（camelCase 偏重）和数值默认值，非数值正确性或接口兼容性需求时不做广泛重构。
- 路径处理使用基于项目根目录的相对路径，遵循现有根路径解析模式。

## 反演流水线 `pipeline/`

- `ZDIPipeline.run()` 是唯一完整的 ZDI 反演流程封装点；CLI（`zdi_runner.py`）和 WebUI（`api/routes/run.py`）均通过它执行。
- `pipeline.py` 中谱线模型分支由 `getattr(par, 'line_model_type', 'voigt')` 控制，默认值保证向后兼容。

## REST API `api/`

- 路由按资源拆分在 `api/routes/`；不得在路由层直接执行重量级计算，长时任务放入后台线程。
- 运行状态机（`idle` / `running` / `done` / `error`）由 `_state_lock` 保护；`_run_lock` 防止并发运行，不得绕过。
- 修改 `PUT /api/config` 保存逻辑时，不得破坏 `observations.files` 的 read-modify-write 保护模式。

## 前端 `frontend/`

- 纯原生 HTML/CSS/JS，无框架；Plotly.js 由 CDN 加载，不引入其他前端依赖。
- `CONFIG_SCHEMA`（`frontend/js/config.js`）驱动配置表单的全部渲染、填充和收集；新增/移除配置项只需修改该数组。
- 条件显隐通过 `_applyAllConditionals()` + `_setupConditionalListeners()` 机制实现（`data-cond-field` / `data-cond-val` 属性）；`select` 或 `checkbox` 控制器变化时自动触发，`_populate()` 后也会触发一次。
- 修改 `frontend/` 后需递增 `index.html` 中对应资源的 `?v=N` 缓存破坏参数，确保浏览器加载最新版本。
- `input[type="number"]` 通过 `color-scheme: dark` 修复深色模式下 spinner 按钮颜色。

## CLI / WebUI 一致性

- CLI 和 WebUI 必须消费等价的配置；任何影响 CLI 行为的变更需同步反映在 WebUI 配置项中。
- 不要在 CLI 路径中添加 WebUI 独有的配置来源，反之亦然。

## 观测数据

- `LSDprof/` 文件的使用和路径映射是明确的配置行为；不要静默改变观测数据的假设（文件格式、Vr 中心偏移等）。
- 新增观测相关字段时要同步更新 `api/routes/observations.py`、`frontend/js/observations.js`，以及 `config.json` 中的 `observations.files` 结构。

## 编辑规范

- 保持编辑最小化和局部化；不做风格性全局重写。
- 旧 JF/Donati 代码的反向工程、输入映射、历史参数说明文档统一放在 `docs/JFcode/` 下；这类文档应优先记录源码位置、交互提示、物理含义和样例值，方便迁移到当前 `config.json` 体系。
- 修改影响用户可见行为（命令行参数、配置键名、API 路径）时，同步更新 `README.md`。
- 验证流程：`pip install -r requirements.txt` → `python -m pytest tests/ -q` → `python app.py` / `python app.py --cli --config config.json`。

