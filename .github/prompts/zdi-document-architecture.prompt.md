---
description: "Explore the ZDIpy_WebUI codebase and write a complete architecture overview with brief module descriptions to .github/instructions/project-architecture.instructions.md"
name: "Document Project Architecture"
argument-hint: "Optional focus area, e.g. 'core only' or 'API and frontend'"
agent: "agent"
tools: [file_search, grep_search, read_file, list_dir, semantic_search, create_file, replace_string_in_file]
---

你是本项目的架构文档员。请系统性地探索 ZDIpy_WebUI 代码库，并在 `.github/instructions/project-architecture.instructions.md` 中写出一份完整的架构总览与简要说明。

## 任务目标

生成（或覆盖）文件：`.github/instructions/project-architecture.instructions.md`

该文件应作为所有贡献者快速理解项目全貌的入口文档，内容需准确、简洁、可维护。

---

## 探索步骤

按以下顺序逐步收集信息，**不得猜测——必须逐文件核实**：

1. **顶层结构**：列出所有顶层目录和关键文件（`app.py`、`zdi_runner.py`、`config_loader.py`、`requirements.txt`）。
2. **入口与启动流程**：阅读 `app.py` 和 `zdi_runner.py`，明确 CLI / WebUI 两种启动路径及参数传递方式。
3. **配置系统**：阅读 `config_loader.py` 和 `config/config.json`，记录所有顶层配置字段及其含义。
4. **科学核心 `core/`**：逐一阅读每个模块（`geometryStellar.py`、`brightnessGeom.py`、`magneticGeom.py`、`lineprofileVoigt.py`、`mainFuncs.py`、`memSimple3.py`、`observations.py`、`readObs.py`），总结职责与对外接口。同时梳理子包 `core/mem/` 和 `core/plotting/` 的模块职责。
5. **流水线 `pipeline/`**：阅读 `pipeline/pipeline.py` 和 `pipeline/result.py`，记录 ZDI 反演流程的调用链。
6. **API 层 `api/`**：阅读 `api/server.py`、`api/models.py` 及所有路由，记录 REST 接口和数据模型。
7. **前端 `frontend/`**：阅读主要 JS 模块（`app.js`、`config.js`、`observations.js`、`results.js`、`run.js`），记录前后端交互方式。
8. **遗留 WebUI `webui_legacy/`**：简要说明与现有 `frontend/` 的关系。
9. **依赖关系**：从 `requirements.txt` 提取关键科学/工程依赖，注明用途。

---

## 输出格式

生成的文件须遵循以下结构（用中文撰写主体内容）：

```markdown
---
description: "ZDIpy_WebUI 完整项目架构总览，供所有贡献者快速理解模块职责、调用链与依赖关系。"
name: "Project Architecture Overview"
applyTo: ["**"]
---

# ZDIpy_WebUI 项目架构总览

> 本文档由 `zdi-document-architecture` prompt 自动生成，如有变更请重新运行并更新。

## 1. 项目简介

[2–3 句话说明项目用途：ZDI 反演、WebUI、CLI、支持的数据格式]

## 2. 顶层目录结构

[用树状列表展示目录，每项附一行职责说明]

## 3. 启动路径

### 3.1 WebUI 模式
[入口 → 配置加载 → API 服务 → 前端交互]

### 3.2 CLI 模式
[入口 → 配置加载 → pipeline 调用]

## 4. 配置系统

[config.json 字段表：字段名 | 类型 | 含义]

## 5. 科学核心模块（`core/`）

[每个模块：名称 | 职责 | 关键类/函数]

### 5.1 几何与网格
### 5.2 亮度与磁场
### 5.3 谱线轮廓
### 5.4 MEM 优化引擎（`core/mem/`）
### 5.5 绘图系统（`core/plotting/`）
### 5.6 观测数据读取

## 6. 反演流水线（`pipeline/`）

[调用链示意：入口 → 各阶段 → 结果输出]

## 7. REST API 层（`api/`）

[路由表：方法 | 路径 | 说明]

## 8. 前端（`frontend/`）

[各 JS 模块职责，前后端通信方式]

## 9. 遗留 WebUI（`webui_legacy/`）

[与现有前端的关系，是否仍在使用]

## 10. 关键依赖

[表格：包名 | 版本约束 | 用途]

## 11. 模块依赖方向

[文字或 ASCII 图说明单向依赖规则：科学引擎 ← pipeline ← API ← frontend]
```

---

## 约束

- 所有内容必须来自实际代码，不得凭印象填写。
- 若某模块内容不确定，注明"待补充"而非编造。
- 文件已存在时直接覆盖（用 `replace_string_in_file` 或 `create_file`）。
- 若用户提供了 `argument-hint` 中的聚焦范围，则仅深入展开对应章节，其余章节保留简短说明。
