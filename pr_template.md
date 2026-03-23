## 新增机器学习关键靶点筛选模块

本 PR 为 TCMDATA 添加了基于机器学习的关键靶点筛选模块，定位于网络药理学分析的下游步骤——从 PPI / WGCNA 候选基因集中，结合表达谱数据筛选关键基因。

### 新增内容

**核心函数**（`R/ml_models.r`、`R/ml_utils.R`）

- `prepare_ml_data()` — 数据预处理，支持三种验证模式：
  - Mode A：全数据交叉验证（无留出集）
  - Mode B：内部训练/测试随机划分
  - Mode C：外部独立验证集
- 六种模型：`ml_lasso()`、`ml_enet()`、`ml_ridge()`、`ml_rf()`、`ml_svm_rfe()`、`ml_xgboost()`
- `run_ml_screening()` — 批量执行多种方法
- `get_ml_consensus()` — 提取被 ≥ k 个模型共同选出的基因
- `select_features()` — 对 Ridge、XGBoost 等保留全部特征的模型进行事后裁剪

**可视化**（`R/ml_plots.R`）

- 各模型的 CV 曲线、系数路径图、重要性条形图
- `plot_ml_roc()` — 多方法 ROC 曲线叠加（三种模式均支持）
- `plot_ml_venn()` — 各方法筛选基因集的韦恩图

**数据与文档**

- `covid19` 数据集（GSE157103，100 样本 × 5000 基因），用于示例与测试
- bookdown 章节 `08-Machine-Learning-analysis.Rmd`，覆盖三种模式的完整使用示例
- 所有导出函数均配有完整 roxygen 文档

### 设计说明

- 所有模型函数返回统一的 S3 对象（`tcm_ml`），具有一致的访问接口（`$selected_features`、`$importance`、`$cv_performance`、`$test_performance`）。
- 超参数调优全部在训练集内部通过交叉验证完成，测试集不参与任何模型选择过程。
- 已兼容 xgboost ≥ 2.1 / 3.x 的 API 变更（`$early_stop$best_iteration`、`$cv_predict$pred`）。

### 测试情况

在 macOS（aarch64，R 4.5，xgboost 3.2.0）上对 Mode A、B、C 进行了冒烟测试，0 error，0 warning。bookdown 章节可正常渲染，GitHub Actions 部署通过。

### 文件变更

- **新增**：`R/ml_models.r`、`R/ml_plots.R`、`R/ml_utils.R`、`data/covid19.rda`、`docs/08-Machine-Learning-analysis.Rmd`、22 个 man pages
- **修改**：`DESCRIPTION`、`NAMESPACE`、`README.md`、`.gitignore`、`docs/_bookdown.yml`
- **重命名**：`docs/08-Other-resources.Rmd` → `docs/09-Other-resources.Rmd`
