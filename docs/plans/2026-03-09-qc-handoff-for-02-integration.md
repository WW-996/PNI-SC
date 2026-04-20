# 胆囊癌 PNI 单细胞项目交接文档：从 `01_qc_preprocess.ipynb` 延续到 `02_integration_annotation.ipynb`

## 1. 项目当前状态

- 项目对象：`3` 位病人、`6` 个样本的胆囊癌单细胞数据。
- 样本设计：每位病人各有 `tumor + adjacent` 配对，`PNI` 分别为 `high / medium / low`。
- 当前已完成模块：`notebooks/01_qc_preprocess.ipynb`
- 当前未开始模块：`02_integration_annotation.ipynb`

## 2. 当前数据与路径约定

### 原始输入

- 原始数据根目录：`data/raw`
- 每个样本当前实际目录命名为连字符格式，例如：
  - `data/raw/P01-T`
  - `data/raw/P01-A`
- `metadata` 中样本命名为下划线格式，例如：
  - `P01_T`
  - `P01_A`
- `01_qc_preprocess.ipynb` 已经补丁支持下划线/连字符自动兼容，不需要再手工改 metadata。

### 关键输入文件

- 样本 metadata：`templates/sample_metadata_template.csv`
- QC 阈值：`templates/sample_qc_thresholds_template.csv`

### 当前 QC 输出对象

- 清洗后的单样本对象索引：`results/tables/qc/qc_clean_object_index.tsv`
- 清洗后的 `h5ad` 对象：`results/intermediate/qc_clean/*.h5ad`
- 非整合合并对象：`results/intermediate/qc_merged_unintegrated.h5ad`
- 批次诊断对象：`results/intermediate/qc_merged_unintegrated_diagnostic.h5ad`

## 3. `01_qc_preprocess.ipynb` 已完成内容

- 按样本读取 `10x` 原始矩阵
- 计算基础 QC 指标：`n_genes_by_counts`、`total_counts`、`percent_mt`
- 按样本独立阈值完成基础 QC 过滤
- 运行 `Scrublet` 完成 doublet 检测与剔除
- 生成批次/样本效应诊断图
- 导出清洗后的单样本对象与结果表

### 当前 notebook 的重要实现约定

- notebook 标题和 markdown 为中文
- 代码里的 warning / exception / 运行反馈使用英文，避免 Windows 编码问题
- `Scrublet` 已配置为 `use_approx_neighbors=False`，因此当前不依赖 `annoy`
- 批次诊断中的 HVG 计算使用 `flavor="seurat"`，避免 `scikit-misc` 依赖

## 4. 当前环境状态

- 推荐内核：`.conda (Python 3.11.15)`
- 当前已补齐并验证可导入的核心依赖：
  - `ipykernel`
  - `numpy`
  - `pandas`
  - `scanpy`
  - `seaborn`
  - `matplotlib`
  - `scipy`
  - `scrublet`
  - `cython`
- 当前 `scrublet` 已实际运行通过
- `annoy` 未安装，但当前 notebook 不需要它

## 5. `01` 阶段的结果摘要

### 样本级 QC 结果

来源：`results/tables/qc/sample_qc_summary.tsv`

| sample_id | tissue_type | pni_level | n_cells_raw | n_cells_after_qc | n_cells_after_doublet_removal | doublet_rate_predicted |
|---|---|---:|---:|---:|---:|---:|
| P01_A | adjacent | high | 8299 | 7761 | 7751 | 0.0013 |
| P01_T | tumor | high | 9619 | 8396 | 8391 | 0.0006 |
| P02_A | adjacent | medium | 6304 | 4843 | 4784 | 0.0122 |
| P02_T | tumor | medium | 6867 | 5475 | 5423 | 0.0095 |
| P03_A | adjacent | low | 8348 | 7542 | 7530 | 0.0016 |
| P03_T | tumor | low | 18395 | 16914 | 16914 | 0.0000 |

### 基础观察

- 所有样本均顺利通过 QC，没有样本触发 `risk_low_cell_count=True`
- 样本间细胞量差异明显，`P03_T` 细胞数显著高于其他样本
- `P02_A` 和 `P02_T` 的 QC 后保留比例相对更低，提示其技术质量或文库特征更特殊
- doublet 预测率整体偏低，但 `P02_A`、`P02_T` 相对更高

### 当前阈值策略

来源：`results/tables/qc/sample_thresholds_applied.tsv`

- 所有样本当前统一采用：
  - `gene_low_cutoff = 300`
  - `umi_low_cutoff = 500`
  - `mt_cutoff = 20`
- `tumor` 当前使用更宽的高端上限：
  - `gene_high_cutoff = 8000`
  - `umi_high_cutoff = 50000`
- `adjacent` 当前使用：
  - `gene_high_cutoff = 7000`
  - `umi_high_cutoff = 40000`

## 6. 批次/样本效应判断

来源图像：

- `results/figures/qc/batch_diagnostic_umap__sample_id.png`
- `results/figures/qc/batch_diagnostic_umap__patient_id.png`
- `results/figures/qc/batch_diagnostic_umap__sequencing_batch.png`
- `results/figures/qc/batch_diagnostic_umap__tissue_type.png`
- `results/figures/qc/technical_metrics_by_sample.png`

### 当前判断

- **不能把当前问题简单表述为“测序批次效应很强”**，因为当前 `sequencing_batch` 只有一个水平：`B1`
- 当前最明显的是：
  - `sample effect`
  - `patient effect`
  - 技术深度/复杂度差异
- `sample_id` 图显示多个区域由单一样本明显主导
- `patient_id` 图显示 `P02` 有较强独立性，提示个体效应较明显
- `tissue_type` 图说明 `tumor vs adjacent` 有真实信号，但并未完全压过样本来源差异
- `technical_metrics_by_sample.png` 显示不同样本在 `n_genes`、`UMI`、`percent_mt` 上差异明显，提示技术异质性客观存在

### 推荐结论表述

- 建议表述为：
  - `There is evident sample/patient-driven separation before integration.`
  - `Tumor-versus-adjacent signal exists, but sample-specific variation remains substantial.`
- 不建议直接写成：
  - `strong sequencing batch effect`

## 7. 后续 `02_integration_annotation.ipynb` 的设计起点

### 目标

- 在尽量保留 `tumor vs adjacent` 以及恶性/基质真实差异的前提下，降低样本/个体来源主导的分离
- 输出一个适合后续聚类、注释、恶性上皮提取和 TME 分析的整合对象

### 推荐输入

- 优先直接读取：`results/intermediate/qc_clean/*.h5ad`
- 通过：`results/tables/qc/qc_clean_object_index.tsv` 建立样本加载顺序

### 推荐在 `02` 中完成的核心任务

- 合并所有清洗后的 `h5ad`
- 比较整合前后表现：
  - 未整合对象
  - `Harmony`（推荐优先尝试）
  - 必要时再看 `BBKNN` 或其他轻量方法
- 用整合前后 UMAP 对比判断：
  - 样本混合是否改善
  - `tumor vs adjacent` 是否被过度抹平
- 选择一个“不过度校正”的方案进入聚类与大类注释

### 对 `02` 的约束

- 不以“把所有样本完全混平”为目标
- 不把病人差异完全视为纯技术噪音
- 要保留后续恶性上皮、CAF、髓系、T/NK 等亚群的真实结构

## 8. 已知技术坑与处理记录

- **路径命名坑**：原始目录为 `P01-T`，metadata 为 `P01_T`
  - 已在 `01` 中兼容处理
- **`scrublet` 依赖坑**：默认安装尝试拉取 `annoy`，在 Windows 上需要 C++ 编译工具
  - 已改为 `use_approx_neighbors=False`
  - 当前 notebook 可正常运行，不依赖 `annoy`
- **`seurat_v3` HVG 依赖坑**：会要求 `scikit-misc`
  - 已改为 `flavor="seurat"`
- **代码中文字符串乱码坑**：Windows 当前链路里代码内中文容易变成 `?`
  - 当前约定为：markdown 中文，代码提示英文

## 9. 下次重新打开 Codex 时的推荐接入方式

可以直接告诉 Codex：

> 请先阅读 `docs/plans/2026-03-09-qc-handoff-for-02-integration.md`，基于当前 `01_qc_preprocess.ipynb` 的结果继续设计并实现 `02_integration_annotation.ipynb`。

如果希望更明确一点，可以补一句：

> `02` 需要围绕样本/病人效应较明显、但 tumor-vs-adjacent 信号仍存在这一现状，优先比较未整合与 Harmony 的效果，避免过度校正。

## 10. 当前最关键的下一步

- 以 `results/intermediate/qc_clean/*.h5ad` 为输入，开始撰写 `02_integration_annotation.ipynb`
- `02` 的第一部分应先做：
  - QC-clean 对象读取
  - 非整合 UMAP 复现
  - `sample_id / patient_id / tissue_type` 上色复核
  - 再进入整合方法比较

