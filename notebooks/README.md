# Notebook 规划

- `01_qc_preprocess.ipynb`：逐样本 QC、doublet 检测与剔除、样本/病人效应诊断、中间结果导出。
- `02_integration_annotation_v3.ipynb`：合并 QC-clean 对象、强制运行 `Harmony` + `Leiden`、完成大类注释、保存 `annotated_broad.h5ad`，并输出 broad-level 图谱。
- `02.2_subclustering_annotation.ipynb`：读取 `annotated_broad.h5ad`，继续进行 epithelial 与 myeloid 的亚群再整合、人工判读和亚群注释导出。
- `02.3_myeloid_subclustering.ipynb`：读取 `annotated_broad.h5ad`，专门对 `broad_celltype == myeloid` 进行亚群重聚类、marker 导出与结果保存。

## 运行前准备

- 填写 `templates/sample_metadata_template.csv`
- 复制并修改 `templates/sample_qc_thresholds_template.csv`
- 将每个样本的 `10x` 原始矩阵放在 `data/raw/<sample_id>/` 或等价可读取路径中
- 先完整运行 `notebooks/01_qc_preprocess.ipynb`
- 确认 `.conda` 环境已安装 `harmonypy`、`igraph`、`leidenalg`

## 推荐内核与环境

- 推荐使用项目内核：`.conda (Python 3.11.15)`
- 当前 `01` / `02` / `02.2` 已验证可导入的核心依赖包括：`ipykernel`、`numpy`、`pandas`、`scanpy`、`seaborn`、`matplotlib`、`scipy`
- `01_qc_preprocess.ipynb` 中 `Scrublet` 已配置为 `use_approx_neighbors=False`，因此当前不要求安装 `annoy`
- `02_integration_annotation_v3.ipynb` 强制运行 `Harmony` 与 `Leiden`，不再提供跳过或回退分支
- `02.2_subclustering_annotation.ipynb` 依赖 `results/intermediate/integration/annotated_broad.h5ad`，运行前需先完成 `02_integration_annotation_v3.ipynb`
- `02.3_myeloid_subclustering.ipynb` 同样依赖 `results/intermediate/integration/annotated_broad.h5ad`，建议在完成 `02` 后直接运行用于髓系精细分析
- 如果 VS Code 提示缺少内核，可直接使用 `D:\Codex\Project\PNI SC\.conda\python.exe`

## 结果输出

- `results/intermediate/qc_clean/*.h5ad`
- `results/tables/qc/*.tsv`
- `results/figures/qc/*.png`
- `results/intermediate/integration/*.h5ad`
- `results/tables/integration/*.tsv`
- `results/figures/integration/*.png`
