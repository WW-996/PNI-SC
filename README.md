# PNI SC

胆囊癌单细胞分析项目骨架，面向 `3` 位病人、`6` 个样本（每位病人 `tumor + adjacent` 配对），主线围绕 `PNI high / medium / low` 的肿瘤微环境重塑。

## 项目结构

- `docs/plans/2026-03-09-gbc-pni-scrna-analysis-design.md`：完整分析蓝图
- `notebooks/01_qc_preprocess.ipynb`：逐样本严格 QC notebook
- `notebooks/README.md`：notebook 运行说明
- `templates/sample_metadata_template.csv`：样本信息模板
- `templates/sample_qc_thresholds_template.csv`：逐样本 QC 阈值模板
- `templates/contrast_plan.tsv`：推荐比较设计模板
- `references/marker_reference.tsv`：细胞注释与机制候选 marker 参考
- `references/signature_reference.tsv`：可迁移到空间组学的 signature 参考

## 推荐使用方式

1. 先填写 `templates/sample_metadata_template.csv`
2. 按样本填写 `templates/sample_qc_thresholds_template.csv`
3. 运行 `notebooks/01_qc_preprocess.ipynb` 完成严格 QC
4. 按 `docs/plans/2026-03-09-gbc-pni-scrna-analysis-design.md` 继续后续分析
5. 用 `references/marker_reference.tsv` 做注释与候选筛选
6. 用 `references/signature_reference.tsv` 生成模块分数并为后续空间组学做映射

## 环境说明

- 推荐使用项目内的 `.conda` 环境运行 notebook
- `01_qc_preprocess.ipynb` 已补齐并验证的核心依赖包括：`ipykernel`、`numpy`、`pandas`、`scanpy`、`seaborn`、`matplotlib`、`scipy`、`scrublet`、`cython`
- 当前 `Scrublet` 在 notebook 中使用 `use_approx_neighbors=False`，因此不依赖 `annoy` 即可运行
- 如需手动指定解释器，可使用 `D:\Codex\Project\PNI SC\.conda\python.exe`

## 当前假设

- 数据类型：`10x scRNA-seq`
- 标签粒度：`PNI` 为病人级标签
- 研究阶段：以机制筛选和假设生成优先，后续兼容空间组学验证

## Quick environment setup (minimal)

```bash
python -m venv .venv
.venv/Scripts/activate
python -m pip install -r requirements.txt
```

For notebook variant generation:

```bash
python tmp/make_no_p02_highlow_notebooks.py
python tmp/make_no_p02_highlow_notebooks.py --root "D:\\path\\to\\PNI SC"
```

You can also override root path with environment variable:

```bash
set PNI_SC_ROOT=D:\\path\\to\\PNI SC
python tmp/make_no_p02_highlow_notebooks.py
```
