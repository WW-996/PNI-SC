# 02_integration_annotation 重构设计

**日期**：2026-03-10
**版本**：v3.0（重构版）
**项目**：胆囊癌 PNI 单细胞分析（PNI SC）
**状态**：已批准，待实现

---

## 背景与目标

在 `01_qc_preprocess.ipynb` 完成逐样本 QC 后，本 notebook 承担从 QC-clean 数据出发的完整"整合 → 诊断 → 大类注释 → 谱系精细注释"全流程。

核心设计原则：
- **半自动化注释**：自动出图，人工判断，手动填写 cluster→亚型映射，再由代码应用
- **所有谱系统一模式**：相同的4-cell结构，降低维护成本
- **小样本保护**：3个患者的数据需要在每步可视化中暴露患者贡献，避免患者特异性假象

---

## Notebook 结构总览

```
Section 0:  配置中心 & 共用工具函数
Section 1:  数据载入与合并
Section 2:  预处理 & Harmony 整合 & 细胞周期评分
Section 3:  整合效果诊断（矫正前后对比 + 患者贡献）
Section 4:  全局聚类 & 大类注释（手动映射）→ 中间检查点保存
─────────────── 谱系精细注释（统一模式）───────────────
Section 5:  Epithelial（恶性亚型细化）
Section 6:  Myeloid（TAM/DC/Neutrophil）
Section 7:  Lymphoid（T_NK + B_cell + plasma_cell）
Section 8:  Fibroblast / CAF
Section 9:  Mesenchymal/Neural-Cres（神经/施旺相关群体排查与细化）
Section 10: Endothelial
────────────────────────────────────────────────────────
Section 11: 汇总注释回写 & 最终导出
```

---

## 各 Section 详细设计

### Section 0 — 配置 & 共用工具函数

**Config 类**（集中管理所有路径与超参数）：
- 路径：`QC_INDEX`、`OUT_DIR`、`FIG_DIR`、`TAB_DIR`
- 参数：`RANDOM_SEED=2026`、`N_HVG=3000`、`N_PCS=30`、`RESOLUTION_MAIN=0.6`、`RESOLUTION_SUBSET=0.4`

**共用工具函数**：

| 函数 | 说明 |
|------|------|
| `plot_markers(adata, marker_dict, title)` | dotplot + violin + feature UMAP 三联图 |
| `plot_patient_contribution(adata, cluster_key)` | 堆叠柱状图，X=cluster，Y=比例，色=patient_id |
| `apply_annotation(adata, mapping, col_name)` | 将 cluster→亚型 映射写入 obs |
| `run_subset_pipeline(adata, batch_key, resolution)` | HVG→PCA→Harmony→Neighbors→UMAP→Leiden |
| `make_annotation_template(adata, cluster_key)` | 自动生成 `{cluster_id: ""}` 模板并打印 |

---

### Section 1 — 数据载入与合并

```
读取 results/tables/qc/qc_clean_object_index.tsv
→ 逐样本载入 h5ad，附加元数据列：
    sample_id / patient_id / tissue_type / pni_level
→ sc.concat(join="outer", fill_value=0, index_unique="_")
→ 保存 adata.layers["counts"]（原始 count 备份）
→ 打印：总细胞数、各样本细胞数分布
```

---

### Section 2 — 预处理 & Harmony 整合 & 细胞周期评分

```
normalize_total(target_sum=1e4) → log1p
→ HVG 筛选（n_top_genes=3000，flavor="seurat"）
→ PCA（n_comps=30，use_highly_variable=True）
→ 细胞周期评分：
    sc.tl.score_genes_cell_cycle(adata, s_genes, g2m_genes)
    （用于后续判断是否需要回归，此处不自动回归，人工决定）
→ Harmony 整合（key=sample_id）→ X_pca_harmony
→ Neighbors（use_rep="X_pca_harmony"）→ UMAP → Leiden(0.6)
```

**注意**：细胞周期评分结果将在 Section 3 中叠加可视化，由分析者判断增殖驱动的cluster是否需要特殊处理。

---

### Section 3 — 整合效果诊断

**矫正前后并排 UMAP（4图）**：
1. `X_pca` colored by `sample_id`（矫正前）
2. `X_pca_harmony` colored by `sample_id`（矫正后）
3. colored by `patient_id`
4. colored by `tissue_type`

**细胞周期叠加 UMAP**：
- colored by `phase`（G1/S/G2M）
- 判断是否存在增殖驱动的大型cluster

**患者贡献柱状图**：
- `plot_patient_contribution(adata, "leiden")`
- 识别患者特异性cluster（>80%来自单个患者 → 注释时标记为`patient_specific`）

---

### Section 4 — 全局聚类 & 大类注释

当前 notebook 已不再采用基于 marker score 的自动大类注释表述，而是以全局 Leiden 聚类、marker 图谱和 UMAP 人工判读为基础，在 `Section 4.2` 中通过 `cluster2celltype` 手动映射写入 `adata.obs["broad_celltype"]`。

```
先绘制全局 marker 图谱与 UMAP，人工确认各 Leiden cluster 的生物学属性
→ 在 Section 4.2 中手动定义 cluster2celltype
→ 写入 adata.obs["broad_celltype"]
→ UMAP：colored by leiden / broad_celltype / pni_level
→ 保存大类注释检查点 annotated_broad.h5ad
```

**当前 `Section 4.2` 已确定的大类映射**：

| Leiden cluster | broad_celltype |
|------|------|
| 0 | T_NK |
| 1 | T_NK |
| 2 | epithelial |
| 3 | myeloid |
| 4 | fibroblast |
| 5 | B_cell |
| 6 | myeloid |
| 7 | endothelial |
| 8 | Mesenchymal/Neural-Cres |
| 9 | T_NK |
| 10 | T_NK |
| 11 | plasma_cell |
| 12 | plasma_cell |
| 13 | myeloid |
| 14 | myeloid |

**当前大类细胞数（后续谱系拆分起点）**：

| broad_celltype | n_cells |
|------|------:|
| T_NK | 20823 |
| myeloid | 9800 |
| epithelial | 6544 |
| fibroblast | 3610 |
| B_cell | 3431 |
| endothelial | 2308 |
| plasma_cell | 2252 |
| Mesenchymal/Neural-Cres | 2076 |

**关于双细胞（Doublet）的保守说明**：

- 大类注释后仍保留检查并剔除明显双细胞簇的机制。
- 但按照 `notebooks/02_integration_annotation_v3.ipynb` 当前 `Section 4.2` 的实际结果，本轮没有任何 cluster 被正式标注为 `Doublet`，实际剔除细胞数为 `0`。
- 因此后续计划中不得再将 `cluster 10` 或其他 cluster 预设为既定双细胞结论。

**Section 4 输出与检查点**：

- 主检查点文件：`results/intermediate/integration/annotated_broad.h5ad`
- 后续所有谱系子集均以 `broad_celltype` 为切分依据，而不再使用 `cell_class_auto`。

---

### Sections 5–10 — 谱系精细注释（统一4-cell模式）

各 Section 的入口以当前 `broad_celltype` 为准：

| Section | 入口 broad_celltype |
|------|------|
| 5 | epithelial |
| 6 | myeloid |
| 7 | T_NK + B_cell + plasma_cell |
| 8 | fibroblast |
| 9 | Mesenchymal/Neural-Cres |
| 10 | endothelial |

每个谱系的结构完全一致：

```python
# ── [Cell A] 载入 + 谱系内再整合 ──────────────────────────
sub = run_subset_pipeline(subsets["Myeloid"], batch_key="sample_id",
                          resolution=Config.RESOLUTION_SUBSET)

MYELOID_MARKERS = {
    "Macrophage_M1":  ["IL1B", "TNF", "CXCL10", "CD80", "HLA-DRA"],
    "Macrophage_M2":  ["CD163", "MRC1", "TGFB1", "CCL18", "FOLR2"],
    "TAM_SPP1":       ["SPP1", "APOE", "C1QA", "TREM2"],
    "cDC1":           ["CLEC9A", "CADM1", "XCR1"],
    "cDC2":           ["ITGAX", "CD1C", "FCER1A", "CLEC10A"],
    "pDC":            ["LILRA4", "IRF7", "CLEC4C"],
    "Neutrophil":     ["S100A8", "S100A9", "FCGR3B", "CSF3R"],
    "Monocyte":       ["CD14", "LYZ", "S100A4", "CCR2"],
}

# ── [Cell B] 自动出图 ─────────────────────────────────────
plot_markers(sub, MYELOID_MARKERS, title="Myeloid Markers")
plot_patient_contribution(sub, "leiden")
# ↑ 在此查看 dotplot / violin / feature UMAP + 患者贡献图

# ── [Cell C] ⬇ 人工填写映射 ──── 暂停点 ──────────────────
make_annotation_template(sub, "leiden")   # 自动打印模板，复制后填写
cluster_annotation = {
    "0": "",   # ← 根据上图填写亚型名称
    "1": "",
    # ...
}

# ── [Cell D] 应用注释 + 保存 + 汇总图 ────────────────────
sub = apply_annotation(sub, cluster_annotation, "cell_subtype")
sc.pl.umap(sub, color=["cell_subtype", "patient_id", "pni_level"])
sub.write_h5ad(Config.OUT_DIR / "subset_myeloid_annotated.h5ad")
```

**各谱系的 marker 字典预设范围**：

| Section | 谱系 | 预期亚型 |
|---------|------|---------|
| 5 | Epithelial | likely_malignant / normal_biliary / malignant_EMT / malignant_PNI |
| 6 | Myeloid | TAM_M1 / TAM_M2 / TAM_SPP1 / cDC1 / cDC2 / pDC / Neutrophil / Monocyte |
| 7 | Lymphoid | CD8_effector / CD8_exhausted / CD4_Tconv / Treg / NK / B_cell / plasma_cell |
| 8 | Fibroblast | myoCAF / iCAF / apCAF / vCAF / normal_fibroblast |
| 9 | Mesenchymal/Neural-Cres | 从 `Mesenchymal/Neural-Cres` 入口进一步排查与细化神经/施旺相关群体，可根据结果分出 Schwann_repair / Schwann_promyelinating / Schwann_PNI / neural_progenitor 或保留更保守标签 |
| 10 | Endothelial | normal_EC / tumor_EC / tip_EC |

---

### Section 11 — 汇总回写 & 最终导出

```
将各谱系 cell_subtype 注释合并回主 adata.obs["cell_subtype"]
→ 双层 UMAP：broad_celltype（大类）+ cell_subtype（亚型）
→ 各谱系细胞数 & 亚型分布汇总表（per sample / per patient）
→ PNI level × cell_subtype 组成比例热图
→ 保存 integrated_final_annotated_v2.h5ad
```

---

## 关键设计决策

| 决策 | 选择 | 理由 |
|------|------|------|
| 人工介入方式 | Notebook 分段执行，手动填写字典 | 最灵活，无额外依赖 |
| Marker 来源 | 各谱系在 notebook 内维护独立字典 | 方便迭代修改 |
| 批次矫正 | Harmony（每级均执行） | 全局和谱系内均需要 |
| 细胞周期 | 评分但不自动回归 | 保留人工判断空间 |
| 患者贡献 | 每个谱系必出柱状图 | 3患者数据的必要保护 |
| 中间保存 | Section 4 后强制保存 `annotated_broad.h5ad` | 避免重跑耗时整合步骤，并固定后续谱系拆分入口 |
| 内皮细胞 | 加入为 Section 10，轻量注释 | TME完整性，轻量处理 |

---

## 输出文件清单

```
results/intermediate/integration/
├── annotated_broad.h5ad                # Section 4 检查点（broad_celltype 入口对象）
├── subset_epithelial_annotated.h5ad
├── subset_myeloid_annotated.h5ad
├── subset_lymphoid_annotated.h5ad
├── subset_fibroblast_annotated.h5ad
├── subset_mesenchymal_neural_annotated.h5ad
├── subset_endothelial_annotated.h5ad
└── integrated_final_annotated_v2.h5ad  # 最终主对象

results/figures/integration/
├── harmony_before_after_umap.png
├── cell_cycle_umap.png
├── patient_contribution_global.png
├── global_umap_annotation.png
└── subset_<lineage>_markers.png × 6

results/tables/integration/
├── lineage_summary.tsv
└── celltype_composition_per_sample.tsv
```

---

## 依赖环境

`.conda` 环境需确认安装：
- `harmonypy`、`leidenalg`、`igraph`（整合必需）
- `scanpy >= 1.9`（`score_genes_cell_cycle` 内置）
