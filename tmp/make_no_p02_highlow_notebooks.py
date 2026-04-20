from __future__ import annotations

import json
from pathlib import Path
from textwrap import dedent


ROOT = Path(r"d:\Codex\Project\PNI SC")
NOTEBOOKS = ROOT / "notebooks"


def to_source(text: str) -> list[str]:
    text = dedent(text).strip("\n") + "\n"
    return text.splitlines(keepends=True)


def load_notebook(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def save_notebook(path: Path, notebook: dict) -> None:
    path.write_text(json.dumps(notebook, ensure_ascii=False, indent=1), encoding="utf-8")


def update_markdown_title(cell: dict, suffix: str) -> None:
    src = "".join(cell.get("source", []))
    cell["source"] = to_source(src.rstrip() + suffix)


SHARED_SETUP = """
# === Section 0: 导入依赖包与全局配置 ===

import os
import random
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from pathlib import Path
import matplotlib.pyplot as plt


class Config:
    BASE_DIR = Path(r"D:\\Codex\\Project\\PNI SC")

    QC_INDEX = BASE_DIR / "results" / "tables" / "qc" / "qc_clean_object_index.tsv"
    OUT_DIR = BASE_DIR / "results" / "intermediate" / "integration_noP02_highlow"
    FIG_DIR = BASE_DIR / "results" / "figures" / "integration_noP02_highlow"
    TAB_DIR = BASE_DIR / "results" / "tables" / "integration_noP02_highlow"

    for d in [OUT_DIR, FIG_DIR, TAB_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    EXCLUDED_PATIENT_IDS = ["P02"]
    TUMOR_ONLY_PLOT_SUFFIX = "tumor_only"
    TUMOR_ONLY_PLOT_LABEL = " (Tumor only)"
    EXPORT_PDF = True
    CHECKPOINT_NAME = "annotated_broad_noP02_highlow.h5ad"

    RANDOM_SEED = 20260310
    DOT_SIZE = 8
    N_HVG = 3000
    N_PCS = 30
    RESOLUTION_MAIN = 0.5


np.random.seed(Config.RANDOM_SEED)
random.seed(Config.RANDOM_SEED)

sc.settings.set_figure_params(dpi=100, facecolor="white", frameon=False, vector_friendly=True)
sc.settings.verbosity = 3

print("Section 0 ready.")
print(f"Excluded patients: {Config.EXCLUDED_PATIENT_IDS}")
print(f"Figure output dir: {Config.FIG_DIR}")


def save_figure(stem, out_dir=None, bbox_inches="tight", dpi=None):
    out_dir = Path(out_dir) if out_dir is not None else Config.FIG_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    stem_path = Path(stem)
    (out_dir / stem_path.parent).mkdir(parents=True, exist_ok=True)
    fig = plt.gcf()
    fig.savefig(out_dir / stem_path.with_suffix(".png"), bbox_inches=bbox_inches, dpi=dpi)
    if Config.EXPORT_PDF:
        fig.savefig(out_dir / stem_path.with_suffix(".pdf"), bbox_inches=bbox_inches, dpi=dpi)


def tumor_only_view(adata):
    if "tissue_type" not in adata.obs.columns:
        return None
    mask = adata.obs["tissue_type"].astype(str).str.lower() == "tumor"
    if int(mask.sum()) == 0:
        return None
    return adata[mask].copy()


def iter_plot_views(adata):
    yield adata, "", ""
    tumor = tumor_only_view(adata)
    if tumor is not None and tumor.n_obs > 0:
        yield tumor, f"_{Config.TUMOR_ONLY_PLOT_SUFFIX}", Config.TUMOR_ONLY_PLOT_LABEL


def plot_markers(adata, marker_dict, title):
    for view, suffix, label in iter_plot_views(adata):
        sc.pl.dotplot(view, marker_dict, groupby="leiden", title=f"{title} - Dotplot{label}", show=False)
        save_figure(f"{title}_dotplot{suffix}")
        plt.show()

        sc.pl.stacked_violin(
            view,
            marker_dict,
            groupby="leiden",
            title=f"{title} - Violin{label}",
            swap_axes=True,
            show=False,
        )
        save_figure(f"{title}_violin{suffix}")
        plt.show()

        rep_markers = [genes[0] for genes in marker_dict.values() if genes]
        if rep_markers:
            sc.pl.umap(
                view,
                color=rep_markers,
                title=[f"{title}{label}: {marker}" for marker in rep_markers],
                ncols=4,
                show=False,
            )
            save_figure(f"{title}_feature_umap{suffix}")
            plt.show()


def plot_patient_contribution(adata, cluster_key, group_key="patient_id", filename_prefix="patient_contribution"):
    for view, suffix, label in iter_plot_views(adata):
        tmp = pd.crosstab(view.obs[cluster_key], view.obs[group_key], normalize="index")
        tmp.plot(kind="bar", stacked=True, figsize=(10, 5), colormap="tab20")
        plt.title(f"{group_key} contribution to {cluster_key}{label}")
        plt.ylabel("Proportion")
        plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", title=group_key)
        plt.tight_layout()
        save_figure(f"{filename_prefix}_{cluster_key}{suffix}")
        plt.show()


def apply_annotation(adata, mapping, col_name):
    adata.obs[col_name] = adata.obs["leiden"].map(mapping).astype("category")
    return adata


def run_subset_pipeline(adata, batch_key, resolution):
    sc.pp.highly_variable_genes(adata, n_top_genes=Config.N_HVG, flavor="seurat", batch_key=batch_key)
    sc.tl.pca(adata, n_comps=Config.N_PCS, use_highly_variable=True, random_state=Config.RANDOM_SEED)
    sce.pp.harmony_integrate(
        adata,
        batch_key,
        basis="X_pca",
        adjusted_basis="X_pca_harmony",
        random_state=Config.RANDOM_SEED,
    )
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", random_state=Config.RANDOM_SEED)
    sc.tl.umap(adata, random_state=Config.RANDOM_SEED)
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden", random_state=Config.RANDOM_SEED)
    return adata


def make_annotation_template(adata, cluster_key):
    clusters = sorted(adata.obs[cluster_key].unique(), key=int)
    template = "cluster_annotation = {\\n"
    for cluster_id in clusters:
        template += f'    "{cluster_id}": "",\\n'
    template += "}"
    print("--- 请复制下方模板并填写注释 ---")
    print(template)
"""


MAIN_LOAD = """
# 1. 读取索引表
df_index = pd.read_csv(Config.QC_INDEX, sep="\\t")
print(f"Index contains {len(df_index)} samples before filtering.")

df_index = df_index[~df_index['patient_id'].isin(Config.EXCLUDED_PATIENT_IDS)].copy()
print(f"Keeping {len(df_index)} samples after excluding {Config.EXCLUDED_PATIENT_IDS}.")

adatas = []
for _, row in df_index.iterrows():
    file_path = row["object_path"]
    sample_id = row["sample_id"]
    print(f"Loading sample: {sample_id} ({file_path})")

    ad = sc.read_h5ad(file_path)
    ad.obs["sample_id"] = sample_id
    ad.obs["patient_id"] = row["patient_id"]
    ad.obs["tissue_type"] = row["tissue_type"]
    ad.obs["pni_level"] = row["pni_level"]
    adatas.append(ad)

adata = sc.concat(adatas, join="outer", fill_value=0, index_unique="_")
adata.layers["counts"] = adata.X.copy()

remaining_patients = set(adata.obs["patient_id"].astype(str).unique())
remaining_pni_levels = set(adata.obs["pni_level"].astype(str).unique())
assert remaining_patients == {'P01', 'P03'}, f"Unexpected patients: {remaining_patients}"
assert remaining_pni_levels == {'high', 'low'}, f"Unexpected PNI levels: {remaining_pni_levels}"

print("=" * 30)
print("Data merge complete after no-P02 filtering.")
print(f"Cells: {adata.n_obs}; genes: {adata.n_vars}")
print(adata.obs["sample_id"].value_counts())
print("=" * 30)

del adatas
"""


MAIN_DIAGNOSTICS = """
import matplotlib.pyplot as plt
import pandas as pd

print("Generating diagnostic plots...")

harmony_neighbors = adata.uns["neighbors"].copy()
harmony_obsp_distances = adata.obsp["distances"].copy()
harmony_obsp_connectivities = adata.obsp["connectivities"].copy()
harmony_umap = adata.obsm["X_umap"].copy()

sc.pp.neighbors(adata, use_rep="X_pca")
sc.tl.umap(adata)
adata.obsm["X_umap_pca"] = adata.obsm["X_umap"].copy()

adata.uns["neighbors"] = harmony_neighbors
adata.obsp["distances"] = harmony_obsp_distances
adata.obsp["connectivities"] = harmony_obsp_connectivities
adata.obsm["X_umap"] = harmony_umap

for view, suffix, label in iter_plot_views(adata):
    fig, axs = plt.subplots(2, 2, figsize=(14, 12))
    sc.pl.embedding(view, basis="X_umap_pca", color="sample_id", title="1. Before Harmony (Colored by Sample)" + label, ax=axs[0, 0], show=False)
    sc.pl.umap(view, color="sample_id", size=Config.DOT_SIZE, title="2. After Harmony (Colored by Sample)" + label, ax=axs[0, 1], show=False)
    sc.pl.umap(view, color="patient_id", size=Config.DOT_SIZE, title="3. After Harmony (Colored by Patient)" + label, ax=axs[1, 0], show=False)
    sc.pl.umap(view, color="tissue_type", size=Config.DOT_SIZE, title="4. After Harmony (Colored by Tissue)" + label, ax=axs[1, 1], show=False)
    plt.tight_layout()
    save_figure(f"harmony_before_after_umap{suffix}")
    plt.show()

    sc.pl.umap(view, color=["leiden", "phase"], ncols=2, title=["Leiden Clusters" + label, "Cell Cycle Phase" + label], show=False)
    save_figure(f"cell_cycle_umap{suffix}")
    plt.show()

plot_patient_contribution(adata, cluster_key="leiden")
"""


MAIN_LINEAGE_MARKERS = """
import matplotlib.pyplot as plt

print("Preparing lineage marker plots...")

custom_markers = {
    "pan_epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "TACSTD2", "MSLN", "MUC1"],
    "pan_fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"],
    "pan_myeloid": ["LYZ", "S100A8", "C1QC", "C1QB", "CTSS", "FCER1G"],
    "pan_T_NK": ["CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R", "NKG7"],
    "B_cell": ["MS4A1", "CD79A", "CD74", "HLA-DRA", "CD37"],
    "plasma_cell": ["MZB1", "JCHAIN", "SDC1", "DERL3", "IGHG1"],
    "pan_endothelial": ["PECAM1", "KDR", "EMCN", "RAMP2", "ESAM"],
    "pan_pericyte": ["RGS5", "MCAM", "CSPG4", "PDGFRB", "NOTCH3"],
    "pan_neural_glial": ["S100B", "PMP22", "MPZ", "PLP1", "SOX10", "NGFR"],
}

valid_markers = {key: [gene for gene in genes if gene in adata.var_names] for key, genes in custom_markers.items()}

for view, suffix, label in iter_plot_views(adata):
    sc.pl.dotplot(view, valid_markers, groupby="leiden", standard_scale="var", color_map="Reds", title=f"Lineage Markers across Leiden Clusters{label}", show=False)
    save_figure(f"custom_lineage_dotplot{suffix}")
    plt.show()

    core_genes = [genes[0] for genes in valid_markers.values() if genes]
    plot_features = ["leiden"] + core_genes
    sc.pl.umap(view, color=plot_features, ncols=3, cmap="Reds", size=Config.DOT_SIZE, title="Lineage marker overview" + label, show=False)
    save_figure(f"custom_lineage_umap{suffix}")
    plt.show()
"""


MAIN_SCHWANN = """
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

target_genes = ["SOX10", "MPZ", "PLP1", "S100B", "NGFR", "L1CAM"]
genes_in_data = [gene for gene in target_genes if gene in adata.var_names]

if genes_in_data:
    for view, suffix, label in iter_plot_views(adata):
        sc.pl.umap(
            view,
            color=genes_in_data,
            cmap="Reds",
            size=50,
            title=[f"{gene} Expression{label}" for gene in genes_in_data],
            show=False,
        )
        save_figure(f"schwann_feature_umap{suffix}")
        plt.show()
"""


MAIN_BROAD_ANNOTATION = """
print("Writing broad annotations and saving checkpoint...")

cluster2celltype = {
    "0": "T_NK",
    "1": "T_NK",
    "2": "epithelial",
    "3": "myeloid",
    "4": "fibroblast",
    "5": "B_cell",
    "6": "myeloid",
    "7": "endothelial",
    "8": "Mesenchymal/Neural-Cres",
    "9": "T_NK",
    "10": "T_NK",
    "11": "plasma_cell",
    "12": "plasma_cell",
    "13": "myeloid",
    "14": "myeloid",
}

adata.obs["broad_celltype"] = adata.obs["leiden"].map(cluster2celltype).astype("category")
adata_clean = adata[adata.obs["broad_celltype"] != "Doublet"].copy()
adata_clean.obs["broad_celltype"] = adata_clean.obs["broad_celltype"].cat.remove_unused_categories()

for view, suffix, label in iter_plot_views(adata_clean):
    sc.pl.umap(view, color="broad_celltype", title="Broad Cell Type Annotation" + label, size=Config.DOT_SIZE, show=False)
    save_figure(f"final_broad_annotation_umap{suffix}")
    plt.show()

save_path = Config.OUT_DIR / Config.CHECKPOINT_NAME
adata_clean.write(save_path)
adata = adata_clean
"""


MAIN_POST_ANNOTATION = """
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc

PATIENT_COL = "patient_id"
SAMPLE_COL = "sample_id"
TISSUE_COL = "tissue_type"
PNI_COL = "pni_level"

fig_out = Config.FIG_DIR / "post_annotation"
fig_out.mkdir(parents=True, exist_ok=True)

custom_markers = {
    "T_NK": ["CD3D", "CD3E", "TRBC1", "TRBC2", "IL7R", "NKG7"],
    "B_cell": ["MS4A1", "CD79A", "CD74", "HLA-DRA", "CD37"],
    "plasma_cell": ["MZB1", "JCHAIN", "DERL3", "IGHG1"],
    "myeloid": ["LYZ", "S100A8", "C1QC", "C1QB", "CTSS", "FCER1G"],
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "TACSTD2", "MSLN", "MUC1"],
    "endothelial": ["PECAM1", "KDR", "EMCN", "RAMP2", "ESAM"],
    "Mesenchymal/Neural-Cres": ["RGS5", "MCAM", "CSPG4", "PDGFRB", "NOTCH3"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"],
}

def canonical_label(value):
    return str(value).strip().lower().replace(" ", "_").replace("-", "_")

observed_categories = [str(cat) for cat in adata.obs["broad_celltype"].cat.categories]
observed_lookup = {canonical_label(cat): cat for cat in observed_categories}
aligned_markers = {}
for marker_label, genes in custom_markers.items():
    canon = canonical_label(marker_label)
    if canon in observed_lookup:
        observed_label = observed_lookup[canon]
        aligned_markers[observed_label] = [gene for gene in genes if gene in adata.var_names]

panel_order = list(aligned_markers.keys())
remaining_categories = [cat for cat in observed_categories if cat not in panel_order]
panel_order.extend(remaining_categories)
adata.obs["broad_celltype"] = pd.Categorical(
    adata.obs["broad_celltype"].astype(str),
    categories=panel_order,
    ordered=True,
)

for view, suffix, label in iter_plot_views(adata):
    sc.tl.dendrogram(view, groupby="broad_celltype")
    sc.pl.dotplot(view, aligned_markers, groupby="broad_celltype", standard_scale="var", color_map="Reds", dendrogram=True, title="Broad cell type marker dotplot" + label, show=False)
    save_figure(f"post_annotation/1_broad_celltype_dotplot{suffix}")
    plt.show()

    for obs_col, stem in [(PATIENT_COL, "2a_proportion_by_patient"), (SAMPLE_COL, "2b_proportion_by_sample"), (TISSUE_COL, "2c_proportion_by_tissue"), (PNI_COL, "2d_proportion_by_pni")]:
        counts = pd.crosstab(view.obs[obs_col], view.obs["broad_celltype"], normalize="index") * 100
        if obs_col == PNI_COL:
            counts = counts.reindex([level for level in ["low", "high"] if level in counts.index])
        counts.plot(kind="bar", stacked=True, figsize=(8, 6), colormap="tab20", edgecolor="white", linewidth=0.5)
        plt.title(f"{stem}{label}")
        plt.tight_layout()
        save_figure(f"post_annotation/{stem}{suffix}")
        plt.show()

    obs_df = pd.crosstab(view.obs["broad_celltype"], view.obs[PNI_COL])
    obs_df = obs_df.loc[obs_df.sum(axis=1) > 0, obs_df.sum(axis=0) > 0]
    expected_df = pd.DataFrame(np.outer(obs_df.sum(axis=1), obs_df.sum(axis=0)) / obs_df.sum().sum(), index=obs_df.index, columns=obs_df.columns)
    roe_df = obs_df / expected_df
    plt.figure(figsize=(8, 6))
    sns.heatmap(roe_df, annot=True, cmap="RdBu_r", center=1, fmt=".2f")
    plt.title("Ro/e heatmap by PNI" + label)
    plt.tight_layout()
    save_figure(f"post_annotation/3_Roe_heatmap{suffix}")
    plt.show()
"""


MAIN_POST_ANNOTATION_GROUPING_UMAPS = """
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd

TISSUE_COL = "tissue_type"
PNI_COL = "pni_level"
AT_LABELS = {"adjacent": "A", "tumor": "T"}
PNI_STATUS_LABELS = {"high": "PNI+", "low": "PNI-"}
TWO_GROUP_COLORS = ["#617fa1", "#d45d5d"]

adata.obs["tissue_group_display"] = adata.obs[TISSUE_COL].astype(str).str.lower().map(AT_LABELS).fillna(adata.obs[TISSUE_COL].astype(str))
adata.obs["tissue_group_display"] = pd.Categorical(adata.obs["tissue_group_display"], categories=["A", "T"], ordered=True)
adata.obs["pni_status_display"] = adata.obs[PNI_COL].astype(str).str.lower().map(PNI_STATUS_LABELS).fillna(adata.obs[PNI_COL].astype(str))
adata.obs["pni_status_display"] = pd.Categorical(adata.obs["pni_status_display"], categories=["PNI-", "PNI+"], ordered=True)

sc.pl.umap(
    adata,
    color="tissue_group_display",
    title="All samples: A/T",
    size=Config.DOT_SIZE,
    palette=TWO_GROUP_COLORS,
    show=False,
)
save_figure("post_annotation/0a_umap_by_AT")
plt.show()

sc.pl.umap(
    adata,
    color="pni_status_display",
    title="All samples: PNI+/-",
    size=Config.DOT_SIZE,
    palette=TWO_GROUP_COLORS,
    show=False,
)
save_figure("post_annotation/0b_umap_by_pni_status")
plt.show()

tumor_all = tumor_only_view(adata)
if tumor_all is not None and tumor_all.n_obs > 0:
    tumor_all.obs["pni_status_display"] = tumor_all.obs[PNI_COL].astype(str).str.lower().map(PNI_STATUS_LABELS).fillna(tumor_all.obs[PNI_COL].astype(str))
    tumor_all.obs["pni_status_display"] = pd.Categorical(tumor_all.obs["pni_status_display"], categories=["PNI-", "PNI+"], ordered=True)
    sc.pl.umap(
        tumor_all,
        color="pni_status_display",
        title="Tumor only: PNI+/-",
        size=Config.DOT_SIZE,
        palette=TWO_GROUP_COLORS,
        show=False,
    )
    save_figure("post_annotation/0c_umap_tumor_only_by_pni_status_tumor_only")
    plt.show()
"""


SUB_LOAD = """
checkpoint_path = Config.OUT_DIR / Config.CHECKPOINT_NAME
if not checkpoint_path.exists():
    raise FileNotFoundError("Missing annotated_broad_noP02_highlow.h5ad.")

adata = sc.read_h5ad(checkpoint_path)
required_obs = ["broad_celltype", "sample_id", "patient_id", "tissue_type", "pni_level"]
missing_obs = [col for col in required_obs if col not in adata.obs.columns]
if missing_obs:
    raise KeyError(f"Missing required obs columns: {missing_obs}")
if "counts" not in adata.layers:
    raise KeyError("Missing adata.layers['counts']")

remaining_patients = set(adata.obs["patient_id"].astype(str).unique())
remaining_pni_levels = set(adata.obs["pni_level"].astype(str).unique())
assert remaining_patients == {'P01', 'P03'}, f"Unexpected patients: {remaining_patients}"
assert remaining_pni_levels == {'high', 'low'}, f"Unexpected PNI levels: {remaining_pni_levels}"
"""


SUB_MYE_OVERVIEW = """
assert "broad_celltype" in adata.obs.columns
assert "counts" in adata.layers
mye = adata[adata.obs["broad_celltype"] == "myeloid"].copy()
mye.X = mye.layers["counts"].copy()
sc.pp.normalize_total(mye, target_sum=1e4)
sc.pp.log1p(mye)
mye = run_subset_pipeline(mye, batch_key="sample_id", resolution=0.2)
for view, suffix, label in iter_plot_views(mye):
    sc.pl.umap(view, color=["leiden", "sample_id", "patient_id", "tissue_type", "pni_level"], ncols=3, size=Config.DOT_SIZE, frameon=False, title=[f"Myeloid overview{label}: leiden", f"Myeloid overview{label}: sample_id", f"Myeloid overview{label}: patient_id", f"Myeloid overview{label}: tissue_type", f"Myeloid overview{label}: pni_level"], show=False)
    save_figure(f"subset_myeloid_umap_overview{suffix}")
    plt.show()
"""


SUB_MYE_MARKERS = """
MYELOID_MARKERS = {
    "Pan_Myeloid": ["LYZ", "TYROBP", "AIF1"],
    "Mono_Inflammatory": ["FCN1", "VCAN", "S100A8", "S100A9", "IL1B"],
    "Neutrophils": ["FCGR3B", "CSF3R", "CXCL8"],
    "C1QC_TAM": ["C1QA", "C1QB", "APOE", "TREM2", "CD163", "MRC1"],
    "SPP1_TAM": ["SPP1", "PLAUR", "CTSD"],
    "cDC2": ["CD1C", "CLEC10A", "FCER1A", "HLA-DRA", "CD74"],
    "pDC": ["LILRA4", "CLEC4C", "IRF7", "GZMB"],
    "Mast_cells": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
}
valid = {key: [gene for gene in genes if gene in mye.var_names] for key, genes in MYELOID_MARKERS.items() if any(gene in mye.var_names for gene in genes)}
plot_markers(mye, valid, title="Myeloid_Markers")
plot_patient_contribution(mye, cluster_key="leiden")
for view, suffix, label in iter_plot_views(mye):
    sc.pl.umap(view, color=["leiden", "tissue_type", "pni_level", "patient_id"], ncols=2, size=Config.DOT_SIZE, frameon=False, title=[f"Myeloid marker overview{label}: leiden", f"Myeloid marker overview{label}: tissue_type", f"Myeloid marker overview{label}: pni_level", f"Myeloid marker overview{label}: patient_id"], show=False)
    save_figure(f"Myeloid_Markers_overview{suffix}")
    plt.show()
"""


SUB_CLUSTER3 = """
from IPython.display import display

def filter_marker_table(df, top_n=25):
    exclude_patterns = (r"^RPL", r"^RPS", r"^MT-", r"^MTRNR", r"^HSP", r"^MALAT1$", r"^NEAT1$", r"^XIST$")
    mask = pd.Series(True, index=df.index)
    for pattern in exclude_patterns:
        mask &= ~df["names"].str.contains(pattern, regex=True, na=False)
    out = df.loc[mask].copy()
    keep_cols = [c for c in ["names", "logfoldchanges", "pvals_adj", "pct_nz_group"] if c in out.columns]
    return out[keep_cols].head(top_n)

mye.obs["leiden"] = mye.obs["leiden"].astype(str).astype("category")
if "3" in list(mye.obs["leiden"].cat.categories):
    sc.tl.rank_genes_groups(mye, groupby="leiden", groups=["3"], reference="rest", method="wilcoxon", pts=True)
    cluster3_vs_rest_filt = filter_marker_table(sc.get.rank_genes_groups_df(mye, group="3"), top_n=30)
    cluster3_vs_rest_filt.to_csv(Config.TAB_DIR / "myeloid_cluster3_vs_rest_markers_noP02_highlow.tsv", sep="\\t", index=False)
    contamination_valid = {
        "Lymphoid_check": [gene for gene in ["IL7R", "LTB", "CD3D", "TRAC", "CD247", "ITK", "BCL11B", "CD96"] if gene in mye.var_names],
        "DC_check": [gene for gene in ["CD1C", "CLEC10A", "FCER1A", "HLA-DRA", "CD74"] if gene in mye.var_names],
        "Myeloid_core": [gene for gene in ["LYZ", "TYROBP", "AIF1", "FCER1G"] if gene in mye.var_names],
    }
    contamination_valid = {k: v for k, v in contamination_valid.items() if v}
    for view, suffix, label in iter_plot_views(mye):
        sc.pl.dotplot(view, contamination_valid, groupby="leiden", standard_scale="var", color_map="Reds", dendrogram=False, title=f"Cluster3_Contamination_Check{label}", show=False)
        save_figure(f"Cluster3_Contamination_Check_dotplot{suffix}")
        plt.show()
        sc.pl.stacked_violin(view, contamination_valid, groupby="leiden", swap_axes=True, title=f"Cluster3_Contamination_Check{label}", show=False)
        save_figure(f"Cluster3_Contamination_Check_violin{suffix}")
        plt.show()
"""


SUB_CLUSTER12 = """
mye.obs["leiden"] = mye.obs["leiden"].astype(str).astype("category")
available_clusters = list(mye.obs["leiden"].cat.categories)
if "1" in available_clusters and "2" in available_clusters:
    sc.tl.rank_genes_groups(mye, groupby="leiden", groups=["1"], reference="2", method="wilcoxon", pts=True)
    filter_marker_table(sc.get.rank_genes_groups_df(mye, group="1"), top_n=30).to_csv(Config.TAB_DIR / "myeloid_cluster1_vs_2_markers_noP02_highlow.tsv", sep="\\t", index=False)
    sc.tl.rank_genes_groups(mye, groupby="leiden", groups=["2"], reference="1", method="wilcoxon", pts=True)
    filter_marker_table(sc.get.rank_genes_groups_df(mye, group="2"), top_n=30).to_csv(Config.TAB_DIR / "myeloid_cluster2_vs_1_markers_noP02_highlow.tsv", sep="\\t", index=False)
    refined_valid = {
        "Cluster1_C1QC_TAM_like": [gene for gene in ["C1QA", "C1QB", "C1QC", "APOE", "TREM2", "CD163", "A2M", "GPNMB", "AXL", "APOC1"] if gene in mye.var_names],
        "Cluster2_Activated_Inflammatory_Myeloid_like": [gene for gene in ["FCN1", "VCAN", "S100A8", "S100A9", "EREG", "AREG", "PLAUR", "TNFAIP6", "TREM1", "AQP9"] if gene in mye.var_names],
    }
    refined_valid = {k: v for k, v in refined_valid.items() if v}
    for view, suffix, label in iter_plot_views(mye):
        sc.pl.dotplot(view, refined_valid, groupby="leiden", standard_scale="var", color_map="Reds", dendrogram=False, title=f"Myeloid_Cluster1_vs_2_Refined_Markers{label}", show=False)
        save_figure(f"Myeloid_Cluster1_vs_2_Refined_Markers_dotplot{suffix}")
        plt.show()
        sc.pl.stacked_violin(view, refined_valid, groupby="leiden", swap_axes=True, title=f"Myeloid_Cluster1_vs_2_Refined_Markers{label}", show=False)
        save_figure(f"Myeloid_Cluster1_vs_2_Refined_Markers_violin{suffix}")
        plt.show()
"""


SUB_FINAL_MYE = """
for view, suffix, label in iter_plot_views(mye):
    sc.pl.umap(view, color=["cell_subtype", "patient_id", "tissue_type", "pni_level"], ncols=2, size=Config.DOT_SIZE, frameon=False, title=[f"Myeloid subtype overview{label}: cell_subtype", f"Myeloid subtype overview{label}: patient_id", f"Myeloid subtype overview{label}: tissue_type", f"Myeloid subtype overview{label}: pni_level"], show=False)
    save_figure(f"Myeloid_Subtype_overview{suffix}")
    plt.show()

FINAL_MYLOID_MARKERS = {
    "Neutrophils_MDSC": ["S100A8", "S100A9", "IL1B", "FCGR3B", "CSF3R", "CXCL8"],
    "FCN1_Monocytes": ["FCN1", "VCAN", "EREG", "AREG"],
    "C1QC_TAMs": ["C1QA", "C1QB", "C1QC", "APOE", "TREM2", "CD163"],
    "pDC": ["LILRA4", "CLEC4C", "IRF7", "GZMB"],
    "Mast_cells": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
}
final_valid = {k: [gene for gene in genes if gene in mye.var_names] for k, genes in FINAL_MYLOID_MARKERS.items() if any(gene in mye.var_names for gene in genes)}
for view, suffix, label in iter_plot_views(mye):
    sc.pl.dotplot(view, final_valid, groupby="cell_subtype", categories_order=subtype_order, standard_scale="var", color_map="Reds", dendrogram=False, title=f"Myeloid_Subtype_Markers_Final{label}", show=False)
    save_figure(f"Myeloid_Subtype_Markers_Final_dotplot{suffix}")
    plt.show()
    sc.pl.stacked_violin(view, final_valid, groupby="cell_subtype", categories_order=subtype_order, swap_axes=True, standard_scale="var", title=f"Myeloid_Subtype_Markers_Final{label}", show=False)
    save_figure(f"Myeloid_Subtype_Markers_Final_violin{suffix}")
    plt.show()
    rep_markers = [genes[0] for genes in final_valid.values() if genes]
    if rep_markers:
        sc.pl.umap(view, color=rep_markers, ncols=4, size=Config.DOT_SIZE, frameon=False, title=[f"Myeloid_Subtype_Markers_Final{label}: {marker}" for marker in rep_markers], show=False)
        save_figure(f"Myeloid_Subtype_Markers_Final_feature_umap{suffix}")
        plt.show()

plot_patient_contribution(mye, cluster_key="cell_subtype", filename_prefix="patient_contribution_myeloid_cell_subtype")
out_path = Config.OUT_DIR / "subset_myeloid_annotated_noP02_highlow.h5ad"
mye.write_h5ad(out_path)
"""


SUB_TARGET = """
target_genes = ["KRT23", "SERPINB3", "CARD14", "KYNU", "CA12", "CMBL", "KRT17", "BFSP2", "TESC", "TBC1D2", "MSLN"]
valid_target_genes = [gene for gene in target_genes if gene in mye.var_names]
if valid_target_genes:
    for view, suffix, label in iter_plot_views(mye):
        sc.pl.umap(view, color=valid_target_genes, cmap="Reds", size=50, title=[f"{gene} Expression{label}" for gene in valid_target_genes], show=False)
        save_figure(f"myeloid_target_gene_umap{suffix}")
        plt.show()
"""


SUB_EPI_OVERVIEW = """
epi = adata[adata.obs["broad_celltype"] == "epithelial"].copy()
epi.X = epi.layers["counts"].copy()
sc.pp.normalize_total(epi, target_sum=1e4)
sc.pp.log1p(epi)
epi = run_subset_pipeline(epi, batch_key="sample_id", resolution=0.4)
for view, suffix, label in iter_plot_views(epi):
    sc.pl.umap(view, color=["leiden", "sample_id", "patient_id", "tissue_type", "pni_level"], ncols=3, size=Config.DOT_SIZE, frameon=False, title=[f"Epithelial overview{label}: leiden", f"Epithelial overview{label}: sample_id", f"Epithelial overview{label}: patient_id", f"Epithelial overview{label}: tissue_type", f"Epithelial overview{label}: pni_level"], show=False)
    save_figure(f"subset_epithelial_umap_overview{suffix}")
    plt.show()
"""


SUB_EPI_MARKERS = """
EPI_MARKERS = {
    "Epithelial_core": ["EPCAM", "TACSTD2", "KRT8", "KRT18", "KRT19"],
    "Biliary_identity": ["KRT19", "KRT7", "SOX9", "MUC1"],
    "Tumor_like": ["MSLN", "MUC1", "EPCAM", "TACSTD2", "LGALS3"],
    "EMT_invasion": ["VIM", "ITGA6", "ITGA3", "MMP7", "FN1"],
    "Proliferation": ["MKI67", "TOP2A", "HMGB2", "TYMS"],
    "PNI_hint": ["NGFR", "L1CAM", "NCAM1", "S100B", "SOX10"],
}
valid = {k: [gene for gene in genes if gene in epi.var_names] for k, genes in EPI_MARKERS.items() if any(gene in epi.var_names for gene in genes)}
plot_markers(epi, valid, title="Epithelial_Markers")
plot_patient_contribution(epi, cluster_key="leiden")
for view, suffix, label in iter_plot_views(epi):
    sc.pl.umap(view, color=["leiden", "tissue_type", "pni_level", "patient_id"], ncols=2, size=Config.DOT_SIZE, frameon=False, title=[f"Epithelial marker overview{label}: leiden", f"Epithelial marker overview{label}: tissue_type", f"Epithelial marker overview{label}: pni_level", f"Epithelial marker overview{label}: patient_id"], show=False)
    save_figure(f"Epithelial_Markers_overview{suffix}")
    plt.show()
"""


SUB_EPI_FINAL = """
if not cluster_annotation:
    raise ValueError("cluster_annotation is empty. Fill it first in Section 5C.")

epi = apply_annotation(epi, cluster_annotation, col_name="cell_subtype")
for view, suffix, label in iter_plot_views(epi):
    sc.pl.umap(view, color=["cell_subtype", "patient_id", "tissue_type", "pni_level"], ncols=2, size=Config.DOT_SIZE, frameon=False, title=[f"Epithelial subtype overview{label}: cell_subtype", f"Epithelial subtype overview{label}: patient_id", f"Epithelial subtype overview{label}: tissue_type", f"Epithelial subtype overview{label}: pni_level"], show=False)
    save_figure(f"Epithelial_Subtype_overview{suffix}")
    plt.show()

out_path = Config.OUT_DIR / "subset_epithelial_annotated_noP02_highlow.h5ad"
epi.write_h5ad(out_path)
"""


def build_main() -> None:
    src = load_notebook(NOTEBOOKS / "02_integration_annotation_v3.ipynb")
    update_markdown_title(src["cells"][0], " (no P02 high/low)")
    src["cells"].insert(
        16,
        {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [],
        },
    )
    src["cells"][2]["source"] = to_source(SHARED_SETUP)
    src["cells"][4]["source"] = to_source(MAIN_LOAD)
    src["cells"][8]["source"] = to_source(MAIN_DIAGNOSTICS)
    src["cells"][10]["source"] = to_source(MAIN_LINEAGE_MARKERS)
    src["cells"][12]["source"] = to_source(MAIN_SCHWANN)
    src["cells"][14]["source"] = to_source(MAIN_BROAD_ANNOTATION)
    src["cells"][16]["source"] = to_source(MAIN_POST_ANNOTATION_GROUPING_UMAPS)
    src["cells"][17]["source"] = to_source(MAIN_POST_ANNOTATION)
    save_notebook(NOTEBOOKS / "02_integration_annotation_noP02_highlow.ipynb", src)


def build_subcluster() -> None:
    src = load_notebook(NOTEBOOKS / "02.2_subclustering_annotation.ipynb")
    update_markdown_title(src["cells"][0], " (no P02 high/low)")
    src["cells"][2]["source"] = to_source(SHARED_SETUP)
    src["cells"][4]["source"] = to_source(SUB_LOAD)
    src["cells"][7]["source"] = to_source(SUB_MYE_OVERVIEW)
    src["cells"][9]["source"] = to_source(SUB_MYE_MARKERS)
    src["cells"][11]["source"] = to_source(SUB_CLUSTER3)
    src["cells"][12]["source"] = to_source(SUB_CLUSTER12)
    src["cells"][17]["source"] = to_source(SUB_FINAL_MYE)
    src["cells"][18]["source"] = to_source(SUB_TARGET)
    src["cells"][21]["source"] = to_source(SUB_EPI_OVERVIEW)
    src["cells"][23]["source"] = to_source(SUB_EPI_MARKERS)
    src["cells"][27]["source"] = to_source(SUB_EPI_FINAL)
    save_notebook(NOTEBOOKS / "02.2_subclustering_annotation_noP02_highlow.ipynb", src)


if __name__ == "__main__":
    build_main()
    build_subcluster()
    print("No-P02 high/low notebooks generated.")
