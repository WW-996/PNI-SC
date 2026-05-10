import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS_DIR = ROOT / "notebooks"
MAIN_NOTEBOOK = NOTEBOOKS_DIR / "02_integration_annotation_noP02_highlow.ipynb"
SUBCLUSTER_NOTEBOOK = NOTEBOOKS_DIR / "02.2_subclustering_annotation_noP02_highlow.ipynb"
MESENCHYMAL_NEURAL_NOTEBOOK = NOTEBOOKS_DIR / "02.4_mesenchymal_neural_subclustering.ipynb"
MESENCHYMAL_NEURAL_NOTEBOOK_V2 = NOTEBOOKS_DIR / "02.4_mesenchymal_neural_subclustering_v2.ipynb"


def load_notebook(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def joined_source(path: Path) -> str:
    notebook = load_notebook(path)
    return "\n".join("".join(cell.get("source", [])) for cell in notebook["cells"])


class NoP02HighLowNotebookTests(unittest.TestCase):
    def test_main_variant_exists(self) -> None:
        self.assertTrue(MAIN_NOTEBOOK.exists(), f"Missing notebook: {MAIN_NOTEBOOK}")

    def test_subcluster_variant_exists(self) -> None:
        self.assertTrue(
            SUBCLUSTER_NOTEBOOK.exists(),
            f"Missing notebook: {SUBCLUSTER_NOTEBOOK}",
        )

    def test_mesenchymal_neural_variant_exists(self) -> None:
        self.assertTrue(
            MESENCHYMAL_NEURAL_NOTEBOOK.exists(),
            f"Missing notebook: {MESENCHYMAL_NEURAL_NOTEBOOK}",
        )

    def test_mesenchymal_neural_v2_variant_exists(self) -> None:
        self.assertTrue(
            MESENCHYMAL_NEURAL_NOTEBOOK_V2.exists(),
            f"Missing notebook: {MESENCHYMAL_NEURAL_NOTEBOOK_V2}",
        )

    def test_main_variant_uses_no_p02_outputs_and_filtering(self) -> None:
        joined = joined_source(MAIN_NOTEBOOK)

        self.assertIn('EXCLUDED_PATIENT_IDS = ["P02"]', joined)
        self.assertIn('EXPORT_PDF = True', joined)
        self.assertIn('integration_noP02_highlow', joined)
        self.assertIn("annotated_broad_noP02_highlow.h5ad", joined)
        self.assertIn("df_index = df_index[~df_index['patient_id'].isin(Config.EXCLUDED_PATIENT_IDS)].copy()", joined)
        self.assertIn("remaining_patients == {'P01', 'P03'}", joined)
        self.assertIn("remaining_pni_levels == {'high', 'low'}", joined)
        self.assertIn("tumor_only", joined)
        self.assertIn("save_figure(", joined)
        self.assertIn('TUMOR_ONLY_PLOT_LABEL = " (Tumor only)"', joined)
        self.assertIn('stem_path.with_suffix(".pdf")', joined)
        self.assertIn('(out_dir / stem_path.parent).mkdir(parents=True, exist_ok=True)', joined)
        self.assertIn('title=f"{title} - Dotplot{label}"', joined)
        self.assertIn('plt.title(f"{group_key} contribution to {cluster_key}{label}")', joined)
        self.assertIn('title=["Leiden Clusters" + label, "Cell Cycle Phase" + label]', joined)
        self.assertIn('title="Broad Cell Type Annotation" + label', joined)
        self.assertIn('title="Broad cell type marker dotplot" + label', joined)
        self.assertIn("def canonical_label(value):", joined)
        self.assertIn('aligned_markers[observed_label] = [gene for gene in genes if gene in adata.var_names]', joined)
        self.assertIn('panel_order = list(aligned_markers.keys())', joined)
        self.assertIn('adata.obs[\"broad_celltype\"] = pd.Categorical(', joined)
        self.assertIn('categories=panel_order,', joined)

    def test_main_variant_has_post_annotation_grouping_umaps(self) -> None:
        joined = joined_source(MAIN_NOTEBOOK)

        self.assertIn('PNI_STATUS_LABELS = {"high": "PNI+", "low": "PNI-"}', joined)
        self.assertIn('AT_LABELS = {"adjacent": "A", "tumor": "T"}', joined)
        self.assertIn('adata.obs["tissue_group_display"] = adata.obs[TISSUE_COL].astype(str).str.lower().map(AT_LABELS).fillna(adata.obs[TISSUE_COL].astype(str))', joined)
        self.assertIn('adata.obs["pni_status_display"] = adata.obs[PNI_COL].astype(str).str.lower().map(PNI_STATUS_LABELS).fillna(adata.obs[PNI_COL].astype(str))', joined)
        self.assertIn('title="All samples: A/T"', joined)
        self.assertIn('title="All samples: PNI+/-"', joined)
        self.assertIn('title="Tumor only: PNI+/-"', joined)
        self.assertIn('save_figure("post_annotation/0a_umap_by_AT")', joined)
        self.assertIn('save_figure("post_annotation/0b_umap_by_pni_status")', joined)
        self.assertIn('save_figure("post_annotation/0c_umap_tumor_only_by_pni_status_tumor_only")', joined)

    def test_subcluster_variant_uses_new_checkpoint_and_tumor_only_plotting(self) -> None:
        joined = joined_source(SUBCLUSTER_NOTEBOOK)

        self.assertIn('EXCLUDED_PATIENT_IDS = ["P02"]', joined)
        self.assertIn('EXPORT_PDF = True', joined)
        self.assertIn('integration_noP02_highlow', joined)
        self.assertIn("annotated_broad_noP02_highlow.h5ad", joined)
        self.assertIn("tumor_only", joined)
        self.assertIn("save_figure(", joined)
        self.assertIn("plot_markers(", joined)
        self.assertIn("plot_patient_contribution(", joined)
        self.assertIn('TUMOR_ONLY_PLOT_LABEL = " (Tumor only)"', joined)
        self.assertIn('stem_path.with_suffix(".pdf")', joined)
        self.assertIn('title=f"{title} - Dotplot{label}"', joined)
        self.assertIn('title=f"Cluster3_Contamination_Check{label}"', joined)
        self.assertIn('title=f"Myeloid_Cluster1_vs_2_Refined_Markers{label}"', joined)

    def test_mesenchymal_neural_variant_uses_expected_no_p02_flow(self) -> None:
        joined = joined_source(MESENCHYMAL_NEURAL_NOTEBOOK)

        self.assertIn('EXCLUDED_PATIENT_IDS = ["P02"]', joined)
        self.assertIn('EXPORT_PDF = True', joined)
        self.assertIn('integration_noP02_highlow', joined)
        self.assertIn("annotated_broad_noP02_highlow.h5ad", joined)
        self.assertIn('mnc = adata[adata.obs["broad_celltype"] == "Mesenchymal/Neural-Cres"].copy()', joined)
        self.assertIn('subset_mesenchymal_neural_annotated_noP02_highlow.h5ad', joined)
        self.assertIn('mesenchymal_neural_cluster_annotation_template_noP02_highlow.tsv', joined)
        self.assertIn('mesenchymal_neural_patient_contribution_noP02_highlow.tsv', joined)
        self.assertIn('mesenchymal_neural_cluster_markers_noP02_highlow.tsv', joined)
        self.assertIn('mesenchymal_neural_subtype_summary_noP02_highlow.tsv', joined)
        self.assertIn('tumor_only', joined)
        self.assertIn('save_figure(', joined)
        self.assertIn("plot_markers(", joined)
        self.assertIn("plot_patient_contribution(", joined)
        self.assertIn('Pericyte_like', joined)
        self.assertIn('Schwann_core', joined)
        self.assertIn('Glial_like', joined)
        self.assertIn('PNI_hint', joined)
        self.assertIn('mesenchymal_pericyte_like', joined)
        self.assertIn('schwann_core', joined)
        self.assertIn('glial_like', joined)
        self.assertIn('undetermined_mesenchymal_neural', joined)
        self.assertIn('sc.tl.rank_genes_groups(mnc, groupby="leiden", method="wilcoxon", pts=True)', joined)
        self.assertIn('SUSPECT_CLUSTER_ID = "8"', joined)
        self.assertIn('DROP_SUSPECT_CLUSTER = False', joined)
        self.assertIn('if Config.SUSPECT_CLUSTER_ID in available_clusters:', joined)
        self.assertIn('sc.tl.rank_genes_groups(mnc, groupby="leiden", groups=[Config.SUSPECT_CLUSTER_ID], reference="rest", method="wilcoxon", pts=True)', joined)
        self.assertIn('mesenchymal_neural_cluster8_vs_rest_markers_noP02_highlow.tsv', joined)
        self.assertIn('Cluster8_Identity_Check', joined)
        self.assertIn('if Config.DROP_SUSPECT_CLUSTER:', joined)
        self.assertIn('mnc = mnc[mnc.obs["leiden"].astype(str) != Config.SUSPECT_CLUSTER_ID].copy()', joined)
        self.assertIn('RECLUSTER_AFTER_DROP = False', joined)
        self.assertIn('mnc_pre_drop = mnc.copy()', joined)
        self.assertIn('reclustered = mnc_pre_drop.copy()', joined)
        self.assertIn('reclustered = reclustered[reclustered.obs["leiden"].astype(str) != Config.SUSPECT_CLUSTER_ID].copy()', joined)
        self.assertIn('reclustered = run_subset_pipeline(reclustered, batch_key="sample_id", resolution=Config.RESOLUTION_MAIN)', joined)
        self.assertIn('mesenchymal_neural_reclustered_cluster_markers_noP02_highlow.tsv', joined)
        self.assertIn('Mesenchymal_Neural_Reclustered_overview', joined)
        self.assertIn('Mesenchymal_Neural_Reclustered_Markers_Final', joined)
        self.assertIn('Mesenchymal_Neural_Reclustered_Markers_Final_dotplot', joined)
        self.assertIn('Mesenchymal_Neural_Reclustered_Markers_Final_violin', joined)
        self.assertIn('Mesenchymal_Neural_Reclustered_Markers_Final_feature_umap', joined)
        self.assertIn('patient_contribution_mesenchymal_neural_reclustered', joined)
        self.assertIn('subset_mesenchymal_neural_reclustered_noP02_highlow.h5ad', joined)

    def test_mesenchymal_neural_v2_removes_cluster8_drop_flow(self) -> None:
        joined = joined_source(MESENCHYMAL_NEURAL_NOTEBOOK_V2)

        self.assertIn('EXCLUDED_PATIENT_IDS = ["P02"]', joined)
        self.assertIn('EXPORT_PDF = True', joined)
        self.assertIn('integration_noP02_highlow', joined)
        self.assertIn("annotated_broad_noP02_highlow.h5ad", joined)
        self.assertIn('mnc = adata[adata.obs["broad_celltype"] == "Mesenchymal/Neural-Cres"].copy()', joined)
        self.assertIn('mesenchymal_neural_patient_contribution_noP02_highlow.tsv', joined)
        self.assertIn('mesenchymal_neural_recluster_cluster_markers_noP02_highlow.tsv', joined)
        self.assertIn('mesenchymal_neural_recluster_annotation_template_noP02_highlow.tsv', joined)
        self.assertIn('mesenchymal_neural_recluster_subtype_summary_noP02_highlow.tsv', joined)
        self.assertIn('mesenchymal_neural_recluster_cell_annotations_noP02_highlow.tsv', joined)
        self.assertIn('subset_mesenchymal_neural_recluster_final_annotated_noP02_highlow.h5ad', joined)
        self.assertIn('tumor_only', joined)
        self.assertIn('save_figure(', joined)
        self.assertIn("plot_markers(", joined)
        self.assertIn("plot_patient_contribution(", joined)
        self.assertIn('SMC_Myofibroblast', joined)
        self.assertIn('Stressed_Cell', joined)
        self.assertIn('Myeloid_Cell', joined)
        self.assertIn('T_NK_Cell', joined)
        self.assertIn('Plasma_Cell', joined)
        self.assertIn('Structural_Fibroblast', joined)
        self.assertIn('Pericyte', joined)
        self.assertIn('Matrix_Fibroblast_mCAF', joined)
        self.assertIn('Schwann_Cell', joined)
        self.assertIn('Immuno_Fibroblast_iCAF', joined)
        self.assertIn('undetermined_mesenchymal_neural', joined)
        self.assertIn('sc.tl.rank_genes_groups(mnc_recluster, groupby="leiden", method="wilcoxon", pts=True)', joined)
        self.assertIn('mnc_recluster.obs["mesenchymal_neural_subtype_recluster"]', joined)

        self.assertNotIn('SUSPECT_CLUSTER_ID = "8"', joined)
        self.assertNotIn('DROP_SUSPECT_CLUSTER = True', joined)
        self.assertNotIn('RECLUSTER_AFTER_DROP = True', joined)
        self.assertNotIn('Section 7A-2', joined)
        self.assertNotIn('Cluster 8 身份核查与可选删除接口', joined)
        self.assertNotIn('if Config.SUSPECT_CLUSTER_ID in available_clusters:', joined)
        self.assertNotIn('cluster8_vs_rest', joined)
        self.assertNotIn('mesenchymal_neural_cluster8_vs_rest_markers_noP02_highlow.tsv', joined)
        self.assertNotIn('Cluster8_Identity_Check', joined)
        self.assertNotIn('if Config.DROP_SUSPECT_CLUSTER:', joined)
        self.assertNotIn('mnc = mnc[mnc.obs["leiden"].astype(str) != Config.SUSPECT_CLUSTER_ID].copy()', joined)
        self.assertNotIn('mnc_pre_drop = mnc.copy()', joined)
        self.assertNotIn('Section 7E', joined)
        self.assertNotIn('reclustered = mnc_pre_drop.copy()', joined)
        self.assertNotIn('reclustered = run_subset_pipeline(reclustered, batch_key="sample_id", resolution=Config.RESOLUTION_MAIN)', joined)
        self.assertNotIn('Mesenchymal_Neural_Reclustered_overview', joined)
        self.assertNotIn('Mesenchymal_Neural_Reclustered_Markers_Final', joined)
        self.assertNotIn('patient_contribution_mesenchymal_neural_reclustered', joined)
        self.assertNotIn('mesenchymal_neural_reclustered_cluster_markers_noP02_highlow.tsv', joined)
        self.assertNotIn('subset_mesenchymal_neural_reclustered_noP02_highlow.h5ad', joined)

    def test_mesenchymal_neural_v2_uses_updated_exploration_markers(self) -> None:
        joined = joined_source(MESENCHYMAL_NEURAL_NOTEBOOK_V2)
        section_7b = joined.split("## Section 7B — Marker 图谱 + 患者贡献（人工判读）", 1)[1]
        section_7b = section_7b.split("## Section 7B-2 — 差异表达辅助判读", 1)[0]

        self.assertIn('"SMC_Myofibroblast": [gene for gene in ["ACTA2", "TAGLN", "MYL9", "FLNA"] if gene in mnc.var_names]', section_7b)
        self.assertIn('"Stressed_Cell": [gene for gene in ["JUNB", "FOS", "EGR1", "JUND", "NR4A1"] if gene in mnc.var_names]', section_7b)
        self.assertIn('"Myeloid_Cell": [gene for gene in ["LYZ", "CD74", "HLA-DRA", "HLA-DRB1", "LCN2"] if gene in mnc.var_names]', section_7b)
        self.assertIn('"T_NK_Cell": [gene for gene in ["PTPRC", "FYN", "CXCR4", "SRGN"] if gene in mnc.var_names]', section_7b)
        self.assertIn('"Plasma_Cell": [gene for gene in ["IGHG1", "IGLC2", "IGHG3", "IGHG4", "IGLC3"] if gene in mnc.var_names]', section_7b)
        self.assertIn('"Structural_Fibroblast": [gene for gene in ["LUM", "DCN", "COL1A2", "COL6A3", "LAMA2"] if gene in mnc.var_names]', section_7b)
        self.assertIn('"Pericyte": [gene for gene in ["RGS5", "NDUFA4L2", "HIGD1B", "SPARC", "COL18A1"] if gene in mnc.var_names]', section_7b)
        self.assertIn('"Matrix_Fibroblast_mCAF": [gene for gene in ["POSTN", "FN1", "VCAN", "BGN", "TIMP1"] if gene in mnc.var_names]', section_7b)
        self.assertIn('"Schwann_Cell": [gene for gene in ["S100B", "PLP1", "PMP22", "NRXN1", "CDH19"] if gene in mnc.var_names]', section_7b)
        self.assertIn('"Immuno_Fibroblast_iCAF": [gene for gene in ["CCL19", "CCL21", "C1S", "CFH", "SERPING1"] if gene in mnc.var_names]', section_7b)

        self.assertNotIn('"vSMC / mural": [gene for gene in ["ACTA2", "TAGLN", "MYH11", "CNN1"] if gene in mnc.var_names]', section_7b)
        self.assertNotIn('"Fibroblast/perivascular fibroblast": [gene for gene in ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"] if gene in mnc.var_names]', section_7b)
        self.assertNotIn('"Schwann/glial-like": [gene for gene in ["SOX10", "S100B", "PLP1", "MPZ", "PMP22", "NGFR"] if gene in mnc.var_names]', section_7b)
        self.assertNotIn('"neural crest progenitor-like": [gene for gene in ["FOXD3", "TFAP2A", "PAX3", "EDNRB", "NES"] if gene in mnc.var_names]', section_7b)

    def test_mesenchymal_neural_v2_exports_recluster_markers_after_section_7b(self) -> None:
        joined = joined_source(MESENCHYMAL_NEURAL_NOTEBOOK_V2)
        section_7b2 = joined.split("## Section 7B-2 — 差异表达辅助判读", 1)[1]
        section_7b2 = section_7b2.split("## Section 7C — 手动填写 cluster→亚型映射（保守版）", 1)[0]

        self.assertIn('mnc_recluster.obs["leiden"] = mnc_recluster.obs["leiden"].astype(str).astype("category")', section_7b2)
        self.assertIn('review_clusters = sorted(mnc_recluster.obs["leiden"].cat.categories, key=int)', section_7b2)
        self.assertIn('recluster_marker_top_dict = {}', section_7b2)
        self.assertIn('recluster_marker_top_dict[cluster_id] = df', section_7b2)
        self.assertIn('display(recluster_marker_top_dict[cluster_id])', section_7b2)
        self.assertIn('mesenchymal_neural_recluster_cluster_markers_noP02_highlow.tsv', section_7b2)

    def test_mesenchymal_neural_v2_uses_updated_cluster_annotations(self) -> None:
        joined = joined_source(MESENCHYMAL_NEURAL_NOTEBOOK_V2)
        section_7c = joined.split("## Section 7C — 手动填写 cluster→亚型映射（保守版）", 1)[1]
        section_7c = section_7c.split("## Section 7D — 应用注释、可视化与结果导出", 1)[0]

        self.assertIn("make_annotation_template(mnc_recluster, 'leiden')", section_7c)
        self.assertIn('"cluster_id": sorted(mnc_recluster.obs["leiden"].astype(str).unique(), key=int)', section_7c)
        self.assertIn('mesenchymal_neural_recluster_annotation_template_noP02_highlow.tsv', section_7c)
        self.assertIn('"0": "SMC_Myofibroblast"', section_7c)
        self.assertIn('"1": "SMC_Myofibroblast"', section_7c)
        self.assertIn('"2": "Pericyte"', section_7c)
        self.assertIn('"3": "SMC_Myofibroblast"', section_7c)
        self.assertIn('"4": "Matrix_Fibroblast_mCAF"', section_7c)
        self.assertIn('"5": "Immuno_Fibroblast_iCAF"', section_7c)
        self.assertIn('"6": "Structural_Fibroblast"', section_7c)
        self.assertIn('"7": "Schwann_Cell"', section_7c)
        self.assertIn('mnc_recluster = apply_annotation(mnc_recluster, recluster_annotation, col_name="mesenchymal_neural_subtype_recluster")', section_7c)
        self.assertIn('allowed_labels": "SMC_Myofibroblast; Structural_Fibroblast; Pericyte; Matrix_Fibroblast_mCAF; Schwann_Cell; Immuno_Fibroblast_iCAF; undetermined_mesenchymal_neural"', section_7c)

    def test_mesenchymal_neural_v2_uses_recluster_marker_analysis_mainline(self) -> None:
        joined = joined_source(MESENCHYMAL_NEURAL_NOTEBOOK_V2)
        section_7b2 = joined.split("## Section 7B-2 — 差异表达辅助判读", 1)[1]
        section_7b2 = section_7b2.split("## Section 7C — 手动填写 cluster→亚型映射（保守版）", 1)[0]

        self.assertIn('RECLUSTER_REVIEW_MARKERS = {', section_7b2)
        self.assertIn('"SMC_Myofibroblast_recluster"', section_7b2)
        self.assertIn('"Structural_Fibroblast_recluster"', section_7b2)
        self.assertIn('"Pericyte_recluster"', section_7b2)
        self.assertIn('"Matrix_Fibroblast_mCAF_recluster"', section_7b2)
        self.assertIn('"Schwann_Cell_recluster"', section_7b2)
        self.assertIn('"Immuno_Fibroblast_iCAF_recluster"', section_7b2)
        self.assertIn('plot_markers(mnc_recluster, recluster_review_markers_valid, title="Mesenchymal_Neural_Recluster_Clusters_Review")', section_7b2)
        self.assertNotIn('NATURAL_CLUSTER_MARKERS = {', section_7b2)
        self.assertNotIn('Mesenchymal_Neural_Natural_Clusters_0_to_9', section_7b2)

    def test_mesenchymal_neural_v2_uses_recluster_outputs_in_final_annotation_panel(self) -> None:
        joined = joined_source(MESENCHYMAL_NEURAL_NOTEBOOK_V2)
        section_7d = joined.split("## Section 7D — 应用注释、可视化与结果导出", 1)[1]

        self.assertIn('for view, suffix, label in iter_plot_views(mnc_recluster):', section_7d)
        self.assertIn('color=["leiden", "mesenchymal_neural_subtype_recluster", "patient_id", "pni_level"]', section_7d)
        self.assertIn('save_figure(f"Mesenchymal_Neural_Recluster_Subtype_overview{suffix}")', section_7d)
        self.assertIn('FINAL_RECLUSTER_MARKERS = {', section_7d)
        self.assertIn('"SMC_Myofibroblast": [gene for gene in ["ACTA2", "TAGLN", "MYL9", "FLNA"] if gene in mnc_recluster.var_names]', section_7d)
        self.assertIn('"Structural_Fibroblast": [gene for gene in ["LUM", "DCN", "COL1A2"] if gene in mnc_recluster.var_names]', section_7d)
        self.assertIn('"Pericyte": [gene for gene in ["RGS5", "NDUFA4L2", "HIGD1B", "COL18A1"] if gene in mnc_recluster.var_names]', section_7d)
        self.assertIn('"Matrix_Fibroblast_mCAF": [gene for gene in ["POSTN", "FN1", "VCAN"] if gene in mnc_recluster.var_names]', section_7d)
        self.assertIn('"Schwann_Cell": [gene for gene in ["S100B", "PLP1", "PMP22", "NRXN1", "CDH19"] if gene in mnc_recluster.var_names]', section_7d)
        self.assertIn('"Immuno_Fibroblast_iCAF": [gene for gene in ["CCL19", "CCL21", "C1S", "CFH", "SERPING1"] if gene in mnc_recluster.var_names]', section_7d)
        self.assertIn('plot_patient_contribution(mnc_recluster, cluster_key="mesenchymal_neural_subtype_recluster", filename_prefix="patient_contribution_mesenchymal_neural_recluster_subtype")', section_7d)
        self.assertIn('subset_mesenchymal_neural_recluster_final_annotated_noP02_highlow.h5ad', section_7d)
        self.assertIn('mesenchymal_neural_recluster_cell_annotations_noP02_highlow.tsv', section_7d)
        self.assertIn('mesenchymal_neural_recluster_subtype_summary_noP02_highlow.tsv', section_7d)

    def test_mesenchymal_neural_v2_adds_post_filter_reclustering_branch(self) -> None:
        joined = joined_source(MESENCHYMAL_NEURAL_NOTEBOOK_V2)

        self.assertIn("## Section 7B-1 — 去除 Cluster 2/3/4 后重新聚类", joined)
        self.assertIn('clusters_to_remove = {"2", "3", "4"}', joined)
        self.assertIn('clusters_to_keep = ["0", "1", "5", "6", "7", "8", "9"]', joined)
        self.assertIn('mnc_recluster = mnc[mnc.obs["leiden"].astype(str).isin(clusters_to_keep)].copy()', joined)
        self.assertIn('print("Removed original clusters:", sorted(clusters_to_remove, key=int))', joined)
        self.assertIn('sc.tl.tsne(mnc_recluster, use_rep="X_pca_harmony", random_state=Config.RANDOM_SEED)', joined)
        self.assertIn('mnc_recluster = run_subset_pipeline(mnc_recluster, batch_key="sample_id", resolution=Config.RESOLUTION_MAIN)', joined)
        self.assertIn('save_figure(f"subset_mesenchymal_neural_recluster_umap_overview{suffix}")', joined)
        self.assertIn('save_figure(f"subset_mesenchymal_neural_recluster_tsne_overview{suffix}")', joined)
        self.assertIn('recluster_marker_df.to_csv(Config.TAB_DIR / "mesenchymal_neural_recluster_cluster_markers_noP02_highlow.tsv"', joined)
        self.assertIn('recluster_template_df.to_csv(Config.TAB_DIR / "mesenchymal_neural_recluster_annotation_template_noP02_highlow.tsv"', joined)
        self.assertIn('recluster_annotation = {', joined)
        self.assertIn('"0": "SMC_Myofibroblast"', joined)
        self.assertIn('"1": "SMC_Myofibroblast"', joined)
        self.assertIn('"2": "Pericyte"', joined)
        self.assertIn('"3": "SMC_Myofibroblast"', joined)
        self.assertIn('"4": "Matrix_Fibroblast_mCAF"', joined)
        self.assertIn('"5": "Immuno_Fibroblast_iCAF"', joined)
        self.assertIn('"6": "Structural_Fibroblast"', joined)
        self.assertIn('"7": "Schwann_Cell"', joined)
        self.assertIn('mnc_recluster = apply_annotation(mnc_recluster, recluster_annotation, col_name="mesenchymal_neural_subtype_recluster")', joined)
        self.assertIn('present_recluster_subtypes = [label for label in recluster_subtype_order if label in set(mnc_recluster.obs["mesenchymal_neural_subtype_recluster"].astype(str))]', joined)
        self.assertIn('save_figure(f"subset_mesenchymal_neural_recluster_annotated_umap_overview{suffix}")', joined)
        self.assertIn('Mesenchymal_Neural_Recluster_Subtype_Markers_Final', joined)
        self.assertIn('plot_patient_contribution(mnc_recluster, cluster_key="mesenchymal_neural_subtype_recluster", filename_prefix="patient_contribution_mesenchymal_neural_recluster_subtype")', joined)
        self.assertIn('recluster_subtype_summary.to_csv(Config.TAB_DIR / "mesenchymal_neural_recluster_subtype_summary_noP02_highlow.tsv"', joined)
        self.assertIn('subset_mesenchymal_neural_recluster_annotated_noP02_highlow.h5ad', joined)

    def test_mesenchymal_neural_v2_recluster_branch_uses_trimmed_marker_panels(self) -> None:
        joined = joined_source(MESENCHYMAL_NEURAL_NOTEBOOK_V2)
        section_7b1 = joined.split("## Section 7B-1 — 去除 Cluster 2/3/4 后重新聚类", 1)[1]
        section_7b1 = section_7b1.split("## Section 7B-2 — 差异表达辅助判读", 1)[0]

        self.assertIn('"Structural_Fibroblast": [gene for gene in ["LUM", "DCN", "COL1A2"] if gene in mnc_recluster.var_names]', section_7b1)
        self.assertIn('"Pericyte": [gene for gene in ["RGS5", "NDUFA4L2", "HIGD1B", "COL18A1"] if gene in mnc_recluster.var_names]', section_7b1)
        self.assertIn('"Matrix_Fibroblast_mCAF": [gene for gene in ["POSTN", "FN1", "VCAN"] if gene in mnc_recluster.var_names]', section_7b1)

        self.assertNotIn('"Stressed_Cell": [gene for gene in ["JUNB", "FOS", "EGR1", "JUND", "NR4A1"] if gene in mnc_recluster.var_names]', section_7b1)
        self.assertNotIn('COL6A3', section_7b1)
        self.assertNotIn('LAMA2', section_7b1)
        self.assertNotIn('SPARC', section_7b1)
        self.assertNotIn('BGN', section_7b1)
        self.assertNotIn('TIMP1', section_7b1)


if __name__ == "__main__":
    unittest.main()
