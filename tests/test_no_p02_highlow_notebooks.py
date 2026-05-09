import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS_DIR = ROOT / "notebooks"
MAIN_NOTEBOOK = NOTEBOOKS_DIR / "02_integration_annotation_noP02_highlow.ipynb"
SUBCLUSTER_NOTEBOOK = NOTEBOOKS_DIR / "02.2_subclustering_annotation_noP02_highlow.ipynb"
MESENCHYMAL_NEURAL_NOTEBOOK = NOTEBOOKS_DIR / "02.4_mesenchymal_neural_subclustering.ipynb"


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
        self.assertIn('mnc.obs["mesenchymal_neural_subtype"]', joined)
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
        self.assertIn('DROP_SUSPECT_CLUSTER = True', joined)
        self.assertIn('if Config.SUSPECT_CLUSTER_ID in available_clusters:', joined)
        self.assertIn('sc.tl.rank_genes_groups(mnc, groupby="leiden", groups=[Config.SUSPECT_CLUSTER_ID], reference="rest", method="wilcoxon", pts=True)', joined)
        self.assertIn('mesenchymal_neural_cluster8_vs_rest_markers_noP02_highlow.tsv', joined)
        self.assertIn('Cluster8_Identity_Check', joined)
        self.assertIn('if Config.DROP_SUSPECT_CLUSTER:', joined)
        self.assertIn('mnc = mnc[mnc.obs["leiden"].astype(str) != Config.SUSPECT_CLUSTER_ID].copy()', joined)
        self.assertIn('RECLUSTER_AFTER_DROP = True', joined)
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


if __name__ == "__main__":
    unittest.main()
