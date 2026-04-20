import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS_DIR = ROOT / "notebooks"
SOURCE_NOTEBOOK = NOTEBOOKS_DIR / "02_integration_annotation_v3.ipynb"
SPLIT_NOTEBOOK = NOTEBOOKS_DIR / "02.2_subclustering_annotation.ipynb"


def load_notebook(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def cell_source(cell: dict) -> str:
    return "".join(cell.get("source", []))


class NotebookSplitTests(unittest.TestCase):
    def test_split_notebook_exists(self) -> None:
        self.assertTrue(SPLIT_NOTEBOOK.exists(), f"Missing split notebook: {SPLIT_NOTEBOOK}")

    def test_source_notebook_stops_after_broad_level_sections(self) -> None:
        nb = load_notebook(SOURCE_NOTEBOOK)
        sources = [cell_source(cell) for cell in nb["cells"]]

        self.assertTrue(any("## Section 4.2" in src for src in sources))
        self.assertTrue(any("### 大类作图" in src for src in sources))
        self.assertFalse(any("## Section 5" in src for src in sources))
        self.assertFalse(any("## Section 6" in src for src in sources))

    def test_split_notebook_starts_from_subclustering(self) -> None:
        nb = load_notebook(SPLIT_NOTEBOOK)
        sources = [cell_source(cell) for cell in nb["cells"]]

        self.assertTrue(any("annotated_broad.h5ad" in src for src in sources))
        self.assertTrue(any("## Section 5" in src for src in sources))
        self.assertTrue(any("## Section 6" in src for src in sources))
        self.assertFalse(any("## Section 1: 数据载入与合并" in src for src in sources))
        self.assertFalse(any("### 大类作图" in src for src in sources))

    def test_split_notebook_checks_required_columns_and_counts_layer(self) -> None:
        nb = load_notebook(SPLIT_NOTEBOOK)
        joined = "\n".join(cell_source(cell) for cell in nb["cells"])

        for required in [
            "broad_celltype",
            "sample_id",
            "patient_id",
            "tissue_type",
            "pni_level",
            "counts",
        ]:
            self.assertIn(required, joined)

    def test_myeloid_section_uses_cluster_annotation_interface(self) -> None:
        nb = load_notebook(SPLIT_NOTEBOOK)
        joined = "\n".join(cell_source(cell) for cell in nb["cells"])

        self.assertIn("cluster_annotation = {", joined)
        self.assertNotIn("rename_dict =", joined)

    def test_myeloid_section_keeps_post_annotation_plotting(self) -> None:
        nb = load_notebook(SPLIT_NOTEBOOK)
        joined = "\n".join(cell_source(cell) for cell in nb["cells"])

        self.assertIn("groupby='cell_subtype'", joined)
        self.assertIn("Myeloid_Subtype_Markers", joined)


if __name__ == "__main__":
    unittest.main()
