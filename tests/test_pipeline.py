"""Tests for hit_validation pipeline."""
import pytest


class TestLigandPreparation:
    def test_import(self):
        from hit_validation.m00_preparation.ligand_preparation import run_ligand_preparation
        assert callable(run_ligand_preparation)

    def test_validate_mol2_import(self):
        from hit_validation.m00_preparation.ligand_preparation import validate_mol2
        assert callable(validate_mol2)


class TestBindingSiteDefinition:
    def test_import(self):
        from hit_validation.m00_preparation.binding_site_definition import run_binding_site_definition
        assert callable(run_binding_site_definition)


class TestGridGeneration:
    def test_import(self):
        from hit_validation.m01_docking.grid_generation import validate_existing_grids
        assert callable(validate_existing_grids)

    def test_missing_grids(self, tmp_path):
        from hit_validation.m01_docking.grid_generation import validate_existing_grids
        assert validate_existing_grids(str(tmp_path), "s.sph", "g.nrg", "g.bmp") is False


class TestDock6Runner:
    def test_import(self):
        from hit_validation.m01_docking.dock6_runner import run_dock6_batch
        assert callable(run_dock6_batch)

    def test_resolve_grid_prefix(self):
        from hit_validation.m01_docking.dock6_runner import resolve_grid_prefix
        result = resolve_grid_prefix("/path/to/grids", "grid.nrg")
        assert result.replace("\\", "/") == "/path/to/grids/grid"

    def test_rigid_template_no_secondary(self):
        from hit_validation.m01_docking.dock6_runner import DOCK6_RIGID_TEMPLATE
        assert "conformer_search_type" in DOCK6_RIGID_TEMPLATE
        assert "rigid" in DOCK6_RIGID_TEMPLATE
        assert "_secondary" not in DOCK6_RIGID_TEMPLATE

    def test_flex_template_no_secondary(self):
        from hit_validation.m01_docking.dock6_runner import DOCK6_FLEX_TEMPLATE
        assert "_secondary" not in DOCK6_FLEX_TEMPLATE


class TestFootprintRescore:
    def test_import(self):
        from hit_validation.m01_docking.footprint_rescore import run_footprint_rescore
        assert callable(run_footprint_rescore)

    def test_template_no_secondary(self):
        from hit_validation.m01_docking.footprint_rescore import FPS_RESCORE_TEMPLATE
        assert "footprint_similarity_score_primary" in FPS_RESCORE_TEMPLATE
        assert "_secondary" not in FPS_RESCORE_TEMPLATE
        assert "gbsa_hawkins" not in FPS_RESCORE_TEMPLATE


class TestGbsaRescore:
    def test_import(self):
        from hit_validation.m01_docking.gbsa_rescore import run_gbsa_rescore
        assert callable(run_gbsa_rescore)

    def test_template_primary(self):
        from hit_validation.m01_docking.gbsa_rescore import GBSA_RESCORE_TEMPLATE
        assert "gbsa_hawkins_score_primary" in GBSA_RESCORE_TEMPLATE
        assert "yes" in GBSA_RESCORE_TEMPLATE.split("gbsa_hawkins_score_primary")[1].split("\n")[0]
        assert "orient_ligand" in GBSA_RESCORE_TEMPLATE
        assert "grid_score_primary" in GBSA_RESCORE_TEMPLATE


class TestScoreCollector:
    def test_import(self):
        from hit_validation.m01_docking.score_collector import run_score_collection
        assert callable(run_score_collection)


class TestPlipAnalysis:
    def test_import(self):
        from hit_validation.m03_interaction_analysis.plip_interaction_analysis import run_plip_analysis
        assert callable(run_plip_analysis)

    def test_batch_import(self):
        from hit_validation.m03_interaction_analysis.plip_interaction_analysis import run_plip_batch_analysis
        assert callable(run_plip_batch_analysis)


class TestFootprintAnalysis:
    def test_import(self):
        from hit_validation.m04_dock6_analysis.footprint_analysis import run_footprint_analysis
        assert callable(run_footprint_analysis)

    def test_zone_definitions(self):
        from hit_validation.m04_dock6_analysis.footprint_analysis import ZONE_DEFINITIONS
        assert "phosphate" in ZONE_DEFINITIONS
        assert "xylose" in ZONE_DEFINITIONS
        assert "uracil" in ZONE_DEFINITIONS
        assert "catalytic" in ZONE_DEFINITIONS
        assert "ARG598" in ZONE_DEFINITIONS["phosphate"]["residues"]
        assert "TRP392" in ZONE_DEFINITIONS["xylose"]["residues"]


class TestDecisionReport:
    def test_import(self):
        from hit_validation.m07_decision_report.decision_report import run_decision_report
        assert callable(run_decision_report)

    def test_classify_molecule_absolute(self):
        from hit_validation.m07_decision_report.decision_report import classify_molecule
        # Absolute thresholds — no UDX reference needed
        result = classify_molecule(
            scores={"Grid_Score": -40.0},
            zone_coverage={"xylose": True, "phosphate": True, "uracil": False,
                           "ribose": False, "catalytic": False},
            zone_energies={"xylose": -5.0, "phosphate": -3.0, "uracil": 0,
                           "ribose": 0, "catalytic": 0},
        )
        assert "recommendation" in result
        assert "zones_covered" in result
        # xylose (-5.0) and phosphate (-3.0) both < -0.5 → 2 zones → "Strong" template
        assert result["n_zones_covered"] >= 2

    def test_classify_no_reference_context(self):
        from hit_validation.m07_decision_report.decision_report import run_decision_report
        import tempfile, os
        import pandas as pd
        # Pipeline works without reference_context (all None)
        with tempfile.TemporaryDirectory() as tmpdir:
            scores_csv = os.path.join(tmpdir, "scores.csv")
            pd.DataFrame([{"Name": "mol1", "Grid_Score": -35.0}]).to_csv(scores_csv, index=False)
            result = run_decision_report(
                scores_csv=scores_csv,
                output_dir=tmpdir,
                reference_scores_path=None,
                reference_plip_path=None,
                reference_footprint_path=None,
            )
            assert result["success"]
            assert result["n_molecules"] == 1