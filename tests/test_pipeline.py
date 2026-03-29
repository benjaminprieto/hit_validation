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

    def test_accepts_rescore_params(self):
        """run_score_collection should accept fps_rescore_dir and gbsa_rescore_dir."""
        import inspect
        from hit_validation.m01_docking.score_collector import run_score_collection
        sig = inspect.signature(run_score_collection)
        assert "fps_rescore_dir" in sig.parameters
        assert "gbsa_rescore_dir" in sig.parameters

    def test_all_poses_with_rescore_merge(self, tmp_path):
        """Test that all_poses CSV includes rescore columns when available."""
        import os
        from hit_validation.m01_docking.score_collector import run_score_collection

        # Create minimal scored mol2 for 01c (docking)
        mol_dir = tmp_path / "docking" / "mol1"
        mol_dir.mkdir(parents=True)
        docking_mol2 = mol_dir / "mol1_scored.mol2"
        docking_mol2.write_text(
            "##########                          Name:                mol1\n"
            "##########                    Grid_Score:          -35.000000\n"
            "##########               Grid_vdw_energy:          -30.000000\n"
            "##########                Grid_es_energy:           -5.000000\n"
            "\n"
            "@<TRIPOS>MOLECULE\nmol1\n 3 2 0 0 0\nSMALL\nGASTEIGER\n\n"
            "@<TRIPOS>ATOM\n"
            "      1 C1        0.0  0.0  0.0 C.3  1 LIG  0.0\n"
            "      2 C2        1.0  0.0  0.0 C.3  1 LIG  0.0\n"
            "      3 C3        2.0  0.0  0.0 C.3  1 LIG  0.0\n"
            "@<TRIPOS>BOND\n"
            "     1     1     2    1\n"
            "     2     2     3    1\n"
        )

        # Create FPS rescore mol2 for 01d
        fps_dir = tmp_path / "fps" / "mol1"
        fps_dir.mkdir(parents=True)
        fps_mol2 = fps_dir / "mol1_fps_scored.mol2"
        fps_mol2.write_text(
            "##########                          Name:                mol1\n"
            "##########                    FPS_Score:          -12.500000\n"
            "##########               FPS_vdw_energy:          -10.000000\n"
            "##########                FPS_es_energy:           -2.500000\n"
            "\n"
            "@<TRIPOS>MOLECULE\nmol1\n 3 2 0 0 0\nSMALL\nGASTEIGER\n\n"
            "@<TRIPOS>ATOM\n"
            "      1 C1        0.0  0.0  0.0 C.3  1 LIG  0.0\n"
            "      2 C2        1.0  0.0  0.0 C.3  1 LIG  0.0\n"
            "      3 C3        2.0  0.0  0.0 C.3  1 LIG  0.0\n"
            "@<TRIPOS>BOND\n"
            "     1     1     2    1\n"
            "     2     2     3    1\n"
        )

        # Create GBSA rescore mol2 for 01f
        gbsa_dir = tmp_path / "gbsa" / "mol1"
        gbsa_dir.mkdir(parents=True)
        gbsa_mol2 = gbsa_dir / "mol1_gbsa_scored.mol2"
        gbsa_mol2.write_text(
            "##########                          Name:                mol1\n"
            "##########                   GBSA_Score:          -22.000000\n"
            "\n"
            "@<TRIPOS>MOLECULE\nmol1\n 3 2 0 0 0\nSMALL\nGASTEIGER\n\n"
            "@<TRIPOS>ATOM\n"
            "      1 C1        0.0  0.0  0.0 C.3  1 LIG  0.0\n"
            "      2 C2        1.0  0.0  0.0 C.3  1 LIG  0.0\n"
            "      3 C3        2.0  0.0  0.0 C.3  1 LIG  0.0\n"
            "@<TRIPOS>BOND\n"
            "     1     1     2    1\n"
            "     2     2     3    1\n"
        )

        out_dir = tmp_path / "output"
        result = run_score_collection(
            docking_dir=str(tmp_path / "docking"),
            output_dir=str(out_dir),
            keep_all_poses=True,
            fps_rescore_dir=str(tmp_path / "fps"),
            gbsa_rescore_dir=str(tmp_path / "gbsa"),
        )
        assert result["success"]

        # Check all_poses CSV has rescore columns
        import pandas as pd
        all_poses = pd.read_csv(result["all_poses_csv"])
        assert "Grid_Score" in all_poses.columns
        assert "FPS_Score" in all_poses.columns or "FPS_vdw_energy" in all_poses.columns
        assert "GBSA_Score" in all_poses.columns
        assert len(all_poses) == 1  # 1 pose for mol1
        assert all_poses.iloc[0]["Grid_Score"] == pytest.approx(-35.0)
        assert all_poses.iloc[0]["GBSA_Score"] == pytest.approx(-22.0)


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

    def test_zone_definitions_empty_default(self):
        """ZONE_DEFINITIONS should be empty — zones come from campaign_config."""
        from hit_validation.m04_dock6_analysis.footprint_analysis import ZONE_DEFINITIONS
        assert ZONE_DEFINITIONS == {}

    def test_run_accepts_zones_param(self):
        """run_footprint_analysis should accept a zones parameter."""
        import inspect
        from hit_validation.m04_dock6_analysis.footprint_analysis import run_footprint_analysis
        sig = inspect.signature(run_footprint_analysis)
        assert "zones" in sig.parameters


class TestDecisionReport:
    def test_import(self):
        from hit_validation.m07_decision_report.decision_report import run_decision_report
        assert callable(run_decision_report)

    def test_classify_molecule_with_dynamic_zones(self):
        from hit_validation.m07_decision_report.decision_report import classify_molecule
        # Dynamic zones from campaign_config — no hardcoded zone names
        zones = {
            "zone_A": {"residues": {"ARG598", "LYS599"}, "label": "Salt bridge pocket"},
            "zone_B": {"residues": {"TRP392", "TRP495"}, "label": "Aromatic pocket"},
        }
        result = classify_molecule(
            scores={"Grid_Score": -40.0},
            zone_coverage={"zone_A": True, "zone_B": True},
            zone_energies={"zone_A": -3.0, "zone_B": -5.0},
            zones=zones,
        )
        assert "recommendation" in result
        assert "zones_covered" in result
        # Both zones < -0.5 -> 2 zones -> "Strong" template
        assert result["n_zones_covered"] >= 2
        assert "Strong" in result["recommendation"]

    def test_classify_no_zones(self):
        """Classification works with no zones defined."""
        from hit_validation.m07_decision_report.decision_report import classify_molecule
        result = classify_molecule(
            scores={"Grid_Score": -35.0},
            zone_coverage={},
            zone_energies={},
        )
        assert result["recommendation"] == "Weak candidate"
        assert result["n_zones_covered"] == 0

    def test_classify_purchase_candidate(self):
        from hit_validation.m07_decision_report.decision_report import classify_molecule
        result = classify_molecule(
            scores={"Grid_Score": -35.0},
            zone_coverage={"zone_A": True},
            zone_energies={"zone_A": -2.0},
        )
        assert result["recommendation"] == "Purchase candidate"

    def test_run_without_reference_context(self):
        from hit_validation.m07_decision_report.decision_report import run_decision_report
        import tempfile, os
        import pandas as pd
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

    def test_pose_selection(self):
        """Test multi-criteria pose selection logic."""
        from hit_validation.m07_decision_report.decision_report import select_best_pose
        import pandas as pd
        poses = pd.DataFrame([
            {"Pose_Index": 0, "Grid_Score": -30.0, "GBSA_Score": -20.0},
            {"Pose_Index": 1, "Grid_Score": -35.0, "GBSA_Score": -25.0},
            {"Pose_Index": 2, "Grid_Score": -32.0, "GBSA_Score": -28.0},
        ])
        result = select_best_pose(poses, {"Grid_Score": 0.5, "GBSA_Score": 0.5})
        assert "best_pose_idx" in result
        assert "composite_score" in result
        assert "top5" in result
        # Pose 1 has best Grid, pose 2 has best GBSA — composite determines winner
        assert result["best_pose_idx"] in [1, 2]

    def test_extract_pose_mol2(self):
        """Test single pose extraction from multi-pose mol2."""
        from hit_validation.m07_decision_report.decision_report import extract_pose_mol2
        import tempfile, os
        mol2_content = """@<TRIPOS>MOLECULE
mol_pose0
 10 10 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C1        0.0000    0.0000    0.0000 C.3     1  LIG1        0.0000
@<TRIPOS>BOND
     1     1     2    1
@<TRIPOS>MOLECULE
mol_pose1
 10 10 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C1        1.0000    1.0000    1.0000 C.3     1  LIG1        0.0000
@<TRIPOS>BOND
     1     1     2    1
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".mol2", delete=False) as f:
            f.write(mol2_content)
            tmp_path = f.name
        try:
            pose0 = extract_pose_mol2(tmp_path, 0)
            assert pose0 is not None
            assert "mol_pose0" in pose0
            pose1 = extract_pose_mol2(tmp_path, 1)
            assert pose1 is not None
            assert "mol_pose1" in pose1
        finally:
            os.unlink(tmp_path)

    def test_no_hardcoded_zone_definitions(self):
        """ZONE_DEFINITIONS should be empty — zones come from campaign_config."""
        from hit_validation.m07_decision_report.decision_report import ZONE_DEFINITIONS
        assert ZONE_DEFINITIONS == {}

    def test_no_hardcoded_key_residues(self):
        """KEY_RESIDUES should not exist as a module-level constant."""
        import hit_validation.m07_decision_report.decision_report as dr
        assert not hasattr(dr, "KEY_RESIDUES")
