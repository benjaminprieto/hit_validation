"""
hit_validation - Screening Hit Validation Pipeline
=================================================================
Rigorous validation of virtual screening hits using DOCK6 with
footprint-based interaction analysis and GB/SA implicit solvation.

Evaluates each hit on its own absolute merits to produce evidence
for purchase decisions and pharmacophore template selection.

Modules:
    m00_preparation        - Ligand preparation (00a), receptor (00b), binding site (00d)
    m01_docking            - DOCK6 engine: grids (01b), docking (01c), footprint (01d),
                             GB/SA (01f), scores (01e)
    m03_interaction_analysis - PLIP interaction analysis (03a)
    m04_dock6_analysis     - Footprint analysis: per-residue energy consensus (04b)
    m07_decision_report    - Decision report with pose selection + mol2 export (07a)
"""
__version__ = "1.0.0"