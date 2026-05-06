"""
m07_decision_report — Decision Report Generation
====================================================
Generate consolidated validation evidence reports for screening hits.
Multi-criteria pose selection, classification, and mol2 export.

Modules:
    decision_report      07a — Pose selection + decision report (HTML + CSV + mol2)
    residue_comparison   07b — Residue-by-residue comparison vs reference (HTML + CSV)
    temporal_analysis    07e — Wrapper around mmpbsa_analysis (per-replica + consolidation)
"""
from .temporal_analysis import run_temporal_analysis

__version__ = "2.1.0"
