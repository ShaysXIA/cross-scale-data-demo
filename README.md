# Case 4.5 Cross-Scale Data Demo

![Case 4.5 SciForge demo](demo.gif)

This repository contains the local SciForge Case 4.5 demo artifacts:

- `demo.gif`: GitHub-previewable demo recording.
- `demo - HD 1080p.mov`: source video used to generate the GIF.
- `data/final/case45_cross_scale_integrated_dataset.jsonl`: final integrated guide-level cross-scale dataset.
- `data/processed/case45_l0_l5/case45_p0_20260708T1015Z/`: raw-derived processed cross-scale data layer.
- `run-output/datasets/`: constructed review dataset.
- `run-output/results/`: source inventory, quality audit, context expansion analysis, and final checklist.
- `run-output/figures/`: generated analysis figures.
- `run-output/reports/`: final demo report.
- `run-output/scripts/`: reproducible stage scripts.

## Final Integrated Dataset

The complete integrated result dataset is stored at:

`data/final/case45_cross_scale_integrated_dataset.jsonl`

Additional formats and metadata:

- `data/final/case45_cross_scale_integrated_dataset.json`: same records as a JSON object.
- `data/final/case45_cross_scale_integrated_dataset.tsv`: flattened table format for spreadsheet review.
- `data/final/case45_cross_scale_integrated_dataset_schema.json`: field-level schema summary.
- `data/final/case45_cross_scale_integrated_dataset_manifest.json`: file checksums and coverage summary.

This final dataset has one row per perturbation guide and integrates:

- L0 perturbation metadata.
- L1 UniProt protein annotation.
- L2 pathway annotation.
- L3 processed guide-vs-control RNA response.
- L4 processed guide-vs-control ADT/protein response.
- L5 DepMap endpoint label and score.
- selected cell counts and provenance links back to `data/processed/`.

Coverage in the integrated file: L0 35/35, L1 27/35, L2 27/35, L3 measured 27/35, L4 measured 35/35, L5 27/35.

## Raw-Derived Processed Layer

The processed data layer is stored at:

`data/processed/case45_l0_l5/case45_p0_20260708T1015Z/`

Key files:

- `processed_manifest.json`: inventory, checksums, generation metadata, and known limitations.
- `raw_matrix_inventory.json`: raw matrix source inventory.
- `metadata_cells.tsv`: parsed cell-level metadata.
- `guide_cell_assignment.tsv`: cell barcode to guide and target assignment.
- `selected_cells.tsv`: selected cells used for the demo run.
- `selected_adt_counts.tsv`: selected ADT count matrix subset.
- `selected_cdna_counts.summary.json` and `selected_cdna_counts_summary.tsv`: selected cDNA matrix summary.
- `l3_guide_rna_response.tsv` / `.json`: guide-level RNA response.
- `l4_guide_adt_response.tsv` / `.json`: guide-level ADT protein response.

`run-output/datasets/case45_review_table.*` is a downstream review table built from the processed layer and annotations; it is not the final integrated dataset or the raw-derived processed data itself.

The `recordings/` workspace is intentionally excluded from version control.
