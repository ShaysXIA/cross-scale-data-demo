#!/usr/bin/env python3
"""Stage 2: reconstruct a review dataset linking L0-L5 fields from source evidence.

Reads from the existing dataset tables under source-data/datasets/ and
the unified atlas, then produces a normalised review table.
"""

import hashlib
import json
import os
import sys
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[4]
SOURCE_ROOT = WORKSPACE / "cross-scale-bench" / "case45-local-demo" / "input" / "source-data"
DATASET_DIR = SOURCE_ROOT / "datasets" / "case45_l0_l5" / "case45_p0_20260708T1015Z"
OUTPUT_DIR = WORKSPACE / "cross-scale-bench" / "case45-local-demo" / "run-output" / "datasets"


def load_json(path):
    with open(path) as f:
        return json.load(f)


def load_table_json(path):
    """Load a table JSON file that has {rows, data} or {records} structure."""
    data = load_json(path)
    if "data" in data:
        return data["data"]
    if "records" in data:
        return data["records"]
    if "per_guide" in data:
        return data["per_guide"]
    return data


def build_l3_by_guide(l3_records):
    """Index L3 records by guide_id."""
    l3_map = {}
    for r in l3_records:
        gid = r["guide_id"]
        top_genes = r.get("top_response_genes", [])
        l3_map[gid] = {
            "L3_status": r.get("status", "missing"),
            "L3_n_cells_guide": r.get("n_cells_guide"),
            "L3_n_cells_control": r.get("n_cells_control"),
            "L3_top_gene_count": r.get("top_gene_count", 0),
            "L3_top_response_genes": [
                {
                    "gene": g["gene"],
                    "log2fc": g["log2fc"],
                    "delta_mean": g["delta_mean"],
                }
                for g in top_genes
            ],
            "L3_measured": r.get("status") == "measured",
            "L3_source_file": r.get("source_cdna_file", ""),
        }
    return l3_map


def build_l4_by_guide(l4_per_guide):
    """Index L4 records by guide_id."""
    l4_map = {}
    for g in l4_per_guide:
        gid = g["guide_id"]
        markers = g.get("markers", {})
        l4_map[gid] = {
            "L4_status": g.get("status", "missing"),
            "L4_n_guide_cells": g.get("n_guide_cells"),
            "L4_n_control_cells": g.get("n_control_cells"),
            "L4_strongest_marker": g.get("strongest_marker"),
            "L4_strongest_delta_mean": g.get("strongest_delta"),
            "L4_strongest_log2fc": g.get("strongest_log2fc"),
            "L4_markers": {
                m: {
                    "guide_mean": v.get("guide_mean"),
                    "control_mean": v.get("control_mean"),
                    "delta_mean": v.get("delta_mean"),
                    "log2fc": v.get("log2fc"),
                }
                for m, v in markers.items()
            },
            "L4_measured": g.get("status") == "measured",
            "L4_measured_marker_count": len(markers),
        }
    return l4_map


def build_l1_l2_l5_maps(l1_rows, l2_rows, l5_rows):
    """Build gene-level maps for L1, L2, L5."""
    l1_map = {}
    for r in l1_rows:
        l1_map[r["target_gene"]] = {
            "L1_uniprot_id": r.get("uniprot_id"),
            "L1_found": r.get("found", False),
        }

    l2_map = {}
    for r in l2_rows:
        l2_map[r["target_gene"]] = {
            "L2_found": r.get("found", False),
            "L2_pathway_count": r.get("pathway_count", 0),
            "L2_top_pathways": r.get("top_pathways", []),
        }

    l5_map = {}
    for r in l5_rows:
        l5_map[r["target_gene"]] = {
            "L5_depmap_score": r.get("score"),
            "L5_dependency_class": r.get("class"),
            "L5_found": r.get("found", False),
        }

    return l1_map, l2_map, l5_map


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # --- Load source tables ---
    l0_rows = load_table_json(DATASET_DIR / "tables" / "l0_perturbation.json")
    l1_rows = load_table_json(DATASET_DIR / "tables" / "l1_uniprot.json")
    l2_rows = load_table_json(DATASET_DIR / "tables" / "l2_reactome.json")
    l3_records = load_table_json(DATASET_DIR / "tables" / "l3_guide_rna_response.json")
    l4_per_guide = load_table_json(DATASET_DIR / "tables" / "l4_adt_summary.json")
    l5_rows = load_table_json(DATASET_DIR / "tables" / "l5_depmap.json")

    # --- Build lookup maps ---
    l3_by_guide = build_l3_by_guide(l3_records)
    l4_by_guide = build_l4_by_guide(l4_per_guide)
    l1_by_gene, l2_by_gene, l5_by_gene = build_l1_l2_l5_maps(l1_rows, l2_rows, l5_rows)

    # --- Build review rows ---
    rows = []
    for l0 in l0_rows:
        guide_id = l0["guide_id"]
        target_gene = l0["target_gene"]

        l1 = l1_by_gene.get(target_gene, {
            "L1_uniprot_id": None,
            "L1_found": False,
        })
        l2 = l2_by_gene.get(target_gene, {
            "L2_found": False,
            "L2_pathway_count": 0,
            "L2_top_pathways": [],
        })
        l3 = l3_by_guide.get(guide_id, {
            "L3_status": "missing",
            "L3_n_cells_guide": None,
            "L3_n_cells_control": None,
            "L3_top_gene_count": 0,
            "L3_top_response_genes": [],
            "L3_measured": False,
            "L3_source_file": "",
        })
        l4 = l4_by_guide.get(guide_id, {
            "L4_status": "missing",
            "L4_n_guide_cells": None,
            "L4_n_control_cells": None,
            "L4_strongest_marker": None,
            "L4_strongest_delta_mean": None,
            "L4_strongest_log2fc": None,
            "L4_markers": {},
            "L4_measured": False,
            "L4_measured_marker_count": 0,
        })
        l5 = l5_by_gene.get(target_gene, {
            "L5_depmap_score": None,
            "L5_dependency_class": None,
            "L5_found": False,
        })

        row = {
            "guide_id": guide_id,
            "target_gene": target_gene,
            "class": l0.get("class", "unknown"),
            "evidence_source": l0.get("evidence_source", ""),
            # L1
            "L1_uniprot_id": l1["L1_uniprot_id"],
            "L1_found": l1["L1_found"],
            # L2
            "L2_found": l2["L2_found"],
            "L2_pathway_count": l2["L2_pathway_count"],
            "L2_top_pathways": l2["L2_top_pathways"],
            # L3
            "L3_status": l3["L3_status"],
            "L3_measured": l3["L3_measured"],
            "L3_n_cells_guide": l3["L3_n_cells_guide"],
            "L3_n_cells_control": l3["L3_n_cells_control"],
            "L3_top_gene_count": l3["L3_top_gene_count"],
            "L3_top_response_genes": l3["L3_top_response_genes"],
            "L3_source_file": l3["L3_source_file"],
            # L4
            "L4_status": l4["L4_status"],
            "L4_measured": l4["L4_measured"],
            "L4_n_guide_cells": l4["L4_n_guide_cells"],
            "L4_n_control_cells": l4["L4_n_control_cells"],
            "L4_strongest_marker": l4["L4_strongest_marker"],
            "L4_strongest_delta_mean": l4["L4_strongest_delta_mean"],
            "L4_strongest_log2fc": l4["L4_strongest_log2fc"],
            "L4_markers": l4["L4_markers"],
            "L4_measured_marker_count": l4["L4_measured_marker_count"],
            # L5
            "L5_depmap_score": l5["L5_depmap_score"],
            "L5_dependency_class": l5["L5_dependency_class"],
            "L5_found": l5["L5_found"],
        }
        rows.append(row)

    # --- Write JSONL ---
    jsonl_path = OUTPUT_DIR / "case45_review_table.jsonl"
    with open(jsonl_path, "w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row, ensure_ascii=False) + "\n")
    print(f"Wrote {jsonl_path}")

    # --- Write TSV ---
    tsv_path = OUTPUT_DIR / "case45_review_table.tsv"
    tsv_columns = [
        "guide_id", "target_gene", "class", "evidence_source",
        "L1_uniprot_id", "L1_found",
        "L2_found", "L2_pathway_count",
        "L3_status", "L3_measured", "L3_n_cells_guide", "L3_n_cells_control",
        "L3_top_gene_count",
        "L4_status", "L4_measured", "L4_n_guide_cells", "L4_n_control_cells",
        "L4_strongest_marker", "L4_strongest_delta_mean", "L4_strongest_log2fc",
        "L4_measured_marker_count",
        "L5_depmap_score", "L5_dependency_class", "L5_found",
    ]
    with open(tsv_path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(tsv_columns) + "\n")
        for row in rows:
            vals = []
            for col in tsv_columns:
                v = row.get(col, "")
                if isinstance(v, (list, dict)):
                    v = json.dumps(v, ensure_ascii=False)
                if v is None:
                    v = ""
                vals.append(str(v))
            fh.write("\t".join(vals) + "\n")
    print(f"Wrote {tsv_path}")

    # --- Compute coverage stats ---
    total = len(rows)
    target_genes = sorted(set(r["target_gene"] for r in rows))
    l0_count = total
    l1_found = sum(1 for r in rows if r["L1_found"])
    l2_found = sum(1 for r in rows if r["L2_found"])
    l3_measured = sum(1 for r in rows if r["L3_measured"])
    l3_missing = sum(1 for r in rows if r["L3_status"] == "missing")
    l4_measured = sum(1 for r in rows if r["L4_measured"])
    l4_missing = sum(1 for r in rows if r["L4_status"] == "missing")
    l5_found = sum(1 for r in rows if r["L5_found"])

    coverage = {
        "L0": {"count": l0_count, "coverage": 1.0},
        "L1": {"count": l1_found, "coverage": round(l1_found / total, 4) if total else 0},
        "L2": {"count": l2_found, "coverage": round(l2_found / total, 4) if total else 0},
        "L3": {"count": l3_measured, "coverage": round(l3_measured / total, 4) if total else 0,
               "missing": l3_missing, "is_measured": True},
        "L4": {"count": l4_measured, "coverage": round(l4_measured / total, 4) if total else 0,
               "missing": l4_missing, "is_measured": True},
        "L5": {"count": l5_found, "coverage": round(l5_found / total, 4) if total else 0},
    }

    # --- Write manifest ---
    manifest = {
        "dataset": "case45_review_table",
        "stage": "stage02",
        "source": "Reconstructed from L0-L5 source tables in case45_p0_20260708T1015Z",
        "files": {
            "jsonl": str(jsonl_path.relative_to(WORKSPACE)),
            "tsv": str(tsv_path.relative_to(WORKSPACE)),
        },
        "stats": {
            "total_rows": total,
            "unique_target_genes": len(target_genes),
            "target_genes": target_genes,
            "layer_coverage": coverage,
        },
        "missingness_policy": {
            "L3_missing": "Guides with status 'missing' or not found in L3 records have null L3 fields",
            "L4_missing": "Guides with status 'missing' or not found in L4 records have null L4 fields",
            "L1_L2_L5_missing": "Genes not found in annotation caches have found=False and null annotation fields",
            "control_guides": "NT guides (class=control_reference) are included with measured L3/L4 but null L5 annotations",
        },
    }
    manifest_path = OUTPUT_DIR / "dataset_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False)
    print(f"Wrote {manifest_path}")

    # --- Summary ---
    print(f"\n=== Stage 2 Review Dataset Summary ===")
    print(f"Total rows: {total}")
    print(f"Unique target genes: {len(target_genes)}")
    print(f"Target genes: {', '.join(target_genes)}")
    print(f"\nLayer coverage:")
    for layer, info in coverage.items():
        is_measured = " (measured)" if info.get("is_measured") else ""
        print(f"  {layer}: {info['count']}/{total}{is_measured}")
    print(f"\nL4 evidence: measured ({l4_measured}/{total} guides have measured ADT protein response)")
    print(f"L4 strongest markers: {', '.join(sorted(set(r['L4_strongest_marker'] for r in rows if r['L4_strongest_marker']))) or 'none'}")
    print(f"\nMissingness: L3_missing={l3_missing}, L4_missing={l4_missing}, L1_not_found={total - l1_found}, L2_not_found={total - l2_found}, L5_not_found={total - l5_found}")


if __name__ == "__main__":
    main()
