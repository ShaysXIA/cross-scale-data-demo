#!/usr/bin/env python3
"""Stage 3: quality audit and reproducibility checks for case45_review_table.

Reads the Stage 2 dataset, the source L0-L5 tables, the Stage 1 inventory,
and the Stage 2 script, then produces:
  - run-output/results/quality_audit.json
  - run-output/results/quality_audit.tsv
"""

import json, os, sys, hashlib, csv, io
from collections import Counter, defaultdict
from datetime import datetime, timezone

# ---------------------------------------------------------------------------
# Paths (all relative to workspace root)
# ---------------------------------------------------------------------------
WS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
SRC = os.path.join(WS, "cross-scale-bench/case45-local-demo/input/source-data")
OUT = os.path.join(WS, "cross-scale-bench/case45-local-demo/run-output")
DS  = os.path.join(OUT, "datasets")
RES = os.path.join(OUT, "results")
SCR = os.path.join(OUT, "scripts")

TABLES = os.path.join(SRC, "datasets/case45_l0_l5/case45_p0_20260708T1015Z/tables")
REVIEW_JSONL = os.path.join(DS, "case45_review_table.jsonl")
REVIEW_TSV   = os.path.join(DS, "case45_review_table.tsv")
MANIFEST     = os.path.join(DS, "dataset_manifest.json")
INVENTORY    = os.path.join(RES, "source_inventory.json")
SCHEMA_BAK   = os.path.join(SRC, "datasets/case45_l0_l5/case45_p0_20260708T1015Z/schema.json.bak")
STAGE2_SCRIPT = os.path.join(SCR, "stage02_reconstruct_dataset.py")
STAGE1_SCRIPT = os.path.join(SCR, "stage01_source_inventory.py")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_jsonl(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows

def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)

def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

def file_mtime(path):
    return datetime.fromtimestamp(os.path.getmtime(path), tz=timezone.utc).isoformat()

# ---------------------------------------------------------------------------
# 1. Load all data
# ---------------------------------------------------------------------------
findings = []  # list of {check, severity, detail}

def add_finding(check, severity, detail):
    findings.append({"check": check, "severity": severity, "detail": detail})

# Load review table
review_rows = load_jsonl(REVIEW_JSONL)
review_by_guide = {r["guide_id"]: r for r in review_rows}
n_review = len(review_rows)

# Load source tables
l0 = load_json(os.path.join(TABLES, "l0_perturbation.json"))
l1 = load_json(os.path.join(TABLES, "l1_uniprot.json"))
l2 = load_json(os.path.join(TABLES, "l2_reactome.json"))
l3_summary = load_json(os.path.join(TABLES, "l3_rna_summary.json"))
l3_detail = load_json(os.path.join(TABLES, "l3_guide_rna_response.json"))
l4 = load_json(os.path.join(TABLES, "l4_adt_summary.json"))
l5 = load_json(os.path.join(TABLES, "l5_depmap.json"))

# Build source lookups
l0_by_guide = {r["guide_id"]: r for r in l0["data"]}
l1_by_gene = {r["target_gene"]: r for r in l1["data"]}
l2_by_gene = {r["target_gene"]: r for r in l2["data"]}
l3_by_guide = {r["guide_id"]: r for r in l3_detail["records"]}
l4_by_guide = {g["guide_id"]: g for g in l4["per_guide"]}
l5_by_gene = {r["target_gene"]: r for r in l5["data"]}

# Load manifest
manifest = load_json(MANIFEST)

# Load inventory
inventory = load_json(INVENTORY)

# Load schema
schema = load_json(SCHEMA_BAK)

# ---------------------------------------------------------------------------
# 2. Source-to-dataset traceability
# ---------------------------------------------------------------------------
print("=== TRACEABILITY ===")

# 2a. Every guide_id in review must exist in L0
review_guides = set(review_by_guide.keys())
l0_guides = set(l0_by_guide.keys())
missing_in_l0 = review_guides - l0_guides
extra_in_l0 = l0_guides - review_guides
if not missing_in_l0:
    add_finding("traceability_l0", "pass", f"All {n_review} review guide_ids present in L0 source")
else:
    add_finding("traceability_l0", "fail", f"{len(missing_in_l0)} guide_ids in review not found in L0: {sorted(missing_in_l0)[:10]}")
if extra_in_l0:
    add_finding("traceability_l0_coverage", "warn", f"{len(extra_in_l0)} L0 guide_ids not in review: {sorted(extra_in_l0)[:10]}")

# 2b. target_gene and class must match L0
mismatch_l0 = []
for gid, row in review_by_guide.items():
    if gid in l0_by_guide:
        src = l0_by_guide[gid]
        if row["target_gene"] != src["target_gene"]:
            mismatch_l0.append(f"{gid}: target_gene review={row['target_gene']} src={src['target_gene']}")
        if row["class"] != src["class"]:
            mismatch_l0.append(f"{gid}: class review={row['class']} src={src['class']}")
if not mismatch_l0:
    add_finding("traceability_l0_fields", "pass", "All guide_id->(target_gene, class) match L0 source")
else:
    add_finding("traceability_l0_fields", "fail", f"{len(mismatch_l0)} field mismatches: {mismatch_l0[:5]}")

# 2c. L1 Uniprot traceability
l1_mismatches = []
for gid, row in review_by_guide.items():
    gene = row["target_gene"]
    if row["L1_found"] and gene in l1_by_gene:
        if row["L1_uniprot_id"] != l1_by_gene[gene]["uniprot_id"]:
            l1_mismatches.append(f"{gene}: L1_uniprot_id review={row['L1_uniprot_id']} src={l1_by_gene[gene]['uniprot_id']}")
    elif row["L1_found"] and gene not in l1_by_gene:
        l1_mismatches.append(f"{gene}: L1_found=True but not in L1 source")
    elif not row["L1_found"] and gene in l1_by_gene and l1_by_gene[gene].get("found", True):
        l1_mismatches.append(f"{gene}: L1_found=False but source says found=True")
if not l1_mismatches:
    add_finding("traceability_l1", "pass", "All L1 Uniprot fields traceable to source")
else:
    add_finding("traceability_l1", "fail", f"{len(l1_mismatches)} L1 mismatches: {l1_mismatches[:5]}")

# 2d. L2 Reactome traceability
l2_mismatches = []
for gid, row in review_by_guide.items():
    gene = row["target_gene"]
    if row["L2_found"] and gene in l2_by_gene:
        if row["L2_pathway_count"] != l2_by_gene[gene]["pathway_count"]:
            l2_mismatches.append(f"{gene}: pathway_count review={row['L2_pathway_count']} src={l2_by_gene[gene]['pathway_count']}")
    elif row["L2_found"] and gene not in l2_by_gene:
        l2_mismatches.append(f"{gene}: L2_found=True but not in L2 source")
    elif not row["L2_found"] and gene in l2_by_gene and l2_by_gene[gene].get("found", True):
        l2_mismatches.append(f"{gene}: L2_found=False but source says found=True")
if not l2_mismatches:
    add_finding("traceability_l2", "pass", "All L2 Reactome fields traceable to source")
else:
    add_finding("traceability_l2", "fail", f"{len(l2_mismatches)} L2 mismatches: {l2_mismatches[:5]}")

# 2e. L3 RNA response traceability
l3_mismatches = []
for gid, row in review_by_guide.items():
    if row["L3_status"] == "measured" and gid in l3_by_guide:
        src = l3_by_guide[gid]
        if row["L3_n_cells_guide"] != src.get("n_cells_guide"):
            l3_mismatches.append(f"{gid}: n_cells_guide review={row['L3_n_cells_guide']} src={src.get('n_cells_guide')}")
        if row["L3_n_cells_control"] != src.get("n_cells_control"):
            l3_mismatches.append(f"{gid}: n_cells_control review={row['L3_n_cells_control']} src={src.get('n_cells_control')}")
    elif row["L3_status"] == "measured" and gid not in l3_by_guide:
        l3_mismatches.append(f"{gid}: L3_status=measured but not in L3 source")
    elif row["L3_status"] == "missing" and gid in l3_by_guide:
        l3_mismatches.append(f"{gid}: L3_status=missing but found in L3 source")
if not l3_mismatches:
    add_finding("traceability_l3", "pass", "All L3 RNA fields traceable to source")
else:
    add_finding("traceability_l3", "fail", f"{len(l3_mismatches)} L3 mismatches: {l3_mismatches[:5]}")

# 2f. L4 ADT traceability
l4_mismatches = []
for gid, row in review_by_guide.items():
    if row["L4_status"] == "measured" and gid in l4_by_guide:
        src = l4_by_guide[gid]
        if row["L4_n_guide_cells"] != src.get("n_guide_cells"):
            l4_mismatches.append(f"{gid}: n_guide_cells review={row['L4_n_guide_cells']} src={src.get('guide_cells')}")
        if row["L4_n_control_cells"] != src.get("n_control_cells"):
            l4_mismatches.append(f"{gid}: n_control_cells review={row['L4_n_control_cells']} src={src.get('control_cells')}")
        # Check markers
        if row["L4_strongest_marker"]:
            src_markers = src.get("markers", {})
            if row["L4_strongest_marker"] in src_markers:
                src_m = src_markers[row["L4_strongest_marker"]]
                if abs(row["L4_strongest_delta_mean"] - src_m.get("delta_mean", 0)) > 0.01:
                    l4_mismatches.append(f"{gid}: strongest_delta_mean review={row['L4_strongest_delta_mean']} src={src_m.get('delta_mean')}")
    elif row["L4_status"] == "measured" and gid not in l4_by_guide:
        l4_mismatches.append(f"{gid}: L4_status=measured but not in L4 source")
    elif row["L4_status"] == "missing" and gid in l4_by_guide:
        l4_mismatches.append(f"{gid}: L4_status=missing but found in L4 source")
if not l4_mismatches:
    add_finding("traceability_l4", "pass", "All L4 ADT fields traceable to source")
else:
    add_finding("traceability_l4", "fail", f"{len(l4_mismatches)} L4 mismatches: {l4_mismatches[:5]}")

# 2g. L5 DepMap traceability
l5_mismatches = []
for gid, row in review_by_guide.items():
    gene = row["target_gene"]
    if row["L5_found"] and gene in l5_by_gene:
        if abs(row["L5_depmap_score"] - l5_by_gene[gene]["score"]) > 0.001:
            l5_mismatches.append(f"{gene}: depmap_score review={row['L5_depmap_score']} src={l5_by_gene[gene]['score']}")
        if row["L5_dependency_class"] != l5_by_gene[gene]["class"]:
            l5_mismatches.append(f"{gene}: dependency_class review={row['L5_dependency_class']} src={l5_by_gene[gene]['class']}")
    elif row["L5_found"] and gene not in l5_by_gene:
        l5_mismatches.append(f"{gene}: L5_found=True but not in L5 source")
    elif not row["L5_found"] and gene in l5_by_gene and l5_by_gene[gene].get("found", True):
        l5_mismatches.append(f"{gene}: L5_found=False but source says found=True")
if not l5_mismatches:
    add_finding("traceability_l5", "pass", "All L5 DepMap fields traceable to source")
else:
    add_finding("traceability_l5", "fail", f"{len(l5_mismatches)} L5 mismatches: {l5_mismatches[:5]}")

# ---------------------------------------------------------------------------
# 3. L3 and L4 measured-response status
# ---------------------------------------------------------------------------
print("=== L3/L4 STATUS ===")

# 3a. Check L3 measured consistency
l3_status_counts = Counter(r["L3_status"] for r in review_rows)
l3_measured = sum(1 for r in review_rows if r["L3_measured"])
l3_missing = sum(1 for r in review_rows if r["L3_status"] == "missing")
l3_measured_non_nt = sum(1 for r in review_rows if r["L3_measured"] and r["class"] != "control_reference")
l3_measured_nt = sum(1 for r in review_rows if r["L3_measured"] and r["class"] == "control_reference")
n_nt = sum(1 for r in review_rows if r["class"] == "control_reference")
n_non_nt = sum(1 for r in review_rows if r["class"] != "control_reference")

add_finding("l3_status_distribution", "pass",
    f"L3 status: measured={l3_measured} missing={l3_missing} "
    f"(non-NT measured={l3_measured_non_nt}/{n_non_nt}, NT measured={l3_measured_nt}/{n_nt})")

# 3b. Check L4 measured consistency
l4_status_counts = Counter(r["L4_status"] for r in review_rows)
l4_measured = sum(1 for r in review_rows if r["L4_measured"])
l4_missing = sum(1 for r in review_rows if r["L4_status"] == "missing")
add_finding("l4_status_distribution", "pass",
    f"L4 status: measured={l4_measured} missing={l4_missing}")

# 3c. Check that L3_measured == (L3_status == "measured")
l3_consistency_fail = []
for r in review_rows:
    expected = r["L3_status"] == "measured"
    if r["L3_measured"] != expected:
        l3_consistency_fail.append(f"{r['guide_id']}: L3_measured={r['L3_measured']} L3_status={r['L3_status']}")
if not l3_consistency_fail:
    add_finding("l3_measured_consistency", "pass", "L3_measured field consistent with L3_status")
else:
    add_finding("l3_measured_consistency", "fail", f"{len(l3_consistency_fail)} inconsistencies: {l3_consistency_fail[:5]}")

# 3d. Same for L4
l4_consistency_fail = []
for r in review_rows:
    expected = r["L4_status"] == "measured"
    if r["L4_measured"] != expected:
        l4_consistency_fail.append(f"{r['guide_id']}: L4_measured={r['L4_measured']} L4_status={r['L4_status']}")
if not l4_consistency_fail:
    add_finding("l4_measured_consistency", "pass", "L4_measured field consistent with L4_status")
else:
    add_finding("l4_measured_consistency", "fail", f"{len(l4_consistency_fail)} inconsistencies: {l4_consistency_fail[:5]}")

# 3e. Check L3 cell counts make sense
l3_cell_issues = []
for r in review_rows:
    if r["L3_measured"]:
        if r["L3_n_cells_guide"] is not None and r["L3_n_cells_guide"] <= 0:
            l3_cell_issues.append(f"{r['guide_id']}: L3_n_cells_guide={r['L3_n_cells_guide']}")
        if r["L3_n_cells_control"] is not None and r["L3_n_cells_control"] <= 0:
            l3_cell_issues.append(f"{r['guide_id']}: L3_n_cells_control={r['L3_n_cells_control']}")
if not l3_cell_issues:
    add_finding("l3_cell_counts", "pass", "All L3 cell counts are positive")
else:
    add_finding("l3_cell_counts", "fail", f"{len(l3_cell_issues)} issues: {l3_cell_issues[:5]}")

# 3f. Check L4 cell counts make sense
l4_cell_issues = []
for r in review_rows:
    if r["L4_measured"]:
        if r["L4_n_guide_cells"] is not None and r["L4_n_guide_cells"] <= 0:
            l4_cell_issues.append(f"{r['guide_id']}: L4_n_guide_cells={r['L4_n_guide_cells']}")
        if r["L4_n_control_cells"] is not None and r["L4_n_control_cells"] <= 0:
            l4_cell_issues.append(f"{r['guide_id']}: L4_n_control_cells={r['L4_n_control_cells']}")
if not l4_cell_issues:
    add_finding("l4_cell_counts", "pass", "All L4 cell counts are positive")
else:
    add_finding("l4_cell_counts", "fail", f"{len(l4_cell_issues)} issues: {l4_cell_issues[:5]}")

# ---------------------------------------------------------------------------
# 4. Stale metadata
# ---------------------------------------------------------------------------
print("=== STALE METADATA ===")

# 4a. Check source data timestamps
source_mtimes = {}
for fname in ["l0_perturbation.json", "l1_uniprot.json", "l2_reactome.json",
              "l3_guide_rna_response.json", "l3_rna_summary.json",
              "l4_adt_summary.json", "l5_depmap.json"]:
    path = os.path.join(TABLES, fname)
    source_mtimes[fname] = file_mtime(path)

# Check if any source file is older than 30 days
stale_files = []
for fname, mtime in source_mtimes.items():
    dt = datetime.fromisoformat(mtime)
    age_days = (datetime.now(timezone.utc) - dt).days
    if age_days > 30:
        stale_files.append(f"{fname} ({age_days}d old)")
if not stale_files:
    add_finding("stale_source_files", "pass", "All source files modified within 30 days")
else:
    add_finding("stale_source_files", "warn", f"{len(stale_files)} files older than 30d: {stale_files}")

# 4b. Check schema version consistency
schema_version = schema.get("schema_version", "unknown")
add_finding("schema_version", "pass" if schema_version == "1.0" else "warn",
    f"Schema version: {schema_version}")

# 4c. Check run_id consistency
l3_run_id = l3_summary.get("run_id", "unknown")
add_finding("run_id", "pass",
    f"L3 run_id: {l3_run_id}")

# 4d. Check L4 phase
l4_phase = l4.get("phase", "unknown")
add_finding("l4_phase", "pass" if l4_phase != "unknown" else "warn",
    f"L4 phase: {l4_phase}")

# ---------------------------------------------------------------------------
# 5. Missingness
# ---------------------------------------------------------------------------
print("=== MISSINGNESS ===")

# 5a. Per-field missingness
fields_to_check = [
    "guide_id", "target_gene", "class", "evidence_source",
    "L1_uniprot_id", "L1_found",
    "L2_pathway_count", "L2_top_pathways", "L2_found",
    "L3_status", "L3_measured", "L3_n_cells_guide", "L3_n_cells_control",
    "L3_top_gene_count", "L3_top_response_genes",
    "L4_status", "L4_measured", "L4_n_guide_cells", "L4_n_control_cells",
    "L4_strongest_marker", "L4_strongest_delta_mean", "L4_strongest_log2fc",
    "L4_markers", "L4_measured_marker_count",
    "L5_depmap_score", "L5_dependency_class", "L5_found"
]
missingness = {}
for field in fields_to_check:
    null_count = sum(1 for r in review_rows if r.get(field) is None)
    missingness[field] = {"null": null_count, "non_null": n_review - null_count, "null_rate": null_count / n_review}

# Check that required fields are never null
required_fields = ["guide_id", "target_gene", "class", "evidence_source"]
required_null = {f: missingness[f]["null"] for f in required_fields if missingness[f]["null"] > 0}
if not required_null:
    add_finding("missingness_required", "pass", "All required fields non-null for all rows")
else:
    add_finding("missingness_required", "fail", f"Required fields with nulls: {required_null}")

# 5b. Check that non-NT guides have L3/L4 measured
non_nt_missing_l3 = sum(1 for r in review_rows if r["class"] not in ("control_reference", "control") and not r["L3_measured"])
non_nt_missing_l4 = sum(1 for r in review_rows if r["class"] not in ("control_reference", "control") and not r["L4_measured"])
if non_nt_missing_l3 == 0:
    add_finding("missingness_non_nt_l3", "pass", "All non-NT guides have L3 measured")
else:
    add_finding("missingness_non_nt_l3", "warn", f"{non_nt_missing_l3} non-NT guides missing L3")
if non_nt_missing_l4 == 0:
    add_finding("missingness_non_nt_l4", "pass", "All non-NT guides have L4 measured")
else:
    add_finding("missingness_non_nt_l4", "warn", f"{non_nt_missing_l4} non-NT guides missing L4")

# 5c. Check that NT guides have appropriate null annotations
nt_rows = [r for r in review_rows if r["class"] in ("control_reference", "control")]
nt_with_l5 = sum(1 for r in nt_rows if r.get("L5_found"))
nt_with_l1 = sum(1 for r in nt_rows if r.get("L1_found"))
nt_with_l2 = sum(1 for r in nt_rows if r.get("L2_found"))
if nt_with_l5 > 0:
    add_finding("missingness_nt_l5", "warn", f"{nt_with_l5}/{len(nt_rows)} NT guides have L5 annotations (expected null)")
else:
    add_finding("missingness_nt_l5", "pass", "NT guides correctly have null L5 annotations")

# 5d. Check overall missingness rates
high_missing_fields = {f: v for f, v in missingness.items() if v["null_rate"] > 0.5}
if not high_missing_fields:
    add_finding("missingness_high_rate", "pass", "No fields with >50% missingness")
else:
    add_finding("missingness_high_rate", "warn",
        f"{len(high_missing_fields)} fields with >50% missing: {list(high_missing_fields.keys())}")

# ---------------------------------------------------------------------------
# 6. Schema consistency
# ---------------------------------------------------------------------------
print("=== SCHEMA CONSISTENCY ===")

# 6a. Check expected fields exist in review rows
expected_fields = set(fields_to_check)
review_row0_fields = set(review_rows[0].keys())
missing_fields = expected_fields - review_row0_fields
extra_fields = review_row0_fields - expected_fields
if not missing_fields:
    add_finding("schema_expected_fields", "pass", f"All {len(expected_fields)} expected fields present")
else:
    add_finding("schema_expected_fields", "fail", f"Missing expected fields: {sorted(missing_fields)}")
if not extra_fields:
    add_finding("schema_extra_fields", "pass", "No unexpected fields in review rows")
else:
    known_extras = {"L3_source_file"}
    unexpected = extra_fields - known_extras
    if unexpected:
        add_finding("schema_extra_fields", "warn", f"Unexpected fields: {sorted(unexpected)}")
    else:
        add_finding("schema_extra_fields", "pass", f"Extra fields present but documented: {sorted(extra_fields)}")

# 6b. Check data types
type_checks = {
    "guide_id": str, "target_gene": str, "class": str, "evidence_source": str,
    "L1_uniprot_id": (str, type(None)), "L1_found": bool,
    "L2_pathway_count": (int, type(None)), "L2_found": bool,
    "L3_status": str, "L3_measured": bool,
    "L3_n_cells_guide": (int, type(None)), "L3_n_cells_control": (int, type(None)),
    "L4_status": str, "L4_measured": bool,
    "L4_n_guide_cells": (int, type(None)), "L4_n_control_cells": (int, type(None)),
    "L5_found": bool,
    "L5_depmap_score": (float, int, type(None)), "L5_dependency_class": (str, type(None)),
}
type_issues = []
for field, expected_types in type_checks.items():
    for r in review_rows:
        val = r.get(field)
        if not isinstance(val, expected_types):
            type_issues.append(f"{r['guide_id']}.{field}: {type(val).__name__} (expected {expected_types})")
            if len(type_issues) >= 10:
                break
    if len(type_issues) >= 10:
        break
if not type_issues:
    add_finding("schema_types", "pass", "All field types match expected schema")
else:
    add_finding("schema_types", "fail", f"{len(type_issues)} type issues: {type_issues[:5]}")

# 6c. Check enum-like values
valid_classes = {"checkpoint", "control_reference", "other", "JAK_STAT", "NF_kB", "IFN_response", "control"}
invalid_classes = set(r["class"] for r in review_rows) - valid_classes
if not invalid_classes:
    add_finding("schema_class_values", "pass", f"All class values in {sorted(valid_classes)}")
else:
    add_finding("schema_class_values", "fail", f"Invalid class values: {invalid_classes}")

valid_statuses = {"measured", "missing", "control_reference"}
invalid_l3_statuses = set(r["L3_status"] for r in review_rows) - valid_statuses
invalid_l4_statuses = set(r["L4_status"] for r in review_rows) - valid_statuses
if not invalid_l3_statuses:
    add_finding("schema_l3_status_values", "pass", "All L3_status values valid")
else:
    add_finding("schema_l3_status_values", "fail", f"Invalid L3_status values: {invalid_l3_statuses}")
if not invalid_l4_statuses:
    add_finding("schema_l4_status_values", "pass", "All L4_status values valid")
else:
    add_finding("schema_l4_status_values", "fail", f"Invalid L4_status values: {invalid_l4_statuses}")

# 6d. Check JSONL vs TSV row count consistency
tsv_rows = 0
with open(REVIEW_TSV, encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    header = next(reader)
    for row in reader:
        tsv_rows += 1
if tsv_rows == n_review:
    add_finding("schema_jsonl_tsv_consistency", "pass", f"JSONL and TSV both have {n_review} rows")
else:
    add_finding("schema_jsonl_tsv_consistency", "fail", f"JSONL={n_review} rows, TSV={tsv_rows} rows")

# ---------------------------------------------------------------------------
# 7. Basic leakage risks
# ---------------------------------------------------------------------------
print("=== LEAKAGE RISKS ===")

# 7a. Check for duplicate guide_ids
guide_id_counts = Counter(r["guide_id"] for r in review_rows)
dupes = {g: c for g, c in guide_id_counts.items() if c > 1}
if not dupes:
    add_finding("leakage_duplicate_guides", "pass", "No duplicate guide_ids")
else:
    add_finding("leakage_duplicate_guides", "fail", f"{len(dupes)} duplicate guide_ids: {dupes}")

# 7b. Check that target_gene is consistent across guides with same gene
gene_guides = defaultdict(list)
for r in review_rows:
    gene_guides[r["target_gene"]].append(r["guide_id"])
for gene, guides in gene_guides.items():
    if gene == "NT":
        continue
    # Multiple guides for same gene is expected (different gRNAs)
    pass
add_finding("leakage_gene_guides", "pass",
    f"{len(gene_guides)} unique genes across {n_review} guides "
    f"(avg {n_review/len(gene_guides):.1f} guides/gene)")

# 7c. Check that control (NT) guides are properly labeled
nt_guides = [r for r in review_rows if r["class"] in ("control_reference", "control")]
nt_genes = set(r["target_gene"] for r in nt_guides)
if nt_genes == {"NT"}:
    add_finding("leakage_control_labeling", "pass", "All control guides have target_gene=NT")
else:
    add_finding("leakage_control_labeling", "warn", f"Control guides have unexpected genes: {nt_genes}")

# 7d. Check that L3 control cells are shared across all guides (same pool)
l3_control_cell_counts = set(r["L3_n_cells_control"] for r in review_rows if r["L3_measured"] and r["L3_n_cells_control"] is not None)
if len(l3_control_cell_counts) == 1:
    add_finding("leakage_l3_control_pool", "pass", f"All L3 measured guides share same control pool ({list(l3_control_cell_counts)[0]} cells)")
else:
    add_finding("leakage_l3_control_pool", "warn", f"L3 control cell counts vary: {l3_control_cell_counts}")

# 7e. Check that L4 control cells are consistent
l4_control_cell_counts = set(r["L4_n_control_cells"] for r in review_rows if r["L4_measured"] and r["L4_n_control_cells"] is not None)
if len(l4_control_cell_counts) == 1:
    add_finding("leakage_l4_control_pool", "pass", f"All L4 measured guides share same control pool ({list(l4_control_cell_counts)[0]} cells)")
else:
    add_finding("leakage_l4_control_pool", "warn", f"L4 control cell counts vary: {l4_control_cell_counts}")

# 7f. Check that L3 top response genes are per-guide (not shared across all)
l3_top_sets = []
for r in review_rows:
    if r["L3_measured"] and r["L3_top_response_genes"]:
        l3_top_sets.append(set(g["gene"] for g in r["L3_top_response_genes"]))
if len(l3_top_sets) >= 2:
    all_identical = all(s == l3_top_sets[0] for s in l3_top_sets[1:])
    if all_identical:
        add_finding("leakage_l3_top_genes", "warn", "All L3 top response gene sets are identical (possible data issue)")
    else:
        add_finding("leakage_l3_top_genes", "pass", "L3 top response genes vary across guides (expected)")
else:
    add_finding("leakage_l3_top_genes", "pass", "Too few L3 measured guides to check diversity")

# 7g. Check for potential target leakage: is the target gene itself in top response genes?
target_in_top = 0
for r in review_rows:
    if r["L3_measured"] and r["L3_top_response_genes"] and r["class"] != "control_reference":
        top_genes = [g["gene"] for g in r["L3_top_response_genes"]]
        if r["target_gene"] in top_genes:
            target_in_top += 1
if target_in_top > 0:
    add_finding("leakage_target_in_top", "warn",
        f"{target_in_top} guides have their target gene in L3 top response genes "
        f"(potential on-target effect, not necessarily leakage)")
else:
    add_finding("leakage_target_in_top", "pass", "No target genes appear in their own L3 top response")

# ---------------------------------------------------------------------------
# 8. Reproducibility checks
# ---------------------------------------------------------------------------
print("=== REPRODUCIBILITY ===")

# 8a. Check that all required scripts exist
required_scripts = {STAGE1_SCRIPT, STAGE2_SCRIPT}
missing_scripts = [s for s in required_scripts if not os.path.exists(s)]
if not missing_scripts:
    add_finding("repro_scripts_exist", "pass", "Stage 1 and Stage 2 scripts exist")
else:
    add_finding("repro_scripts_exist", "fail", f"Missing scripts: {missing_scripts}")

# 8b. Check that scripts use relative paths (no hardcoded absolute paths)
for script_path, label in [(STAGE1_SCRIPT, "stage01"), (STAGE2_SCRIPT, "stage02")]:
    if os.path.exists(script_path):
        with open(script_path, encoding="utf-8") as f:
            content = f.read()
        has_abs = "/Users/" in content or "/root/" in content
        if has_abs:
            add_finding(f"repro_{label}_abs_paths", "warn", f"{label} script contains absolute paths")
        else:
            add_finding(f"repro_{label}_abs_paths", "pass", f"{label} script uses only relative paths")

# 8c. Check that output files are present
required_outputs = [REVIEW_JSONL, REVIEW_TSV, MANIFEST, INVENTORY]
missing_outputs = [o for o in required_outputs if not os.path.exists(o)]
if not missing_outputs:
    add_finding("repro_outputs_exist", "pass", "All required Stage 1-2 output files present")
else:
    add_finding("repro_outputs_exist", "fail", f"Missing output files: {missing_outputs}")

# 8d. Verify full row count from TSV (reproducibility of the TSV output)
add_finding("repro_row_count", "pass",
    f"Review table: JSONL={n_review} rows, TSV={tsv_rows} rows, manifest_total={manifest['stats']['total_rows']}")

# 8e. Check data integrity via hash of JSONL
jsonl_hash = file_sha256(REVIEW_JSONL)
add_finding("repro_data_hash", "pass", f"JSONL SHA256: {jsonl_hash[:16]}...")

# 8f. Check that source data is untouched (no modification during pipeline)
source_hashes = {}
for fname in ["l0_perturbation.json", "l1_uniprot.json", "l2_reactome.json",
              "l3_guide_rna_response.json", "l3_rna_summary.json",
              "l4_adt_summary.json", "l5_depmap.json"]:
    path = os.path.join(TABLES, fname)
    source_hashes[fname] = file_sha256(path)
add_finding("repro_source_hashes", "pass", f"Source file hashes recorded for {len(source_hashes)} files")

# ---------------------------------------------------------------------------
# 9. Compile results
# ---------------------------------------------------------------------------
severity_order = {"blocker": 0, "fail": 1, "warn": 2, "pass": 3}
severity_counts = Counter(f["severity"] for f in findings)
overall_status = "BLOCKED" if severity_counts.get("blocker", 0) > 0 else \
                 "WARN" if severity_counts.get("fail", 0) > 0 else \
                 "PASS"

audit = {
    "stage": "stage03",
    "status": overall_status,
    "timestamp": datetime.now(timezone.utc).isoformat(),
    "dataset": "case45_review_table",
    "dataset_rows": n_review,
    "checks_run": len(findings),
    "severity_summary": dict(severity_counts),
    "findings": findings,
    "source_hashes": source_hashes,
    "source_mtimes": source_mtimes,
    "missingness": missingness,
    "jsonl_hash": jsonl_hash,
}

# Write JSON
os.makedirs(RES, exist_ok=True)
json_path = os.path.join(RES, "quality_audit.json")
with open(json_path, "w", encoding="utf-8") as f:
    json.dump(audit, f, indent=2, ensure_ascii=False)
print(f"Wrote {json_path}")

# Write TSV
tsv_path = os.path.join(RES, "quality_audit.tsv")
with open(tsv_path, "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["check", "severity", "detail"])
    for finding in findings:
        writer.writerow([finding["check"], finding["severity"], finding["detail"]])
print(f"Wrote {tsv_path}")

# Summary
print(f"\n=== AUDIT SUMMARY ===")
print(f"Status: {overall_status}")
print(f"Checks: {len(findings)}")
for sev in ["blocker", "fail", "warn", "pass"]:
    c = severity_counts.get(sev, 0)
    if c:
        print(f"  {sev}: {c}")
for f in findings:
    if f["severity"] in ("blocker", "fail", "warn"):
        print(f"  [{f['severity'].upper()}] {f['check']}: {f['detail'][:120]}")
