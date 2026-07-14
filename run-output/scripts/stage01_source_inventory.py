#!/usr/bin/env python3
"""Stage 1 source inventory for case45-local-demo.

Reads CASE.md and all files under the source-data symlink, classifies inputs
into categories, and writes source_inventory.json and source_inventory.tsv
into run-output/results/.
"""

import hashlib
import json
import os
import sys
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[4]
SOURCE_ROOT = WORKSPACE / "cross-scale-bench" / "case45-local-demo" / "input" / "source-data"
OUTPUT_DIR = WORKSPACE / "cross-scale-bench" / "case45-local-demo" / "run-output" / "results"

CATEGORY_RULES = {
    "CASE.md": "contract",
    "atlas/": "atlas",
    "data/raw/": "raw",
    "data/external/": "external-cache",
    "data/processed/": "processed",
    "data/fixtures/": "fixture",
    "datasets/": "dataset",
    "figures/": "figure",
    "reports/": "report",
    "results/": "result",
    "scripts/": "script",
}


def classify(relpath: str) -> str:
    for prefix, category in CATEGORY_RULES.items():
        if relpath == prefix.rstrip("/") or relpath.startswith(prefix):
            return category
    return "unknown"


def file_info(filepath: Path, root: Path) -> dict:
    rel = str(filepath.relative_to(root))
    st = filepath.stat()
    sha = hashlib.sha256(filepath.read_bytes()).hexdigest()
    return {
        "relative_path": rel,
        "category": classify(rel),
        "size_bytes": st.st_size,
        "sha256": sha,
    }


def summarize_case(case_path: Path) -> dict:
    text = case_path.read_text(encoding="utf-8")
    lines = text.splitlines()
    return {
        "path": str(case_path.relative_to(WORKSPACE)),
        "lines": len(lines),
        "bytes": len(text.encode("utf-8")),
        "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "first_20_lines": lines[:20],
        "last_20_lines": lines[-20:],
    }


def main():
    if not SOURCE_ROOT.exists():
        print(f"ERROR: source root not found: {SOURCE_ROOT}", file=sys.stderr)
        sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    case_path = SOURCE_ROOT / "CASE.md"
    case_summary = summarize_case(case_path) if case_path.exists() else None

    files = []
    for fp in sorted(SOURCE_ROOT.rglob("*")):
        if fp.is_file() and ".git" not in fp.parts:
            files.append(file_info(fp, SOURCE_ROOT))

    categories = {}
    for f in files:
        cat = f["category"]
        if cat not in categories:
            categories[cat] = {"count": 0, "total_bytes": 0, "files": []}
        categories[cat]["count"] += 1
        categories[cat]["total_bytes"] += f["size_bytes"]
        categories[cat]["files"].append(f["relative_path"])

    total_files = len(files)
    total_bytes = sum(f["size_bytes"] for f in files)

    contract_summary = (
        "Case 4.5: Cross-Scale Cell Atlas from Multi-Database Integration. "
        "Anchored on Papalexi/Satija ECCITE-seq (GSE153056) with strict raw GEO route. "
        "L0-L5 layers: perturbation, target, pathway, RNA response, ADT protein response, "
        "endpoint phenotype (DepMap/Achilles THP-1). L6 excluded. "
        "P0-P6 construct the dataset; P7 performs automated exploratory analysis and visualization. "
        "External caches for L1/L2/L5 from UniProt, Reactome/GO, and DepMap/Achilles."
    )

    inventory = {
        "stage": "stage01",
        "status": "complete",
        "workspace": str(WORKSPACE),
        "source_root": str(SOURCE_ROOT),
        "case_contract": {
            "summary": contract_summary,
            "manifest": case_summary,
        },
        "inventory": {
            "total_files": total_files,
            "total_bytes": total_bytes,
            "categories": {cat: {"count": v["count"], "total_bytes": v["total_bytes"]}
                           for cat, v in sorted(categories.items())},
        },
        "files": files,
    }

    json_path = OUTPUT_DIR / "source_inventory.json"
    json_path.write_text(json.dumps(inventory, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"Wrote {json_path}")

    tsv_path = OUTPUT_DIR / "source_inventory.tsv"
    with open(tsv_path, "w", encoding="utf-8") as fh:
        fh.write("relative_path\tcategory\tsize_bytes\tsha256\n")
        for fi in files:
            fh.write(f"{fi['relative_path']}\t{fi['category']}\t{fi['size_bytes']}\t{fi['sha256']}\n")
    print(f"Wrote {tsv_path}")

    print(f"\nStage 1 source inventory complete: {total_files} files, {total_bytes:,} bytes")
    for cat, v in sorted(categories.items()):
        print(f"  {cat}: {v['count']} files, {v['total_bytes']:,} bytes")


if __name__ == "__main__":
    main()
