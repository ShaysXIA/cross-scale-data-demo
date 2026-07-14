#!/usr/bin/env python3
"""
Stage 4: Context-expansion and causal-readiness analysis.

Reads the L0-L5 review table (case45_review_table.jsonl) and produces:
  - run-output/results/context_expansion_analysis.json
  - run-output/results/context_expansion_analysis.tsv
"""

import json
import os
import hashlib
from datetime import datetime, timezone
from collections import Counter, defaultdict

# --- Paths ---
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATASET_PATH = os.path.join(ROOT, "datasets", "case45_review_table.jsonl")
RESULTS_DIR = os.path.join(ROOT, "results")
JSON_OUT = os.path.join(RESULTS_DIR, "context_expansion_analysis.json")
TSV_OUT = os.path.join(RESULTS_DIR, "context_expansion_analysis.tsv")

os.makedirs(RESULTS_DIR, exist_ok=True)

# --- Load data ---
with open(DATASET_PATH) as f:
    rows = [json.loads(line) for line in f if line.strip()]

N = len(rows)

# --- Helper: classify rows ---
def is_control(r):
    return r["class"] == "control" or r["L3_status"] == "control_reference"

def is_perturbation(r):
    return not is_control(r)

controls = [r for r in rows if is_control(r)]
perturbations = [r for r in rows if is_perturbation(r)]

# ============================================================
# 1. Context-expansion value analysis
# ============================================================

# 1a. L1 UniProt target annotation value
l1_found = sum(1 for r in rows if r["L1_found"])
l1_unique_ids = len(set(r["L1_uniprot_id"] for r in rows if r["L1_uniprot_id"]))
l1_multi_gene_per_id = defaultdict(list)
for r in rows:
    l1_multi_gene_per_id[r["L1_uniprot_id"]].append(r["target_gene"])

l1_analysis = {
    "total_rows": N,
    "l1_found_count": l1_found,
    "l1_found_rate": round(l1_found / N, 4),
    "l1_unique_uniprot_ids": l1_unique_ids,
    "l1_value": "UniProt IDs provide reviewed protein-level annotation (accession, function, domains) "
                "that gene symbols alone cannot. This enables cross-reference with protein databases, "
                "structural biology, and drug-target resources.",
    "l1_limitation": "All targets are human-reviewed Swiss-Prot entries; isoforms and PTMs are not captured here."
}

# 1b. L2 Reactome/GO pathway value
l2_found = sum(1 for r in rows if r["L2_found"])
l2_pathway_counts = [r["L2_pathway_count"] for r in rows]
l2_zero_pathways = sum(1 for c in l2_pathway_counts if c == 0)

# Collect all unique pathways
all_pathways = set()
for r in rows:
    for p in r.get("L2_top_pathways", []):
        all_pathways.add(p)

# Per-class pathway enrichment
class_pathways = defaultdict(Counter)
for r in rows:
    for p in r.get("L2_top_pathways", []):
        class_pathways[r["class"]][p] += 1

l2_analysis = {
    "l2_found_count": l2_found,
    "l2_found_rate": round(l2_found / N, 4),
    "l2_pathway_count_min": min(l2_pathway_counts),
    "l2_pathway_count_max": max(l2_pathway_counts),
    "l2_pathway_count_mean": round(sum(l2_pathway_counts) / N, 2),
    "l2_zero_pathway_guides": l2_zero_pathways,
    "l2_unique_pathways": len(all_pathways),
    "l2_value": "Reactome/GO pathways place individual gene perturbations into a mechanistic context. "
                "They reveal whether a perturbation targets a hub in signaling cascades (e.g., JAK-STAT, "
                "NF-kB, immune checkpoint) or a peripheral node. Pathway annotations enable grouping "
                "perturbations by biological process rather than by gene family alone.",
    "l2_limitation": "Pathway annotations are curated human knowledge with varying granularity; "
                     "top-level pathways (e.g., 'Immune System') are broad and less informative. "
                     "Zero-pathway guides exist."
}

# 1c. L3 RNA response value
l3_measured = sum(1 for r in rows if r["L3_measured"] and not is_control(r))
l3_perturbation_rows = [r for r in perturbations if r["L3_measured"]]
l3_n_cells_guide = [r["L3_n_cells_guide"] for r in l3_perturbation_rows]
l3_n_cells_control = [r["L3_n_cells_control"] for r in l3_perturbation_rows]
l3_top_gene_counts = [r["L3_top_gene_count"] for r in l3_perturbation_rows]

# Collect top response genes across all perturbations
l3_all_response_genes = Counter()
for r in l3_perturbation_rows:
    for g in r.get("L3_top_response_genes", []):
        l3_all_response_genes[g["gene"]] += 1

l3_analysis = {
    "l3_measured_perturbations": l3_measured,
    "l3_control_references": sum(1 for r in rows if r["L3_status"] == "control_reference"),
    "l3_n_cells_guide_min": min(l3_n_cells_guide) if l3_n_cells_guide else None,
    "l3_n_cells_guide_max": max(l3_n_cells_guide) if l3_n_cells_guide else None,
    "l3_n_cells_guide_mean": round(sum(l3_n_cells_guide) / len(l3_n_cells_guide), 1) if l3_n_cells_guide else None,
    "l3_n_cells_control_total": l3_n_cells_control[0] if l3_n_cells_control else None,
    "l3_top_gene_count_mean": round(sum(l3_top_gene_counts) / len(l3_top_gene_counts), 1) if l3_top_gene_counts else None,
    "l3_recurrent_response_genes": l3_all_response_genes.most_common(10),
    "l3_value": "scRNA-seq provides transcriptome-wide readout of perturbation effects. "
               "It captures both direct target gene knockdown and indirect transcriptional "
               "responses. The guide vs. control comparison identifies differentially expressed "
               "genes that may mediate the perturbation-to-phenotype path.",
    "l3_limitation": "Single time point; dynamic responses missed. Cell counts vary per guide. "
                     "Top-N gene lists may miss subtle but important changes."
}

# 1d. L4 ADT protein response value
l4_measured = sum(1 for r in rows if r["L4_measured"])
l4_perturbation_rows = [r for r in perturbations if r["L4_measured"]]
l4_control_rows = [r for r in controls if r["L4_measured"]]

# ADT marker distribution
l4_marker_stats = Counter()
for r in rows:
    for marker_name in r.get("L4_markers", {}):
        l4_marker_stats[marker_name] += 1

l4_delta_means = [r["L4_strongest_delta_mean"] for r in l4_perturbation_rows
                  if r["L4_strongest_delta_mean"] is not None]
l4_log2fcs = [r["L4_strongest_log2fc"] for r in l4_perturbation_rows
              if r["L4_strongest_log2fc"] is not None]

# Per-class ADT marker shifts
l4_class_delta = defaultdict(list)
for r in perturbations:
    if r["L4_strongest_delta_mean"] is not None:
        l4_class_delta[r["class"]].append(r["L4_strongest_delta_mean"])

l4_analysis = {
    "l4_measured_count": l4_measured,
    "l4_measured_rate": round(l4_measured / N, 4),
    "l4_adt_markers_available": len(l4_marker_stats),
    "l4_adt_markers": dict(l4_marker_stats.most_common()),
    "l4_delta_mean_min": round(min(l4_delta_means), 2) if l4_delta_means else None,
    "l4_delta_mean_max": round(max(l4_delta_means), 2) if l4_delta_means else None,
    "l4_delta_mean_mean": round(sum(l4_delta_means) / len(l4_delta_means), 2) if l4_delta_means else None,
    "l4_log2fc_min": round(min(l4_log2fcs), 4) if l4_log2fcs else None,
    "l4_log2fc_max": round(max(l4_log2fcs), 4) if l4_log2fcs else None,
    "l4_log2fc_mean": round(sum(l4_log2fcs) / len(l4_log2fcs), 4) if l4_log2fcs else None,
    "l4_per_class_delta_mean": {k: round(sum(v) / len(v), 2) for k, v in l4_class_delta.items()},
    "l4_value": "ADT protein measurements capture surface protein abundance changes, "
               "providing a functional readout closer to phenotype than transcriptomics. "
               "Markers include immune checkpoint proteins (PDL1, PDL2, CD86) directly "
               "relevant to the checkpoint biology under study. Protein-level changes "
               "may diverge from RNA due to post-transcriptional regulation.",
    "l4_limitation": "Limited to 4 ADT markers (PDL1, PDL2, CD86, CD274). "
                     "Does not capture intracellular protein changes or secreted factors."
}

# 1e. L5 DepMap endpoint value
l5_found = sum(1 for r in rows if r["L5_found"])
l5_scores = [r["L5_depmap_score"] for r in rows if r["L5_depmap_score"] is not None]
l5_classes = Counter(r["L5_dependency_class"] for r in rows)

# Per-class DepMap score
l5_class_scores = defaultdict(list)
for r in rows:
    if r["L5_depmap_score"] is not None:
        l5_class_scores[r["class"]].append(r["L5_depmap_score"])

l5_analysis = {
    "l5_found_count": l5_found,
    "l5_found_rate": round(l5_found / N, 4),
    "l5_missing_count": N - l5_found,
    "l5_depmap_score_min": round(min(l5_scores), 4) if l5_scores else None,
    "l5_depmap_score_max": round(max(l5_scores), 4) if l5_scores else None,
    "l5_depmap_score_mean": round(sum(l5_scores) / len(l5_scores), 4) if l5_scores else None,
    "l5_dependency_class_distribution": dict(l5_classes),
    "l5_per_class_score_mean": {k: round(sum(v) / len(v), 4) for k, v in l5_class_scores.items()},
    "l5_value": "DepMap/Achilles gene effect scores provide a cancer-relevant endpoint phenotype "
               "in the matched THP-1 cell line. This anchors the molecular perturbation layers "
               "(L3, L4) to a disease-relevant fitness outcome. The score quantifies how essential "
               "each gene is for THP-1 proliferation/survival.",
    "l5_limitation": "8 of 35 targets missing from DepMap. Scores reflect pooled CRISPR screen "
                     "in a single cell line, not individual-level clinical outcome. "
                     "Gene effect in THP-1 may not generalize to other AML contexts."
}

# ============================================================
# 2. Cross-layer consistency analysis
# ============================================================

# 2a. L3-L4 concordance: compare direction of strongest RNA response vs ADT response
# We don't have a direct RNA equivalent of ADT markers, but we can check if L3 and L4 both show activity
l3_l4_concordance = []
for r in perturbations:
    l3_active = r["L3_measured"] and r["L3_top_gene_count"] > 0
    l4_active = r["L4_measured"] and r["L4_strongest_delta_mean"] is not None
    l3_l4_concordance.append({
        "guide_id": r["guide_id"],
        "target_gene": r["target_gene"],
        "class": r["class"],
        "l3_active": l3_active,
        "l4_active": l4_active,
        "concordant": l3_active == l4_active
    })

concordant_count = sum(1 for c in l3_l4_concordance if c["concordant"])

l3_l4_consistency = {
    "total_perturbations": len(l3_l4_concordance),
    "both_active": sum(1 for c in l3_l4_concordance if c["l3_active"] and c["l4_active"]),
    "neither_active": sum(1 for c in l3_l4_concordance if not c["l3_active"] and not c["l4_active"]),
    "l3_only": sum(1 for c in l3_l4_concordance if c["l3_active"] and not c["l4_active"]),
    "l4_only": sum(1 for c in l3_l4_concordance if not c["l3_active"] and c["l4_active"]),
    "concordant_count": concordant_count,
    "concordant_rate": round(concordant_count / len(l3_l4_concordance), 4) if l3_l4_concordance else 0,
    "interpretation": "L3 and L4 are concordant when both show measurable perturbation response. "
                      "Divergence between RNA and protein levels may reflect post-transcriptional "
                      "regulation, protein half-life differences, or ADT marker panel limitations."
}

# 2b. L0-L5 alignment: does perturbation target show dependency?
l0_l5_alignment = []
for r in perturbations:
    score = r["L5_depmap_score"]
    dep_class = r["L5_dependency_class"]
    l0_l5_alignment.append({
        "guide_id": r["guide_id"],
        "target_gene": r["target_gene"],
        "class": r["class"],
        "depmap_score": score,
        "dependency_class": dep_class,
        "has_dependency": dep_class in ("strong_dependency", "moderate_dependency") if dep_class != "missing" else None
    })

l0_l5_consistency = {
    "total_perturbations": len(l0_l5_alignment),
    "with_depmap_data": sum(1 for a in l0_l5_alignment if a["depmap_score"] is not None),
    "moderate_dependency": sum(1 for a in l0_l5_alignment if a["dependency_class"] == "moderate_dependency"),
    "weak_or_no_dependency": sum(1 for a in l0_l5_alignment if a["dependency_class"] == "weak_or_no_dependency"),
    "missing": sum(1 for a in l0_l5_alignment if a["dependency_class"] == "missing"),
    "interpretation": "Most perturbation targets show weak or no dependency in THP-1 DepMap. "
                      "This suggests either (a) the perturbation effect is context-dependent and "
                      "not captured by pooled CRISPR fitness screens, (b) the gene is not essential "
                      "for THP-1 proliferation, or (c) the ECCITE-seq perturbation (single time point) "
                      "does not translate to the long-term fitness endpoint measured by DepMap."
}

# 2c. L2-L3 pathway-to-response alignment
# For each perturbation class, check if top pathways align with observed response
l2_l3_pathway_alignment = {}
for cls in set(r["class"] for r in rows):
    cls_rows = [r for r in rows if r["class"] == cls]
    cls_pathways = Counter()
    for r in cls_rows:
        for p in r.get("L2_top_pathways", []):
            cls_pathways[p] += 1
    cls_genes = Counter()
    for r in cls_rows:
        if r["L3_measured"]:
            for g in r.get("L3_top_response_genes", []):
                cls_genes[g["gene"]] += 1
    l2_l3_pathway_alignment[cls] = {
        "top_pathways": cls_pathways.most_common(5),
        "top_response_genes": cls_genes.most_common(5),
        "pathway_count": len(cls_pathways),
        "gene_count": len(cls_genes)
    }

# ============================================================
# 3. Causal-inference quality view
# ============================================================

causal_quality = {
    "treatment_definition": {
        "layer": "L0",
        "source": "GSE153056 ECCITE-seq CRISPR guide metadata",
        "quality": "well-defined",
        "details": "Each guide targets a specific gene via CRISPR interference. "
                   "Guide identity is confirmed by sgRNA barcode. 35 guides targeting "
                   "genes in checkpoint, JAK-STAT, IFN-response, NF-kB, and control classes.",
        "concerns": [
            "Guide efficiency varies across targets; no on-target cleavage validation data included.",
            "Off-target effects are not quantified in this dataset.",
            "Single guide per target; no independent replicates with different sgRNAs."
        ]
    },
    "outcome_definition": {
        "layer": "L5",
        "source": "DepMap/Achilles CRISPR GeneEffect (THP-1, ACH-000146)",
        "quality": "moderate",
        "details": "DepMap gene effect score represents the dependency of THP-1 cell proliferation "
                   "on each gene, measured via pooled CRISPR screen. Score < -0.5 = strong dependency, "
                   "-0.5 to -0.1 = moderate, > -0.1 = weak or no dependency.",
        "concerns": [
            "8 of 35 targets (22.9%) missing from DepMap; cannot assess endpoint for these genes.",
            "Pooled CRISPR screen outcome is a population-level fitness measure, not individual-cell phenotype.",
            "Single cell line (THP-1); generalizability to other AML subtypes or primary cells is unknown.",
            "Gene effect score is a composite across multiple sgRNAs; may mask heterogeneous effects."
        ]
    },
    "mediator_layers": {
        "L3_RNA": {
            "layer": "L3",
            "type": "transcriptomic mediator",
            "measurement": "scRNA-seq differential expression (guide vs control)",
            "quality": "well-defined",
            "concerns": [
                "Single time point post-perturbation; temporal dynamics unknown.",
                "Top-12 genes per guide may miss indirect downstream mediators.",
                "Cell counts vary across guides (min 142, max varies); statistical power uneven."
            ]
        },
        "L4_ADT": {
            "layer": "L4",
            "type": "protein surface marker mediator",
            "measurement": "ADT (antibody-derived tag) protein abundance",
            "quality": "moderate",
            "concerns": [
                "Only 4 ADT markers (PDL1, PDL2, CD86, CD274); limited proteomic coverage.",
                "Surface proteins only; intracellular signaling proteins not captured.",
                "ADT signal may be affected by antibody binding competition or epitope masking."
            ]
        }
    },
    "confounding_limitations": [
        {
            "type": "batch_effects",
            "severity": "moderate",
            "detail": "ECCITE-seq data from a single experiment batch; no batch correction needed "
                      "but also no cross-batch validation. Any technical artifact in this batch "
                      "cannot be distinguished from biological signal."
        },
        {
            "type": "cell_cycle",
            "severity": "low",
            "detail": "Cell cycle state may confound RNA expression but is unlikely to be "
                      "systematically biased across guide targets in this randomized design."
        },
        {
            "type": "guide_efficiency_variation",
            "severity": "moderate",
            "detail": "CRISPRi knockdown efficiency varies by target gene and guide sequence. "
                      "Without protein-level validation of knockdown, weaker responses may "
                      "reflect poor knockdown rather than biological irrelevance."
        },
        {
            "type": "selection_bias",
            "severity": "moderate",
            "detail": "Only 35 guides across 5 functional classes were selected. This is not a "
                      "genome-wide screen. Results may not generalize to untested genes. "
                      "Class selection (checkpoint, JAK-STAT, IFN, NF-kB) biases toward "
                      "immune signaling pathways."
        },
        {
            "type": "leakage",
            "severity": "low",
            "detail": "Control guides are used for L3 differential expression normalization. "
                      "If control guides themselves have off-target effects, the reference "
                      "baseline is compromised. 8 control guides provide reasonable reference."
        },
        {
            "type": "temporal_confounding",
            "severity": "high",
            "detail": "L3/L4 measurements are from a single post-perturbation time point. "
                      "The DepMap L5 endpoint is a long-term proliferation assay (days to weeks). "
                      "The temporal gap between acute molecular response (L3/L4) and chronic "
                      "fitness endpoint (L5) is a major limitation for causal inference."
        },
        {
            "type": "cell_line_specificity",
            "severity": "high",
            "detail": "All measurements are in THP-1, a single AML cell line. "
                      "The molecular response to a perturbation may be entirely different "
                      "in primary cells, other AML subtypes, or in vivo contexts."
        }
    ],
    "causal_claim_support": {
        "supported_claims": [
            "Descriptive: Gene X perturbation in THP-1 is associated with RNA expression changes in genes A, B, C.",
            "Descriptive: Gene X perturbation alters surface protein abundance of immune checkpoint markers.",
            "Associative: Genes with stronger L4 ADT response tend to have different L5 DepMap scores.",
            "Hypothesis-generating: Targets with large L3/L4 responses and moderate L5 dependency "
            "are candidates for functional validation."
        ],
        "unsupported_claims": [
            "Causal: Perturbation of gene X causes DepMap dependency via L3/L4 molecular mediators.",
            "Mediation: The effect of L0 perturbation on L5 phenotype is mediated through L3 RNA changes.",
            "Generalizable: Findings in THP-1 apply to other AML contexts or primary patient samples.",
            "Clinical: Gene X is a therapeutic target in AML."
        ],
        "overall_assessment": "The dataset supports hypothesis generation and associative analysis. "
                             "It does NOT support causal claims due to: (1) single time point, "
                             "(2) lack of mediation analysis design, (3) temporal mismatch between "
                             "L3/L4 (acute) and L5 (chronic), (4) single cell line, (5) small N (35 guides). "
                             "The data is best used to nominate targets for follow-up validation "
                             "with orthogonal methods (e.g., time-course, independent cell lines, "
                             "in vivo models)."
    }
}

# ============================================================
# 4. Context-expansion summary
# ============================================================

# Quantify "value added" per layer
layer_value = {
    "L0_perturbation": {
        "description": "CRISPR guide identity and target gene. The anchor of the dataset.",
        "unique_contribution": "Defines the experimental perturbation; all other layers are responses to L0.",
        "rows_covered": N,
        "completeness": 1.0
    },
    "L1_target": {
        "description": "UniProt protein annotation of target genes.",
        "unique_contribution": "Links gene symbols to reviewed protein entries with functional domains, "
                              "enabling protein-level queries and cross-database integration.",
        "rows_covered": l1_found,
        "completeness": round(l1_found / N, 4)
    },
    "L2_pathway": {
        "description": "Reactome/GO pathway membership of target proteins.",
        "unique_contribution": "Places individual perturbations into biological pathway context. "
                              "Enables grouping by mechanism rather than gene family.",
        "rows_covered": l2_found,
        "completeness": round(l2_found / N, 4)
    },
    "L3_RNA": {
        "description": "scRNA-seq transcriptomic response to perturbation.",
        "unique_contribution": "Transcriptome-wide readout of perturbation effects. "
                              "Identifies downstream genes and pathways affected by each perturbation.",
        "rows_covered": l3_measured + sum(1 for r in rows if r["L3_status"] == "control_reference"),
        "completeness": round((l3_measured + sum(1 for r in rows if r["L3_status"] == "control_reference")) / N, 4)
    },
    "L4_ADT": {
        "description": "ADT protein surface marker response to perturbation.",
        "unique_contribution": "Protein-level readout of immune checkpoint markers. "
                              "Captures functional changes closer to phenotype than transcriptomics.",
        "rows_covered": l4_measured,
        "completeness": round(l4_measured / N, 4)
    },
    "L5_endpoint": {
        "description": "DepMap/Achilles gene effect score in THP-1.",
        "unique_contribution": "Disease-relevant fitness endpoint. "
                              "Anchors molecular perturbations to a cancer dependency phenotype.",
        "rows_covered": l5_found,
        "completeness": round(l5_found / N, 4)
    }
}

# ============================================================
# 5. Build output
# ============================================================

timestamp = datetime.now(timezone.utc).isoformat()

output = {
    "stage": "stage04",
    "status": "complete",
    "timestamp": timestamp,
    "dataset": "case45_review_table",
    "dataset_rows": N,
    "context_expansion_value": {
        "layer_value_added": layer_value,
        "l1_analysis": l1_analysis,
        "l2_analysis": l2_analysis,
        "l3_analysis": l3_analysis,
        "l4_analysis": l4_analysis,
        "l5_analysis": l5_analysis,
        "overall_assessment": (
            "The L0-L5 cross-scale dataset adds substantial value beyond simple layer coverage. "
            "Key contributions: (1) L1/L2 provide protein-level annotation and pathway context "
            "that transforms a list of gene perturbations into a mechanistically organized resource. "
            "(2) L3/L4 provide dual-modal (RNA + protein) molecular response data that enables "
            "comparison of transcriptomic and proteomic perturbation effects. "
            "(3) L5 anchors the molecular responses to a disease-relevant cancer dependency endpoint, "
            "enabling association between acute molecular changes and chronic fitness outcomes. "
            "(4) The cross-layer structure enables multi-scale queries: from pathway (L2) to "
            "molecular response (L3/L4) to phenotype (L5). "
            "Limitations include: single time point, single cell line, 4-plex ADT panel, "
            "22.9% missing L5 data, and small N (35 guides)."
        )
    },
    "cross_layer_consistency": {
        "l3_l4_consistency": l3_l4_consistency,
        "l0_l5_alignment": l0_l5_consistency,
        "l2_l3_pathway_alignment": l2_l3_pathway_alignment
    },
    "causal_inference_quality": causal_quality,
    "input_hashes": {
        "case45_review_table_jsonl": hashlib.sha256(
            open(DATASET_PATH, "rb").read()
        ).hexdigest()
    }
}

# Write JSON
with open(JSON_OUT, "w") as f:
    json.dump(output, f, indent=2, default=str)

# Write TSV
tsv_rows = []

# Header
tsv_rows.append([
    "section", "subsection", "metric", "value", "detail"
])

# Layer value
for layer_id, layer_data in layer_value.items():
    tsv_rows.append([
        "layer_value", layer_id, "description", layer_data["description"], ""
    ])
    tsv_rows.append([
        "layer_value", layer_id, "unique_contribution", layer_data["unique_contribution"], ""
    ])
    tsv_rows.append([
        "layer_value", layer_id, "rows_covered", str(layer_data["rows_covered"]), ""
    ])
    tsv_rows.append([
        "layer_value", layer_id, "completeness", str(layer_data["completeness"]), ""
    ])

# L1-L5 analysis summaries
for layer_key in ["l1_analysis", "l2_analysis", "l3_analysis", "l4_analysis", "l5_analysis"]:
    layer_data = output["context_expansion_value"][layer_key]
    for k, v in layer_data.items():
        if k in ("l3_recurrent_response_genes", "l4_adt_markers", "l4_per_class_delta_mean",
                 "l5_dependency_class_distribution", "l5_per_class_score_mean"):
            v = json.dumps(v)
        tsv_rows.append([
            layer_key, "", k, str(v), ""
        ])

# L3-L4 concordance
tsv_rows.append([
    "cross_layer", "l3_l4", "concordant_rate", str(l3_l4_consistency["concordant_rate"]),
    l3_l4_consistency["interpretation"]
])
tsv_rows.append([
    "cross_layer", "l3_l4", "both_active", str(l3_l4_consistency["both_active"]), ""
])
tsv_rows.append([
    "cross_layer", "l3_l4", "neither_active", str(l3_l4_consistency["neither_active"]), ""
])

# L0-L5 alignment
tsv_rows.append([
    "cross_layer", "l0_l5", "moderate_dependency", str(l0_l5_consistency["moderate_dependency"]),
    l0_l5_consistency["interpretation"]
])
tsv_rows.append([
    "cross_layer", "l0_l5", "weak_or_no_dependency", str(l0_l5_consistency["weak_or_no_dependency"]), ""
])
tsv_rows.append([
    "cross_layer", "l0_l5", "missing", str(l0_l5_consistency["missing"]), ""
])

# Causal quality
for section, data in causal_quality.items():
    if isinstance(data, dict):
        for k, v in data.items():
            if isinstance(v, (str, int, float, bool)):
                tsv_rows.append([
                    "causal_quality", section, k, str(v), ""
                ])
            elif isinstance(v, list):
                tsv_rows.append([
                    "causal_quality", section, k, str(len(v)), json.dumps(v)
                ])
            elif isinstance(v, dict):
                for k2, v2 in v.items():
                    tsv_rows.append([
                        "causal_quality", f"{section}/{k}", k2, str(v2), ""
                    ])

# Write TSV
with open(TSV_OUT, "w") as f:
    for row in tsv_rows:
        f.write("\t".join(str(c).replace("\t", " ").replace("\n", " ") for c in row) + "\n")

print(f"JSON written to: {JSON_OUT}")
print(f"TSV written to: {TSV_OUT}")
print(f"Dataset rows analyzed: {N}")
print(f"L3-L4 concordance rate: {l3_l4_consistency['concordant_rate']}")
print(f"L5 missing rate: {(N - l5_found) / N:.2%}")
print("Stage 4: COMPLETE")
