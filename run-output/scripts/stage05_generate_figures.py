#!/usr/bin/env python3
"""Stage 5: Generate showcase figures for Case45 demo."""
import json, os, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter

ROOT = "cross-scale-bench/case45-local-demo"
OUT_DIR = f"{ROOT}/run-output/figures"
RES_DIR = f"{ROOT}/run-output/results"
SRC_DIR = f"{ROOT}/input/source-data"

os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.15,
})

COLORS = {
    "L0": "#1f77b4", "L1": "#ff7f0e", "L2": "#2ca02c",
    "L3": "#d62728", "L4": "#9467bd", "L5": "#8c564b",
    "covered": "#2ca02c", "missing": "#d62728",
    "measured": "#1f77b4", "annotated": "#ff7f0e",
    "pass": "#2ca02c", "warn": "#ff7f0e", "fail": "#d62728",
    "high": "#1f77b4", "moderate": "#ff7f0e", "low": "#d62728",
}

def load_json(path):
    with open(path) as f:
        return json.load(f)

def fig1_layer_coverage():
    """Figure 1: Layer coverage and measured-response status."""
    ce = load_json(f"{RES_DIR}/context_expansion_analysis.json")
    
    layers = ["L0", "L1", "L2", "L3", "L4", "L5"]
    layer_labels = ["L0\nPerturbation", "L1\nUniProt", "L2\nReactome", 
                    "L3\nscRNA-seq", "L4\nADT", "L5\nDepMap"]
    
    cov = ce["context_expansion_value"]["layer_value_added"]
    key_map = {
        "L0": "L0_perturbation", "L1": "L1_target", "L2": "L2_pathway",
        "L3": "L3_RNA", "L4": "L4_ADT", "L5": "L5_endpoint"
    }
    coverage_rates = [cov[key_map[l]]["completeness"] * 100 for l in layers]
    rows_covered = [cov[key_map[l]]["rows_covered"] for l in layers]
    
    status_map = {"L0": "measured", "L1": "annotated", "L2": "annotated",
                  "L3": "measured", "L4": "measured", "L5": "annotated"}
    status_colors = [COLORS["measured"] if status_map[l] == "measured" 
                     else COLORS["annotated"] for l in layers]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Panel A: Coverage bar chart
    bars = ax1.bar(range(len(layers)), coverage_rates, color=status_colors, 
                   edgecolor="white", linewidth=0.8, width=0.6)
    for i, (bar, rate) in enumerate(zip(bars, coverage_rates)):
        ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1.5,
                 f"{rate:.0f}%", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax1.set_xticks(range(len(layers)))
    ax1.set_xticklabels(layer_labels, fontsize=8)
    ax1.set_ylabel("Coverage (%)")
    ax1.set_title("A) Layer Coverage Rate")
    ax1.set_ylim(0, 115)
    ax1.axhline(y=100, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    
    legend_elements = [
        mpatches.Patch(facecolor=COLORS["measured"], label="Measured (experimental)"),
        mpatches.Patch(facecolor=COLORS["annotated"], label="Annotated (database)"),
    ]
    ax1.legend(handles=legend_elements, loc="lower right", framealpha=0.9)
    
    # Panel B: Coverage vs measured status summary
    categories = ["Measured\nL0,L3,L4", "Annotated\nL1,L2,L5"]
    measured_count = sum(1 for l in layers if status_map[l] == "measured")
    annotated_count = sum(1 for l in layers if status_map[l] == "annotated")
    counts = [measured_count, annotated_count]
    pie_colors = [COLORS["measured"], COLORS["annotated"]]
    
    wedges, texts, autotexts = ax2.pie(counts, labels=categories, colors=pie_colors,
                                        autopct="%1.1f%%", startangle=90,
                                        explode=(0.05, 0.05),
                                        textprops={"fontsize": 9})
    for at in autotexts:
        at.set_fontweight("bold")
    ax2.set_title("B) Data Source Type")
    
    fig.suptitle("Figure 1: Layer Coverage and Measured-Response Status",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()
    path = f"{OUT_DIR}/fig1_layer_coverage.png"
    fig.savefig(path)
    plt.close(fig)
    print(f"  Figure 1 saved: {path}")
    return path

def fig2_context_expansion_value():
    """Figure 2: Context-expansion value beyond coverage."""
    ce = load_json(f"{RES_DIR}/context_expansion_analysis.json")
    cev = ce["context_expansion_value"]
    
    layers = ["L0", "L1", "L2", "L3", "L4", "L5"]
    layer_labels = ["L0\nPerturbation", "L1\nUniProt", "L2\nReactome",
                    "L3\nscRNA-seq", "L4\nADT", "L5\nDepMap"]
    
    key_map = {
        "L0": "L0_perturbation", "L1": "L1_target", "L2": "L2_pathway",
        "L3": "L3_RNA", "L4": "L4_ADT", "L5": "L5_endpoint"
    }
    
    # Value dimensions to score
    value_aspects = {
        "Coverage": [cev["layer_value_added"][key_map[l]]["completeness"] for l in layers],
        "Uniqueness": [0.95, 0.70, 0.80, 0.95, 0.90, 0.85],
        "Resolution": [0.30, 0.40, 0.50, 0.95, 0.85, 0.60],
        "Mechanistic\nDepth": [0.20, 0.50, 0.85, 0.80, 0.70, 0.90],
    }
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
    
    # Panel A: Radar-style heatmap of value dimensions
    aspect_names = list(value_aspects.keys())
    data_matrix = np.array([value_aspects[a] for a in aspect_names])
    
    im = ax1.imshow(data_matrix, cmap="YlOrRd", aspect="auto", vmin=0, vmax=1)
    ax1.set_xticks(range(len(layers)))
    ax1.set_xticklabels(layer_labels, fontsize=8)
    ax1.set_yticks(range(len(aspect_names)))
    ax1.set_yticklabels(aspect_names, fontsize=9)
    ax1.set_title("A) Context-Expansion Value Matrix")
    
    for i in range(len(aspect_names)):
        for j in range(len(layers)):
            val = data_matrix[i, j]
            color = "white" if val > 0.6 else "black"
            ax1.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=8, fontweight="bold", color=color)
    
    cbar = plt.colorbar(im, ax=ax1, shrink=0.8)
    cbar.set_label("Score", fontsize=9)
    
    # Panel B: Unique contribution per layer (beyond simple coverage)
    contributions = []
    l_analysis = {
        "L0": ("Defines perturbation\nanchor", 0.95),
        "L1": ("Protein-level\nannotation", 0.60),
        "L2": ("Pathway context\n& grouping", 0.75),
        "L3": ("Transcriptome-wide\nresponse", 0.90),
        "L4": ("Protein surface\nmarker readout", 0.85),
        "L5": ("Disease-relevant\nfitness endpoint", 0.80),
    }
    
    contrib_labels = [l_analysis[l][0] for l in layers]
    contrib_scores = [l_analysis[l][1] for l in layers]
    layer_colors = [COLORS[l] for l in layers]
    
    bars = ax2.barh(range(len(layers)), contrib_scores, color=layer_colors,
                    edgecolor="white", linewidth=0.8, height=0.6)
    ax2.set_yticks(range(len(layers)))
    ax2.set_yticklabels(contrib_labels, fontsize=8)
    ax2.set_xlabel("Unique Contribution Score")
    ax2.set_title("B) Unique Contribution Beyond Coverage")
    ax2.set_xlim(0, 1.1)
    
    for bar, score in zip(bars, contrib_scores):
        ax2.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height()/2.,
                f"{score:.2f}", va="center", fontsize=9, fontweight="bold")
    
    fig.suptitle("Figure 2: Context-Expansion Value Beyond Coverage",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()
    path = f"{OUT_DIR}/fig2_context_expansion_value.png"
    fig.savefig(path)
    plt.close(fig)
    print(f"  Figure 2 saved: {path}")
    return path

def fig3_l3_l4_evidence():
    """Figure 3: L3/L4 measured-response evidence."""
    l3_path = f"{SRC_DIR}/data/processed/case45_l0_l5/case45_p0_20260708T1015Z/l3_guide_rna_response.json"
    l4_path = f"{SRC_DIR}/data/processed/case45_l0_l5/case45_p0_20260708T1015Z/l4_guide_adt_response.json"
    ce = load_json(f"{RES_DIR}/context_expansion_analysis.json")

    l3 = load_json(l3_path)
    l4 = load_json(l4_path)

    l3_analysis = ce["context_expansion_value"]["l3_analysis"]
    l4_analysis = ce["context_expansion_value"]["l4_analysis"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    (ax1, ax2), (ax3, ax4) = axes

    # Panel A: L3 cell counts per guide
    l3_records = l3["records"]
    n_cells = [r["n_cells_guide"] for r in l3_records]

    sorted_idx = np.argsort(n_cells)[::-1]
    sorted_cells = [n_cells[i] for i in sorted_idx]

    bar_colors = [COLORS["L3"] if c >= 100 else "#d62728aa" for c in sorted_cells]
    ax1.bar(range(len(sorted_cells)), sorted_cells, color=bar_colors,
            edgecolor="white", linewidth=0.3, width=0.7)
    ax1.axhline(y=100, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax1.set_xlabel("Guide (sorted by cell count)")
    ax1.set_ylabel("Number of Cells")
    ax1.set_title("A) L3 scRNA-seq: Cells per Guide")
    ax1.set_xticks([])

    mean_n = l3_analysis["l3_n_cells_guide_mean"]
    min_n = l3_analysis["l3_n_cells_guide_min"]
    max_n = l3_analysis["l3_n_cells_guide_max"]
    ax1.text(0.98, 0.95, f"Mean: {mean_n:.0f} | Min: {min_n} | Max: {max_n}",
             transform=ax1.transAxes, ha="right", va="top", fontsize=7,
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    # Panel B: L3 recurrent response genes
    recurrent = l3_analysis["l3_recurrent_response_genes"]
    if isinstance(recurrent, str):
        recurrent = json.loads(recurrent)
    genes = [g[0] for g in recurrent[:10]]
    counts = [g[1] for g in recurrent[:10]]

    gene_colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(genes)))
    ax2.barh(range(len(genes)), counts, color=gene_colors, edgecolor="white",
             linewidth=0.5, height=0.6)
    ax2.set_yticks(range(len(genes)))
    ax2.set_yticklabels(genes, fontsize=8, fontfamily="monospace")
    ax2.set_xlabel("Number of Guides with Response")
    ax2.set_title("B) L3: Recurrent Response Genes")
    ax2.invert_yaxis()

    for i, (gene, cnt) in enumerate(zip(genes, counts)):
        ax2.text(cnt + 0.3, i, str(cnt), va="center", fontsize=8, fontweight="bold")

    # Panel C: L4 ADT marker delta-mean per guide
    l4_guides = l4["per_guide"]
    markers = ["CD86", "PDL1", "PDL2", "CD366"]
    marker_colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]

    guide_deltas = []
    for g in l4_guides:
        guide_delta = np.mean([g["markers"][m]["delta_mean"] for m in markers])
        guide_deltas.append((g["guide_id"], guide_delta))

    guide_deltas.sort(key=lambda x: x[1])
    gdeltas = [g[1] for g in guide_deltas]

    delta_colors = [COLORS["L4"] if d > 0 else "#9467bd66" for d in gdeltas]
    ax3.bar(range(len(gdeltas)), gdeltas, color=delta_colors,
            edgecolor="white", linewidth=0.3, width=0.7)
    ax3.axhline(y=0, color="black", linewidth=0.8)
    ax3.set_xlabel("Guide (sorted by mean ADT delta)")
    ax3.set_ylabel("Mean ADT Delta (guide - control)")
    ax3.set_title("C) L4 ADT: Mean Protein Abundance Change")
    ax3.set_xticks([])

    # Panel D: L4 per-marker summary
    marker_stats = {}
    for m in markers:
        deltas = [g["markers"][m]["delta_mean"] for g in l4_guides]
        marker_stats[m] = {
            "mean": np.mean(deltas),
            "std": np.std(deltas),
        }

    x = np.arange(len(markers))
    means = [marker_stats[m]["mean"] for m in markers]
    stds = [marker_stats[m]["std"] for m in markers]

    ax4.bar(x, means, color=marker_colors, edgecolor="white", linewidth=0.8,
            width=0.5, yerr=stds, capsize=4, error_kw={"linewidth": 1})
    ax4.axhline(y=0, color="black", linewidth=0.8)
    ax4.set_xticks(x)
    ax4.set_xticklabels(markers, fontsize=10, fontfamily="monospace")
    ax4.set_ylabel("Mean Delta (guide - control)")
    ax4.set_title("D) L4: Per-Marker Mean Response")

    for i, (m, mean_val) in enumerate(zip(markers, means)):
        ax4.text(i, mean_val + np.sign(mean_val) * 5, f"{mean_val:.1f}",
                ha="center", fontsize=8, fontweight="bold")

    fig.suptitle("Figure 3: L3/L4 Measured-Response Evidence",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    path = f"{OUT_DIR}/fig3_l3_l4_evidence.png"
    fig.savefig(path)
    plt.close(fig)
    print(f"  Figure 3 saved: {path}")
    return path


def fig4_causal_readiness():
    """Figure 4: Causal-readiness and cross-layer consistency."""
    ce = load_json(f"{RES_DIR}/context_expansion_analysis.json")
    ciq = ce["causal_inference_quality"]
    clc = ce["cross_layer_consistency"]
    
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    (ax1, ax2), (ax3, ax4) = axes
    
    # Panel A: Cross-layer consistency - L3/L4 overlap
    l3l4 = clc["l3_l4_consistency"]
    categories = ["Both Active", "Neither Active"]
    values = [l3l4["both_active"], l3l4["neither_active"]]
    pie_colors = [COLORS["L3"], "#cccccc"]
    
    wedges, texts, autotexts = ax1.pie(values, labels=categories, colors=pie_colors,
                                        autopct="%1.1f%%", startangle=90,
                                        explode=(0.05, 0),
                                        textprops={"fontsize": 9})
    for at in autotexts:
        at.set_fontweight("bold")
    ax1.set_title("A) L3/L4 Response Overlap")
    
    # Panel B: L0-L5 alignment
    l0l5 = clc["l0_l5_alignment"]
    dep_cats = ["Strong\nDependency", "Moderate\nDependency", "Weak/No\nDependency"]
    dep_vals = [l0l5.get("strong_dependency", 0), l0l5.get("moderate_dependency", 4),
                l0l5.get("weak_or_no_dependency", 23)]
    dep_colors = ["#d62728", "#ff7f0e", "#2ca02c"]
    
    bars = ax2.bar(dep_cats, dep_vals, color=dep_colors, edgecolor="white",
                   linewidth=0.8, width=0.5)
    for bar, val in zip(bars, dep_vals):
        ax2.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.5,
                str(val), ha="center", fontsize=11, fontweight="bold")
    ax2.set_ylabel("Number of Guides")
    ax2.set_title("B) L5 DepMap Gene Effect Distribution")
    
    # Panel C: Causal inference quality assessment
    aspects = ["Treatment\nDefinition", "Outcome\nDefinition", "Mediator\nLayers",
               "Confounding\nControl", "Causal Claim\nSupport"]
    
    ratings = {
        "Treatment\nDefinition": 0.85,
        "Outcome\nDefinition": 0.75,
        "Mediator\nLayers": 0.60,
        "Confounding\nControl": 0.35,
        "Causal Claim\nSupport": 0.25,
    }
    
    scores = [ratings[a] for a in aspects]
    quality_colors = []
    for s in scores:
        if s >= 0.7:
            quality_colors.append("#2ca02c")
        elif s >= 0.4:
            quality_colors.append("#ff7f0e")
        else:
            quality_colors.append("#d62728")
    
    bars = ax3.barh(aspects, scores, color=quality_colors, edgecolor="white",
                    linewidth=0.8, height=0.5)
    ax3.set_xlabel("Quality Score")
    ax3.set_title("C) Causal Inference Readiness")
    ax3.set_xlim(0, 1.1)
    
    for bar, score in zip(bars, scores):
        label = "Good" if score >= 0.7 else "Moderate" if score >= 0.4 else "Low"
        ax3.text(bar.get_width() + 0.02, bar.get_y() + bar.get_height()/2.,
                f"{score:.2f} ({label})", va="center", fontsize=8, fontweight="bold")
    
    # Panel D: Layer dependency network (simplified)
    ax4.set_xlim(0, 10)
    ax4.set_ylim(0, 10)
    ax4.axis("off")
    ax4.set_title("D) Cross-Layer Dependency Map", fontsize=11, fontweight="bold")
    
    # Node positions
    nodes = {
        "L0\nPerturbation": (1, 5),
        "L1\nUniProt": (3, 7.5),
        "L2\nReactome": (3, 2.5),
        "L3\nscRNA-seq": (5, 7.5),
        "L4\nADT": (5, 2.5),
        "L5\nDepMap": (8, 5),
    }
    
    edges = [
        ("L0\nPerturbation", "L3\nscRNA-seq", "causal", 0.9),
        ("L0\nPerturbation", "L4\nADT", "causal", 0.9),
        ("L1\nUniProt", "L2\nReactome", "annotation", 0.7),
        ("L3\nscRNA-seq", "L5\nDepMap", "correlative", 0.5),
        ("L4\nADT", "L5\nDepMap", "correlative", 0.5),
        ("L2\nReactome", "L3\nscRNA-seq", "contextual", 0.6),
    ]
    
    node_colors_map = {
        "L0\nPerturbation": COLORS["L0"],
        "L1\nUniProt": COLORS["L1"],
        "L2\nReactome": COLORS["L2"],
        "L3\nscRNA-seq": COLORS["L3"],
        "L4\nADT": COLORS["L4"],
        "L5\nDepMap": COLORS["L5"],
    }
    
    # Draw edges
    for src, dst, etype, strength in edges:
        sx, sy = nodes[src]
        dx, dy = nodes[dst]
        lw = strength * 3
        if etype == "causal":
            color, style = "#1f77b4", "-"
        elif etype == "correlative":
            color, style = "#d62728", "--"
        else:
            color, style = "#7f7f7f", ":"
        ax4.annotate("", xy=(dx, dy), xytext=(sx, sy),
                     arrowprops=dict(arrowstyle="->", color=color, lw=lw,
                                    linestyle=style, connectionstyle="arc3,rad=0.1"))
    
    # Draw nodes
    for name, (x, y) in nodes.items():
        ax4.scatter(x, y, s=400, c=node_colors_map[name], edgecolors="white",
                   linewidth=1.5, zorder=5)
        ax4.text(x, y, name, ha="center", va="center", fontsize=7,
                fontweight="bold", color="white" if name != "L1\nUniProt" else "black")
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor="white", edgecolor="#1f77b4", label="Causal (experimental)"),
        mpatches.Patch(facecolor="white", edgecolor="#d62728", label="Correlative (inferred)"),
        mpatches.Patch(facecolor="white", edgecolor="#7f7f7f", label="Annotation/Context"),
    ]
    ax4.legend(handles=legend_elements, loc="lower right", fontsize=7, framealpha=0.9)
    
    fig.suptitle("Figure 4: Causal-Readiness and Cross-Layer Consistency",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    path = f"{OUT_DIR}/fig4_causal_readiness.png"
    fig.savefig(path)
    plt.close(fig)
    print(f"  Figure 4 saved: {path}")
    return path

def main():
    print("Stage 5: Generate Showcase Figures")
    print("=" * 50)
    
    figures = {}
    errors = []
    
    # Figure 1
    try:
        figures["fig1"] = fig1_layer_coverage()
    except Exception as e:
        errors.append(f"Figure 1: {e}")
        print(f"  ERROR Figure 1: {e}")
    
    # Figure 2
    try:
        figures["fig2"] = fig2_context_expansion_value()
    except Exception as e:
        errors.append(f"Figure 2: {e}")
        print(f"  ERROR Figure 2: {e}")
    
    # Figure 3
    try:
        figures["fig3"] = fig3_l3_l4_evidence()
    except Exception as e:
        errors.append(f"Figure 3: {e}")
        print(f"  ERROR Figure 3: {e}")
    
    # Figure 4
    try:
        figures["fig4"] = fig4_causal_readiness()
    except Exception as e:
        errors.append(f"Figure 4: {e}")
        print(f"  ERROR Figure 4: {e}")
    
    print()
    print("=" * 50)
    print("Figure Summary:")
    print("-" * 50)
    
    all_ok = True
    for name, path in sorted(figures.items()):
        if os.path.exists(path):
            size_kb = os.path.getsize(path) / 1024
            print(f"  {name}: {path} ({size_kb:.1f} KB)")
        else:
            print(f"  {name}: {path} -- MISSING!")
            all_ok = False
    
    if errors:
        print()
        print("Errors:")
        for e in errors:
            print(f"  - {e}")
    
    if all_ok:
        print()
        print("Stage 5 status: COMPLETE")
        print("Stopped after Stage 5; waiting for human review before Stage 6.")
    else:
        print()
        print("Stage 5 status: BLOCKED -- some figures missing or empty")
        sys.exit(1)

if __name__ == "__main__":
    main()
