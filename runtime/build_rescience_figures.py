"""Render five retained analyses directly from verified replay outputs.

No scientific computation or source-data changes. SVGs are editable deliverables;
PNG previews belong in the caller's task-marked scratch folder.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
FULL = ROOT / "rescience_offline_candidate_20260927/full_replay/runs"
INPUTS = {
    "core": FULL / "core_seed1/result.json",
    "loo": FULL / "loo_frozen/result.json",
    "gsea": FULL / "gsea_current/result.json",
    "itpr3": ROOT / "results/itpr3_claim_verification.csv",
}
EXPECTED = {
    "core": "a52f57efdc4c82a08264ce4b17ce9c3eed45ff079bec1d2ed21df5d93c1cb221",
    "loo": "29e349d0ef9198d15a1b085e07cb78790ce3e471e321e1b74cb708677c1a3aa6",
    "gsea": "3bba0847db786c4a591d73263c2d95355768146ebdd83600374a7241c2d2d790",
}
MODELS = ["CEP290_LCA", "P23H", "rd10"]
LABELS = ["Human CEP290-LCA", "P23H rat", "rd10 mouse"]
COLOURS = ["#0072B2", "#D55E00", "#009E73"]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def make(args):
    assert args.preview_dir.is_dir()
    assert (args.preview_dir / ".codex-task-temp").is_file()
    assert not args.out.exists(), "Refusing to overwrite a prior figure revision"
    hashes = {k: sha(p) for k, p in INPUTS.items()}
    for key, expected in EXPECTED.items():
        assert hashes[key] == expected, f"Changed input: {key}"
    core, loo, gsea = [json.loads(INPUTS[k].read_text(encoding="utf-8"))
                       for k in ["core", "loo", "gsea"]]
    with INPUTS["itpr3"].open(encoding="utf-8-sig", newline="") as stream:
        itpr3 = list(csv.DictReader(stream))
    assert [r["species"] for r in itpr3] == ["human", "rat", "mouse"]
    args.out.mkdir(parents=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none", "svg.hashsalt": "e071-rescience-v1"})
    figures = []

    def save(fig, stem, source_values):
        svg = args.out / (stem + ".svg")
        fig.savefig(svg, metadata={"Date": None}, bbox_inches="tight")
        preview = args.preview_dir / (stem + ".png")
        fig.savefig(preview, dpi=160, bbox_inches="tight")
        plt.close(fig)
        figures.append({"svg": svg.name, "sha256": sha(svg),
                        "preview": str(preview), "plotted_values": source_values})

    fig, ax = plt.subplots(figsize=(7.8, 4.8), layout="constrained")
    vals = [float(r["logFC"]) for r in itpr3]
    ax.bar(LABELS, vals, color=COLOURS, width=.55)
    ax.axhline(0, color="#444444", lw=.8)
    ax.set(ylabel="Deposited log2 fold change (HAL − DMSO)",
           title="Deposited ITPR3 / Itpr3 effects", ylim=(-1.5, 1.35))
    for i, (value, row) in enumerate(zip(vals, itpr3)):
        text = f"p = {float(row['nominal_p']):.3g}\nadj. p = {float(row['adjusted_p']):.3g}"
        ax.text(i, value + (.1 if value >= 0 else -.1), text,
                ha="center", va="bottom" if value >= 0 else "top")
    save(fig, "figure1_deposited_itpr3", itpr3)

    counts = {m: {k: core["models"][m][k] for k in ["nominal_p05_fc12", "fdr_selected"]}
              for m in MODELS}
    fig, ax = plt.subplots(figsize=(7.8, 4.7), layout="constrained")
    x = np.arange(3)
    for shift, key, label, colour in [(-.18, "nominal_p05_fc12", "Nominal p < 0.05", "#0072B2"),
                                    (.18, "fdr_selected", "Adjusted p < 0.05", "#E69F00")]:
        bar = ax.bar(x + shift, [counts[m][key] for m in MODELS], .36,
                     label=label, color=colour)
        ax.bar_label(bar, padding=3)
    ax.set_xticks(x, LABELS)
    ax.set(ylabel="Selected genes", ylim=(0, 4300),
           title="Genes meeting nominal and adjusted thresholds")
    ax.legend(frameon=False)
    fig.text(.5, -.025, "Both gates also require |log2 fold change| ≥ log2(1.2).",
             ha="center", fontsize=9)
    save(fig, "figure2_gene_thresholds", counts)

    pairs = ["CEP290_LCA__P23H", "CEP290_LCA__rd10", "P23H__rd10"]
    pair_labels = ["Human–P23H", "Human–rd10", "P23H–rd10"]
    records = [loo["cross_model_loo"][p] for p in pairs]
    full = np.array([r["full"]["spearman"] for r in records])
    lower = np.array([r["min"] for r in records])
    upper = np.array([r["max"] for r in records])
    assert np.all(lower <= full) and np.all(full <= upper)
    fig, ax = plt.subplots(figsize=(7.8, 4.6), layout="constrained")
    ax.errorbar(x, full, yerr=np.stack([full-lower, upper-full]), fmt="o",
                color="#0072B2", capsize=6, lw=1.6)
    ax.axhline(0, color="#777777", lw=.8, linestyle="--")
    ax.set_xticks(x, pair_labels)
    ax.set(ylabel="Spearman correlation", ylim=(-.32, .49),
           title="Cross-model correlations under sample deletion")
    fig.text(.5, -.025, "Points: full-sample contrasts. Bars: tested state range, not a confidence interval.",
             ha="center", fontsize=9)
    save(fig, "figure3_deletion_ranges", {
        p: {"full": r["full"]["spearman"], "min": r["min"], "max": r["max"],
            "state_pairs": r["all_state_pairs"]} for p, r in zip(pairs, records)})

    overlap = core["author_overlap"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), layout="constrained")
    ax.errorbar([0], [overlap["null_ge2_mean"]], yerr=[overlap["null_ge2_sd"]],
                fmt="o", capsize=7, color="#0072B2", label="Random-set mean ± SD")
    ax.scatter([.18], [overlap["observed_ge2"]], marker="x", s=75, color="#D55E00",
               label="Observed overlap")
    ax.set(xlim=(-.5,.65), ylim=(0,80), ylabel="Genes present in at least two lists",
           title="Shared-gene count and random-set null")
    ax.set_xticks([.09], ["Source-like reconstruction"])
    ax.legend(frameon=False, loc="lower center")
    ax.text(.03,.97, f"10,000 draws; Monte Carlo p = {overlap['p_ge2']:.3f}",
            transform=ax.transAxes, va="top", fontsize=9)
    save(fig, "figure4_overlap_null", overlap)

    terms = ["Oxidative phosphorylation", "Ribosome", "Circadian entrainment"]
    rows = [r for rows in gsea["highlighted"].values() for r in rows]
    selected = []
    for term in terms:
        match = [r for r in rows if r["term"] == term]
        assert len(match) == 1
        selected.append(match[0])
    fig, ax = plt.subplots(figsize=(8.2, 5.0), layout="constrained")
    for i, (model, label, colour) in enumerate(zip(MODELS, LABELS, COLOURS)):
        values = [r["NES_"+model] for r in selected]
        positions = x + (i-1)*.25
        ax.bar(positions, values, width=.24, label=label, color=colour)
        for position, value, row in zip(positions, values, selected):
            if row["FDR_"+model] < .05:
                ax.text(position, value+.08 if value >= 0 else value-.08, "*",
                        ha="center", va="bottom" if value >= 0 else "top", fontsize=14)
    ax.axhline(0, color="#777777", lw=.8)
    ax.set_xticks(x, ["Oxidative\nphosphorylation", "Ribosome", "Circadian\nentrainment"])
    ax.set(ylabel="Normalized enrichment score", ylim=(-1.25, 4.25),
           title="Ranked enrichment in three source-highlighted pathways")
    ax.legend(frameon=False, ncol=3, loc="upper center", fontsize=9)
    fig.text(.5,-.025, "* Returned GSEA FDR < 0.05; absence of a star is not evidence of no biological effect.",
             ha="center", fontsize=9)
    save(fig, "figure5_ranked_pathways", selected)
    assert all(sha(INPUTS[k]) == v for k,v in hashes.items())
    (args.out / "FIGURE_PROVENANCE.json").write_text(json.dumps({
        "scope": "Five original analysis displays retained; no scientific rerun or input change.",
        "input_hashes": hashes, "input_paths": {k:str(p) for k,p in INPUTS.items()},
        "figures": figures, "visual_review": "PENDING", "independent_review": "PENDING",
    }, indent=2, allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"figures_created": len(figures), "out": str(args.out),
                      "inputs_unchanged": True}))

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--preview-dir", type=Path, required=True)
    make(p.parse_args())
