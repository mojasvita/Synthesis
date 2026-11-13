#!/usr/bin/env python3
"""Make a conservative species-tree figure for one query gene.

Inputs (relative to project directory):
  reference_data/phylogeny/species_tree.nwk
  results/<FBgn>/locus_states.tsv
  results/<FBgn>/spaln/spaln_run_summary.tsv  (optional)

Outputs:
  results/<FBgn>/locus_states_with_spaln.tsv
  results/<FBgn>/species_status.tsv
  results/<FBgn>/species_tree.png / .pdf / .svg

Spaln results are supplementary evidence: no hit is NOT evidence of absence,
and any alignment from Spaln is NOT automatically verified homology.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def load_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=fields, extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def normalize_species(value):
    return (value or "").strip().replace(" ", "_").replace("-", "_").lower()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("gene_id", help="Gene to plot, e.g. FBgn0035694")
    ap.add_argument(
        "--tree", type=Path, default=Path("reference_data/phylogeny/species_tree.nwk")
    )
    ap.add_argument("--results-dir", type=Path, default=Path("results"))
    args = ap.parse_args()

    try:
        from Bio import Phylo
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError as exc:
        ap.error(
            f"Missing plotting dependency: {exc}. Install: conda install -c conda-forge biopython matplotlib"
        )

    root = args.results_dir / args.gene_id
    table_path = root / "locus_states.tsv"
    if not table_path.is_file() or not args.tree.is_file():
        ap.error(f"Input missing: {table_path} or {args.tree}")
    records = load_tsv(table_path)
    records = [r for r in records if r.get("gene_id") == args.gene_id]
    if not records:
        ap.error(f"No rows for gene {args.gene_id} in {table_path}")

    # Link Spaln results to individual locus IDs, not just species.
    spaln_path = root / "spaln" / "spaln_run_summary.tsv"
    spaln_rows = load_tsv(spaln_path) if spaln_path.is_file() else []
    spaln_map = {(r["species"], r["locus_id"]): r for r in spaln_rows}
    augmented = []
    for original in records:
        row = dict(original)
        r = spaln_map.get((row["species"], row["locus_id"]))
        status = "not_tested"
        if r:
            if r["run_status"] != "completed_check_alignment":
                status = "run_error"
            else:
                output = Path(r["alignment_or_error"])
                # Summary often stores a project-relative path; also support absolute paths.
                if not output.is_file():
                    status = "alignment_file_missing"
                else:
                    content = output.read_text(encoding="utf-8", errors="replace")
                    status = (
                        "candidate_alignment_unvalidated"
                        if "ALIGNMENT" in content and "Score" in content
                        else "no_alignment_reported"
                    )
        row["spaln_status"] = status
        augmented.append(row)
    out_augmented = root / "locus_states_with_spaln.tsv"
    write_tsv(out_augmented, augmented, list(records[0]) + ["spaln_status"])

    # Collapse loci to species. One strong syntenic protein match is enough to
    # mark evidence at species level; NOT enough to prove one-to-one orthology.
    by_species = defaultdict(list)
    for r in augmented:
        by_species[normalize_species(r["species"])].append(r)
    species_data = {}
    for norm, loci in by_species.items():
        strong = sum(r["blast_status"] == "annotated_homolog_supported" for r in loci)
        status = "coding_homolog_supported" if strong else "uncertain"
        species_data[norm] = {
            "species": loci[0]["species"],
            "state": status,
            "n_loci": len(loci),
            "n_supported_loci": strong,
            "blast_statuses": ",".join(sorted(set(r["blast_status"] for r in loci))),
            "spaln_statuses": ",".join(sorted(set(r["spaln_status"] for r in loci))),
        }

    tree = Phylo.read(str(args.tree), "newick")

    # Удаляем технические листья из дерева только для визуализации
    for leaf_name in ("14", "11"):
        matches = [tip for tip in tree.get_terminals() if tip.name == leaf_name]
        for tip in matches:
            tree.prune(tip)

    # Удаляем оставшиеся терминальные узлы без названия
    for tip in list(tree.get_terminals()):
        if not tip.name:
            tree.prune(tip)

    terminals = tree.get_terminals()

    print(f"Species tree: {len(terminals)} terminal taxa")

    leaf_names = [normalize_species(t.name) for t in terminals]
    if len(set(leaf_names)) != len(leaf_names):
        ap.error("Ambiguous/duplicate species names after normalization in Newick tips")
    matched = set(leaf_names) & set(species_data)
    if not matched:
        ap.error(
            "No tree tips matched TSV species. Inspect Newick names and TSV species column."
        )

    # Include all Newick terminal species, including those not assessed.
    status_rows = []
    for tip in terminals:
        key = normalize_species(tip.name)
        if key == "drosophila_melanogaster":
            state, label = "query_reference", "Query species (D. melanogaster)"
        elif key in species_data:
            state, label = species_data[key]["state"], species_data[key][
                "state"
            ].replace("_", " ")
        else:
            state, label = "not_assessed", "Not assessed"
        data = species_data.get(key, {})
        status_rows.append(
            {
                "species": tip.name,
                "tree_state": state,
                "n_loci": data.get("n_loci", 0),
                "n_supported_loci": data.get("n_supported_loci", 0),
                "blast_statuses": data.get("blast_statuses", ""),
                "spaln_statuses": data.get("spaln_statuses", ""),
            }
        )
    out_status = root / "species_status.tsv"
    write_tsv(out_status, status_rows, list(status_rows[0]))

    palette = {
        "coding_homolog_supported": "#28866a",
        "uncertain": "#e5a037",
        "query_reference": "#3d73ba",
        "not_assessed": "#a0a0a0",
    }
    state_by_tip = {
        normalize_species(row["species"]): row["tree_state"] for row in status_rows
    }

    # A cladogram presents only topology (edge lengths are not interpreted as time).
    # Tree y-order is preserved; source branch lengths are intentionally ignored.
    ymap = {tip: i for i, tip in enumerate(reversed(terminals))}
    xmap = {}

    def place(clade, depth):
        xmap[clade] = depth
        if clade not in ymap:
            for child in clade.clades:
                place(child, depth + 1)
            ymap[clade] = (
                min(ymap[ch] for ch in clade.clades)
                + max(ymap[ch] for ch in clade.clades)
            ) / 2

    place(tree.root, 0)
    maxdepth = max(xmap[tip] for tip in terminals)
    fig, ax = plt.subplots(figsize=(14, max(6.0, len(terminals) * 0.40 + 1.5)))
    for clade in tree.find_clades(order="preorder"):
        if not clade.clades:
            continue
        low = min(ymap[ch] for ch in clade.clades)
        high = max(ymap[ch] for ch in clade.clades)
        ax.plot(
            [xmap[clade], xmap[clade]], [low, high], color="#686868", lw=1.2, zorder=1
        )
        for child in clade.clades:
            ax.plot(
                [xmap[clade], xmap[child]],
                [ymap[child], ymap[child]],
                color="#686868",
                lw=1.2,
                zorder=1,
            )
    for tip in terminals:
        state = state_by_tip[normalize_species(tip.name)]
        ax.scatter(xmap[tip] + 0.06, ymap[tip], s=54, color=palette[state], zorder=5)
        nice_name = tip.name.replace("_", " ")
        data = species_data.get(normalize_species(tip.name))
        suffix = ""
        if data and data["n_loci"] > 1:
            suffix = f"  [{data['n_loci']} loci]"
        ax.text(
            maxdepth + 0.45,
            ymap[tip],
            nice_name + suffix,
            fontsize=9.3,
            va="center",
            color="#222222",
        )
    ax.set_xlim(-0.25, maxdepth + 8.5)
    ax.set_ylim(-0.9, len(terminals) - 0.1)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(
        f"{args.gene_id}: syntenic coding-homology evidence on the species tree",
        fontsize=14,
        loc="left",
        pad=18,
    )
    legend_items = [
        ("coding_homolog_supported", "Annotated coding homolog supported"),
        ("uncertain", "Uncertain (no supported protein homology)"),
        ("query_reference", "D. melanogaster query"),
        ("not_assessed", "Not assessed"),
    ]
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=palette[k],
            markeredgecolor="none",
            markersize=9,
            label=label,
        )
        for k, label in legend_items
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=2,
        frameon=False,
        fontsize=9,
        bbox_to_anchor=(0.5, 0.02),
    )
    fig.text(
        0.09,
        0.095,
        "Cladogram topology only (not time-scaled). Spaln candidates are unvalidated; "
        "no hit ≠ gene absence.",
        fontsize=8.4,
        color="#555555",
    )
    fig.subplots_adjust(left=0.05, right=0.97, top=0.91, bottom=0.15)
    for extension in ("png", "pdf", "svg"):
        out = root / f"species_tree.{extension}"
        fig.savefig(out, dpi=250, bbox_inches="tight")
        print("Saved:", out)
    plt.close(fig)
    print("Saved:", out_augmented)
    print("Saved:", out_status)
    print(f"Matched assessed species to tree: {len(matched)}/{len(species_data)}")
    missing = sorted(set(species_data) - set(leaf_names))
    if missing:
        print("WARNING: Assessed species absent from tree:", ", ".join(missing))
    print("Tree state counts:")
    from collections import Counter

    for name, count in Counter(r["tree_state"] for r in status_rows).items():
        print(f"  {name}: {count}")


if __name__ == "__main__":
    main()
