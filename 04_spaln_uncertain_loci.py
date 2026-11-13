#!/usr/bin/env python3
"""Run Spaln for syntenic loci lacking supported protein homology.

Coordinates in locus_states.tsv are treated as 0-based, half-open BED intervals.
The extracted sequence is in reference (forward) orientation; Spaln tests both.
"""

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path

DEFAULT_HAL2FASTA = Path(
    "/home/fv_kalashik/bioinformatics/tools/cactus-bin-v3.3.0/bin/hal2fasta"
)
TARGET_STATUSES = {
    "annotated_without_supported_homology",
    "annotated_without_blast_hit",
}


def resolve_program(name, fallback=None):
    found = shutil.which(name)
    if found:
        return found
    if fallback and fallback.is_file():
        return str(fallback)
    raise FileNotFoundError(f"Cannot find '{name}'. Install it / add to PATH.")


def read_fasta_bases(path):
    sequence = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                sequence.append(line.strip())
    return "".join(sequence)


def run_and_capture(command, output_path, error_path, cwd=None):
    with output_path.open("w", encoding="utf-8") as out, error_path.open(
        "w", encoding="utf-8"
    ) as err:
        completed = subprocess.run(
            command, stdout=out, stderr=err, cwd=cwd, check=False
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {completed.returncode}). See {error_path}: {' '.join(map(str, command))}"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene_id", help="e.g. FBgn0035694")
    parser.add_argument(
        "--flank",
        type=int,
        default=10_000,
        help="Bases to extend on EACH side (default: 10000; --flank 10 for exactly 10 bp)",
    )
    parser.add_argument("--hal", type=Path, default=Path("Acalyptratae20.hal"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    parser.add_argument(
        "--hal2fasta", default=None, help="Path to hal2fasta binary, if not on PATH"
    )
    parser.add_argument(
        "--spaln", default=None, help="Path to spaln binary, if not on PATH"
    )
    args = parser.parse_args()

    if args.flank < 0:
        parser.error("--flank must be nonnegative")
    gene_dir = args.results_dir / args.gene_id
    table_path = gene_dir / "locus_states.tsv"

    # 01_prepare_query.py normally writes query_proteins.faa. Keep a
    # fallback for older runs that used dmel_protein.faa.
    protein = gene_dir / "query_proteins.faa"
    if not protein.is_file():
        protein = gene_dir / "dmel_protein.faa"
    if not protein.is_file():
        parser.error(
            f"Protein FASTA not found for {args.gene_id}: "
            f"expected {gene_dir / 'query_proteins.faa'} "
            f"or {gene_dir / 'dmel_protein.faa'}"
        )

    for p in (table_path, protein, args.hal):
        if not p.is_file():
            parser.error(f"Missing file: {p}")
    if not read_fasta_bases(protein):
        parser.error(f"Protein FASTA has no sequence: {protein}")
    hal2fasta = args.hal2fasta or resolve_program("hal2fasta", DEFAULT_HAL2FASTA)
    spaln = args.spaln or resolve_program("spaln")

    with table_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    selected = [
        r
        for r in rows
        if r.get("gene_id") == args.gene_id and r.get("blast_status") in TARGET_STATUSES
    ]

    if not selected:
        print(f"No rows with target BLAST statuses for {args.gene_id}")
        print("Target statuses: " + ", ".join(sorted(TARGET_STATUSES)))
        other = sorted({r.get("blast_status", "") for r in rows if r.get("gene_id") == args.gene_id})
        print("Statuses found: " + (", ".join(other) if other else "none"))
        print("NOTE: no_coding_annotation_needs_genewise is not automatically analyzed by this script.")
        return

    out_root = args.results_dir / args.gene_id / "spaln"
    out_root.mkdir(parents=True, exist_ok=True)
    summary = []
    for row in selected:
        species, locus, chrom = row["species"], row["locus_id"], row["chromosome"]
        start, end = int(row["start"]), int(row["end"])
        if start < 0 or end <= start:
            print(
                f"WARN invalid BED interval: {species} {locus}: {start}:{end}",
                file=sys.stderr,
            )
            summary.append((species, locus, "invalid_coordinates", "", "", "", ""))
            continue

        # BED / HAL genomic coordinates are zero-based, end-exclusive.
        extracted_start = max(0, start - args.flank)
        requested_end = end + args.flank
        requested_length = requested_end - extracted_start
        dest = out_root / f"{species}_{locus}"
        dest.mkdir(parents=True, exist_ok=True)
        dna = dest / "genomic_window.fasta"
        hal_err = dest / "hal2fasta.stderr.txt"
        spaln_out = dest / "spaln_alignment.txt"
        spaln_err = dest / "spaln.stderr.txt"
        coords = dest / "window_coordinates.tsv"
        print(
            f"[{species} / {locus}] {chrom}:{extracted_start}-{requested_end} ({requested_length} bp requested)",
            flush=True,
        )

        try:
            extract_command = [
                str(hal2fasta),
                str(args.hal.resolve()),
                species,
                "--sequence",
                chrom,
                "--start",
                str(extracted_start),
                "--length",
                str(requested_length),
                "--upper",
            ]
            run_and_capture(extract_command, dna, hal_err)
            sequence = read_fasta_bases(dna)
            if len(sequence) != requested_length:
                raise ValueError(
                    f"Extracted {len(sequence)} bp, expected {requested_length} bp"
                )
            if sequence.upper().count("N") == len(sequence):
                raise ValueError("Extracted sequence consists entirely of Ns")
            with coords.open("w", encoding="utf-8") as fh:
                fh.write(
                    "species\tlocus_id\tchromosome\twindow_start_0based\twindow_end_exclusive\twindow_length_bp\tflank_bp\n"
                )
                fh.write(
                    f"{species}\t{locus}\t{chrom}\t{extracted_start}\t{requested_end}\t{len(sequence)}\t{args.flank}\n"
                )

            # Spaln mode Q2: local genomic segment against query protein, HSP-assisted.
            # O1 prints the alignment; S3 examines both DNA strands.
            command = [
                str(spaln),
                "-Q2",
                "-O1",
                "-S3",
                str(dna.resolve()),
                str(protein.resolve()),
            ]
            run_and_capture(command, spaln_out, spaln_err)
            summary.append(
                (
                    species,
                    locus,
                    "completed_check_alignment",
                    chrom,
                    str(extracted_start),
                    str(requested_end),
                    str(spaln_out),
                )
            )
            print(f"  -> Saved: {spaln_out}")
        except (OSError, RuntimeError, ValueError) as exc:
            print(f"  ERROR: {exc}", file=sys.stderr)
            summary.append(
                (
                    species,
                    locus,
                    "error",
                    chrom,
                    str(extracted_start),
                    str(requested_end),
                    str(exc),
                )
            )

    summary_path = out_root / "spaln_run_summary.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "species",
                "locus_id",
                "run_status",
                "chromosome",
                "window_start_0based",
                "window_end_exclusive",
                "alignment_or_error",
            ]
        )
        writer.writerows(summary)
    print(f"\nSummary: {summary_path}")
    print(
        "NOTE: completed_check_alignment means the command ran, NOT that homology was proven."
    )
    if any(r[2] == "error" for r in summary):
        sys.exit(1)


if __name__ == "__main__":
    main()
