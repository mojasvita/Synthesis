#!/usr/bin/env python3
"""Step 3. Link syntenic protein-coding genes to proteins using GFF3 CDS
protein_id attributes, run BLASTP and write locus_states.tsv.

Requires syntenic_loci.tsv from step 02 and query_proteins.faa from step 01.
No whole-proteome search and no GeneWise inference are performed here.
"""

import argparse
import csv
import gzip
import subprocess
from collections import defaultdict
from pathlib import Path

EXTRA = [
    "target_proteins",
    "best_hit",
    "identity_pct",
    "query_coverage_pct",
    "evalue",
    "blast_status",
    "blast_note",
]


def open_text(path):
    return (
        gzip.open(path, "rt")
        if str(path).endswith(".gz")
        else open(path, encoding="utf-8")
    )


def find_file(root, names):
    base = root / "ncbi_dataset" / "data"
    if not base.is_dir():
        return None
    for name in names:
        paths = sorted(base.rglob(name))
        if paths:
            return paths[0]
    return None


def parse_attributes(text):
    result = {}
    for field in text.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            result[key] = value
    return result


def protein_ids_for_genes(gff_path, wanted_gene_ids):
    """Return {gene-ID: set(protein_id)} using gene -> mRNA -> CDS Parent links.

    Also accept a direct CDS Parent=gene-ID relationship. Two passes avoid
    reliance on record order; do not use unreliable FASTA description matching.
    """
    wanted = set(wanted_gene_ids)
    transcript_to_gene = {}
    with open_text(gff_path) as inp:
        for line in inp:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"mRNA", "transcript"}:
                continue
            attrs = parse_attributes(fields[8])
            transcript = attrs.get("ID")
            for parent in attrs.get("Parent", "").split(","):
                if parent in wanted and transcript:
                    transcript_to_gene[transcript] = parent

    by_gene = defaultdict(set)
    with open_text(gff_path) as inp:
        for line in inp:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = parse_attributes(fields[8])
            protein_id = attrs.get("protein_id")
            if not protein_id:
                continue
            for parent in attrs.get("Parent", "").split(","):
                gene_id = parent if parent in wanted else transcript_to_gene.get(parent)
                if gene_id:
                    by_gene[gene_id].add(protein_id)
    return by_gene


def proteins_by_accession(faa_path, accessions):
    """Read only requested FASTA IDs, matching accession as first header token."""
    wanted = set(accessions)
    result = {}
    current = None
    sequence = []
    with open_text(faa_path) as inp:
        for line in inp:
            if line.startswith(">"):
                if current is not None:
                    result[current] = "".join(sequence)
                accession = line[1:].split()[0]
                current = accession if accession in wanted else None
                sequence = []
            elif current is not None:
                sequence.append(line.strip())
        if current is not None:
            result[current] = "".join(sequence)
    return result


def write_fasta(path, records):
    with path.open("w", encoding="utf-8") as out:
        for accession, sequence in records:
            out.write(f">{accession}\n")
            for i in range(0, len(sequence), 80):
                out.write(sequence[i : i + 80] + "\n")


def best_blast_hit(blast_tsv):
    """Best single HSP by E-value, then bitscore. Coverage is query residues."""
    best = None
    with blast_tsv.open(encoding="utf-8") as inp:
        for line in inp:
            f = line.rstrip("\n").split("\t")
            if len(f) != 8:
                continue
            # qseqid sseqid pident qstart qend qlen evalue bitscore
            try:
                qstart, qend, qlen = int(f[3]), int(f[4]), int(f[5])
                identity, evalue, bitscore = float(f[2]), float(f[6]), float(f[7])
            except ValueError:
                continue
            if qlen <= 0:
                continue
            coverage = 100.0 * (abs(qend - qstart) + 1) / qlen
            hit = (evalue, -bitscore, f[1], identity, coverage)
            if best is None or hit < best:
                best = hit
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene_id")
    parser.add_argument("--out", default="results")
    parser.add_argument("--reference", default="reference_data")
    parser.add_argument("--blastp", default="blastp")
    parser.add_argument("--evalue", type=float, default=0.05)
    parser.add_argument(
        "--min-qcov", type=float, default=0, help="Minimum query coverage in percent"
    )
    args = parser.parse_args()

    out = Path(args.out) / args.gene_id
    source = out / "syntenic_loci.tsv"
    query = out / "query_proteins.faa"
    if not source.is_file():
        parser.error(f"Missing {source}; run 02_find_syntenic_loci.py first")
    if not query.is_file():
        parser.error(f"Missing {query}; run 01_prepare_query.py first")

    with source.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames
        if not fields or not {
            "species",
            "status",
            "coding_gene_ids",
            "locus_id",
        } <= set(fields):
            parser.error("Unexpected columns in syntenic_loci.tsv")
        columns = list(fields) + [f for f in EXTRA if f not in fields]
        rows = list(reader)

    # Process each species once (not one GFF/FASTA pass per locus).
    species = sorted(
        {r["species"] for r in rows if r["status"] == "annotated_locus_needs_blastp"}
    )
    for org in species:
        related_rows = [
            r
            for r in rows
            if r["species"] == org and r["status"] == "annotated_locus_needs_blastp"
        ]
        base = Path(args.reference) / org
        gff = find_file(
            base, ("genomic.gff", "genomic.gff3", "genomic.gff.gz", "genomic.gff3.gz")
        )
        faa = find_file(base, ("protein.faa", "protein.faa.gz"))
        if not gff or not faa:
            for row in related_rows:
                row.update(
                    blast_status="annotation_or_protein_fasta_missing",
                    blast_note="Required GFF3 or protein FASTA not found",
                )
            continue

        wanted = {
            gid.strip()
            for row in related_rows
            for gid in row["coding_gene_ids"].split(",")
            if gid.strip()
        }
        gene_to_protein = protein_ids_for_genes(gff, wanted)
        accessions = (
            set().union(*(gene_to_protein.values())) if gene_to_protein else set()
        )
        sequences = proteins_by_accession(faa, accessions)

        for row in related_rows:
            ids = [
                gid.strip() for gid in row["coding_gene_ids"].split(",") if gid.strip()
            ]
            protein_ids = sorted(
                set().union(*(gene_to_protein.get(gid, set()) for gid in ids))
            )
            if not protein_ids:
                row.update(
                    blast_status="protein_mapping_unresolved",
                    blast_note="No protein_id linked through GFF3 gene -> mRNA -> CDS",
                )
                continue
            found = [(pid, sequences[pid]) for pid in protein_ids if pid in sequences]
            if not found:
                row.update(
                    blast_status="protein_mapping_unresolved",
                    blast_note="GFF3 protein_id not present in protein.faa: "
                    + ",".join(protein_ids),
                )
                continue

            target = out / f"{row['locus_id']}_proteins.faa"
            blast_out = out / f"{row['locus_id']}_blastp.tsv"
            write_fasta(target, found)
            row["target_proteins"] = ",".join(pid for pid, _ in found)
            subprocess.run(
                [
                    args.blastp,
                    "-query",
                    str(query),
                    "-subject",
                    str(target),
                    "-outfmt",
                    "6 qseqid sseqid pident qstart qend qlen evalue bitscore",
                    "-out",
                    str(blast_out),
                ],
                check=True,
            )
            hit = best_blast_hit(blast_out)
            if hit is None:
                row["blast_status"] = "annotated_without_blast_hit"
            else:
                row.update(
                    best_hit=hit[2],
                    identity_pct=f"{hit[3]:.2f}",
                    query_coverage_pct=f"{hit[4]:.1f}",
                    evalue=f"{hit[0]:.3g}",
                )
                row["blast_status"] = (
                    "annotated_homolog_supported"
                    if hit[0] < args.evalue and hit[4] >= args.min_qcov
                    else "annotated_without_supported_homology"
                )
            row["blast_note"] = (
                "BLAST E-value from proteins linked to the syntenic locus, "
                "not a proteome-wide search; query coverage is for one HSP"
            )
            print(org, row["locus_id"], row["blast_status"], row.get("best_hit", ""))
    for row in rows:
        if row["status"] != "annotated_locus_needs_blastp":
            row["blast_status"] = "not_applicable"

    dest = out / "locus_states.tsv"
    with dest.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=columns, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
    print("Saved:", dest)


if __name__ == "__main__":
    main()
