
#!/usr/bin/env python3

import csv
import gzip
import re
from pathlib import Path
from collections import defaultdict, Counter

# ============================================================
# 1. PATHS AND SPECIES
# ============================================================

BASE = Path("/home/e.nikolaeva/orfan_gene/drosophila_proteomes/proteomes_1isoform")
RESULTS = BASE / "OrthoFinder" / "Results_Nov12_1"

PROTEOME = BASE / "D_melanogaster.faa"

OUTDIR = BASE / "Dmel_origin_analysis"
OUTDIR.mkdir(exist_ok=True)

DMEL = "D_melanogaster"
DSIM = "D_simulans"
DSEC = "D_sechellia"

TRIO = {DMEL, DSIM, DSEC}


# ============================================================
# 2. HELPERS
# ============================================================

def read_fasta_ids(path):
    """Return protein IDs from a FASTA file."""
    ids = set()

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])

    return ids


def read_tsv(path):
    """Read a tab-separated file."""
    with open(path, newline="") as f:
        return list(csv.reader(f, delimiter="\t"))


def split_protein_ids(cell):
    """Split a cell containing comma-separated protein IDs."""
    if not cell or cell in {"", "-", "NA"}:
        return []

    return [
        x.strip()
        for x in cell.split(",")
        if x.strip()
    ]


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def find_file(filename):
    """Search for a file under BASE."""
    matches = list(BASE.rglob(filename))
    return matches[0] if matches else None


def read_diamond_no_hit_ids():
    """
    Read the corrected DIAMOND summary.
    Returns:
      - set of D. melanogaster proteins with zero cross-species hits
      - set of proteins for which DIAMOND status is known
    """
    path = find_file("unassigned_diamond_summary.tsv")

    no_hits = set()
    known = set()

    if path is None:
        print("WARNING: DIAMOND summary not found.")
        return no_hits, known

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            query_species = row.get("query_species", "")
            query_protein = row.get("query_protein", "")
            query = row.get("query", "")

            if query_species == DMEL and query_protein:
                protein = query_protein
            elif query.startswith(DMEL + "|"):
                protein = query.split("|", 1)[1]
            else:
                continue

            try:
                n_species = int(row.get("target_species_count", ""))
            except (ValueError, TypeError):
                continue

            known.add(protein)

            if n_species == 0:
                no_hits.add(protein)

    return no_hits, known


# ============================================================
# 3. READ D. MELANOGASTER PROTEOME
# ============================================================

print("Reading D. melanogaster proteome...")

dmel_proteins = read_fasta_ids(PROTEOME)

print("Proteins in D. melanogaster proteome:", len(dmel_proteins))


# ============================================================
# 4. READ ORTHOGROUPS
# ============================================================

print("Reading Orthogroups.tsv...")

og_file = RESULTS / "Orthogroups" / "Orthogroups.tsv"
og_table = read_tsv(og_file)

header = og_table[0]

species_columns = {
    species: i
    for i, species in enumerate(header)
}

required_species = [DMEL, DSIM, DSEC]

for species in required_species:
    if species not in species_columns:
        raise ValueError(
            f"Species column not found in Orthogroups.tsv: {species}"
        )

# protein -> orthogroup(s)
protein_to_ogs = defaultdict(list)

# orthogroup -> species -> protein IDs
og_data = {}

for row in og_table[1:]:
    if not row:
        continue

    og = row[0]
    species_proteins = {}

    for species, col in species_columns.items():
        if col >= len(row):
            continue

        ids = split_protein_ids(row[col])

        if ids:
            species_proteins[species] = ids

        for protein in ids:
            if species == DMEL:
                protein_to_ogs[protein].append(og)

    og_data[og] = species_proteins


# ============================================================
# 5. READ UNASSIGNED GENES
# ============================================================

print("Reading Orthogroups_UnassignedGenes.tsv...")

unassigned_file = (
    RESULTS / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv"
)

unassigned_table = read_tsv(unassigned_file)

unassigned_dmel = set()

unassigned_header = unassigned_table[0]

if DMEL not in unassigned_header:
    raise ValueError(
        "D_melanogaster column not found in UnassignedGenes.tsv"
    )

dmel_col = unassigned_header.index(DMEL)

for row in unassigned_table[1:]:
    if dmel_col >= len(row):
        continue

    for protein in split_protein_ids(row[dmel_col]):
        unassigned_dmel.add(protein)



# ============================================================
# 6. READ ORTHOLOG TABLES
# ============================================================

print("Reading pairwise ortholog tables...")

ortholog_dmel = set()

ortholog_dir = RESULTS / "Orthologs"

for path in list(ortholog_dir.glob("*.tsv")) + list(ortholog_dir.glob("*.tsv.gz")):

    opener = gzip.open if path.suffix == ".gz" else open

    with opener(path, "rt") as f:
        for line in f:
            # Split line into fields and then into individual IDs.
            tokens = re.split(r"[,;\s]+", line.strip())

            for token in tokens:
                if token in dmel_proteins:
                    ortholog_dmel.add(token)


# ============================================================
# 7. READ DUPLICATION EVENTS
# ============================================================

print("Reading duplication events...")

dup_file = RESULTS / "Gene_Duplication_Events" / "Duplications.tsv"

duplication_dmel = set()

if dup_file.exists():
    with open(dup_file) as f:
        for line in f:
            tokens = re.split(r"[,;\s]+", line.strip())

            for token in tokens:
                if token in dmel_proteins:
                    duplication_dmel.add(token)

# ============================================================
# 8. TECHNICAL CLASSIFICATION OF ALL D. MELANOGASTER PROTEINS
# ============================================================

print("Classifying D. melanogaster proteins...")

all_rows = []

for protein in sorted(dmel_proteins):

    in_og = protein in protein_to_ogs
    in_unassigned = protein in unassigned_dmel
    has_ortholog = protein in ortholog_dmel
    has_duplication = protein in duplication_dmel

    ogs = protein_to_ogs.get(protein, [])

    species_sets = []
    for og in ogs:
        species_sets.append(set(og_data[og].keys()))

    if in_unassigned:
        classification = "UNASSIGNED_CANDIDATE"
    elif has_ortholog and has_duplication:
        classification = "ORTHOLOGS_AND_DUPLICATION"
    elif has_ortholog:
        classification = "ORTHOLOGS_DETECTED"
    elif has_duplication:
        classification = "DUPLICATION_DETECTED"
    elif in_og:
        classification = "HOMOLOGY_GROUP_ONLY"
    else:
        classification = "NO_ORTHOGROUP_OR_ORTHOLOG_RECORD"

    all_rows.append([
        protein,
        classification,
        ";".join(ogs),
        int(in_og),
        int(in_unassigned),
        int(has_ortholog),
        int(has_duplication),
    ])


write_tsv(
    OUTDIR / "Dmel_all_genes.tsv",
    [
        "protein_id",
        "technical_classification",
        "orthogroups",
        "in_orthogroup",
        "unassigned",
        "ortholog_record",
        "duplication_record",
    ],
    all_rows
)


# ============================================================
# 9. SUMMARY OF CLASSIFICATION
# ============================================================

counts = Counter(row[1] for row in all_rows)

summary_rows = [
    [category, count]
    for category, count in sorted(counts.items())
]

write_tsv(
    OUTDIR / "Dmel_origin_summary.tsv",
    ["technical_classification", "protein_count"],
    summary_rows
)


# ============================================================
# 10. ORTHOGROUPS RESTRICTED TO THE THREE SPECIES
# ============================================================

print("Finding orthogroups restricted to the Drosophila trio...")

trio_og_rows = []
trio_gene_rows = []
all_three_rows = []

for og, species_proteins in og_data.items():

    present_species = set(species_proteins.keys())

    # Must contain D. melanogaster and at least one close relative.
    if DMEL not in present_species:
        continue

    if not (present_species & {DSIM, DSEC}):
        continue

    # No species outside the trio may be present.
    if not present_species.issubset(TRIO):
        continue

    dmel_ids = species_proteins.get(DMEL, [])
    dsim_ids = species_proteins.get(DSIM, [])
    dsec_ids = species_proteins.get(DSEC, [])

    trio_og_rows.append([
        og,
        ",".join(sorted(present_species)),
        len(dmel_ids),
        len(dsim_ids),
        len(dsec_ids),
        ",".join(dmel_ids),
        ",".join(dsim_ids),
        ",".join(dsec_ids),
    ])

    for protein in dmel_ids:
        trio_gene_rows.append([
            protein,
            og,
            ",".join(sorted(present_species)),
            ",".join(dsim_ids),
            ",".join(dsec_ids),
            int(bool(dsim_ids)),
            int(bool(dsec_ids)),
            int(bool(dsim_ids) and bool(dsec_ids)),
        ])

    # Strict subset: all three species represented.
    if {DMEL, DSIM, DSEC}.issubset(present_species):
        all_three_rows.append([
            og,
            ",".join(dmel_ids),
            ",".join(dsim_ids),
            ",".join(dsec_ids),
        ])


write_tsv(
    OUTDIR / "trio_restricted_orthogroups.tsv",
    [
        "orthogroup",
        "species_present",
        "Dmel_protein_count",
        "Dsim_protein_count",
        "Dsec_protein_count",
        "Dmel_protein_ids",
        "Dsim_protein_ids",
        "Dsec_protein_ids",
    ],
    trio_og_rows
)

write_tsv(
    OUTDIR / "Dmel_trio_restricted_genes.tsv",
    [
        "Dmel_protein_id",
        "orthogroup",
        "species_present",
        "Dsim_protein_ids",
        "Dsec_protein_ids",
        "present_in_Dsim",
        "present_in_Dsec",
        "present_in_all_three",
    ],
    trio_gene_rows
)

write_tsv(
    OUTDIR / "trio_shared_all_three.tsv",
    [
        "orthogroup",
        "Dmel_protein_ids",
        "Dsim_protein_ids",
        "Dsec_protein_ids",
    ],
    all_three_rows
)


# ============================================================
# 11. D. MELANOGASTER-ONLY ORTHOGROUPS
# ============================================================

print("Finding D. melanogaster-only orthogroups...")

dmel_only_rows = []

for og, species_proteins in og_data.items():

    present_species = set(species_proteins.keys())

    if present_species == {DMEL}:
        for protein in species_proteins.get(DMEL, []):
            dmel_only_rows.append([
                protein,
                og,
                len(species_proteins.get(DMEL, [])),
            ])


write_tsv(
    OUTDIR / "Dmel_only_orthogroup_genes.tsv",
    [
        "Dmel_protein_id",
        "orthogroup",
        "Dmel_proteins_in_orthogroup",
    ],
    dmel_only_rows
)


# ============================================================
# 12. DIAMOND: D. MELANOGASTER UNASSIGNED PROTEINS WITHOUT
#     CROSS-SPECIES HITS
# ============================================================

print("Reading DIAMOND summary...")

diamond_no_hits, diamond_known = read_diamond_no_hit_ids()

diamond_rows = []

for protein in sorted(unassigned_dmel):

    if protein in diamond_no_hits:
        status = "NO_CROSS_SPECIES_DIAMOND_HITS"
    elif protein in diamond_known:
        status = "CROSS_SPECIES_DIAMOND_HIT"
    else:
        status = "DIAMOND_STATUS_NOT_AVAILABLE"

    diamond_rows.append([protein, status])


write_tsv(
    OUTDIR / "Dmel_unassigned_diamond_status.tsv",
    ["Dmel_protein_id", "diamond_status"],
    diamond_rows
)


# ============================================================
# 13. COMBINED SPECIFICITY CANDIDATES
# ============================================================

specific_rows = []

# Candidates in Dmel-only orthogroups
for protein, og, n in dmel_only_rows:
    specific_rows.append([
        protein,
        "Dmel_only_orthogroup",
        og,
        "No other species represented in this orthogroup",
    ])

# Unassigned Dmel proteins with no cross-species DIAMOND hit
for protein in sorted(unassigned_dmel & diamond_no_hits):
    specific_rows.append([
        protein,
        "Unassigned_no_cross_species_DIAMOND_hit",
        "",
        "No cross-species hit in the targeted DIAMOND analysis",
    ])


write_tsv(
    OUTDIR / "Dmel_specific_candidates.tsv",
    [
        "Dmel_protein_id",
        "candidate_type",
        "orthogroup",
        "evidence",
    ],
    specific_rows
)


# ============================================================
# 14. SIMPLE ID LISTS
# ============================================================

def write_ids(filename, ids):
    with open(OUTDIR / filename, "w") as f:
        for protein in sorted(set(ids)):
            f.write(protein + "\n")


write_ids("Dmel_all_protein_ids.txt", dmel_proteins)

write_ids(
    "Dmel_trio_restricted_protein_ids.txt",
    [row[0] for row in trio_gene_rows]
)

write_ids(
    "Dmel_only_orthogroup_protein_ids.txt",
    [row[0] for row in dmel_only_rows]
)

write_ids(
    "Dmel_specific_candidate_ids.txt",
    [row[0] for row in specific_rows]
)


# ============================================================
# 15. FINAL REPORT
# ============================================================

print("\n========== FINISHED ==========")
print("D. melanogaster proteins:", len(dmel_proteins))

print("\nTechnical classification:")
for category, count in sorted(counts.items()):
    print(f"{category}: {count}")

print("\nDrosophila trio:")
print("Trio-restricted orthogroups:", len(trio_og_rows))
print("Dmel proteins in trio-restricted orthogroups:", len(trio_gene_rows))
print("Orthogroups represented in all three species:", len(all_three_rows))

print("\nD. melanogaster-only:")
print("Proteins in Dmel-only orthogroups:", len(dmel_only_rows))
print("Unassigned proteins:", len(unassigned_dmel))
print("Unassigned proteins without cross-species DIAMOND hits:",
      len(unassigned_dmel & diamond_no_hits))

print("\nOutput directory:", OUTDIR)
