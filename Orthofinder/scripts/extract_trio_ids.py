
#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

BASE = Path("/home/e.nikolaeva/orfan_gene/drosophila_proteomes/proteomes_1isoform")
RESULTS = BASE / "OrthoFinder/Results_Nov12_1"
OUT = BASE / "Dmel_origin_analysis"
OUT.mkdir(exist_ok=True)

DMEL = "D_melanogaster"
DSIM = "D_simulans"
DSEC = "D_sechellia"

TRIO = {DMEL, DSIM, DSEC}

INPUT = RESULTS / "Orthogroups/Orthogroups.tsv"


def parse_ids(value):
    if not value or value.strip() in {"-", "NA"}:
        return []
    return [x.strip() for x in value.split(",") if x.strip()]


def save_tsv(filename, header, rows):
    with open(OUT / filename, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def save_ids(filename, ids):
    with open(OUT / filename, "w") as f:
        for protein in sorted(ids):
            f.write(protein + "\n")


# Read orthogroups
print("Reading orthogroups...")

groups = []

with open(INPUT) as f:
    reader = csv.DictReader(f, delimiter="\t")

    species_names = [
        x for x in reader.fieldnames if x != "Orthogroup"
    ]

    drosophila_species = [
        x for x in species_names if x.startswith("D_")
    ]

    for row in reader:
        data = {
            species: parse_ids(row.get(species, ""))
            for species in species_names
        }

        present = {
            species for species, ids in data.items() if ids
        }

        groups.append({
            "og": row["Orthogroup"],
            "data": data,
            "present": present
        })


# Define categories
dmel_groups = []
dmel_dsim = []
dmel_dsec = []
all_three = []
trio_restricted = []
drosophila_restricted = []

for g in groups:
    present = g["present"]

    if DMEL not in present:
        continue

    dmel_groups.append(g)

    if DSIM in present:
        dmel_dsim.append(g)

    if DSEC in present:
        dmel_dsec.append(g)

    if TRIO.issubset(present):
        all_three.append(g)

    if TRIO.issubset(present) and present.issubset(TRIO):
        trio_restricted.append(g)

    # Only Drosophila genus species, no other genera
    if present.issubset(set(drosophila_species)):
        drosophila_restricted.append(g)


# Write detailed orthogroup tables
def make_rows(group_list, species_list):
    rows = []

    for g in group_list:
        data = g["data"]

        rows.append([
            g["og"],
            len(g["present"]),
            ",".join(sorted(g["present"])),
            *[
                len(data[s]) for s in species_list
            ],
            *[
                ",".join(data[s]) for s in species_list
            ]
        ])

    return rows


trio_header = [
    "orthogroup",
    "species_count",
    "species_present"
]

trio_header += [s + "_protein_count" for s in [DMEL, DSIM, DSEC]]
trio_header += [s + "_protein_ids" for s in [DMEL, DSIM, DSEC]]

save_tsv(
    "trio_all_details.tsv",
    trio_header,
    make_rows(
        [
            g for g in groups
            if DMEL in g["present"]
            and bool(g["present"] & {DSIM, DSEC})
        ],
        [DMEL, DSIM, DSEC]
    )
)

save_tsv(
    "trio_all_shared.tsv",
    trio_header,
    make_rows(all_three, [DMEL, DSIM, DSEC])
)

save_tsv(
    "trio_restricted.tsv",
    trio_header,
    make_rows(trio_restricted, [DMEL, DSIM, DSEC])
)


# Drosophila genus-restricted groups
drosophila_header = [
    "orthogroup",
    "species_count",
    "species_present",
]

drosophila_header += [
    s + "_protein_count" for s in drosophila_species
]

drosophila_header += [
    s + "_protein_ids" for s in drosophila_species
]

save_tsv(
    "drosophila_genus_restricted.tsv",
    drosophila_header,
    make_rows(drosophila_restricted, drosophila_species)
)


# Collect unique protein IDs per species and category
def collect_ids(group_list, species):
    result = set()

    for g in group_list:
        result.update(g["data"][species])

    return result


for species in [DMEL, DSIM, DSEC]:
    save_ids(
        f"{species}_all_three_ids.txt",
        collect_ids(all_three, species)
    )

    save_ids(
        f"{species}_trio_restricted_ids.txt",
        collect_ids(trio_restricted, species)
    )

    save_ids(
        f"{species}_all_dmel_associated_ids.txt",
        collect_ids(dmel_groups, species)
    )


# Summary
summary_rows = []

categories = {
    "Dmel_all_orthogroups": dmel_groups,
    "Dmel_plus_Dsim": dmel_dsim,
    "Dmel_plus_Dsec": dmel_dsec,
    "All_three_species": all_three,
    "Trio_restricted": trio_restricted,
    "Drosophila_genus_restricted": drosophila_restricted,
}

for category, group_list in categories.items():
    row = [
        category,
        len(group_list)
    ]

    for species in [DMEL, DSIM, DSEC]:
        row.append(
            len(collect_ids(group_list, species))
        )

    summary_rows.append(row)


save_tsv(
    "trio_summary.tsv",
    [
        "category",
        "orthogroup_count",
        "Dmel_unique_proteins",
        "Dsim_unique_proteins",
        "Dsec_unique_proteins"
    ],
    summary_rows
)


# Final report
print("\n========== FINISHED ==========")

for category, group_list in categories.items():
    print(f"{category}: {len(group_list)} orthogroups")

    for species in [DMEL, DSIM, DSEC]:
        print(
            f"  {species}: "
            f"{len(collect_ids(group_list, species))} unique proteins"
        )

print("\nDrosophila species in dataset:")
print(", ".join(drosophila_species))

print("\nOutput:", OUT)
