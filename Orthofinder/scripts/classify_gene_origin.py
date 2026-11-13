import csv
import gzip
import os
import re
from collections import defaultdict, Counter

BASE = "."
OG_FILE = os.path.join(BASE, "Orthogroups", "Orthogroups.tsv")
UNASSIGNED_FILE = os.path.join(
    BASE, "Orthogroups", "Orthogroups_UnassignedGenes.tsv"
)
DUP_FILE = os.path.join(
    BASE, "Gene_Duplication_Events", "Duplications.tsv"
)
ORTHOLOG_DIR = os.path.join(BASE, "Orthologs")

# species -> protein -> information
genes = defaultdict(lambda: defaultdict(lambda: {
    "orthogroups": set(),
    "unassigned": False,
    "orthology": False,
    "partners": set(),
    "duplication": False,
    "duplication_events": set(),
}))

species_list = []
og_species = defaultdict(set)

def split_ids(value):
    if not value:
        return []
    return [x.strip() for x in value.split(",") if x.strip()]

def add_gene(species, protein):
    if species and protein:
        genes[species][protein]  # initialize record

def read_orthogroups():
    with open(OG_FILE, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        species = [x for x in reader.fieldnames if x != "Orthogroup"]
        for sp in species:
            if sp not in species_list:
                species_list.append(sp)

        for row in reader:
            og = row["Orthogroup"]
            for sp in species:
                for protein in split_ids(row.get(sp, "")):
                    add_gene(sp, protein)
                    genes[sp][protein]["orthogroups"].add(og)
                    og_species[og].add(sp)

def read_unassigned():
    with open(UNASSIGNED_FILE, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            for sp in species_list:
                for protein in split_ids(row.get(sp, "")):
                    add_gene(sp, protein)
                    genes[sp][protein]["unassigned"] = True

def read_orthologs():
    files = [
        x for x in os.listdir(ORTHOLOG_DIR)
        if x.endswith(".tsv.gz")
    ]

    for filename in files:
        source_species = filename[:-7]

        path = os.path.join(ORTHOLOG_DIR, filename)
        with gzip.open(path, "rt") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)

            for row in reader:
                if len(row) < 4:
                    continue

                target_species = row[1].strip()
                source_ids = split_ids(row[2])
                target_ids = split_ids(row[3])

                for protein in source_ids:
                    add_gene(source_species, protein)
                    genes[source_species][protein]["orthology"] = True
                    for partner in target_ids:
                        genes[source_species][protein]["partners"].add(
                            f"{target_species}:{partner}"
                        )

                for protein in target_ids:
                    add_gene(target_species, protein)
                    genes[target_species][protein]["orthology"] = True
                    for partner in source_ids:
                        genes[target_species][protein]["partners"].add(
                            f"{source_species}:{partner}"
                        )

def parse_prefixed_gene(value):
    """Разбирает вид_идентификатор белка по известным названиям видов."""
    value = value.strip()
    for sp in sorted(species_list, key=len, reverse=True):
        prefix = sp + "_"
        if value.startswith(prefix):
            return sp, value[len(prefix):]
    return None, None

def read_duplications():
    with open(DUP_FILE, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            og = row.get("Orthogroup", "")
            event = "|".join([
                og,
                row.get("Species Tree Node", ""),
                row.get("Gene Tree Node", ""),
                row.get("Type", ""),
            ])

            for column in ("Genes 1", "Genes 2"):
                for item in split_ids(row.get(column, "")):
                    sp, protein = parse_prefixed_gene(item)
                    if sp and protein:
                        add_gene(sp, protein)
                        genes[sp][protein]["duplication"] = True
                        genes[sp][protein]["duplication_events"].add(event)

def get_distribution(og):
    return sorted(og_species.get(og, set()))


def classify(info):
    if info["unassigned"]:
        return "UNASSIGNED_CANDIDATE"
    if info["orthology"] and info["duplication"]:
        return "ORTHOLOGS_AND_DUPLICATION"
    if info["orthology"]:
        return "ORTHOLOGS_DETECTED"
    if info["duplication"]:
        return "DUPLICATION_DETECTED"
    if info["orthogroups"]:
        return "HOMOLOGY_GROUP_ONLY"
    return "NO_ORTHOFINDER_GROUP"

def main():
    print("Reading orthogroups...")
    read_orthogroups()

    print("Reading unassigned proteins...")
    read_unassigned()

    print("Reading ortholog tables...")
    read_orthologs()

    print("Reading duplication events...")
    read_duplications()

    # Cache distribution for every orthogroup
    all_ogs = set()
    for sp in genes:
        for protein, info in genes[sp].items():
            all_ogs.update(info["orthogroups"])

    distributions = {
        og: get_distribution(og)
        for og in all_ogs
    }

    output_rows = []
    unassigned_rows = []

    for sp in sorted(genes):
        for protein, info in sorted(genes[sp].items()):
            ogs = sorted(info["orthogroups"])
            category = classify(info)

            distribution = sorted({
                species
                for og in ogs
                for species in distributions.get(og, [])
            })

            row = {
                "species": sp,
                "protein_id": protein,
                "orthogroups": ";".join(ogs),
                "orthology_detected": int(info["orthology"]),
                "ortholog_partner_count": len(info["partners"]),
                "duplication_detected": int(info["duplication"]),
                "duplication_event_count": len(info["duplication_events"]),
                "unassigned": int(info["unassigned"]),
                "orthogroup_species_count": len(distribution),
                "orthogroup_species": ";".join(distribution),
                "classification": category,
            }
            output_rows.append(row)

            if info["unassigned"]:
                unassigned_rows.append(row.copy())

    columns = [
        "species",
        "protein_id",
        "orthogroups",
        "orthology_detected",
        "ortholog_partner_count",
        "duplication_detected",
        "duplication_event_count",
        "unassigned",
        "orthogroup_species_count",
        "orthogroup_species",
        "classification",
    ]

    def write_table(path, rows):
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=columns, delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)

    write_table("gene_origin_classification.tsv", output_rows)
    write_table("unassigned_candidates.tsv", unassigned_rows)

    # Summary by species and classification
    with open("gene_origin_summary.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "species", "classification", "protein_count"
        ])

        for sp in sorted(genes):
            counts = Counter(
                classify(info) for info in genes[sp].values()
            )
            for category in sorted(counts):
                writer.writerow([sp, category, counts[category]])

    totals = Counter(
        classify(info)
        for sp in genes
        for info in genes[sp].values()
    )

    print("\n=== DONE ===")
    print("Proteins:", len(output_rows))
    print("Unassigned candidates:", len(unassigned_rows))
    print("\nClassification totals:")
    for category, count in sorted(totals.items()):
        print(f"{category}\t{count}")

    print("\nCreated:")
    print("  gene_origin_classification.tsv")
    print("  gene_origin_summary.tsv")
    print("  unassigned_candidates.tsv")

if __name__ == "__main__":
    main()
