import csv
import os
from collections import defaultdict

PROTEOME_DIR = (
    "/home/e.nikolaeva/orfan_gene/"
    "drosophila_proteomes/proteomes_1isoform"
)

CANDIDATES = "unassigned_candidates.tsv"

def read_fasta(path):
    header = None
    seq = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header.split()[0], "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)

        if header is not None:
            yield header.split()[0], "".join(seq)

# Собираем ID кандидатов по видам
candidate_ids = defaultdict(set)

with open(CANDIDATES) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        candidate_ids[row["species"]].add(row["protein_id"])

queries = []
database = []

for filename in os.listdir(PROTEOME_DIR):
    if not filename.endswith(".faa"):
        continue

    species = filename[:-4]
    path = os.path.join(PROTEOME_DIR, filename)

    print("Reading:", species)

    for protein_id, sequence in read_fasta(path):
        if not sequence:
            continue

        full_id = f"{species}|{protein_id}"

        # Все белки идут в базу
        database.append((full_id, sequence))

        # Нераспределённые белки идут в запрос
        if protein_id in candidate_ids[species]:
            queries.append((full_id, sequence))

def write_fasta(path, records):
    with open(path, "w") as f:
        for name, sequence in records:
            f.write(f">{name}\n")
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + "\n")

write_fasta("unassigned_queries.faa", queries)
write_fasta("all_proteomes_db.faa", database)

print("\nQueries:", len(queries))
print("Database proteins:", len(database))
print("Created unassigned_queries.faa")