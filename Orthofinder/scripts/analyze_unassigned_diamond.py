import csv
from collections import defaultdict

HITS_FILE = "unassigned_diamond_hits.tsv"
QUERIES_FILE = "unassigned_queries.faa"

# Таксономические группы для первичного анализа.
def taxonomic_group(species):
    if species.startswith("D_"):
        return "Drosophila"
    if species == "Scaptodrosophila_lebanonensis":
        return "Scaptodrosophila"
    if species in {
        "Aedes_aegypti",
        "Anopheles_gambiae",
        "Bactrocera_dorsalis",
        "Bactrocera_latrifrons",
        "Musca_domestica",
    }:
        return "Other_Diptera"
    if species in {
        "Apis_florea",
        "Apis_mellifera",
        "Bombyx_mori",
        "Monomorium_pharaonis",
    }:
        return "Other_insects"
    if species == "Daphnia_magna":
        return "Daphnia_outgroup"
    return "Other"


def parse_id(full_id):
    if "|" not in full_id:
        return "UNKNOWN", full_id
    return full_id.split("|", 1)


# Считываем все запросы, чтобы не потерять те, для которых нет совпадений.
query_ids = []
with open(QUERIES_FILE) as f:
    for line in f:
        if line.startswith(">"):
            query_ids.append(line[1:].strip().split()[0])

# Считываем DIAMOND.
hits_by_query = defaultdict(list)
all_hits = []

with open(HITS_FILE) as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 8:
            continue

        qid, sid = parts[0], parts[1]
        pident = float(parts[2])
        length = int(parts[3])
        qlen = int(parts[4])
        slen = int(parts[5])
        evalue = float(parts[6])
        bitscore = float(parts[7])

        qspecies, qprotein = parse_id(qid)
        sspecies, sprotein = parse_id(sid)

        # Исключаем совпадение белка с самим собой.
        if qid == sid:
            continue

        qcov = 100 * length / qlen if qlen else 0
        scov = 100 * length / slen if slen else 0

        record = {
            "query": qid,
            "query_species": qspecies,
            "query_protein": qprotein,
            "target": sid,
            "target_species": sspecies,
            "target_protein": sprotein,
            "group": taxonomic_group(sspecies),
            "identity": pident,
            "alignment_length": length,
            "query_length": qlen,
            "target_length": slen,
            "query_coverage": qcov,
            "target_coverage": scov,
            "evalue": evalue,
            "bitscore": bitscore,
        }

        hits_by_query[qid].append(record)
        all_hits.append(record)


# 1. Подробная таблица межвидовых совпадений.
hit_columns = [
    "query", "query_species", "query_protein",
    "target", "target_species", "target_protein",
    "group", "identity", "alignment_length",
    "query_length", "target_length",
    "query_coverage", "target_coverage",
    "evalue", "bitscore"
]

with open("unassigned_cross_species_hits.tsv", "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=hit_columns, delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(all_hits)


# 2. Сводка по каждому кандидату.
summary_columns = [
    "query", "query_species", "query_protein",
    "cross_species_hit_count",
    "target_species_count",
    "taxonomic_groups",
    "best_hit_species",
    "best_hit_group",
    "best_hit_identity",
    "best_hit_query_coverage",
    "best_hit_target_coverage",
    "best_hit_evalue",
    "best_hit_bitscore",
    "best_hit_protein"
]

summary = []

for qid in query_ids:
    qspecies, qprotein = parse_id(qid)
    hits = hits_by_query.get(qid, [])

    # Лучшее совпадение по bitscore.
    best = max(hits, key=lambda x: x["bitscore"]) if hits else None

    target_species = sorted({
        h["target_species"] for h in hits
    })

    groups = sorted({
        h["group"] for h in hits
    })

    row = {
        "query": qid,
        "query_species": qspecies,
        "query_protein": qprotein,
        "cross_species_hit_count": len(hits),
        "target_species_count": len(target_species),
        "taxonomic_groups": ";".join(groups),
        "best_hit_species": best["target_species"] if best else "",
        "best_hit_group": best["group"] if best else "",
        "best_hit_identity": best["identity"] if best else "",
        "best_hit_query_coverage": best["query_coverage"] if best else "",
        "best_hit_target_coverage": best["target_coverage"] if best else "",
        "best_hit_evalue": best["evalue"] if best else "",
        "best_hit_bitscore": best["bitscore"] if best else "",
        "best_hit_protein": best["target_protein"] if best else "",
    }
    summary.append(row)

with open("unassigned_diamond_summary.tsv", "w", newline="") as f:
    writer = csv.DictWriter(
        f, fieldnames=summary_columns, delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(summary)


# 3. Запросы без межвидовых совпадений.
no_cross_hits = [
    qid for qid in query_ids
    if len(hits_by_query.get(qid, [])) == 0
]

with open("unassigned_no_cross_species_hits.txt", "w") as f:
    f.write("\n".join(no_cross_hits))


# 4. Общая статистика.
print("=== DIAMOND ANALYSIS ===")
print("Queries in FASTA:", len(query_ids))
print("Queries with cross-species hits:", sum(
    1 for q in query_ids if hits_by_query.get(q)
))
print("Queries without cross-species hits:", len(no_cross_hits))
print("Cross-species alignments:", len(all_hits))

print("\n=== TAXONOMIC GROUPS ===")
group_counts = defaultdict(int)
for hit in all_hits:
    group_counts[hit["group"]] += 1

for group, count in sorted(group_counts.items()):
    print(f"{group}: {count}")

print("\nCreated:")
print("  unassigned_cross_species_hits.tsv")
print("  unassigned_diamond_summary.tsv")
print("  unassigned_no_cross_species_hits.txt")
