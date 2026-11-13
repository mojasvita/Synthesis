
import csv
from pathlib import Path
from Bio import SeqIO

# =====================================================
# PATHS
# =====================================================

BASE = Path("/home/e.nikolaeva/orfan_gene/drosophila_proteomes/proteomes_1isoform")
OUT = BASE / "Dmel_origin_analysis"

FASTA_FILE = BASE / "D_melanogaster.faa"

TRIO_FILE = OUT / "trio_restricted.tsv"
DIAMOND_FILE = OUT / "Dmel_unassigned_diamond_status.tsv"


# =====================================================
# READ PROTEIN SEQUENCES
# =====================================================

print("Reading D. melanogaster proteome...")

sequences = {}

for record in SeqIO.parse(FASTA_FILE, "fasta"):
    sequences[record.id] = record

print(f"Proteins in original FASTA: {len(sequences)}")


# =====================================================
# READ TRIO-RESTRICTED IDS
# =====================================================

print("\nReading trio-restricted IDs...")

trio_ids = set()

with open(TRIO_FILE) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        ids = row["D_melanogaster_protein_ids"].strip()

        if ids:
            for protein_id in ids.split(","):
                protein_id = protein_id.strip()
                if protein_id:
                    trio_ids.add(protein_id)

print(f"Unique trio-restricted Dmel IDs: {len(trio_ids)}")


# =====================================================
# READ UNASSIGNED DIAMOND STATUS IDS
# =====================================================

print("\nReading unassigned DIAMOND status IDs...")

diamond_ids = set()

with open(DIAMOND_FILE) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        protein_id = row["Dmel_protein_id"].strip()

        if protein_id:
            diamond_ids.add(protein_id)

print(f"Unique unassigned Dmel IDs: {len(diamond_ids)}")


# =====================================================
# EXTRACT SEQUENCES
# =====================================================

def extract_sequences(ids, output_file):
    found = []
    missing = []

    for protein_id in sorted(ids):
        if protein_id in sequences:
            found.append(sequences[protein_id])
        else:
            missing.append(protein_id)

    SeqIO.write(found, output_file, "fasta")

    return len(found), missing


trio_output = OUT / "Dmel_trio_restricted.faa"
diamond_output = OUT / "Dmel_unassigned_diamond_status.faa"

trio_found, trio_missing = extract_sequences(trio_ids, trio_output)
diamond_found, diamond_missing = extract_sequences(diamond_ids, diamond_output)


# =====================================================
# SAVE MISSING IDS
# =====================================================

with open(OUT / "trio_missing_ids.txt", "w") as f:
    f.write("\n".join(sorted(trio_missing)))

with open(OUT / "diamond_missing_ids.txt", "w") as f:
    f.write("\n".join(sorted(diamond_missing)))


# =====================================================
# CHECK OVERLAP
# =====================================================

overlap = trio_ids & diamond_ids

# =====================================================
# SUMMARY
# =====================================================

summary = [
    {
        "dataset": "trio_restricted",
        "input_unique_ids": len(trio_ids),
        "sequences_found": trio_found,
        "ids_missing": len(trio_missing),
        "output_file": trio_output.name
    },
    {
        "dataset": "unassigned_diamond_status",
        "input_unique_ids": len(diamond_ids),
        "sequences_found": diamond_found,
        "ids_missing": len(diamond_missing),
        "output_file": diamond_output.name
    }
]

summary_file = OUT / "sequence_extraction_summary.tsv"

with open(summary_file, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=summary[0].keys(),
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(summary)


# =====================================================
# REPORT
# =====================================================

print("\n========== EXTRACTION FINISHED ==========")

print("\nTRIO RESTRICTED:")
print(f"  Input IDs:       {len(trio_ids)}")
print(f"  Sequences found: {trio_found}")
print(f"  Missing IDs:     {len(trio_missing)}")
print(f"  Output:          {trio_output}")

print("\nUNASSIGNED DIAMOND STATUS:")
print(f"  Input IDs:       {len(diamond_ids)}")
print(f"  Sequences found: {diamond_found}")
print(f"  Missing IDs:     {len(diamond_missing)}")
print(f"  Output:          {diamond_output}")

print(f"\nIDs present in both datasets: {len(overlap)}")

print("\nSummary:", summary_file)
