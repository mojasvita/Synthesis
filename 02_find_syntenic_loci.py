#!/usr/bin/env python3
"""Step 2: HAL liftover -> BED fragments -> coding-annotation overlap per locus.

This script does NOT infer non-genic hits. It only records coding annotation overlap.
"""
import argparse
import csv
import gzip
import subprocess
from pathlib import Path

COLS = ['gene_id', 'species', 'locus_id', 'chromosome', 'start', 'end',
        'aligned_bp', 'coding_gene_ids', 'coding_gene_tags', 'overlap_bp',
        'status', 'note']


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def read_bed(path):
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split('\t')
        if len(parts) >= 3:
            try:
                a, b = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            if b > a:
                rows.append((parts[0], a, b, parts[5] if len(parts) > 5 else '.'))
    return sorted(rows)


def merged_segments(intervals):
    merged = []
    for a, b in sorted(intervals):
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    return merged


def group_bed(fragments, max_gap):
    groups = []
    for chrom, a, b, strand in fragments:
        if not groups or groups[-1]['chrom'] != chrom or groups[-1]['strand'] != strand or a - groups[-1]['end'] > max_gap:
            groups.append(dict(chrom=chrom, start=a, end=b, strand=strand, segments=[(a, b)]))
        else:
            groups[-1]['end'] = max(groups[-1]['end'], b)
            groups[-1]['segments'].append((a, b))
    return groups


def protein_coding_genes(gff, chromosomes):
    genes = []
    with open_text(gff) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9 or parts[2] != 'gene' or parts[0] not in chromosomes:
                continue
            attrs = dict(x.split('=', 1) for x in parts[8].split(';') if '=' in x)
            if attrs.get('gene_biotype') == 'protein_coding':
                genes.append((parts[0], int(parts[3]) - 1, int(parts[4]),
                              attrs.get('ID', ''), attrs.get('locus_tag', '')))
    return genes


def overlap(segments, start, end):
    return sum(max(0, min(b, end) - max(a, start)) for a, b in merged_segments(segments))


def find_gff(root):
    base = root / 'ncbi_dataset' / 'data'
    for name in ('genomic.gff', 'genomic.gff3', 'genomic.gff.gz', 'genomic.gff3.gz'):
        paths = sorted(base.rglob(name)) if base.is_dir() else []
        if paths:
            return paths[0]
    return None


def seq_sizes_from_gff(gff):
    """Read ALL sequence-region declarations; NCBI may interleave them with features."""
    sizes = {}
    with open_text(gff) as fh:
        for line in fh:
            if line.startswith('##sequence-region '):
                cols = line.split()
                if len(cols) >= 4:
                    sizes[cols[1]] = int(cols[3])
    return sizes


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('gene_id')
    p.add_argument('--hal', default='Acalyptratae20.hal')
    p.add_argument('--assemblies', default='supplementary_file_genome_assembly.csv')
    p.add_argument('--reference', default='reference_data')
    p.add_argument('--out', default='results')
    p.add_argument('--hal-bin', default='~/bioinformatics/tools/cactus-bin-v3.3.0/bin')
    p.add_argument('--species', nargs='*', help='simulans yakuba erecta (default: 19 author Cactus outgroups)')
    p.add_argument('--max-gap', type=int, default=2000, help='Group nearby HAL fragments (heuristic)')
    args = p.parse_args()
    out = Path(args.out) / args.gene_id
    query_bed = out / 'query.bed'
    if not query_bed.is_file():
        p.error(f'Missing {query_bed}; run 01_prepare_query.py first')
    with open(args.assemblies, newline='', encoding='utf-8-sig') as f:
        all_rows = list(csv.DictReader(f))[:20]  # 20 species in original Cactus subset, including Dmel
    species = [r['Species'].replace(' ', '_') for r in all_rows if r['Species'] != 'Drosophila melanogaster']
    if args.species:
        wanted = {s.lower() for s in args.species}
        species = [s for s in species if s.lower() in wanted or s.lower().split('_')[-1] in wanted]
        if not species:
            p.error('No matching species in CSV')
    halbin = Path(args.hal_bin).expanduser()
    liftover = str(halbin / 'halLiftover')
    halstats = str(halbin / 'halStats')
    hal2fasta = str(halbin / 'hal2fasta')
    rows = []
    for org in species:
        print(f'[{org}]', flush=True)
        raw = out / f'{org}_liftover.bed'
        try:
            subprocess.run([liftover, args.hal, 'Drosophila_melanogaster',
                            str(query_bed), org, str(raw)], check=True, capture_output=True, text=True)
            loci = group_bed(read_bed(raw), args.max_gap)
            if not loci:
                rows.append(dict(gene_id=args.gene_id, species=org, status='unaligned',
                                 note='No aligned HAL segments'))
                continue
            gff = find_gff(Path(args.reference) / org)
            gff_sizes = seq_sizes_from_gff(gff) if gff else {}
            hal_out = subprocess.check_output([halstats, '--chromSizes', org, args.hal], text=True)
            hal_sizes = dict((v[0], int(v[1])) for l in hal_out.splitlines() if len(v := l.split()) == 2)
            genes = protein_coding_genes(gff, {loc['chrom'] for loc in loci}) if gff else []
            for idx, loc in enumerate(loci, 1):
                row = {k: '' for k in COLS}
                row.update(gene_id=args.gene_id, species=org, locus_id=f'{org}_L{idx}',
                           chromosome=loc['chrom'], start=loc['start'], end=loc['end'],
                           aligned_bp=sum(b-a for a, b in merged_segments(loc['segments'])))
                chrom = loc['chrom']
                if chrom not in hal_sizes or loc['end'] > hal_sizes[chrom]:
                    row.update(status='invalid_hal_coordinates', note='Chromosome absent or interval outside HAL')
                elif not gff:
                    row.update(status='annotation_missing', note='GFF3 not downloaded')
                elif chrom not in gff_sizes:
                    row.update(status='assembly_unverified', note='Chromosome not found in GFF sequence-region headers')
                elif abs(gff_sizes[chrom] - hal_sizes[chrom]) > 1:
                    row.update(status='assembly_mismatch', note=f'GFF size={gff_sizes[chrom]}, HAL size={hal_sizes[chrom]}')
                else:
                    hits = []
                    for gc, ga, gb, gid, tag in genes:
                        if gc == chrom:
                            bp = overlap(loc['segments'], ga, gb)
                            if bp:
                                hits.append((bp, gid, tag))
                    hits.sort(reverse=True)
                    row['coding_gene_ids'] = ','.join(h[1] for h in hits)
                    row['coding_gene_tags'] = ','.join(h[2] for h in hits)
                    row['overlap_bp'] = sum(h[0] for h in hits)
                    row['status'] = 'annotated_locus_needs_blastp' if hits else 'no_coding_annotation_needs_genewise'
                    if gff_sizes[chrom] != hal_sizes[chrom]:
                        row['note'] = f'WARNING: GFF/HAL length differs by 1 bp ({gff_sizes[chrom]}/{hal_sizes[chrom]}); verify assembly'
                # Extract DNA even when annotation unavailable, if HAL coordinates valid.
                if chrom in hal_sizes and loc['end'] <= hal_sizes[chrom]:
                    dna_path = out / f'{org}_L{idx}_dna.fna'
                    with dna_path.open('w') as stream:
                        subprocess.run([hal2fasta, args.hal, org, '--sequence', chrom,
                                        '--start', str(loc['start']), '--length', str(loc['end'] - loc['start'])],
                                       check=True, stdout=stream)
                rows.append(row)
                print(' ', row['locus_id'], row['status'], row['coding_gene_ids'])
        except (OSError, subprocess.CalledProcessError, ValueError) as exc:
            rows.append(dict(gene_id=args.gene_id, species=org, status='error', note=str(exc)[:250]))
            print('  ERROR:', exc)
    dest = out / 'syntenic_loci.tsv'
    with dest.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, COLS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    print('Saved:', dest)


if __name__ == '__main__':
    main()
