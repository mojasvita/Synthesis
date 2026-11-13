#!/usr/bin/env python3
"""Step 1: FBgn -> BED, DNA from HAL, proteins from FlyBase. Python stdlib only."""
import argparse
import gzip
import re
import subprocess
from pathlib import Path


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def gene_coordinates(gff, gene_id):
    hits = []
    with open_text(gff) as source:
        for line in source:
            if line.startswith('#') or f'ID={gene_id}' not in line:
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) == 9 and fields[2] == 'gene':
                attrs = dict(item.split('=', 1) for item in fields[8].split(';') if '=' in item)
                if attrs.get('ID') == gene_id:
                    hits.append((fields[0], int(fields[3]) - 1, int(fields[4]), fields[6]))
    if len(hits) != 1:
        raise ValueError(f'Expected one GFF gene ID={gene_id}; found {len(hits)}')
    return hits[0]


def extract_proteins(source_path, gene_id, output_path):
    n = 0
    write = False
    # Match exact FBgn token in FlyBase parent field; preserve original FASTA headers.
    pattern = re.compile(r'parent=[^;\n]*?(?<!\w)' + re.escape(gene_id) + r'(?!\w)')
    with open_text(source_path) as source, output_path.open('w') as out:
        for line in source:
            if line.startswith('>'):
                write = bool(pattern.search(line))
                if write:
                    n += 1
            if write:
                out.write(line)
    if n == 0:
        output_path.unlink(missing_ok=True)
        raise ValueError(f'No FlyBase protein FASTA record for {gene_id}')
    return n


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gene_id')
    parser.add_argument('--hal', default='Acalyptratae20.hal')
    parser.add_argument('--reference', default='reference_data')
    parser.add_argument('--out', default='results')
    parser.add_argument('--hal-bin', default='~/bioinformatics/tools/cactus-bin-v3.3.0/bin')
    args = parser.parse_args()
    if not re.fullmatch(r'FBgn\d+', args.gene_id):
        parser.error('Expected FBgn ID')
    ref = Path(args.reference) / 'Drosophila_melanogaster'
    gff = ref / 'dmel-all-r6.36.gff.gz'
    proteins = ref / 'dmel-all-translation-r6.36.fasta.gz'
    output = Path(args.out) / args.gene_id
    output.mkdir(parents=True, exist_ok=True)
    chrom, start, end, strand = gene_coordinates(gff, args.gene_id)
    # Verify query interval lies within the HAL chromosome; one-base release
    # differences in GFF ##sequence-region are warnings, not hard failures.
    halstats = Path(args.hal_bin).expanduser() / 'halStats'
    sizes = dict((fields[0], int(fields[1])) for line in subprocess.check_output(
        [str(halstats), '--chromSizes', 'Drosophila_melanogaster', args.hal],
        text=True).splitlines() if len(fields := line.split()) == 2)
    if chrom not in sizes or end > sizes[chrom]:
        raise ValueError(f'Query interval {chrom}:{start}-{end} is outside HAL chromosome')
    with open_text(gff) as gff_source:
        for line in gff_source:
            if line.startswith('##sequence-region '):
                fields = line.split()
                if fields[1] == chrom and int(fields[3]) != sizes[chrom]:
                    print(f'WARNING: GFF chromosome length={fields[3]}, HAL={sizes[chrom]}; '
                          'check assembly compatibility (continuing because query fits HAL)')
                    break
            if not line.startswith('#'):
                break
    # BED is 0-based, end-exclusive; GFF is 1-based, end-inclusive.
    (output / 'query.bed').write_text(
        f'{chrom}\t{start}\t{end}\t{args.gene_id}\t0\t{strand}\n')
    n = extract_proteins(proteins, args.gene_id, output / 'query_proteins.faa')
    hal2fasta = Path(args.hal_bin).expanduser() / 'hal2fasta'
    with (output / 'query_dna.fna').open('w') as dna:
        subprocess.run([str(hal2fasta), args.hal, 'Drosophila_melanogaster',
                        '--sequence', chrom, '--start', str(start),
                        '--length', str(end - start)], check=True, stdout=dna)
    print(f'Query: {args.gene_id}, {chrom}:{start}-{end} ({strand}), BED coordinates')
    print(f'FlyBase proteins: {n}; saved in {output}')


if __name__ == '__main__':
    main()
