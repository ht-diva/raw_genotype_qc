"""Rename BIM IDs without changing row order or allele order; fail on collisions."""
import csv
import sys
from collections import Counter
from pathlib import Path


def main(bim, chromosome_map, output, mapping):
    chrom = dict(line.split() for line in Path(chromosome_map).read_text().splitlines() if line.strip())
    rows = [line.split() for line in Path(bim).read_text().splitlines() if line.strip()]
    if not rows or any(len(row) != 6 for row in rows):
        raise ValueError('Expected a nonempty, six-column BIM')
    ids = [f'{chrom[r[0]]}:{int(r[3])}:{r[5]}:{r[4]}' for r in rows]
    counts = Counter(ids)
    Path(mapping).parent.mkdir(parents=True, exist_ok=True)
    with open(mapping, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['original_row', 'original_id', 'original_chr', 'chr', 'pos', 'provisional_ref', 'alt', 'working_id', 'id_count'])
        for i, (row, vid) in enumerate(zip(rows, ids), 1):
            writer.writerow([i, row[1], row[0], chrom[row[0]], row[3], row[5], row[4], vid, counts[vid]])
    duplicates = sum(n > 1 for n in counts.values())
    if duplicates:
        raise ValueError(f'{duplicates} coordinate/allele ID collisions; inspect {mapping}. Resolve duplicates before rerunning; none were silently dropped.')
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, 'w') as handle:
        for row, vid in zip(rows, ids):
            row[1] = vid
            handle.write('\t'.join(row) + '\n')


if __name__ == '__main__':
    main(*sys.argv[1:])