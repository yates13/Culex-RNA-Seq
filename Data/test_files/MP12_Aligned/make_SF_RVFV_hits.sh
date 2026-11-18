#!/bin/bash
# Extract and document RVFV-aligned reads from SF-RNA-D3-R3 (reads are NOT removed from the data)
# Run from MP12_Aligned/ with the hisat2_env_fixed conda env active.
set -euo pipefail

IN=trimmed.SF-RNA-D3-R3-Cxt_.sam
OUT=SF_RVFV_hits
command -v samtools >/dev/null || { echo "samtools not found - activate hisat2_env_fixed first"; exit 1; }
[ -s "$IN" ] || { echo "$IN not found - run this from MP12_Aligned/"; exit 1; }
mkdir -p "$OUT"

# 1. SAM of mapped records (with header)
samtools view -h -F 4 "$IN" > "$OUT/SF_RVFV_hits.sam"

# 2. TSV + annotated FASTA, with duplicate pairs grouped into fragment IDs
samtools view -F 4 "$IN" > "$OUT/.mapped.tmp"
python3 - "$OUT" <<'PY'
import sys, collections
out = sys.argv[1]
seg = {'NC_014395.1': 'S', 'NC_014396.1': 'M', 'NC_014397.1': 'L'}
pairs = collections.OrderedDict()
for line in open(f'{out}/.mapped.tmp'):
    f = line.rstrip('\r\n').split('\t')
    if len(f) < 10:
        continue
    pairs.setdefault(f[0], []).append(f)
# fragment = same segment + same mate start positions + same CIGARs
groups = collections.OrderedDict()
for name, rs in pairs.items():
    key = (rs[0][2],) + tuple(sorted(((r[3], r[5]) for r in rs), key=lambda x: int(x[0])))
    groups.setdefault(key, []).append(name)
fid = {}
for i, (key, names) in enumerate(sorted(groups.items(), key=lambda x: (-len(x[1]), x[0]))):
    for n in names:
        fid[n] = f'F{i+1:02d}'
with open(f'{out}/SF_RVFV_hits.tsv', 'w') as t, open(f'{out}/SF_RVFV_hits.fa', 'w') as fa:
    t.write('read_name\tmate\tflag\tsegment\trefseq\tpos\tmapq\tcigar\tfragment_id\tsequence\n')
    for name, rs in pairs.items():
        for r in rs:
            mate = 'R1' if int(r[1]) & 64 else 'R2'
            s = seg.get(r[2], r[2])
            t.write('\t'.join([name, mate, r[1], s, r[2], r[3], r[4], r[5], fid[name], r[9]]) + '\n')
            fa.write(f'>{name}/{mate[1]} {s}:{r[3]} {r[5]} {fid[name]}\n{r[9]}\n')
with open(f'{out}/SF_RVFV_fragments.tsv', 'w') as g:
    g.write('fragment_id\tsegment\tmate_starts\tcigars\tn_pairs\ttiles\n')
    for key, names in sorted(groups.items(), key=lambda x: (-len(x[1]), x[0])):
        g.write('\t'.join([fid[names[0]], seg.get(key[0], key[0]),
                           '/'.join(p for p, c in key[1:]), '/'.join(c for p, c in key[1:]),
                           str(len(names)), ','.join(n.split(':')[4] for n in names)]) + '\n')
print(f'{len(pairs)} pairs, {sum(len(v) for v in pairs.values())} records, {len(groups)} distinct fragments')
PY
rm "$OUT/.mapped.tmp"

# 3. README write-up
cat > "$OUT/README.md" <<'EOF'
# RVFV-aligned reads in sample SF-RNA-D3-R3-Cxt

**Location:** `/scratch/alpine/c832500103@colostate.edu/RNA_Seq/Data/test_files/MP12_Aligned/SF_RVFV_hits/`
**Source alignment:** `../trimmed.SF-RNA-D3-R3-Cxt_.sam` (HISAT2, trimmed reads vs. RVFV reference)
**Date checked:** 2026-11-18

## Summary

Sample SF-RNA-D3-R3 is not expected to contain Rift Valley fever virus, but **22 read pairs (44 alignment records)** in it align to the RVFV reference. The reads have **not been removed** from the data. They are extracted here for documentation.

The alignments are high quality, so these are real RVFV sequences, not spurious alignments of mosquito reads. They collapse to only **12 distinct fragments** because of PCR duplication. That pattern fits low-level carry-over from an RVFV-positive library (index hopping or cross-contamination) better than a true infection. This is a working interpretation and has not been confirmed yet (see "Checks still to do").

## Files

| File | Contents |
|---|---|
| `SF_RVFV_hits.sam` | The 44 mapped records, with the SAM header (`samtools view -h -F 4`) |
| `SF_RVFV_hits.tsv` | One row per record: read name, mate, flag, segment, RefSeq, position, MAPQ, CIGAR, fragment ID, sequence |
| `SF_RVFV_hits.fa` | The same reads as FASTA. The headers give the segment, position, CIGAR and fragment ID (BLAST-ready) |
| `SF_RVFV_fragments.tsv` | One row per distinct fragment, with the number of duplicate pairs and the tiles they came from |
| `../make_SF_RVFV_hits.sh` | The script that generated all of the above |

Reference segments: `NC_014395.1` = S, `NC_014396.1` = M, `NC_014397.1` = L (RVFV RefSeq, ZH-548 lineage, the parent strain of MP-12).

## What the alignments look like

- **All 44 records have MAPQ 60, and every CIGAR is a full match** (112M to 150M) with no soft-clipping or indels.
- **All 22 pairs are properly paired in the correct orientation** (flags 99/147 or 83/163).
- **All three genome segments are hit**: S in 7 pairs, M in 7 pairs and L in 8 pairs. Spurious alignments of host or low-complexity reads would not normally produce clean, full-length, properly paired hits spread across S, M and L.
- **All reads come from the same flowcell and lane** (`LH00401:707:23WW2FLT4`, lane 6).

## Fragments (duplicates collapsed)

Pairs with the same segment, the same mate start positions and the same sequence were grouped as one fragment. The duplicates come from different tiles, including tiles on both flowcell surfaces (1xxx and 2xxx), so they are PCR duplicates rather than optical duplicates.

| Fragment | Segment | Mate starts | Read length | Pairs | Tiles |
|---|---|---|---|---|---|
| F01 | L | 5187 / 5206 | 150 | 5 | 1180, 1386, 1464, 2180, 2390 |
| F02 | L | 4969 / 4969 | 129 | 4 | 1311, 1493, 2266, 2294 |
| F03 | S | 229 / 229 | 122 | 2 | 2164, 2311 |
| F04 | S | 314 / 392 | 150 | 2 | 1215, 1404 |
| F05 | M | 992 / 1004 | 150 | 2 | 1343, 2263 * |
| F06 | S | 522 / 576 | 150 | 1 | 2220 |
| F07 | S | 532 / 532 | 112 | 1 | 2196 |
| F08 | M | 2430 / 2430 | 142 | 1 | 1294 |
| F09 | M | 2585 / 2662 | 150 | 1 | 2202 |
| F10 | M | 3503 / 3503 | 146 | 1 | 1439 |
| F11 | L | 1559 / 1559 | 126 | 1 | 2288 |
| F12 | L | 2440 / 2491 | 150 | 1 | 1187 |

\* In the F05 copy from tile 2263, R1 differs from the other copy by one base (C→G near the 3′ end). This is consistent with a sequencing or PCR error.

## Interpretation

1. **The reads are real RVFV sequence**, as the alignment quality above shows.
2. **The sample is not expected to be infected, and the signal is very low.** There are 22 pairs, only 12 of them unique, against the sample's full library depth. A real infection would normally produce many more independent fragments spread along each segment.
3. **The heavy duplication of a few fragments** (F01 ×5, F02 ×4) suggests that these molecules were abundant in a *source* library, and that a small number of them were mis-assigned to SF. The most likely route is index hopping or barcode bleed on the shared lane. Low-level cross-contamination during library preparation is also possible.
4. **The reads were left in the dataset.** At this level they have no meaningful effect on host (Culex) gene quantification. They define the background level of RVFV signal to expect in uninfected samples on this run.

## Checks still to do

Run these from `MP12_Aligned/`:

```bash
# Is the RVFV-positive MP library on the same flowcell/lane?
grep -v '^@' trimmed.MP-RNA-D1-R2-Cxt_.sam | head -1 | cut -f1

# Do the exact SF fragment sequences also appear in the MP library?
tail -n +2 SF_RVFV_hits/SF_RVFV_hits.tsv | cut -f10 | sort -u > SF_RVFV_hits/SF_RVFV_seqs.txt
samtools view -F 4 trimmed.MP-RNA-D1-R2-Cxt_.sam | cut -f10 | grep -c -F -x -f SF_RVFV_hits/SF_RVFV_seqs.txt

# How many RVFV reads does MP have? (to estimate the bleed rate: SF hits / MP hits)
samtools view -c -F 0x904 trimmed.MP-RNA-D1-R2-Cxt_.sam
```

If the MP library shares lane 6 and contains these identical fragments, and SF/MP is in the typical index-hopping range (~0.1–1%), carry-over is confirmed. For an independent check of identity, the reads in `SF_RVFV_hits.fa` can be run through NCBI blastn.
EOF

echo "Done. Files in $OUT/:"
ls -l "$OUT"
