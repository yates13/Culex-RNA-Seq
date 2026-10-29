# Culex tarsalis RNA-seq: host / RVFV / unaligned-reads pipeline

RNA-seq processing for *Culex tarsalis* (CtarK1) on CSU/CURC Alpine (SLURM).
The main goal of the HISAT2 stage is to split every library into three datasets:

1. **Reads that aligned to the Culex genome** (host)
2. **Reads that aligned to Rift Valley fever virus (RVFV)**
3. **Reads that aligned to neither**, which go on to the ssRNA+ alignment and
   BLAST against insect-specific viruses

Aligning to the host first and keeping what does not align helps rule out
endogenous viral elements (EVEs) integrated into the mosquito genome.

## Directory layout

```
RNA_Seq/
  Raw_Files/        raw FASTQs from the sequencing core
  Genome/
    Index_Genome/hisatIndex/   Cxt genome HISAT2 index (prefix: CtarK1_index)
    Index_RVFV/                RVFV HISAT2 index (prefix: RVFV_index)
  Scripts/          this pipeline (see below)
  Data/             all outputs, one folder per step
  slurmlogs/        outlog/ and errlog/ (must exist before submitting)
```

## Scripts

| Script | Step | Notes |
|---|---|---|
| `01b_runFastp_RNA_HTML.sh` | Trim reads with fastp, HTML QC reports | Output FASTQs are named `trimmed.<sample>_R1_001.fastq.gz` / `_R2_001.fastq.gz` |
| `02_RNA_Cxthisatarray_unaligned.sh` | Two-stage HISAT2 array: Cxt genome, then RVFV | Produces the three datasets above. **Current script.** |
| `02_RNA_AedesHisatArray.sh` | Original single-stage HISAT2 array | Legacy name; it aligned Culex data. Replaced by the script above. |
| `02b_RerunFailed.sh` | Re-run failed array tasks | *(confirm description)* |
| `04_GTFnew_runAeHTseqArray.sh` | Gene counts with HTSeq against the updated GTF | *(confirm description)* |
| `05_GTF_MosquitoRescue_v2.sh` | GTF-related step | *(confirm description)* |
| `DiffBindFunction_ME_d3_BF.R` | CUT&RUN differential binding (not part of the RNA-seq flow) | |
| `RVFV_IndexNotes.txt` | How the RVFV index was built | Same steps as below |
| `old_scripts/` | Superseded versions | |

## Workflow

1. **Trim** all libraries (`01b`). Make sure every library, including the remaining SF
   libraries, is trimmed before building the sample list.
2. **Direct-alignment stats** from the earlier alignments are collected first, then the
   old `.sam` files are deleted (they are being redone with this approach).
3. **Indexes**: Cxt genome index and RVFV index in `Genome/` (see below).
4. **Two-stage alignment** (`02_RNA_Cxthisatarray_unaligned.sh`):
   - Stage 1: trimmed reads to the Cxt genome. Writes the host alignments, and the pairs that
     did not align (`--un-conc-gz`).
   - Stage 2: those unaligned pairs to RVFV. Writes the RVFV alignments, and the pairs that
     aligned to neither.
5. **ssRNA+ re-alignment**: re-align the "neither" reads to the ssRNA+ strand index
   using the earlier ssRNA+ scripts.
6. **Report** results in the Excel summary as before.
7. **Decide** whether to repeat with a specific virus family (e.g. Flaviviridae).

## Inputs required before running the HISAT2 array

The script checks 1 to 4 and exits if any are missing.

1. Trimmed paired-end FASTQs in `Data/slurm_script/trimmed/`
   (`trimmed.<sample>_R1_001.fastq.gz` and matching `_R2_001`).
2. Sample list `CulexTrimmed.txt` in that folder (see "Setting up the array").
3. Cxt genome index: `Genome/Index_Genome/hisatIndex/CtarK1_index.*.ht2` (8 files).
4. RVFV index: `Genome/Index_RVFV/RVFV_index.*.ht2` (8 files).
5. A working conda environment with `hisat2` (and `samtools` for downstream checks).
6. `slurmlogs/outlog` and `slurmlogs/errlog` folders. SLURM cannot write logs if they are missing.

## Building the RVFV index (one time)

RVFV has three segments, so all three go into one FASTA. Reference: RefSeq
`NC_014397` (L), `NC_014396` (M), `NC_014395` (S).

```bash
mkdir -p Genome/Index_RVFV && cd Genome/Index_RVFV

# download all three segments as one FASTA
curl -o RVFV_all.fasta \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_014397,NC_014396,NC_014395&rettype=fasta&retmode=text"

# check: 3 sequences, ~6.4 kb (L), ~3.9 kb (M), ~1.7 kb (S)
grep ">" RVFV_all.fasta
awk '/^>/{if(n)print n,l; n=$1; l=0; next}{l+=length($0)}END{print n,l}' RVFV_all.fasta

# build (hisat2 environment active); takes seconds
hisat2-build -p 4 RVFV_all.fasta RVFV_index
ls RVFV_index.*.ht2      # expect 8 files
```

Do not run `cat ... > RVFV_all.fasta` after the download: the `>` empties the file first and
wipes the good copy.

## Setting up the array with CulexTrimmed.txt

`CulexTrimmed.txt` lists one R1 filename per line, one line per sample. The array task ID is
the 0-indexed line number: task 0 runs line 1, task 1 runs line 2, and so on.

```bash
cd Data/slurm_script/trimmed
ls | grep -v / | grep "_R1_001.fastq.gz" | grep "trimmed" > CulexTrimmed.txt
wc -l CulexTrimmed.txt        # N lines, so the array range is 0 to N-1
```

Regenerate the list whenever trimmed files are added or removed. The list is alphabetical,
so adding samples can shift which task ID maps to which sample.

## Running

Test one sample first, check its outputs, then submit the rest:

```bash
sbatch --array=0 02_RNA_Cxthisatarray_unaligned.sh            # test sample
sbatch --array=1-23 02_RNA_Cxthisatarray_unaligned.sh         # rest (24 samples)
```

Per-sample run time has been roughly 10 minutes to 1.5 hours at 16 CPUs for the earlier
alignments. Use a time limit of about 6 hours at 8 CPUs, since a job killed at its limit loses
all its work.

## HISAT2 settings

| Option | Why |
|---|---|
| `--phred33` | Standard Illumina quality encoding |
| `--rna-strandness RF` | Stranded library (read 1 reverse, read 2 forward) |
| `--un-conc-gz <prefix>_R%.fastq.gz` | Writes pairs that did not align concordantly; `%` becomes 1 and 2 |
| `--summary-file` | Per-sample alignment statistics for the Excel report |
| `--no-spliced-alignment` (stage 2) | Viral RNA is not spliced |

## Outputs

Under `Data/<script name>/`. When submitted with `sbatch` the script name resolves to
`slurm_script`, so outputs land in `Data/slurm_script/`.

| Folder | Contents |
|---|---|
| `1_CxtAligned/` | Host alignments (SAM) |
| `1_CxtUnaligned/` | Pairs not aligned to Cxt; input to stage 2. Keep until counts are checked. |
| `2_RVFVAligned/` | RVFV alignments (SAM) |
| `3_NeitherAligned/` | Pairs aligned to neither: the dataset for ssRNA+ re-alignment / BLAST |
| `Summaries/` | `<sample>cxt_summary.txt` and `<sample>rvfv_summary.txt` |

## Verifying a run

- The RVFV summary's total pairs should equal the number of pairs unaligned to Cxt.
- Count read pairs in a FASTQ as lines / 4. R1 and R2 must have the same count.
- An uninfected sample (e.g. BF) should have no reads aligned to RVFV, so
  `3_NeitherAligned/` should hold the same reads as `1_CxtUnaligned/`.
- HISAT2 with multiple threads can write reads in a different order, so compare sorted
  read names, not the files byte by byte.
- By default HISAT2 also writes unaligned reads into the SAM, so a SAM can be much larger than
  its aligned reads. Use `--no-unal` to drop them.

## Known caveats

- **Concordant vs. unaligned.** `--un-conc-gz` writes pairs that did not align as a proper
  pair. Pairs with only one host-aligned mate may pass into the "neither" set. If that matters,
  extract pairs where both mates are unmapped from a BAM instead (`samtools view -f 12`).
- **Disk space.** Cxt SAMs run to hundreds of GB per sample. Convert to BAM, or add `--no-unal`,
  before running all samples. A gzip "Disk quota exceeded" error means scratch is full and
  the partial outputs are corrupt; delete them and rerun.
- **Scratch is purged every 90 days.** Move conda environments and anything to keep to
  `/projects`. A conda install whose Python has been purged fails with
  "bad interpreter: No such file or directory". Environment yml exports and rebuild scripts live
  in `/projects/<user>/CondaEnvs`.
- **Home directory quota is small (2 GB).** Write SLURM logs and working files to scratch.

## Software

fastp (trimming), HISAT2 2.2.2 (alignment), samtools (BAM conversion and checks),
HTSeq (counting), SLURM on CURC Alpine.

