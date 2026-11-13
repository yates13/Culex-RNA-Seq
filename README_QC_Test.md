# Verifying a two-stage HISAT2 pipeline (host, then virus)

`verify_all.sh` checks that the outputs of a two-stage HISAT2 run are complete and consistent with each other, before they are used downstream. It was written to verify files received from a colleague, where the pipeline itself was not run by the person checking.

## The pipeline being checked

```
trimmed reads
   |  stage 1: align to host genome (HISAT2, --un-conc-gz keeps unaligned pairs)
   v
CXT_Aligned/*.sam           all stage-1 records (aligned and unaligned)
CXT_Unaligned/*.fastq.gz    pairs that did not align concordantly to the host
   |  stage 2: align those pairs to a virus index (HISAT2, --un-conc-gz)
   v
MP12_Aligned/*.sam          all stage-2 records
Neither_Aligned/*.fastq.gz  pairs that aligned concordantly to neither
```

The `Neither` FASTQs are the end product (for example, input for ssRNA+ virus alignment or BLAST). The goal is to confirm they contain exactly the reads they should.

## Expected layout

```
test_files/
  CXT_Aligned/trimmed.<sample>.sam
  CXT_Unaligned/trimmed.<sample>unaligned_R*.fastq.gz
  MP12_Aligned/trimmed.<sample>.sam
  Neither_Aligned/trimmed.<sample>unaligned_R*.fastq.gz
```

`<sample>` is a prefix such as `BF-RNA-D1-R1-Cxt_` (it includes the trailing `Cxt_`). Each sample is given a role:

- `control`: must not align to the virus (mock-infected or uninfected)
- `virus`: must align to the virus (infected)

**File names can mislead.** In the data this was written for, some stage-2 SAMs carried a Cxt-style name. Trust the header (`@SQ` count, `@PG` command line) over the file name; check C1 does this for you.

## Running

```bash
cd <results directory>
T=/path/to/test_files W=/path/to/work sbatch verify_all.sh
# other samples:
SAMPLE_LIST="X-Cxt_:control Y-Cxt_:virus" T=... W=... sbatch verify_all.sh
```

- Output goes to `verify_all.<jobid>.out`. Each check prints `[PASS]`, `[FAIL]` or `[SKIP]` with its number; the last line is a summary, and the exit code is non-zero if any check failed.
- The `#SBATCH` lines (partition, qos, time, memory) are for one particular cluster; edit them for yours. Without SLURM, `bash verify_all.sh` works the same way.
- Each SAM is read three times, and the read-name lists in `$W` can reach tens of GB for large samples. Delete `$W` when finished.

## Checks

| ID | Tests | Expected | A failure suggests |
|---|---|---|---|
| C1 | Stage-1 SAM has more `@SQ` references than the stage-2 SAM | pass | files swapped or mislabelled |
| C2 | Stage-1 SAM ends with a newline and a full record | pass | file cut off during write or transfer |
| C3 | Stage-2 SAM ends with a newline and a full record | pass | file cut off during write or transfer |
| C4 | Stage-1 SAM has properly paired reads | > 0 | wrong file, or empty alignment |
| C5 | Stage-1 unaligned FASTQ has reads | > 0 | wrong file |
| C6 | No read is both properly paired in stage 1 and in the unaligned FASTQ | 0 | unaligned FASTQ contains host-aligned pairs |
| C7 | Unaligned FASTQ equals (all stage-1 SAM reads minus properly paired) | 0 differences | SAM and FASTQ from different runs, or SAM incomplete |
| C8 | Every `Neither` read is in the stage-1 unaligned FASTQ | 0 outside | stage 2 read the wrong input |
| C9 | No `Neither` read is properly paired in stage 1 | 0 | host-aligned pairs leaked through |
| C10 | Stage-2 SAM covers exactly the stage-1 unaligned reads | 0 differences | stage-2 SAM truncated or from another input |
| C11 | `Neither` equals (stage-1 unaligned minus stage-2 properly paired) | 0 differences | `Neither` not what stage 2 should output |
| C12 | Every virus-mapped read is in the stage-1 unaligned FASTQ | 0 outside | stage 2 ran on the wrong reads |
| C13 | No virus-mapped read remains in `Neither` | 0 | virus reads were not removed |
| C14 | Role-specific: `control` has `Neither` identical to unaligned; `virus` has `Neither` a strict subset | pass | no reads removed for a `virus` sample, or reads removed from a `control` |
| C15 | Role-specific: `control` has 0 virus-mapped reads; `virus` has more than 0 | pass | contamination in a control, or no virus signal |

C10 to C13 are skipped (and reported as `[SKIP]`) if the stage-2 SAM is missing.

The script also prints an `[info]` line: the number of unaligned pairs with a host-mapped mate (one mate mapped, or both mapped non-concordantly). With `--un-conc-gz`, such pairs stay in the unaligned and `Neither` FASTQs. If only pairs with both mates unmapped are wanted, filter them out downstream.

## Limitations

- Checks compare read names, not sequences. They do not confirm the reads themselves are intact.
- The two mates are not cross-checked. If the stage-1 unaligned FASTQ and the `Neither` FASTQ are different mates (R1 and R2), the names still match, but nothing verifies the other mate's file.
- There is no raw input to compare against, so the script cannot prove stage 1 saw every trimmed read. Compare the total read-name count of the stage-1 SAM against the trimmed FASTQ or the HISAT2 summary file.
- Name-level checks can miss a truncation that damages only the final pair. C2 and C3 cover that case.
- Both SAMs must contain unaligned records (no `--no-unal`), because C7 and C10 depend on them.
- `proper` means SAM flag `0x2` is set on any record of the pair; `mapped` means flag `0x4` is unset on any record.
