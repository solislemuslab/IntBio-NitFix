
# fastp Trimming QC Check

This folder contains the quality-control summary for the `fastp` trimming step of the August 2025 `symbiosis_sorted` reads.

The goal of this QC step was to check whether Ryan's recommended `fastp` trimming settings improved or maintained read quality before mapping the reads to the symbiosis island reference.

## Folder Contents

| File | What it contains | How to use it |
|---|---|---|
| `README.md` | This explanation file. | Use this to understand the QC table, figures, and log file. |
| `fastp_logs_ryan_full_summary.tsv` | Per-sample QC summary extracted from 1,116 `fastp` JSON reports. | Main table for reporting read counts, read retention, Q20/Q30 rates, and GC content before/after trimming. |
| `fastp_qc_summary_run.log` | Terminal output from the QC-summary script. | Reproducibility record showing that 1,116 JSON files were processed and recording the main average QC values. |
| `fastp_quality_before_trimming_all_samples_lightpurple_mean_blue.svg` | Per-base quality profiles before trimming for all samples. | Use to show raw read quality before trimming. |
| `fastp_quality_after_trimming_all_samples_lightpurple_mean_blue.svg` | Per-base quality profiles after trimming for all samples. | Recommended report figure for showing quality after trimming. |
| `fastp_quality_mean_only.svg` | Mean per-base quality profiles before and after trimming, without individual sample lines. | Optional simpler figure if the all-sample figure is too busy. |
| `fastp_quality_before_after_and_retention.svg` | Earlier summary figure showing before/after Q20/Q30 and read retention. | Optional internal summary; not the preferred final figure. |
| `fastp_quality_all_samples_gray_mean_blue.svg` | Earlier version of the after-trimming all-sample figure with gray sample lines. | Older version; replaced by the clearer light-purple figure. |

## Figures 

1. `fastp_quality_before_trimming_all_samples_lightpurple_mean_blue.svg`
2. `fastp_quality_after_trimming_all_samples_lightpurple_mean_blue.svg`

In both figures, each light-purple line represents one sample and the blue line represents the mean across all 1,116 samples. Read 1 corresponds to the forward reads (`R1`) and Read 2 corresponds to the reverse reads (`R2`). The red dashed horizontal lines mark the Q20 and Q30 quality thresholds.

These figures are useful because they show both the average read quality and the variation among samples. They also make it clear that most bases remain above the Q30 threshold across read positions.

## Summary Table Columns

Each row in `fastp_logs_ryan_full_summary.tsv` represents one sample.

| Column | Description |
|---|---|
| `sample` | Sample name used throughout the pipeline. |
| `before_total_reads` | Total number of reads before trimming. This includes both paired-end reads. |
| `after_total_reads` | Total number of reads remaining after `fastp` trimming and filtering. |
| `reads_retained_percent` | Percentage of reads retained after trimming, calculated as `(after_total_reads / before_total_reads) x 100`. |
| `before_q20_rate` | Fraction of bases with Phred quality score >= Q20 before trimming. Q20 is approximately 99% base-call accuracy. |
| `after_q20_rate` | Fraction of bases with Phred quality score >= Q20 after trimming. |
| `before_q30_rate` | Fraction of bases with Phred quality score >= Q30 before trimming. Q30 is approximately 99.9% base-call accuracy. |
| `after_q30_rate` | Fraction of bases with Phred quality score >= Q30 after trimming. |
| `before_gc_content` | Fraction of bases that are G or C before trimming. For example, `0.58` means 58% GC content. |
| `after_gc_content` | Fraction of bases that are G or C after trimming. This checks whether trimming strongly changed base composition. |

## What Q20 and Q30 Mean

FASTQ files store sequencing quality as characters, not normal numbers. For example, a quality line may look like:

```text
IIIIIIIIIIIIIIIIIIII
```

Each character represents the confidence of one base call. These characters are converted into Phred quality scores.

| Score | Approximate error rate | Approximate accuracy |
|---|---:|---:|
| Q20 | 1 error per 100 bases | 99% |
| Q30 | 1 error per 1,000 bases | 99.9% |

Therefore, Q20 and Q30 rates report the fraction of bases that meet or exceed those quality thresholds.

## What GC Content Means

GC content is the fraction of bases that are either `G` or `C`:

```text
GC content = (G + C) / total bases
```

GC content is not expected to improve after trimming. Instead, it is used as a bias check. If GC content changes strongly after trimming, it may suggest that trimming removed some sequence types more than others.

## QC Run Log

The file `fastp_qc_summary_run.log` records the terminal output from the QC-summary script. It confirms that 1,116 `fastp` JSON files were found and processed, and that the summary table and QC figure were written successfully. The log also records the main average QC metrics used in the report, including mean read retention, mean Q20 before/after trimming, and mean Q30 before/after trimming.

## Main Result

The trimming QC summarized 1,116 `fastp` JSON reports.

| Metric | Result |
|---|---:|
| Samples summarized | 1,116 |
| Mean reads retained after trimming | 95.84% |
| Mean Q20 before trimming | 98.73% |
| Mean Q20 after trimming | 99.65% |
| Mean Q30 before trimming | 95.39% |
| Mean Q30 after trimming | 97.49% |

The raw reads were already high quality before trimming. After `fastp` trimming, Q20 and Q30 rates increased modestly while most reads were retained. GC content remained similar before and after trimming, indicating that trimming improved read quality without strongly changing base composition.


The raw reads were already high quality before trimming. After `fastp` trimming, Q20 and Q30 rates increased modestly while most reads were retained. GC content remained similar before and after trimming, indicating that trimming improved read quality without strongly changing base composition.

