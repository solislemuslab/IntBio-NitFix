# fastp Trimming QC Summary

This folder contains quality-control outputs summarizing the `fastp` trimming step for the August 2025 `symbiosis_sorted` reads.

The trimming step followed Ryan's recommended `fastp` settings and was applied to all paired-end samples before mapping to the symbiosis island reference.

## Files

| File | Description |
|---|---|
| `fastp_logs_ryan_full_summary.tsv` | Per-sample summary table extracted from the `fastp` JSON reports. |
| `fastp_quality_before_after_and_retention.svg` | QC figure showing read quality before/after trimming and read-retention summary. |
| `fastp_qc_summary_run.log` | Log file from the QC-summary script. |

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

