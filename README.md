# exploreatacseq
Exploratory analysis of multiple ATAC-seq datasets

## Installation

```r
library(devtools)
install_github("anthony-aylward/exploreatacseq")
```

## Example usage

Here is an example JSON file that is formatted for use with `exploreatacseq`
(see [tssenrich](https://github.com/anthony-aylward/tssenrich) for computation
of TSS enrichment values):

```json
{
  "SAMN10079665": {
    "untreated": {
      "reads": "SAMN10079665_untreated.sort.filt.rmdup.bam",
      "tssenrich": 7.62,
      "peaks": "SAMN10079665_untreated_peaks.narrowPeak"
    },
    "dex": {
      "reads": "SAMN10079665_dex.sort.filt.rmdup.bam",
      "tssenrich": 6.67,
      "peaks": "SAMN10079665_dex_peaks.narrowPeak"
    }
  },
  "SAMN09767462": {
    "untreated": {
      "reads": "SAMN09767462_untreated.sort.filt.rmdup.bam",
      "tssenrich": 5.37,
      "peaks": "SAMN09767462_untreated_peaks.narrowPeak"
    },
    "dex": {
      "reads": "SAMN09767462_dex.sort.filt.rmdup.bam",
      "tssenrich": 5.16,
      "peaks": "SAMN09767462_dex_peaks.narrowPeak"
    }
  },
  "AFA3256": {
    "untreated": {
      "reads": "AFA3256_untreated.sort.filt.rmdup.bam",
      "tssenrich": 5.48,
      "peaks": "AFA3256_untreated_peaks.narrowPeak"
    },
    "dex": {
      "reads": "AFA3256_dex.sort.filt.rmdup.bam",
      "tssenrich": 5.88,
      "peaks": "AFA3256_dex_peaks.narrowPeak"
    }
  }
}
```

Assuming the above file is at `atacseq.json`, here is a simple example
application of `exploreatacseq` in R:
```r
library(exploreatacseq)

explore(
  "atacseq.json",
  "output_prefix",
  treatment_groups = list(c("untreated", "dex")),
  cores = 2
)
```