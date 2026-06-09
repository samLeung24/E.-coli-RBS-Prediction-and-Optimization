# Reproducibility Notes

Run scripts from the repository root unless a script says otherwise. Several
historical scripts use relative paths to files that were stored next to this
repository in the original lab workspace.

## R Packages

Install the packages listed in `DESCRIPTION`, including Bioconductor packages:

```r
install.packages(c(
  "cluster",
  "dplyr",
  "ggplot2",
  "ggseqlogo",
  "gridExtra",
  "seqinr",
  "stringr",
  "tidyverse"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("Biostrings", "seqLogo"))
```

## External Tools

Some workflow stages require command-line tools outside R:

- Prodigal
- ViennaRNA/RNAfold
- MXfold2

Install these separately and make sure their executables are available on
`PATH` before running scripts that call them.

## Expected Local Inputs

The historical workflow expects local copies of files such as:

- E. coli K-12 MG1655 genome FASTA files
- gene annotation tables
- Prodigal prediction outputs
- expression and abundance tables
- RNAfold/MXfold2 secondary-structure outputs

Place raw inputs under `data/raw/` and generated intermediates under
`data/processed/` when modernizing the workflow. The existing scripts still
contain older relative paths and may need path cleanup before they are fully
portable.
