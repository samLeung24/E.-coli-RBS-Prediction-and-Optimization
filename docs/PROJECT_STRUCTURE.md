# Project Structure

This repository is organized as a source-first GitHub project.

## Tracked Directories

- `scripts/`: R scripts used across the historical analysis workflow.
- `docs/`: project notes, setup details, and reproducibility guidance.
- `data/raw/README.md`: documents expected raw inputs without committing them.
- `data/processed/README.md`: documents generated/intermediate tables.
- `results/figures/README.md`: documents generated visual outputs.
- `results/tables/README.md`: documents exported result tables.

## Local-Only Artifacts

The original workspace contains many large or generated artifacts: genome FASTA
files, Prodigal output, RNAfold/MXfold2 output, CSV intermediates, plots, PDFs,
and RStudio local state. These are intentionally ignored by `.gitignore`.

Keep these files locally or archive them in a data repository/release asset if
they need to be preserved. Avoid committing them directly to the main source
repository unless they are small, stable, and required for examples or tests.
