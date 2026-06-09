# E. coli RBS Prediction and Optimization

This repository contains analysis scripts for identifying and evaluating
ribosome binding site (RBS) features in *Escherichia coli* K-12 MG1655.
The workflow combines gene-start prediction, Shine-Dalgarno motif matching,
upstream sequence extraction, RNA secondary-structure analysis, and exploratory
models relating sequence/structure features to expression data.

The project was originally developed as an exploratory research workspace. It
has been reorganized so GitHub tracks source code and documentation, while raw
genomes, generated FASTA files, intermediate CSV tables, figures, and local
RStudio state remain outside version control.

## Repository Layout

```text
.
|-- scripts/          R analysis and workflow scripts
|-- data/
|   |-- raw/          raw external inputs, not tracked by git
|   `-- processed/    generated/intermediate tables, not tracked by git
|-- results/
|   |-- figures/      generated plots and PDFs, not tracked by git
|   `-- tables/       exported result tables, not tracked by git
|-- docs/             project notes and reproducibility guidance
|-- DESCRIPTION       R dependency metadata
`-- README.md
```

## Main Workflow

The scripts are historical analysis steps rather than one consolidated command.
Run them from the repository root so relative file paths resolve consistently.
See `docs/SCRIPT_INDEX.md` for a fuller script-by-script index.

1. Predict candidate genes and start sites with Prodigal:
   `scripts/bulkProdigal.R`, `scripts/individualProdigal.R`
2. Extract upstream regions and RBS candidates:
   `scripts/locateUpstreams.R`, `scripts/extractRBS.R`
3. Match Shine-Dalgarno motifs:
   `scripts/vMatch.R`, `scripts/vMatch1.R`, `scripts/vMatch2.R`,
   `scripts/vMatch3.R`, `scripts/vMatch4.R`
4. Analyze RNA secondary structure and RBS pairing:
   `scripts/bashOfRNAfold.sh`, `scripts/calculateProportionPaired.R`,
   `scripts/calculateProportionPairedMXFOLD2.R`,
   `scripts/calculateRBSProportionPairedMXFOLD2.R`,
   `scripts/furtherAnnotationOnSS.R`
5. Explore expression, tRNA, folding-energy, and modeling relationships:
   `scripts/prelimAnalysis.R`, `scripts/RFPrep_and_SSAnalysis.R`,
   `scripts/modellingEnergy.R`, `scripts/tRNA_analysis.R`

## Dependencies

R package dependencies are listed in `DESCRIPTION`. Several scripts also expect
external command-line biology tools and local datasets, including:

- Prodigal for gene prediction
- ViennaRNA/RNAfold for folding predictions
- MXfold2 outputs for some secondary-structure analyses
- E. coli K-12 MG1655 genome and annotation files
- expression and abundance tables used in the historical analysis

See `docs/REPRODUCIBILITY.md` for the expected local inputs and setup notes.

## Version-Control Policy

Git should track:

- R source scripts
- documentation
- small project metadata/configuration files

Git should not track:

- raw genome/transcriptome files
- generated FASTA/CSV/GFF/DBN/fold files
- plots, PDFs, and model output artifacts
- `.RData`, `.Rhistory`, `.Rproj.user`, `.DS_Store`, and other local state

If a generated artifact is needed for publication or a stable release, place it
under `results/` and explicitly unignore that file in `.gitignore`.
