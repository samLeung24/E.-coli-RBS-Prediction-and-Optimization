# Script Index

Run scripts from the repository root unless noted otherwise.

## Gene Prediction and Sequence Extraction

- `scripts/bulkProdigal.R`: bulk Prodigal workflow over genome inputs.
- `scripts/individualProdigal.R`: per-gene Prodigal workflow.
- `scripts/gff_to_csv.sh`: extracts CDS records from individual GFF files.
- `scripts/locateUpstreams.R`: historical upstream/RNAfold workflow script.
- `scripts/extractRBS.R`: extracts RBS-region sequence data.

## Shine-Dalgarno Matching

- `scripts/vMatch.R`: motif matching workflow.
- `scripts/vMatch1.R`: motif matching variant.
- `scripts/vMatch2.R`: motif matching variant.
- `scripts/vMatch3.R`: motif matching variant.
- `scripts/vMatch4.R`: motif matching variant.

## Secondary Structure

- `scripts/bashOfRNAfold.sh`: generates RNAfold command batches.
- `scripts/calculateProportionPaired.R`: calculates paired-region proportions.
- `scripts/calculateProportionPairedMXFOLD2.R`: MXfold2 paired-region analysis.
- `scripts/calculateRBSProportionPairedMXFOLD2.R`: RBS-specific MXfold2 analysis.
- `scripts/calculateRBS-rRNAcomplementarity.R`: RBS/rRNA complementarity scoring.
- `scripts/furtherAnnotationOnSS.R`: additional secondary-structure annotation.

## Exploratory Analysis and Modeling

- `scripts/prelimAnalysis.R`: preliminary correlations and regressions.
- `scripts/RFPrep_and_SSAnalysis.R`: feature preparation and secondary-structure analysis.
- `scripts/modellingEnergy.R`: folding-energy modeling.
- `scripts/Clustering.R`: clustering exploration.
- `scripts/RNASeq_analysis.R`: RNA-seq expression analysis.
- `scripts/seqLogo.R`: sequence-logo generation.
- `scripts/tRNA_analysis.R`: tRNA/codon availability analysis.
- `scripts/Termination.R`: termination-region exploration.
- `scripts/BSub.R`: *B. subtilis* comparison workflow.
- `scripts/bsub_expression_exploration.R`: *B. subtilis* expression exploration.
