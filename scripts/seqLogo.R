### Faire un logo de séquence consensuelle
library(Biostrings)
library(seqLogo)
library(ggseqlogo)
library(ggplot2)
library(seqinr)

fasta <- readDNAStringSet("~/Desktop/projectRBS/eschColi_K_12_MG1655-mature-tRNAs.fa")
RBS <- DNAStringSet()
for (i in 1:length(fasta)) {
  seq_length <- length(fasta[[i]])
  RBS[[i]] <- substring(fasta[[i]],first = seq_length-49 , last = seq_length)
}
RBS_seq <- as.character(RBS)

### Automatic formation of PWM and Logo
p1 <- ggseqlogo(RBS_seq,method = "bits")
p2 <- ggseqlogo(RBS_seq,method = "prob")
gridExtra::grid.arrange(p1, p2)
