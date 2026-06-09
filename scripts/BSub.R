library(ggplot2)
library(dplyr)

bsub_genes <- readDNAStringSet("~/Desktop/projectRBS/All-genes-of-B.-subtilis-subtilis-168.fasta")
bsub_genome_f <- readDNAStringSet("~/Desktop/projectRBS/GCF_000009045.1_ASM904v1_genomic.fna")
bsub_genome_r <- reverse(complement(bsub_genome_f))
bsub_annotation <- read.table("~/Desktop/projectRBS/bsub168_genomic.gff",skip = 7, sep = "\t")
bsub_abundance <- read.table("~/Desktop/projectRBS/224308-WHOLE_ORGANISM-integrated.txt")
  for (i in 1:length(bsub_abundance$V1)) {
    bsub_abundance$Locus[i] <- strsplit(bsub_abundance$V2[i],"\\.")[[1]][2]
  }
bsub_info <- strsplit(bsub_annotation$V9,";")

ggplot(data = bsub_abundance, mapping = aes(x = log(V3+0.0001))) +
  geom_density()+
  theme_bw()+
  labs(x = "Abundance")

bsub <- data.frame(matrix(nrow = length(bsub_genes),ncol = 3))
names(bsub) <- c("Locus","Locus.tag","ORF")

# Extract gene information
for (i in 1:length(bsub_genes)) {
  print(i)
  info <- strsplit(names(bsub_genes),"\\|")[[i]]
  bsub$Locus[i] <- info[1]
  bsub$Locus.tag[i] <- info[2]
  bsub$ORF[i] <- as.character(bsub_genes[[i]])
}

## Search for genes in the genome and extract 30 nt upstream
for (i in 1:length(bsub_genes)) {
  print(i)
  ORF <- substr(bsub$ORF[i],1,50)
  if (grepl(pattern = ORF,bsub_genome_f)) {
    index <- regexpr(ORF,bsub_genome_f)[1]
    bsub$Upstream30[i] <- substr(bsub_genome_f,index-30,index-1)
  } else if (grepl(pattern = ORF,bsub_genome_r)) {
    index <- regexpr(ORF,bsub_genome_r)[1]
    bsub$Upstream30[i] <- substr(bsub_genome_r,index-30,index-1)
  } else bsub$Upstream30[i] <- NA
}

## Calculate GC content
for (i in 1:length(bsub$ORF)) {
  gc_count <- sum(strsplit(bsub$ORF[i], "")[[1]] %in% c("G", "C"))
  bsub$gc_content[i] <- (gc_count / nchar(bsub$ORF[i]))
}

## Combine first 30 nt of the ORF and 30 nt upstream

for (i in 1:length(bsub$ORF)) {
  locus <- bsub$Locus[i]
  header <- paste0(">",locus)
  seq <- paste0(substr(bsub$Upstream30[i],1,30),substr(bsub$ORF[i],1,30))
  seq <- as.character(RNAString(DNAString(seq)))
  cat(header, file = "BSUB_3030.fasta", sep = "\n", append = TRUE);
  cat(seq, file = "BSUB_3030.fasta", sep = "\n", append = TRUE)
}

### Shell Script for Mxfold2 prediction ###

### Read Mxfold2 prediction
mxf <- readLines("~/Desktop/projectRBS/BSUB_3030.predict")
for (i in 1:length(bsub$Locus)){
  info = strsplit(mxf[3*i], " ")
  bsub$SecondaryStructure[i] = info [[1]][1]
  bsub$FoldingScore[i] = gsub("\\(|\\)","",info[[1]][2])
}
