library(ggplot2)

bsub_genes <- readDNAStringSet("~/Desktop/projectRBS/All-genes-of-B.-subtilis-subtilis-168.fasta")
bsub_genome_f <- readDNAStringSet("~/Desktop/projectRBS/GCF_000009045.1_ASM904v1_genomic.fna")
bsub_genome_r <- complement(bsub_genome_f)
bsub_annotation <- read.table("~/Desktop/projectRBS/bsub168_genomic.gff",skip = 7, sep = "\t")
bsub_abundance <- read.table("~/Desktop/projectRBS/224308-WHOLE_ORGANISM-integrated.txt")
bsub_info <- strsplit(bsub_annotation$V9,";")

ggplot(data = bsub_abundance, mapping = aes(x = log(V3+0.0001))) +
  geom_density()+
  theme_bw()+
  labs(x = "Abundance")

table <- data.frame(matrix(nrow = length(bsub_genes),ncol = 3))
names(table) <- c("Locus","Locus.tag","ORF")

for (i in 1:length(bsub_genes)) {
  info <- strsplit(names(bsub_genes),"\\|")[[i]]
  table$Locus[i] <- info[1]
  table$Locus.tag[i] <- info[2]
  table$ORF[i] <- as.character(bsub_genes[[i]])
}
