### Construct a linking sequence
library(ggplot2)
library(gridExtra)

table <- read.csv("forRF8.csv")
subseq <- readRNAStringSet("30nt_upstream_30nt_ORF.RNA.fasta")

Genome <- readDNAStringSet("../MG1655_complete_genome.fasta")
rRNA_start <- 4166659
rRNA_end <- 4168200
rRNA <- Genome[[1]][rRNA_start:rRNA_end]
rRNA <- RNAStringSet(rRNA)

mRNA_name <- table$Locus
mRNA_sequence <- table$ORF
sc <- substr(table$ORF,1,3)
### Extract annotated secondary structure of mRNAs

dir <- "~/Desktop/projectRBS/30_30_structural_annotation/"
for (i in 1:length(table$Locus)) {
  fileName <- paste0(dir,table$gene_id[i],"|",table$Locus[i],".st")
  info <- readLines(fileName)
  info_anno <- info[6]
  table$SS_info[i] <- info_anno
}

write.csv(table,"forRF8.csv",row.names = FALSE)

for (i in 1:length(subseq)) {
  dir <- "~/Desktop/projectRBS/30_30_rRNA_fastas/"
  fileName <- names(subseq[i])
  cat(paste0(">",fileName), file = paste0(dir,fileName,".fasta"),sep = "\n")
  cat(as.character(subseq(subseq[[i]],1,30)), file = paste0(dir,fileName,".fasta"), append = TRUE,sep = "\n")
  cat(as.character(subseq(rRNA,1534,1542)),file = paste0(dir,fileName,".fasta"), append = TRUE)
}

hybrid <- readLines("../30_30_rRNA_RNAfold.dbn")
nGenes <- length(hybrid)/2
hybrid_info <- data.frame(Locus =rep("",nGenes,),
                          structure = rep("",nGenes),
                          start = rep(0,nGenes),
                          end = rep(0,nGenes),
                          energy = rep(0,nGenes))
for (i in 1:nGenes) {
  print(i)
  name <- hybrid[2*i-1]
  input_string <- hybrid[[2*i]]
  pattern1 <- "^([\\.\\(\\)&]+)"
  pattern2 <- "\\s*(\\d+,\\d+)\\s*(?=\\s*:)"
  pattern3 <- "\\(([^)]+)\\)$"

  # Extracting
  hybrid_info$Locus[i] <- strsplit(name,"\\|")[[1]][2]
  hybrid_info$structure[i] <- str_match(input_string, pattern1)[, 1]
  sep <- str_match(input_string, pattern2)[, 1]
  sep <- trimws(sep)
  hybrid_info$start[i] <- strsplit(sep,",")[[1]][1]
  hybrid_info$end[i] <- strsplit(sep,",")[[1]][2]
  hybrid_info$energy[i] <- gsub("[())]","",str_match(input_string, pattern3)[, 1])
}

p1 <- ggplot(data = hybrid_info,aes(x = (30-as.numeric(end)))) + geom_density()
p2 <- ggplot(data = hybrid_info,aes(x = (30-as.numeric(end)))) + geom_density()
p3 <- ggplot(data = hybrid_info,aes(x = (30-as.numeric(end)))) + geom_density()
p4 <- ggplot(data = hybrid_info,aes(x = (30-as.numeric(end)))) + geom_density()
p5 <- ggplot(data = hybrid_info,aes(x = (30-as.numeric(end)))) + geom_density()
p6 <- ggplot(data = hybrid_info,aes(x = (30-as.numeric(end)))) + geom_density()
p7 <- ggplot(data = hybrid_info,aes(x = (30-as.numeric(end)))) + geom_density()
p8 <- ggplot(data = hybrid_info,aes(x = (30-as.numeric(end)))) + geom_density()
p9 <- ggplot(data = hybrid_info,aes(x = (30-as.numeric(end)))) + geom_density()

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)
