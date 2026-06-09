library(ggplot2)
library(Biostrings)
library(stringr)
library(dplyr)

### Retrieve expression level data
exptab <- read.table("../511145-WHOLE_ORGANISM-integrated.txt",sep = "\t")
ratiotab <- read.csv("../edited_processed_sequences_pax_dna_prot_Ecoli-K12_20230118.csv")
seqtab <- read.csv("../gene_abundances.csv")
names(seqtab)[1] <- "Locus"
for (i in 1:length(exptab$V2)) {
  exptab$Locus.tag[i] <- strsplit(exptab$V2[i],".", fixed =TRUE)[[1]][2]
}
temptab <- cbind(gene_id=ratiotab$GeneID,Locus.tag=ratiotab$Locus.tag)
finaltab <- merge(temptab,exptab[,c("V3","Locus.tag")],by = "Locus.tag")
## <Abandoned> positives <- read.csv("bulkProdigal.csv")
## <Abandoned> positives <- positives[complete.cases(positives),]
## <Abandoned> finaltab <- finaltab[finaltab$gene_id %in% positives$gene_id,]
colnames(finaltab) <- c("Locus.tag", "gene_id", "expression_level")

### Import secondary structure predictions (Mxfold2)
sectab <- read.csv("WholeProportionPairedMxfold2.csv")
table <- merge(sectab,finaltab, by = "gene_id")
## Abandoned table <- merge(table,positives,by = "gene_id")
table$loged_el <- log(table$expression_level+0.0001,base = 10)
table <- table[is.finite(table$loged_el),]

## <Abandoned> plot(table$spacer_length,table$loged_el,main = "Expression Level v. Spacer Length (30)")
## <Abandoned> abline(lm(table$loged_el~table$spacer_length),col = "red")
## <Abandoned> summary(lm(table$spacer_length~table$loged_el))


par(mfrow=c(2,2))

plot(table$X30nt_score,table$loged_el,main = "Expression Level v. Folding Score (30)")
abline(lm(table$loged_el~table$X30nt_score),col = "red")
summary(lm(table$X30nt_score~table$loged_el))

plot(table$X40nt_score,table$loged_el,main = "Expression Level v. Folding Score (40)")
abline(lm(table$loged_el~table$X40nt_score),col = "red")
summary(lm(table$X40nt_score~table$loged_el))

plot(table$X50nt_score,table$loged_el,main = "Expression Level v. Folding Score (50)")
abline(lm(table$loged_el~table$X50nt_score),col = "red")
summary(lm(table$X50nt_score~table$loged_el))

par(mfrow=c(2,2))
plot(lm(table$loged_el~table$X30nt_score))

table <- table[,-2]
moreInfo <- ratiotab[,c("Locus.tag","Locus","dna_seqs","X5utr_seqs","cai")]
table <- merge(table,moreInfo,by = "Locus.tag")
names(ratiotab)

names(table)[15] <- "ORF"
names(table)[16] <- "RBS"
table$gc_content <- str_count(tolower(table$ORF), "g|c")/nchar(table$ORF)

table <- merge(table,seqtab,by = "Locus")

table$SC_pairing <- c()
table$start_codon <- c()
for (i in 1:length(table$gene_id)) {
  SC <- substr(table$X30nt[i],1,3)
  table$start_codon[i] <- SC
  table$SC_pairing[i] <- (sum(strsplit(SC, "")[[1]] == "(")+ sum(strsplit(SC, "")[[1]] == ")"))/3
}

## Remove duplications
unique_table <- table %>% distinct()
write.csv(unique_table,"forRF7.csv",row.names = FALSE)
