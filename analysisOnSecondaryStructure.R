library(ggplot2)
### Retrieve expression level data
exptab <- read.table("../511145-WHOLE_ORGANISM-integrated.txt",sep = "\t")
ratiotab <- read.csv("../edited_processed_sequences_pax_dna_prot_Ecoli-K12_20230118.csv",skip = 1)
for (i in 1:length(exptab$V2)) {
  exptab$Locus.tag[i] <- strsplit(exptab$V2[i],".", fixed =TRUE)[[1]][2]
}
temptab <- cbind(gene_id=ratiotab$gene_id,Locus.tag=ratiotab$Locus.tag)
finaltab <- merge(temptab,exptab[,c("V3","Locus.tag")],by = "Locus.tag")  
positives <- read.csv("bulkProdigal.csv")
positives <- positives[complete.cases(positives),]
finaltab <- finaltab[finaltab$gene_id %in% positives$gene_id,]
colnames(finaltab) <- c("Locus.tag", "gene_id", "expression_level")

### Import secondary structure predictions (Mxfold2)
sectab <- read.csv("proportionPairedMxfold2.csv")
table <- merge(sectab,finaltab, by = "gene_id")
table <- merge(table,positives,by = "gene_id")
table$loged_el <- log(table$expression_level,base = 10)
table <- table[is.finite(table$loged_el),]

plot(table$spacer_length,table$loged_el,main = "Expression Level v. Spacer Length (30)")
abline(lm(table$loged_el~table$spacer_length),col = "red")
summary(lm(table$spacer_length~table$loged_el))


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

Coefficients <- vector(length = 100,mode = "numeric")

for (i in 1:100) {
  R = somekindoffunction(i)
  Coefficients[i] <- R
}

print(Coefficients)
write.csv(Coefficients,"a.txt")

