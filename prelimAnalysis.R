library(ggplot2)

### Correlations between GC-content, length of spacer, and length of SD
table <- read.csv("bulkProdigal.csv")
table2 <- read.csv("../processed_sequences_pax_dna_prot_Ecoli-K12_20230118.csv")

## Merging the abundance column with the SD Prediction Table
tab_need <- cbind(table2$GeneID,table2$abundance_top_ratio)
colnames(tab_need) <- c("gene_id","abundance")
tab_merged <- merge(table,tab_need,by="gene_id",all.x=TRUE)
tab_merged <- tab_merged[complete.cases(tab_merged),]
tab_merged$sd_length <- rep(NA,nrow(tab_merged))
for (i in 1:nrow(tab_merged)){
  tab_merged$sd_length[i] <- nchar(tab_merged$sd_seq[i])
}

class(tab_merged[,9])

sd_length <- rep(NA,nrow(table))
sd_seq <- table$sd_seq
table$spacer_length[is.na(table$spacer_length)] <- 0
for (i in 1:nrow(table)) {
  if (is.na(sd_seq[i])){
    sd_length[i] <- 0
  } else {
    sd_length[i] <- nchar(sd_seq[i])
  }
}

table$sd_length <- sd_length

for (i in 1:nrow(table)) {
  if (is.na(table$sd_seq[i])){
    table$sd_seq[i] <- 0
  } else {
    table$sd_seq[i] <- 1
  }
}

table$sd_seq <- as.numeric(table$sd_seq)

## Calculate the association between GC-content and SD occurence
cor(table$gc_content, table$sd_seq)
model <- lm(sd_seq ~ gc_content, data = table)
summary(model)
plot(table$gc_content, table$sd_seq)

# Checking for collinearity, not suitable for multi-variable
cor(tab_merged$sd_length,tab_merged$spacer_length)
plot(tab_merged$sd_length,tab_merged$spacer_length)
abline(lm(tab_merged$sd_length ~ tab_merged$spacer_length))
summary(lm(tab_merged$sd_length ~ tab_merged$spacer_length))

## Calculate the association spacer_length and abundance
cor(tab_merged$spacer_length,y=tab_merged$abundance)
model <- lm(abundance ~ spacer_length, data = tab_merged)
plot(tab_merged$spacer_length,y=tab_merged$abundance)
abline(model)
summary(model)

## Calculate the associationsd_length and abundance
cor(tab_merged$sd_length,y=tab_merged$abundance)
model <- lm(abundance ~ sd_length, data = tab_merged)
plot(tab_merged$sd_length,y=tab_merged$abundance)
abline(model)
summary(model)




