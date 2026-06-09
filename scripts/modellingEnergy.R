## Energy Model
library(stringr)
library(ggplot2)

table <- read.csv("forRF8.csv")
nRow <- length(table$Locus)

## Algorith for spacing compensation
calculate_spacing_free_energy <- function(s) {
  ## reference
  if (s>=5) deltaG = 0.048*(s-5)^2+0.24*(s-5)
  else deltaG = 12.2/(1+exp(2.5*(s-5+2)))^3
  return(deltaG)
}

sc_tab <- table[,c("Locus","ORF")]
aug <- -2.47
gug <- 0.62
uug <- 1.23
sc_tab$sc <- as.character(RNAStringSet(DNAStringSet(substr(sc_tab$ORF,1,3))))

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

mRNA <- readLines("../30_30_RNAfold.dbn")
mRNA_info <- data.frame(Locus =rep("",nGenes),
                        structure = rep("",nGenes),
                        energy = rep(0,nGenes))
for (i in 1:nGenes) {
  print(i)
  name <- mRNA[3*i-2]
  input_structure <- mRNA[[3*i]]
  mRNA_info$Locus[i] <- strsplit(name,"\\|")[[1]][2]
  mRNA_info$structure[i] <- strsplit(input_structure," ")[[1]][1]
  mRNA_info$energy[i] <- gsub("[)]","",strsplit(input_structure,"\\s\\(")[[1]][2])
}

energy_tab <- data.frame(Locus = table$Locus,
                         hybrid = rep(0,nRow),
                         mRNA = rep(0,nRow),
                         start = rep(0,nRow),
                         spacing = rep(0,nRow),
                         total = rep(0,nRow))
for (i in 1:nRow){
  energy_tab$hybrid[i] <- as.numeric(hybrid_info[hybrid_info$Locus == energy_tab$Locus[i],"energy"])
  energy_tab$mRNA[i] <- as.numeric(mRNA_info[mRNA_info$Locus == energy_tab$Locus[i],"energy"])
  sc <- sc_tab[sc_tab$Locus == energy_tab$Locus[i],"sc"][1]
  if (sc == "AUG") {
    energy_tab$start[i] <- aug
  } else if (sc == "GUG")  {
    energy_tab$start[i] <- gug
  } else {
    energy_tab$start[i] <- uug
  }
  spacer_length <- 30-as.numeric(hybrid_info[hybrid_info$Locus == energy_tab$Locus[i],"end"])
  energy_tab$spacer_length[i] <- spacer_length
  energy_tab$spacing[i] <- calculate_spacing_free_energy(s = spacer_length)
  energy_tab$total <- energy_tab$hybrid+energy_tab$start+energy_tab$spacing-energy_tab$mRNA
}

table$energy_total <- energy_tab$total

## Convert SC to corresponding energy
for (i in 1:nRow) {
  sc <- sc_tab$sc[i]
  if (sc == "AUG") {
    table$sc_energy[i] <- aug
  } else if (sc == "GUG")  {
    table$sc_energy[i] <- gug
  } else {
    table$sc_energy[i] <- uug
  }
}

## Convert SC to corresponding cai
sc_freq <- table(sc_tab$sc)
for (i in 1:nRow) {
  sc <- sc_tab$sc[i]
  if (sc == "AUG") {
    table$sc_cai[i] <- 1
  } else if (sc == "GUG")  {
    table$sc_cai[i] <- 0.496
  } else if (sc == "UUG") table$sc_cai[i] <- 0.106
}

write.csv(table,"forRF9.csv",row.names = FALSE)
