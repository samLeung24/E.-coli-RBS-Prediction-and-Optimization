library(ggplot2)

table <- read.csv("~/Desktop/projectRBS/full_tRNA/gene_abundances.csv")
## Filter undetected mappings
table <- table[table$Gene.ID!=".",]
info <- readLines("../eschColi_K_12_MG1655-mature-tRNAs.fa")
new <- data.frame(vector("character"),vector("character"),vector("character"))
names(new) <- c("Start","End","Codon")

extract_info <- function(input_string) {
  # Extract the start position
  number_pattern <- "chr:(\\d+)-"
  number_matches <- regexpr(number_pattern, input_string, perl = TRUE)
  if (number_matches != -1) {
    number_match <- regmatches(input_string, number_matches)[[1]]
    extracted_number <- as.numeric(gsub("chr:|-", "", number_match))
  } else {
    extracted_number <- NA
  }

  # Extract the end
  number_pattern2 <- "-(\\d+)\\s*\\([+-]"
  number_matches2 <- regexpr(number_pattern2, input_string, perl = TRUE)
  if (number_matches2 != -1) {
    number_match2 <- regmatches(input_string, number_matches2)[[1]]
    extracted_number_dash <- as.numeric(gsub("-|\\s*\\([+-]", "", number_match2))
  } else {
    extracted_number_dash <- NA
  }

  # Extract the codon pattern within the brackets
  codon_pattern <- "\\((\\w{3})\\)"
  codon_matches <- regexpr(codon_pattern, input_string, perl = TRUE)
  if (codon_matches != -1) {
    codon_match <- regmatches(input_string, codon_matches)[[1]]
    extracted_codon <- gsub("[()]", "", codon_match)
  } else {
    extracted_codon <- NA
  }

  # Return a list with the extracted number and codon pattern
  return(c(extracted_number, extracted_number_dash, extracted_codon))
}

convert_codon <- function(anticodon) {
  base1 <- substr(anticodon,1,1)
  base2 <- substr(anticodon,2,2)
  base3 <- substr(anticodon,3,3)

  # Replace each nucleotide base in the anti-codon with its complement using the base-pairing rule
  base1_complement <- if (base1 == "A") {
    "U"
  } else if (base1 == "T") {
    "A"
  } else if (base1 == "C") {
    "G"
  } else {
    "C"
  }

  base2_complement <- if (base2 == "A") {
    "U"
  } else if (base2 == "T") {
    "A"
  } else if (base2 == "C") {
    "G"
  } else {
    "C"
  }

  base3_complement <- if (base3 == "A") {
    "U"
  } else if (base3 == "T") {
    "A"
  } else if (base3 == "C") {
    "G"
  } else {
    "C"
  }

  # Combine the complementary bases to form the codon
  codon <- paste0(base3_complement, base2_complement, base1_complement)

  # Return the codon
  return(codon)
}

for (i in 1:(length(info)/3)) {
  new[i,] <- extract_info(info[3*i-2])
}

codons <- c("UUU", "UUC", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG", "AUU",
            "AUC", "AUA", "AUG", "GUU", "GUC", "GUA", "GUG", "UCU", "UCC",
            "UCA", "UCG", "CCU", "CCC", "CCA", "CCG", "ACU", "ACC", "ACA",
            "ACG", "GCU", "GCC", "GCA", "GCG", "UAU", "UAC", "UAA", "UAG",
            "CAU", "CAC", "CAA", "CAG", "AAU", "AAC", "AAA", "AAG", "GAU",
            "GAC", "GAA", "GAG", "UGU", "UGC", "UGA", "UGG", "CGU", "CGC",
            "CGA", "CGG", "AGU", "AGC", "AGA", "AGG", "GGU", "GGC", "GGA", "GGG")

amino_acids <- c("Phe", "Phe", "Leu", "Leu", "Leu", "Leu", "Leu", "Leu",
                 "Ile", "Ile", "Ile", "Met", "Val", "Val", "Val", "Val",
                 "Ser", "Ser", "Ser", "Ser", "Pro", "Pro", "Pro", "Pro",
                 "Thr", "Thr", "Thr", "Thr", "Ala", "Ala", "Ala", "Ala",
                 "Tyr", "Tyr", "Stop", "Stop", "His", "His", "Gln", "Gln",
                 "Asn", "Asn", "Lys", "Lys", "Asp", "Asp", "Glu", "Glu",
                 "Cys", "Cys", "Stop", "Trp", "Arg", "Arg", "Arg", "Arg",
                 "Ser", "Ser", "Arg", "Arg", "Gly", "Gly", "Gly", "Gly")

codon_table <- data.frame(codons, amino_acids)


merged_tab <- merge(new,table,by = "End", all = FALSE)
merged_tab <- merged_tab[,-which(names(merged_tab) == "Start.y")]
names(merged_tab)[2] <- "Start"
positive_tab <- merged_tab[complete.cases(merged_tab),]
positive_tab$loggedTPM <- log(positive_tab$TPM+0.00001,base=10)
positive_tab$percentTPM <- positive_tab$TPM/10000
colnames(positive_tab)[3] <- "Anticodon"

for (i in 1:length(positive_tab$Anticodon)) {
  positive_tab$Codon[i] <- convert_codon(positive_tab$Anticodon[i])
}

for (j in 1:length(positive_tab$Anticodon)) {
  positive_tab$AA[j] <- codon_table[codon_table$codons==positive_tab$Codon[j],][2]
}

positive_tab$AA <- as.character(positive_tab$AA)


freq <- as.data.frame(table(positive_tab$AA))

ggplot(data = freq, aes(x = Var1,y = Freq))+
  geom_col(fill = "#4285F4")+
  theme_classic()+
  labs(x = "Amino Acids", y = "No. of tRNA")+
  theme(plot.title=element_text(color="Black",size=20),
        plot.subtitle = element_text(color="Black",size=14),
        plot.caption = element_text(color="Black",size=12),
        axis.text = element_text(color="Black",size=14),
        axis.title = element_text(color="Black",size=16))

aa_sum <- aggregate(TPM ~ AA, data = positive_tab, FUN = sum)
aa_sum$"#/1000" <- aa_sum$TPM/1000
aa_sum$logged <- log(aa_sum$TPM,base = 2)

ggplot(data = aa_sum, aes(x = reorder(AA,-`#/1000`),y = `#/1000`))+
  geom_col(fill = "#4285F4")+
  theme_classic()+
  labs(x = "Amino Acids", y = "Codon Abundance (logged)")+
  theme(plot.title=element_text(color="Black",size=20),
        plot.subtitle = element_text(color="Black",size=14),
        plot.caption = element_text(color="Black",size=12),
        axis.text = element_text(color="Black",size=12),
        axis.text.x = element_text(color = "Black", size = 12, angle = 45, vjust = 0.7 ),
        axis.title = element_text(color="Black",size=14))

codon_availability <- aggregate(TPM~Codon,data = positive_tab, sum)


write.csv(positive_tab,"../full_tRNA_abundance.csv",row.names = FALSE)
write.csv(codon_availability,"../codon_availability.csv",row.names = FALSE)
