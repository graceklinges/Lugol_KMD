library("phyloseq")
library("microbiome")

rm(list=ls())

setwd("~/Lugol_KMD")

load(file = "ps_unpruned.RData") #unpruned, renamed
ntaxa(ps) #130081

mapfile = "LuKMD_metadata.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map
ps

tax.clean <- data.frame(tax_table(ps))
#write.table(tax.clean, "taxonomy_unrare_fullseq.txt", sep = "\t")

taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps)))
tax.clean <- data.frame(tax_table(ps))
#write.table(tax.clean, "taxonomy_unrare.txt", sep = "\t")

new_names <- paste0("ASV_", seq(ntaxa(ps)))

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){

   if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == "Family_"){
    family <- paste("", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  }
  if (tax.clean[i,6] == "MD3-55"){
    tax.clean[i, 6:7] <- "Aquarickettsia"
  }
}

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == ""){
    family <- paste0("Unclassified_",new_names[i], sep = "") 
    tax.clean[i, 5:7] <- family
  }
}


tax_table(ps) <- as.matrix(tax.clean)
save(ps, file = "ps_unpruned.RData")
#CLR transform using microbiome package for beta div analyses
#CLR transform applies a pseudocount of min(relative abundance)/2 
#to exact zero relative abundance entries in OTU table before taking logs.

rm(list=ls())
load(file = "ps_pruned.RData") #pruned
ntaxa(ps) #102028

mapfile = "LuKMD_metadata.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map
ps

tax.clean <- data.frame(tax_table(ps))
#write.table(tax.clean, "taxonomy_unrare_fullseq.txt", sep = "\t")

taxa_names(ps) <- paste0("ASV_", seq(ntaxa(ps)))
tax.clean <- data.frame(tax_table(ps))
#write.table(tax.clean, "taxonomy_unrare.txt", sep = "\t")

new_names <- paste0("ASV_", seq(ntaxa(ps)))

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == "Family_"){
    family <- paste("", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  }
  if (tax.clean[i,6] == "MD3-55"){
    tax.clean[i, 6:7] <- "Aquarickettsia"
  }
}

for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,6] == ""){
    family <- paste0("Unclassified_",new_names[i], sep = "") 
    tax.clean[i, 5:7] <- family
  }
}

tax_table(ps) <- as.matrix(tax.clean)
save(ps, file = "ps_pruned.RData")
#load("ps_pruned.RData")

summary(sample_sums(ps)) #min 561, median 109194, mean 126192, max 590327

#relative abundance transform
ps_rel <- microbiome::transform(ps, "compositional")
save(ps_rel, file = "ps_rel.RData") #rel abundance, pruned

clr <- microbiome::transform(ps, 'clr')
save(clr, file = "ps_clr.RData") #clr transformed, pruned
