library(dada2); packageVersion("dada2")
library(phyloseq)

#runs need to have algorithm run separately
setwd("~/Desktop/LugolKMD/") #400 samples, DADA2 processing and chimera removal done separately 
rm(list=ls())
## starting from sequence tables for lane 1 and lane 2 saved 
st.1<-readRDS("seqtab_Lugol.rds")
st.2<-readRDS("seqtab-KMD.rds")

# Inspect distribution of read lengths
table(nchar(getSequences(st.1)))
table(nchar(getSequences(st.2)))

#merge the two sequence tables for lane 1 and 2
st.all<-mergeSequenceTables(st.1,st.2)
table(nchar(getSequences(st.all)))
hist(nchar(getSequences(st.all)))

seqtab_trimmed <- st.all 
ncol(seqtab_trimmed) #138596

names_seqtab <- row.names(seqtab_trimmed)
write.table(names_seqtab, file = "names_LugolKMD_combined.txt", sep="\t")

# Save chimera-free ASV table as downstream tasks may cause R to crash
saveRDS(seqtab_trimmed, "seqtab_chimerafree_merged.rds")

# Assign taxonomy based on silva reference database at species (100%) level, you must have the appropriate Silva database downloaded
tax_silva <- assignTaxonomy(seqtab_trimmed, "~/Desktop/Databases/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)

library(phyloseq)
# Export sequence table with genus and species assignments as phyloseq objects
ps <- phyloseq(otu_table(seqtab_trimmed, taxa_are_rows=FALSE), tax_table(tax_silva))

# Save as RDS objects
saveRDS(ps, file = "ps_KMD_Lugol.rds")
