library("dada2")
library("seqinr")
library("biomformat")
library("data.table")
library("ggplot2")
library("ggthemes")
library("ampvis2")
library("cowplot")
library("phyloseq")

rm(list=ls())
setwd("~/Lugol_KMD/") 

# functions
fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

# phyloseq object output from Dada2
ps <- readRDS("ps_KMD_Lugol.rds")
seqtab <- readRDS("seqtab_chimerafree_merged.rds")
# sequences object made from Dada2

# exporting sequences to a fasta file for import into qiime
uniqueSeqs <- as.list(colnames(seqtab))
write.fasta(uniqueSeqs, uniqueSeqs, "uniqueSeqs.fasta")

# import metadata and merge into phyloseq object
mapfile = "LuKMD_metadata.txt"
map = import_qiime_sample_data(mapfile)
sample_data(ps) <- map

# export taxonomy to import into qiime2
tax <-as(tax_table(ps),"matrix")
tax_cols <- c("Kingdom", "Phylum", "Class", "Order","Family","Genus", "Species")
tax <-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, "taxonomy.txt", quote=FALSE, col.names=FALSE, sep="\t")

# summary of data
#ps
#summary(sample_data(ps))

ntaxa(ps) #138596
nsamples(ps) #400
rank_names(ps)
sample_names(ps)[1:5]
sample_variables(ps)

# remove mitochondria and chloroplasts, is.na important becuase if not included
# this command will also remove all Family = NA or Order = NA
ps_full <- ps
ps_with_mito = subset_taxa(ps, (Order!="Chloroplast") | is.na(Order))
ps_no_mito = subset_taxa(ps_with_mito, (Family!="Mitochondria") | is.na(Family))
ps_no_Eukaryota = subset_taxa(ps_no_mito, (Kingdom!="Eukaryota") | is.na(Kingdom))

ntaxa(ps_full) - ntaxa(ps_with_mito) #7475 taxa lost
ntaxa(ps_with_mito) - ntaxa(ps_no_mito)#lost 1040 taxa
ntaxa(ps_no_mito) - ntaxa(ps_no_Eukaryota)#0 taxa lost

ps = ps_no_Eukaryota
ntaxa(ps) #130081
save(ps, file = "ps_unpruned.RData") #dataset with no pruning but contaminants removed

summary(taxa_sums(ps)) # 18 was first quantile 
# Filter on prevalence or total counts
pst = fast_melt(ps)
write.table(pst, file = "pst_LuKMD.txt", sep = "\t")

#remove ASVs with total counts in the bottom quartile
prevdt = pst[, list(Prevalence = sum(count > 0), 
                    TotalCounts = sum(count)),
             by = taxaID]
keepTaxa = prevdt[(Prevalence >=0 & TotalCounts >15), taxaID]
ps_pruned = prune_taxa(keepTaxa,ps)
ps_pruned
sample_sums(ps_pruned)
ntaxa(ps_pruned) #102028
min(sample_sums(ps_pruned)) #561 woof

summary_tab <- data.frame(init=(as.matrix(sample_sums(ps_full)))[,1],
                          chloros_removed=(as.matrix(sample_sums(ps_with_mito)))[,1],
                          mitos_removed=(as.matrix(sample_sums(ps_no_mito)))[,1],
                          pruned=(as.matrix(sample_sums(ps_pruned)))[,1])
write.table(summary_tab, "reads_lost_phyloseq.txt", quote=FALSE, col.names=TRUE, sep="\t")

#going to use my pruning step
ps <- ps_pruned
save(ps, file = "ps_pruned.RData")
