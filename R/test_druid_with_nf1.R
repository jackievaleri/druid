# Test concoct with diff exp genes from nf1
setwd("~/Documents/Lab/druid/R")

# import statements
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
require(biomaRt)
library(biomaRt)
#devtools::install_github("diogocamacho/druid")
library(DRUID)
#devtools::install_github("diogocamacho/cauldron")

# load in diff exp genes
diff_exp = read.csv("../../nf1/out/all_significant_4_30_19.csv")
lg2fc = diff_exp[,3]
pval = diff_exp[,6]
dge_mat = cbind(lg2fc, pval)

# take off the decimal points from ensembl ids
gene_names = diff_exp[,1]
num_genes = length(gene_names)
new_names = {num_genes * 0}

# strip period
for (i in 1:num_genes) {
  current_line = as.character(gene_names[i])
  gene_split = strsplit(current_line, "[.]")[[1]]
  # note the [.] instead of . because the split is a regexp
  new_names[i] = gene_split[1]
}

# convert ensembl IDs using biomart
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# get gene symbols instead of ensembl ids
map = getBM(attributes= c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=new_names, mart=mart)
map = map[match(new_names, map$ensembl_gene_id),] # reorder to match input
entrez_nf1_genes = map$hgnc_symbol

# try druid
druid_results = DRUID::concoct(dge_matrix = dge_mat, entrez = entrez_nf1_genes)
druid_results = DRUID::concoct(dge_matrix = dge_mat, entrez = entrez_nf1_genes, druid_direction = "pos")

# ERRORS: 
# RUN with CMAP
# make clear you choose a number
# download cauldron too

# Selection: 1
#:: Running DRUID on CMAP data ::
# Generating query vector...
#Error in druid_geneset(dge_matrix = dge_matrix, desired_effect = druid_direction,  : 
       #                  No genes macthed in drug profiles.

#Selection: 6
#:: Running DRUID on all data sets ::
# !!Warning: depending on processor speed, this could take a while. Go get a coffee or something.

#CMAP data...
#Generating query vector...
#Error in druid_geneset(dge_matrix = dge_matrix, desired_effect = druid_direction,  : 
 #                        No genes macthed in drug profiles.


## Okay.... maybe I need all genes and druid can decide what it needs

# load in all exp genes
diff_exp = read.csv("../../nf1/out/all_genes_4_30_19.csv")
lg2fc = diff_exp[,3]
pval = diff_exp[,6]
dge_mat = cbind(lg2fc, pval)

# take off the decimal points
gene_names = diff_exp[,1]
num_genes = length(gene_names)
new_names = {num_genes * 0}

# strip period
for (i in 1:num_genes) {
  current_line = as.character(gene_names[i])
  gene_split = strsplit(current_line, "[.]")[[1]]
  # note the [.] instead of . because the split is a regexp
  new_names[i] = gene_split[1]
}

# get gene symbols instead of ensembl ids
map = getBM(attributes= c("ensembl_gene_id", "hgnc_symbol"), filters="ensembl_gene_id", values=new_names, mart=mart)
map = map[match(new_names, map$ensembl_gene_id),] # reorder to match input
entrez_nf1_genes = map$hgnc_symbol

# test druid
druid_results = DRUID::concoct(dge_matrix = dge_mat, entrez = entrez_nf1_genes)
druid_results = DRUID::concoct(dge_matrix = dge_mat, entrez = entrez_nf1_genes, druid_direction = "pos")
druid_results = DRUID::concoct(dge_matrix = dge_mat, entrez = entrez_nf1_genes, pvalue_thr = 0.25)
druid_results = DRUID::concoct(dge_matrix = dge_mat, entrez = entrez_nf1_genes, pvalue_thr = 1, druid_direction = "pos")
