# Test concoct with diff exp genes from nf1

setwd("~/Documents/Lab/druid/R")

# load in diff exp genes
diff_exp = read.csv("../../nf1/out/all_significant_4_30_19.csv")
lg2fc = diff_exp[,3]
pval = diff_exp[,6]
dge_mat = cbind(lg2fc, pval)

# convert ensembl IDs to entrez names
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
require(biomaRt)

devtools::install_github("diogocamacho/druid")
library(DRUID)
devtools::install_github("diogocamacho/cauldron")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

listAttributes(mart)
listFilters()

#getBM(attributes= "hgnc_symbol", filters="ensembl_gene_id", values=gene_names, mart=mart)

# take off the decimal points
gene_names = diff_exp[,1]
num_genes = length(gene_names)
new_names = {num_genes * 0}

for (i in 1:num_genes) {
  current_line = as.character(gene_names[i])
  gene_split = strsplit(current_line, "[.]")[[1]]
  # note the [.] instead of . because the split is a regexp
  print(gene_split[1])
  new_names[i] = gene_split[1]
}

entrez_gene_names = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart=mart)
entrez_nf1_map=  entrez_gene_names[which(entrez_gene_names$ensembl_gene_id %in% new_names),]
# okay that gives me the genes but not in the right order.....
entrez_nf1_genes = {num_genes * 0}
for (i in 1:num_genes) {
  curr_gene = new_names[i]
  index_curr_gene = which(entrez_nf1_map$ensembl_gene_id == curr_gene)
  entrez_nf1_genes[i] = entrez_nf1_map[index_curr_gene,2]
  print(entrez_nf1_genes[i])
}

druid_results = DRUID::concoct(dge_matrix = dge_mat, entrez = entrez_nf1_genes)
druid_results = DRUID::concoct(dge_matrix = dge_mat, entrez = entrez_nf1_genes, druid_direction = "pos")


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


## Okay.... maybe I need all genes

# load in all exp genes
diff_exp = read.csv("../../nf1/out/all_genes_4_30_19.csv")
lg2fc = diff_exp[,3]
pval = diff_exp[,6]
dge_mat = cbind(lg2fc, pval)

# take off the decimal points
gene_names = diff_exp[,1]
num_genes = length(gene_names)
new_names = {num_genes * 0}

for (i in 1:num_genes) {
  current_line = as.character(gene_names[i])
  gene_split = strsplit(current_line, "[.]")[[1]]
  # note the [.] instead of . because the split is a regexp
  #print(gene_split[1])
  new_names[i] = gene_split[1]
}

# get entrez id names
entrez_gene_names = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart=mart)
entrez_nf1_map=  entrez_gene_names[which(entrez_gene_names$ensembl_gene_id %in% new_names),]
# okay that gives me the genes but not in the right order.....
entrez_nf1_genes = {num_genes * 0}
for (i in 1:num_genes) {
  curr_gene = new_names[i]
  index_curr_gene = which(entrez_nf1_map$ensembl_gene_id == curr_gene)
  print(curr_gene)
  replacement = entrez_nf1_map[index_curr_gene,2]
  print(entrez_nf1_map[index_curr_gene,2])
  if (length(replacement) > 0) {
    if (replacement != "") {
    entrez_nf1_genes[i] = entrez_nf1_map[index_curr_gene,2]
    }
  }
}

# test druid
druid_results = DRUID::concoct(dge_matrix = dge_mat, entrez = entrez_nf1_genes)
druid_results = DRUID::concoct(dge_matrix = dge_mat, entrez = entrez_nf1_genes, druid_direction = "pos")
druid_results = DRUID::concoct(dge_matrix = dge_mat, entrez = entrez_nf1_genes, pvalue_thr = 0.25)
druid_results = DRUID::concoct(dge_matrix = dge_mat, entrez = entrez_nf1_genes, pvalue_thr = 1, druid_direction = "pos")

