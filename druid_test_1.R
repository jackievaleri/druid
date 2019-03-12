## testing with version 1.1.0

devtools::install_github("diogocamacho/druid", force=TRUE)
library(DRUID)

# diogo installation example
gset <- unique(gsub(" down", "", gsub(" up", "", sample(colnames(DRUID::cmap_druid$tfidf), 100))))

query_matrix <- matrix(1, ncol = 2, nrow = length(gset))
query_matrix[, 2] <- 0
query_matrix[sample(x = seq(1, length(gset)), size = 0.25 * length(gset)), 1] <- -1

tfidf_mat = DRUID::cmap_druid$tfidf
print(class(tfidf_mat))

example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = tfidf_mat, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)
# problems: cant leave out crossproduct without getting "Error in concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf,  : 
# TF-IDF matrix is of wrong type. Please revise.""
# also, example concoct wrapper needs to be changed to leave out crossproduct in README

print(example_druid)

example_druid <- example_druid %>% 
  tibble::add_column(., drug_name = DRUID::cmap_druid$drugs$name, .before = 1) %>%
  tibble::add_column(., concentration = DRUID::cmap_druid$drugs$concentration, .before = 2) %>%
  tibble::add_column(., cell_line = DRUID::cmap_druid$drugs$cell_line, .before = 3)

print(example_druid)

library(ggplot2)
example_druid %>% dplyr::filter(., cosine_similarity == 0) %>% ggplot() + geom_point(aes(x = drug_name, y = druid_score, color = cell_line), alpha = 0.5) + facet_grid(. ~ cell_line, scales = "free") + theme_bw() + theme(axis.text.x = element_blank())

# test dge_matrix argument
query_matrix_empty <- matrix(0, ncol = 2, nrow = length(gset))
query_matrix_empty[, 2] <- 0
example_druid <- concoct(dge_matrix = query_matrix_empty, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)

query_matrix_ones <- matrix(1, ncol = 2, nrow = length(gset))
query_matrix_ones[, 2] <- 0
example_druid <- concoct(dge_matrix = query_matrix_ones, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)
print(example_druid)

# test tfidf_matrix argument
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = "", tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)

# test crossproduct argument with empty string
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = "", num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)
# test crossproduct argument with empty vector
empty_cross_test <- vector(mode="numeric", length=0)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = empty_cross_test, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)

# test num_random argument
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 0, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = -1, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 1000.0, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 1000.5, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)

# test druid_direction argument
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "pos", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "fake", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)

# test fold_thr argument
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 1, pvalue_thr = 0.05, entrez = gset)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = -1, pvalue_thr = 0.05, entrez = gset)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0.5, pvalue_thr = 0.05, entrez = gset)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 2, pvalue_thr = 0.05, entrez = gset)

# test pvalue_thr argument
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.01, entrez = gset)
print(example_druid)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.001, entrez = gset)
print(example_druid)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.00001, entrez = gset)
print(example_druid)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 1, entrez = gset)
print(example_druid)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0, entrez = gset)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = -1, entrez = gset)

# test entrez argument
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = "")
gset = as.integer(gset)
example_druid <- concoct(dge_matrix = query_matrix, tfidf_matrix = DRUID::cmap_druid$tfidf, tfidf_crossproduct = DRUID::cmap_druid$cpm, num_random = 10000, druid_direction = "neg", fold_thr = 0, pvalue_thr = 0.05, entrez = gset)
