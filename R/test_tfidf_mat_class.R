devtools::install_github("diogocamacho/druid", force=TRUE)
library(DRUID)

tfidf_mat = DRUID::cmap_druid$tfidf

print(class(tfidf_mat))
print(class(tfidf_mat) != 'dgCMatrix') # results in FALSE
print(class(tfidf_mat) != 'matrix') # results in TRUE - this is bad
print((class(tfidf_mat) != "dgCMatrix") | (class(tfidf_mat) != "matrix")) # results in TRUE - also bad
print((class(tfidf_mat) != "dgCMatrix") & (class(tfidf_mat) != "matrix")) # results in FALSE - this is what we want
