# Load, Handle, & Stack clinical+sample files of ALL cancers chosen

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
source("data_loaders.R")

CANCERS <- c("COAD","READ", "GBM","LGG", "KIRC","KIRP", "STAD","ESCA", "LUAD","LUSC", "BRCA")

## Stacking clinical:
meta_clin <- NULL
for(cancer_k in CANCERS) {
  tmp <- load_clinical(cancer=cancer_k)
  meta_clin <- rbind(meta_clin, tmp)
}
rm(tmp, cancer_k)

meta_clin <- subset(meta_clin, ! is.na(FGA))

meta_clin$Macro <- NA
meta_clin$Macro[meta_clin$TCGA %in% c("GBM","LGG")] <- "GBMLGG"
meta_clin$Macro[meta_clin$TCGA %in% c("COAD","READ")] <- "COADREAD"
meta_clin$Macro[meta_clin$TCGA %in% c("LUSC","LUAD")] <- "LUNG"
meta_clin$Macro[meta_clin$TCGA %in% c("STAD","ESCA")] <- "STES"
meta_clin$Macro[meta_clin$TCGA %in% c("KIRP","KIRC")] <- "KIR"
meta_clin$Macro[meta_clin$TCGA == "BRCA"] <- "BRCA"

table(meta_clin$Macro, useNA="always")
table(substr(meta_clin$sample, 14, 15)) #check tumor sample types: 7313    9  363 

## Stacking expression via joining with pre-processing:
mega_expr <- load_expr_as_df(CANCERS[1]) #init

for(cancer_k in CANCERS[-1]) {
  tmp <- load_expr_as_df(cancer=cancer_k)
  mega_expr <- merge(mega_expr, tmp, by="gene_id")
}
rm(tmp, cancer_k)

rownames(mega_expr) <- mega_expr$gene_id
mega_expr$gene_id <- NULL
mega_expr <- data.matrix(mega_expr)

## Post Processing:
## For the very few (n=2) repeat samples at **vial** level, keep 1 only:
sum(duplicated(substr(colnames(mega_expr), 1, 16)))
mega_expr <- mega_expr[ , ! duplicated(substr(colnames(mega_expr), 1, 16))]
colnames(mega_expr) <- substr(colnames(mega_expr), 1, 16)

mega_expr <- mega_expr[ ,colnames(mega_expr) %in% meta_clin$sample]
meta_clin <- subset(meta_clin, sample %in% colnames(mega_expr))

stopifnot(identical(meta_clin$sample, colnames(mega_expr)))

## Export:
# save(list=c("meta_clin","mega_expr"), file="../results/240217a_mega.RData", compress=TRUE)

## Export for python ML/DL:
mega_expr <- as.data.frame(t(mega_expr))
mega_expr$sample <- rownames(mega_expr)
rownames(mega_expr) <- NULL

mega_expr <- merge(meta_clin[,c("sample","TCGA","normFga","highFga")], mega_expr, by="sample")
# write.csv(mega_expr, file="../results/240217b_df4ml.csv", row.names=FALSE, quote=FALSE)
