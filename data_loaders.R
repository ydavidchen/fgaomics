# Raw Data Loaders: cBioPortal clinical & Synapse samples+RNAseqV2
# Use in conjunction with utils.R
# Ref. gdc.cancer.gov/resources-tcga-users/tcga-code-tables

library(data.table)

CANCERS <- c(
  "COAD","READ", 
  "GBM","LGG", 
  "KIRC","KIRP", 
  "STAD","ESCA", 
  "LUAD","LUSC", 
  
  "LIHC", "PAAD",
  "BRCA", "OV",
  "BLCA",
  "HNSC",
  "PRAD",
  "SKCM"
)

COLS_CBIO <- c("Patient ID","Sample ID", 
               "TCGA PanCanAtlas Cancer Type Acronym",
               "Cancer Type", "Cancer Type Detailed", "Subtype",
               "Diagnosis Age", "Sex",
               "Mutation Count","TMB (nonsynonymous)","Fraction Genome Altered")

load_clinical <- function(cancer) {
  #'@description Loads & joins cBioPortal & Synapse sample files
  
  cbio <- read.csv(paste0(DIR_CBIO, tolower(cancer), "_tcga_pan_can_atlas_2018_clinical_data.tsv"), sep="\t", check.names=FALSE)
  cbio <- cbio[ , COLS_CBIO]
  
  samples <- read.csv(paste0(DIR_SYN, cancer, "/", "nationwidechildrens.org_", cancer, "_bio.sample.tsv"), sep="\t")
  samples <- samples[ , c("sample","sample_type")]
  samples$`Sample ID` <- substr(samples$sample, 1, 15)
  
  res <- merge(samples, cbio, by="Sample ID")
  colnames(res)[colnames(res)=="TCGA PanCanAtlas Cancer Type Acronym"] <- "TCGA"
  colnames(res)[colnames(res)=="Diagnosis Age"] <- "AgeAtDx"
  colnames(res)[colnames(res)=="Fraction Genome Altered"] <- "FGA"
  colnames(res)[colnames(res)=="TMB (nonsynonymous)"] <- "TMB"
  
  res$normFga <- MinMaxScale(res$FGA)
  res$normTmb <- MinMaxScale(res$TMB)
  res$highFga <- res$normFga > median(res$normFga, na.rm=TRUE)
  res$highTmb <- res$normTmb > median(res$normTmb, na.rm=TRUE)
  
  return(res)
}


load_expr_as_df <- function(cancer, prefix="/unc.edu_", suffix="_IlluminaHiSeq_RNASeqV2.geneExp.tsv", 
                            excludeUnchar=TRUE, normalize=TRUE, propFilter=0.10, df=TRUE) {
  #'@description Loads & prelim-processes Synapse RNAseqV2
  expr <- fread(paste0(DIR_SYN, cancer, prefix, cancer, suffix), data.table=FALSE)
  expr <- subset(expr, ! grepl("?", gene_id, fixed=TRUE))
  
  if(excludeUnchar) {
    expr <- subset(expr, ! grepl("C\\d+orf", gene_id, ignore.case=FALSE))
    expr <- subset(expr, ! grepl("^LOC\\d+", gene_id, ignore.case=FALSE))
  }
  
  rownames(expr) <- expr$gene_id
  expr$gene_id <- NULL
  expr <- expr[ , substr(colnames(expr), 14, 15) %in% c("01","02","06")]
  expr <- data.matrix(expr)
  
  if(propFilter > 0) expr <- select_most_var(expr, round(nrow(expr)*(1-propFilter)))
  
  if(normalize) expr <- t(scale(t(expr))) #column-wise
  
  if(df) {
    expr <- as.data.frame(expr)
    expr$gene_id <- rownames(expr)
    rownames(expr) <- NULL
    expr <- expr[ , c(ncol(expr),1:(ncol(expr)-1))]
  }
  
  return(expr)
}
