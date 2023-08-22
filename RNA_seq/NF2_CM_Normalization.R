# Set working environment -------------------------------------------------
cat('\014')
rm(list = ls())
setwd("C:/Workspace/R/Project/NF2_CM")

# Load library ------------------------------------------------------------
library(readr)
library(tidyverse)
library(janitor)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Load count matrix -------------------------------------------------------
RCs_iPSC_Wi_HiN <- read_csv("Rawdata/RCs_iPSC_Wi_HiN.csv")
rc_data <- as.data.frame(clean_names(RCs_iPSC_Wi_HiN))
head(rc_data)
rownames(rc_data) <- rc_data[["x1"]]
cln_data <- rc_data[, -which(names(rc_data) == "x1")]

rc_cols <- c("Wi-1", "Wi-2", "Wi-3",
             "HiN-1", "H11-2", "H11-3")
colnames(cln_data) <- rc_cols
head(cln_data)

group_info <- c(rep('WT', 3), rep('KO', 3))
info_df <- data.frame(sample_info = colnames(cln_data), group_info)
info_df$group_info <- factor(info_df$group_info, levels = c("WT", "KO"))
info_df

# DESeq2 ------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = cln_data,
                               colData = info_df,
                               design = ~ group_info)
keep <- rowSums(counts(dds)) >=  ncol(dds)
dds <- dds[keep, ]
ddsDE <- DESeq(dds)
ncs <- as.data.frame(counts(ddsDE, normalized = T))
anno_ncs <- ncs %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = rownames(ncs), keytype = "ENSEMBL", column = "SYMBOL"))
head(anno_ncs)
#write.csv(anno_ncs, "Export/NCS_NF2_CM_iPSC.csv")
vsd_ncs <- as.data.frame(assay(vst(dds, blind = F)))
anno_vsd <- vsd_ncs %>%
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = rownames(vsd_ncs), keytype = "ENSEMBL", column = "SYMBOL"))
head(anno_vsd)
#write.csv(anno_vsd, "Export/VST_NF2_CM_iPSC.csv")
resultsNames(ddsDE)
res_DE <- results(ddsDE, name = 'group_info_KO_vs_WT')
anno_res <- as.data.frame(res_DE) %>%
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = rownames(res_DE), keytype = "ENSEMBL", column = "SYMBOL"))
head(anno_res)
#write.csv(anno_res, "Export/RES_NF2_CM_iPSC.csv")
anno_res[anno_res$hgnc_symbol %in% "NF2", ]

lfc_DE <- lfcShrink(ddsDE, coef = 'group_info_KO_vs_WT', type = 'apeglm')
anno_lfc <- as.data.frame(lfc_DE) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = rownames(lfc_DE), keytype = "ENSEMBL", column = "SYMBOL"))
head(anno_lfc)
write.csv(anno_lfc, "Export/LFC_NF2_CM_iPSC.csv")
anno_lfc[anno_lfc$hgnc_symbol %in% "NF2", ]
