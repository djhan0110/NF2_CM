# Set working environment -------------------------------------------------
rm(list = ls())
cat("\014")
setwd("C:/Workspace/R/NF2_CM/RNA-Seq")

# Load packages -----------------------------------------------------------
library(readr)
library(tidyverse)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggforce)
library(PCAtools)

# Load count matrix -------------------------------------------------------
NF2_Cardio_RCs <- read.csv("00_Raw_data/NF2_Cardio_RCs.txt", sep = "\t", skip = 1)
clean_df <- janitor::clean_names(NF2_Cardio_RCs)

select_df <- clean_df %>%
  dplyr::select(matches("12|geneid"))
select_df <- select_df %>%
  rename_with(~ gsub("x_home_djhan0110_rna_seq_nf2_cardio_03_hisat2_|_s1_l001_r1_bam", "", .)) 

select_df <- as.data.frame(select_df[, order(as.numeric(gsub("\\D", "", colnames(select_df))))]) %>% 
  column_to_rownames("geneid") %>% 
  rename("124" = "Ctrl-1", "125" = "Ctrl-2", "126" = "Ctrl-3",
         "127" = "KO-1", "128" = "KO-2", "129" = "KO-3")

genotype <- factor(rep(c("Ctrl", "KO"), each = 3)) %>% relevel(ref = "Ctrl")

info_df <- data.frame(genotype)
rownames(info_df) <- colnames(select_df)
info_df

select_df[rownames(select_df) %in% "ENSG00000186575", ]

# DESeq2 -------------------------------------------------------------
ddsm <- DESeqDataSetFromMatrix(countData = select_df,
                               colData = info_df,
                               design = ~ genotype)

keep_rows <- rowSums(counts(ddsm)) >= ncol(ddsm)
keep_ddsm <- ddsm[keep_rows, ]
dds <- DESeq(keep_ddsm)
vst_dds <- vst(dds, blind = F)
annotation_vst <- as.data.frame(assay(vst_dds)) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = row.names(vst_dds), keytype = "ENSEMBL", column = "SYMBOL"))

write.csv(annotation_vst, "01_Export/NF2_Cardio_VST.csv")

resultsNames(dds)
resultsNames(dds)[2]

LFC <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm")
annotation_LFC <- as.data.frame(LFC) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = row.names(LFC), keytype = "ENSEMBL", column = "SYMBOL"))
annotation_LFC[annotation_LFC$hgnc_symbol %in% "NF2", ]

write.csv(annotation_LFC, "01_Export/NF2_Cardio_LFC.csv")

# Principal component analysis --------------------------------------------
mapping_vst <- annotation_vst

annotation_symbol <- mapIds(org.Hs.eg.db, keys = rownames(mapping_vst), keytype = "ENSEMBL",column = "SYMBOL")
rownames(mapping_vst) <- make.unique(ifelse(is.na(annotation_symbol), rownames(mapping_vst), annotation_symbol))
mapping_vst <- mapping_vst %>% 
  dplyr::select(-hgnc_symbol)

pca_vst <- pca(mapping_vst, metadata = info_df, removeVar = 0.1)

scree_plot <- screeplot(pca_vst, axisLabSize = 10, titleLabSize = 10)
scree_plot

bi_plot <- biplot(pca_vst, showLoadings = TRUE,
                  labSize = 5, pointSize = 5, sizeLoadingsNames = 5)

bi_plot

pairs_plot <- pairsplot(pca_vst)
pairs_plot

loadings_plot <- plotloadings(pca_vst, labSize = 2)
loadings_plot

precomp_vst <- prcomp(t(as.matrix(assay(vst_dds))))
pc_eigenvalues <- precomp_vst$sdev^2
eigenvalues_df <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))

# Scree plot --------------------------------------------------------------
eigenvalues_df %>%
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) +
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

# Loadings plot -----------------------------------------------------------
annotation_loadings <- data.frame(precomp_vst$rotation) %>% 
  mutate(hgnc_symbol = mapIds(org.Hs.eg.db, keys = rownames(precomp_vst$rotation), keytype = "ENSEMBL", column = "SYMBOL"))
topN = 10
topN_genes <- na.omit(annotation_loadings) %>% 
  dplyr::select(PC1, PC2, hgnc_symbol) %>% 
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>%
  group_by(PC) %>%
  arrange(desc(abs(loading))) %>% 
  dplyr::slice(1:topN) %>% 
  pull(hgnc_symbol) %>%
  unique() 
topN_loadings <- annotation_loadings %>% 
  filter(hgnc_symbol %in% topN_genes)

loadings_plot <- ggplot(topN_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "cm")),
               color = "brown", linewidth = 0.2, alpha = 0.5) +
  geom_text_repel(aes(x = PC1, y = PC2, label = hgnc_symbol),
                  color = ifelse(topN_loadings$hgnc_symbol %in% c(
                    "PLN",
                    "TECRL",
                    "MYL7",
                    "ACTN2",
                    "SMYD1",
                    "TTN",
                    "MYH6",
                    "MYL4",
                    "HSPB7",
                    "MYBPC3"), "red", "black"),
                  min.segment.length = unit(100, "lines"),
                  force = 2, box.padding = 0.1, max.overlaps = 100,
                  fontface = "bold.italic", color = "black", size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02)) +
  scale_y_continuous(expand = c(0.02, 0.02)) +
  labs(x = "PC1", y = "PC2") +
  theme_classic() +
  theme(text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(color = "black", face = "bold", angle = 90),
        axis.text = element_text(color = "black", face = "bold"))
loadings_plot
ggsave("02_Figure/NF2_Cardio_loadings.tiff", loadings_plot, units = "in", width = 2.5, height = 2.5, device = "tiff", dpi = 300)

# PCA plot ----------------------------------------------------------------
PCA_values <- plotPCA(vst_dds, intgroup = "genotype", returnData = T)
percentVar <- round(100 * attr(PCA_values, "percentVar"))

pca_plot <- ggplot(PCA_values, aes(PC1, PC2, fill = genotype)) +
  geom_point(size = 2, stroke = 0.3, shape = 21) +
  geom_mark_ellipse(alpha = 0.1, linewidth = NA, expand = unit(1, "mm")) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        text = element_text(color = "black", face = "bold"),
        axis.text = element_text(color = "black", face = "bold"),
        axis.ticks.length = unit(0.1, "cm"),
        axis.ticks = element_line(linewidth = 0.5, colour = "black"),
        legend.title = element_blank(),
        legend.margin = margin(-2, 0, -5, 0),
        legend.key.size = unit(0.25, "cm"),
        legend.position = "top") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  labs(x = paste0("PC1: ",percentVar[1],"% variance"),
       y = paste0("PC2: ",percentVar[2],"% variance")) 
pca_plot

ggsave("02_Figure/NF2_Cardio_PCA.tiff", pca_plot, units = "in", width = 2, height = 2, device = "tiff", dpi = 300)

