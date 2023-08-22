# Set working environment -------------------------------------------------
rm(list = ls())
current_path <- getwd()
#setwd("C:/Workspace/R/Project/NF2_CM")

# Load libraries ----------------------------------------------------------
library(readr)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(ggplot2)
library(ComplexHeatmap)
library(colorRamp2)
library(seriation)
library(gridExtra)

# Load result -------------------------------------------------------------
GOIDs <- read.table(paste0(current_path, "/RNA_seq/Rawdata/GOIDs_NF2_CM_iPSC.txt"), header = T, colClasses = "character")
GOIDs$IDs <- paste0("GO:", GOIDs$IDs)
anno_GO <- AnnotationDbi::select(GO.db, keys = GOIDs$ID, keytype = "GOID",  columns = "TERM") %>% 
  mutate(TERM = str_to_sentence(TERM))
ret_genes <- AnnotationDbi::select(org.Hs.eg.db, keys = anno_GO$GOID, keytype = "GOALL",  columns=c("ENSEMBL", "SYMBOL"))
GO_genes <- left_join(ret_genes, anno_GO, by = c("GOALL" = "GOID"))

LFC_NF2_CM_iPSC <- read_csv(paste0(current_path, "/RNA_seq/Export/LFC_NF2_CM_iPSC.csv"))
LFC <- as.data.frame(janitor::clean_names(LFC_NF2_CM_iPSC)) %>% 
  mutate(nlogpadj = -log10(padj),
         attr = ifelse(padj < 0.05 & log2fold_change > 1, "up",
                       ifelse(padj < 0.05 & log2fold_change < -1, "down", "ns")))
filtered_LFC <- LFC %>% 
  filter(attr != "ns") %>% na.omit()
table(LFC$attr)
table(filtered_LFC$attr)

VST_NF2_CM_iPSC <- read_csv(paste0(current_path, "/RNA_seq/Export/VST_NF2_CM_iPSC.csv"))
VST <- as.data.frame(janitor::clean_names(VST_NF2_CM_iPSC)) %>% 
  left_join(LFC %>% dplyr::select(x1, log2fold_change, padj, nlogpadj, attr), by = "x1") %>% 
  na.omit()

merged_data <- left_join(GO_genes, filtered_LFC %>% dplyr::select(x1, attr), by = c("ENSEMBL" = "x1")) %>% 
  na.omit()

GO_list <- split(merged_data, merged_data$GOALL)
cln_list <- lapply(GO_list, function(ls) ls[!duplicated(ls$SYMBOL), ])
comb_df <- do.call(rbind, cln_list) 
uniq_GO <- unique(comb_df$GOALL)
uniq_symbol <- unique(comb_df$SYMBOL)
filtered_VST <- VST %>% 
  filter(hgnc_symbol %in% uniq_symbol)
row.names(filtered_VST) <- filtered_VST$x1

for (i in 1:length(uniq_GO)) {
  current_GO <- as.character(uniq_GO[i])
  current_SYMBOL <- GO_list[[current_GO]]$SYMBOL
  filtered_VST <- filtered_VST %>% 
    mutate(!!current_GO := ifelse((log2fold_change >= 1) & (filtered_VST$hgnc_symbol %in% current_SYMBOL), "up",
                                  ifelse((log2fold_change <= -1) & (filtered_VST$hgnc_symbol %in% current_SYMBOL), "down", NA)))
  }

sample_list <- colnames(VST)[grep("^wi_.*|^h11_.*", colnames(VST))]
sample_order <- unique(sample_list)
sample_list <- factor(sample_list, levels = sample_order)

hm1_data <- filtered_VST %>%
  dplyr::select(matches("^wi_.*|^h11_.*"))
hm1_data <- as.matrix(hm1_data)
row_hm1 <- filtered_VST$x1 %in% rownames(hm1_data)
rownames(hm1_data) <- ifelse(row_hm1, filtered_VST$hgnc_symbol, rownames(hm1_data))
hm1_order <- c(grep("^wi_.*", colnames(hm1_data), value = TRUE),
               grep("^h11_.*", colnames(hm1_data), value = TRUE))


hm2_data <- filtered_VST %>%
  dplyr::select(matches("^GO:"))
row_hm2 <- filtered_VST$x1 %in% rownames(hm2_data)
rownames(hm2_data) <- ifelse(row_hm2, filtered_VST$hgnc_symbol, rownames(hm2_data))
hm2_data <- as.matrix(hm2_data)
hm2_na <- colSums(!is.na(hm2_data))
hm2_order <- order(hm2_na, decreasing = TRUE)

row_order = seriate(dist(hm1_data), method = "TSP")

hm_1 <- Heatmap(hm1_data, name = "VST",
                show_row_dend = FALSE, row_order = get_order(row_order),
                row_names_gp = gpar(fontface = 'bold.italic'), row_names_side = "left",
                column_title = c("1" = "WT", "2" = "NF2 KO"),
                column_title_gp = gpar(fontface = 'bold'),
                column_order = hm1_order, column_km = 2,
                column_names_gp = gpar(fontface = 'bold'),
                col = colorRamp2(c(min(hm1_data), mean(hm1_data), max(hm1_data)), c('steelblue', 'white', 'firebrick')),
                heatmap_width = unit(2.5, "in"),
                border = TRUE)

hm_1
hm_2 <- Heatmap(hm2_data,
                show_row_dend = FALSE, row_order = get_order(row_order), show_row_names = FALSE, 
                column_order = hm2_order,
                column_names_gp = gpar(fontface = 'bold'),
                col = c("up" = "firebrick", "down" = "steelblue"), na_col = "white",
                rect_gp = gpar(col = "white", lty = 1),
                heatmap_legend_param = list(
                  title = "DEG", at = c("up", "down"), labels = c("up" = "Up", "down" = "Down")),
                heatmap_width = unit(1, "in"),
                border = TRUE)
hm_2 
draw(hm_1 + hm_2, auto_adjust = FALSE, ht_gap = unit(c(1, 0), "mm"))

# 
# 
# anno_GO[1,2]
# 
# GO_text <- character()
# for (i in 1:nrow(anno_GO)) {
#   loop_text <- paste(anno_GO[i,1] , anno_GO[i,2])
#   GO_text  %>% rbind(loop_text)
# }
# 
# grid.newpage()
# anno_txt <- textbox_grob(anno_GO, gp = gpar(fontsize = 8))
# grid.draw(anno_txt)

tiff(file = paste0(current_path, "/RNA_seq/Figure/NF2_CM_GO_DEG_HM.tiff"), res = 300, width = 4.3, height = 4.3, units = "in")
draw(hm_1 + hm_2, auto_adjust = FALSE, ht_gap = unit(c(1, 0), "mm"))
dev.off()

