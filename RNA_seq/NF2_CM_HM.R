# Set working environment -------------------------------------------------
rm(list = ls())
setwd("C:/Workspace/R/Project/NF2_CM")

# Load libraries ----------------------------------------------------------
library(readr)
library(readxl)
library(tidyverse)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)
library(ComplexHeatmap)
library(colorRamp2)

# Import normalized count numbers -----------------------------------------
LFC_NF2_iPSC <- janitor::clean_names(na.omit(read_csv("Export/LFC_NF2_CM_iPSC.csv"))) %>% 
  mutate(signif = ifelse(padj < 0.05 & log2fold_change >= 1, "up",
                          ifelse(padj < 0.05 & log2fold_change <= -1, "down", "ns")))
VST_NF2_iPSC <- janitor::clean_names(na.omit(read_csv("Export/VST_NF2_CM_iPSC.csv")))
merged_df <- merge(LFC_NF2_iPSC, VST_NF2_iPSC, by = "hgnc_symbol") 

CM_markers <- janitor::clean_names(na.omit(read_excel("Rawdata/CM_markers.xlsx")))
duplicated_markers <- CM_markers[duplicated(CM_markers) | duplicated(CM_markers, fromLast = TRUE), ]
duplicated_markers

ensembl_ids <- CM_markers %>% 
  mutate(ensembl_id = mapIds(org.Hs.eg.db, keys = CM_markers$target, keytype = "SYMBOL", column = "ENSEMBL"))

filtered_df <- left_join(ensembl_ids, merged_df, by = c("ensembl_id" = "x1.x"))

# Signaling pathways ------------------------------------------------------
GOs <- data.frame(ID = c("0051897", "0051898", "0006915"))
GOs$ID <- gsub("^GO:", "", GOs$ID)
GOs$ID <- paste0("GO:", GOs$ID)      

# Retrieve GO terms and descriptions for all GO IDs -----------------------
GO_list <- lapply(GOs$ID, function(GO_ID) {
  terms <- AnnotationDbi::select(GO.db, keys = GO_ID, keytype = "GOID", columns = c("GOID", "TERM"))
  terms$TERM <- paste0(toupper(substring(terms$TERM, 1, 1)), substring(terms$TERM, 2))
  return(terms)
})
GO_bind <- do.call(rbind, GO_list)
print(GO_bind)

anno_df <- AnnotationDbi::select(org.Hs.eg.db, keys = GO_bind$GOID, keytype = "GOALL", columns=c("ENSEMBL", "SYMBOL"))
joined_df <- anno_df %>% 
  left_join(GO_bind, by = c("GOALL" = "GOID")) %>% 
  left_join(merged_df, by = c("SYMBOL" = "hgnc_symbol")) %>% 
  janitor::clean_names() %>% 
  na.omit

sig_df <- joined_df %>% 
  filter(signif != "ns")

split_df <- split(sig_df, sig_df$goall)
GO_1 <- split_df[[1]][!duplicated(split_df[[1]][[5]]), ]
GO_2 <- split_df[[2]][!duplicated(split_df[[2]][[5]]), ]


bind_df <- rbind(GO_1, GO_2) %>% 
  mutate(go_1 = ifelse(goall == GOs[1, 1] & signif == "up", "up",
                      ifelse(goall == GOs[1, 1] & signif == "down", "down", "ns")),
         go_2 = ifelse(goall == GOs[2, 1] & signif == "up", "up",
                      ifelse(goall == GOs[2, 1] & signif == "down", "down", "ns")))

table(bind_df$go_1)
table(bind_df$go_2)

count_mat <- as.matrix(bind_df[, grepl("i_psc", colnames(bind_df))])
row.names(count_mat) <- bind_df$symbol
count_clabs <- c("Wi-1", "Wi-2", "Wi-3", "H11-1", "H11-2", "H11-3")

GO_mat <- as.matrix(bind_df[, grepl("go_", colnames(bind_df))])
row.names(GO_mat) <- bind_df$symbol
GO_clabs <- c("Positive regulation", "Negative regulation")
GO_clabs <- sapply(GO_clabs, function(name) paste(strwrap(name, width = 20), collapse = "\n"))
GO_col <- c("ns" = "gray95", "down" = "blue", "up" = "red")

hm_1 <- Heatmap(count_mat, name = "VST", col = colorRamp2(c(min(count_mat), mean(count_mat), max(count_mat)),
                                                         c("blue", "white", "red")), border = T,
                width = 5,
                column_dend_height = unit(3, "mm"),
                
                column_title = NULL,
                column_names_gp = gpar(fontface = 'bold', fontsize = 8),
                column_labels = count_clabs,
                column_split = factor(rep(c("Wi", "H11"), each = 3), levels = c( "H11", "Wi")),
                show_row_dend = F,
                show_row_names = F,
                heatmap_legend_param = list(title_position = 'topcenter',
                                            grid_width = unit(3, 'mm'),
                                            title_gp = gpar(fontface = 'bold', fontsize = 5),
                                            labels_gp = gpar(fontface = 'bold', fontsize = 5),
                                            border = "black"))
hm_1

hm_2 <- Heatmap(GO_mat, name = "DEG", border = T, col = GO_col,
                width = 1,
                show_row_names = T,
                row_names_side = 'right',
                row_names_gp = gpar(fontface = 'bold.italic', fontsize = 6),
                column_labels = GO_clabs,
                column_names_centered = TRUE,
                column_names_gp = gpar(fontface = 'bold', fontsize = 8),
                rect_gp = gpar(col = "white", lty = 1),
                heatmap_legend_param = list(title_position = 'topcenter',
                                            grid_width = unit(3, 'mm'),
                                            title_gp = gpar(fontface = 'bold', fontsize = 5),
                                            labels_gp = gpar(fontface = 'bold', fontsize = 5),
                                            at = c('up', 'down', 'ns'),
                                            legend_gp = gpar(fill = 1:3),
                                            border = "black"))
hm_2

# Set the file path for the TIFF image
output_file <- "Figure/HM_iPSC_AKT.tiff"

# Open the TIFF graphics device
tiff(output_file, width = 5, height = 5, units = "in", res = 300)

# Plot the heatmap
draw(hm_1 + hm_2, heatmap_legend_side = "right", ht_gap = unit(1, "mm"), auto_adjust = F)

# Close the graphics device
dev.off()

