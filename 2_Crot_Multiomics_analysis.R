# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        Crotonylation Project                       -----
# -----                                                                    -----
# -----                           Chatterjee Lab                           -----
# -----                         University of Iowa                         -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
# Budhaditya Basu
# 06/07/2024
library(Seurat)
library(Signac)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(scCustomize)
library(paletteer)
library(ggthemes)
library(ggrastr)
library(patchwork)
library(ggrepel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)
library(TFBSTools)
library(UpSetR)
# Significant DARs condition: abs(avg_log2FC) > 0.2 & p.adj < 0.05
# Significant DEGs condition: abs(avg_log2FC) > 0.2 & p.adj < 0.05
# tssRegion=c(-2000, 2000) # Promoter
# TF motif enrichment only in the Promoter region
merged <- readRDS("/home/bbasu/hpchome/R_codes/utsav/integrated_crotonate_celltypes_labeled.rds")
head(merged@meta.data)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 1: Calculate Differentially accessible peaks across all cluster
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DefaultAssay(merged) <- "peaks"
DA <- list()
Idents(merged) <- merged$celltype
for(i in levels(merged)){
  cluster = i
  message("Calculating DA for ",cluster)
  Idents(merged) <- merged$celltype.treatment
  cluster.marker <- FindMarkers(merged, ident.1 = paste0(cluster,"_Crotonate"), 
                                ident.2 = paste0(cluster, "_Saline"), 
                                min.pct=0.05,
                                test.use = 'LR',
                                latent.vars = 'nCount_peaks') # It takes long time!
  DA[[cluster]] <- cluster.marker
}

saveRDS(DA, file = "/home/bbasu/hpchome/R_codes/utsav/Crot_Differential_peaks.rds")
DA <- readRDS("/home/bbasu/hpchome/R_codes/utsav/Crot_Differential_peaks.rds")

#===============================================================================
# Annotate the DARs per cell type
#===============================================================================
DefaultAssay(merged) <- "peaks"
celltype.DAR.list <- list()
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
for (df in names(DA)) {
  tryCatch({
    DARs = DA[[df]]
    message("Doing analysis for ", df)
    # sc-ATAC
    # Run BH procedure
    DARs$p.adj <- p.adjust(DARs$p_val, method = "BH")
    cellType.DARs <- DARs %>%
      dplyr::mutate(celltype = df) %>%
      dplyr::select(-p_val_adj)
    # dplyr::filter(abs(avg_log2FC) > 0.2 & p.adj < 0.05) # Take both upregulated and downregulated peaks
    # find the closest gene to each of these peaks using the ClosestFeature function
    cf <- ClosestFeature(merged, regions = rownames(cellType.DARs)) 
    cellType.DARs <- cbind(cellType.DARs, gene=cf$gene_name, genomicRegion = cf$gene_biotype, 
                           type = cf$type, distance=cf$distance)
    #All Peak Annotation
    All.gr <- StringToGRanges(rownames(cellType.DARs))
    All.peakAnno <- annotatePeak(All.gr, tssRegion=c(-2000, 2000), 
                                 TxDb = txdb, annoDb="org.Mm.eg.db")
    peakAnno.info <- All.peakAnno@anno
    # remove unwanted strings from the annotation column
    peakAnno.info$annotation <- gsub("\\([^()]*\\)", "", peakAnno.info$annotation) 
    peak.info <- data.frame(peaks = rownames(cellType.DARs),
                            peakType = peakAnno.info$annotation, 
                            distanceToTSS = peakAnno.info$distanceToTSS,
                            SYMBOL = peakAnno.info$SYMBOL)
    # Merge peak annotation with the cellType.DARs
    cellType.DARs <- cbind(cellType.DARs, peak.info)
    
    # Store
    celltype.DAR.list[[df]] <- cellType.DARs
  })
}

head(celltype.DAR.list$CA1)
saveRDS(celltype.DAR.list, file = "/home/bbasu/hpchome/R_codes/utsav/Crot_cellType_DARs_annotated.rds")
celltype.DAR.list <- readRDS("/home/bbasu/hpchome/R_codes/utsav/Crot_cellType_DARs_annotated.rds")
# Store in spreadsheet
write.xlsx(celltype.DAR.list, file = "/home/bbasu/hpchome/R_codes/utsav/Differentially_accessible_regions_annotated.xlsx", 
           rowNames = F)

#===============================================================================
# Peak type composition barplot (Upregulated)
#===============================================================================
# Calculate proportion of each peak type which are more accessible (Upregulated)
peak_composition <- data.frame()
num_dars_prom <- list()
for (df in names(celltype.DAR.list)) {
  message("Doing analysis for ", df)
  cellType.DARs = celltype.DAR.list[[df]] %>%
    dplyr::filter(avg_log2FC > 0.2 & p.adj < 0.05) # take positive log2FC values
  # Store number of DARs 
  num_dars_prom[[df]] <- cellType.DARs %>%
    dplyr::filter(str_detect(peakType, "Promoter"))%>%
    distinct(gene)%>%
    nrow()
  # peak_types <- data.frame(prop.table(table(cellType.DARs$peakType)))
  peak_types <- data.frame(table(cellType.DARs$peakType))
  peak_types$celltype <- df
  #Merge
  peak_composition <- rbind(peak_composition, peak_types)
}

peak_composition_filtered <- peak_composition %>%
  dplyr::filter(celltype %in% c("CA1", "CA3", "Dentate Gyrus", "SUB"))

# barplot1 <- ggbarplot(peak_composition_filtered, x = "celltype", y = "Freq", fill = "Var1", #color = "Var1",
#                       palette = alpha(c("plum1", "orange", "skyblue", "palegreen",
#                                         "khaki", "turquoise", "violetred1"), 1), 
#                       xlab = "", 
#                       ylab = "Proportion", label = F)+
#   labs_pubr()+
#   theme(legend.position = "right",
#         legend.title = element_blank())+
#   rotate_x_text(45)

# for (cluster in unique(peak_composition_filtered$celltype)) {
#   cluster_index = which(unique(peak_composition_filtered$celltype) == cluster)
#   barplot1 <- barplot1 + annotate("text", x = cluster_index, y = 0.05,
#                                   label = as.character(num_dars_prom[[cluster]]),
#                                   size = 5)
# }

barplot1 <- ggbarplot(peak_composition_filtered, x = "celltype", y = "Freq", fill = "Var1", #color = "Var1",
                      palette = alpha(c("plum1", "orange", "skyblue", "palegreen",
                                        "khaki", "turquoise", "violetred1"), 1), 
                      xlab = "", 
                      ylab = "Number of peaks", label = F)+
  labs_pubr()+
  theme(legend.position = "right",
        legend.title = element_blank())+
  rotate_x_text(45)



ggsave(filename = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/DARs_annotation_hippocampal_neurons.pdf",
       barplot1,
       height = 6,
       width = 8,
       units = "in",
       dpi = 600)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 2: Calculate Differentially expressed genes across all cluster
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Calculate DEGs across all clusters
DefaultAssay(merged) <- "RNA"
DEG <- list()
Idents(merged) <- merged$celltype
for(i in levels(merged)){
  cluster = i
  message("Calculating DEG for ",cluster)
  Idents(merged) <- merged$celltype.treatment
  cluster.marker <- FindMarkers(merged, ident.1 = paste0(cluster,"_Crotonate"), 
                                ident.2 = paste0(cluster, "_Saline"), 
                                min.pct=0.2, logfc.threshold=0, 
                                test.use = "wilcox") 
  DEG[[cluster]] <- cluster.marker
}
saveRDS(DEG, file = "/home/bbasu/hpchome/R_codes/utsav/DEG_Crotonylation_list_new.rds")
#Import DEGs of all clusters
DEG <- readRDS("/home/bbasu/hpchome/R_codes/utsav/DEG_Crotonylation_list_new.rds")
# Store in spreadsheet
write.xlsx(DEG, file = "/home/bbasu/hpchome/R_codes/utsav/DEG_complete_snRNA.xlsx", 
           rowNames = T)


#===============================================================================
# Upset plot for CA1 CA3 SUB, DG Upregulated only
#===============================================================================
DAR_DEG_overlapped <- list()

for (df in names(celltype.DAR.list)) {
  message("Doing analysis for ", df)
  cellType.DARs = celltype.DAR.list[[df]]
  cellType.DEGs = DEG[[df]]
  # Filter the upregulated(more accessible peaks)
  cellType.DARs <- cellType.DARs %>%
    dplyr::filter(avg_log2FC > 0.2 & p.adj < 0.05)
  # Filter those peaks which are in the promoter region
  DAR_prom <- cellType.DARs %>%
    dplyr::filter(str_detect(peakType, "Promoter"))
  DAR_genebody <- cellType.DARs %>%
    dplyr::filter(!peakType %in% c("Distal Intergenic","Downstream "))
  # scRNA-Seq
  gene.list <- cellType.DEGs %>%
    dplyr::filter(avg_log2FC > 0.20 & p_val_adj < 0.05)%>%
    tibble::rownames_to_column(var = "gene")
  # Overlap 
  overlap_genome <- data.frame("Genomewide" = intersect(unique(cellType.DARs$gene), gene.list$gene)) # For overlapping all accessible peaks
  overlap_prom <- data.frame("Promoter" = intersect(unique(DAR_prom$gene), gene.list$gene)) # For overlapping all accessible promoter peaks
  overlap_genebody <-  data.frame("Genebody" = intersect(unique(DAR_genebody$gene), gene.list$gene)) # For overlapping all accessible promoter and gene body peaks
  overlap_combined <- gdata::cbindX(overlap_prom, 
                                    overlap_genebody, 
                                    overlap_genome)
  # DAR_DEG_overlapped[[df]] <- overlap_combined
  DAR_DEG_overlapped[[df]] <- overlap_genebody$Genebody
}


DAR_DEG_overlapped_filter <- DAR_DEG_overlapped[c("CA1", "CA3", "Dentate Gyrus", "SUB")]

DAR_DEG_overlapped_filter$SUB
write.xlsx(DAR_DEG_overlapped_filter, file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Upregulated_DAR_DEG_overlap_prom_genebody_genomewide.xlsx",
           rowNames = F)



#===============================================================================
# Upset Plot (Upregulated)
#===============================================================================

?upset
pdf("/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/UpregulatedDARgenebody_DEGs_Upset.pdf",
    height = 8, width = 14, onefile = FALSE)
upset_plot <- print(upset(fromList(DAR_DEG_overlapped_filter), 
                          sets.bar.color = c("#2CA02CFF", "#D62728FF", 
                                             "#9467BDFF", "#1F77B4FF"),
                          order.by = "freq",
                          text.scale = 2,
                          mainbar.y.label = "Unique and overlapping genes",
                          sets.x.label = "DAR_DEG overlapped genes",
                          point.size = 4,
                          set_size.show = TRUE))
dev.off()


#===============================================================================
# Perform DARprom/DARgenebody/DARgenomewide and DEG overlapped (Upregulated)
#===============================================================================
DAR_DEG_overlapped <- list()

for (df in names(celltype.DAR.list)) {
  message("Doing analysis for ", df)
  cellType.DARs = celltype.DAR.list[[df]]
  cellType.DEGs = DEG[[df]]
  # Filter the upregulated(more accessible peaks)
  cellType.DARs <- cellType.DARs %>%
    dplyr::filter(avg_log2FC > 0.2 & p.adj < 0.05)
  # Filter those peaks which are in the promoter region
  DAR_prom <- cellType.DARs %>%
    dplyr::filter(str_detect(peakType, "Promoter"))
  DAR_genebody <- cellType.DARs %>%
    dplyr::filter(!peakType %in% c("Distal Intergenic","Downstream "))
  # scRNA-Seq
  gene.list <- cellType.DEGs %>%
    dplyr::filter(avg_log2FC > 0.20 & p_val_adj < 0.05)%>%
    tibble::rownames_to_column(var = "gene")
  # Overlap 
  overlap_genome <- data.frame("Genomewide" = intersect(unique(cellType.DARs$gene), gene.list$gene)) # For overlapping all accessible peaks
  overlap_prom <- data.frame("Promoter" = intersect(unique(DAR_prom$gene), gene.list$gene)) # For overlapping all accessible promoter peaks
  overlap_genebody <-  data.frame("Genebody" = intersect(unique(DAR_genebody$gene), gene.list$gene)) # For overlapping all accessible promoter and gene body peaks
  overlap_combined <- gdata::cbindX(overlap_prom, 
                                    overlap_genebody, 
                                    overlap_genome)
  DAR_DEG_overlapped[[df]] <- overlap_combined
}


DAR_DEG_overlapped_filter <- DAR_DEG_overlapped[c("CA1", "CA3", "Dentate Gyrus", "SUB")]

DAR_DEG_overlapped_filter$CA1


#===============================================================================
# Perform GO enrichment for CA1 CA3 SUB, DG Upregulated only
#===============================================================================

ego_result <- list()
cnet_plot <- list()

for (df in names(DAR_DEG_overlapped_filter)) {
  tryCatch({
    message("Doing analysis for ", df)
    genes = DAR_DEG_overlapped_filter[[df]]$Genebody # Or use Promoter/Genebody
    cellType.DEGs = DEG[[df]]
    # scRNA-Seq
    data <- cellType.DEGs %>%
      dplyr::filter(abs(avg_log2FC) > 0.2 & p_val_adj < 0.05)%>%
      tibble::rownames_to_column(var = "gene_name")
    
    ego <- enrichGO(gene          = genes,
                    OrgDb         = org.Mm.eg.db, # or Org.Hs.eg.db
                    ont           = "MF",
                    #one of “BP”, “MF”, “CC” or “ALL”
                    pAdjustMethod = "BH",
                    #one of “bonferroni”, “BH”, “BY”, “fdr”, “none”
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    keyType = "SYMBOL",
                    #“ENSEMBL”, “ENTREZID”, “SYMBOL”
                    readable      = TRUE)
    ego_result[[df]] <- ego@result
    
    # Make CNET plot
    fc <- data %>%
      dplyr::select(gene_name, avg_log2FC)%>%
      dplyr::filter(gene_name %in% genes)%>%
      deframe()
    
    cnet_plot[[df]] <- cnetplot(ego, circular = TRUE, colorEdge = TRUE, foldChange = fc,
                                showCategory = 5,
                                cex_category = 3.5,
                                cex_gene = 3.5,
                                cex_label_category = 3.5,
                                cex_label_gene = 3.5)+
      paletteer::scale_color_paletteer_c("ggthemes::Red-Blue Diverging", direction = -1,
                                         limits = c(-2,2))+
      theme(legend.text = element_text(size = 16, face= "bold"),
            legend.title = element_text(size = 18, face = "bold"),
            legend.key.size = unit(2, 'cm'), 
            legend.key.height = unit(2, 'cm'), 
            legend.key.width = unit(2, 'cm'))+
      labs(color = "log2FC") 
  }, error=function(e){})
}

write.xlsx(ego_result, file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Upregulated_DARgenebody_DEG_overlapped_MF_enrichment.xlsx", 
           rowNames = F)

#===============================================================================
#  CNET plots for DARgenebody DEG (upregulated)
#===============================================================================
pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Upreg_MF_DARgenebody_CA1_CNET_plot.pdf", 
    height = 26, width = 31)
cnet_plot$CA1
dev.off()

pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Upreg_MF_DARgenebody_CA3_CNET_plot.pdf", 
    height = 26, width = 31)
cnet_plot$CA3
dev.off()

pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Upreg_MF_DARgenebody_SUB_CNET_plot.pdf", 
    height = 26, width = 31)
cnet_plot$SUB
dev.off()

pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Upreg_MF_DARgenebody_DG_CNET_plot.pdf", 
    height = 26, width = 31)
cnet_plot$`Dentate Gyrus`
dev.off()

#===============================================================================
# Upset plot for CA1 CA3 SUB, DG Downregulated only [DARgenebody]
#===============================================================================

DAR_DEG_overlapped <- list()


for (df in names(celltype.DAR.list)) {
  message("Doing analysis for ", df)
  cellType.DARs = celltype.DAR.list[[df]]
  cellType.DEGs = DEG[[df]]
  # Filter the downregulated(less accessible peaks)
  cellType.DARs <- cellType.DARs %>%
    dplyr::filter(avg_log2FC < -0.2 & p.adj < 0.05)
  # Filter those peaks which are in the promoter region
  DAR_prom <- cellType.DARs %>%
    dplyr::filter(str_detect(peakType, "Promoter"))
  # Filter those peaks which are in the promoter plus genebody region
  DAR_genebody <- cellType.DARs %>%
    dplyr::filter(!peakType %in% c("Distal Intergenic","Downstream "))
  # scRNA-Seq
  gene.list <- cellType.DEGs %>%
    dplyr::filter(avg_log2FC < -0.20 & p_val_adj < 0.05)%>%
    tibble::rownames_to_column(var = "gene")
  # Overlap 
  overlap_genome <- data.frame("Genomewide" = intersect(unique(cellType.DARs$gene), gene.list$gene)) # For overlapping all accessible peaks
  overlap_prom <- data.frame("Promoter" = intersect(unique(DAR_prom$gene), gene.list$gene)) # For overlapping all accessible promoter peaks
  overlap_genebody <-  data.frame("Genebody" = intersect(unique(DAR_genebody$gene), gene.list$gene)) # For overlapping all accessible promoter and gene body peaks
  overlap_combined <- gdata::cbindX(overlap_prom, 
                                    overlap_genebody, 
                                    overlap_genome)
  # DAR_DEG_overlapped[[df]] <- overlap_combined
  DAR_DEG_overlapped[[df]] <- overlap_genebody$Genebody
}


DAR_DEG_overlapped_filter <- DAR_DEG_overlapped[c("CA1", "CA3", "Dentate Gyrus", "SUB")]

DAR_DEG_overlapped_filter$SUB
write.xlsx(DAR_DEG_overlapped_filter, file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Downregulated_DAR_DEG_overlap_prom_genebody_genomewide.xlsx",
           rowNames = F)




#===============================================================================
# Upset Plot (Downregulated only)
#===============================================================================

?upset
pdf("/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Downregulated_DARgenebody_DEGs_Upset.pdf",
    height = 8, width = 14, onefile = FALSE)
upset_plot <- print(upset(fromList(DAR_DEG_overlapped_filter), 
                          sets.bar.color = c("#2CA02CFF", "#D62728FF", 
                                             "#9467BDFF", "#1F77B4FF"),
                          order.by = "freq",
                          text.scale = 2,
                          mainbar.y.label = "Unique and overlapping genes",
                          sets.x.label = "DAR_DEG overlapped genes",
                          point.size = 4,
                          set_size.show = TRUE))
dev.off()

#===============================================================================
# Perform GO enrichment for CA1 CA3 SUB, DG Downregulated only
#===============================================================================

DAR_DEG_overlapped <- list()


for (df in names(celltype.DAR.list)) {
  message("Doing analysis for ", df)
  cellType.DARs = celltype.DAR.list[[df]]
  cellType.DEGs = DEG[[df]]
  # Filter the downregulated(less accessible peaks)
  cellType.DARs <- cellType.DARs %>%
    dplyr::filter(avg_log2FC < -0.2 & p.adj < 0.05)
  # Filter those peaks which are in the promoter region
  DAR_prom <- cellType.DARs %>%
    dplyr::filter(str_detect(peakType, "Promoter"))
  # Filter those peaks which are in the promoter plus genebody region
  DAR_genebody <- cellType.DARs %>%
    dplyr::filter(!peakType %in% c("Distal Intergenic","Downstream "))
  # scRNA-Seq
  gene.list <- cellType.DEGs %>%
    dplyr::filter(avg_log2FC < -0.20 & p_val_adj < 0.05)%>%
    tibble::rownames_to_column(var = "gene")
  # Overlap 
  overlap_genome <- data.frame("Genomewide" = intersect(unique(cellType.DARs$gene), gene.list$gene)) # For overlapping all accessible peaks
  overlap_prom <- data.frame("Promoter" = intersect(unique(DAR_prom$gene), gene.list$gene)) # For overlapping all accessible promoter peaks
  overlap_genebody <-  data.frame("Genebody" = intersect(unique(DAR_genebody$gene), gene.list$gene)) # For overlapping all accessible promoter and gene body peaks
  overlap_combined <- gdata::cbindX(overlap_prom, 
                                    overlap_genebody, 
                                    overlap_genome)
  DAR_DEG_overlapped[[df]] <- overlap_combined
}


DAR_DEG_overlapped_filter <- DAR_DEG_overlapped[c("CA1", "CA3", "Dentate Gyrus", "SUB")]




ego_result <- list()
cnet_plot <- list()

for (df in names(DAR_DEG_overlapped_filter)) {
  tryCatch({
    message("Doing analysis for ", df)
    genes = DAR_DEG_overlapped_filter[[df]]$Genebody # Use Promoter/Genebody
    cellType.DEGs = DEG[[df]]
    # scRNA-Seq gene names and fold change values
    data <- cellType.DEGs %>%
      dplyr::filter(abs(avg_log2FC) > 0.2 & p_val_adj < 0.05)%>%
      tibble::rownames_to_column(var = "gene_name")
    
    ego <- enrichGO(gene          = genes,
                    OrgDb         = org.Mm.eg.db, # or Org.Hs.eg.db
                    ont           = "MF",
                    #one of “BP”, “MF”, “CC” or “ALL”
                    pAdjustMethod = "BH",
                    #one of “bonferroni”, “BH”, “BY”, “fdr”, “none”
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    keyType = "SYMBOL",
                    #“ENSEMBL”, “ENTREZID”, “SYMBOL”
                    readable      = TRUE)
    ego_result[[df]] <- ego@result
    
    # Make CNET plot
    fc <- data %>%
      dplyr::select(gene_name, avg_log2FC)%>%
      dplyr::filter(gene_name %in% genes)%>%
      deframe()
    
    cnet_plot[[df]] <- cnetplot(ego, circular = TRUE, colorEdge = TRUE, foldChange = fc,
                                showCategory = 5,
                                cex_category = 3.5,
                                cex_gene = 3.5,
                                cex_label_category = 3.5,
                                cex_label_gene = 3.5)+
      paletteer::scale_color_paletteer_c("ggthemes::Red-Blue Diverging", direction = -1,
                                         limits = c(-2,2))+
      theme(legend.text = element_text(size = 16, face= "bold"),
            legend.title = element_text(size = 18, face = "bold"),
            legend.key.size = unit(2, 'cm'), 
            legend.key.height = unit(2, 'cm'), 
            legend.key.width = unit(2, 'cm'))+
      labs(color = "log2FC") 
  }, error=function(e){})
}


write.xlsx(ego_result, file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Downregulated_DARgenebody_DEG_overlapped_MF_enrichment.xlsx", 
           rowNames = F)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  CNET plots for DARgenebody DEG downregulated
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Downreg_MF_DARgenebody_CA1_CNET_plot.pdf", 
#     height = 26, width = 31)
# cnet_plot$CA1
# dev.off()

pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Downreg_MF_DARgenebody_CA3_CNET_plot.pdf", 
    height = 26, width = 31)
cnet_plot$CA3
dev.off()

# pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Downreg_MF_DARgenebody_SUB_CNET_plot.pdf", 
#     height = 26, width = 31)
# cnet_plot$SUB
# dev.off()

pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Downreg_MF_DARgenebody_DG_CNET_plot.pdf", 
    height = 26, width = 31)
cnet_plot$`Dentate Gyrus`
dev.off()


#===============================================================================
# Compare DA and DE by scatter plot (Scatter plot DARgenebody)
#===============================================================================
lfc.thresh <- 0.2
padj.thresh <- 0.05
scatter.plot <- list()
DEG_DAR_merge <- list()

for (df in names(DEG)){
  tryCatch({
    message("Doing analysis for ", df)
    dars = celltype.DAR.list[[df]] %>%
      dplyr::filter(!peakType %in% c("Distal Intergenic","Downstream "))%>%
      dplyr::filter(abs(avg_log2FC) > 0.2 & p.adj < 0.05)%>%
      dplyr::select(gene, avg_log2FC, p.adj, peaks, celltype, distanceToTSS, peakType)%>%
      dplyr::rename(ATAC_LFC = avg_log2FC,
                    ATAC_BH = p.adj)
    degs = DEG[[df]] %>%
      tibble::rownames_to_column(var = "gene")%>%
      dplyr::filter(abs(avg_log2FC) > 0.2 & p_val_adj < 0.05)%>%
      dplyr::select(gene, avg_log2FC, p_val_adj)%>%
      dplyr::rename(RNA_LFC = avg_log2FC,
                    RNA_BH = p_val_adj)
    
    # Merge together
    res <- merge(degs, dars)
    # Create signed log10(FDR) for visualization.
    res$atac.slfdr = -log10(res$ATAC_BH)*sign(res$ATAC_LFC)
    res$rna.slfdr = -log10(res$RNA_BH)*sign(res$RNA_LFC)
    # Add color column
    res$color <- rep("notSig", nrow(res))
    res[which((abs(res$RNA_LFC) > lfc.thresh & res$RNA_BH < padj.thresh) | (abs(res$ATAC_LFC) > lfc.thresh & res$ATAC_BH < padj.thresh)), "color"] <- "sigInOnlyOne"
    res[which(res$RNA_LFC > lfc.thresh & res$RNA_BH < padj.thresh & res$ATAC_LFC > lfc.thresh & res$ATAC_BH < padj.thresh), "color"] <- "UpATAC+UpRNA"
    res[which(res$RNA_LFC < -lfc.thresh & res$RNA_BH < padj.thresh & res$ATAC_LFC < -lfc.thresh & res$ATAC_BH < padj.thresh), "color"] <- "DownATAC+DownRNA"
    res[which(res$RNA_LFC > lfc.thresh & res$RNA_BH < padj.thresh & res$ATAC_LFC < -lfc.thresh & res$ATAC_BH < padj.thresh), "color"] <- "Discordant.change"
    res[which(res$RNA_LFC < -lfc.thresh & res$RNA_BH < padj.thresh & res$ATAC_LFC > lfc.thresh & res$ATAC_BH < padj.thresh), "color"] <- "Discordant.change"
    # Calculate limits for square visual
    # lim <- max(max(abs(res$rna.slfdr)), max(abs(res$atac.slfdr)))
    
    # Calculate the number of concordant increase and decreased unique genes
    concord_up <- res %>%
      dplyr::filter(color == "UpATAC+UpRNA")%>%
      distinct(gene, .keep_all = TRUE)%>%
      nrow()
    concord_down <- res %>%
      dplyr::filter(color == "DownATAC+DownRNA")%>%
      distinct(gene, .keep_all = TRUE)%>%
      nrow()
    
    cor_test <- cor.test(res$rna.slfdr, res$atac.slfdr, method = "pearson")
    # Plot data
    p <- ggplot(res,
                aes(x = atac.slfdr,y = rna.slfdr,fill = color,
                    label = gene,
                    size = abs(ATAC_LFC*RNA_LFC))) +
      theme_base() +
      geom_point(alpha = 0.9, shape = 21, stroke = 0.1, colour = "black")+
      geom_vline(aes(xintercept = 0)) +
      geom_hline(aes(yintercept = 0)) +
      geom_vline(xintercept = 0.83,linetype = 2,color = "gray") +
      geom_vline(xintercept = -0.83,linetype = 2,color = "gray") +
      geom_hline(yintercept = 0.83,linetype = 2,color = "gray") +
      geom_hline(yintercept = -0.83,linetype = 2,color = "gray") +
      scale_fill_manual(values = c("notSig" = "gray90",
                                   "sigInOnlyOne" = "gray30",
                                   "UpATAC+UpRNA" = "red",
                                   "DownATAC+DownRNA" = "blue",
                                   "Discordant.change"= "gray90"),
                        labels = c("Discordant change",
                                   "Concordant decrease",
                                   "Concordant Increase",
                                   "Significant in one assay",
                                   "Not Significant",
                                   "Concordant Increase"))+
      # scale_x_continuous(limits = c(-lim, lim)) +
      # scale_y_continuous(limits = c(-lim, lim))+
      ggtitle(paste0("DARs intersecting DEGs for ", df))+
      labs(x = expression("Signed -log"[10]*"(FDR) :  ATAC"), 
           y = expression("Signed -log"[10]*"(FDR) :  RNA"), 
           fill = "Direction of LFC")
    
    # if(df %in% c("SUB")){
    #   # Add ggplot2 layers
    #   p1 <- p +
    #     annotate("text",x=20,y=-40, label=paste0( "R: ",round(cor_test$estimate,4)))+
    #     annotate("text",x=20,y=-50, label=paste0( "P Value: ",round(cor_test$p.value,4)))+
    #     annotate("text", x = 30, y = 30,
    #              label = as.character(paste0("Concordant \nIncrease: ", concord_up)),
    #              size = 4)+
    #     annotate("text", x = -10, y = -20,
    #              label = as.character(paste0("Concordant \nDecrease: ", concord_down)),
    #              size = 4)
    # } 
    # else{
    #   p1 <- p +
    #     annotate("text",x=30,y=-40, label=paste0( "R: ",round(cor_test$estimate,4)))+
    #     annotate("text",x=30,y=-50, label=paste0( "P Value: ",round(cor_test$p.value,4)))+
    #     annotate("text", x = 50, y = 100,
    #              label = as.character(paste0("Concordant \nIncrease: ", concord_up)),
    #              size = 4)+
    #     annotate("text", x = -30, y = -100,
    #              label = as.character(paste0("Concordant \nDecrease: ", concord_down)),
    #              size = 4)
    # }
    # +geom_label_repel(
    #   data = res[which(res$color %in% c("UpATAC+UpRNA", "DownATAC+DownRNA")), ],
    #   inherit.aes = T,
    #   fill = 'white',
    #   color = 'black',
    #   max.overlaps = 20,
    #   size = 2,
    #   force = 5)
    DEG_DAR_merge[[df]] <- res
    scatter.plot[[df]] <- p
  }, error=function(e){})
}


saveRDS(DEG_DAR_merge, file = "/home/bbasu/hpchome/R_codes/utsav/DAR_DEGs_Overlapped_complete.rds")
DEG_DAR_merge <- readRDS(file = "/home/bbasu/hpchome/R_codes/utsav/DAR_DEGs_Overlapped_complete.rds")
# Store in spreadsheet
write.xlsx(DEG_DAR_merge, file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/DAR_DEGs_Overlapped_complete.xlsx", 
           rowNames = F)
png(filename = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Overlap_scatter_corr.png",
    height = 12, width = 20, units = "in", res = 600)
(scatter.plot$CA1+scatter.plot$CA3)/(scatter.plot$SUB+scatter.plot$`Dentate Gyrus`)
dev.off()

tiff(filename = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Overlap_scatter_corr.tiff",
     height = 12, width = 20, units = "in", res = 600)
(scatter.plot$CA1+scatter.plot$CA3)/(scatter.plot$SUB+scatter.plot$`Dentate Gyrus`)
dev.off()

# Save as pdf
pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Overlap_scatter_corr.pdf",
    height = 12, width = 18)
(scatter.plot$CA1+scatter.plot$CA3)/(scatter.plot$SUB+scatter.plot$`Dentate Gyrus`)
dev.off()


#===============================================================================
# Volcano plot for all clusters (DARgenebody)
#===============================================================================
DEG_DAR_merge <- readRDS(file = "/home/bbasu/hpchome/R_codes/utsav/DAR_DEGs_Overlapped_complete.rds")
colors <- paletteer::paletteer_d("ggsci::category20_d3")%>%
  head(12)
plots <- list()
for(i in names(DEG)){
  cluster = i
  message("Plotting for ",cluster)
  cluster_index = which(names(DEG) == cluster)
  #Data wrangling
  cluster.marker <- DEG[[cluster]]
  res <- DEG_DAR_merge[[cluster]]
  # Calculate the number of concordant increase and decreased unique genes
  concord_up <- res %>%
    dplyr::filter(color == "UpATAC+UpRNA")%>%
    distinct(gene)
  concord_down <- res %>%
    dplyr::filter(color == "DownATAC+DownRNA")%>%
    distinct(gene)
  
  data <- data.frame(gene = row.names(cluster.marker),
                     pval = -log10(cluster.marker$p_val_adj), 
                     lfc = cluster.marker$avg_log2FC)
  
  data <- mutate(data, color = case_when(data$gene %in% concord_up$gene ~ "Increased",
                                         data$gene %in% concord_down$gene ~ "Decreased", 
                                         !(data$gene %in% c(concord_up$gene, concord_down$gene))~"nonsignificant"))
  
  # Make a basic ggplot2 object with x-y values
  vol <- ggplot(data, aes(x = lfc, y = pval, color = color))+
    geom_point(data = subset(data, data$color %in% c("nonsignificant")),
               size = 1, na.rm = T) +
    geom_point(data = subset(data, data$color %in% c("Increased", "Decreased")),
               size = 1, na.rm = T) +
    scale_color_manual(name = "Directionality",
                       values = c(Increased = colors[cluster_index], 
                                  Decreased = colors[cluster_index], 
                                  nonsignificant = "gray90")) +
    theme_base() + # change overall theme
    theme(legend.position = "none") + # change the legend
    xlim(-2,2)+
    scale_y_continuous(breaks = seq(0, 300, by=150), limits=c(0,300))+
    theme(plot.title = element_text(size = 10))+
    ggtitle(cluster)+
    annotate("text", x = -1.4, y = 250, label = capture.output(writeLines(paste0("Up :", nrow(concord_up)))), 
             size = 4,
             fontface =1.2)+
    annotate("text", x = -1.4, y = 230, label = capture.output(writeLines(paste0("Down :", nrow(concord_down)))), 
             size = 4,
             fontface =1.2)
  
  if(cluster %in% c("CA1", "Inhibitory neurons")){
    # Add ggplot2 layers
    p1 <- vol+
      xlab("")+
      ylab(expression(-log[10]~"(Adj P Value)"))
  } 
  else if(cluster %in% c("APC", "Cajal-Retzius", "Excitatory cortical neurons")){
    # Add ggplot2 layers
    p1 <- vol +
      xlab(expression(log[2]~"(Fold Change)"))+
      ylab("")
  } 
  else if(cluster == "OPC"){
    # Add ggplot2 layers
    p1 <- vol +
      xlab(expression(log[2]~"(Fold Change)"))+
      ylab(expression(-log[10]~"(Adj P Value)"))
  }
  else{
    # Add ggplot2 layers
    p1 <- vol +
      xlab("")+
      ylab("")
    
  }
  plots[[cluster]] <- p1
}
#saveRDS(plots, file = "/home/bbasu/hpchome/R_codes/utsav/Plots_Crotonylation_list.rds")

pdf("/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/Volcano_Crotonylation.pdf", height = 10, width = 12)
plots$CA1+plots$CA3+plots$SUB+plots$`Dentate Gyrus`+plots$`Inhibitory neurons`+plots$Oligodendrocytes+plots$Astrocytes+plots$Microglia+plots$OPC+plots$`Excitatory cortical neurons`+plots$APC+plots$`Cajal-Retzius`
dev.off()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step 3: TF motif enrichment analysis only in the promoter region only from ATAC peaks
# [tssRegion=c(-2000, 2000) # Promoter] 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#===============================================================================
# Compute Motif enrichment analysis
#===============================================================================

DefaultAssay(merged) <- "peaks"

# first compute the GC content for each peak
merged <- RegionStats(merged, genome = BSgenome.Mmusculus.UCSC.mm10)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

# add motif information
merged <- AddMotifs(
  object = merged,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm) # Takes 5 mins

# Motif Enrichment only in the promoter region
#===============================================================================
DEG_DAR_merge <- list()
padj.thresh <- 0.05
lfc.thresh <- 0.2

for (df in names(DEG)){
  tryCatch({
    message("Doing analysis for ", df)
    dars = celltype.DAR.list[[df]] %>%
      dplyr::filter(str_detect(peakType, "Promoter"))%>%
      # dplyr::filter(abs(avg_log2FC) > 0.2 & p.adj < 0.05)%>%
      dplyr::select(gene, avg_log2FC, p.adj, peaks, celltype, distanceToTSS, peakType)%>%
      dplyr::rename(ATAC_LFC = avg_log2FC,
                    ATAC_BH = p.adj)
    degs = DEG[[df]] %>%
      tibble::rownames_to_column(var = "gene")%>%
      # dplyr::filter(abs(avg_log2FC) > 0.2 & p_val_adj < 0.05)%>%
      dplyr::select(gene, avg_log2FC, p_val_adj)%>%
      dplyr::rename(RNA_LFC = avg_log2FC,
                    RNA_BH = p_val_adj)
    
    # Merge together
    res <- merge(degs, dars)
    # Create signed log10(FDR) for visualization.
    res$atac.slfdr = -log10(res$ATAC_BH)*sign(res$ATAC_LFC)
    res$rna.slfdr = -log10(res$RNA_BH)*sign(res$RNA_LFC)
    # Add color column
    res$color <- rep("notSig", nrow(res))
    res[which((abs(res$RNA_LFC) > lfc.thresh & res$RNA_BH < padj.thresh) | (abs(res$ATAC_LFC) > lfc.thresh & res$ATAC_BH < padj.thresh)), "color"] <- "sigInOnlyOne"
    res[which(res$RNA_LFC > lfc.thresh & res$RNA_BH < padj.thresh & res$ATAC_LFC > lfc.thresh & res$ATAC_BH < padj.thresh), "color"] <- "UpATAC+UpRNA"
    res[which(res$RNA_LFC < -lfc.thresh & res$RNA_BH < padj.thresh & res$ATAC_LFC < -lfc.thresh & res$ATAC_BH < padj.thresh), "color"] <- "DownATAC+DownRNA"
    res[which(res$RNA_LFC > lfc.thresh & res$RNA_BH < padj.thresh & res$ATAC_LFC < -lfc.thresh & res$ATAC_BH < padj.thresh), "color"] <- "Discordant.change"
    res[which(res$RNA_LFC < -lfc.thresh & res$RNA_BH < padj.thresh & res$ATAC_LFC > lfc.thresh & res$ATAC_BH < padj.thresh), "color"] <- "Discordant.change"
    
    # Store the data
    DEG_DAR_merge[[df]] <- res
    
  })
}


enriched.motif <- list()

for (cluster in c("CA1", "CA3", "SUB", "Dentate Gyrus")) {
  tryCatch({
    message("Performing motif enrichment for ",cluster)
    
    res <- DEG_DAR_merge[[cluster]]
    
    shared.up.peaks <- res %>%
      dplyr::filter(color == "UpATAC+UpRNA")
    # shared.down.peaks <- res %>%
    #   dplyr::filter(color == "DownATAC+DownRNA")
    
    # Enriched in Up DARs+DEGs
    up_enriched_motifs <- FindMotifs(
      object = merged,
      features = shared.up.peaks$peaks,
      assay = "peaks")
    
    up_enriched_motifs$LFC <- log2(up_enriched_motifs$fold.enrichment)
    up_enriched_motifs <- up_enriched_motifs[,c("motif.name", "LFC", "p.adjust")]
    enriched_motifs <- up_enriched_motifs %>%
      dplyr::filter(p.adjust < 0.05)
    enriched.motif[[cluster]] <- enriched_motifs
    # Enriched in Down DARs+DEGs
    # down_enriched_motifs <- FindMotifs(
    #   object = merged,
    #   features = shared.down.peaks$peaks,
    #   assay = "peaks")
    # 
    # down_enriched_motifs$LFC <- log2(down_enriched_motifs$fold.enrichment)
    # down_enriched_motifs <- down_enriched_motifs[,c("motif.name",  "LFC", "p.adjust")]
    # 
    # # Merge motifs
    # enriched_motifs <- rbind(up_enriched_motifs, down_enriched_motifs)
    # enriched_motifs <- enriched_motifs %>%
    #   dplyr::filter(p.adjust < 0.05)
    
    # enriched.motif[[cluster]] <- enriched_motifs
  }, error=function(e){})
}

# Save enriched motifs
saveRDS(enriched.motif, file = "/home/bbasu/hpchome/R_codes/utsav/TF.enriched.motifs_list_DARprom_DEG_overlapped.rds")

enriched.motif <- readRDS("/home/bbasu/hpchome/R_codes/utsav/TF.enriched.motifs_list_DARprom_DEG_overlapped.rds")
# Store in spreadsheet
write.xlsx(enriched.motif, file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/TF.enriched.motifs_DARprom_DEG_overlapped.xlsx", 
           rowNames = F)


# Motif Plot
#===============================================================================
# Select top15 CA1 motifs
ca1_top15_motifs <- enriched.motif$CA1 %>%
  dplyr::arrange(p.adjust)%>%
  dplyr::select(motif.name)%>%
  head(15)%>% deframe()

# Select top15 CA3 motifs
ca3_top15_motifs <- enriched.motif$CA3 %>%
  dplyr::arrange(p.adjust)%>%
  dplyr::select(motif.name)%>%
  head(15)%>% deframe()

# Select top15 SUB motifs
sub_top15_motifs <- enriched.motif$SUB %>%
  dplyr::arrange(p.adjust)%>%
  dplyr::select(motif.name)%>%
  head(15)%>% deframe()


pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/CA1_promoter_top15_motifs_DARprom_DEG_overlapped.pdf",
    height = 10, width = 12)
MotifPlot(
  object = merged,
  motifs = ca1_top15_motifs)
dev.off()

pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/CA3_promoter_top15_motifs_DARprom_DEG_overlapped.pdf",
    height = 10, width = 12)
MotifPlot(
  object = merged,
  motifs = ca3_top15_motifs)
dev.off()

pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/SUB_promoter_top15_motifs_DARprom_DEG_overlapped.pdf",
    height = 10, width = 12)
MotifPlot(
  object = merged,
  motifs = sub_top15_motifs)
dev.off()

#===============================================================================
# Heatmap of top15 transcription factors enrichment
#===============================================================================
library(purrr)
library(ComplexHeatmap)
library(tidyverse)
enriched.motif <- readRDS("/home/bbasu/hpchome/R_codes/utsav/TF.enriched.motifs_list_DARprom_DEG_overlapped.rds")

#===============================================================================
# Plot only p Value
ca1_top15_motifs <- enriched.motif$CA1 %>%
  dplyr::arrange(p.adjust)%>%
  dplyr::mutate(ca1_es = -log10(p.adjust))%>%
  dplyr::select(motif.name, ca1_es)%>%
  head(15)

ca3_top15_motifs <- enriched.motif$CA3 %>%
  dplyr::arrange(p.adjust)%>%
  dplyr::mutate(ca3_es = -log10(p.adjust))%>%
  dplyr::select(motif.name, ca3_es)%>%
  head(15)

sub_top15_motifs <- enriched.motif$SUB %>%
  dplyr::arrange(p.adjust)%>%
  dplyr::mutate(sub_es = -log10(p.adjust))%>%
  dplyr::select(motif.name, sub_es)%>%
  head(15)

#===============================================================================

motifs_list <- list(ca1_top15_motifs, ca3_top15_motifs, sub_top15_motifs)

merged_motifs <- purrr::reduce(motifs_list, full_join, by = "motif.name")

merged_motifs[is.na(merged_motifs)] <- 0

mat <- merged_motifs %>%
  column_to_rownames(var = "motif.name")



mat <- as.matrix(mat[,c(1:3)])
# Scale the expression value
# scaled_mat <- t(scale(t(mat)))

my_palette <- paletteer::paletteer_c("grDevices::Blues 3", 30, direction = -1)
# my_palette <- paletteer::paletteer_c("grDevices::Blue-Red 2", 30)
# my_palette <- paletteer::paletteer_c("grDevices::Tropic", 30)

ha = HeatmapAnnotation(df = data.frame(Celltypes = c("CA1","CA3",
                                                     "SUB")),
                       col = list(Celltypes = c("CA1" = "#6FCEB4",
                                                "CA3" = "#DC0085",
                                                "SUB" = "#FFBB78")),
                       annotation_legend_param = list(title = "TF Motif Enrichment",
                                                      title_gp = gpar(fontsize = 10,
                                                                      fontface = "bold")),
                       annotation_name_side = "left",
                       simple_anno_size = unit(2, "mm"))

ht1 <- Heatmap(mat,cluster_columns = F, cluster_rows = T,
               width = unit(3, "cm"), 
               height = unit(15, "cm"),
               show_column_names = F,
               show_row_names = T,
               col = my_palette,
               bottom_annotation = ha,
               name = "-log10(P Value)",
               heatmap_legend_param = list(
                 title_gp = gpar(fontsize = 10,
                                 fontface = "bold")
               ))

ht1


pdf(file = "/home/bbasu/hpchome/R_codes/utsav/Multiomics_crotonylation_data/TF_motif_enrichment_heatmap_DARprom_DEG_overlapped.pdf",
    height = 8, width = 8)
ht1
dev.off()

