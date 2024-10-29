# ------------------------------------------------------------------------------
# -----                                                                    -----
# -----                        Crotonylation Project                       -----
# -----                                                                    -----
# -----                           Chatterjee Lab                           -----
# -----                         University of Iowa                         -----
# -----                                                                    -----
# ------------------------------------------------------------------------------
# Budhaditya Basu
# 06/03/2024
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(tidyverse)
library(patchwork)
library(scCustomize)
library(ggpubr)

#############################SALINE 1 ##########################################
#Load the RNA and ATAC data
counts_saline1 <- Read10X_h5("/home/bbasu/LSS/lss_schatterj/Multiomic_crotonylation/Saline1/outs/filtered_feature_bc_matrix.h5")
metadata_saline1 <- read.csv('/home/bbasu/LSS/lss_schatterj/Multiomic_crotonylation/Saline1/outs/per_barcode_metrics.csv', 
                             header = TRUE, row.names = 1, stringsAsFactors = FALSE)
fragpath_saline1 <- "/home/bbasu/LSS/lss_schatterj/Multiomic_crotonylation/Saline1/outs/atac_fragments.tsv.gz"

#Create a Seurat object containing the RNA data
saline1 <- CreateSeuratObject(
  counts = counts_saline1$`Gene Expression`,
  assay = "RNA",
  project = "Saline1",
  meta.data = metadata_saline1
)
head(saline1@meta.data)
grep ("^mt-", rownames(saline1[["RNA"]]), value = T)

saline1[["percent.mt"]] <- PercentageFeatureSet(saline1, pattern = "^mt-")
## Now add in the ATAC-seq data
# Only use peaks in standard chromosomes
atac_counts <- counts_saline1$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Get gene annotations for mm10
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

saline1[['ATAC']] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = fragpath_saline1,
  min.cells = 10,
  annotation = annotations
)

DefaultAssay(saline1) <- 'ATAC'
saline1 <- NucleosomeSignal(saline1)
saline1 <- TSSEnrichment(saline1)
saline1$pct_reads_in_peaks <- saline1$atac_peak_region_fragments / saline1$atac_fragments * 100
saline1$blacklist_fraction <- FractionCountsInRegion(saline1, assay = 'ATAC', regions = blacklist_mm10)


## Peak calling using MACS2
peaks <- CallPeaks(saline1, 
                   assay = 'ATAC',
                   macs2.path = "/home/bbasu/.local/bin/macs2")
peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
peaks <- subsetByOverlaps(x = peaks,
                          ranges = blacklist_mm10, 
                          invert = TRUE)
# quantify counts in each peak
macs_count_saline1 <- FeatureMatrix(fragments = Fragments(saline1),
                                    features = peaks,
                                    cells = colnames(saline1))

saline1[['peaks']] <- CreateChromatinAssay(
  counts = macs_count_saline1,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = fragpath_saline1,
  min.cells = 10,
  annotation = annotations
)
## Save the sample-specific Seurat object 
saveRDS(saline1, file = "/home/bbasu/hpchome/R_codes/utsav/saline1.rds")

#############################SALINE 2 ##########################################
#Load the RNA and ATAC data
counts_saline2 <- Read10X_h5("/home/bbasu/LSS/lss_schatterj/Multiomic_crotonylation/Saline2/outs/filtered_feature_bc_matrix.h5")
metadata_saline2 <- read.csv('/home/bbasu/LSS/lss_schatterj/Multiomic_crotonylation/Saline2/outs/per_barcode_metrics.csv', 
                             header = TRUE, row.names = 1, stringsAsFactors = FALSE)
fragpath_saline2 <- "/home/bbasu/LSS/lss_schatterj/Multiomic_crotonylation/Saline2/outs/atac_fragments.tsv.gz"

#Create a Seurat object containing the RNA data
saline2 <- CreateSeuratObject(
  counts = counts_saline2$`Gene Expression`,
  assay = "RNA",
  project = "Saline2",
  meta.data = metadata_saline2
)
head(saline2@meta.data)
grep ("^mt-", rownames(saline2[["RNA"]]), value = T)

saline2[["percent.mt"]] <- PercentageFeatureSet(saline2, pattern = "^mt-")
## Now add in the ATAC-seq data
# Only use peaks in standard chromosomes
atac_counts <- counts_saline2$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Get gene annotations for mm10
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "mm10"

saline2[['ATAC']] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = fragpath_saline2,
  min.cells = 10,
  annotation = annotations
)

DefaultAssay(saline2) <- 'ATAC'
saline2 <- NucleosomeSignal(saline2)
saline2 <- TSSEnrichment(saline2)
saline2$pct_reads_in_peaks <- saline2$atac_peak_region_fragments / saline2$atac_fragments * 100
saline2$blacklist_fraction <- FractionCountsInRegion(saline2, assay = 'ATAC', regions = blacklist_mm10)

## Peak calling using MACS2
peaks <- CallPeaks(saline2, 
                   assay = 'ATAC',
                   macs2.path = "/home/bbasu/.local/bin/macs2")
peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
peaks <- subsetByOverlaps(x = peaks,
                          ranges = blacklist_mm10, 
                          invert = TRUE)
# quantify counts in each peak
macs_count_saline2 <- FeatureMatrix(fragments = Fragments(saline2),
                                    features = peaks,
                                    cells = colnames(saline2))

saline2[['peaks']] <- CreateChromatinAssay(
  counts = macs_count_saline2,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = fragpath_saline2,
  min.cells = 10,
  annotation = annotations
)
## Save the sample-specific Seurat object 
saveRDS(saline2, file = "/home/bbasu/hpchome/R_codes/utsav/saline2.rds")


#############################CROTONATE 1 #######################################
#Load the RNA and ATAC data
counts_crotonate1 <- Read10X_h5("/home/bbasu/LSS/lss_schatterj/Multiomic_crotonylation/Crotonate1/outs/filtered_feature_bc_matrix.h5")
metadata_crotonate1 <- read.csv('/home/bbasu/LSS/lss_schatterj/Multiomic_crotonylation/Crotonate1/outs/per_barcode_metrics.csv', 
                                header = TRUE, row.names = 1, stringsAsFactors = FALSE)
fragpath_crotonate1 <- "/home/bbasu/LSS/lss_schatterj/Multiomic_crotonylation/Crotonate1/outs/atac_fragments.tsv.gz"

#Create a Seurat object containing the RNA data
crotonate1 <- CreateSeuratObject(
  counts = counts_crotonate1$`Gene Expression`,
  assay = "RNA",
  project = "Crotonate1",
  meta.data = metadata_crotonate1
)
head(crotonate1@meta.data)
grep ("^mt-", rownames(crotonate1[["RNA"]]), value = T)

crotonate1[["percent.mt"]] <- PercentageFeatureSet(crotonate1, pattern = "^mt-")
## Now add in the ATAC-seq data
# Only use peaks in standard chromosomes
atac_counts <- counts_crotonate1$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Get gene annotations for mm10
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

crotonate1[['ATAC']] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = fragpath_crotonate1,
  min.cells = 10,
  annotation = annotations
)

DefaultAssay(crotonate1) <- 'ATAC'
crotonate1 <- NucleosomeSignal(crotonate1)
crotonate1 <- TSSEnrichment(crotonate1)
crotonate1$pct_reads_in_peaks <- crotonate1$atac_peak_region_fragments / crotonate1$atac_fragments * 100
crotonate1$blacklist_fraction <- FractionCountsInRegion(crotonate1, assay = 'ATAC', regions = blacklist_mm10)


## Peak calling using MACS2
peaks <- CallPeaks(crotonate1, 
                   assay = 'ATAC',
                   macs2.path = "/home/bbasu/.local/bin/macs2")
peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
peaks <- subsetByOverlaps(x = peaks,
                          ranges = blacklist_mm10, 
                          invert = TRUE)
# quantify counts in each peak
macs_count_crotonate1 <- FeatureMatrix(fragments = Fragments(crotonate1),
                                       features = peaks,
                                       cells = colnames(crotonate1))

crotonate1[['peaks']] <- CreateChromatinAssay(
  counts = macs_count_crotonate1,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = fragpath_crotonate1,
  min.cells = 10,
  annotation = annotations
)
## Save the sample-specific Seurat object 
saveRDS(crotonate1, file = "/home/bbasu/hpchome/R_codes/utsav/crotonate1.rds")

#############################CROTONATE 2 #######################################
#Load the RNA and ATAC data
counts_crotonate2 <- Read10X_h5("/home/bbasu/LSS/lss_schatterj/Multiomic_crotonylation/Crotonate2/outs/filtered_feature_bc_matrix.h5")
metadata_crotonate2 <- read.csv('/home/bbasu/LSS/lss_schatterj/Multiomic_crotonylation/Crotonate2/outs/per_barcode_metrics.csv', 
                                header = TRUE, row.names = 1, stringsAsFactors = FALSE)
fragpath_crotonate2 <- "/home/bbasu/LSS/lss_schatterj/Multiomic_crotonylation/Crotonate2/outs/atac_fragments.tsv.gz"

#Create a Seurat object containing the RNA data
crotonate2 <- CreateSeuratObject(
  counts = counts_crotonate2$`Gene Expression`,
  assay = "RNA",
  project = "Crotonate2",
  meta.data = metadata_crotonate2
)
head(crotonate2@meta.data)
grep ("^mt-", rownames(crotonate2[["RNA"]]), value = T)

crotonate2[["percent.mt"]] <- PercentageFeatureSet(crotonate2, pattern = "^mt-")
## Now add in the ATAC-seq data
# Only use peaks in standard chromosomes
atac_counts <- counts_crotonate2$Peaks
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Get gene annotations for mm10
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "mm10"

crotonate2[['ATAC']] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = fragpath_crotonate2,
  min.cells = 10,
  annotation = annotations
)

DefaultAssay(crotonate2) <- 'ATAC'
crotonate2 <- NucleosomeSignal(crotonate2)
crotonate2 <- TSSEnrichment(crotonate2)
crotonate2$pct_reads_in_peaks <- crotonate2$atac_peak_region_fragments / crotonate2$atac_fragments * 100
crotonate2$blacklist_fraction <- FractionCountsInRegion(crotonate2, assay = 'ATAC', regions = blacklist_mm10)


## Peak calling using MACS2
peaks <- CallPeaks(crotonate2, 
                   assay = 'ATAC',
                   macs2.path = "/home/bbasu/.local/bin/macs2")
peaks <- keepStandardChromosomes(peaks, pruning.mode = 'coarse')
peaks <- subsetByOverlaps(x = peaks,
                          ranges = blacklist_mm10, 
                          invert = TRUE)
# quantify counts in each peak
macs_count_crotonate2 <- FeatureMatrix(fragments = Fragments(crotonate2),
                                       features = peaks,
                                       cells = colnames(crotonate2))

crotonate2[['peaks']] <- CreateChromatinAssay(
  counts = macs_count_crotonate2,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = fragpath_crotonate2,
  min.cells = 10,
  annotation = annotations
)
## Save the sample-specific Seurat object 
saveRDS(crotonate2, file = "/home/bbasu/hpchome/R_codes/utsav/crotonate2.rds")

##################################Integration ##################################
#Import RDS object
library(Seurat)
library(Signac)
library(tidyverse)
library(ggpubr)

saline1 <- readRDS("/home/bbasu/hpchome/R_codes/utsav/saline1.rds")
saline2 <- readRDS("/home/bbasu/hpchome/R_codes/utsav/saline2.rds")
crotonate1 <- readRDS("/home/bbasu/hpchome/R_codes/utsav/crotonate1.rds")
crotonate2 <- readRDS("/home/bbasu/hpchome/R_codes/utsav/crotonate2.rds")

#Assign condition before integration
saline1$treatment <- "Saline"
saline2$treatment <- "Saline"
crotonate1$treatment <- "Crotonate"
crotonate2$treatment <- "Crotonate"

head(saline1@meta.data)
head(saline2@meta.data)
head(crotonate1@meta.data)
head(crotonate2@meta.data)

DefaultAssay(saline1) <- 'RNA'
DefaultAssay(saline2) <- 'RNA'
DefaultAssay(crotonate1) <- 'RNA'
DefaultAssay(crotonate2) <- 'RNA'

# Filter out low quality cells
# RNA-seq: nCount_RNA, percent.mt
# ATAC-seq: nCount_ATAC, nucleosome_signal, TSS.enrichment
saline1 <- subset(
  x = saline1,
  subset = nCount_RNA > 2e2 & 
    nCount_RNA < 5e4 &
    nCount_ATAC > 2e2 &
    nCount_ATAC < 1e5 &
    percent.mt < 5 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3 &
    pct_reads_in_peaks > 15
)

saline2 <- subset(
  x = saline2,
  subset = nCount_RNA > 2e2 & 
    nCount_RNA < 5e4 &
    nCount_ATAC > 2e2 &
    nCount_ATAC < 1e5 &
    percent.mt < 5 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3 &
    pct_reads_in_peaks > 15
)

crotonate1 <- subset(
  x = crotonate1,
  subset = nCount_RNA > 2e2 & 
    nCount_RNA < 5e4 &
    nCount_ATAC > 2e2 &
    nCount_ATAC < 1e5 &
    percent.mt < 5 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3 &
    pct_reads_in_peaks > 15
)

crotonate2 <- subset(
  x = crotonate2,
  subset = nCount_RNA > 2e2 & 
    nCount_RNA < 5e4 &
    nCount_ATAC > 2e2 &
    nCount_ATAC < 1e5 &
    percent.mt < 5 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3 &
    pct_reads_in_peaks > 15
)

# Find most variable features in each cell
saline1 <- NormalizeData(saline1) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
saline2 <- NormalizeData(saline2) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
crotonate1 <- NormalizeData(crotonate1) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
crotonate2 <- NormalizeData(crotonate2) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

anchor_features <- SelectIntegrationFeatures(object.list = list(saline1, saline2,
                                                         crotonate1, crotonate2))
anchors <- FindIntegrationAnchors(object.list = list(saline1, saline2,
                                                     crotonate1, crotonate2), 
                                  anchor.features = anchor_features)
merged <- IntegrateData(anchorset = anchors)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- RunUMAP(merged, reduction = "pca", dims = 1:30)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)
# Save the seurat object
saveRDS(merged, file = "/home/bbasu/hpchome/R_codes/utsav/integrated_crotonate.rds")
merged <- readRDS("/home/bbasu/hpchome/R_codes/utsav/integrated_crotonate.rds")

head(merged@meta.data)
DimPlot_scCustom(merged)
DimPlot_scCustom(merged, group.by = "orig.ident")
DimPlot_scCustom(merged, group.by = "treatment")
######## Cell types annotation #################################################
#InN:Gad2, Gad1
#DG: Prox1, Ahcyl2
#Oligo: Mog, Prr5l
#CA1: Gm10754, Tbc1d1
#Astrocytes:Rgs20
#CA3:Nectin3, Il16
#SUB:Fn1, Stac
#Microglia:Inpp5d, C1qa
#OPC: Pdgfra
#APC: Igfbpl1
#ExN:Tshz2
#Cajal-Retzius:Trp73
#Endothelial: Cp
#CA2: Amigo2
#Ependymal cells: Dnah7b
#Assessing the expression of biomarkers in specific clusters
DefaultAssay(merged) <- "RNA"

#############Plot 1.Stacked Violin Plot All cluster markers#####################
colors <- paletteer::paletteer_d("ggsci::category20_d3")%>%
  head(12)
pdf(file = "/home/bbasu/hpchome/R_codes/utsav/plots/1_stacked_violin_markers_all_clusters.pdf")
Stacked_VlnPlot(merged, features = c("Prox1", "Mog", "Gm10754", "Nectin3", "Fn1", "Rgs20", "Gad2", "Inpp5d", 
                                     "Pdgfra", "Tshz2", "Igfbpl1","Trp73"),
                x_lab_rotate = 45, colors_use = colors) + ggtitle("Cluster Identification")
dev.off()
DimPlot_scCustom(merged)
FeaturePlot_scCustom(merged, features = c("Prox1"), reduction = "umap") #DG
FeaturePlot_scCustom(merged, features = c("Mog"), reduction = "umap") 
FeaturePlot_scCustom(merged, features = c("Gm10754"), reduction = "umap")
FeaturePlot_scCustom(merged, features = c("Nectin3"), reduction = "umap")
FeaturePlot_scCustom(merged, features = c("Fn1"), reduction = "umap")
FeaturePlot_scCustom(merged, features = c("Rgs20"), reduction = "umap")
FeaturePlot_scCustom(merged, features = c("Gad2"), reduction = "umap")
FeaturePlot_scCustom(merged, features = c("Inpp5d"), reduction = "umap")
FeaturePlot_scCustom(merged, features = c("Pdgfra"), reduction = "umap")
FeaturePlot_scCustom(merged, features = c("Tshz2"), reduction = "umap")
FeaturePlot_scCustom(merged, features = c("Igfbpl1"), reduction = "umap")
FeaturePlot_scCustom(merged, features = c("Trp73"), reduction = "umap")


merged <- RenameIdents(merged, `0` = "Dentate Gyrus", `1` = "Dentate Gyrus", `2` = "Oligodendrocytes", `3` = "CA1", `4` = "CA1", 
                       `5` = "CA3", `6` = "SUB", `7` = "Astrocytes", `8` = "Inhibitory neurons", `9` = "Microglia", 
                       `10` = "Oligodendrocytes", `11` = "OPC", `12` = "Excitatory cortical neurons", `13` = "Inhibitory neurons", `14` = "CA3", 
                       `15` = "CA1", `16` = "APC", `17` = "Excitatory cortical neurons", `18` = "Excitatory cortical neurons", `19` = "Inhibitory neurons", 
                       `20` = "CA3", `21` = "CA1", `22` = "Excitatory cortical neurons", `23` = "CA3", `24` = "Inhibitory neurons", `25` = "Astrocytes", 
                       `26` = "CA1", `27` = "Oligodendrocytes", `28` = "OPC", `29` = "Dentate Gyrus", `30` = "OPC", `31` = "Cajal-Retzius", `32` = "Dentate Gyrus", 
                       `33` = "Dentate Gyrus", `34` = "CA1", `35` = "Microglia", `36` = "OPC", `37` = "Microglia")
DimPlot_scCustom(merged)
merged$celltype <- Idents(merged)
merged$celltype.treatment <- paste(Idents(merged), merged$treatment, sep = "_")
head(merged@meta.data)

############################ Plot 2. Cell Number per treatment##################
cell.number <- as.data.frame(table(merged@meta.data$orig.ident))

pdf("/home/bbasu/hpchome/R_codes/utsav/plots/2_QualityControl_highQuality_cells.pdf")
ggbarplot(cell.number, x="Var1", y="Freq", fill = "Var1", color = "Var1", palette = "npg",
          sort.by.groups=FALSE,
          label = T, xlab = "", ylab = "Cell Number") + theme(legend.position="none") 
dev.off()
############################ Plot 3. UMAP Plot##################################
colors <- paletteer::paletteer_d("ggsci::category20_d3")%>%
  head(4)
p1 <- DimPlot_scCustom(merged, reduction = 'umap', group.by = "treatment", figure_plot = TRUE)
p2 <- DimPlot_scCustom(merged, reduction = 'umap', group.by = "orig.ident",
                       colors_use = colors, figure_plot = TRUE)

p1|p2

# Compare the dimplots
pdf("/home/bbasu/hpchome/R_codes/utsav/plots/3_Dimplot_grouped.pdf", height = 6, width = 15)
p1|p2
dev.off()

############################ Plot 4. UMAP Plot celltypes########################
length(unique(merged$celltype))
colors <- paletteer::paletteer_d("ggsci::category20_d3")%>%
  head(12)
pdf("/home/bbasu/hpchome/R_codes/utsav/plots/4_UMAP_plot_celltypes_labelled.pdf", height = 8, width = 10)
DimPlot_scCustom(merged, reduction = 'umap',
                 colors_use = colors, label = FALSE, figure_plot = TRUE)
dev.off()

## Create a gene activity matrix 
# peaks-seq
DefaultAssay(merged) <- "peaks"
merged <- RunTFIDF(merged)
merged <- FindTopFeatures(merged, min.cutoff = 5) #min.cutoff to ‘q75’ to use the top 25% all peaks
merged <- RunSVD(merged)

gene.activities <- GeneActivity(merged, extend.upstream = 2000, extend.downstream = 2000, 
                                biotypes = "protein_coding")
merged[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)
merged <- NormalizeData(merged, assay = 'GeneActivity')

DefaultAssay(merged) <- "GeneActivity"

###########################Plot 5. Gene Activity Plot###########################
pdf("/home/bbasu/hpchome/R_codes/utsav/plots/5_Gene_activity_chromatin_accessibility.pdf", height = 8, width = 10)
DotPlot_scCustom(merged, features = c("Prox1", "Mog", "Gm10754", "Nectin3", "Fn1", "Rgs20", "Gad2", "Inpp5d", 
                                      "Pdgfra", "Tshz2", "Igfbpl1","Trp73"),
                 flip_axes = TRUE, x_lab_rotate = TRUE)
dev.off()

################################################################################
# Save the seurat object
merged <- JoinLayers(merged)
saveRDS(merged, file = "/home/bbasu/hpchome/R_codes/utsav/integrated_crotonate_celltypes_labeled.rds")

# Combine all data into a single seurat object
# It takes long time!!
# merged <- merge(saline1, y = c(saline2, crotonate1, crotonate2),
#                 add.cell.id = c("saline1", "saline2", "crotonate1", "crotonate2"),
#                 project = "Crotonate_treatment")
# saveRDS(merged, file = "/home/bbasu/hpchome/R_codes/utsav/merged_crotonate.rds")
