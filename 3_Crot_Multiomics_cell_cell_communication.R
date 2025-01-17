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
# Summary: Ligand Receptor Analysis
devtools::install_github("jinworks/CellChat")
library(Seurat)
library(CellChat)
library(ggpubr)
library(patchwork)
# Data input & processing and initialization of CellChat object
#===============================================================================
merged <- readRDS("/home/bbasu/hpchome/R_codes/utsav/integrated_crotonate_celltypes_labeled.rds")

# Saline group
# Subset the data based on treatment
DefaultAssay(merged)
saline <- subset(merged, subset = treatment == "Saline")
head(saline@meta.data)

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

#===============================================================================
# Step: 1 Process the saline data and create cellchat object
#===============================================================================
Idents(saline) <- saline$celltype
saline.input <- saline[["RNA"]]$data # normalized data matrix
saline.labels <- Idents(saline)
# create a dataframe of the cell labels
saline.meta <- data.frame(labels = saline.labels, row.names = names(saline.labels))
# Create cell-chat object
saline.cellChat <- createCellChat(object = saline, group.by = "ident", assay = "RNA")
# set the used database in the object
saline.cellChat@DB <- CellChatDB

#===============================================================================
# Step: 2 Preprocessing the expression data for cell-cell communication analysis
#===============================================================================
# subset the expression data of signaling genes for saving computation cost
saline.cellChat <- subsetData(saline.cellChat) # This step is necessary even if using the whole database
saline.cellChat <- identifyOverExpressedGenes(saline.cellChat, thresh.fc = 0.1, thresh.p = 0.05)
saline.cellChat <- identifyOverExpressedInteractions(saline.cellChat)
# The number of highly variable ligand-receptor pairs used for signaling inference is 1153

#===============================================================================
# Step: 3 Compute the communication probability and infer cellular communication network
#===============================================================================
saline.cellChat <- computeCommunProb(saline.cellChat, type = "triMean") # Takes time!
saline.cellChat <- filterCommunication(saline.cellChat, min.cells = 10)
# Infer the cell-cell communication at a signaling pathway level
# Inferred intercellular communication network of each ligand-receptor pair in the slot ‘net’ 
# each signaling pathway is stored in the slot ‘netP’
saline.cellChat <- computeCommunProbPathway(saline.cellChat)

#===============================================================================
# Step: 4 Calculate the aggregated cell-cell communication network
#===============================================================================
# by counting the number of links or summarizing the communication probability
saline.cellChat <- aggregateNet(saline.cellChat)
saline.cellChat <- netAnalysis_computeCentrality(saline.cellChat)
saveRDS(saline.cellChat, file = "/home/bbasu/hpchome/R_codes/utsav/saline.cellchat.rds")
saline.cellChat <- readRDS(file = "/home/bbasu/hpchome/R_codes/utsav/saline.cellchat.rds")

################################################################################
# Crotonate group
################################################################################
crotonate <- subset(merged, subset = treatment == "Crotonate")
head(crotonate@meta.data)
#===============================================================================
# Step: 1 Process the crotonate data and create cellchat object
#===============================================================================
Idents(crotonate) <- crotonate$celltype
crotonate.input <- crotonate[["RNA"]]$data # normalized data matrix
crotonate.labels <- Idents(crotonate)
# create a dataframe of the cell labels
crotonate.meta <- data.frame(labels = crotonate.labels, row.names = names(crotonate.labels))
# Create cell-chat object
crotonate.cellChat <- createCellChat(object = crotonate, group.by = "ident", assay = "RNA")
# set the used database in the object
crotonate.cellChat@DB <- CellChatDB

#===============================================================================
# Step: 2 Preprocessing the expression data for cell-cell communication analysis
#===============================================================================
# subset the expression data of signaling genes for saving computation cost
crotonate.cellChat <- subsetData(crotonate.cellChat) # This step is necessary even if using the whole database
crotonate.cellChat <- identifyOverExpressedGenes(crotonate.cellChat, thresh.fc = 0.1, thresh.p = 0.05)
crotonate.cellChat <- identifyOverExpressedInteractions(crotonate.cellChat)
# The number of highly variable ligand-receptor pairs used for signaling inference is 1108

#===============================================================================
# Step: 3 Compute the communication probability and infer cellular communication network
#===============================================================================
crotonate.cellChat <- computeCommunProb(crotonate.cellChat, type = "triMean")
crotonate.cellChat <- filterCommunication(crotonate.cellChat, min.cells = 10)
# Infer the cell-cell communication at a signaling pathway level
# Inferred intercellular communication network of each ligand-receptor pair in the slot ‘net’ 
# each signaling pathway is stored in the slot ‘netP’
crotonate.cellChat <- computeCommunProbPathway(crotonate.cellChat)

#===============================================================================
# Step: 4 Calculate the aggregated cell-cell communication network
#===============================================================================
# by counting the number of links or summarizing the communication probability
crotonate.cellChat <- aggregateNet(crotonate.cellChat)
crotonate.cellChat <- netAnalysis_computeCentrality(crotonate.cellChat)
saveRDS(crotonate.cellChat, file = "/home/bbasu/hpchome/R_codes/utsav/crotonate.cellchat.rds")
crotonate.cellChat <- readRDS(file = "/home/bbasu/hpchome/R_codes/utsav/crotonate.cellchat.rds")
#===============================================================================
#===============================================================================
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load the cellchat object
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saline.cellChat <- readRDS(file = "/home/bbasu/hpchome/R_codes/utsav/saline.cellchat.rds")
crotonate.cellChat <- readRDS(file = "/home/bbasu/hpchome/R_codes/utsav/crotonate.cellchat.rds")
object.list <- list(Saline = saline.cellChat, Crotonate = crotonate.cellChat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#===============================================================================
# Visualize
#===============================================================================

# Compare the total number of interactions and interaction strength
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), width = 0.4,
                           color.use = c("khaki3", "navyblue"))+labs_pubr()+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))

gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",width = 0.4,
                           color.use = c("khaki3", "navyblue"))+labs_pubr()+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))
gg1 + gg2

pdf(file = "/home/bbasu/hpchome/R_codes/utsav/plots/cellchat_interaction.pdf",
    height = 7, width = 12)
gg1+gg2
dev.off()


#  Identify the signaling changes of specific cell populations
gg1 <- netAnalysis_signalingChanges_scatter(cellchat,
                                            color.use = c("#9467BDFF", "#9467BDFF", "#9467BDFF"),
                                            idents.use = c("SUB"),
                                            dot.size = 8, label.size = 10)+
  labs(title = "Signalling changes of SUB (Crot. vs Saline)")+
  theme(axis.text = element_text(size = 25))+
  xlim(-2,14)+ylim(-2,8)+
  labs_pubr()

SUB_differential_signalling <- gg1$data
write.csv(SUB_differential_signalling, file = "/home/bbasu/hpchome/R_codes/utsav/SUB_differential_signalling.csv")

gg2 <- netAnalysis_signalingChanges_scatter(cellchat,
                                            color.use = c("#2CA02CFF", "#2CA02CFF", "#2CA02CFF"),
                                            idents.use = c("CA1"),
                                            dot.size = 8, label.size = 10)+
  labs(title = "Signalling changes of CA1 (Crot. vs Saline)")+
  theme(axis.text = element_text(size = 25))+
  xlim(-2,14)+ylim(-2,8)+
  labs_pubr()

CA1_differential_signalling <- gg2$data
write.csv(CA1_differential_signalling, file = "/home/bbasu/hpchome/R_codes/utsav/CA1_differential_signalling.csv")

gg3 <- netAnalysis_signalingChanges_scatter(cellchat,
                                            color.use = c("#D62728FF","#D62728FF","#D62728FF"),
                                            idents.use = c("CA3"),
                                            dot.size = 8, label.size = 10)+
  labs(title = "Signalling changes of CA3 (Crot. vs Saline)")+
  theme(axis.text = element_text(size = 25))+
  xlim(-2,14)+ylim(-2,8)+
  labs_pubr()

CA3_differential_signalling <- gg3$data
write.csv(CA3_differential_signalling, file = "/home/bbasu/hpchome/R_codes/utsav/CA3_differential_signalling.csv")

gg4 <- netAnalysis_signalingChanges_scatter(cellchat,
                                            color.use = c("#1F77B4FF","#1F77B4FF","#1F77B4FF"),
                                            idents.use = c("Dentate Gyrus"),
                                            dot.size = 8, label.size = 10)+
  labs(title = "Signalling changes of DG (Crot. vs Saline)")+
  theme(axis.text = element_text(size = 25))+
  xlim(-2,14)+ylim(-2,8)+
  labs_pubr()

DG_differential_signalling <- gg4$data
write.csv(DG_differential_signalling, file = "/home/bbasu/hpchome/R_codes/utsav/DG_differential_signalling.csv")


pdf(file = "/home/bbasu/hpchome/R_codes/utsav/plots/Signalling_changes.pdf",
    height = 8, width = 10)
gg1
gg2
gg3
gg4
dev.off()

# Show Glutamate signaling network
pathways.show <- c("Glutamate") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

pdf(file = "/home/bbasu/hpchome/R_codes/utsav/plots/Glutamate_Signalling_changes_Crot_vs_Saline_4celltypes.pdf",
    height = 6, width = 12)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                      thresh = 0.05,
                      sources.use = c("CA1", "CA3", "Dentate Gyrus", "SUB"),
                      targets.use = c("CA1", "CA3", "Dentate Gyrus", "SUB"),
                      remove.isolate = TRUE,
                      color.use = c(
                        "Dentate Gyrus"="#1F77B4FF",
                        "Oligodendrocytes"="#FF7F0EFF",
                        "CA1" = "#2CA02CFF",
                        "CA3"= "#D62728FF",
                        "SUB"="#9467BDFF",
                        "Astrocytes"="#8C564BFF",
                        "Inhibitory neurons"="#E377C2FF",
                        "Microglia"="#7F7F7FFF",
                        "OPC"="#BCBD22FF",
                        "Excitatory cortical neurons"="#17BECFFF",
                        "APC"="#AEC7E8FF",
                        "Cajal-Retzius"="#FFBB78FF"),
                      label.edge = TRUE)
}

#===============================================================================
# Identify the up-regulated and down-regulated signaling ligand-receptor pairs
#===============================================================================
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Crotonate"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, thresh.pc = 0.1, 
                                       thresh.fc = 0.1,thresh.p = 0.05, 
                                       group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in Crotonate
net.up <- subsetCommunication(cellchat, net = net, datasets = "Crotonate",
                              ligand.logFC = 0.1, receptor.logFC = 0.1)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Saline",
                                ligand.logFC = -0.1, receptor.logFC = -0.1)


pdf(file = "/home/bbasu/hpchome/R_codes/utsav/plots/Upregulated_signalling_DG-CA3.pdf",
    height = 8, width = 15)
netVisual_chord_gene(object.list[[2]], 
                     sources.use = c("Dentate Gyrus"), 
                     targets.use = c("CA3"),
                     slot.name = 'net',
                     thresh = 0.05,
                     net = net.up, lab.cex = 0.8, 
                     small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]),
                     color.use = colors)
dev.off()

pdf(file = "/home/bbasu/hpchome/R_codes/utsav/plots/Upregulated_signalling_CA3_CA1.pdf",
    height = 8, width = 15)
netVisual_chord_gene(object.list[[2]], 
                     sources.use = c("CA3"), 
                     targets.use = c("CA1"),
                     slot.name = 'net',
                     thresh = 0.05,
                     net = net.up, lab.cex = 0.8, 
                     small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]),
                     color.use = colors)
dev.off()

pdf(file = "/home/bbasu/hpchome/R_codes/utsav/plots/Upregulated_signalling_CA1_SUB.pdf",
    height = 8, width = 15)
netVisual_chord_gene(object.list[[2]], 
                     sources.use = c("CA1"), 
                     targets.use = c("SUB"),
                     slot.name = 'net',
                     thresh = 0.05,
                     net = net.up, lab.cex = 0.8, 
                     small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]),
                     color.use = colors)
dev.off()

# pdf(file = "/home/bbasu/hpchome/R_codes/utsav/plots/Upregulated_signalling_wordcloud.pdf",
#     height = 8, width = 15)
# computeEnrichmentScore(net.up,
#                        species = 'mouse',
#                        variable.both = TRUE)
# dev.off()




