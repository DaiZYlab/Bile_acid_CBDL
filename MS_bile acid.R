library(Seurat)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

Sham.data <- Read10X(data.dir = "~/R/BDL/Sham/filtered_feature_bc_matrix/")
Sham <- CreateSeuratObject(counts = Sham.data, project = "Sham", min.cells = 5, min.features = 100)

Sham[["percent.mt"]] <- PercentageFeatureSet(object = Sham, pattern = "^mt-")
VlnPlot(Sham, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Sham <- subset(Sham, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10)
Sham <- NormalizeData(object =Sham, normalization.method = "LogNormalize", scale.factor = 10000)

Sham <- FindVariableFeatures(object = Sham, selection.method = "vst", nfeatures = 2000)

Sham <- ScaleData(Sham, features = rownames(Sham), verbose = FALSE)
ElbowPlot(CBDL3W.singlet, 30)

Sham <- RunPCA(Sham, features = VariableFeatures(object = Sham))

#Clustering
ElbowPlot(Sham, 30)

Sham  <- RunTSNE(Sham , dims = 1:30, method = "FIt-SNE")
Sham <- RunUMAP(Sham, reduction = "pca", dims = 1:30)
DimPlot(Sham, reduction = "umap", label = TRUE)
DimPlot(Sham, reduction = "tsne", label = TRUE)
saveRDS(Sham, file = "~/R/BDL/MS/Sham.rds")


CBDL3W.data <- Read10X(data.dir = "~/R/BDL/CBDL3W/filtered_feature_bc_matrix/")
rownames(x = CBDL3W.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalB", replacement = "", 
                                                        x = rownames(x = CBDL3W.data[["Antibody Capture"]]))
CBDL3W <- CreateSeuratObject(counts = CBDL3W.data[["Gene Expression"]], project = "CBDL3W", min.cells = 5, min.features = 100)
CBDL3W<- NormalizeData(CBDL3W)
CBDL3W[["HTO"]] <- CreateAssayObject(CBDL3W.data[["Antibody Capture"]][, colnames(x = CBDL3W)])
CBDL3W <- NormalizeData(CBDL3W, assay = "HTO", normalization.method = "CLR")
FeatureScatter(CBDL3W, feature1 = "hto_3W-1", feature2 = "hto_3W-2", pt.size = 1)

CBDL3W <- HTODemux(CBDL3W, assay = "HTO", positive.quantile = 0.99)
table(CBDL3W$HTO_classification.global)
Idents(CBDL3W) <- "HTO_maxID"
RidgePlot(CBDL3W, assay = "HTO", features = rownames(CBDL3W[["HTO"]])[1:3], ncol = 3)
Idents(CBDL3W) <- "HTO_classification.global"
VlnPlot(CBDL3W, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

CBDL3W.singlet <- subset(CBDL3W, idents = "Singlet")
table(CBDL3W.singlet$HTO_maxID)

CBDL3W.singlet[["percent.mt"]] <- PercentageFeatureSet(object = CBDL3W.singlet, pattern = "^mt-")
VlnPlot(CBDL3W.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
CBDL3W.singlet  <- subset(CBDL3W.singlet, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10)

DefaultAssay(CBDL3W.singlet) <- "RNA"
CBDL3W.singlet <- NormalizeData(CBDL3W.singlet, normalization.method = "LogNormalize", scale.factor = 10000)
CBDL3W.singlet  <- FindVariableFeatures(CBDL3W.singlet, election.method = "vst", nfeatures = 2000)

CBDL3W.singlet <- ScaleData(CBDL3W.singlet, verbose = FALSE )
CBDL3W.singlet<- RunPCA(object = CBDL3W.singlet, npcs = 30, verbose = FALSE, approx=FALSE)

ElbowPlot(CBDL3W.singlet, 30)

CBDL3W.singlet <- RunTSNE(CBDL3W.singlet, dims = 1:18, method = "FIt-SNE")
CBDL3W.singlet <- RunUMAP(CBDL3W.singlet, reduction = "pca", dims = 1:18)
CBDL3W.singlet<- FindNeighbors(CBDL3W.singlet, reduction = "pca", dims = 1:18)

CBDL3W.singlet<- FindClusters(CBDL3W.singlet, resolution = 0.5)

DimPlot(CBDL3W.singlet, reduction = "tsne", split.by= "HTO_maxID", label = TRUE, pt.size = 0.1)
DimPlot(CBDL3W.singlet, reduction = "umap", split.by= "HTO_maxID", label = TRUE, pt.size = 1)

DimPlot(CBDL3W.singlet, reduction = "tsne", label = TRUE, pt.size = 0.1)
DimPlot(CBDL3W.singlet, reduction = "umap", label = TRUE, pt.size = 1)

saveRDS(CBDL3W.singlet, file = "~/R/BDL/MS/CBDL3W.singlet.rds")

# Only merge Sham and new 3 weeks
# integrate 
ifnb.list = c(CBDL3W.singlet,Sham) 
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

merge.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
CBDL3<- IntegrateData(anchorset = merge.anchors)
DefaultAssay(object = CBDL3) <- "integrated"
CBDL3 <- FindVariableFeatures(object = CBDL3, selection.method = "vst", nfeatures = 2000)
CBDL3 <- ScaleData(object = CBDL3, verbose = FALSE)
CBDL3<- RunPCA(object = CBDL3, npcs = 40, verbose = FALSE)
ElbowPlot(CBDL3, 40)
CBDL3 <- RunUMAP(CBDL3, reduction = "pca", dims = 1:32)
CBDL3 <- RunTSNE(CBDL3, dims = 1:32, method = "FIt-SNE")
CBDL3 <- FindNeighbors(CBDL3, reduction = "pca", dims = 1:32)
CBDL3 <- FindClusters(CBDL3, resolution = 0.6)
DimPlot(CBDL3, reduction = "tsne", split.by = "orig.ident", label = TRUE, pt.size = 0.1)
DimPlot(CBDL3, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.1)

DimPlot(CBDL3, reduction = "umap",  label = TRUE, pt.size = 0.1)
DimPlot(CBDL3, reduction = "tsne",  pt.size = 0.1)
DimPlot(CBDL3, reduction = "tsne",  pt.size = 0.1, label = TRUE)
DimPlot(CBDL3, reduction = "tsne", split.by = "orig.ident",  pt.size = 0.1)
DimPlot(CBDL3, reduction = "tsne", group.by = "orig.ident",  pt.size = 0.1)


saveRDS(CBDL3, file = "~/R/BDL/MS/CBDL3.rds")

CBDL3$orig.ident <- factor(CBDL3$orig.ident, levels = c("Sham","CBDL3W"))

DefaultAssay(object = CBDL3) <- "RNA"

library(UCell)
markers <- list()
markers$EC <- c("Cdh5", "Pecam1")
markers$SMC <- c("Acta2",	"Myh11", "Cnn1", "Itga7", "Ntrk3")
markers$Pericyte <- c("Notch3",	"Pdgfrb", "Cspg4", "Lamc3", "Trpc6", "Higd1b", "Fam162b")
markers$Fibroblast <- c("Fn1",	"Col1a1")
markers$LEC <- c("Prox1", "Ccl21a", "Hoxd3", "Reln", "Ptx3", "Nts")
markers$PMN <- c("S100a9", "S100a8", "Mmp9", "Gm5416", "Asprv1", "Stfa2")
markers$AT2 <- c("Sftpc", "Sftpb", "Lamp3")
markers$AT1 <- c("Aqp5", "Hopx", "Ager", "Rtkn2")
markers$Treg <- c("Foxp3", "Ctla4","Il2ra", "Folr4", "Icos", "Tnfrsf18", "Tnfrsf4")
markers$Cd4_T <- c("Cd40lg", "Cd4","Trat1", "Nsg2", "Lef1", "Cd3e", "Ccr7")
markers$Cd8_T <- c("Cd8a", "Cd8b1","Cd3e", "Cd3d", "Trac", "Cd27", "Ccr7")
markers$NK <- c("Klrc1","Klrb1c", "Klre1", "Ncr1", "Xcl1", "Tbx21")
markers$Bcell <- c("Ighm", "Cd19", "Pax5", "Cd79a", "Igkc")
markers$DC <- c("Xcr1", "Gucy2c", "Sept3", "Gpr82", "flt3", "Cd209a", "Cd209d", "Clec10a", "Mgl2", "Ccl17")
markers$iMac <- c("C1qa", "Fcrls", "Folr2", "Ccl12", "Ms4a7", "C1qc", "C1qb")
markers$aMac <- c("Olr1", "F7", "Marco", "Atp6v0d2", "Gpnmb", "Pparg", "Ear2")
markers$inMono <- c("Ifitm6", "Fcnb", "Gpr141", "Gm9733", "Mnda", "Tifab", "Mmp8")
markers$Master <- c("Tpsab1", "Mrgprb2", "Mrgprb1", "Tpsb2", "Cpa3")
markers$AF1 <- c("Wnt2", "Meox2","Slc38a5", "Tcf21", "Col13a1", "Figf", "Enpep", "Adh1")
markers$AdvFib <- c("Aspn", "Serpinf1","Pi16", "Sfrp2")
markers$AF2 <- c("Abca8a", "Agtr2", "Lum", "Scara5", "Entpd2", "Clec3b", "Htra3")
markers$Basal <- c("Krt5", "Pitx1","Sprr2a3", "Capns2", "Pkp1", "Sprr2b", "Trp63")
markers$Cilliated <- c("Ankrd65", "Gm867","Aoc1", "Ldlrad1", "Fam216b", "Cdhr3", "Barx2")
markers$secretory <- c("Gabrp", "Gm15883", "Slc16a11", "Pon1", "Gsta3", "Scga1a1", "Scgb3a2")
markers$Mesothelial <- c("Lrm4", "Upk3b","Wt1", "Slurp1", "Tmem151a", "Aldh1a2", "Cldn15")
markers$Platelet <- c("Gp5", "Tubb1","Gp6", "Ly6g6f", "Mpl", "F2rl2", "Gfi1b")

CBDL3<- AddModuleScore_UCell(CBDL3, features = markers)
signature.names <- paste0(names(markers), "_UCell")
VlnPlot(CBDL3, features = c("EC_UCell", "LEC_UCell"), group.by = "seurat_clusters")
VlnPlot(CBDL3, features = c("inMono_UCell",  "iMac_UCell", "aMac_UCell", "NK_UCell"), group.by = "seurat_clusters")
VlnPlot(CBDL3, features = c( "Bcell_UCell", "PMN_UCell", "DC_UCell",  "Treg_UCell", "Cd8_T_UCell", "Cd4_T_UCell"), group.by = "seurat_clusters")
VlnPlot(CBDL3, features = c("AT1_UCell",  "AT2_UCell", "Cilliated_UCell", "secretory_UCell"), group.by = "seurat_clusters")
VlnPlot(CBDL3, features = c("AF1_UCell", "AdvFib_UCell", "AF2_UCell", "Pericyte_UCell","SMC_UCell"), group.by = "seurat_clusters")
VlnPlot(CBDL3, features = c("Master_UCell", "Mesothelial_UCell", "Platelet_UCell"), group.by = "seurat_clusters")

markers.to.plot <- c("Cdh5", "Acta2", "Notch3","Cd68", "Aqp5", "Sftpb", "Prox1","Chst1", "Car4", "Col1a1", "Scgb1a1", "Slc6a2", "Selp", "Cd3g", "Myh11", "Fn1", "Top2a", "Klrc1", "Ighm", "Cd14", "Cd3e", "Cd7", "Pdgfrb", "Sftpc", "Cspg4", "Rgs5", "Des", "Fabp4", "S100a9", "Cyb561a3", "Cpa3", "Gata3", "Olr1", "Ccl17","Dynlrb2", "Fcnb", "C1qc", "Mal", "Gzma", "Ager", "Dcn", "Msln", "Ctsw", "Icos", "Cxcr4", "Cd8a", "Scgb3a1", "Bmp5", "Lum","Nkd2", "Moxd1", "Vwf", "Gja5","Serpinf1", "Apoe", "Aspn", "Gpc3")
DotPlot(CBDL3, features = rev(markers.to.plot), dot.scale = 8) + RotatedAxis()

CBDL3.markers <- FindAllMarkers(CBDL3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CBDL3.markers, file = "~/R/BDL/MS/CBDL3.markers.csv")

CBDL3 <- RenameIdents(CBDL3, '0' = "EC", '1' = "PMN", '2' = "Mac",  '3' = "Mac", '4' = "EC", '5' = "Fib", '6' = "EC", '7' = "Tcell", '8' = "EC", '9' = "AT2", '10' =  "NK", '11' = "Bcell", '12' = "Tcell", '13' = "Pericyte", '14' = "PMN", '15' = "EC", '16' = "AT1", '17' = "DC", '18' = "Tcell", '19' = "Mac", '20' ="Cilliated_Secretory", '21' ="SMC", '22' = "LEC", '23'= "Bcell", '24'= "Master", '25'= "AT2", '26'= "EC", '27'= "Meso")

CBDL3$MergeCellType <- Idents(CBDL3)
markers.to.plot <- c("Cdh5","Pecam1",  "S100a9", "S100a8", "Cd68", "Mrc1", "Col1a2",  "Col1a1","Cd3e", "Cd3g","Slc34a2", "Sftpb","Gzma","Klrb1c","Ighm", "Igkc", "Pdgfrb", "Notch3","Aqp5","Ager", "Cd83", "Tnip3", "Cyp2f2", "Scgb3a2", "Acta2", "Myh11",  "Prox1", "Ccl21a", "Cpa3","Cd200r3", "Upk3b", "Wt1")

#merge celltype
DefaultAssay(object = CBDL3) <- "integrated"
CBDL3 <- FindClusters(CBDL3, resolution = 0.6)
CBDL3$MergeCellType <- Idents(CBDL3)
DimPlot(CBDL3, reduction = "umap", label = TRUE, pt.size = 0.1)
Idents(CBDL3) <- "CellType"
DimPlot(CBDL3, reduction = "umap", label = TRUE, pt.size = 0.1)
Idents(CBDL3) <- "MergeCellType"
DimPlot(CBDL3, reduction = "umap", label = TRUE, pt.size = 0.1)
DimPlot(CBDL3, reduction = "umap", label = TRUE, split.by = "orig.ident", pt.size = 0.1)

saveRDS(CBDL3, file = "~/R/BDL/MS/CBDL3.rds")

DefaultAssay(object = CBDL3) <- "RNA"

top20 <- CBDL3.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(CBDL3, features = top20$gene) + NoLegend()

#proportion change, merge cell types
CBDL3.propotion <- prop.table(table(Idents(CBDL3),CBDL3$orig.ident))
write.csv(CBDL3.propotion, file = "~/R/BDL/MS/CBDL3.propotion.csv")
table(Idents(CBDL3),CBDL3$orig.ident)
cell.prop<-as.data.frame(prop.table(table(Idents(CBDL3), CBDL3$orig.ident)))
colnames(cell.prop)<-c("cluster","orig.ident","proportion")
ggplot(cell.prop,aes(orig.ident,proportion,fill=cluster)) + geom_bar(stat="identity",position="fill") + ggtitle("") + theme_bw() + theme(axis.ticks.length=unit(0.5,'cm')) + guides(fill=guide_legend(title=NULL))

library("scProportionTest")
prop_test <- sc_utils(CBDL3)
prop_test <- permutation_test(prop_test, cluster_identity = "MergeCellType", sample_1 = "Sham", sample_2 = "CBDL3W", sample_identity = "orig.ident")
permutation_plot(prop_test)

Idents(CBDL3) <- "MergeCellType"

#DEGs for all cell types
Idents(CBDL3) <- "orig.ident"
CBDL3$celltype.test <- paste(Idents(CBDL3), CBDL3$MergeCellType, sep = "_")
Idents(CBDL3) <- "celltype.test"

# Get unique cell types from CBDL3$MergeCellType
cell_types <- unique(CBDL3$MergeCellType)

# Create an empty list to store results for each cell type
results_list <- list()

# Loop through each cell type and perform differential expression analysis
for (cell_type in cell_types) {
  celltype_Sham <- paste("Sham", cell_type, sep = "_")
  celltype_CBDL3W <- paste("CBDL3W", cell_type, sep = "_")
  
  # Perform differential expression analysis
  results <- FindMarkers(
    CBDL3, 
    ident.1 = celltype_CBDL3W, 
    ident.2 = celltype_Sham, 
  )
  
  # Store the results in the list with the cell type as the name
  results_list[[cell_type]] <- results
}

# Loop through each cell type's differential expression results
#for (cell_type in names(results_list)) {
  # Exclude mitochondrial genes
  non_mito_genes <- !grepl("^mt-", rownames(results_list[[cell_type]]), ignore.case = TRUE)
  results_list[[cell_type]] <- results_list[[cell_type]][non_mito_genes, ]
#}

# Now, results_list contains DE analysis for each cell type without mitochondrial genes

library(scRNAtoolVis)
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")


results_list <- lapply(results_list, function(cell_type) {
cell_type$gene <- rownames(cell_type)
return(cell_type)
}) # put the rownames (gene symbol) to a new column gene to avoid the change of gene name after combination into one data frame.

# Ensure cell_types array is defined and matches the names in results_list
cell_types <- names(results_list)
# Add a cell type column to each data frame and combine them into one
degs_combined <- do.call(rbind, lapply(cell_types, function(cell_type) {
  df <- results_list[[cell_type]]
  df$CellType <- cell_type
  return(df)
}))

degs_combined$cluster <- degs_combined$CellType # need cluster column
write.csv(degs_combined, file = "~/R/BDL/MS/degs_combined.csv")

# Assuming degs_combined is your data frame and it has a column named 'cluster'
cell_types <- unique(degs_combined$cluster)

# Split cell types into two groups (first 6 and next 6)
first_six_cell_types <- cell_types[1:6]
second_six_cell_types <- cell_types[7:12]

# Subset data for the first six cell types
degs_first_six <- subset(degs_combined, cluster %in% first_six_cell_types)

# Plot using jjVolcano for the first six cell types
jjVolcano(diffData = degs_first_six, fontface = 'italic', log2FC.cutoff = 0.585)

# Subset data for the next six cell types
degs_second_six <- subset(degs_combined, cluster %in% second_six_cell_types)

# Plot using jjVolcano for the next six cell types
jjVolcano(diffData = degs_second_six, fontface = 'italic', log2FC.cutoff = 0.585)

mygene <- c('Nppa', "Nppb", "Vegfa", "Vegfc", "Postn", "Cxcl12")


#DEGs EC
head(CBDL3)
Idents(CBDL3) <- "orig.ident"
CBDL3$celltype.stim <- paste(Idents(CBDL3), CBDL3$MergeCellType, sep = "_")

Idents(CBDL3) <- "celltype.stim"

CBDL3.response_AT2 <- FindMarkers(CBDL3, ident.1 = "Sham_AT2", ident.2 = "CBDL3W_AT2")
write.csv(CBDL3.response_AT2, file = "~/R/BDL/CBDL3.response_AT2.csv")

# AT1 AT2
Idents(CBDL3) <- "MergeCellType"

CBDL3.AT <- subset(CBDL3, idents = c("AT1","AT2"))

DefaultAssay(object = CBDL3.AT) <- "integrated"
CBDL3.AT <- FindVariableFeatures(object = CBDL3.AT, selection.method = "vst", nfeatures = 2000)
CBDL3.AT <- ScaleData(object =CBDL3.AT, verbose = FALSE)
CBDL3.AT<- RunPCA(object = CBDL3.AT, npcs = 10, verbose = FALSE)
ElbowPlot(CBDL3.AT, 10)
CBDL3.AT <- RunUMAP(CBDL3.AT, reduction = "pca", dims = 1:4)
CBDL3.AT <- RunTSNE(CBDL3.AT, dims = 1:4, method = "FIt-SNE")
CBDL3.AT <- FindNeighbors(CBDL3.AT, reduction = "pca", dims = 1:4)
CBDL3.AT <- FindClusters(CBDL3.AT, resolution = 0.1)
DimPlot(CBDL3.AT, reduction = "umap", label = TRUE)
DimPlot(CBDL3.AT, reduction = "tsne", label = TRUE)
saveRDS(CBDL3.AT, file = "~/R/BDL/MS/CBDL3.AT.rds")
DefaultAssay(object = CBDL3.AT) <- "RNA"
CBDL3.AT.markers <- FindAllMarkers(CBDL3.AT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CBDL3.AT.markers, file = "~/R/BDL/MS/CBDL3.AT.markers.csv")

DefaultAssay(object = CBDL3.AT) <- "RNA"

VlnPlot(CBDL3.AT, features = c("Sftpc", "Sftpb", "Aqp5", "Ager","Krt8", "Cdkn1a", "Ly6a", "Soat1", "Cav1"), ncol = 2, pt.size = 0)
markers.to.plot <- c("Sftpb", "Sftpc","Aqp5","Ager", "Cav1", "Krt8", "Cdkn1a", "Ly6a", "Soat1")
DotPlot(CBDL3.AT, features = rev(markers.to.plot),  dot.scale = 8) + RotatedAxis()
CBDL3.AT <- RenameIdents(CBDL3.AT, '0' = "AT2", '1' = "AT1", '2' = "AT1/2", '3' = "AT1/2" )
DimPlot(CBDL3.AT, reduction = "tsne", split.by = "orig.ident", label = TRUE)
DimPlot(CBDL3.AT, reduction = "umap", split.by = "orig.ident", label = TRUE)

saveRDS(CBDL3.AT, file = "~/R/BDL/CBDL3.AT.rds")
CBDL3.AT$ATCelltype <- Idents(CBDL3.AT)  

top20 <- CBDL3.AT.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(CBDL3.AT, features = top20$gene) + NoLegend()


#proportion change, merge cell types
CBDL3.AT.propotion <- prop.table(table(Idents(CBDL3.AT),CBDL3.AT$orig.ident))
write.csv(CBDL3.AT.propotion, file = "~/R/BDL/MS/CBDL3.AT.propotion.csv")
table(Idents(CBDL3.AT),CBDL3.AT$orig.ident)
cell.prop<-as.data.frame(prop.table(table(Idents(CBDL3.AT), CBDL3.AT$orig.ident)))
colnames(cell.prop)<-c("cluster","orig.ident","proportion")
ggplot(cell.prop,aes(orig.ident,proportion,fill=cluster)) + geom_bar(stat="identity",position="fill") + ggtitle("") + theme_bw() + theme(axis.ticks.length=unit(0.5,'cm')) + guides(fill=guide_legend(title=NULL))

cell.num <- as.data.frame(table(Idents(CBDL3.AT), CBDL3.AT$orig.ident))
colnames(cell.num)<-c("cluster","orig.ident","cell.number")
ggplot(cell.num,aes(orig.ident,cell.number,fill=cluster)) + geom_bar(stat="identity") + ggtitle("") + theme_bw() + theme(axis.ticks.length=unit(0.5,'cm')) + guides(fill=guide_legend(title=NULL))

#DEGs of AT2 cells
Idents(CBDL3.AT) <- "MergeCellType"
CBDL3.AT$celltype.orig.ident <- paste(Idents(CBDL3.AT), CBDL3.AT$orig.ident, sep = "_")
Idents(CBDL3.AT) <- "celltype.orig.ident"
AT2.response <- FindMarkers(CBDL3.AT, ident.1 = "AT2_CBDL3W", ident.2 = "AT2_Sham")
head(AT2.response, n = 40)
write.csv(AT2.response, file = "~/R/BDL/MS/AT2.response.csv")

#Identify differential expressed genes across conditions
FeaturePlot(CBDL3.AT, features = c("Cebpa", "Tead1","Sftpc", "Hopx", "Gdf15", "Mfap4"), split.by = "orig.ident", 
            cols = c("grey", "red"))

VlnPlot(CBDL3.AT, features = c("Cebpa", "Cxcl12", "Tead1", "S100a9","Sftpc","Mfap4", "Epas1"  ), split.by = "orig.ident", group.by = "ATCelltype",
        pt.size = 0, combine = TRUE, ncol = 2)
VlnPlot(CBDL3.AT, features = c("Cebpa", "Irx1", "Irx2", "Atf4", "Foxp1", "Etv5"), split.by = "orig.ident", group.by = "ATCelltype",
        pt.size = 0.1, combine = TRUE, ncol = 2)
VlnPlot(CBDL3.AT, features = c("Sftpc", "Sftpb"), split.by = "orig.ident", group.by = "ATCelltype",
        pt.size = 0.1, combine = TRUE, ncol = 1)
Idents(CBDL3.AT) <- "orig.ident"

CBDL3.AT.averages <- AverageExpression(CBDL3.AT)
write.csv(CBDL3.AT.averages[["RNA"]], file = "~/R/BDL/CBDL3.AT.averages.csv")
Idents(CBDL3.AT) <- "CellType"


# Filter DEGs between "Sham_AT2" and "CBDL3W_AT2"
significant_DEGs <- subset(AT2.response, p_val_adj < 0.05 & (avg_log2FC > 0.585 | avg_log2FC < -0.585))

write.csv(significant_DEGs, file = "~/R/BDL/MS/significant_DEGs_AT2.csv")

significant_DEGs$abs_log2FC <- abs(significant_DEGs$avg_log2FC)  # Add a column with absolute log2 fold change
significant_DEGs_no_mt <- subset(significant_DEGs, !grepl("^mt-", rownames(significant_DEGs)))
write.csv(significant_DEGs_no_mt, file = "~/R/BDL/MS/significant_DEGs_no_mt_AT2.csv")

# Scale data for top genes
Idents(CBDL3.AT) <- "MergeCellType"
AT2_CBDL <- subset(CBDL3.AT, subset = MergeCellType == "AT2")

AT2_CBDL <- ScaleData(AT2_CBDL, features = rownames(significant_DEGs_no_mt))

Idents(AT2_CBDL) <- "celltype.orig.ident"

top_genes <- head(rownames(significant_DEGs_no_mt), 100)

DoHeatmap(AT2_CBDL, features = top_genes) + NoLegend()

#SCENIC

#if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
#devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE, force = TRUE)
getwd()
setwd("~/R/BDL/SCENIC")

library(SCopeLoomR)
exprMat  <-  as.matrix(CBDL3.AT@assays$RNA@data)
#dim(exprMat)

head(CBDL3.AT)
cellInfo <-  CBDL3.AT@meta.data[,c(2,3,1,16)]
colnames(cellInfo)=c( 'nGene' ,'nUMI','Group','CellType')
head(cellInfo)
table(cellInfo$CellType)
saveRDS(cellInfo, file="int/cellInfo.Rds")

library(SCENIC)
org <- "mgi"
#dir.create("cisTarget_databases"); setwd("cisTarget_databases")
dbDir <- "~/cisTarget_databases"
myDatasetTitle <- "SCENIC"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=1) 

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=10*.03*ncol(exprMat),
                           minSamples=ncol(exprMat)*.03)
#genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=10*.03*ncol(exprMat),#minSamples=ncol(exprMat)*.03), TF 177, it did not finish for 5 days
#genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=20*.05*ncol(exprMat),#minSamples=ncol(exprMat)*.05), TF 32, finally 2 TF
#genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=15*.05*ncol(exprMat),#minSamples=ncol(exprMat)*.05), TF 56, finally 4 TF, 666 genes
#genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,minCountsPerGene=10*.05*ncol(exprMat),#minSamples=ncol(exprMat)*.05), TF 94, finally 26 TF, 1102 genes

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered) 

rm(exprMat) # close exprMat

runCorrelation(exprMat_filtered, scenicOptions)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

exprMat_filtered <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered, scenicOptions) # take days to run
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
exprMat_log <- log2(exprMat+1)# if not log normalized.
rm(exprMat) # close exprMat

#scenicOptions <- readRDS("int/scenicOptions.Rds")

scenicOptions@settings$verbose <- TRUE
#scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

#scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
scenicOptions@settings$dbs <- scenicOptions@settings$dbs
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
#scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget"))
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions) 
library(foreach)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")


aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
savedSelections <- shiny::runApp(aucellApp)# determine the threshold for each TF

# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC")

scenicOptions@fileNames$output["loomFile",] <- "SCENIC1.loom"
export2loom(scenicOptions, exprMat)

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

#Projection the AUC and TF expression onto t-SNEs
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log) # default t-SNE
savedSelections <- shiny::runApp(aucellApp)
print(tsneFileName(scenicOptions))
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Tead1", "Klf6", "Cebpa")],], plots="Expression")

#par(bg = "black")
par(mfrow=c(1,2))

regulonNames <- c( "Tead1", "Cebpa")
cellCol <- plotEmb_rgb(scenicOptions,regulonNames, aucType = "Binary",aucMaxContrast = 1,offColor = "lightgray",showLegend = F)

regulonNames <- list(red=c("Sox10", "Sox8"),
                     green=c("Irf1"),
                     blue=c( "Tef"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")

#Average Regulon Activity by cluster
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")

#monocyle3
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
cds <- as.cell_data_set(CBDL3.AT)
#cds <- preprocess_cds(cds, num_dim = 3)
#cds <- reduce_dimension(cds)
#cds <- cluster_cells(cds, resolution = 0.000198)
cds <- cluster_cells(cds, resolution = 0.000185)
cds <- learn_graph(cds)
plot_cells(cds)

cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

max.AT2 <- which.max(unlist(FetchData(CBDL3.AT, "Sftpc")))
max.AT2 <- colnames(CBDL3.AT)[max.AT2]
cds <- order_cells(cds, root_cells = max.AT2)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = TRUE)

#Volcano plot for AT2 DEGs

df = read.table(file = "~/R/BDL/MS/AT2.response.txt", header = T,row.names = 1, sep = '\t')
de <- df[complete.cases(df), ]

# The basic scatter plot: x is "log2FoldChange", y is "pvalue"
ggplot(data=de, aes(x=avg_log2FC, y=p_val_adj)) + geom_point()

# Convert directly in the aes()
p <- ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val_adj))) + geom_point()
# Add more simple "theme"
p <- ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val_adj))) + geom_point() + theme_minimal()

# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.3, 0.3), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$avg_log2FC > 0.6 & de$p_val_adj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$avg_log2FC < -0.6 & de$p_val_adj < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("green", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("green", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de$delabel <- row.names(de)

de$delabel[de$diffexpressed != "NO"] <- de$delabel[de$diffexpressed != "NO"]

ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()

# Finally, we can organize the labels nicely using the "ggrepel" package and the geom_text_repel() function
# load library
library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=de, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("green", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
ggsave("~/R/BDL/MS/DEG_AT2.pdf", width = 10, height = 8)

#pathway AT2_DEG
df_AT2_KEGG = read.table(file = "~/R/BDL/MS/AT2_KEGG_up.txt", header = T, sep = '\t')
p <- ggplot(df_AT2_KEGG, aes(x=Fold.Enrichment, y=Up_KEGG)) + geom_point(aes(size=Count, color=FDR)) + scale_color_gradient(low = "red", high = "blue")
ggplot(df_AT2_KEGG, aes(x=Fold.Enrichment, y=Up_KEGG)) + geom_point(aes(size=Count, color=FDR)) + scale_color_gradient(low = "red", high = "blue")

ggsave("~/R/BDL/MS/KEGG_up_AT2.pdf", width = 8, height = 8, dpi=600)

df_AT2_KEGG_Dn = read.table(file = "~/R/BDL/MS/AT2_Dn_KEGG.txt", header = T, sep = '\t')
p <- ggplot(df_AT2_KEGG_Dn, aes(x=Fold.Enrichment, y=Dn_KEGG)) + geom_point(aes(size=Count, color=FDR)) + scale_color_gradient(low = "red", high = "blue")
ggplot(df_AT2_KEGG_Dn, aes(x=Fold.Enrichment, y=Dn_KEGG)) + geom_point(aes(size=Count, color=FDR)) + scale_color_gradient(low = "red", high = "blue")

ggsave("~/R/BDL/MS/KEGG_Dn_AT2.pdf", width = 8, height = 8, dpi=600)

