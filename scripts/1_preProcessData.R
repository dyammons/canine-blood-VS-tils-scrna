#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")

#set output params
outName <- "tils_naive6_w_all_OS_PBMC"
subName <- "allCells"

#load in annotated tumor data aligned using cellranger 6.1.2
seu.obj.tumor <- readRDS(file = "../input/canine_naive_n6_annotated.rds")

#retrieve object aligned using cellranger 6.1.2
seu.obj.pbmc <- readRDS(file = "../input/220926_cfam_hVoWadj_res1.2_dims45_S3.rds")

#load in annotated object associated with the primary publication
reference <- readRDS(file = "../input/final_dataSet_HvO.rds")

#correct barcode discprencies to enable transfer of cell type labels
rownames(reference@meta.data) <- ifelse(grepl("_8", rownames(reference@meta.data)),paste0(substr(rownames(reference@meta.data), 1, nchar(rownames(reference@meta.data))-2),"_111") ,rownames(reference@meta.data))
rownames(reference@meta.data) <- ifelse(grepl("_9", rownames(reference@meta.data)),paste0(substr(rownames(reference@meta.data), 1, nchar(rownames(reference@meta.data))-2),"_8") ,rownames(reference@meta.data))
rownames(reference@meta.data) <- ifelse(grepl("_10", rownames(reference@meta.data)),paste0(substr(rownames(reference@meta.data), 1, nchar(rownames(reference@meta.data))-3),"_9") ,rownames(reference@meta.data))
rownames(reference@meta.data) <- ifelse(grepl("_111", rownames(reference@meta.data)),paste0(substr(rownames(reference@meta.data), 1, nchar(rownames(reference@meta.data))-4),"_10") ,rownames(reference@meta.data))

ct.l3 <- reference$celltype.l3
names(ct.l3) <- rownames(reference@meta.data)
seu.obj.pbmc <- AddMetaData(seu.obj.pbmc, metadata = ct.l3, col.name = "celltype.l3_pbmc")

ct.l1 <- reference$celltype.l1
names(ct.l3) <- rownames(reference@meta.data)
seu.obj.pbmc <- AddMetaData(seu.obj.pbmc, metadata = ct.l1, col.name = "celltype.l1_pbmc")

#select cells to integrate from tumor (only immune cells from tumor)
Idents(seu.obj.tumor) <- "celltype.l3"
table(seu.obj.tumor@active.ident)
seu.obj.tils <- subset(seu.obj.tumor,
                          idents = c("B cell" ,"CD4_act" , "CD4_fh" , "CD4_naive" ,"CD4_reg" ,
                                     "CD8_SPP1_hi" ,"CD8_eff" , "CD8_ex","IFN-TAM" , "CD4+_TIM",
                                     "NK","Neutrophil","Plasma cell","CD4-_TIM","T_IFN",
                                     "T_cycling","cDC1","cDC2","mregDC","pDC","preDC")
                         )

seu.obj.tils <- dataVisUMAP(seu.obj = seu.obj.tils, outDir = "../output/s3/", outName = "20231014_tils_2500", final.dims = 40, final.res = 0.6, stashID = "clusterID_2", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.6, n.neighbors = 75, assay = "integrated", saveRDS = F,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


### Fig 1a-1: create raw UMAP of tumor cells used in analysis
pi <- DimPlot(seu.obj.tils, 
              reduction = "umap", 
              group.by = "celltype.l1",
             # cols = c("#D4F3A3","#C8C3E2","#0066A5","#AC0535","#148B7D","#F17D00","hotpink"),
              pt.size = 0.1,
              label = F,
              label.box = F,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi) + NoLegend() + theme(axis.title = element_blank(), 
                                                 panel.border = element_blank(), 
                                                 plot.margin = unit(c(0, 0, 0, 0), "cm"),
                                                plot.title = element_text(size= 26)) + ggtitle("Tumor infiltrating immune cells (TILs)")
ggsave(paste("../output/", outName,"/",subName, "/tils_rawUMAP.png", sep = ""), width = 7, height = 7)

#select cells to integrate from pbmcs (all pbmcs from OS dogs)
seu.obj.pbmc.os <- subset(seu.obj.pbmc, 
                          subset = cellSource == "Osteosarcoma"
                         )

seu.obj.pbmc.os <- dataVisUMAP(seu.obj = seu.obj.pbmc.os, outDir = "../output/s3/", outName = "20231014_pbmc_os_2500", final.dims = 40, final.res = 0.2, stashID = "clusterID_2", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.6, n.neighbors = 75, assay = "integrated", saveRDS = F,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


### Fig 1a-2: create raw UMAP of pbmcs used in analysis
pi <- DimPlot(seu.obj.pbmc.os, 
              reduction = "umap", 
              group.by = "clusterID_2",
              pt.size = 0.1,
              label = F,
              label.box = F,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi) + NoLegend() + theme(axis.title = element_blank(), 
                                                 panel.border = element_blank(), 
                                                 plot.margin = unit(c(0, 0, 0, 0), "cm"),
                                                plot.title = element_text(size= 26)) + ggtitle("Circulating immune cells")
ggsave(paste("../output/", outName,"/",subName, "/pbmcOS_rawUMAP.png", sep = ""), width = 7, height = 7)


#integrate the data
seu.sub.list <- c(SplitObject(seu.obj.pbmc.os, split.by = "orig.ident"),SplitObject(seu.obj.tils, split.by = "name"))
seu.obj <- indReClus(seu.obj = NULL, outDir = "../output/s2/", subName = "bloodANDtils", preSub = T, seu.list = seu.sub.list, nfeatures = 2500, vars.to.regress = "percent.mt")

# seu.obj <- readRDS(file = "../output/s2/bloodANDtils_S2.rds")
clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = "bloodANDtils", test_dims = c(50,45,40,35,30,25), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "20230505_bloodANDtils_2500", final.dims = 40, final.res = 0.6, stashID = "clusterID_2", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.6, n.neighbors = 75, assay = "integrated", saveRDS = F,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#remove platelet contaminated cells
seu.obj.sub <- subset(seu.obj, invert = T,
                 subset = clusterID_2 == "24" | clusterID_2 == "17")
seu.obj.sub$clusterID_prePal_filter <- seu.obj.sub$clusterID_2
seu.obj.sub$clusterID_2 <- NULL
seu.obj.sub$clusterID <- NULL

#repeat integration with pal conaminated cells removed
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "../output/s2/", subName = "20230505_bloodANDtils_QCfiltered_2500", preSub = T, nfeatures = 2500,
                      vars.to.regress = "percent.mt"
                       )

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "20230505_bloodANDtils_QCfiltered_2500_QCfiltered_2500", final.dims = 40, final.res = 1.2, stashID = "clusterID", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.5, n.neighbors = 60, assay = "integrated", saveRDS = F,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#clean metadata a bit and save the object for down stream use
seu.obj@meta.data <- seu.obj@meta.data[,!grepl("DF|pANN", colnames(seu.obj@meta.data))]
saveRDS(seu.obj, "../output/s3/20230505_bloodANDtils_QCfiltered_2500_QCfiltered_2500_res1.2_dims40_dist0.5_neigh60_S3.rds")
