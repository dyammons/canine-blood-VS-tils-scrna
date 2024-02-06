#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")

### Analysis note: 
# This script loads in the previously processed Seurat object (./output/s3/230816_duod_h3c4_NoIntrons_res1.3_dims40_dist0.3_neigh50_S3.rds)
# then subsets on T cells and generates all figures assocaited with Figure 3

############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin all cells analysis   ######## <<<<<<<<<<<<<<
############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in data
seu.obj <- readRDS("./output/s3/20230505_bloodANDtils_QCfiltered_2500_QCfiltered_2500_res1.2_dims40_dist0.5_neigh60_S3.rds")
seu.obj$conSense <- ifelse(seu.obj$cellSource == "TILs", seu.obj$celltype.l3,seu.obj$celltype.l3_pbmc)
outName <- "tils_naive6_w_all_OS_PBMC"
subName <- "allCells"

Idents(seu.obj) <- "orig.ident"
seu.obj <- RenameIdents(seu.obj, c("run_count_osa_pbmc_1" = "bl_1", "run_count_osa_pbmc_2" = "bl_2", 
                                   "run_count_osa_pbmc_3" = "bl_3", "run_count_osa_pbmc_4" = "bl_4",
                                   "run_count_osa_pbmc_5" = "bl_5","run_count_pbmc_tp1_pt1" = "bl_6",
                                   "run_count_pbmc_tp1_pt2" = "bl_7", "run_count_pbmc_tp1_pt3" = "bl_8", 
                                   "run_count_pbmc_tp1_pt4" = "bl_9", "run_count_pbmc_tp1_pt5" = "bl_10",
                                   "run_count_tumor_no_tx_1_1" = "tils_1","run_count_tumor_no_tx_1_2" = "tils_1",
                                   "run_count_tumor_no_tx_2_1" = "tils_2", "run_count_tumor_no_tx_2_2" = "tils_2", 
                                   "run_count_tumor_no_tx_4" = "tils_3", "run_count_tumor_no_tx_5" = "tils_4",
                                   "run_count_no_tx_tumor_6" = "tils_5","run_count_tumor_no_tx_7" = "tils_6")
                       )
seu.obj$name <- Idents(seu.obj)


### Fig extra - create raw UMAP
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID",
              #cols = levels(seu.obj$dcColz),
              pt.size = 0.5,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- cusLabels(plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8, labCol = "black") + NoLegend() + theme(axis.title = element_blank(),
                                                                                                                    panel.border = element_blank())
ggsave(paste("./output/", outName, "/", subName, "/", "rawUMAP.png", sep = ""), width = 7, height = 7)

#get label cpnsensous
seu.obj$conSense <- ifelse(seu.obj$cellSource == "TILs", seu.obj$celltype.l3,seu.obj$celltype.l3_pbmc)

# ### DO NOT RUN -- export the seu.obj$conSense metadata slot into a .csv to save as metadata (cellTypez.csv is provided through GitHub)
# Idents(seu.obj) <- "conSense"
# df <- as.data.frame(c(c("CD8_SPP1_hi", "CD8+ Memory","CD8+ Effector", "CD8+ gd T cell","CD8+ Naive", "NK cell", "NK T cell", "NK cell","CD8_eff","CD8_ex", "NK"),
#                 c("CD4+ Naive", "CD4+ T reg", "CD4+ TCM", "CD4+ TEM", 
#                                     "CD4+ TEM, Th1-like","CD4+ TEM, Th17-like", "CD4+ TEM, Th2-like",
#                                    "CD4_act","CD4_ex","CD4_naive","CD4_reg"),
#                 c("Plasmacytoid DC","Pre-DC" ,"Unclassified DC","cDC1","cDC2",
#                                      "migDC","pDC", "preDC", "Myeloid cDC2","Myeloid cDC1" ),
#                 c("IFN-TAM","M-MDSC" ,"Monocyte, CD4+","Monocyte, CD4-","Monocyte, IFN signature",
#                                      "TIM"),
#                 c("Activated B cell", "B cell", "Class switched B cell", "Immature B cell", 
#                                     "Naive B cell","Plasma cell"),
#                 c("Neutrophil", "PMN-MDSC")
#                ))

# colnames(df) <- "conSense"

# df$major <- c(rep("cyto",length(c("CD8_SPP1_hi", "CD8+ Memory","CD8+ Effector", "CD8+ gd T cell","CD8+ Naive", "NK cell", "NK T cell", "NK cell","CD8_eff","CD8_ex", "NK"))),
#                     rep("helper",length(c("CD4+ Naive", "CD4+ T reg", "CD4+ TCM", "CD4+ TEM", 
#                                     "CD4+ TEM, Th1-like","CD4+ TEM, Th17-like", "CD4+ TEM, Th2-like",
#                                    "CD4_act","CD4_ex","CD4_naive","CD4_reg"))),
#                     rep("dc",length(c("Plasmacytoid DC","Pre-DC" ,"Unclassified DC","cDC1","cDC2",
#                                      "migDC","pDC", "preDC", "Myeloid cDC2","Myeloid cDC1"))),
#                     rep("mono",length(c("IFN-TAM","M-MDSC" ,"Monocyte, CD4+","Monocyte, CD4-","Monocyte, IFN signature",
#                                      "TIM"))),
#                     rep("bcell",length(c("Activated B cell", "B cell", "Class switched B cell", "Immature B cell", 
#                                     "Naive B cell","Plasma cell"))),
#                     rep("neut",length(c("Neutrophil", "PMN-MDSC")))
#                    )


# write.csv(df, "./cellTypez.csv")

#load in the cell type consensous classifications
seu.obj <- loadMeta(seu.obj = seu.obj, seu.file = NULL, metaFile = "./cellTypez.csv", groupBy = "conSense", metaAdd = "major", 
                     save = FALSE, outName = "", header = TRUE
                    )


### Fig extra - create raw UMAP
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "major",
              #cols = levels(seu.obj$dcColz),
              split.by = "major",
              pt.size = 0.5,
              label = F,
              ncol = 4,
              label.box = F,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi)
ggsave(paste("./output/", outName, "/", subName, "/", "rawUMAP2.png", sep = ""), width = 7, height = 7)

#clean up the dataset to remove cells that do not have a counterpart in both datasets
Idents(seu.obj) <- "major"
seu.obj <- subset(seu.obj, invert = T,
                          ident = c("Eosinophil", "DN T cell","gd T cell", "CD4+, IFN signature","Basophil", "CD34+ Unclassified", "T_IFN")
                         )
table(seu.obj$major)
dim(seu.obj)

#rename to match cell type across datasets
seu.obj$major <- droplevels(seu.obj$major)
seu.obj <- RenameIdents(seu.obj, c("Cycling T cell" = "T_cycling", "T_cycling" = "T_cycling")
                       )
seu.obj$major <- Idents(seu.obj)

#exclude any remaining cells that do not have an ID associated with the barcode
seu.obj$major <- as.factor(ifelse(is.na(seu.obj$major),"NA",as.character(seu.obj$major)))
Idents(seu.obj) <- "major"
seu.obj <- subset(seu.obj,invert = TRUE,
                  ident = "NA"
                 )
table(seu.obj$major)
dim(seu.obj)


### Fig extra - create raw UMAP split by major ID to ensure the above cade ran as desired
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "major",
              #cols = levels(seu.obj$dcColz),
              split.by = "major",
              pt.size = 0.5,
              label = F,
              ncol = 4,
              label.box = F,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi)
ggsave(paste("./output/", outName, "/", subName, "/", "rawUMAP3.png", sep = ""), width = 10, height = 10)


### Fig 1b - Create raw UMAP
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "major",
              cols = c("#D4F3A3","#C8C3E2","#0066A5","#AC0535","#148B7D","#F17D00","hotpink"),
              pt.size = 0.5,
              label = F,
              label.box = F,
              shuffle = TRUE
) + NoLegend()
pi <- formatUMAP(plot = pi, smallAxes = T) 
ggsave(paste("./output/", outName, "/", subName, "/majorID_UMAP_tiny.png", sep = ""), width = 7, height = 7)

#clean and reorder levels 
seu.obj$major <- droplevels(as.factor(seu.obj$major))
Idents(seu.obj) <- "major"
seu.obj <- RenameIdents(seu.obj, c("T_cycling" = "Cycling T cell", "bcell" = "B cell" ,   "cyto"  = "CD8 T cell"   , "dc"  = "Dendritic cell"  ,    "helper"  = "CD4 T cell" , "mono" = "Monocyte",    "neut" = "Neutrophil")
                       )
seu.obj$major <- Idents(seu.obj)

seu.obj$major <- factor(seu.obj$major, levels = levels(seu.obj$major)[c(5,3,1,2,6,4,7)])
seu.obj$majorID_sub <- seu.obj$major
seu.obj$major <- droplevels(as.factor(seu.obj$major))


### Fig 1c - key feature plots
features <- c("CD3G","CD8A","CD4",  "S100A12","DLA-DRA",
              "FLT3", "ANPEP", "MS4A1","JCHAIN","TOP2A"
             )
colorz <- "black"
fig1b <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol = 5, features = features, color = colorz, order = F) 
ggsave(paste("./output/", outName, "/", subName, "featPlots.png", sep = ""), width = 15, height = 6)


### Fig 1d - barchart comparing cell type proptions
p <- skewPlot(seu.obj, groupBy = "majorID_sub")
ggsave(paste("./output/", outName, "/", subName, "/", "barchart.png", sep = ""), width = 6, height = 4)

#store gene list associated with platelets (as presented in https://doi.org/10.5281/zenodo.7884518)
pal_feats = c('TIMP1', 'NAA10', 'ENSCAFG00000037735', 'GP6', 'SEC11C', 'FTL', 'NRGN', 'ACOT7', 'VCL', 'RSU1', 'ITGB1', 'H3-3A', 'RABGAP1L', 'SELP', 'SH3GLB1', 'ACTB', 'ENSCAFG00000008221', 'TLN1', 'GSN', 'AMD1', 'TREM2', 'SH3BGRL2', 'MYH9', 'PLEK', 'ENSCAFG00000042554', 'RAP1B', 'ENSCAFG00000004260', 'NAP1L1', 'PPBP', 'RASA3', 'ITGA2B', 'EIF1', 'ACTG1', 'C9H17orf64', 'JMJD6', 'CCL14', 'GNG11', 'IGF2BP3', 'TBXAS1', 'VDAC3', 'MARCHF2', 'TPM4', 'TKT', 'FTH1.1', 'FERMT3', 'RTN3', 'PRKAR2B', 'SVIP', 'ENSCAFG00000030286', 'ADA', 'MYL9', 'TUBB1', 'TUBA1B', 'METTL7A', 'THBS1', 'SERF2', 'PIF1', 'B2M', 'GAS2L1', 'YWHAH', 'HPSE', 'ATG3', 'ENSCAFG00000015217', 'ITGA6','RGS18', 'SUB1', 'LGALS1', 'CFL1', 'BIN2', 'CAT', 'RGS10', 'MGST3', 'TMBIM6', 'PFN1', 'CD63', 'RALBP1', 'GNAS', 'SEPTIN7', 'TPT1', 'UBB', 'ATF4', 'BBLN', 'MTDH', 'ENSCAFG00000017655','FYB1', 'ENO1', 'GABARAP', 'SSR4', 'MSN', 'ENSCAFG00000011134', 'ENSCAFG00000046637', 'COX8A', 'DLA-64', 'CD47', 'VASP', 'DYNLRB1', 'DLA88', 'SMDT1', 'ATP5PF','ELOB', 'ENSCAFG00000029155', 'ARPC3', 'VPS28', 'LRRFIP1', 'SRP14', 'ABRACL', 'ENSCAFG00000043577', 'ENSCAFG00000042598')


### Complete pseudobulk DGE by each major cell type
createPB(seu.obj = seu.obj, groupBy = "majorID_sub", comp = "cellSource", biologicalRep = "name", lowFilter = T, dwnSam = F, min.cell = 5,
         clusters = NULL, outDir = paste0("./output/", outName, "/", subName, "/pseudoBulk/"), grepTerm = "tils", grepLabel = c("TILs", "Blood")
)

pseudoDEG(metaPWD = paste0("./output/", outName, "/", subName, "/pseudoBulk/majorID_sub_deg_metaData.csv"),
          padj_cutoff = 0.05, lfcCut = 0.58, outDir = paste0("./output/", outName, "/", subName, "/pseudoBulk/"), 
          outName = outName, 
          idents.1_NAME = "TILs", idents.2_NAME = "Blood",
          inDir = paste0("./output/", outName, "/", subName, "/pseudoBulk/"), title = "All cells", 
          filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)

#load in the DGE results for each major cell type (exluding cycling cells) to plot and extract tissue gene signature 
pwds <- lapply(levels(seu.obj$majorID_sub)[c(1,2,4:7)], function(x){paste0("./output/", outName, "/", subName, "/pseudoBulk/", x, "/", outName, "_cluster_", x, "_all_genes.csv")})
df.list <- lapply(pwds, read.csv, header = T)

feats.list <- lapply(df.list, function(x){feats <- x %>% filter(log2FoldChange > 0) %>% select(gene)})
tumor.sig <- Reduce(intersect, list(feats.list[1][[1]]$gene, feats.list[2][[1]]$gene, feats.list[3][[1]]$gene, feats.list[4][[1]]$gene, feats.list[5][[1]]$gene, feats.list[5][[1]]$gene, feats.list[6][[1]]$gene))

#pathway GO
p <- plotGSEA(geneList = tumor.sig, category = "C5", species = "dog", termsTOplot = 10, upOnly = T, 
                     pvalueCutoff = 0.05, subcategory = NULL
                    ) + theme(axis.title=element_text(size = 16))

p <- p + scale_x_continuous(limits = c(-12,ceiling(max(p$data$x_axis)*1.05)), breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),name = "log10(p.adj)") + ggtitle("Gene ontology") + theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste0("./output/", outName, "/", subName, "/", "gseaPlot_1.png"), width =7, height = 7)


#pathway REACTOME
p <- plotGSEA(geneList = tumor.sig, category = "C2", species = "dog", termsTOplot = 10, upOnly = T,
                     pvalueCutoff = 0.05, subcategory = "CP:REACTOME"
                    ) + theme(axis.title=element_text(size = 16))

p <- p + scale_x_continuous(limits = c(-20,ceiling(max(p$data$x_axis)*1.05)), breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),name = "log10(p.adj)") + ggtitle("Reactome") + theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste0("./output/", outName, "/", subName, "/", "gseaPlot_2.png"), width =7, height = 7)


#also extract the conserved downregulated features
feats.list2 <- lapply(df.list, function(x){feats <- x %>% filter(log2FoldChange < 0) %>% select(gene)})
tumor.sig2 <- Reduce(intersect, list(feats.list2[1][[1]]$gene, feats.list2[2][[1]]$gene, feats.list2[3][[1]]$gene, feats.list2[4][[1]]$gene, feats.list2[5][[1]]$gene, feats.list2[5][[1]]$gene, feats.list2[6][[1]]$gene))

#stash the signature for later use
tumor.sig <- c(tumor.sig, tumor.sig2)


### Fig supp 1a - upset plot demonstrating the conserved gene signature
library(UpSetR)
upSet.df <- as.data.frame(unique(c(feats.list[1][[1]]$gene,feats.list[2][[1]]$gene,feats.list[3][[1]]$gene,feats.list[4][[1]]$gene,feats.list[5][[1]]$gene,feats.list[5][[1]]$gene,feats.list[6][[1]]$gene)))
colnames(upSet.df) <- "gene"
upSet.df[ ,levels(seu.obj$majorID_sub)[c(1,2,4:7)][1]] <- as.integer(ifelse(upSet.df$gene %in% feats.list[1][[1]]$gene, 1, 0))
upSet.df[ ,levels(seu.obj$majorID_sub)[c(1,2,4:7)][2]] <- as.integer(ifelse(upSet.df$gene %in% feats.list[2][[1]]$gene, 1, 0))
upSet.df[ ,levels(seu.obj$majorID_sub)[c(1,2,4:7)][3]] <- as.integer(ifelse(upSet.df$gene %in% feats.list[3][[1]]$gene, 1, 0))
upSet.df[ ,levels(seu.obj$majorID_sub)[c(1,2,4:7)][4]] <- as.integer(ifelse(upSet.df$gene %in% feats.list[4][[1]]$gene, 1, 0))
upSet.df[ ,levels(seu.obj$majorID_sub)[c(1,2,4:7)][5]] <- as.integer(ifelse(upSet.df$gene %in% feats.list[5][[1]]$gene, 1, 0))
upSet.df[ ,levels(seu.obj$majorID_sub)[c(1,2,4:7)][6]] <- as.integer(ifelse(upSet.df$gene %in% feats.list[6][[1]]$gene, 1, 0))

png(file = paste0("./output/", outName, "/", subName, "/",subName, "_upSet.png"), width=4000, height=2000, res=400)
par(mfcol=c(1,1))     
p <- upset(upSet.df, sets = colnames(upSet.df)[2:7], , cutoff = 10,  nintersects = 60)
p
dev.off()

### Complete pseudobulk DGE by each major cell type w/filter
createPB(seu.obj = seu.obj, groupBy = "majorID_sub", comp = "cellSource", biologicalRep = "name", min.cell = 5,
                     clusters = NULL, outDir = paste0("./output/", outName, "/", subName, "/pseudoBulk/"), grepTerm = "tils", grepLabel = c("TILs", "Blood"), featsTOexclude = c(pal_feats,tumor.sig), lowFilter = T, dwnSam =F
                    )

p_volc <- pseudoDEG(metaPWD = paste0("./output/", outName, "/", subName, "/pseudoBulk/majorID_sub_deg_metaData.csv"),
          padj_cutoff = 0.01, lfcCut = 0.58, outDir = paste0("./output/", outName, "/", subName, "/pseudoBulk/"), outName = paste0(subName, "_FILTERED"), idents.1_NAME = "TILs", idents.2_NAME = "Blood",
          inDir = paste0("./output/", outName, "/", subName, "/pseudoBulk/"), title = "TILS vs Blood", fromFile = T, meta = NULL, pbj = NULL, returnVolc = T, paired = F, pairBy = "", 
          minimalOuts = F, saveSigRes = T, filterTerm = "^ENSCAF", addLabs = NULL, mkDir = T
                     )


########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end all cells analysis   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


