#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")

### Analysis note: 

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin DC preprocessing   ######## <<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set output params
outName <- "dc"

#subset on the DCs & complete subset analysis
Idents(seu.obj) <- "conSense"
seu.obj.sub <- subset(seu.obj, idents = c(
    "Plasmacytoid DC","Pre-DC" ,"Unclassified DC","cDC1","cDC2",
    "migDC","pDC", "preDC", "Myeloid cDC2","Myeloid cDC1"
))
table(seu.obj.sub$conSense)
table(seu.obj.sub$name) > 50
seu.sub.list <- SplitObject(seu.obj.sub, split.by = "name")
seu.sub.list <- seu.sub.list[-c(1,2,3,8)] #remove samples < 50 cells
seu.obj <- indReClus(
    seu.list = seu.sub.list, seu.obj = NULL, 
    outDir = "./output/s2/", subName = "20230507_DC_bloodANDtils", 
    preSub = T, vars.to.regress = "percent.mt", k = 50
)
DefaultAssay(seu.obj) <- "integrated"
clusTree(
    seu.obj = seu.obj, 
    dout = "./output/clustree/", outName = "20230507_DC_bloodANDtils", 
    test_dims = c(40,35,30), algorithm = 3, prefix = "integrated_snn_res."
)
seu.obj <- dataVisUMAP(
    seu.obj = seu.obj, outDir = "./output/s3/", 
    outName = "20230507_DC_bloodANDtils", 
    final.dims = 30, final.res = 0.1, algorithm = 3,
    stashID = "clusterID_sub", prefix = "integrated_snn_res.", 
    min.dist = 0.6, n.neighbors = 75, 
    assay = "integrated", saveRDS = T
)

###################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin DC analysis   ######## <<<<<<<<<<<<<<
###################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#load in processed data
seu.obj <- readRDS(
    "../output/s3/20230507_DC_bloodANDtils_res0.1_dims30_dist0.6_neigh75_S3.rds"
)
seu.obj$cellSource <- ifelse(grepl("tils",seu.obj$name),"Tumor","Blood")
outName <- "dc"

#Check consenseous
seu.obj$ct <- ifelse(
    seu.obj$cellSource == "Tumor", 
    paste0("TIL;", seu.obj$celltype.l3), 
    paste0("Blood;",seu.obj$celltype.l3_pbmc)
)
#Plot % in blood
table(seu.obj$ct, seu.obj$clusterID_sub) %>% 
    melt() %>% separate(Var.1, sep = ";", c("source", "ct")) %>% 
    filter(source == "Blood") %>%
    group_by(Var.2) %>% 
    mutate(pct = value/sum(value)) %>% 
    filter(value == max(value))
#Plot % in TIL
table(seu.obj$ct, seu.obj$clusterID_sub) %>% 
    melt() %>% 
    separate(Var.1, sep = ";", c("source", "ct")) %>% 
    filter(source == "TIL") %>% 
    group_by(Var.2) %>% 
    mutate(pct = value/sum(value)) %>% 
    filter(value == max(value))
#Plot % in both
Idents(seu.obj) <- "cellSource"
seu.obj.ds <- subset(seu.obj, downsample = min(table(seu.obj$cellSource)))
table(seu.obj.ds$ct, seu.obj.ds$clusterID_sub) %>% 
    sweep(., 2, colSums(.), `/`) %>% 
    melt() %>% separate(Var.1, sep = ";", c("source", "ct")) %>% 
    group_by(source, Var.2) %>% 
    filter(value == max(value))

#Set metadata
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c(
    "0" = "#6D0026", "1" = "#A0060A", "2" = "#C12000", "3" = "#AC0535",
    "4" = "#7D0025"
))
seu.obj$sub_colz <- Idents(seu.obj)
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c(
    "0" = "cDC2 (c0)", "1" = "cDC1 (c1)", "2" = "pDC (c2)", "3" = "mregDC (c3)",
    "4" = "preDC (c4)"
))
seu.obj$majorID_sub <- Idents(seu.obj)
seu.obj$majorID_sub <- factor(seu.obj$majorID_sub, levels = c(
    "cDC2 (c0)","cDC1 (c1)", "pDC (c2)","mregDC (c3)", "preDC (c4)"
)[c(5,2,1,4,3)])

#Export annotations
write.csv(seu.obj@meta.data["majorID_sub"], file = paste0("../output/annotations/", outName, ".csv"))

### Supp data - Generate violin plots for each cluster
vilnPlots(
    seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, 
    outDir = paste0("../output/viln/", outName, "/"), 
    outName = outName, returnViln = F,
    outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
)

### Fig sup: Create violin plots for key feats
features = c(
    "IL3RA", "PGLYRP2", "DDR2", 
    "CADM1","DNASE1L3", "CLEC1B", 
    "CD1C","PID1","CD300H",
    "IL21R","IL4I1", "IDO1",
    "IGF1","RARRES2", "IGHM","IGKC"
)

pi <- VlnPlot(
    object = seu.obj,
    pt.size = 0,
    same.y.lims = F,
    group.by = "majorID_sub",
    combine = T,
    cols = levels(seu.obj$sub_colz),
    stack = T,
    fill.by = "ident",
    flip = T,
    features = rev(features)
) + 
    NoLegend() + 
    theme(
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank()
    )

ggsave(paste0("../output/", outName, "/", outName, "_selectViln.png"), width = 5, height = 6)

### Fig 5a - UMAP by clusterID_sub
pi <- DimPlot(
    seu.obj, 
    reduction = "umap", 
    group.by = "clusterID_sub",
    cols = levels(seu.obj$sub_colz),
    pt.size = 0.5,
    label = TRUE,
    label.box = TRUE,
    shuffle = TRUE
)
pi <- cusLabels(
    plot = pi, shape = 21, size = 10, textSize = 6, 
    alpha = 0.8, labCol = "white"
) + 
    NoLegend() + 
    theme(
        axis.title = element_blank(),
        panel.border = element_blank()
    )
ggsave(paste("../output/", outName, "/", "rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig 5b - skew plot for abundance analysis
p <- skewPlot(
    seu.obj, groupBy = "majorID_sub", yAxisLabel = "Percent of DCs",
    dout = paste0("../output/", outName), outName = outName, 
    sampleRep = "name", grepTerm = "tils", grepRes = c("Tumor","Blood")
)
ggsave(paste0("../output/", outName, "/", "skewPlot.png"), width = 6, height = 4)


# ### Fig 4b - skew plot for abundance analysis
# metadata <- seu.obj@meta.data
# expression <- as.data.frame(t(seu.obj@assays$RNA@data)) #use log noralized count
# expression$majorID_sub <- seu.obj@meta.data[rownames(expression), ]$majorID_sub
# clusAvg_expression <- expression %>% 
#     group_by(majorID_sub) %>% 
#     summarise(across(where(is.numeric), mean)) %>% 
#     column_to_rownames(var = "majorID_sub")

# #complete herichcial cluster and spin clades to get in visually pleasing presnetion
# M <- cor(t(clusAvg_expression), method = "pearson")

# hc <- hclust(as.dist(M), method = "complete")

# M <- M[rev(rownames(M)), ]
# M[row(M) + col(M) > nrow(M) + 1] <- NA
# melted_cormat <- reshape2::melt(M, na.rm = TRUE)
# #Plot heatmap
# ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
#     geom_tile(color = "white")+
#     scale_fill_gradient(
#         low = "white", high = "red", limit = c(0.85, 1), space = "Lab", 
#         name = "Pearson\ncorrelation"
#     ) +
#     theme_minimal() + 
#     coord_fixed() + 
#     geom_text(aes(Var2, Var1, label = round(value, 2)), color = "black", size = 4) +
#     theme(
#         axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
#         axis.text.y = element_text(size = 12),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         axis.ticks = element_blank()
#     )
# ggsave(paste0("../output/", outName, "/", "coor.png"), width = 7, height = 7)


### Fig 5c - DGE analysis

#load in the tumor and pal signatures to exlude from DE analysis
pal_feats = c(
    "TIMP1", "NAA10", "ENSCAFG00000037735", "GP6", "SEC11C", "FTL", 
    "NRGN", "ACOT7", "VCL", "RSU1", "ITGB1", "H3-3A", "RABGAP1L", 
    "SELP", "SH3GLB1", "ACTB", "ENSCAFG00000008221", "TLN1", "GSN", 
    "AMD1", "TREM2", "SH3BGRL2", "MYH9", "PLEK", "ENSCAFG00000042554", 
    "RAP1B", "ENSCAFG00000004260", "NAP1L1", "PPBP", "RASA3", "ITGA2B", 
    "EIF1", "ACTG1", "C9H17orf64", "JMJD6", "CCL14", "GNG11", "IGF2BP3", 
    "TBXAS1", "VDAC3", "MARCHF2", "TPM4", "TKT", "FTH1.1", "FERMT3", 
    "RTN3", "PRKAR2B", "SVIP", "ENSCAFG00000030286", "ADA", "MYL9", 
    "TUBB1", "TUBA1B", "METTL7A", "THBS1", "SERF2", "PIF1", "B2M", 
    "GAS2L1", "YWHAH", "HPSE", "ATG3", "ENSCAFG00000015217", "ITGA6", 
    "RGS18", "SUB1", "LGALS1", "CFL1", "BIN2", "CAT", "RGS10", "MGST3", 
    "TMBIM6", "PFN1", "CD63", "RALBP1", "GNAS", "SEPTIN7", "TPT1", 
    "UBB", "ATF4", "BBLN", "MTDH", "ENSCAFG00000017655", "FYB1", 
    "ENO1", "GABARAP", "SSR4", "MSN", "ENSCAFG00000011134", "ENSCAFG00000046637", 
    "COX8A", "DLA-64", "CD47", "VASP", "DYNLRB1", "DLA88", "SMDT1", 
    "ATP5PF", "ELOB", "ENSCAFG00000029155", "ARPC3", "VPS28", "LRRFIP1", 
    "SRP14", "ABRACL", "ENSCAFG00000043577", "ENSCAFG00000042598"
)
tumor.sig <- read.csv("./metaData/tumorSig.csv", header = T)$x

createPB(
    seu.obj = seu.obj, groupBy = "allCells", 
    comp = "cellSource", biologicalRep = "name",
    outDir = paste0("../output/", outName, "/pseudoBulk/"), 
    grepTerm = "tils", grepLabel = c("Tumor", "Blood"), 
    featsTOexclude = c(pal_feats,tumor.sig), 
    lowFilter = T, dwnSam = F
)
p_volc <- pseudoDEG(
    inDir = paste0("../output/", outName, "/pseudoBulk/"), 
    metaPWD = paste0("../output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"),
    outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = outName,
    padj_cutoff = 0.01, lfcCut = 0.58, strict_lfc = T,
    idents.1_NAME = "Tumor", idents.2_NAME = "Blood", title = "TILS vs Blood", 
    fromFile = T, returnVolc = T, filterTerm = "^ENSCAF", mkDir = T
)
p <- prettyVolc(
    plot = p_volc[[1]], rightLab = "Tumor", leftLab = "Blood", 
    rightCol = "red", leftCol = "blue", arrowz = T
) + 
    theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_blank(),
        legend.position = c(0.10, 0.85),
    )
ggsave(paste0("../output/", outName, "/", "volcPlot.png"), width = 7, height = 7)


### Fig 5d - GO GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", 
                           outName, "_cluster_allCells_all_genes.csv"),
    geneList = NULL, category = "C5", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = F, pvalueCutoff = 0.05, subcategory = NULL,
    lolli = T,
    saveRes = paste0("../output/", outName, "/c5_", outName, "_res.csv")
) + 
    theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)
    )
p <- p + scale_x_continuous(
    limits = c(-25,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Gene ontology") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_1.png"), width = 7, height = 7)

### Fig 5e - Reactome GSEA of DEGs

p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", 
                           outName, "_cluster_allCells_all_genes.csv"),
    geneList = NULL, category = "C2", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "CP:REACTOME", 
    lolli = T,
    saveRes = paste0("../output/", outName, "/c2_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-20,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Reactome") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_2.png"), width = 7, height = 7)

### Fig Sup - ImmuneSigDB GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", 
                           outName, "_cluster_allCells_all_genes.csv"),
    geneList = NULL, category = "C7", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
    lolli = T, filterTerm = "_DC_|BMDC",
    saveRes = paste0("../output/", outName, "/c7_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-60,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("ImmuneSigDB") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 0.9)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_3.png"), width = 7, height = 7)

### Fig 5f - Split UMAP of selected DEGs
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$cellSource)))

features <- c("IL16", "IRF4", "CXCR4", "IL18", "IL22RA2", "IL4I1", "CCR7", "CD274", "IDO1")

p <- FeaturePlot(
    seu.obj.sub, features = features, pt.size = 0.1, split.by = "cellSource", 
    order = T, cols = c("lightgrey","darkblue"), by.col = F
) + 
    labs(x = "UMAP1", y = "UMAP2") & 
    theme(
        axis.text = element_blank(),
        axis.title.y.right = element_text(size = 11),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(size=11),
        title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) 
ggsave(paste("../output/", outName, "/", "splitFeats.png", sep = ""), width = 12, height = 4)


### Fig sup - Split dot plot of select DEGs
res.df <- read.csv(paste0("../output/", outName, "/pseudoBulk/allCells/", 
                          outName, "_cluster_allCells_all_genes.csv"))
res.df <- res.df[!grepl("^ENS", res.df$gene), ]
geneList_UP <- res.df %>% 
    filter(padj < 0.1) %>% 
    filter(log2FoldChange > 1) %>% 
    pull(gene)
geneList_DWN <- res.df %>% 
    filter(padj < 0.1) %>% 
    filter(log2FoldChange < -1) %>% 
    pull(gene)

seu.obj$cellSource <- as.factor(seu.obj$cellSource)
p <- splitDot(
    seu.obj = seu.obj, groupBy = "majorID_sub", splitBy = "cellSource", buffer = 150,
    namedColz = setNames(c("#F8766D", "#00BFC4"),  c("Blood", "Tumor")), 
    geneList_UP = unique(c(geneList_UP[1:20], "CXCR4","VEGFA", "IL1A","IL1B",
                           "CCR7", "CD274","IDO1")), 
    geneList_DWN = unique(c(geneList_DWN[1:20], "IL16", "IRF4")), 
    geneColz = c("red", "blue")
)
p <- p +
    theme(
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.justification = 'center'
    )
ggsave(plot = p, paste0("../output/", outName, "/", outName, "_splitDot.png"), width = 13, height = 5)


## too few DCs in TME to do within cluster de analysis

#rename metadata
seu.obj$celltype <- seu.obj$majorID_sub
seu.obj$clusterID <- seu.obj$clusterID_sub

# Export data for UCSC cell browser
ExportToCB_cus(
    seu.obj = seu.obj, dataset.name = "DC", outDir = "../output/cb_input/", 
    markers = "../output/supplementalData/supplemental_data_6.csv", 
    metaSlots = c(
        "cluster","gene","avg_log2FC","p_val_adj"
    ),
    reduction = "umap",  
    colsTOkeep = c(
        "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase",
        "name", "clusterID", "celltype"
    ),
    skipEXPR = F, test = F,
    feats = c(
        "IL3RA", "PGLYRP2", "DDR2", "CADM1","DNASE1L3", "CLEC1B", 
        "CD1C","PID1","CD300H", "IL21R","IL4I1", "IDO1",
        "IGF1","RARRES2", "IGHM","IGKC"
    )
)

#################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end DC analysis   ######## <<<<<<<<<<<<<<
#################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

