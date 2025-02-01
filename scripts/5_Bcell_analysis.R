#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")

### Analysis note: 

############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin B cell preprocessing   ######## <<<<<<<<<<<<<<
############################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#set output params
outName <- "bcell"

#subset on the B cells & complete subset analysis
table(seu.obj$conSense)[grepl("B",names(table(seu.obj$conSense)))]
Idents(seu.obj) <- "conSense"
seu.obj.sub <- subset(seu.obj, ident = c(
    "Activated B cell", "B cell", "Class switched B cell", "Immature B cell", 
    "Naive B cell","Plasma cell"
))

table(seu.obj.sub$conSense)
min(table(seu.obj.sub$name)) > 100
seu.sub.list <- SplitObject(seu.obj.sub, split.by = "name")
#Remove samples that have too few cells to be included
seu.sub.list <- seu.sub.list[-c(3,5,12,14,16)]
seu.obj <- indReClus(
    seu.list = seu.sub.list, seu.obj = NULL, outDir = "./output/s2/", 
    subName = "20230507_bcell_bloodANDtils", preSub = T, 
    vars.to.regress = "percent.mt", nfeatures = 2000
)
clusTree(
    seu.obj = seu.obj, dout = "./output/clustree/", 
    outName = "20230507_bcell_bloodANDtils", 
    test_dims = c(40,35,30), algorithm = 3, prefix = "integrated_snn_res."
)
seu.obj <- dataVisUMAP(
    seu.obj = seu.obj, outDir = "./output/s3/", 
    outName = "20230507_bcell_bloodANDtils", 
    final.dims = 30, final.res = 0.4, algorithm = 3,
    min.dist = 0.5, n.neighbors = 60, 
    stashID = "clusterID_sub", prefix = "integrated_snn_res.", 
    assay = "integrated", saveRDS = T
)

########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin B cell analysis   ######## <<<<<<<<<<<<<<
########################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#Load in processed data
seu.obj <- readRDS(
    "../output/s3/20230507_bcell_bloodANDtils_res0.4_dims30_dist0.5_neigh60_S3.rds"
)
seu.obj$cellSource <- ifelse(grepl("tils",seu.obj$name),"Tumor","Blood")
outName <- "bcell"

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
colz.df <- read.csv("./metaData/majorGroups.csv")
colz.df <- colz.df[colz.df$majorID2 == "b" | colz.df$majorID2 == "plasma", ]
table(seu.obj$clusterID_sub, seu.obj$conSense)
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c(
    "0" = "#F8766D", "1" = "#A3A500", "2" = "#00BF7D", "3" = "#00B0F6",
    "4" = "#E76BF3"
))
seu.obj$sub_colz <- Idents(seu.obj)
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c(
    "0" = "Naive (c0)", "1" = "Class_switched (c1)", "2" = "Plasma_cell (c2)", 
    "3" = "Immature (c3)", "4" = "Plasma_cell (c2)"
))
seu.obj$majorID_sub <- Idents(seu.obj)
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c(
    "0" = "0", "1" = "1", "2" = "2", "3" = "3", "4" = "2"
))
seu.obj$clusID_new <- Idents(seu.obj)

#Export annotations
write.csv(seu.obj@meta.data["majorID_sub"], file = paste0("../output/annotations/", outName, ".csv"))

### Supp data - Generate violin plots for each cluster
vilnPlots(
    seu.obj = seu.obj, groupBy = "majorID_sub", numOfFeats = 24, 
    outDir = paste0("../output/viln/", outName, "/"), 
    outName = outName, returnViln = F,
    outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS")
)

modulez <- list(
  `Immature` = c("SYT1", "PAX5", "VPREB3", "ERC2", "TMTC2", "KLHL14", "F8", 
                 "TEX9", "TDRP", "ADGRF1"),
  `Naive` = c("TNFRSF13C", "BANK1", "HTR1F", "PAX5", "EBF1", "BTLA", "NRIP1", 
              "ADAM9"),
  `Class switched` = c("TNFRSF13C", "GOLM1", "BANK1", "BTLA", "EBF1", 
                       "DYNC1I1", "MTMR2", "PAX5"),
  `Activated` = c("IGKC", "CACNB2", "PAX5", "TNFRSF13C", "IGHM", "RASGRF2", 
                  "AOX2", "BCAR3", "ADAM32"),
  `Plasma` = c("JCHAIN", "MZB1", "TXNDC5", "LMAN1", "FKBP11", "LAP3", "DERL3", 
               "CCR10", "MKI67", "TNFRSF13B")
)
#run module score
seu.obj <- AddModuleScore(
    seu.obj, features = modulez, name = "_score"
)
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)
#plot the results of enrichment scores
features <- names(modulez)
ecScores <- majorDot(
    seu.obj = seu.obj, groupBy = "clusID_new", features = rev(features)
) + 
    coord_flip() + 
    theme(
        plot.margin = margin(3, 0, 3, 0, "pt"),
        axis.text.y=element_text(size=10),
        axis.title = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        axis.text.x = element_text(angle=0, hjust = 0.5)
    ) + 
    scale_y_discrete(position = "right") + 
    scale_colour_continuous(name="Module score", type = "viridis")
ggsave(paste("../output/", outName, "/", outName, "_dotPlot_ecScores.png", sep = ""), width = 6, height = 4)

#plot indivdual members of each term
modulez <- c(list("Module score" = names(modulez)), modulez)
labelz <- as.data.frame(names(modulez))
colnames(labelz) <- "labz"
labelz$modLen <- unname(unlist(lapply(modulez, length)))
cntr <- 0
plots <- lapply(modulez, function(x){
    cntr <<- cntr+1
    labz.df <- labelz[cntr,]

    majorDot(
        seu.obj = seu.obj, groupBy = "clusID_new", features = rev(unname(unlist(x)))
    ) + 
        theme(
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "right",
            legend.direction = "vertical",
            axis.title = element_blank(),
            plot.margin = margin(3, 0, 3, 0, "pt")
        ) + 
        scale_colour_distiller(palette = "RdYlBu", name='Average\nexpression', limits = c(-2.5,2.5)) + 
        geom_text(
            data = labz.df, aes(label = labz, y = 4.85, x = (modLen+1)/2),
            angle = 270, vjust = 0.5, hjust=0.5, size = 12*0.36
        ) + 
        coord_flip(ylim = c(1, 4.75), clip = "off") + 
        annotate(
            "segment", x = -Inf, y = 4.5, xend = Inf, yend = 4.5, 
            lineend = "round", linejoin = "bevel", linetype ="solid", 
            colour = "grey70", alpha = 0.7, size = 0.5
        )
})
#get all the plots together
    patch <- area()
    nrow <- length(modulez)
    ncol <- 1
    counter=0
    for (i in 1:nrow) {
        for (x in 1:ncol) {
                patch <- append(patch, area(t = i, l = x, b = i, r = x))
        }
    }
#change to the color of the module scores for visual distinction & plot final
plots$`Module score` <- plots$`Module score` + 
    theme(axis.text.x = element_text(angle=0, hjust = 0.5)) + 
    scale_y_discrete(position = "right") + 
    scale_colour_viridis() + 
    guides(color = guide_colorbar(title = 'Module\nscore'), limits = c(-2.5,2.5))
p <- Reduce( `+`, plots ) +
    plot_layout(
        guides = "collect", design = patch, 
        height = unname(unlist(lapply(modulez, length)))/sum(unname(unlist(lapply(modulez, length))))
    )
ggsave(paste("../output/", outName, "/", outName, "_dotPlot_ecScores_2.png", sep = ""), width = 2, height = 7, scale = 2)




### Fig 4a - UMAP by clusterID_sub
pi <- DimPlot(
    seu.obj, 
    reduction = "umap", 
    group.by = "clusID_new",
    cols = colz.df$colour[c(1,4,5,2)],
    pt.size = 0.5,
    label = TRUE,
    label.box = TRUE,
    shuffle = TRUE
)
pi <- cusLabels(
    plot = pi, shape = 21, size = 10, textSize = 6, alpha = 0.8, 
    labCol = colz.df$labCol[c(1,4,5,2)]
) + 
    NoLegend() + 
    theme(
        axis.title = element_blank(),
        panel.border = element_blank()
    )
ggsave(paste("../output/", outName, "/", "rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig 4b - skew plot for abundance analysis
p <- skewPlot(
    seu.obj, groupBy = "majorID_sub", yAxisLabel = "Percent of B cells",
    dout = paste0("../output/", outName), outName = outName, 
    sampleRep = "name", grepTerm = "tils", grepRes = c("Tumor","Blood")
)
ggsave(paste0("../output/", outName, "/", "skewPlot.png"), width = 6, height = 4)

### Fig supp - calc sample distances to support seperation of B and plasam cells
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@data))
expression$majorID_sub <- seu.obj@meta.data[rownames(expression), ]$majorID_sub
clusAvg_expression <- expression %>% 
    group_by(majorID_sub) %>% 
    summarise(across(where(is.numeric), mean)) %>% 
    column_to_rownames(var = "majorID_sub")

M <- as.matrix(dist(clusAvg_expression, method = "euclidean"))
M <- M[rev(rownames(M)), ]
M[row(M) + col(M) > nrow(M) + 1] <- NA
melted_cormat <- reshape2::melt(M, na.rm = TRUE)
#Plot heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "grey")+
    scale_fill_gradient(
        low = "white", high = "red", limit = c(0, 25), space = "Lab", 
        name = "Euclidean\ndistance"
    ) +
    theme_minimal() + 
    coord_fixed() + 
    geom_text(aes(Var2, Var1, label = round(value, 2)), color = "black", size = 4) +
    theme(
        axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()
    ) +
    NoLegend() +
    ggtitle("Euclidean distance")
ggsave(paste0("../output/", outName, "/", "coor.png"), width = 5, height = 5)

### Fig 4c - DGE analysis

#exclude plasma cells from DGE analysis
seu.obj.sub <- subset(
    seu.obj, invert = T, subset = majorID_sub == "Plasma_cell (c2)"
)

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
    seu.obj = seu.obj.sub, groupBy = "allCells", 
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

### Fig 4d - GO GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", 
                           outName, "_cluster_allCells_all_genes.csv"),
    geneList = NULL, category = "C5", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = NULL,
    lolli = T,
    saveRes = paste0("../output/", outName, "/c5_", outName, "_res.csv")
) + 
    theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)
    )
p <- p + scale_x_continuous(
    limits = c(-50,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Gene ontology") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.07),
        legend.background = element_rect(fill = 'transparent', colour = NA),
        legend.key = element_rect(fill = 'transparent', colour = NA)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_1.png"), width = 7, height = 7)

### Fig 4e - Reactome GSEA of DEGs
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
    limits = c(-35,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Reactome") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_2.png"), width = 7, height = 7)

### Fig extra - ImmuneSigDB GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", 
                           outName, "_cluster_allCells_all_genes.csv"),
    geneList = NULL, category = "C7", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
    lolli = T, filterTerm = "BCELL|PLASMA_CELL|B_CELL",
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
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_3.png"), width = 7, height = 7)

### Fig extra - Cell type (C8) GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/pseudoBulk/allCells/", 
                           outName, "_cluster_allCells_all_genes.csv"),
    geneList = NULL, category = "C8", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, 
    lolli = T,
    saveRes = paste0("../output/", outName, "/c8_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-85,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Cell type") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_4.png"), width = 7, height = 7)


### Fig 2f - Split UMAP of selected DEGs
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$cellSource)))
features <- c("LYZ", "BTLA", "ADD3", "IGHM", "FOS", "FOSB", "FAS", "IGFLR1", "EGFL7")
p <- FeaturePlot(
    seu.obj.sub, features = features, pt.size = 0.1, split.by = "cellSource", 
    order = T, cols = c("lightgrey", "darkblue"), by.col = F
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


res.df <- read.csv(file = paste0("../output/", outName, "/pseudoBulk/allCells/",
                                 outName, "_cluster_allCells_all_genes.csv"))
res.df <- res.df[!grepl("^ENS", res.df$gene), ]
geneList_UP <- res.df[res.df$padj < 0.1 & res.df$log2FoldChange > 1, ]$gene
geneList_DWN <- res.df[res.df$padj < 0.1 & res.df$log2FoldChange < -1, ]$gene
seu.obj$cellSource <- as.factor(seu.obj$cellSource)
p <- splitDot(
    seu.obj = seu.obj, groupBy = "majorID_sub", splitBy = "cellSource", buffer = 150,
    namedColz = setNames(c("#F8766D", "#00BFC4"),  c("Blood", "Tumor")), 
    geneList_UP = unique(c(geneList_UP[1:20], "FOS", "FOSB", "FAS", "IGFLR1", "EGFL7")), 
    geneList_DWN = unique(c(geneList_DWN[1:20], "LYZ", "BTLA", "ADD3", "IGHM")), 
    geneColz = c("red", "blue")
)
p <- p +
    theme(
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.justification = 'center'
    )
ggsave(plot = p, paste0("../output/", outName, "/", outName, "_splitDot.png"), width = 12, height = 5)


## too few B cells in TME to do within cluster de analysis
## Repeat analysis on only plasma cells
outName <- "plasma"
seu.obj.sub <- subset(
    seu.obj, invert = F, subset = majorID_sub == "Plasma_cell (c2)"
)

createPB(
    seu.obj = seu.obj.sub, groupBy = "allCells", 
    comp = "cellSource", biologicalRep = "name",
    outDir = paste0("../output/", outName, "/pseudoBulk/"), 
    grepTerm = "tils", grepLabel = c("TILs", "Blood"), 
    featsTOexclude = c(pal_feats,tumor.sig), 
    lowFilter = T, dwnSam = F
)
p_volc <- pseudoDEG(
    inDir = paste0("../output/", outName, "/pseudoBulk/"), 
    metaPWD = paste0("../output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"),
    outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = outName,
    padj_cutoff = 0.01, lfcCut = 0.58, strict_lfc = T,
    idents.1_NAME = "TILs", idents.2_NAME = "Blood", title = "TILS vs Blood", 
    fromFile = T, returnVolc = T, filterTerm = "^ENSCAF", mkDir = T
)
p <- prettyVolc(
    plot = p_volc[[1]], rightLab = NULL, leftLab = NULL, 
    rightCol = "red", leftCol = "blue", arrowz = F
) + 
    theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_blank()
    )
ggsave(paste0("../output/", outName, "/", "volcPlot.png"), width = 7, height = 7)

### Fig 2d - GO GSEA of DEGs
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
    limits = c(-80,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Gene ontology") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_5.png"), width = 7, height = 7)

### Fig 2e - Reactome GSEA of DEGs -- no enrichment

### Fig 2f - Split UMAP of selected DEGs
Idents(seu.obj) <- "cellSource"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$cellSource)))
features <- c("LYZ", "BTLA", "ADD3", "IGHM", "FOS", "FOSB", "FAS", "IGFLR1", "EGFL7")
p <- FeaturePlot(
    seu.obj.sub, features = features, pt.size = 0.1, split.by = "cellSource", 
    order = T, cols = c("lightgrey", "darkblue"), by.col = F
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

#plasma vs other B cell
res <- btwnClusDEG(
  seu.obj = seu.obj, groupBy = "majorID_sub", 
  idents.1 = "Plasma_cell (c2)", 
  bioRep = "name", padj_cutoff = 0.05, lfcCut = 1, topn=c(20,20),
  minCells = 5, outDir = paste0("../output/", outName, "/"), 
  title = "Wald",
  idents.1_NAME = "Plasma", idents.2_NAME = "Other B cell", 
  strict_lfc = T,
  returnVolc = F, doLinDEG = F, 
  paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 12, 
  dwnCol = "blue", stblCol = "grey",upCol = "red", labSize = 3
)

### Fig extra - ImmuneSigDB GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/", 
                           "Plasma_vs_Other_B_cell_all_genes.csv"),
    geneList = NULL, category = "C7", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
    lolli = T, filterTerm = "BCELL|PLASMA_CELL|B_CELL",
    saveRes = paste0("../output/", outName, "/c7_", outName, "_plasmaVother.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-90,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("ImmuneSigDB") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_6.png"), width = 7, height = 7)

### Fig extra - Cell type (C8) GSEA of DEGs
p <- plotGSEA(
    pwdTOgeneList = paste0("../output/", outName, "/", 
                           "Plasma_vs_Other_B_cell_all_genes.csv"),
    geneList = NULL, category = "C8", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05,
    lolli = T,
    saveRes = paste0("../output/", outName, "/c7_", outName, "_plasmaVother.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-60,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "-log10(p.adj)"
) + 
    ggtitle("Cell type") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_7.png"), width = 7, height = 7)

#rename metadata
seu.obj$celltype <- seu.obj$majorID_sub
seu.obj$clusterID <- seu.obj$clusID_new

# Export data for UCSC cell browser
ExportToCB_cus(
    seu.obj = seu.obj, dataset.name = "Bcell", outDir = "../output/cb_input/", 
    markers = "../output/supplementalData/supplemental_data_5.csv", 
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
        "ERC2", "PAX5", "BANK1", "CCR10", "BTLA", "MZB1", "DYNC1I1", 
        "DERL3", "VPREB3", "MKI67", "EBF1", "BTLA", "FKBP11", "LAP3", 
        "RASGRF2", "LMAN1", "TNFRSF13B", "PAX5", "JCHAIN", "GOLM1"
    )
)

######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end B cell analysis   ######## <<<<<<<<<<<<<<
######################################## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


