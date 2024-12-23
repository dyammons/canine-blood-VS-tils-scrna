#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")
outName <- "hu_blood"
# load10x(
#     din = "../external_data/blood/", dout = "../output/s1_hu_blood/", readCnts = T,
#     outName = outName, testQC = F, nFeature_RNA_high = 5500, nFeature_RNA_low = 200, 
#     nCount_RNA_high = 100000, nCount_RNA_low = 100, percent.mt_high = 15, 
#     mt_pattern = "^MT-", removeDubs = FALSE, 
#     removeRBC_pal = FALSE, pal_feats = NULL, isolatePalRBC = FALSE
# )

files <- list.files(
    path = "../output/s1_hu_blood", pattern = ".rds", all.files = FALSE,
    full.names = TRUE
)
seu.obj.list <- mapply(readRDS, files)
names(seu.obj.list) <- paste0("H", 2:7)
blood_feats <- unique(unlist(lapply(seu.obj.list, rownames)))

seu.tumor <- readRDS("../output/s3/human_naive_n6_res0.8_dims45_dist0.3_neigh30_S3.rds")
seu.tumor <- loadMeta(
    seu.obj = seu.tumor, metaFile = "./metaData/all_celltype_human.csv", 
    groupBy = "clusterID", metaAdd = "celltype.l1"
)
seu.tumor <- loadMeta(
    seu.obj = seu.tumor, metaFile = "./metaData/all_celltype_human.csv", 
    groupBy = "clusterID", metaAdd = "celltype.l2"
)

Idents(seu.tumor) <- "orig.ident"
seu.tumor <- RenameIdents(seu.tumor, c(
    "OS_1" = "GSM4952363", "OS_2" = "GSM4952364", 
    "OS_3" = "GSM4952365", "OS_4" = "GSM5155198",
    "OS_5" = "GSM5155199", "OS_6" = "GSM5155200"
))
seu.tumor$GSM_ID <- Idents(seu.tumor)

write.csv(seu.tumor@meta.data[c("GSM_ID", "celltype.l1", "celltype.l2")], file = "./metaData/human_os_anno.csv")

seu.tumor <- subset(seu.tumor, subset = celltype.l2 %in% c(
    "CD8 T cell", "CD4 T cell", "CD14_TIM", "Plasma cell", "IFN-TAM", "NK", 
    "cDC2", "B cell", "T_IFN", "pDC"
))
tumor_feats <- rownames(seu.tumor)
seu.tumor.list <- SplitObject(seu.tumor, split.by = "orig.ident")

shared_feats <- intersect(blood_feats, tumor_feats)
seu.list <- lapply(c(seu.obj.list, seu.tumor.list), function(seu){
    cnts <- seu@assays$RNA@counts
    cnts <- cnts[rownames(cnts) %in% shared_feats, ]
    seu.obj <- CreateSeuratObject(
        cnts, project = "human_analysis", assay = "RNA",
        min.cells = 0, min.features = 0, names.field = 1,
        names.delim = "_", meta.data = seu@meta.data
    )
})

seu.obj <- indReClus(
    seu.obj = NULL, outDir = "../output/s2_human/", subName = "hu_bloodANDtils",
    preSub = T, seu.list = seu.list, nfeatures = 2500, 
    vars.to.regress = "percent.mt"
)
seu.obj <- dataVisUMAP(
    seu.obj = seu.obj, outDir = "../output/s3/", 
    outName = "20241208_HU_bloodANDtils_2500", final.dims = 40, final.res = 0.6,
    stashID = "clusterID_new", algorithm = 3, prefix = "integrated_snn_res.", 
    min.dist = 0.6, n.neighbors = 75, assay = "integrated", saveRDS = T,
)
#Load in the integrated data
seu.obj <- readRDS("../output/s3/20241208_HU_bloodANDtils_2500_res0.6_dims40_dist0.6_neigh75_S3.rds")
seu.obj$cellSource <- ifelse(grepl("H", seu.obj$orig.ident), "blood", "tumor")
outName <- "human"
singleR(
    seu.obj, outName = "celltypeID", clusters = "clusterID_new", 
    outDir = paste0("../output/", outName, "/")
)

#Identify genes to use for comparative analysis
#load in the tumor and pal signatures to exlude from DE analysis
genez <- read.csv("./metaData/canine_human_convert.csv")
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

genez <- genez[! genez$input_gene %in% c(pal_feats,tumor.sig), ]



# B cells
seu.obj.sub <- subset(seu.obj, subset = clusterID_new == 5)
seu.obj.sub$allCells <- "allCells"
featsTOexclude <- rownames(seu.obj.sub)[! rownames(seu.obj.sub) %in% genez$ortholog_gene]

createPB(
    seu.obj = seu.obj.sub, groupBy = "allCells", min.cell = 20,
    comp = "cellSource", biologicalRep = "orig.ident",
    outDir = paste0("../output/", outName, "/pseudoBulk/"), 
    grepTerm = "OS", grepLabel = c("TILs", "Blood"), 
    featsTOexclude = featsTOexclude, 
    lowFilter = T, dwnSam = F
)
p_volc <- pseudoDEG(
    inDir = paste0("../output/", outName, "/pseudoBulk/"), 
    metaPWD = paste0("../output/", outName, "/pseudoBulk/allCells_deg_metaData.csv"),
    outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = outName,
    padj_cutoff = 0.01, lfcCut = 0.58, strict_lfc = T, saveFullRes = T,
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



cnts <- read.csv("../output/bcell/pseudoBulk/allCells_pb_matrix.csv")
cnts <- cnts %>%
    left_join(genez, by = "X") %>%
    filter(! ortholog_gene %in% featsTOexclude) %>%
    na.omit() %>%
    rownames_to_column() %>%
    select(-X, -input_gene, -rowname, -index) %>%
    column_to_rownames("ortholog_gene")
write.csv(cnts, file = "../output/bcell/pseudoBulk/allCells_pb_matrix.csv")

p_volc <- pseudoDEG(
    inDir = paste0("../output/bcell/pseudoBulk/"), 
    metaPWD = paste0("../output/bcell/pseudoBulk/allCells_deg_metaData.csv"),
    outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = "bcell_hu",
    padj_cutoff = 0.01, lfcCut = 0.58, strict_lfc = T, saveFullRes = T,
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


can_dat <- read.csv(paste0("../output/bcell/pseudoBulk/allCells/",
                   "bcell_hu_cluster_allCells_all_genes.csv"))

hu_dat <- read.csv(paste0("../output/human/pseudoBulk/allCells/",
                          "human_cluster_allCells_all_genes.csv"))

genez <- intersect(can_dat$gene, hu_dat$gene)

p <- crossSpeciesDEG(
    pwdTOspecies1 = paste0("../output/bcell/pseudoBulk/allCells/",
                           "bcell_hu_cluster_allCells_all_genes.csv"), 
    pwdTOspecies2 = paste0("../output/human/pseudoBulk/allCells/",
                           "human_cluster_allCells_all_genes.csv"), 
    species1 = "Canine", species2 = "Human", 
    nlab = 10, nlab_conflict = 5, nlab_axis = 2, 
    overlapGenes = overlapGenes, colUp = "red",
    colDwn = "blue", cONtrast = c("Tumor","Blood"), seed = 12,
    hjustvar = c(1,0.5,0,0.5,1,0,1,0),
    vjustvar = c(0.5,1,0.5,0,0,1,1,0),
    outDir = "../output/", saveGeneList = F
)
ggsave(paste("../output/", outName, "/", outName, "_pSigned_cfibVend.png", sep = ""), width = 7,height=5)


### Fig 2d - GO GSEA of DEGs
p <- plotGSEA(
#     pwdTOgeneList = paste0("../output/human/pseudoBulk/allCells/",
#                           "human_cluster_allCells_all_genes.csv"),
    geneList = genez, category = "C5", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = NULL, 
    lolli = T,
    saveRes = paste0("../output/", outName, "/c5_", outName, "_res.csv")
) + 
    theme(axis.title = element_text(size = 16))
p <- p + scale_x_continuous(
    limits = c(-80,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "log10(p.adj)"
) + 
    ggtitle("Gene ontology") + 
    theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste0("../output/", outName, "/", "gseaPlot_1.png"), width = 7, height = 7)

### Fig 2e - Reactome GSEA of DEGs
p <- plotGSEA(
#     pwdTOgeneList = paste0("../output/human/pseudoBulk/allCells/",
#                           "human_cluster_allCells_all_genes.csv"),
    geneList = genez, category = "C2", species = "dog", termsTOplot = 10, 
    upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "CP:REACTOME", 
    lolli = T,
    saveRes = paste0("../output/", outName, "/c2_", outName, "_res.csv")
) + theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12)
)
p <- p + scale_x_continuous(
    limits = c(-12,ceiling(max(p$data$x_axis)*1.05)), 
    breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
    name = "log10(p.adj)"
) + 
    ggtitle("Reactome") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_2.png"), width = 7, height = 7)

### Fig extra - ImmuneSigDB GSEA of DEGs
p <- plotGSEA(
#     pwdTOgeneList = paste0("../output/human/pseudoBulk/allCells/",
#                           "human_cluster_allCells_all_genes.csv"),
    geneList = genez, category = "C7", species = "dog", termsTOplot = 10, 
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
    name = "log10(p.adj)"
) + 
    ggtitle("ImmuneSigDB") + 
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = c(0.05, 1.05)
    )
ggsave(paste0("../output/", outName, "/", "gseaPlot_3.png"), width = 7, height = 7)
