#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")
outName <- "hu_blood_os"
# load10x(
#    din = "../external_data/human_OSLui_pbmcLee/", dout = "../output/s1_hu/",
#    outName = outName, testQC = F, nFeature_RNA_high = 10000, nFeature_RNA_low = 200, 
#    nCount_RNA_high = 100000, nCount_RNA_low = 100, percent.mt_high = 15, 
#    mt_pattern = "^MT-", removeDubs = FALSE, 
#    removeRBC_pal = FALSE, pal_feats = NULL, isolatePalRBC = FALSE
# )
# seu.obj <- sctIntegrate(
#     din = "../output/s1_hu/", outName = outName, dout = "../output/s2_human/",
#     vars.to.regress = "percent.mt", nfeatures = 2000
# )

# seu.obj <- dataVisUMAP(
#     seu.obj = seu.obj, outDir = "../output/s3/", outName = outName, 
#     final.dims = 45, final.res = 0.8, stashID = "clusterID", algorithm = 3, 
#     prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, 
#     assay = "integrated", saveRDS = T
# )


#Load in the processed data
seu.obj <- readRDS("../output/s3/hu_blood_os_res0.8_dims45_dist0.3_neigh30_S3.rds")

#Transfer predefined annotations
pbmc_anno <- read.csv("./metaData/pbmc_meta.csv") %>%
    separate(barcode, into = c("barcode", "sample"), sep = "-") %>%

    mutate(
        sample = recode(sample, 
                        "5" = "GSM4509015", "13" = "GSM4509023", 
                        "14" = "GSM4509024", "19" = "GSM4509029"
                       ),
        mergeVal = paste0(barcode, "-", sample),
        cellSource = "Blood"
    ) %>% 
    select(barcode, sample, mergeVal, cellSource, celltype)
colnames(pbmc_anno) <-   c("X", "GSM_ID", "mergeVal", "cellSource", "celltype.l1")
pbmc_anno$celltype.l2 <- pbmc_anno$celltype.l1
os_anno <- read.csv("./metaData/human_os_anno.csv") %>%
    mutate(
        mergeVal = paste0(substr(X, 1, 17), GSM_ID),
        cellSource = "Tumor"
    ) %>%
    select(X, GSM_ID, mergeVal, cellSource, celltype.l1, celltype.l2)
annos <- rbind(pbmc_anno, os_anno)

seu.obj@meta.data <- seu.obj@meta.data %>%
    rownames_to_column() %>%
    mutate(
        mergeVal = paste0(substr(rowname, 1, 17), orig.ident)
    ) %>%
    left_join(annos, by = "mergeVal") %>%
    column_to_rownames()

#Remove unannotated and cells without counterpart
seu.obj <- subset(seu.obj, invert = TRUE, 
                  cells = colnames(seu.obj)[is.na(seu.obj$celltype.l2)]
                 )
seu.obj$celltype.l2 <- as.factor(seu.obj$celltype.l2)
seu.obj <- subset(
    seu.obj, invert = TRUE, subset = celltype.l2 %in% c(
        "Cycling osteoblast", "Endothelial_cell", "FABP5_TAM", "Fibroblast",
        "Mast cell", "Mature_OC", "NR4A3_TAM", "Osteoblast", "Osteoblast_cycling", 
        "Platelet", "RBC", "TXNIP_TAM", "Uncategorized1", "Uncategorized2", 
        "preOC", "Plasma cell", "IFN-TAM" #exclude Plasma cell b/c none in human tumor
    )
)
seu.obj$celltype.l2 <- droplevels(seu.obj$celltype.l2)
Idents(seu.obj) <- "celltype.l2"
seu.obj <- RenameIdents(seu.obj, c(
    "B cell" = "B cell", "B cell, IgG+" = "B cell", "B cell, IgG-" = "B cell",

    "CD14_TIM" = "Monocyte", "classical Monocyte" = "Monocyte", 
    "intermediate Monocyte" = "Monocyte", "nonclassical Monocyte" = "Monocyte",

    "CD4 T cell" = "CD4", "CD4, EM-like" = "CD4", "CD4, non-EM-like" = "CD4", 

    "CD8 T cell" = "CD8", "CD8, EM-like" = "CD8", "CD8, non-EM-like" = "CD8", 
    "NK" = "CD8", "NK cell" = "CD8", "T_IFN" = "CD8", 

    "DC" = "DC", "cDC2" = "DC", "pDC" = "DC"
))
seu.obj$majorID_sub <- Idents(seu.obj)

gc()

seu.sub.list <- SplitObject(seu.obj, split.by = "orig.ident")
seu.obj <- indReClus(
    seu.list = seu.sub.list, seu.obj = NULL, 
    outDir = "../output/s2/", subName = "20241226_focused_hu", 
    preSub = T, vars.to.regress = "percent.mt"
)
seu.obj <- dataVisUMAP(
    seu.obj = seu.obj, outDir = "../output/s3/", 
    outName = "20241226_focused_hu", 
    final.dims = 30, final.res = 0.1, algorithm = 3,
    stashID = "clusterID", prefix = "integrated_snn_res.", 
    min.dist = 0.3, n.neighbors = 30, 
    assay = "integrated", saveRDS = T
)


#Load in the processed data
seu.obj <- readRDS("../output/s3/20241226_focused_hu_res0.1_dims30_dist0.3_neigh30_S3.rds")
table(seu.obj$celltype.l2, seu.obj$clusterID)

features <- c(
        "CD3G", "CD8A", "GZMA", "CD4", 
        "S100A12", "DLA-DRA","CD177", "CD68", 
        "FLT3", "MS4A1", "JCHAIN", "TOP2A"
)
p <- prettyFeats(
        seu.obj = seu.obj, nrow = 3, ncol = 4, features = features, color = "black"
) 
ggsave(paste0("../output/", outName, "/featPlots.png"), width = 12, height = 9)


## Fig 1c - UMAP colorized by majorID of integrated dataset
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "celltype.l2",
#               cols = c("#D4F3A3","#C8C3E2","#0066A5","#AC0535","#148B7D","#F17D00","hotpink"),
              pt.size = 0.5,
              label = T,
              label.box = T,
              shuffle = TRUE
) + NoLegend()
pi <- formatUMAP(plot = pi, smallAxes = T) 
ggsave(paste("../output/", outName, "/majorID_UMAP_tiny.png", sep = ""), width = 7, height = 7)



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

tumor_pal_sig <- genez[genez$input_gene %in% c(pal_feats,tumor.sig), ]$ortholog_gene

featsTOexclude <- c(
    rownames(seu.obj)[! rownames(seu.obj) %in% genez$ortholog_gene], #only retain 1:1
    tumor_pal_sig
)

# Complete DE analysis in each subcluser - human data
createPB(
    seu.obj = seu.obj, groupBy = "majorID_sub", comp = "cellSource", 
    biologicalRep = "orig.ident", outDir = paste0("../output/", outName, "/pseudoBulk/"), 
    min.cell = 25, grepTerm = "GSM450", grepLabel = c("Blood", "Tumor"), 
    featsTOexclude = featsTOexclude, lowFilter = T, dwnSam = F
)
pseudoDEG(
    inDir = paste0("../output/", outName, "/pseudoBulk/"), 
    metaPWD = paste0("../output/", outName, "/pseudoBulk/",
                     "majorID_sub_deg_metaData.csv"),
    outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = outName,
    padj_cutoff = 0.01, lfcCut = 0.58,  idents.1_NAME = "Tumor", 
    idents.2_NAME = "Blood", title = "Tumor vs Blood", fromFile = T, 
    returnVolc = F, filterTerm = "^ENS", mkDir = T, saveFullRes = T
)


### Evalutate exhaustion gene signature enrichment in human data
modulez <- list(
    "Exhausted" = c("BTLA", "CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT")
)
#run module score
seu.obj <- AddModuleScore(
    seu.obj, features = modulez, name = "_score"
)
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)
p <- VlnPlot(
    seu.obj, group.by = "majorID_sub", 
    features = "Exhausted", 
    split.by = "cellSource"
)
p <- p +
    theme(
        axis.title = element_blank()
    ) +
    ggtitle("Exhaustion enrichment score")
ggsave(paste0("../output/", outName, "/", outName, "viln.png"), width = 6, height = 3, scale = 1.25)

#Run stats on the enrichment scores
stat_df <- seu.obj@meta.data %>%
    select(orig.ident, majorID_sub, cellSource, Exhausted) %>%
    group_by(orig.ident, majorID_sub) %>%
    mutate(
        sum_score = mean(Exhausted)
    ) %>%
    ungroup() %>%
    distinct(sum_score, .keep_all = TRUE)
statz <- compare_means(sum_score ~ cellSource, group.by = "majorID_sub", stat_df) %>%
    select(-p.format, -.y.) %>%
    mutate(
        p.signif = symnum(p.adj, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                          symbols = c("****", "***", "**", "*", "ns"))
    )
statz
stat_df %>%
    group_by(cellSource, majorID_sub) %>%
    summarize(MEAN = mean(sum_score))

# n.s. in CD8, but * sig in CD4s

### Load in  canine data and run DE with same feats used for human analysis
seu.obj.can <- readRDS("../output/s3/bloodANDtils_filtered_allCells_S3.rds")
Idents(seu.obj.can) <- "orig.ident"
seu.obj.can <- RenameIdents(seu.obj.can, c(
        "run_count_osa_pbmc_1" = "bl_1", 
        "run_count_osa_pbmc_2" = "bl_2", 
        "run_count_osa_pbmc_3" = "bl_3", 
        "run_count_osa_pbmc_4" = "bl_4",
        "run_count_osa_pbmc_5" = "bl_5",
        "run_count_pbmc_tp1_pt1" = "bl_6",
        "run_count_pbmc_tp1_pt2" = "bl_7", 
        "run_count_pbmc_tp1_pt3" = "bl_8", 
        "run_count_pbmc_tp1_pt4" = "bl_9", 
        "run_count_pbmc_tp1_pt5" = "bl_10",
        "run_count_tumor_no_tx_1_1" = "tils_1",
        "run_count_tumor_no_tx_1_2" = "tils_1",
        "run_count_tumor_no_tx_2_1" = "tils_2", 
        "run_count_tumor_no_tx_2_2" = "tils_2", 
        "run_count_tumor_no_tx_4" = "tils_3", 
        "run_count_tumor_no_tx_5" = "tils_4",
        "run_count_no_tx_tumor_6" = "tils_5",
        "run_count_tumor_no_tx_7" = "tils_6"
))
seu.obj.can$name <- Idents(seu.obj.can)
#pull in subcluster analysis annotations -- store as celltype.l3
files <- list.files(path = "../output/annotations/", pattern = ".csv", 
                    all.files = FALSE, full.names = T)
df.list <- lapply(files, read.csv, header = T)
anno.df <- do.call(rbind, df.list)
seu.obj.can <- AddMetaData(
    seu.obj.can, setNames(anno.df$majorID_sub, anno.df$X), col.name = "celltype.l3"
)
seu.obj.can$major <- as.factor(ifelse(
    seu.obj.can$celltype.l3 == "Plasma_cell (c2)",
    "Plasma",
    as.character(seu.obj.can$major)
))
Idents(seu.obj.can) <- "major"
seu.obj.can <- RenameIdents(seu.obj.can, c(
    "CD4 T cell" = "CD4",
    "CD8 T cell" = "CD8",
    "Dendritic cell" = "DC",
    "B cell" = "B cell", 
    "Monocyte" = "Monocyte", 
    "Neutrophil" = "Exlude", 
    "Plasma" = "Exlude"
))
seu.obj.can$major <- Idents(seu.obj.can)

# Find the genes used in the human analysis and only use those in the canine
featsUseInHuman <- rownames(seu.obj)[! rownames(seu.obj) %in% featsTOexclude]
featsUseInHuman <- genez[genez$ortholog_gene %in% featsUseInHuman, ]$input_gene

featsTOexclude <- c(
    rownames(seu.obj.can)[! rownames(seu.obj.can) %in% featsUseInHuman], #only use genes evalutated in human
    tumor_pal_sig
)

# Complete DE analysis in each subcluser - canine data using same filtering
createPB(
    seu.obj = seu.obj.can, groupBy = "major", comp = "cellSource", 
    biologicalRep = "name", outDir = paste0("../output/", outName, "/canine/pseudoBulk/"), 
    min.cell = 25, grepTerm = "tils", grepLabel = c("Tumor", "Blood"), 
    featsTOexclude = featsTOexclude, lowFilter = T, dwnSam = F
)
pseudoDEG(
    inDir = paste0("../output/", outName, "/canine/pseudoBulk/"), 
    metaPWD = paste0("../output/", outName, "/canine/pseudoBulk/",
                     "major_deg_metaData.csv"),
    outDir = paste0("../output/", outName, "/canine/pseudoBulk/"), outName = outName,
    padj_cutoff = 0.01, lfcCut = 0.58,  idents.1_NAME = "Tumor", 
    idents.2_NAME = "Blood", title = "Tumor vs Blood", fromFile = T, 
    returnVolc = F, filterTerm = "^ENS", mkDir = T, saveFullRes = T
)

# Convert gene symbols to human
lapply(levels(seu.obj.can$major), function(x){
    res <- read.csv(paste0("../output/", outName, "/canine/pseudoBulk/", x, "/", 
                            outName, "_cluster_", x, "_all_genes.csv"))
    res <- res %>%
        left_join(genez, by = c("gene" = "input_gene")) %>%
        select(-index, -X)
    res$gene <- res$ortholog_gene
    res$ortholog_gene <- NULL
    write.csv(res, file = paste0("../output/", outName, "/canine/pseudoBulk/", x, "/", 
                                 outName, "_cluster_", x, "_all_genes.csv")
             )    
})

# Plot DE results comparing expression in the two specices
pis <- lapply(c("CD4", "CD8", "DC", "Monocyte"), function(x){
    p <- crossSpeciesDEG(
        pwdTOspecies1 = paste0("../output/", outName, "/canine/pseudoBulk/", x, "/", 
                                 outName, "_cluster_", x, "_full_all_genes.csv"), 
        pwdTOspecies2 = paste0("../output/hu_blood_os/pseudoBulk/", x, "/",
                               "hu_blood_os_cluster_", x, "_full_all_genes.csv"), 
        species1 = "Canine", species2 = "Human", 
        nlab = 10, nlab_conflict = 5, nlab_axis = 2, 
        overlapGenes = overlapGenes, colUp = "red",
        colDwn = "blue", cONtrast = c("Tumor","Blood"), seed = 12,
        hjustvar = c(1,0.5,0,0.5,1,0,1,0),
        vjustvar = c(0.5,1,0.5,0,0,1,1,0), outName = x, 
        outDir = paste0("../output/", outName, "/"), saveGeneList = T, useFC = T
    )
    ggsave(paste0("../output/", outName, "/", outName, "_pSigned_", gsub(" ", "_", x), ".png"), width = 7, height = 5)
    return(p)
})

#Add titles
titlez <- c("CD4 T cell", "CD8 T cell", "Dendritic cell", "Monocyte")
pis <- lapply(1:4, function(i){
    pis[[i]] + 
        ggtitle(titlez[i]) +
        theme(
            plot.title = element_text(size = 18, hjust = 0.05),
            axis.title = element_text(size = 18)
        )
})

pi <- Reduce("+", pis) +
    patchwork::plot_layout(ncol = 2, nrow = 2, axis_titles = "collect")
ggsave(paste0("../output/", outName, "/", outName, "_pSigned_compiled.png"), width = 7, height = 5, scale = 2)

cts <- c("CD4", "CD8", "DC", "B cell", "Monocyte")
queryTerms <- c(
    "TCELL|THYMOCYTE|THYMUS", "CD8_T_CELL|T_CELL|THYMOCYTE|THYMUS",
    "_DC_|BMDC", "BCELL|PLASMA_CELL|B_CELL", "MONOCYTE|MACROPHAGE"
)
lapply(1:5, function(i){
    ct <- cts[i]
    queryTerm <- queryTerms[i]
    dat <- as.data.frame(pis[[i]]$data)
    shared_genez <- dat[dat$species == "Both" & dat$direction == "up", ]$gene
    # GO
    p <- plotGSEA(
        geneList = shared_genez, category = "C5", species = "human", termsTOplot = 10, 
        upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = NULL, 
        lolli = T,
        saveRes = paste0("../output/", outName, "/c5_", gsub(" ", "_", ct), "_res.csv")
    )
        p <- p + scale_x_continuous(
            limits = c(-80,ceiling(max(p$data$x_axis)*1.05)), 
            breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
            name = "-log10(p.adj)"
        ) + 
            ggtitle("Gene ontology") + 
            theme(plot.title = element_text(size = 20, hjust = 0.5))
        ggsave(paste0("../output/", outName, "/", "gseaPlot_c5", gsub(" ", "_", ct), ".png"), width = 7, height = 7)
    # Reactome
    p <- plotGSEA(
        geneList = shared_genez, category = "C2", species = "human", termsTOplot = 10, 
        upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "CP:REACTOME", 
        lolli = T,
        saveRes = paste0("../output/", outName, "/c2_", gsub(" ", "_", ct), "_res.csv")
    )
    # Nothing enriched in B cells so skip
    if(ct != "B cell"){
        p <- p + scale_x_continuous(
            limits = c(-12,ceiling(max(p$data$x_axis)*1.05)), 
            breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
            name = "-log10(p.adj)"
        ) + 
            ggtitle("Reactome") + 
            theme(
                plot.title = element_text(size = 20, hjust = 0.5),
                legend.position = c(0.05, 1.05),
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 12)
            )
        ggsave(paste0("../output/", outName, "/", "gseaPlot_c2", gsub(" ", "_", ct), ".png"), width = 7, height = 7)
    }
    # ImmuneSigDB
    p <- plotGSEA(
        geneList = shared_genez, category = "C7", species = "human", termsTOplot = 10, 
        upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
        lolli = T, filterTerm = queryTerm,
        saveRes = paste0("../output/", outName, "/c7_", gsub(" ", "_", ct), "_res.csv")
    )
        p <- p + scale_x_continuous(
            limits = c(-60,ceiling(max(p$data$x_axis)*1.05)), 
            breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
            name = "-log10(p.adj)"
        ) + 
            ggtitle("ImmuneSigDB") + 
            theme(
                plot.title = element_text(size = 20, hjust = 0.5),
                legend.position = c(0.05, 1.05),
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 12)
            )
        ggsave(paste0("../output/", outName, "/", "gseaPlot_c7", gsub(" ", "_", ct), ".png"), width = 7, height = 7)
    # Run ImmuneSigDB using human DE results
    p <- plotGSEA(
        pwdTOgeneList = paste0("../output/hu_blood_os/pseudoBulk/", ct, "/",
                               "hu_blood_os_cluster_", ct, "_all_genes.csv"),
        category = "C7", species = "human", termsTOplot = 10, 
        upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
        lolli = T, filterTerm = queryTerm,
        saveRes = paste0("../output/", outName, "/c7_", gsub(" ", "_", ct), "_human_res.csv")
    )
    ggsave(paste0("../output/", outName, "/human_", gsub(" ", "_", ct), "_c7.png"), width = 7, height = 7)
    # Run Reactome using human DE results    
    p <- plotGSEA(
        pwdTOgeneList = paste0("../output/hu_blood_os/pseudoBulk/", ct, "/",
                               "hu_blood_os_cluster_", ct, "_all_genes.csv"),
        category = "C2", species = "human", termsTOplot = 10, 
        upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "CP:REACTOME", 
        lolli = T,
        saveRes = paste0("../output/", outName, "/c2_", gsub(" ", "_", ct), "_human_res.csv")
    )
    ggsave(paste0("../output/", outName, "/human_", gsub(" ", "_", ct), "_c2.png"), width = 7, height = 7)
})

### Clean gene classifications for supplemental data
pValFilter <- 0.01
df.list <- lapply(c("CD4", "CD8", "DC", "Monocyte"), function(x){
    pValFilter <- -log10(pValFilter)
    df <- read.csv(paste0("../output/hu_blood_os/", x, 
                          "_Canine_v_Human_Tumor_Blood.csv")) %>%
        mutate(
            significance_classification = case_when(
                abs(signed_pVal_spec1) > pValFilter & abs(signed_pVal_spec2) > pValFilter ~ "Both",
                abs(signed_pVal_spec1) > pValFilter & abs(signed_pVal_spec2) < pValFilter ~ "Canine Only",
                abs(signed_pVal_spec1) < pValFilter & abs(signed_pVal_spec2) > pValFilter ~ "Human Only",
                abs(signed_pVal_spec1) < pValFilter & abs(signed_pVal_spec2) < pValFilter ~ "n.s."
            ),
            FC_classification = case_when(
                direction == "up" ~ "FC > 0 in both",
                direction == "dwn" ~ "FC < 0 in both",
                direction == "axis1" ~ "FC > 0 in canine",
                direction == "axis2" ~ "FC > 0 in human",
                direction == "axis3" ~ "FC < 0 in canine",
                direction == "axis4" ~ "FC > 0 in human",
                direction == "conflict1" ~ "FC > 0 in canine AND FC < 0 in human",
                direction == "conflict2" ~ "FC < 0 in canine AND FC > 0 in human",
                direction == "n.s." ~ "Not DE in either species"
            ),
            quadrant = case_when(
                direction == "up" ~ "Top right (up in both)",
                direction == "dwn" ~ "Bottom left (down in both)",
                direction == "axis1" ~ "Right axis (Up in canine)",
                direction == "axis2" ~ "Up axis (Up in human)",
                direction == "axis3" ~ "Left axis (Down in canine)",
                direction == "axis4" ~ "Down axis (Down in human)",
                direction == "conflict1" ~ "Bottom right (up in canine, down in human)",
                direction == "conflict2" ~ "Top left (down in canine, up in human)",
                direction == "n.s." ~ "Not DE in either species"
            ),
            final_classification = case_when(
                direction == "up" ~ "Conserved",
                direction == "dwn" ~ "Conserved",
                direction == "axis1" ~ "Ambiguous",
                direction == "axis2" ~ "Ambiguous",
                direction == "axis3" ~ "Ambiguous",
                direction == "axis4" ~ "Ambiguous",
                direction == "conflict1" ~ "Divergent",
                direction == "conflict2" ~ "Divergent",
                direction == "n.s." ~ "Conserved (n.s.)"
            )
        ) %>%
        select(gene, log2FoldChange_spec1, signed_pVal_spec1, 
               log2FoldChange_spec2, signed_pVal_spec2, significance_classification,
               FC_classification, quadrant, final_classification
              )
    colnames(df) <- c("gene", "log2FoldChange_Canine", "signed_pVal_Canine", 
                      "log2FoldChange_Human", "signed_pVal_Human", 
                      "significance_classification", "FC_classification", "quadrant", 
                      "final_classification")
    df$cell_type <- x
    return(df)
})
write.csv(do.call(rbind, df.list), "../output/supplementalData/human_canine_comp.csv", row.names = F)


# ### Core cross species TIMs signature
# #Compare to human BrCA TEMo vs Mono
# hu_degs <- read.csv("./metaData/brCA_TEMo_VS_Mo.csv") %>%
#     column_to_rownames(var = "gene.name")

# hu_degs <- read.csv("./metaData/lung_TIM_vs_Mono_Zilionis.csv") %>%
#     column_to_rownames(var = "gene")

# hu_degs <- hu_degs[hu_degs$FDR < 0.01, ]

# x <- "Monocyte"
# degs_can <- read.csv(paste0("../output/", outName, "/canine/pseudoBulk/", x, "/", 
#                          outName, "_cluster_", x, "_all_genes.csv"))
# degs_hu <- read.csv(paste0("../output/hu_blood_os/pseudoBulk/", x, "/",
#                        "hu_blood_os_cluster_", x, "_all_genes.csv"))

# intersect(degs_can$gene, degs_hu$gene)
# intersect(degs_can$gene, rownames(hu_degs))

# intersect(degs_hu$gene, rownames(hu_degs))


# genez <- orthogene::convert_orthologs(
#     gene_df = rownames(hu_degs), gene_output = "columns", 
#     input_species = "human", output_species = "dog", 
#     non121_strategy = "drop_both_species"
# )
# genez <- genez %>% 
#     left_join(rownames_to_column(hu_degs), by = c("input_gene" = "rowname"))

# degs_dog <- read.csv(paste0(
#     "../output/", outName, "/pseudoBulk/allCells/", outName, 
#     "_cluster_allCells_all_genes.csv"
# ))
# degs_dog <- degs_dog[degs_dog$log2FoldChange > 0, ]
# genez_up <- genez[genez$FDR < 0.01 & genez$fold_change > 1.5, ]$ortholog_gene
# overlap_tims <- intersect(genez_up, degs_dog$gene)

# hu_only <- genez_up[! genez_up %in% overlap_tims]
# can_only <- degs_dog$gene
# can_only <- genez_up[! can_only %in% overlap_tims]


df.res <- read.csv("../output/supplementalData/human_canine_comp.csv")

# Extract overlap percentages for manuscript
df.res_summary <- df.res %>%
    filter(significance_classification != "n.s.") %>%
    group_by(cell_type, final_classification) %>%
    summarize(
        COUNT = n()
    ) %>%
    group_by(cell_type) %>%
    mutate(
        PCT = round(COUNT / sum(COUNT), 2)
    ) %>%
    filter(final_classification == "Divergent") 

diverg_gene <- df.res %>%
    filter(cell_type != "B cell" & final_classification == "Divergent") %>%
    pull(gene)

ambig_gene <- df.res %>%
    filter(cell_type != "B cell" & final_classification == "Ambiguous") %>%
    pull(gene)

ambig_gene <- df.res %>%
    filter(cell_type != "B cell" & final_classification == "Conserved") %>%
    pull(gene)

dput(setNames(df.res_summary$PCT, df.res_summary$cell_type))

cts <- c("CD4", "CD8", "DC", "Monocyte")
queryTerms <- c(
#     "TCELL|THYMOCYTE|THYMUS", "CD8_T_CELL|T_CELL|THYMOCYTE|THYMUS",
    "_CD4_", "_CD8_",
    "_DC_|BMDC", "BCELL|PLASMA_CELL|B_CELL", "MONOCYTE|MACROPHAGE"
)

lapply(1:5, function(i){

    # Strict conserved classifcation
#     gene_res <- list(
#         "Human" = df.res[df.res$signed_pVal_Human > 2 & 
#                          df.res$cell_type == x, ]$gene,
#         "Canine" = df.res[df.res$signed_pVal_Canine > 2 & 
#                          df.res$cell_type == x, ]$gene
#     )
    
    # Lenient conserved classification
    gene_res <- list(
        "Human" = df.res[df.res$log2FoldChange_Human > 0 & 
                         df.res$significance_classification != "n.s." & 
                         df.res$cell_type == cts[i], ]$gene,
        "Canine" = df.res[df.res$log2FoldChange_Canine > 0 & 
                         df.res$significance_classification != "n.s." & 
                         df.res$cell_type == cts[i], ]$gene
    )

    #Make the plot
    VennDiagram::venn.diagram(
        x = gene_res,
        filename = paste0("../output/", outName, "/", cts[i], "_venn.png"),
        output = TRUE,
        imagetype = "png",
        fontfamily = "sans",
        cat.fontfamily = "sans",
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.dist = c(0.05, 0.05),
        margin = 0.05,
        disable.logging = TRUE,
        height = 2250, 
        width = 2250
    )

    gene_res$Core <- intersect(gene_res[[1]], gene_res[[2]])

    #Run Reactome for Canine, Human, and Core - results to be used for heatmap
    lapply(1:3, function(ii){
        res <- gene_res[[ii]]
        saveName <- names(gene_res)[ii]
        p <- plotGSEA(
            geneList = res, category = "C7", species = "human", termsTOplot = 10, 
            upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "IMMUNESIGDB", 
            lolli = T, filterTerm = queryTerms[i],
            saveRes = paste0("../output/", outName, "/c7_", gsub(" ", "_", cts[i]), "_", saveName, "_res.csv")
        )#         p <- plotGSEA(
#             geneList = res, category = "C2", species = "human", termsTOplot = 10, 
#             upOnly = T, trunkTerm = T, pvalueCutoff = 0.05, subcategory = "CP:REACTOME", 
#             lolli = T,
#             saveRes = paste0("../output/", outName, "/c2_", gsub(" ", "_", cts[i]), "_", saveName, "_res.csv")
#         )
        if(is.null(p)) {
            message("No pathways enriched....")
        } else {
            p <- p + scale_x_continuous(
                limits = c(-12,ceiling(max(p$data$x_axis)*1.05)), 
                breaks = c(0,ceiling(max(p$data$x_axis)*1.05)/2,ceiling(max(p$data$x_axis)*1.05)),
                name = "-log10(p.adj)"
            ) + 
                ggtitle("Reactome") + 
                theme(
                    plot.title = element_text(size = 20, hjust = 0.5),
                    legend.position = c(0.05, 1.05),
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 12)
                )
            ggsave(paste0("../output/", outName, "/", "gseaPlot_c7", gsub(" ", "_", cts[i]), "_", saveName, ".png"), width = 7, height = 7)
        }
    })
})

#Organize gsea results
df.list <- lapply(c("CD4", "CD8", "DC", "Monocyte"), function(ct){
    res.df.list <- lapply(c("Canine", "Human"), function(gset){
        inFile <- paste0("../output/", outName, "/c7_", ct, "_", gset, "_res.csv")

        if(file.exists(inFile)){
            df <- read.csv(inFile, header = T)
            df$CellType <- ct
            df$Database <- "C7"
            df$gset <- gset
            message(paste("Loading", ct, gset))
            return(df)
        }
    })
    res.df <- do.call(rbind, res.df.list)
    return(res.df)
})
res.df <- do.call(rbind, df.list)

#data cleaning and conversion to matrix
df.list <- lapply(c("CD4", "CD8", "DC", "Monocyte"), function(ct){
    df <- res.df %>% 
        filter(CellType == ct)
    
    #fix labeling error - all get the GSE* prefix to correct downstream error
    df$Description <- df$ID
    df$Description <- unlist(lapply(df$Description, function(x){
            cutz <- unlist(gregexpr(pattern = '_', x))
            cut_loc <- max(cutz[cutz < 60])
            if(nchar(x) > 60){
                paste0(substr(x, 1, cut_loc - 1), "\n", substr(x, cut_loc + 1, nchar(x)))
            } else {
                x
            }
    }))
    df$Description <- gsub("_", " ", df$Description)

    df <- df %>%
        mutate(
            row = row_number()
        ) %>% 
        select(CellType, Description, x_axis, row, gset) %>%
        pivot_wider(names_from = gset, values_from = x_axis) %>% 
        as.data.frame()

    mat <- aggregate(
        df[ , 4:ncol(df)], by = df[2], 
        FUN = function(x){sum(as.numeric(x), na.rm = TRUE)}
    ) %>%
        column_to_rownames("Description") %>%
        as.matrix()
    
    mat <- mat[rev(order(rowSums(mat))), ]
    row_split <- unlist(unname(apply(mat, MARGIN = 1, function(x){
        case_when(
            sum(x > 0) > 1 ~ "both",
            sum(x > 0) == 1 & x[1] > 0 ~ "canine",
            sum(x > 0) == 1 & x[2] > 0 ~ "human"
    )})))
    position <- c(
        which(row_split == "both")[1:10],
        which(row_split == "canine")[1:5],
        which(row_split == "human")[1:5]
    )

    labelz <- rownames(mat)[c(position)]
    labelz <- gsub("^\\S+\\s*", "", labelz)

    ra_right <- rowAnnotation(key_feats = anno_mark(
        at = position, labels = labelz, 
        labels_gp = gpar(fontsize = 9)
    ))
    
    #create the heatmap
    ht <- Heatmap(
        mat,
        name = "-log10(FDR)",
        column_title = ct,
        col = circlize::colorRamp2(c(0, 15), c("#EEEEEE", "red")),
        cluster_rows = F,
        row_title_gp = gpar(fontsize = 24),
        row_title = NULL,
        border = TRUE,
        show_row_names = F,
        right_annotation = ra_right,
        row_split = row_split,
        cluster_columns = F,
        show_column_names = T,
        heatmap_legend_param = list(
            direction = "horizontal",
            legend_width = unit(6, "cm")
        )
    )
    png(file = paste0("../output/", outName, "/",  outName, "_", ct, "_heatTest.png"), width = 2500, height = 3000, res = 400)
    par(mfcol = c(1, 1))   
    draw(
        ht, padding = unit(c(2, 2, 2, 2), "mm"), heatmap_legend_side = "top"
    )
    dev.off()

    #clean data for source
    mat_dat <- arrange(cbind(as.data.frame(mat), row_split), row_split)
    mat_dat$cell_type <- ct
    return(mat_dat)

#     as.data.frame(table(row_split)) %>%
#         mutate(
#             tot = sum(Freq),
#             PCT = round((Freq / sum(Freq)) * 100, 2)
#         )
    
})
mat_dat <- do.call(rbind, df.list)

write.csv(mat_dat, "../output/supplementalData/supplemental_data_10.csv")

df_hu <- read.csv("../output/hu_blood_os/c7_Monocyte_Human_res.csv")
df_cn <- read.csv("../output/hu_blood_os/c7_Monocyte_Canine_res.csv")
total_pathways <- length(unique(c(df_cn$ID, df_hu$ID)))
overlap_pathways <- intersect(df_cn$ID, df_hu$ID)
can_pathways <- length(df_cn$ID[ ! df_cn$ID %in% overlap_pathways])
hu_pathways <-  length(df_hu$ID[ ! df_hu$ID %in% overlap_pathways])
num_overlap_pathways <- length(intersect(df_cn$ID, df_hu$ID))


cross_validate <- function(
    ct = NULL,
    disc_features = NULL

) {
    geneListConserved <- df.res[df.res$cell_type == ct & 
           df.res$final_classification == "Conserved" &
           df.res$significance_classification != "n.s.", ]$gene
    message(paste0("Features discussed that were conserved: ",
                   paste(disc_features[disc_features %in% geneListConserved], 
                         collapse = ", ")))
    message(paste0("Features discussed that were NOT conserved: ",
                   paste(disc_features[! disc_features %in% geneListConserved], 
                         collapse = ", ")))
            
    geneListDivergent <- df.res[df.res$cell_type == ct & 
           df.res$final_classification == "Divergent" &
           df.res$significance_classification != "n.s.", ]$gene
    message(paste0("Features discussed that were divergent: ",
                   paste(disc_features[disc_features %in% geneListDivergent], 
                         collapse = ", ")))
    #plot human data
    seu.obj$cellSource <- as.factor(seu.obj$cellSource)
    seu.obj.sub <- subset(seu.obj, subset = majorID_sub == ct)
    pi <- VlnPlot(
        object = seu.obj.sub,
        pt.size = 0,
        same.y.lims = F,
        group.by = "cellSource",
        combine = T,
        add.noise = T, 
        cols = c("#F8766D", "#00BFC4"),
        stack = T,
        fill.by = "ident",
        flip = T,
        features = disc_features
    ) + 
        NoLegend() + 
        theme(
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5)
        ) +
        ggtitle("Human")
    pi <- pi + geom_jitter(mapping = aes(color = ident), alpha = 0.5, data = pi$data)
    ggsave(plot = pi, paste0("../output/", outName, "/", ct, "_human_splitViln.png"), width = 5, height = 9)

    #plot canine data
    seu.obj.can.sub <- subset(seu.obj.can, subset = major == ct)
    pi <- VlnPlot(
        object = seu.obj.can.sub,
        pt.size = 0,
        same.y.lims = F,
        group.by = "cellSource",
        combine = T,
        add.noise = T, 
        cols = c("#F8766D", "#00BFC4"),
        stack = T,
        fill.by = "ident",
        flip = T,
        features = disc_features
    ) + 
        NoLegend() + 
        theme(
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5)
        ) +
        ggtitle("Canine")
    pi <- pi + geom_jitter(mapping = aes(color = ident), alpha = 0.5, data = pi$data)
    ggsave(plot = pi, paste0("../output/", outName, "/", ct, "_canine_splitViln.png"), width = 5, height = 9)
}

Idents(seu.obj.can) <- "cellSource"
seu.obj.can <- RenameIdents(seu.obj.can, c(
    "Peripheral" = "Blood",
    "TILs" = "Tumor"
))
seu.obj.can$cellSource <- Idents(seu.obj.can)
seu.obj.can$cellSource <- as.factor(seu.obj.can$cellSource)

cross_validate(
    ct = "CD4",
    disc_features = c("SELL", "LEF1", "TMEM154", "TNFRSF18", "TNFRSF4", 
                      "HAVCR1", "CXCL13", "PDCD1", "LAG3","IL4I1")
)

cross_validate(
    ct = "CD8",
    disc_features = c("TRIM22", "LEF1", "PTGDR", "CD53", "CX3CR1", "HAVCR2", 
                      "TNFSF9", "LAG3","TNFRSF18")
)

cross_validate(
    ct = "Monocyte",
    disc_features = c("LTF", "IL16","CCR2", "IL1A", "OSM", "CD274", 
                      "PTGES", "C1QC", "CXCL10", "CXCL16", "CCL19", "CCL5", 
                      "CCL7", "CCL8")
)

cross_validate(
    ct = "DC",
    disc_features = c("IL16", "IRF4", "CXCR4", "VEGFA", "IL1A", "IL4I1", "CCR7",
                      "CD274", "IDO1")
)

cross_validate(
    ct = "B cell",
    disc_features = c("LYZ", "BTLA", "ADD3", "IGHM", "FOS", "FOSB", "FAS", 
                      "IGFLR1", "EGFL7")
)

# Features discussed that were conserved: TNFRSF18, TNFRSF4, CXCL13, PDCD1, LAG3, IL4I1
# Features discussed that were NOT conserved: SELL, LEF1, TMEM154, HAVCR1
# Features discussed that were divergent: 

# Features discussed that were conserved: HAVCR2, LAG3, TNFRSF18
# Features discussed that were NOT conserved: LEF1, PTGDR, CX3CR1, TNFSF9
# Features discussed that were divergent: TRIM22, CD53

# Features discussed that were conserved: OSM, C1QC, CXCL10, CXCL16, CCL8
# Features discussed that were NOT conserved: LTF, IL16, IL1A, CD274, CCL19, CCL7
# Features discussed that were divergent: CCR2, PTGES, CCL5

# Features discussed that were conserved: IL4I1, CCR7, CD274, IDO1
# Features discussed that were NOT conserved: IL16, IRF4, CXCR4, VEGFA, IL1A
# Features discussed that were divergent: 

# Features discussed that were conserved: FOS, FOSB, FAS, IGFLR1
# Features discussed that were NOT conserved: LYZ, ADD3, IGHM, EGFL7
# Features discussed that were divergent: BTLA