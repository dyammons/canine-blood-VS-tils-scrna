#!/usr/bin/Rscript

### Author: Dylan Ammons
### Analysis code for [title here]

### Purpose: 
# Load in the integrated object generated in `1_preProcessData.R` and complete 
# abundance analysis and provide overview of cell type classifications.

### Analysis approach (3 steps):

## (1) Data preparation for analysis - Load in data and clean up metadata to 
# enable plotting

## (2) Plot primary figures (Figure 1) - Complete visualization and run 
# differential abundance analysis

## (3) Complete differential expression (DE) analysis and make supplemental 
# figures - first pass completes DE within each major cell type using all 
# features. Results from DE are then loaded in and used to identify DEGs that
# are conserved across all cell types. These genes are considered the tissue
# signatures and are then excluded during a second round of DE analysis within
# each cell type.

############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   begin all cells analysis   ######## <<<<<<<<<<<<<<
############################################# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


### (1) Data preparation for analysis

## Load custom functions, packages, and set output path
source("/pl/active/dow_lab/dylan/repos/scrna-seq/analysis-code/customFunctions.R")
library(UpSetR)
outName <- "allCells"

## Load in data preprocessed and clean up metadata
seu.obj <- readRDS("../output/s3/20230505_bloodANDtils_QCfiltered_2500_QCfiltered_2500_res1.2_dims40_dist0.5_neigh60_S3.rds")
seu.obj$conSense <- ifelse(
        seu.obj$cellSource == "TILs", 
        seu.obj$celltype.l3,
        seu.obj$celltype.l3_pbmc
)
Idents(seu.obj) <- "orig.ident"
seu.obj <- RenameIdents(seu.obj, c(
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
seu.obj$name <- Idents(seu.obj)
#Load in the cell type consensus classifications
seu.obj <- loadMeta(
        seu.obj = seu.obj, seu.file = NULL, metaFile = "../metaData/cellTypez.csv", 
        groupBy = "conSense", metaAdd = "major", save = FALSE, outName = "", header = TRUE
)

## Exclude cells that do not have a counterpart in both datasets
Idents(seu.obj) <- "major"
seu.obj <- subset(seu.obj, invert = T, ident = c(
        "Eosinophil", "DN T cell","gd T cell", "CD4+, IFN signature",
        "Basophil", "CD34+ Unclassified", "T_IFN"
        )
)
#Rename to match cell type across datasets
seu.obj$major <- droplevels(seu.obj$major)
seu.obj <- RenameIdents(seu.obj, c("Cycling T cell" = "T_cycling"))
seu.obj$major <- Idents(seu.obj)
#Exclude any remaining cells that do not have an ID associated with the barcode
seu.obj$major <- as.factor(ifelse(is.na(seu.obj$major), "NA", as.character(seu.obj$major)))
Idents(seu.obj) <- "major"
seu.obj <- subset(seu.obj,invert = TRUE, ident = "NA")
#Rename and reorder major cell type levels
seu.obj$major <- droplevels(as.factor(seu.obj$major))
Idents(seu.obj) <- "major"
seu.obj <- RenameIdents(seu.obj, c(
        "T_cycling" = "Cycling T cell", "bcell" = "B cell" ,
        "cyto"  = "CD8 T cell", "dc"  = "Dendritic cell",
        "helper"  = "CD4 T cell", "mono" = "Monocyte",
        "neut" = "Neutrophil"
        )
)
seu.obj$major <- Idents(seu.obj)
seu.obj$major <- factor(seu.obj$major, levels = levels(seu.obj$major)[c(5,3,1,2,6,4,7)])
seu.obj$majorID_sub <- seu.obj$major
seu.obj$major <- droplevels(as.factor(seu.obj$major))


### (2) Plot primary figures (Figure 1)

## Fig extra - UMAP of unsupervised clustering results using integrated dataset
pi <- DimPlot(
        seu.obj, 
        reduction = "umap", 
        group.by = "clusterID",
        pt.size = 0.5,
        label = TRUE,
        label.box = TRUE,
        shuffle = TRUE
)
pi <- cusLabels(
        plot = pi, shape = 21, size = 10, textSize = 6, 
        alpha = 0.8, labCol = "black"
        ) + 
        NoLegend() + 
        theme(
                axis.title = element_blank(),
                panel.border = element_blank()
        )
ggsave(paste("../output/", outName, "/rawUMAP.png", sep = ""), width = 7, height = 7)

## Fig 1c - UMAP colorized by majorID of integrated dataset
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
ggsave(paste("../output/", outName, "/majorID_UMAP_tiny.png", sep = ""), width = 7, height = 7)

## Fig 1d - feature plots of cell type defining genes
features <- c(
        "CD3G", "CD8A", "GZMA", "CD4", 
        "S100A12", "DLA-DRA","ANPEP", "CD68", 
        "FLT3", "MS4A1", "JCHAIN", "TOP2A"
)
p <- prettyFeats(
        seu.obj = seu.obj, nrow = 3, ncol = 4, features = features, color = "black"
) 
ggsave(paste("../output/", outName, "/featPlots.png", sep = ""), width = 12, height = 9)

## Fig 1e - boxplots comparing cell type proportions
p <- skewPlot(
        seu.obj, groupBy = "major", sigTextSpacing = 0.035, 
        yAxisLabel = "Percent of all cells",
        dout = paste0("../output/", outName), outName = outName
)
ggsave(paste0("../output/", outName, "/barchart.png"), width = 8, height = 3)


### (3) Complete DE and make supplemental figures

#Store gene list associated with platelets
#Gene lists from https://doi.org/10.5281/zenodo.7884518
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

## Complete pseudobulk DGE by each major cell type
createPB(
        seu.obj = seu.obj, groupBy = "majorID_sub", comp = "cellSource", 
        biologicalRep = "orig.ident", lowFilter = T, dwnSam = F, min.cell = 5,
        outDir = paste0("../output/", outName, "/pseudoBulk/"), 
        grepTerm = "tumor", grepLabel = c("TILs", "Blood")
)
pseudoDEG(
        metaPWD = paste0("../output/", outName, "/pseudoBulk/majorID_sub_deg_metaData.csv"),
        padj_cutoff = 0.05, lfcCut = 0.58, 
        outDir = paste0("../output/", outName, "/pseudoBulk/"), outName = outName, 
        idents.1_NAME = "TILs", idents.2_NAME = "Blood",
        inDir = paste0("../output/", outName, "/pseudoBulk/"), title = "All cells", 
        filterTerm = "ZZZZ", addLabs = NULL, mkDir = T
)

## Load in DE results for each major cell type (excluding cycling cells)
pwds <- lapply(levels(seu.obj$majorID_sub)[c(1,2,4:7)], function(x){paste0("../output/", outName, "/pseudoBulk/", x, "/", outName, "_cluster_", x, "_all_genes.csv")})
df.list <- lapply(pwds, read.csv, header = T)

## Derive tissue gene signatures
#Identify genes that are upregulated in TILs across all cell types 
feats.list <- lapply(df.list, function(x){
        feats <- x %>% filter(log2FoldChange > 0) %>% pull(gene)
})
tumor.sig <- Reduce(intersect, do.call(c, feat.list))
#Identify genes upregulated in blood
feats.list2 <- lapply(df.list, function(x){
        feats <- x %>% filter(log2FoldChange < 0) %>% pull(gene)
})
blood.sig <- Reduce(intersect, do.call(c, feat.list2))
#Stash the signature for later use
tissue.sig <- c(tumor.sig, blood.sig)
write.csv(tissue.sig, "./metaData/tumorSig.csv", row.names = F)

## Complete GSEA of tissue signatures
#Gene ontology
p <- plotGSEA(
        geneList = tumor.sig, category = "C5", species = "dog", termsTOplot = 10, 
        upOnly = T, pvalueCutoff = 0.05, subcategory = NULL
        ) + theme(axis.title = element_text(size = 16))
p <- p + 
        scale_x_continuous(
                limits = c(-12, ceiling(max(p$data$x_axis)*1.05)), 
                breaks = c(0, ceiling(max(p$data$x_axis)*1.05)/2, ceiling(max(p$data$x_axis)*1.05)),
                name = "log10(p.adj)"
        ) + 
        ggtitle("Gene ontology") + 
        theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste0("../output/", outName, "/gseaPlot_1.png"), width =7, height = 7)
#Reactome
p <- plotGSEA(
        geneList = tumor.sig, category = "C2", species = "dog", termsTOplot = 10, 
        upOnly = T, pvalueCutoff = 0.05, subcategory = "CP:REACTOME"
        ) + theme(axis.title=element_text(size = 16))
p <- p + 
        scale_x_continuous(
                limits = c(-20, ceiling(max(p$data$x_axis)*1.05)), 
                breaks = c(0, ceiling(max(p$data$x_axis)*1.05)/2, ceiling(max(p$data$x_axis)*1.05)),
                name = "log10(p.adj)"
        ) + 
        ggtitle("Reactome") + 
        theme(plot.title = element_text(size = 20, hjust = 0.5))
ggsave(paste0("../output/", outName, "/gseaPlot_2.png"), width = 7, height = 7)

## Fig supp 1a - upset plot demonstrating the conserved gene signature
#Format the input data
upSet.df <- as.data.frame(unique(tumor.sig))
colnames(upSet.df) <- "gene"
for (i in 1:6){
        upSet.df[ , levels(seu.obj$majorID_sub)[c(1,2,4:7)][i]] <- 
                as.integer(ifelse(upSet.df$gene %in% feats.list[i], 1, 0))
}
#Generate the UpSet plot
png(file = paste0("../output/", outName, "/upSet.png"), width=4000, height=2000, res=400)
par(mfcol=c(1,1))     
upset(upSet.df, sets = colnames(upSet.df)[2:7], , cutoff = 10,  nintersects = 60)
dev.off()

## Complete pseudobulk DE by each major cell type with tissue sig genes excluded
createPB(
        seu.obj = seu.obj, groupBy = "majorID_sub", comp = "cellSource", 
        biologicalRep = "orig.ident", min.cell = 5, clusters = NULL, 
        outDir = paste0("../output/", outName, "/pseudoBulk/"), 
        grepTerm = "tumor", grepLabel = c("TILs", "Blood"), 
        featsTOexclude = c(pal_feats, tissue.sig), lowFilter = T, dwnSam = F
)
p_volc <- pseudoDEG(
        metaPWD = paste0("../output/", outName, "/pseudoBulk/majorID_sub_deg_metaData.csv"),
        padj_cutoff = 0.01, lfcCut = 0.58, 
        outDir = paste0("../output/", outName, "/pseudoBulk/"), 
        outName = paste0(outName, "_FILTERED"), 
        idents.1_NAME = "TILs", idents.2_NAME = "Blood",
        inDir = paste0("../output/", outName, "/pseudoBulk/"), 
        title = "TILS vs Blood", 
        fromFile = T, meta = NULL, pbj = NULL, returnVolc = T, 
        paired = F, pairBy = "", 
        minimalOuts = F, saveSigRes = T, filterTerm = "^ENSCAF", addLabs = NULL, 
        mkDir = T
)
saveRDS(seu.obj, "../output/s3/bloodANDtils_filtered_allCells_S3.rds")

########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#######   end all cells analysis   ######## <<<<<<<<<<<<<<<<<
########################################### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


