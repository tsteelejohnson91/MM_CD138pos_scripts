
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggplot2)
library(infercnv)
`%notin%` <- Negate(`%in%`)

releveler = function(list){
  old_elem = as.numeric(as.character(unique(list)))
  old_elem = as.character(sort(old_elem))
  converter = c(1:length(old_elem))
  names(converter) = as.character(old_elem)
  new_elem = as.factor(converter[as.character(list)])
  return(new_elem)
}

FindClusterMarkers = function(X,g,n){
  MarkerTable = list()
  for(grp in unique(g)){
    Pvals = apply(X,2,function(x) getPval(x,g==grp))
    Log2FCs = apply(X,2,function(x) getLog2FC(x,g==grp))
    FDRs = p.adjust(Pvals,method="BH")
    NegLog10_Pvals = -log10(Pvals)
    NegLog10_FDRs = -log10(FDRs)
    FDR_Ranks = rank(FDRs)
    Log2FC_Ranks = rank(-Log2FCs)
    FDR_Log2FC_Ranksums = FDR_Ranks+Log2FC_Ranks
    MarkerTable[[grp]] = data.frame(Gene = names(Pvals),Cluster=rep(grp,length(Pvals)),Log2FC=Log2FCs,Pval=Pvals,FDR=FDRs,NegLog10_Pval=NegLog10_Pvals,
                                    NegLog10_FDR=NegLog10_FDRs,Log2FC_Rank=Log2FC_Ranks,FDR_Rank=FDR_Ranks,FDR_Log2FC_Ranksum=FDR_Log2FC_Ranksums)
    MarkerTable[[grp]] = arrange(MarkerTable[[grp]],FDR_Log2FC_Ranksum) %>% top_n(-n)
  }
  MarkerTable = do.call(rbind,MarkerTable)
  return(MarkerTable)
}

npc.seurat = readRDS("~/Desktop/Genentech_project/npc_seurat.rds")
qc_cutoffs <- read.csv("~/Desktop/Genentech_project/QCcutoffs.csv")
atac_qc_cutoffs <- readRDS("~/Desktop/Genentech_project/atacQCcutoffs.rds")
rownames(qc_cutoffs) = qc_cutoffs$sample
patID = "0661_1274_5961619"
inp.dir = "~/Desktop/MM_multiome/data/SMM/"
out.dir = paste0("~/Desktop/MM_multiome/indivSampleAnalysis/SMM/Sample_",patID,"/")
system(command=paste0("mkdir ",out.dir))

# Loading counts
counts <- Read10X(data.dir = paste0(inp.dir,"Sample_",patID,"/filtered_feature_bc_matrix/"))
# Initialize the Seurat object with the raw (non-normalized data).

# Create Seurat object
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
#newseqlevelsStyle <- mapSeqlevels(seqlevels(annotation), "UCSC")
#annotation <- renameSeqlevels(annotation, newseqlevelsStyle)
genome(annotation) <- "hg38"

# Cleaning counts object
bleh = counts$Peaks[sub("[:].*","",rownames(counts$Peaks)) %in% seqlevels(annotation),]

mm_cd138pos <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)
mm_cd138pos[["percent.mt"]] <- PercentageFeatureSet(mm_cd138pos, pattern = "^MT-")

frag.file <- paste0(inp.dir,"Sample_",patID,"/atac_fragments.tsv.gz")
mm_cd138pos[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  annotation = annotation
)

# Viewing the RNA and ATAC counts
#VlnPlot(mm_cd138pos, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
DefaultAssay(mm_cd138pos) <- "ATAC"
mm_cd138pos <- NucleosomeSignal(mm_cd138pos)
mm_cd138pos <- TSSEnrichment(mm_cd138pos)
if(paste0("Sample_",patID) %notin% rownames(atac_qc_cutoffs)){
  plot(VlnPlot(object = mm_cd138pos,features = c("nFeature_RNA", "nFeature_ATAC","percent.mt","TSS.enrichment","nucleosome_signal"),ncol = 5))
  pdf(paste0(out.dir,"MultiomeQCpre.pdf"),width=20,height=4)
  VlnPlot(object = mm_cd138pos,features = c("nFeature_RNA", "nFeature_ATAC","percent.mt","TSS.enrichment","nucleosome_signal"),ncol = 5)
  dev.off()
  atac_qc_cutoffs[paste0("Sample_",patID),"nFeature_ATAC_max"] <- readline(prompt = "Enter nFeatureATAC max cutoff: ")
  atac_qc_cutoffs[paste0("Sample_",patID),"nFeature_ATAC_min"] <- readline(prompt = "Enter nFeatureATAC min cutoff: ")
  atac_qc_cutoffs[paste0("Sample_",patID),"nucleosome_signal_max"] <- readline(prompt = "Enter nucleosome signal max cutoff: ")
  atac_qc_cutoffs[paste0("Sample_",patID),"nucleosome_signal_min"] <- readline(prompt = "Enter nucleosome signal min cutoff: ")
  atac_qc_cutoffs[paste0("Sample_",patID),"TSS.enrichment_max"] <- readline(prompt = "Enter TSS enrichment max cutoff: ")
  atac_qc_cutoffs[paste0("Sample_",patID),"TSS.enrichment_min"] <- readline(prompt = "Enter TSS enrichment min cutoff: ")
}
rnames_tmp = rownames(atac_qc_cutoffs)
atac_qc_cutoffs <- as.data.frame(lapply(atac_qc_cutoffs,as.numeric))
rownames(atac_qc_cutoffs) <- rnames_tmp
saveRDS(atac_qc_cutoffs,file="~/Desktop/Genentech_project/atacQCcutoffs.rds")

# Removing low quality cells
# Mitochondrial content mt% 
mm_cd138pos <- subset(
  x = mm_cd138pos,
  subset = nFeature_RNA < qc_cutoffs[paste0("Sample_",patID),"nFeature_RNA_max"] &
    nFeature_RNA > qc_cutoffs[paste0("Sample_",patID),"nFeature_RNA_min"] &
    percent.mt < qc_cutoffs[paste0("Sample_",patID),"percentMT_max"] &
    nFeature_ATAC < as.numeric(atac_qc_cutoffs[paste0("Sample_",patID),"nFeature_ATAC_max"]) &
    nFeature_ATAC > as.numeric(atac_qc_cutoffs[paste0("Sample_",patID),"nFeature_ATAC_min"]) &
    nucleosome_signal < as.numeric(atac_qc_cutoffs[paste0("Sample_",patID),"nucleosome_signal_max"]) &
    nucleosome_signal > as.numeric(atac_qc_cutoffs[paste0("Sample_",patID),"nucleosome_signal_min"]) &
    TSS.enrichment < as.numeric(atac_qc_cutoffs[paste0("Sample_",patID),"TSS.enrichment_max"]) &
    TSS.enrichment > as.numeric(atac_qc_cutoffs[paste0("Sample_",patID),"TSS.enrichment_min"])
)

# After plot
VlnPlot(object = mm_cd138pos,features = c("nFeature_RNA", "nFeature_ATAC","percent.mt","TSS.enrichment","nucleosome_signal"),ncol = 5)
pdf(paste0(out.dir,"MultiomeQCpost.pdf"),width=16,height=4)
VlnPlot(object = mm_cd138pos,features = c("nFeature_RNA", "nFeature_ATAC","percent.mt","TSS.enrichment","nucleosome_signal"),ncol = 5)
dev.off()

# RNA analysis
DefaultAssay(mm_cd138pos) <- "RNA"
mm_cd138pos <- SCTransform(mm_cd138pos, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(mm_cd138pos) <- "ATAC"
mm_cd138pos <- RunTFIDF(mm_cd138pos)
mm_cd138pos <- FindTopFeatures(mm_cd138pos, min.cutoff = 'q0')
mm_cd138pos <- RunSVD(mm_cd138pos)
mm_cd138pos <- RunUMAP(mm_cd138pos, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

mm_cd138pos <- FindMultiModalNeighbors(mm_cd138pos, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
mm_cd138pos <- RunUMAP(mm_cd138pos, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mm_cd138pos <- FindClusters(mm_cd138pos, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(mm_cd138pos, reduction = "umap.rna", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(mm_cd138pos, reduction = "umap.atac", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(mm_cd138pos, reduction = "wnn.umap", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN")
pdf(paste0(out.dir,"Multiome_scatter_allcells.pdf"),width=10,height=4)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Remove non-plasma cells
tmp = data.frame(SDC1 = as.numeric(mm_cd138pos@assays$RNA@counts["SDC1",]),
                 Cluster = mm_cd138pos$seurat_clusters)
summary = tmp %>%
  group_by(Cluster) %>%
  summarize(quant10 = quantile(SDC1, probs = 0.10),
            quant25 = quantile(SDC1, probs = 0.25), 
            quant50 = quantile(SDC1, probs = 0.50),
            quant75 = quantile(SDC1, probs = 0.75),
            quant90 = quantile(SDC1, probs = 0.90))

Cluster_remove = as.numeric(which(summary$quant75==0))

mm_cd138pos = subset(mm_cd138pos, subset = seurat_clusters %in% summary$Cluster[summary$quant75>0])

# RNA analysis
DefaultAssay(mm_cd138pos) <- "RNA"
mm_cd138pos <- SCTransform(mm_cd138pos, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(mm_cd138pos) <- "ATAC"
mm_cd138pos <- RunTFIDF(mm_cd138pos)
mm_cd138pos <- FindTopFeatures(mm_cd138pos, min.cutoff = 'q0')
mm_cd138pos <- RunSVD(mm_cd138pos)
mm_cd138pos <- RunUMAP(mm_cd138pos, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

mm_cd138pos <- FindMultiModalNeighbors(mm_cd138pos, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
mm_cd138pos <- RunUMAP(mm_cd138pos, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mm_cd138pos <- FindClusters(mm_cd138pos, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# Remove non-plasma cells
tmp = data.frame(SDC1 = as.numeric(mm_cd138pos@assays$RNA@counts["SDC1",]),
                 Cluster = mm_cd138pos$seurat_clusters)
summary = tmp %>%
  group_by(Cluster) %>%
  summarize(quant10 = quantile(SDC1, probs = 0.10),
            quant25 = quantile(SDC1, probs = 0.25), 
            quant50 = quantile(SDC1, probs = 0.50),
            quant75 = quantile(SDC1, probs = 0.75),
            quant90 = quantile(SDC1, probs = 0.90))

Cluster_remove = as.numeric(which(summary$quant75==0))
mm_cd138pos = subset(mm_cd138pos, subset = seurat_clusters %in% summary$Cluster[summary$quant75>0])

# releveling the clusters 1-k
mm_cd138pos$seurat_clusters_orig = mm_cd138pos$seurat_clusters
mm_cd138pos$seurat_clusters = releveler(mm_cd138pos$seurat_clusters_orig)
Idents(mm_cd138pos) = mm_cd138pos$seurat_clusters
p1 <- DimPlot(mm_cd138pos, reduction = "umap.rna", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(mm_cd138pos, reduction = "umap.atac", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(mm_cd138pos, reduction = "wnn.umap", label = TRUE, label.size = 5, repel = TRUE) + ggtitle("WNN")
pdf(paste0(out.dir,"Multiome_scatter_plasma.pdf"),width=10,height=4)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(paste0(out.dir,"Multiome_PHF19_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_PHF19")
dev.off()
pdf(paste0(out.dir,"Multiome_NSD2_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_NSD2")
dev.off()
pdf(paste0(out.dir,"Multiome_CCND1_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_CCND1")
dev.off()
pdf(paste0(out.dir,"Multiome_CCND3_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_CCND3")
dev.off()
pdf(paste0(out.dir,"Multiome_MAF_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_MAF")
dev.off()
pdf(paste0(out.dir,"Multiome_MAFB_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_MAFB")
dev.off()
pdf(paste0(out.dir,"Multiome_AHCYL1_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_AHCYL1")
dev.off()
pdf(paste0(out.dir,"Multiome_CKS1B_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_CKS1B")
dev.off()
pdf(paste0(out.dir,"Multiome_TP53_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_TP53")
dev.off()
pdf(paste0(out.dir,"Multiome_RB1_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_RB1")
dev.off()
pdf(paste0(out.dir,"Multiome_MYC_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_MYC")
dev.off()
pdf(paste0(out.dir,"Multiome_CDKN2C_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_CDKN2C")
dev.off()
pdf(paste0(out.dir,"Multiome_FGFR3_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_FGFR3")
dev.off()
pdf(paste0(out.dir,"Multiome_ITGB7_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_ITGB7")
dev.off()
pdf(paste0(out.dir,"Multiome_CCND2_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_CCND2")
dev.off()
pdf(paste0(out.dir,"Multiome_TNFRSF17_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_TNFRSF17")
dev.off()
pdf(paste0(out.dir,"Multiome_FCRL5_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_FCRL5")
dev.off()
pdf(paste0(out.dir,"Multiome_CD274_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_CD274")
dev.off()
pdf(paste0(out.dir,"Multiome_GPRC5D_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_GPRC5D")
dev.off()
pdf(paste0(out.dir,"Multiome_LAG3_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_LAG3")
dev.off()
pdf(paste0(out.dir,"Multiome_IL15_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_IL15")
dev.off()
pdf(paste0(out.dir,"Multiome_IL2RA_scatter_plasma.pdf"),width=4,height=4)
FeaturePlot(mm_cd138pos,reduction = "wnn.umap",features="sct_IL2RA")
dev.off()

tmp = data.frame(CD56 = as.numeric(mm_cd138pos@assays$RNA@counts["NCAM1",]),
                 Cluster = mm_cd138pos$seurat_clusters)
summaryCD56 = tmp %>%
  group_by(Cluster) %>%
  summarize(quant10 = quantile(CD56, probs = 0.10),
            quant25 = quantile(CD56, probs = 0.25), 
            quant50 = quantile(CD56, probs = 0.50),
            quant75 = quantile(CD56, probs = 0.75),
            quant90 = quantile(CD56, probs = 0.90))

tmp = data.frame(CD45 = as.numeric(mm_cd138pos@assays$RNA@counts["PTPRC",]),
                 Cluster = mm_cd138pos$seurat_clusters)
summaryCD45 = tmp %>%
  group_by(Cluster) %>%
  summarize(quant10 = quantile(CD45, probs = 0.10),
            quant25 = quantile(CD45, probs = 0.25), 
            quant50 = quantile(CD45, probs = 0.50),
            quant75 = quantile(CD45, probs = 0.75),
            quant90 = quantile(CD45, probs = 0.90))


deg_obj_tmp = mm_cd138pos
DefaultAssay(deg_obj_tmp) = "RNA"
deg_obj_tmp <- NormalizeData(object=deg_obj_tmp, normalization.method = "LogNormalize",
                             scale.factor = 1e4)
cluster.markers = FindAllMarkers(deg_obj_tmp,min.pct=0.25,min.diff.pct=0.25)
write.csv(cluster.markers,paste0(out.dir,"cluster_markers.csv"))
gne.targets = c("TNFRSF17","FCRL5","CD274","GPRC5D","LAG3","IL15","IL2RA")
gne.markers = FindAllMarkers(deg_obj_tmp,features=gne.targets,logfc.threshold=0,min.pct=0,min.cells.feature=0,min.cells.group=0,return.thresh=1)
write.csv(gne.markers,paste0(out.dir,"gne_markers.csv"))

# Infer CNV analysis
mm_cd138pos2 = subset(mm_cd138pos,subset = nFeature_RNA > quantile(mm_cd138pos$nFeature_RNA,0.5))
tumor_counts = mm_cd138pos2@assays$RNA@counts
npc_counts = npc.seurat@assays$RNA@counts
npc_counts = npc_counts[,npc_counts["CKS1B",]==0]
#npc_counts = npc.seurat@assays$RNA@counts[,npc.seurat$new.ident=="RRMM-11"]
int_feats = intersect(rownames(tumor_counts),rownames(npc_counts))
counts = cbind(npc_counts[int_feats,],tumor_counts[int_feats,])
cluster_var = factor(c(rep("NPC",dim(npc_counts)[2]),paste0("Cluster_",mm_cd138pos2$seurat_clusters)),
                     levels = c("NPC",paste0("Cluster_",1:length(unique(mm_cd138pos2$seurat_clusters)))))
counts = counts[,order(cluster_var)]
meta = data.frame(cluster = cluster_var[order(cluster_var)],row.names=colnames(counts))

# inferCNV with multiome clusters
out.dir2 = paste0(out.dir,"inferCNV/")
system(command=paste0("mkdir ",out.dir2))
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix=counts,
  annotations_file=meta,
  delim="\t",
  gene_order_file="~/Desktop/grants/LLS_grant/genome_annotations/grch38_gene_pos.txt",
  ref_group_names=c("NPC"),
  max_cells_per_group = 500)

copy_infercnv_obj = infercnv_obj
new_gene_order = data.frame()
for (chr_name in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) {
  new_gene_order = rbind(new_gene_order, infercnv_obj@gene_order[which(infercnv_obj@gene_order[["chr"]] == chr_name) , , drop=FALSE])
}
names(new_gene_order) <- c("chr", "start", "stop")
copy_infercnv_obj@gene_order = new_gene_order
copy_infercnv_obj@expr.data = infercnv_obj@expr.data[rownames(new_gene_order), , drop=FALSE]

infercnv_obj_default = infercnv::run(
  copy_infercnv_obj,
  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=out.dir2,
  cluster_by_groups=TRUE, 
  plot_steps=FALSE,
  denoise=TRUE,
  HMM=FALSE,
  no_prelim_plot=TRUE,
  png_res=300
)

# inferCNV with inferCNV clustering
meta2 = meta
meta2$cluster = as.character(meta2$cluster)
meta2$cluster[meta2$cluster!="NPC"] = "Myeloma"
meta2$cluster = factor(meta2$cluster,levels=c("NPC","Myeloma"))
out.dir3 = paste0(out.dir,"inferCNV2/")
system(command=paste0("mkdir ",out.dir3))
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix=counts,
  annotations_file=meta2,
  delim="\t",
  gene_order_file="~/Desktop/grants/LLS_grant/genome_annotations/grch38_gene_pos.txt",
  ref_group_names=c("NPC"),
  max_cells_per_group = 1000)

copy_infercnv_obj = infercnv_obj
new_gene_order = data.frame()
for (chr_name in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) {
  new_gene_order = rbind(new_gene_order, infercnv_obj@gene_order[which(infercnv_obj@gene_order[["chr"]] == chr_name) , , drop=FALSE])
}
names(new_gene_order) <- c("chr", "start", "stop")
copy_infercnv_obj@gene_order = new_gene_order
copy_infercnv_obj@expr.data = infercnv_obj@expr.data[rownames(new_gene_order), , drop=FALSE]

infercnv_obj_default = infercnv::run(
  copy_infercnv_obj,
  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=out.dir3,
  cluster_by_groups=FALSE, 
  plot_steps=FALSE,
  denoise=TRUE,
  HMM=FALSE,
  no_prelim_plot=TRUE,
  png_res=300
)




# Testing
if(FALSE){
infercnv_obj_default = infercnv::run(
  analysis_mode = "subcluster",
  copy_infercnv_obj,
  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=out.dir3,
  cluster_by_groups=FALSE,
  hclust_method='ward.D2',
  tumor_subcluster_partition_method="random_trees",
  tumor_subcluster_pval=0.05,
  plot_steps=FALSE,
  denoise=TRUE,
  HMM=TRUE,
  no_prelim_plot=TRUE,
  png_res=300
)
}


