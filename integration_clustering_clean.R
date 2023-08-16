###############################################
# This script integrates all of the RNA from
# the multiome samples together into a single
# representation.
# Travis S. Johnson
# February 7, 2022

library(dplyr)
library(Seurat)
library(ggplot2)

releveler = function(list){
  old_elem = as.numeric(as.character(unique(list)))
  old_elem = as.character(sort(old_elem))
  converter = c(1:length(old_elem))
  names(converter) = as.character(old_elem)
  new_elem = as.factor(converter[as.character(list)])
  return(new_elem)
}
vers.id = "v3_2"
data.path <- "/N/project/phi_ingest_mm/Brian_Walker/Multiome_analysis/"
proc.path <- "/N/u/johnstrs/Carbonate/Multiome_data/"

qc_cutoffs <- read.csv(paste0(proc.path,"QCcutoffs.csv"))
rownames(qc_cutoffs) = qc_cutoffs$sample

sample_info = read.csv(paste0(proc.path,"sample_info.csv"))
samples = sample_info$samples
diagnosis = sample_info$diagnosis

sample.list <- list()
QC_metrics <- matrix(NA,length(samples),3)
rownames(QC_metrics) = samples
colnames(QC_metrics) = c("preQC_cells","postQC_cells","plasma_cells")
QC_metrics = as.data.frame(QC_metrics)

options(warn = -1)
for(s in samples){
message(s)
tmp <- Read10X(data.dir=paste0(data.path,s,"/outs/filtered_feature_bc_matrix"))
tmp <- CreateSeuratObject(counts=tmp$`Gene Expression`, min.cells = 10, min.features=100)
tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")

QC_metrics[s,"preQC_cells"] = dim(tmp)[2]
pdf(paste0(proc.path,s,"_preQC_violinplot_v3.pdf"),width=21,height=7)
plot(VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()

tmp <- subset(
  x = tmp,
  subset = nFeature_RNA < qc_cutoffs[s,"nFeature_RNA_max"] &
           nFeature_RNA > qc_cutoffs[s,"nFeature_RNA_min"] &
           percent.mt < qc_cutoffs[s,"percentMT_max"]
  )

QC_metrics[s,"postQC_cells"] = dim(tmp)[2]
pdf(paste0(proc.path,s,"_postQC_violinplot_v3.pdf"),width=21,height=7)
plot(VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()

sample.list[[s]] <- tmp
rm(tmp)
message(paste0(dim(sample.list[[s]])[1]," ",dim(sample.list[[s]])[2]))
}

##################################################
# Integration

for (s in samples){
  message(s)
  sample.list[[s]]$orig.ident <- rep(s,dim(sample.list[[s]])[2])
  sample.list[[s]] <- FindVariableFeatures(sample.list[[s]])
}
features <- SelectIntegrationFeatures(object.list = sample.list)
sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = sample.list, reduction = "rpca",
                                  dims = 1:50)
samples.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
samples.integrated <- ScaleData(samples.integrated, verbose = FALSE)
samples.integrated <- RunPCA(samples.integrated, verbose = FALSE)
samples.integrated <- RunUMAP(samples.integrated, dims = 1:50)
samples.integrated <- FindNeighbors(samples.integrated)
samples.integrated <- FindClusters(samples.integrated)
rm(sample.list,features,anchors)

###############################################
# Removing non-plasma cells

tmp = data.frame(SDC1 = as.numeric(samples.integrated@assays$RNA["SDC1",]),
                 Cluster = samples.integrated$seurat_clusters)
summary = tmp %>%
  group_by(Cluster) %>%
  summarize(quant10 = quantile(SDC1, probs = 0.10),
            quant25 = quantile(SDC1, probs = 0.25), 
            quant50 = quantile(SDC1, probs = 0.50),
            quant75 = quantile(SDC1, probs = 0.75),
            quant90 = quantile(SDC1, probs = 0.90))

samples.integrated <- subset(
  x = samples.integrated,
  subset = seurat_clusters %in% summary$Cluster[summary$quant90>0]
)

###############################################
# Reintegrating the samples

sample.list = SplitObject(samples.integrated, split.by = "orig.ident")
rm(samples.integrated)
for (s in samples){
  message(s)
  sample.list[[s]]$orig.ident <- rep(s,dim(sample.list[[s]])[2])
  sample.list[[s]] <- FindVariableFeatures(sample.list[[s]])
}
features <- SelectIntegrationFeatures(object.list = sample.list)
sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = sample.list, reduction = "rpca",
                                  dims = 1:50)
samples.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
samples.integrated <- ScaleData(samples.integrated, verbose = FALSE)
samples.integrated <- RunPCA(samples.integrated, verbose = FALSE)
samples.integrated <- RunUMAP(samples.integrated, dims = 1:50)
samples.integrated <- FindNeighbors(samples.integrated)
samples.integrated <- FindClusters(samples.integrated)
rm(sample.list,features,anchors)

tmp = data.frame(SDC1 = as.numeric(samples.integrated@assays$RNA["SDC1",]),
                 Cluster = samples.integrated$seurat_clusters)
summary2 = tmp %>%
  group_by(Cluster) %>%
  summarize(quant10 = quantile(SDC1, probs = 0.10),
            quant25 = quantile(SDC1, probs = 0.25), 
            quant50 = quantile(SDC1, probs = 0.50),
            quant75 = quantile(SDC1, probs = 0.75),
            quant90 = quantile(SDC1, probs = 0.90))

samples.integrated <- subset(
  x = samples.integrated,
  subset = seurat_clusters %in% summary2$Cluster[summary2$quant75>0]
)

mapper = diagnosis
names(mapper) = samples
samples.integrated$diagnosis = mapper[samples.integrated$orig.ident]
samples.integrated$seurat_clusters_orig = samples.integrated$seurat_clusters
samples.integrated$seurat_clusters = releveler(samples.integrated$seurat_clusters_orig)
Idents(samples.integrated) = samples.integrated$seurat_clusters

diagTable = table(samples.integrated$seurat_clusters,samples.integrated$diagnosis)
diagTable = as.data.frame.matrix(diagTable)
diagTable = diagTable[,c("SMM","NDMM","RRMM")]
write.csv(diagTable,file=paste0(proc.path,"diagnosis_table_counts_",vers.id,".csv"),quote=FALSE)
diagTable = diagTable/rowSums(diagTable)
write.csv(diagTable,file=paste0(proc.path,"diagnosis_table_percent_",vers.id,".csv"),quote=FALSE)

for (s in samples){QC_metrics[s,"plasma_cells"] = sum(samples.integrated$orig.ident==s)}
write.csv(QC_metrics,file=paste0(proc.path,"qc_metrics_",vers.id,".csv"))

saveRDS(samples.integrated,file=paste0(proc.path,"samples_integrated_",vers.id,".rds"))
