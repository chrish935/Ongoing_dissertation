## Made by Christopher M. Horn, MS
## Kielian Lab data
## Analyzing orthopedic model scRNA-seq data over the course of multiple days (D03, D07, D14)
## Created: 2020-05-12
## Updated: 2022-02-03




# Setting up environment -----

## Load in packages
packages <- c('tidyverse',
              'Seurat',
              'patchwork',
              'SingleR',
              'SingleCellExperiment',
              'sctransform',
              'MAST',
              'plotly',
              'ggsci',
              'readxl',
              'fgsea',
              'data.table',
              'dyno',
              'clusterExperiment',
              'pheatmap',
              'mgcv',
              'slingshot',
              'tradeSeq',
              'UpSetR',
              'msigdbr',
              'celldex',
              'escape',
              'dittoSeq',
              'here')

invisible(lapply(packages, library, character.only = T))

rm(packages)

## Set options
options(future.globals.maxSize = 4000 * 1024^2)
set.seed(12345)



# Initializing objects -----

## Setting up Seurat objects
## Load in the data
d3.data <- Read10X(data.dir = here('D3', 'Data'))
d7.data <- Read10X(data.dir = here('D7', 'Data'))
d14.data <- Read10X(data.dir = here('D14', 'Data'))

## Initialize Seurat object w/raw data
d3.ortho <- CreateSeuratObject(counts = d3.data, project = 'D03', min.cells = 3, min.features = 200)
d7.ortho <- CreateSeuratObject(counts = d7.data, project = 'D07', min.cells = 3, min.features = 200)
d14.ortho <- CreateSeuratObject(counts = d14.data, project = 'D14', min.cells = 3, min.features = 200)

## Pre-processing and QC
d3.ortho[['percent.mt']] <- PercentageFeatureSet(d3.ortho, pattern = '^mt-')
d3.ortho <- subset(d3.ortho, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
d7.ortho[['percent.mt']] <- PercentageFeatureSet(d7.ortho, pattern = '^mt-')
d7.ortho <- subset(d7.ortho, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
d14.ortho[['percent.mt']] <- PercentageFeatureSet(d14.ortho, pattern = '^mt-')
d14.ortho <- subset(d14.ortho, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 5)

## Create timing metadata
sample.id_1 <- rep('D03', length(d3.ortho@meta.data$orig.ident))
sample.id_2 <- rep('D07', length(d7.ortho@meta.data$orig.ident))
sample.id_3 <- rep('D14', length(d14.ortho@meta.data$orig.ident))
names(sample.id_1) <- rownames(d3.ortho@meta.data)
names(sample.id_2) <- rownames(d7.ortho@meta.data)
names(sample.id_3) <- rownames(d14.ortho@meta.data)
d3.ortho <- AddMetaData(d3.ortho, sample.id_1, col.name = 'Sample.ID')
d7.ortho <- AddMetaData(d7.ortho, sample.id_2, col.name = 'Sample.ID')
d14.ortho <- AddMetaData(d14.ortho, sample.id_3, col.name = 'Sample.ID')

## Remove temp objects
rm(d3.data,
   d7.data,
   d14.data,
   sample.id_1,
   sample.id_2,
   sample.id_3)

gc()



# Annotating cell types -----

## Annotate the cells w/SingleR
immgen.se <- ImmGenData()

d3.ortho_sce <- as.SingleCellExperiment(d3.ortho)
d7.ortho_sce <- as.SingleCellExperiment(d7.ortho)
d14.ortho_sce <- as.SingleCellExperiment(d14.ortho)
commonGenes.1 <- intersect(rownames(d3.ortho_sce), rownames(immgen.se))
commonGenes.2 <- intersect(rownames(d7.ortho_sce), rownames(immgen.se))
commonGenes.3 <- intersect(rownames(d14.ortho_sce), rownames(immgen.se))
immgen.se.1 <- immgen.se[commonGenes.1,]
immgen.se.2 <- immgen.se[commonGenes.2,]
immgen.se.3 <- immgen.se[commonGenes.3,]
d3.ortho_sce <- d3.ortho_sce[commonGenes.1,]
d7.ortho_sce <- d7.ortho_sce[commonGenes.2,]
d14.ortho_sce <- d14.ortho_sce[commonGenes.3,]
pred.d3.ortho <- SingleR(test = d3.ortho_sce, ref = immgen.se.1, labels = immgen.se.1$label.main)
pred.d7.ortho <- SingleR(test = d7.ortho_sce, ref = immgen.se.2, labels = immgen.se.2$label.main)
pred.d14.ortho <- SingleR(test = d14.ortho_sce, ref = immgen.se.3, labels = immgen.se.3$label.main)
d3.ortho[['celltype']] <- pred.d3.ortho$pruned.labels
d7.ortho[['celltype']] <- pred.d7.ortho$pruned.labels
d14.ortho[['celltype']] <- pred.d14.ortho$pruned.labels

## Remove temp objects
rm(immgen.se,
   immgen.se.1,
   immgen.se.2,
   immgen.se.3,
   d3.ortho_sce,
   d7.ortho_sce,
   d14.ortho_sce,
   commonGenes.1,
   commonGenes.2,
   commonGenes.3,
   pred.d3.ortho,
   pred.d7.ortho,
   pred.d14.ortho)

gc()



# Object integration -----

## Integrating all cells from all days
ortho.list <- c(d3.ortho, d7.ortho, d14.ortho)
names(ortho.list) <- c('D03', 'D07', 'D14')
for (i in 1:length(ortho.list)) {ortho.list[[i]] <- SCTransform(ortho.list[[i]], verbose = T)}
ortho.features <- SelectIntegrationFeatures(object.list = ortho.list, nfeatures = 3000)
ortho.list <- PrepSCTIntegration(object.list = ortho.list, anchor.features = ortho.features, verbose = T)
ortho.anchors <- FindIntegrationAnchors(object.list = ortho.list, normalization.method = 'SCT', anchor.features = ortho.features, verbose = T)
ortho.integrated <- IntegrateData(anchorset = ortho.anchors, normalization.method = 'SCT', verbose = T)
ortho.integrated <- RunPCA(ortho.integrated, verbose = T)
ortho.integrated <- RunUMAP(ortho.integrated, dims = 1:30)
ortho.integrated <- FindNeighbors(ortho.integrated, dims = 1:30)
ortho.integrated <- FindClusters(ortho.integrated, resolution = 0.5)



## Output celltype composition of each cluster
for (i in seq_along(levels(ortho.integrated@meta.data$seurat_clusters))) {
  if(i != length(ortho.integrated@meta.data$seurat_clusters)) {
    ii <- ortho.integrated@meta.data %>% filter(seurat_clusters == i-1) %>% select(celltype) %>% table()
    write.csv(ii, here('Integrated', 'Cluster Composition', 'Ortho Celltype', paste0('cluster_', i-1, '.csv')))
  }
}

rm(i,
   ii)

gc()



## Plotting
ortho.int_pca.time <- DimPlot(ortho.integrated, split.by = 'Sample.ID', reduction = 'pca') +
  theme_classic() +
  NoLegend() + ggtitle('Split by Sample ID') + labs(x = 'PC 1', y = 'PC 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ortho.int_pca.clus <- DimPlot(ortho.integrated, reduction = 'pca') +
  scale_color_d3(palette = 'category20b') +
  theme_classic() +
  ggtitle('', subtitle = 'Colored by Cluster') +
  labs(x = 'PC 1', y = 'PC 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ortho.int_umap.time <- DimPlot(ortho.integrated, split.by = 'Sample.ID') +
  theme_classic() +
  NoLegend() +
  ggtitle('Split by Sample ID') +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ortho.int_umap.clus <- DimPlot(ortho.integrated) +
  scale_color_d3(palette = 'category20b') +
  theme_classic() +
  ggtitle('', subtitle = 'Colored by Cluster') +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

## Remove temp objects
rm(ortho.list,
   ortho.features,
   ortho.anchors,
   d3.ortho,
   d7.ortho,
   d14.ortho)

gc()



write.csv(ortho.integrated@meta.data, here('Integrated', 'ortho.int_metadata.csv'))



# Renaming/plotting -----

## Renaming clusters
# new.cluster.ids <- c('Granulocytes I', 'Granulocytes II', 'Granulocytes III', 'Granulocytes IV', 'Granulocytes V', 'Granulocytes VI', 'Granulocytes VII', '7', 'Granulocytes VIII', 'Granulocytes IX', 'Granulocytes X', '11') # Preparing new cluster labels
new.cluster.ids <- c('G1',
                     'G2',
                     'G3',
                     'G4',
                     'G5',
                     'G6',
                     'G7',
                     'C7',
                     'G8',
                     'G9',
                     'G10',
                     'C11') # Preparing new cluster labels
names(new.cluster.ids) <- levels(ortho.integrated)
ortho.integrated <- RenameIdents(ortho.integrated, new.cluster.ids)

rm(new.cluster.ids)

gc()

## Re-Plotting w/new labels
ortho.int_pca.time <- DimPlot(ortho.integrated, split.by = 'Sample.ID', reduction = 'pca') +
  theme_classic() +
  ggtitle('Split by Sample ID') +
  labs(x = 'PC 1', y = 'PC 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ortho.int_pca.clus <- DimPlot(ortho.integrated, reduction = 'pca') +
  theme_classic() +
  ggtitle('Colored by Cluster') +
  labs(x = 'PC 1', y = 'PC 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ortho.int_umap.time <- DimPlot(ortho.integrated, split.by = 'Sample.ID') +
  theme_classic() +
  ggtitle('Split by Sample ID') +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ortho.int_umap.clus <- DimPlot(ortho.integrated) +
  theme_classic() +
  ggtitle('Colored by Cluster') +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())



## Subsetting out the granulocytes
# neutro.integrated <- subset(ortho.integrated, idents = c('Granulocytes I', 'Granulocytes II', 'Granulocytes III', 'Granulocytes IV', 'Granulocytes V', 'Granulocytes VI', 'Granulocytes VII', 'Granulocytes VIII', 'Granulocytes IX', 'Granulocytes X'))
neutro.integrated <- subset(ortho.integrated, idents = c('G1',
                                                         'G2',
                                                         'G3',
                                                         'G4',
                                                         'G5',
                                                         'G6',
                                                         'G7',
                                                         'G8',
                                                         'G9',
                                                         'G10'))



## Output number of cells from each Sample.ID in each cluster
gran_clust.comp <- table(Idents(neutro.integrated), neutro.integrated$Sample.ID) # easier way

write.csv(gran_clust.comp, here('Integrated', 'Cluster Composition', 'Neutro SampleID', 'gran_clust comp.csv'))

rm(gran_clust.comp)

gc()



# Trajectory Analysis -----

# Finding trajectories =====

## Inferring trajectories
object_counts <- Matrix::t(as(as.matrix(neutro.integrated@assays$RNA@counts), 'sparseMatrix'))
object_expression <- Matrix::t(as(as.matrix(neutro.integrated@assays$RNA@data), 'sparseMatrix'))
neutro.integrated_dyn <- wrap_expression(
  counts = object_counts, 
  expression = object_expression
)

rm(object_counts, object_expression)

## Add a dimensionality reduction
neutro.integrated_dimred <- dyndimred::dimred_umap(neutro.integrated_dyn$expression)

## Infer the trajectory
neutro.integrated_model <- infer_trajectory(neutro.integrated_dyn, ti_slingshot(), verbose = T)

## Plot trajectory & pseudotime
neutro.integrated_milestone.umap <- plot_dimred(neutro.integrated_model, label_milestones = T, dimred = neutro.integrated_dimred, hex_cells = F) +
  theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = 'UMAP 1', y = 'UMAP 2') ## Check this prior to rooting
neutro.integrated_milestone.pca <- plot_dimred(neutro.integrated_model, label_milestones = T, dimred = 'pca', hex_cells = F) +
  theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = 'PC 1', y = 'PC 2') ## Check this prior to rooting

neutro.integrated_traj.umap <- plot_dimred(neutro.integrated_model, dimred = neutro.integrated_dimred, grouping = neutro.integrated@active.ident, color_density = 'grouping', hex_cells = F) +
  theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = 'UMAP 1', y = 'UMAP 2')
neutro.integrated_traj.pca <- plot_dimred(neutro.integrated_model, dimred = 'pca', grouping = neutro.integrated@active.ident, color_density = 'grouping', hex_cells = F) +
  theme_classic() +
  theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = 'PC 1', y = 'PC 2')
neutro.integrated_traj2.umap <- plot_dimred(neutro.integrated_model, dimred = neutro.integrated_dimred, grouping = neutro.integrated@meta.data$Sample.ID, color_density = 'grouping', hex_cells = F) +
  theme_classic() +
  theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = 'UMAP 1', y = 'UMAP 2')
neutro.integrated_traj2.pca <- plot_dimred(neutro.integrated_model, dimred = 'pca', grouping = neutro.integrated@meta.data$Sample.ID, color_density = 'grouping', hex_cells = F) +
  theme_classic() +
  theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(x = 'PC 1', y = 'PC 2')
neutro.integrated_pseudo.umap <- plot_dimred(neutro.integrated_model, "pseudotime", pseudotime = calculate_pseudotime(neutro.integrated_model), dimred = neutro.integrated_dimred, hex_cells = F) +
  theme_classic() +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
neutro.integrated_pseudo.pca <- plot_dimred(neutro.integrated_model, "pseudotime", pseudotime = calculate_pseudotime(neutro.integrated_model), dimred = 'pca', hex_cells = F) +
  theme_classic() +
  labs(x = 'PC 1', y = 'PC 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

## Root trajectory if necessary
neutro.integrated_model <- add_root(neutro.integrated_model, root_milestone_id = "4")

## Simplify trajectory
simp <- simplify_trajectory(neutro.integrated_model)
simp_lab <- simp %>% label_milestones(c('1' = 'Stuck',
                                        '3' = 'Bifurcation',
                                        '4' = 'Mature',
                                        '5' = 'Immature'))



# Plotting gene expression over pseudotime =====

# Slingshot #####

## Get trajectory & clustering information for lineage acquisition (pseudotime)
expression <- Matrix::t(as(as.matrix(neutro.integrated@assays$RNA@data), 'sparseMatrix'))
ndim <- 20L
max_clusters <- min(nrow(expression)-1, 10)
pca <- irlba::prcomp_irlba(expression, n = ndim)
  
  # Select optimal number of dimensions if ndim is large enough
  if (ndim > 3) {
    # This code is adapted from the expermclust() function in TSCAN
    # The only difference is in how PCA is performed
    # (they specify scale. = TRUE and we leave it as FALSE)
    x <- 1:ndim
    optpoint1 <- which.min(sapply(2:10, function(i) {
      x2 <- pmax(0, x - i)
      sum(lm(pca$sdev[1:ndim] ~ x + x2)$residuals^2 * rep(1:2,each = 10))
    }))
    
    # This is a simple method for finding the "elbow" of a curve, from
    # https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
    x <- cbind(1:ndim, pca$sdev[1:ndim])
    line <- x[c(1, nrow(x)),]
    proj <- princurve::project_to_curve(x, line)
    optpoint2 <- which.max(proj$dist_ind)-1
    
    # We will take more than 3 PCs only if both methods recommend it
    optpoint <- max(c(min(c(optpoint1, optpoint2)), 3))
  }

dimred <- pca$x[, seq_len(optpoint)]
rownames(dimred) <- rownames(expression)

clusterings <- lapply(3:max_clusters, function(K){
     cluster::pam(dimred, K) # we generally prefer PAM as a more robust alternative to k-means
  })
wh.cl <- which.max(sapply(clusterings, function(x){ x$silinfo$avg.width })) + 1
labels <- clusterings[[min(c(wh.cl, 8))]]$clustering

rm(ndim,
   optpoint,
   optpoint1,
   optpoint2,
   x,
   pca,
   proj,
   line,
   expression,
   clusterings,
   wh.cl,
   max_clusters)

gc()



## Find trajectory
lineages <- getLineages(dimred, labels, start.clus = '4')
sling <- getCurves(lineages,
                   shrink = 1L,
                   reweight = T,
                   reassign = T,
                   thresh = 0.001,
                   maxit = 10L,
                   stretch = 2L,
                   smoother = 'smooth.spline',
                   shrink.method = 'cosine') ## Parameters taken from dynverse ti_slingshot code

rm(dimred, labels, lineages)

gc()



# TradeSeq #####

## Find optimal number of knots to use for GAM fitting
counts <- neutro.integrated@assays$RNA@counts
filt <- rowSums(counts > 1) >= 120 # Filtering transcripts (genes w/count of at least two in at least 120 different cells)
counts.filt <- counts[filt, ]
rm(filt)

icMat <- evaluateK(counts = as.matrix(counts.filt), sds = sling, k = 3:15, nGenes = 300, verbose = T, plot = T)

## Fit GAM
sce <- fitGAM(as.matrix(counts.filt), sling, nknots = 10) # SCE method
converge <- table(rowData(sce)$tradeSeq$converged) # Check whether genes converged

# List method
# control <- gam.control()
# control$maxit <- 1000 # Set maximum number of iterations to 1,000
# gamList <- fitGAM(counts = as.matrix(counts.filt), pseudotime = slingPseudotime(sling, na = F), cellWeights = slingCurveWeights(sling), control = control, sce = F)
# pvalLineage <- getSmootherPvalues(gamList)
# statLineage <- getSmootherTestStats(gamList)

## Testing
assoRes <- associationTest(sce, lineages = T) # Testing whether genes are significantly changed along pseudotime (independent lineages)
startRes <- startVsEndTest(sce, lineages = T) # Testing whether genes are significantly changed between the start & end of pseudotime (independent lineages)
endRes <- diffEndTest(sce) # Testing whether genes are significantly changed at end of pseudotime
patternRes <- patternTest(sce) # Testing whether genes have significantly different expression patterns throughout pseudotime
earlyDE.plot <- plotGeneCount(curve = sling, counts = as.matrix(counts.filt), clusters = apply(slingClusterLabels(sling), 1, which.max), models = sce) # Visualize where knots are
earlyDERes <- earlyDETest(sce, knots = c(1, 3))

## Identifying significant changes in assoRes
assoRes.sig.lin1 <- rownames(assoRes)[which(p.adjust(assoRes$pvalue_1, "fdr") <= 0.05)]
assoRes.sig.lin2 <- rownames(assoRes)[which(p.adjust(assoRes$pvalue_2, "fdr") <= 0.05)]

assoRes.sig.upset <- upset(fromList(list('Lineage 1' = assoRes.sig.lin1 , 'Lineage 2' = assoRes.sig.lin2)))

yhatSmooth.1 <- predictSmooth(sce, gene = assoRes.sig.lin1, nPoints = 50)
heatSmooth.1 <- pheatmap(t(scale(t(yhatSmooth.1[, 1:50]))), cluster_cols = F, show_rownames = F, show_colnames = F)

yhatSmooth.2 <- predictSmooth(sce, gene = assoRes.sig.lin2, nPoints = 50)
heatSmooth.2 <- pheatmap(t(scale(t(yhatSmooth.2[, 1:50]))), cluster_cols = F, show_rownames = F, show_colnames = F)

## Combining End & Pattern testing
patternRes$Gene <- rownames(patternRes)
patternRes$pattern <- patternRes$waldStat
patternRes.comp <- patternRes[, c('Gene', 'pattern')]

endRes$Gene <- rownames(endRes)
endRes$end <- endRes$waldStat
endRes.comp <- endRes[, c('Gene', 'end')]

compare <- merge(patternRes.comp, endRes.comp, by = 'Gene', all = F)
compare$transientScore <- rank(-compare$end, ties.method = 'min')^2 + rank(compare$pattern, ties.method = 'random')^2

rm(patternRes.comp, endRes.comp)

## Plotting
## Expression vs Pseudotime
plotSmoothers(sce, assays(sce)$counts, gene = 'Ltf')

## End vs Pattern
ggplot(compare, aes(x = log(pattern), y = log(end))) +
  geom_point(aes(col = transientScore)) +
  labs(x = 'patternTest Wald Statistic (log scale)',
       y = 'diffEndTest Wald Statistic (log scale)') +
  scale_color_continuous(low = 'grey', high = 'blue') +
  theme_classic()

## Clustering genes with similar expression along pseudotime
nPointsClus <- 50
clusPat <- clusterExpressionPatterns(sce, nPoints = nPointsClus, genes = rownames(as.matrix(counts.filt)))
clusterLabels <- primaryCluster(clusPat$rsec)
cUniq <- unique(clusterLabels)
cUniq <- cUniq[!cUniq == -1] # remove unclustered genes

gene <- rownames(as.matrix(counts.filt))
Gene.Cluster <- as.data.frame(gene)
Gene.Cluster$cluster <- primaryCluster(clusPat$rsec)
Gene.Cluster <- arrange(Gene.Cluster, cluster)
write.csv(Gene.Cluster, here('Integrated', 'TradeSeq', 'Real time', 'gene cluster', 'full_cluster.csv'))

for (i in seq_along(unique(Gene.Cluster$cluster))) {
  if(i != length(Gene.Cluster$cluster)) {
    clust <- Gene.Cluster %>% filter(Gene.Cluster$cluster %in% c(unique(Gene.Cluster$cluster)[i]))
    write.csv(clust, here('Integrated', 'TradeSeq', 'Real time', 'gene cluster', paste0('cluster_', i-1, '.csv')))
  }
}

rm(clust,
   i)

gc()

for (xx in sort(cUniq)[1:6]) {
  cId <- which(clusterLabels == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
    y = rep(range(clusPat$yhatScaled[cId, ]),
    nPointsClus / 2)),
    aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0('Cluster ', xx),  x = 'Pseudotime', y = 'Normalized expression') +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 2),
                                  y = clusPat$yhatScaled[geneId, ],
                                  lineage = rep(0:1, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = FALSE) +
    scale_color_manual(values = c('#440154FF', '#FDE725FF'),
                       breaks = c('0', '1'))
  assign(paste0('c', xx), p)
}

rm(geneId,
   ii,
   xx,
   p,
   cId,
   nPointsClus,
   cUniq,
   clusterLabels)

gc()



# TradeSeq GSEA #####

## fgsea
geneSets <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
# geneSets <- msigdbr(species = 'Mus musculus', category = 'C2', subcategory = 'CP:KEGG') # KEGG
# geneSets <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark
geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
stats.lin1 <- assoRes$waldStat_1
stats.lin2 <- assoRes$waldStat_2
names(stats.lin1) <- toupper(rownames(assoRes))
names(stats.lin2) <- toupper(rownames(assoRes))
eaRes.lin1 <- fgsea(pathways = m_list, stats = stats.lin1, eps = 0.0, scoreType = 'pos', minSize = 10)
eaRes.lin2 <- fgsea(pathways = m_list, stats = stats.lin2, eps = 0.0, scoreType = 'pos', minSize = 10)
eaRes.lin1 <- filter(eaRes.lin1, padj <= 0.05)
eaRes.lin2 <- filter(eaRes.lin2, padj <= 0.05)
eaRes.lin1_unique <- anti_join(eaRes.lin1, eaRes.lin2, by = 'pathway')
eaRes.lin2_unique <- anti_join(eaRes.lin2, eaRes.lin1, by = 'pathway')
fwrite(eaRes.lin1, file = here('Integrated', 'Pseudotime', 'Real time', 'TradeSeq', 'fgsea', 'GO_lineage 1 enrichment analysis.tsv'), sep="\t", sep2=c("", " ", ""))
fwrite(eaRes.lin2, file = here('Integrated', 'Pseudotime', 'Real time', 'TradeSeq', 'fgsea', 'GO_lineage 2 enrichment analysis.tsv'), sep="\t", sep2=c("", " ", ""))
fwrite(eaRes.lin1_unique, file = here('Integrated', 'Pseudotime', 'Real time', 'TradeSeq', 'fgsea', 'GO_unique lineage 1 enrichment analysis.tsv'), sep="\t", sep2=c("", " ", ""))
fwrite(eaRes.lin2_unique, file = here('Integrated', 'Pseudotime', 'Real time', 'TradeSeq', 'fgsea', 'GO_unique lineage 2 enrichment analysis.tsv'), sep="\t", sep2=c("", " ", ""))

rm(geneSets,
   m_list,
   stats.lin1,
   stats.lin2)

gc()



# Cluster-level DE -----

## Gene DE Analysis
DefaultAssay(neutro.integrated) <- 'RNA'
neutro.integrated <- NormalizeData(neutro.integrated, verbose = T) # Counts are normalized by dividing the counts for a given feature by the total counts per cell, multiplying by a scale factor (default == 10,000), and then taking the natural log using log1p()
neutro.markers <- FindAllMarkers(neutro.integrated, min.pct = 0, logfc.threshold = 0, test.use = 'MAST')
neutro.filtered <- filter(neutro.markers, neutro.markers$p_val_adj <= 0.05)
write.csv(neutro.filtered, here('Integrated', 'DE', 'full DE_gran.csv'))

for (i in seq_along(levels(neutro.filtered$cluster))) {
  if(i != length(neutro.filtered$cluster)) {
    DE <- neutro.filtered %>% filter(neutro.filtered$cluster %in% c(levels(neutro.filtered$cluster)[i]))
    write.csv(DE, here('Integrated', 'DE', paste0('DE_gran', i, '.csv')))
  }
}

rm(DE,
   i)

gc()



markers.define <- subset(neutro.markers, avg_logFC > 1 & p_val_adj < 0.05)$gene

for (i in seq_along(levels(neutro.markers$cluster))) {
  if(i != length(neutro.markers$cluster)) {
    DE <- neutro.markers %>% filter(neutro.markers$cluster %in% c(levels(neutro.markers$cluster)[i])) %>% filter(avg_logFC > 1 & p_val_adj < 0.05)
    write.csv(DE, here('Integrated', 'DE', paste0('markers.define_gran', i, '.csv')))
  }
}

rm(DE,
   i)

gc()

## Within cluster DE
g1 <- subset(neutro.integrated, idents = c('G1'))
g2 <- subset(neutro.integrated, idents = c('G2'))
g3 <- subset(neutro.integrated, idents = c('G3'))
g4 <- subset(neutro.integrated, idents = c('G4'))
g5 <- subset(neutro.integrated, idents = c('G5'))
g6 <- subset(neutro.integrated, idents = c('G6'))
g7 <- subset(neutro.integrated, idents = c('G7'))
g8 <- subset(neutro.integrated, idents = c('G8'))
g9 <- subset(neutro.integrated, idents = c('G9'))
g10 <- subset(neutro.integrated, idents = c('G10'))

g.list <- list(g1,
               g2,
               g3,
               g4,
               g5,
               g6,
               g7,
               g8,
               g9,
               g10)

rm(g1,
   g2,
   g3,
   g4,
   g5,
   g6,
   g7,
   g8,
   g9,
   g10)

gc()

for(i in 1:length(g.list)){
  DefaultAssay(g.list[[i]]) <- 'RNA'
  g.list[[i]] <- NormalizeData(g.list[[i]], verbose = T)
  g.list[[i]] <- ScaleData(g.list[[i]], verbose = T)
  g.list[[i]]$cell_origin <- paste(Idents(g.list[[i]]), g.list[[i]]$Sample.ID, sep = '_')
  g.list[[i]]$cell <- Idents(g.list[[i]])
  Idents(g.list[[i]]) <- 'cell_origin'
}

for(i in 1:length(g.list)){
  DE <- FindAllMarkers(g.list[[i]], min.pct = 0, logfc.threshold = 0, test.use = 'MAST')
  DE <- DE %>%
    filter(p_val_adj <= 0.05) %>%
    arrange(desc(avg_log2FC))
  write.csv(DE, here('Integrated', 'DE', 'within cluster', paste0('gran_', i, '.csv')))
}

rm(i,
   g.list,
   DE)

gc()



## fgsea
# geneSets <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
geneSets <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'MF') # GO:MF
# geneSets <- msigdbr(species = 'Mus musculus', category = 'C2') # C2
# geneSets <- msigdbr(species = 'Mus musculus', category = 'C7') # C7
# geneSets <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark
geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
for (i in seq_along(levels(neutro.markers$cluster))) {
  if(i != length(neutro.markers$cluster)) {
    glist <- neutro.markers %>% filter(neutro.markers$cluster %in% c(levels(neutro.markers$cluster)[i]))
    stats <- glist$avg_logFC
    names(stats) <- toupper(glist$gene)
    stats <- sort(stats, decreasing = T)
    eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
    eaRes <- filter(eaRes, padj <= 0.25 & NES > 0)
    eaRes <- arrange(eaRes, desc(NES))
    # fwrite(eaRes, file = here('Integrated', 'DE', 'fgsea', 'GO', paste0('GO_eaRes.gran', i, '.tsv')), sep="\t", sep2=c("", " ", ""))
    fwrite(eaRes, file = here('Integrated', 'DE', 'fgsea', 'GO', paste0('GO.MF_eaRes.gran', i, '.tsv')), sep="\t", sep2=c("", " ", ""))
    # fwrite(eaRes, file = here('Integrated', 'DE', 'fgsea', 'C2', paste0('C2_eaRes.gran', i, '.tsv')), sep="\t", sep2=c("", " ", ""))
    # fwrite(eaRes, file = here('Integrated', 'DE', 'fgsea', 'C7', paste0('C7_eaRes.gran', i, '.tsv')), sep="\t", sep2=c("", " ", ""))
    # fwrite(eaRes, file = here('Integrated', 'DE', 'fgsea', 'Hallmark', paste0('HM_eaRes.gran', i, '.tsv')), sep="\t", sep2=c("", " ", ""))
  }
}

rm(geneSets,
   m_list,
   i,
   stats,
   glist)

gc()



## Performing GSEA on average normalized counts
GO.set <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
HM.set <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark
CP.set <- msigdbr(species = 'Mus musculus', category = 'C2') %>% filter(gs_subcat != 'CGP') # Canonical Pathways

geneSet_list <- list(GO.set, HM.set, CP.set)
neutro.avg_data <- as.data.frame(neutro.avg@assays$RNA@data)

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in 1:ncol(neutro.avg_data)) {
    stats <- neutro.avg_data[, i]
    names(stats) <- toupper(rownames(neutro.avg_data))
    stats <- sort(stats, decreasing = T)
    eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
    eaRes <- filter(eaRes, padj <= 0.25)
    eaRes <- arrange(eaRes, desc(NES))
    fwrite(eaRes, file = here('Integrated', 'DE', 'fgsea', 'GSEA_Counts', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes.', ii, '.tsv')), sep="\t", sep2=c("", " ", ""))
  }
}

rm(geneSet_list,
   m_list,
   stats,
   eaRes,
   GO.set,
   HM.set,
   CP.set,
   i,
   ii,
   neutro.avg_data)



# Single-cell GSEA ------

GS <- getGeneSets(library = 'H', species = 'Mus musculus')
# GS <- getGeneSets(library = 'C5', species = 'Mus musculus')
ES <- enrichIt(neutro.integrated, gene.sets = GS, groups = 1000, cores = 2)
neutro.integrated <- AddMetaData(neutro.integrated, ES)

colors <- colorRampPalette(c("#FF4B20",
                             "#FFB433",
                             "#C6FDEC",
                             "#7AC5FF",
                             "#0348A6"))

multi_dittoPlot(neutro.integrated, vars = c('HALLMARK_APOPTOSIS',
                                            'HALLMARK_FATTY_ACID_METABOLISM',
                                            'HALLMARK_GLYCOLYSIS',
                                            'HALLMARK_HEME_METABOLISM',
                                            'HALLMARK_HYPOXIA',
                                            'HALLMARK_INFLAMMATORY_RESPONSE',
                                            'HALLMARK_MTORC1_SIGNALING',
                                            'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                                            'HALLMARK_PEROXISOME',
                                            'HALLMARK_PI3K_AKT_MTOR_SIGNALING',
                                            'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY'), 
                group.by = 'seurat_clusters', plots = c('jitter', 'vlnplot', 'boxplot'), 
                ylab = 'Enrichment Scores', 
                theme = theme_classic() +
                  theme(plot.title = element_text(size = 10)))



ES2 <- data.frame(neutro.integrated[[]], Idents(neutro.integrated))
colnames(ES2)[ncol(ES2)] <- 'cluster'
ridgeEnrichment(ES2, gene.set = 'HALLMARK_HYPOXIA', group = 'seurat_clusters', add.rug = TRUE) # Add 'facet = 'Sample.ID'' to split by day



PCA <- performPCA(enriched = ES2, groups = c('seurat_clusters', 'Sample.ID'))
pcaEnrichment(PCA, PCx = 'PC1', PCy = 'PC2', contours = TRUE)
pcaEnrichment(PCA, PCx = 'PC1', PCy = 'PC2', contours = FALSE, facet = 'seurat_clusters')

output <- getSignificance(ES2, group = 'seurat_clusters', fit = 'ANOVA') # Can use linear.model, T.test, or ANOVA



# Complex Heatmap ------

## Average single cell data & pull out normalized counts
## Counts are normalized by dividing the counts for a given feature by the total counts per cell, multiplying by a scale factor (default == 10,000), and then taking the natural log using log1p()
neutro.avg <- AverageExpression(neutro.integrated, return.seurat = T)
avg.norm_counts <- neutro.avg@assays$RNA@data
avg.norm_counts <- as.data.frame(avg.norm_counts)
avg.norm_counts <- rownames_to_column(avg.norm_counts, var = 'gene')

## Filter averaged normalized counts for genes of interest & convert back into a matrix
gene_list <- c('Gpi1',
               'Tpi1',
               'Gapdh',
               'Pgk1',
               'Pgam1',
               'Eno1',
               'Pkm',
               'G6pdx',
               'Pgls',
               'Pgd',
               'Tkt',
               'Taldo1')
heatmap.data <- avg.norm_counts %>%
  filter(gene %in% gene_list)
heatmap.data <- column_to_rownames(heatmap.data, var = 'gene')
heatmap.data <- as.matrix(heatmap.data)

## Create vectors for column/row annotations
mdsc_pmn <- c('MDSC',
              'MDSC',
              'MDSC',
              'MDSC',
              'MDSC',
              'MDSC',
              'MDSC',
              'PMN',
              'PMN',
              'PMN')
gene_anno <- c('PPP',
               'Glycolysis',
               'PPP',
               'Glycolysis',
               'Glycolysis',
               'Glycolysis',
               'Glycolysis',
               'PPP',
               'PPP',
               'PPP',
               'Glycolysis',
               'Glycolysis') # Check to make sure annotation order matches how the genes appear in the matrix

which(heatmap.data == max(heatmap.data), arr.ind = TRUE) # Use this to find the max value within the matrix so that you can set your upper bound

col_fun = circlize::colorRamp2(c(0, 4), c("white", "blue4")) # Check to make sure that the upper bound of this range isn't cutting off the expression of some genes in the matrix
ComplexHeatmap::Heatmap(heatmap.data,
                        name = 'Average\nNormalized\nCounts',
                        col = col_fun,
                        rect_gp = grid::gpar(col = 'black', lwd = 2),
                        column_split = mdsc_pmn,
                        cluster_columns = T,
                        cluster_column_slices = F,
                        column_gap = grid::unit(5, 'mm'),
                        row_split = gene_anno,
                        cluster_rows = T,
                        cluster_row_slices = F,
                        row_gap = grid::unit(5, 'mm'),
                        show_parent_dend_line = F,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid::grid.text(sprintf("%.1f", heatmap.data[i, j]), x, y, gp = grid::gpar(fontsize = 10))
                        })


# Notes -----

# If you'd like to do everything in three dimensions, use RumUMAP(object, dims = 1:n, n.components = 3L). Bringing this into
# Slingshot requires you to create objects that contain the cell embeddings (object@reductions$umap@cell.embeddings) & the
# clustering info (object@active.iden (or any grouping you desire)). You can then supply these to slingshot as you normally
# would. If you would like to visualize, you could do the following (requires rgl package):

# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# plot3d.SlingshotDataSet(sds)
# plot3d(rd, col = gg_color_hue(10)[cl], aspect = 'iso', add = T)


# I'm not sure if using three dimensions has any clear benefit over using two, but it may.



# Test area -----

## Monocle trajectories

cds <- as.cell_data_set(neutro.integrated)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

max.camp <- which.max(unlist(FetchData(integrated.sub, "Camp")))
max.camp <- colnames(integrated.sub)[max.camp]
cds <- order_cells(cds, root_cells = max.camp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)

# Set the assay back as 'integrated'
integrated.sub <- as.Seurat(cds, assay = "integrated")
FeaturePlot(integrated.sub, "monocle3_pseudotime")



