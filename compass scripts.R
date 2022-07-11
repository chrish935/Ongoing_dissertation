## Made by Christopher M. Horn, MS
## Kielian Lab data
## Analyzing orthopedic model scRNA-seq compass results
## Created: 2021-08-05
## Updated: 2021-08-05


## Computing maturity score based on Seurat's AddModuleScore function

object = combined.express
features = list(gene_list)
pool = rownames(combined.express)
nbin = 24
ctrl = 100
k = FALSE
name = 'path_score'
seed = 1

cluster.length <- length(x = features)
assay.data <- object
data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
data.avg <- data.avg[order(data.avg)]
data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                                n = nbin,
                                labels = FALSE,
                                right = FALSE)

names(x = data.cut) <- names(x = data.avg)
ctrl.use <- vector(mode = "list", length = cluster.length)

for (i in 1:cluster.length) {
  features.use <- features[[i]]
  for (j in 1:length(x = features.use)) {
    ctrl.use[[i]] <- c(ctrl.use[[i]],
                       names(x = sample(x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
                                        size = ctrl,
                                        replace = FALSE)))
  }
}

ctrl.use <- lapply(X = ctrl.use, FUN = unique)
ctrl.scores <- matrix(data = numeric(length = 1L),
                      nrow = length(x = ctrl.use),
                      ncol = ncol(x = object))

for (i in 1:length(ctrl.use)) {
  features.use <- ctrl.use[[i]]
  ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use,])
}

features.scores <- matrix(data = numeric(length = 1L),
                          nrow = cluster.length,
                          ncol = ncol(x = object))

for (i in 1:cluster.length) {
  features.use <- features[[i]]
  data.use <- assay.data[features.use, , drop = FALSE]
  features.scores[i, ] <- Matrix::colMeans(x = data.use)
}

features.scores.use <- features.scores - ctrl.scores
rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
features.scores.use <- as.data.frame(x = t(x = features.scores.use))
rownames(x = features.scores.use) <- colnames(x = object)

rm(ctrl.scores, ctrl.use, data.use, features, features.scores, cluster.length, ctrl, data.avg, data.cut, features.use, i, j, k, name, nbin, pool, seed, assay.data, object)



## PCA of compass meta rxn consistencies

meta_rxn_consist <- read.csv('/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/compass meta rxn consistencies.csv', header = T)
meta_rxn_consist <- column_to_rownames(meta_rxn_consist, var = 'X')
meta_rxn_consist <- as.data.frame(t(meta_rxn_consist))
meta_rxn_consist.pca <- prcomp(meta_rxn_consist, scale. = T)
meta_rxn_consist.pca_coord <- meta_rxn_consist.pca$x[ , 1:3]
meta_rxn_consist.pca_load <- meta_rxn_consist.pca$rotation[ , 1:3]

write.csv(meta_rxn_consist.pca_coord, '/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/meta_rxn_pca_coord.csv')
write.csv(meta_rxn_consist.pca_load, '/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/meta_rxn_pca_loadings.csv')

met_activity <- read.csv('/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/metabolic_activity_calculations.csv', header = T) # PCA coords have already been added as columns
met_activity <- select(met_activity, -1)

# PCA plots

p1 <- ggplot(met_activity, aes(x = PC1, y = PC2)) +
  geom_point(shape = 21, color = 'black', size = 3, aes(fill = cell_type)) +
  theme_classic() +
  labs(x = 'PC1 (47.92%)', y = 'PC2 (13.16%)') +
  scale_fill_manual(values = c('tomato', 'slateblue1')) +
  guides(fill = guide_legend(title = 'Cell type'))

p2 <- ggplot(met_activity, aes(x = PC1, y = PC3)) +
  geom_point(shape = 21, color = 'black', size = 3, aes(fill = cell_type)) +
  theme_classic() +
  labs(x = 'PC1 (47.92%)', y = 'PC3 (4.23%)') +
  scale_fill_manual(values = c('tomato', 'slateblue1')) +
  guides(fill = guide_legend(title = 'Cell type'))

p3 <- ggplot(met_activity, aes(x = PC2, y = PC3)) +
  geom_point(shape = 21, color = 'black', size = 3, aes(fill = cell_type)) +
  theme_classic() +
  labs(x = 'PC2 (13.16%)', y = 'PC3 (4.23%)') +
  scale_fill_manual(values = c('tomato', 'slateblue1')) +
  guides(fill = guide_legend(title = 'Cell type'))

meta_rxn_consist.pca_plot <- ggpubr::ggarrange(p1, p2, p3, common.legend = T, legend = 'bottom', nrow = 1, ncol = 3)

rm(p1, p2, p3)

# PCs vs metabolic activity plots

p1 <- ggplot(met_activity, aes(x = PC1, y = metabolic_activity)) +
  geom_point(aes(shape = cell_type, color = Seurat_maturity_score), size = 3) +
  geom_smooth(method = 'lm', color = 'black') +
  theme_classic() +
  labs(x = 'PC1 (47.92%)', y = 'Metabolic activity') +
  scale_fill_manual(values = c('tomato', 'slateblue1')) +
  guides(shape = guide_legend(title = 'Cell type')) +
  viridis::scale_color_viridis(option = 'inferno', 'Maturity score') +
  annotate('text', x = 0, y = 0.12, label = paste('Pearson correlation:', round(cor(met_activity$metabolic_activity, met_activity$PC1), 3)))

p2 <- ggplot(met_activity, aes(x = PC2, y = metabolic_activity)) +
  geom_point(aes(shape = cell_type, color = Seurat_maturity_score), size = 3) +
  geom_smooth(method = 'lm', color = 'black') +
  theme_classic() +
  labs(x = 'PC2 (13.16%)', y = 'Metabolic activity') +
  scale_fill_manual(values = c('tomato', 'slateblue1')) +
  guides(shape = guide_legend(title = 'Cell type')) +
  viridis::scale_color_viridis(option = 'inferno', 'Maturity score') +
  annotate('text', x = 0, y = 0.12, label = paste('Pearson correlation:', round(cor(met_activity$metabolic_activity, met_activity$PC2), 3)))

p3 <- ggplot(met_activity, aes(x = PC3, y = metabolic_activity)) +
  geom_point(aes(shape = cell_type, color = Seurat_maturity_score), size = 3) +
  geom_smooth(method = 'lm', color = 'black') +
  theme_classic() +
  labs(x = 'PC3 (4.23%)', y = 'Metabolic activity') +
  scale_fill_manual(values = c('tomato', 'slateblue1')) +
  guides(shape = guide_legend(title = 'Cell type')) +
  viridis::scale_color_viridis(option = 'inferno', 'Maturity score') +
  annotate('text', x = -10, y = 0.12, label = paste('Pearson correlation:', round(cor(met_activity$metabolic_activity, met_activity$PC3), 3)))

meta_rxn_consist.pca2_plot <- ggpubr::ggarrange(p1, p2, p3, common.legend = T, legend = 'bottom', nrow = 1, ncol = 3)

rm(p1, p2, p3)

rm(met_activity, meta_rxn_consist, meta_rxn_consist.pca, meta_rxn_consist.pca_coord, meta_rxn_consist.pca_load, meta_rxn_consist.pca_plot, meta_rxn_consist.pca2_plot)



## Signed -log x spearman plots

meta_rxn_consist <- read.csv('/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/compass meta rxn consistencies.csv', header = T)
meta_rxn_consist <- column_to_rownames(meta_rxn_consist, var = 'X')
meta_rxn_consist <- as.data.frame(t(meta_rxn_consist))

met_activity <- read.csv('/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/metabolic_activity_calculations.csv', header = T) # PCA coords have already been added as columns
met_activity <- select(met_activity, -1)

patho_cor <- as.data.frame(cor(meta_rxn_consist[-1], met_activity$Seurat_patho_score, method = 'spearman'))
mature_cor <- as.data.frame(cor(meta_rxn_consist[-1], met_activity$Seurat_maturity_score, method = 'spearman'))

patho_cor <- rownames_to_column(patho_cor, var = 'rowname')
mature_cor <- rownames_to_column(mature_cor, var = 'rowname')

colnames(patho_cor)[2] <- 'patho_spearman'
colnames(mature_cor)[2] <- 'mature_spearman'

combine.cor <- merge(patho_cor, mature_cor, by = 'rowname')
rm(patho_cor, mature_cor)

compass.result <- read.csv('/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/compass analysis results.csv', header = T)
colnames(compass.result)[1] <- 'rowname'
compass.result <- select(compass.result, c(1, 4:8))
compass.result <- compass.result %>%
  mutate(neg_log_p = -log10(adjusted_pval))
compass.result <- compass.result %>%
  mutate(signed_neg_log_p = neg_log_p*sign(cohens_d))

compass.cor <- merge(compass.result, combine.cor, by = 'rowname')
rm(compass.result, met_activity, combine.cor)

write.csv(compass.cor, '/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/meta_rxn cor patho mature.csv')

compass.cor <- read.csv('/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/meta_rxn cor patho mature.csv', header = T) # Added column to highlight subsystems of interest

ggplot(compass.cor %>% arrange(label_point), aes(x = signed_neg_log_p, y = patho_spearman)) +
  geom_point(aes(fill = label_point), shape = 21, color = 'black', size = 3) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_classic() +
  labs(x = bquote(Signed~-log[10](BH-adjusted~Wilcoxon~rank~sum~p)), y = 'Spearman correlation with pathogenicity score') +
  scale_fill_manual(values = c('gray', 'red')) +
  Seurat::NoLegend()

ggplot(compass.cor %>% arrange(label_point), aes(x = signed_neg_log_p, y = mature_spearman)) +
  geom_point(aes(fill = label_point), shape = 21, color = 'black', size = 3) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_classic() +
  labs(x = bquote(Signed~-log[10](BH-adjusted~Wilcoxon~rank~sum~p)), y = 'Spearman correlation with maturity score') +
  scale_fill_manual(values = c('gray', 'red')) +
  Seurat::NoLegend()

rm(compass.cor, meta_rxn_consist)



## Generating marker gene x meta rxn heatmap

meta_rxn_consist <- read.csv('/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/compass meta rxn consistencies.csv', header = T)
meta_rxn_consist <- column_to_rownames(meta_rxn_consist, var = 'X')
meta_rxn_consist <- as.data.frame(t(meta_rxn_consist))
meta_rxn_consist <- rownames_to_column(meta_rxn_consist, var = 'rowname')

expression <- read.csv('/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/combined_express.csv', header = T)
colnames(expression)[1] <- 'rowname'
expression <- as.data.frame(expression)
expression <- expression %>%
  select(c('rowname', 'Il1b', 'Clec4e',	'Junb',	'Ctsd',	'Wfdc17',	'Il1f9',	'Pla2g7',	'Arg2',	'Cd84',	'Lcn2',	'Prdx5',	'Ngp',	'Camp',	'Ltf',	'Arhgdib',	'Anxa1',	'Plbd1',	'Tkt',	'Aldh2',	'Ly6c2',	'Adpgk',	'Cd177'))

combined.consist_express <- merge(expression, meta_rxn_consist, by = 'rowname')
combined.consist_express <- column_to_rownames(combined.consist_express, var = 'rowname')

meta_rxn_consist <- column_to_rownames(meta_rxn_consist, var = 'rowname')
meta_rxn_consist <- as.data.frame(t(meta_rxn_consist))

cor.table <- data.frame(reaction = rownames(meta_rxn_consist))
for (i in 1:22) {
  cor.result <- cor(combined.consist_express[23:ncol(combined.consist_express)], combined.consist_express[i], method = 'spearman')
  cor.table[ , ncol(cor.table) + 1] <- cor.result
  colnames(cor.table)[i+1] <- colnames(combined.consist_express)[i]
}

rm(cor.result, i, combined.consist_express, expression, meta_rxn_consist)

write.csv(cor.table, '/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/maturity genes cor compass heatmap.csv')

compass.result <- read.csv('/Volumes/External/Research/PhD/Kielian_Lab/Data/Compass_metabolism/Integrated_ortho/compass analysis results.csv', header = T)
compass.result <- compass.result %>%
  filter(adjusted_pval < 0.1)

cor.table <- column_to_rownames(cor.table, var = 'reaction')
cor.table <- cor.table %>%
  filter(rownames(.) %in% compass.result$X)

heatmap.data <- as.matrix(cor.table)

rm(cor.table, compass.result)

ha <- ComplexHeatmap::HeatmapAnnotation('Gene markers' = c(rep('MDSC', 9), rep('PMN', 13)), col = list('Gene markers' = c('MDSC' = 'tomato', 'PMN' = 'slateblue1')), gp = grid::gpar(col = "black"), show_annotation_name = F)
col_fun = circlize::colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red'))
ComplexHeatmap::Heatmap(heatmap.data,
                        name = 'Spearman correlation',
                        col = col_fun,
                        column_names_side = 'top',
                        cluster_columns = F,
                        cluster_column_slices = F,
                        cluster_rows = T,
                        cluster_row_slices = F,
                        show_parent_dend_line = F,
                        show_row_dend = F,
                        show_column_dend = F,
                        show_row_names = F,
                        column_names_rot = 45,
                        border = T,
                        top_annotation = ha)













