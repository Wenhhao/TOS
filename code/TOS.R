# 0.Set environment ----------
#rm(list = ls())
# gc()
pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)
library(ggrepel)

# 1.Read data and preprocess ---------
mk <- rio::import('subcell_location_MitochondriaOnly566.xlsx', sheet = 1)
df0 <- rio::import('TOS_IDH2_matrix.tsv')
df_label <- rio::import('TOS_IDH2_label.xlsx')
df <- df0[df0$Patient %in% df_label$Patient, ]
df <- rbind(df0[1, ], df)
df %<>% select(MS_file_name:Type, all_of(intersect(mk$Uniprot, colnames(df))))
colnames(df)[1:5]
dfprot <- data.frame(Protein = colnames(df)[-(1:4)],
                     Gene = unname(unlist(df[1, -(1:4)])))
df_info <- df %>% select(MS_file_name:Type) %>% slice(-1)
df <- df_info %>% left_join(df)
df %<>% mutate_at(vars(-(MS_file_name:Type)), as.numeric) # "NA" -> NA_double caused warnings
mat <- df %>% column_to_rownames('MS_file_name') %>% select(-(Patient:Type)) %>% t()


# 2.Missing ratio of proteins (by groups; and summarized) --------
missing_values <- apply(mat, 1, function(row) sum(is.na(row))/length(colnames(mat)))
result_all <- data.frame(Protein = paste(rownames(mat)), Missing_Values_all = missing_values)

mat_copy <-as.data.frame(mat)
colnames(mat_copy) <- df_info$TYPE
mat_FA <- mat_copy[ , colnames(mat_copy)=='FA']
rownames(mat_FA) <- rownames(mat_copy)
missing_values_FA <- apply(mat_FA, 1, function(row) sum(is.na(row))/length(colnames(mat_FA)))
result_FA <- data.frame(Protein = paste(rownames(mat_FA)), Missing_Values_FA = missing_values_FA)

mat_HCC <- mat_copy[ , colnames(mat_copy)=='HCC']
rownames(mat_HCC) <- rownames(mat_copy)
missing_values_HCC <- apply(mat_HCC, 1, function(row) sum(is.na(row))/length(colnames(mat_HCC)))
result_HCC <- data.frame(Protein = paste(rownames(mat_HCC)), Missing_Values_HCC = missing_values_HCC)

mat_HCA <- mat_copy[ , colnames(mat_copy)=='HCA']
rownames(mat_HCA) <- rownames(mat_copy)
missing_values_HCA <- apply(mat_HCA, 1, function(row) sum(is.na(row))/length(colnames(mat_HCA)))
result_HCA <- data.frame(Protein = paste(rownames(mat_HCA)), Missing_Values_HCA = missing_values_HCA)

mat_FTC <- mat_copy[ , colnames(mat_copy)=='FTC']
rownames(mat_FTC) <- rownames(mat_copy)
missing_values_FTC <- apply(mat_FTC, 1, function(row) sum(is.na(row))/length(colnames(mat_FTC)))
result_FTC <- data.frame(Protein = paste(rownames(mat_FTC)), missing_values_FTC = missing_values_FTC)

rlt_summary <- cbind(result_all, result_FA, result_FTC, result_HCA, result_HCC)
rlt_summary <- rlt_summary %>% select(-contains('Protein'))
df_missing <- rlt_summary



# 3.Fold-change (omit NA; and impute min0.8 NA) --------
nafill <- log2(0.8) + min(mat_copy, na.rm = T)

## 3.1 HCC / FTC --------
### 3.1.2 impute min0.8 NA------
mat_HCC_min <- mat_HCC
mat_HCC_min[is.na(mat_HCC_min)] <- nafill
mat_FTC_min <- mat_FTC
mat_FTC_min[is.na(mat_FTC_min)] <- nafill
FC_min_HCC_FTC <- c()
for (i in 1:nrow(mat_copy)) {
  cat(i, '...\r')
  x <- unlist(mat_HCC_min[i,])
  y <- unlist(mat_FTC_min[i,])
  FC_min_HCC_FTC[i] <- mean(t(2^mat_HCC_min[i,]), na.rm = T)/mean(t(2^mat_FTC_min[i,]), na.rm = T)
}

## 3.2 HCA / FA --------
### 3.2.2 impute min0.8 NA------
mat_HCA_min <- mat_HCA
mat_HCA_min[is.na(mat_HCA_min)] <- nafill
mat_FA_min <- mat_FA
mat_FA_min[is.na(mat_FA_min)] <- nafill
FC_min_HCA_FA <- c()
for (i in 1:nrow(mat_copy)) {
  cat(i, '...\r')
  x <- unlist(mat_HCA_min[i,])
  y <- unlist(mat_FA_min[i,])
  FC_min_HCA_FA[i] <- mean(t(2^mat_HCA_min[i,]), na.rm = T)/mean(t(2^mat_FA_min[i,]), na.rm = T)
}

FC_rlt_summary <- cbind(FC_min_HCA_FA, FC_min_HCC_FTC)
rownames(FC_rlt_summary) <- rownames(mat_copy)
df_fc <- FC_rlt_summary %>% as.data.frame()


# 4.Differential expression analysis (omit NA; and impute min0.8 NA) --------
## 4.1 Welch's t-test (HCC vs. FTC; and HCA vs. FA) --------
### HCC vs. FTC ---------
T_min_HCC_FTC <- sapply(1:nrow(mat_HCC_min), function(i){
  x <- unlist(mat_HCC_min[i, ])
  y <- unlist(mat_FTC_min[i, ])
  tryCatch(t.test(x, y, var.equal = F, paired = F)$p.value, error = function(e) {NA})
})

### HCA vs. FA ---------
T_min_HCA_FA <- sapply(1:nrow(mat_HCA_min), function(i){
  x <- unlist(mat_HCA_min[i, ])
  y <- unlist(mat_FA_min[i, ])
  tryCatch(t.test(x, y, var.equal = F, paired = F)$p.value, error = function(e) {NA})
})

## 4.2 Wilconxon's sum-rank test (HCC vs. FTC; and HCA vs. FA) --------
### HCC vs. FTC ---------
identical(rownames(mat_HCC_min), rownames(mat_FTC_min)) # TRUE
U_min_HCC_FTC <- sapply(1:nrow(mat_HCC_min), function(i){
  x <- unlist(mat_HCC_min[i, ])
  y <- unlist(mat_FTC_min[i, ])
  tryCatch(wilcox.test(x, y, exact = F)$p.value, error = function(e) {NA})
})

### HCA vs. FA ---------
identical(rownames(mat_HCA_min), rownames(mat_FA_min)) # TRUE
U_min_HCA_FA <- sapply(1:nrow(mat_HCA_min), function(i){
  x <- unlist(mat_HCA_min[i, ])
  y <- unlist(mat_FA_min[i, ])
  tryCatch(wilcox.test(x, y, exact = F)$p.value, error = function(e) {NA})
})

## 4.3 limma moderate t-test (HCC vs. FTC; and HCA vs. FA; NA impute 0) --------
library(limma)
### HCC vs. FTC ---------
Mat <- cbind(mat_HCC_min, mat_FTC_min)
Type <- factor(c(rep('HCC', ncol(mat_HCC_min)), rep('FTC', ncol(mat_FTC_min))), levels = c('HCC', 'FTC'))
design <- model.matrix(~0 + Type)
contrast <- makeContrasts(Type = 'TypeHCC-TypeFTC', levels = design)

fit1 <- lmFit(Mat, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)
res_min_HCC_FTC <- topTable(fit3, coef = 'Type', number = Inf)
colnames(res_min_HCC_FTC) %<>% str_c('limma_min_', ., '_HCC_FTC')

### HCA vs. FA ---------
Mat <- cbind(mat_HCA_min, mat_FA_min)
Type <- factor(c(rep('HCA', ncol(mat_HCA_min)), rep('FA', ncol(mat_FA_min))), levels = c('HCA', 'FA'))
design <- model.matrix(~0 + Type)
contrast <- makeContrasts(Type = 'TypeHCA-TypeFA', levels = design)

fit1 <- lmFit(Mat, design)
fit2 <- contrasts.fit(fit1, contrasts = contrast)
fit3 <- eBayes(fit2)
res_min_HCA_FA <- topTable(fit3, coef = 'Type', number = Inf)
colnames(res_min_HCA_FA) %<>% str_c('limma_min_', ., '_HCA_FA')


P_rlt_summary <- cbind(T_min_HCC_FTC, T_min_HCA_FA, U_min_HCC_FTC, U_min_HCA_FA, res_min_HCC_FTC[rownames(mat_copy), ], res_min_HCA_FA[rownames(mat_copy), ])
rownames(P_rlt_summary) <- rownames(mat_copy)
df_p <- P_rlt_summary %>% as.data.frame()
df_rlt <- cbind(df_missing, df_fc, df_p) %>% rownames_to_column('Protein')
df_rlt %<>% mutate(
  adj.T_min_HCC_FTC = p.adjust(T_min_HCC_FTC, method = 'BH'),
  adj.T_min_HCA_FA = p.adjust(T_min_HCA_FA, method = 'BH'),
  adj.U_min_HCC_FTC = p.adjust(U_min_HCC_FTC, method = 'BH'),
  adj.U_min_HCA_FA = p.adjust(U_min_HCA_FA, method = 'BH')
)
df_rlt %<>% inner_join(dfprot, .)
df_rlt %>% rio::export('DEA_386_20240628.xlsx')




# 5.Select TOS (thyroid oncocytic score) markers --------
## 5.1 considering multiple conditions --------
# cutoff of < 50% missing ratio
df_rlt <- rio::import('DEA_386_20240628.xlsx')
X <- df_rlt %>% filter(if_all(matches('^Missing_Values'), function(x) x <= 0.4))

# select proteins with all FC > 4 and all adjusted P < 0.05
Y <- X %>% filter(if_all(matches('^FC'), function(x) {x > 4}) | if_all(matches('^FC'), function(x) {x < 1/4})) %>% 
  filter(if_all(matches('adj\\.'), function(x) x < 0.01))

Y %<>% add_column(
  adj.P.max = Y %>% select(matches('adj\\.')) %>% apply(1, max),
  FC.min = Y %>% select(matches('^FC')) %>% apply(1, min),
  Missing.max = Y %>% select(matches('^Missing')) %>% apply(1, max),
  .after = 'Gene'
) %>% 
  arrange(desc(FC.min))
Y %>% rio::export('TOS_124markers_v1_20240628.xlsx')




### volcano plot ----
df_volcano <- df_rlt %>% add_column(
  adj.P.max = df_rlt %>% select(matches('adj\\.')) %>% apply(1, max, na.rm = T),
  FC.min = df_rlt %>% select(matches('^FC')) %>% apply(1, min, na.rm = T),
  Missing.max = df_rlt %>% select(matches('^Missing')) %>% apply(1, max, na.rm = T),
  .after = 'Gene'
) 
df_volcano$IsSignificant <- df_volcano$Protein %in% Y$Protein

plot_vol <- ggplot()+
  aes(x = log2(FC.min), y = -log10(adj.P.max)) +
  geom_point(data = df_volcano %>% filter(!IsSignificant), aes(size = abs(log2(FC.min))), alpha = 0.3, color = '#AAAAAA') +
  geom_point(data = df_volcano %>% filter(IsSignificant, FC.min > 1), aes(size = abs(log2(FC.min)), color = -log10(adj.P.max) * sign(log2(FC.min))), alpha = 0.95) +
  scale_color_gradient(low = '#EEEEEE', high = '#FE5936', breaks = c(0, 15)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = log2(4)) +
  geom_vline(xintercept = -log2(4)) +
  labs(x = 'Log2(FC.min)', y = '-Log10(adj.P.max)', subtitle = 'HCC/FTC & HCA/FA') +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  geom_label_repel(
    data = df_volcano %>% filter(IsSignificant, str_detect(Gene, '^IDH')),
    aes(label = Gene), seed = 1,
    nudge_x = 4,
    point.padding = 0.5,# max.overlaps = 100,
    size = 3,
    segment.size = 0.5,
    direction = "y"
  )
ggsave('TOS_124markers_v1_20240628.pdf', plot_vol, width = 6, height = 4.5)




### boxplot ----
idh <- Y %>% filter(str_detect(Gene, '^IDH')) %>% add_column(label = str_glue("{.$Gene} ({.$Protein})"), .before = 1)
df_box <- df[-1, ] %>% select(TYPE, all_of(idh$Protein)) %>% 
  mutate(TYPE = factor(TYPE, levels = c('FA', 'FTC', 'HCA', 'HCC')))


df_box3 <- df_box %>%
  pivot_longer(cols = -TYPE, names_to = 'Protein', values_to = 'Log2Intensity') %>% 
  inner_join(idh)

plot_box3 <- ggplot(df_box3) +
  facet_wrap(facets = ~label, nrow = 1, ncol = 5)+
  aes(x = TYPE, y = Log2Intensity) +
  geom_jitter(size = 0.1, width = 0.3, color = '#AAAAAA') +
  geom_boxplot(outlier.shape = NA, fill = NA, color = '#000000') +
  labs(x = '', y = 'Log2Intensity') +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  ggsignif::geom_signif(
    comparisons = list(c('HCA', 'FA'), c('HCC', 'FTC')),
    step_increase = 0.1,
    test = wilcox.test)

ggsave('TOS_124markers_v1_3box_20240628_addWilcox.pdf', plot_box3, width = 12, height = 4)



# a simple pathway enrichment
library(clusterProfiler)
library(org.Hs.eg.db)

GO <- enrichGO(Y$Gene, OrgDb = org.Hs.eg.db,
               ont = 'ALL', keyType = 'SYMBOL', readable = T,
               minGSSize = 10, maxGSSize = 500,
               pvalueCutoff = 0.05, pAdjustMethod = 'BH')
GO_simplify <- clusterProfiler::simplify(GO, cutoff = 0.7, by = 'p.adjust', select_fun = min)
dfgo <- GO[]
dfgo_simplify <- GO_simplify[]

pdf('TOS_GO_enrichment_top5_bubble_124marker_20240628.pdf', width = 6.5, height = 9)
enrichplot::dotplot(GO, showCategory = 5, split = 'ONTOLOGY', title = 'Gene Ontology') + facet_grid(ONTOLOGY~., scale='free')
enrichplot::dotplot(GO_simplify, showCategory = 5, split = 'ONTOLOGY', title = 'Gene Ontology (simplified)') + facet_grid(ONTOLOGY~., scale='free')
similar_matrix <- dfgo_simplify %>%
  filter(ONTOLOGY == 'BP') %>% pull(ID) %>% 
  simplifyEnrichment::GO_similarity(ont = 'BP', db = 'org.Hs.eg.db')
GO_simplify_sim <- simplifyEnrichment::simplifyGO(similar_matrix)
graphics.off()



pacman::p_unload(pacman::p_loaded(), character.only = T)
library(magrittr)
library(tidyverse)


## 5.2 refer to TDS --------
TOS <- function(X, log2trans = F, naimp = NULL, center = T, zscore = T){
  # each column is a gene
  if(log2trans){
    X %<>% log2()
  }
  if(!is.null(naimp)){
    X[is.na(X)] <- naimp
  }
  if(center){
    X <- scale(t(X), center = apply(X, 1, median), scale = F) %>% t() # centered at the median
  }
  TOS <- rowMeans(X)
  # TOS <- apply(X, 1, median)
  if(zscore){
    TOS <- scale(TOS)
  }
  data.frame(TOS = TOS, row.names = rownames(X)) %>% return()
}

Y <- rio::import('TOS_124markers_v1_20240628.xlsx')
prot_TOS <- Y$Protein



X <- t(mat[prot_TOS, ])
TOS <- TOS(X, log2trans = F, naimp = min(X, na.rm = T)+log2(0.8), center = F, zscore = T)
TOS %<>% rownames_to_column('MS_file_name') %>% left_join(df_info)

p <- ggplot(TOS) +
  aes(x = TYPE, y = TOS) +
  geom_violin(fill = '#11244655') +
  geom_boxplot(color = '#112446', width = 0.1) +
  theme_bw()
ggsave('TOS_boxplot_124marker_20240628.pdf', p, width = 5, height = 5)


### Heatmap -----------------------------------------------------------------
df2 <- df[, Y$Protein]
rownames(df2) <- df$MS_file_name
df2 <- 2^df2
df2 <- add_column(df2, Patient = df$Patient, TYPE = df$TYPE, .before = 1)
df2_average <- df2 %>% group_by(Patient, TYPE) %>% 
  summarise_all(mean, na.rm = T) %>% 
  ungroup()

TOS_average <- TOS[, c(1:4)] %>% column_to_rownames('MS_file_name') %>% 
  group_by(Patient, TYPE) %>% 
  summarise_all(mean, na.rm = T) %>% 
  ungroup()

library(pheatmap)
library(RColorBrewer)


ann_col <- df_label[, c(1,2,4,5)] %>% 
  inner_join(TOS_average) %>% 
  select(TOS, TYPE, Gender, Age, everything()) %>% 
  column_to_rownames('Patient')

ann_col$TYPE <- factor(ann_col$TYPE, levels = c('FA', 'FTC', 'HCA', 'HCC'))
ann_col %<>% arrange(TYPE, TOS)
ann_clr <- list(
  TYPE = c(HCC = "#EB9D15",HCA = "#F0C116",FTC = "#0E71EB",FA = "#5A9DE0"),
  TOS = brewer.pal(11, 'PRGn'),
  Gender = c('F' = '#F19695', 'M' = '#A3C9DD')
)
mycolor <- colorRampPalette(rev(brewer.pal(9, 'Spectral')))(50)
heat_X <- df2_average %>% column_to_rownames('Patient') %>% select(-TYPE) %>% 
  log2() %>% t() %>% .[, rownames(ann_col)]
clust_col <- hclust(dist(t(heat_X), method = 'euclidean'), method = 'ward.D2') %>%
  as.dendrogram() %>%
  reorder(wts = 1:ncol(heat_X)) %>% # reorder to keep original order as much as possible
  as.hclust()
plot_heat <- pheatmap(heat_X,
                      scale = "row",
                      color = mycolor,
                      fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                      border_color = F,
                      annotation_col = ann_col,
                      annotation_colors = ann_clr,
                      cluster_rows = T, cluster_cols = clust_col,
                      clustering_distance_rows = 'euclidean',
                      clustering_distance_cols = 'euclidean',
                      clustering_method = 'ward.D2',
                      show_rownames = T, show_colnames = T
)
ggsave(filename = 'TOS_124markers_heatmap.pdf', plot_heat, width = 8, height = 12)
ggsave(filename = 'TOS_124markers_heatmap_wide.pdf', plot_heat, width = 25, height = 12)




plot_heatv2 <- pheatmap(heat_X,
                        scale = "row",
                        color = mycolor,
                        fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                        border_color = F,
                        annotation_col = ann_col,
                        annotation_colors = ann_clr,
                        cluster_rows = T, cluster_cols = F,
                        clustering_distance_rows = 'euclidean',
                        clustering_distance_cols = 'euclidean',
                        clustering_method = 'ward.D2',
                        show_rownames = T, show_colnames = T
)
ggsave(filename = 'TOS_124markers_heatmapv2.pdf', plot_heatv2, width = 8, height = 12)
ggsave(filename = 'TOS_124markers_heatmapv2_wide.pdf', plot_heatv2, width = 25, height = 12)


### metascape ---------------------------------------------------------------
df_meta <- rio::import('metascape_result.xlsx', sheet = 2)
df_meta <- df_meta[str_detect(df_meta$GroupID, 'Summary'), ]
df_meta_split <- df_meta %>% separate_rows(Symbols)
# str_subset(df_meta_split$Symbols, '\\W')
df_meta1 <- rio::import('metascape_result.xlsx', sheet = 1) 
df_meta1 %<>% select(MyList, `Gene Symbol`) %>%
  inner_join(df_meta_split, by = c('Gene Symbol' = 'Symbols')) %>% 
  distinct(MyList, .keep_all = T) %>% 
  arrange(GroupID)
df_meta1 %<>% mutate(No = as.numeric(str_extract(GroupID, '[0-9]+')), .before = 1) %>% 
  arrange(No)

ann_row <- df_meta1[, c('MyList', 'Description')] %>% column_to_rownames('MyList')
ann_row$Description <- factor(ann_row$Description, levels = unique(ann_row$Description))
ann_clr <- list(
  TYPE = c(HCC = "#EB9D15",HCA = "#F0C116",FTC = "#0E71EB",FA = "#5A9DE0"),
  Age = brewer.pal(9, 'Greens'),
  TOS = brewer.pal(11, 'PRGn'),
  Gender = c('F' = '#F19695', 'M' = '#A3C9DD'),
  Description = c('#5FBFF9', '#F9A6A6', '#C9F9A6', '#A6CDF9', '#FFA500', '#FFD700', '#FF6347'#, '#FFE4B5', '#228B22'
                  ) %>% setNames(levels(ann_row$Description))
)
heat_X <- df2_average %>% column_to_rownames('Patient') %>% select(-TYPE) %>% 
  log2() %>% t() %>% .[rownames(ann_row), rownames(ann_col)]
clust_col <- hclust(dist(t(heat_X), method = 'euclidean'), method = 'ward.D2') %>%
  as.dendrogram() %>%
  reorder(wts = 1:ncol(heat_X)) %>% # reorder to keep original order as much as possible
  as.hclust()
plot_heat <- pheatmap(heat_X,
                      scale = "row",
                      color = mycolor,
                      fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                      border_color = F,
                      annotation_col = ann_col,
                      annotation_row = ann_row,
                      annotation_colors = ann_clr,
                      cluster_rows = T, cluster_cols = clust_col,
                      clustering_distance_rows = 'euclidean',
                      clustering_distance_cols = 'euclidean',
                      clustering_method = 'ward.D2',
                      show_rownames = T, show_colnames = T
)
ggsave(filename = 'TOS_124markers_pathway_heatmap.pdf', plot_heat, width = 8, height = 12)

plot_heatv2 <- pheatmap(heat_X,
                        scale = "row",
                        color = mycolor,
                        fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                        border_color = F,
                        annotation_col = ann_col,
                        annotation_row = ann_row,
                        annotation_colors = ann_clr,
                        cluster_rows = F, cluster_cols = F,
                        clustering_distance_rows = 'euclidean',
                        clustering_distance_cols = 'euclidean',
                        clustering_method = 'ward.D2',
                        show_rownames = T, show_colnames = T
)
ggsave(filename = 'TOS_124markers_pathway_heatmapv2.pdf', plot_heatv2, width = 8, height = 6)

heat2_col_order <- data.frame(Patient = clust_col$labels[clust_col$order],
                              heat1_order = 1:length(clust_col$labels)) %>% 
  inner_join(ann_col %>% rownames_to_column('Patient')) %>% 
  arrange(TYPE, heat1_order)
plot_heatv3 <- pheatmap(heat_X[, heat2_col_order$Patient],
                        scale = "row",
                        color = mycolor,
                        fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                        border_color = F,
                        annotation_col = ann_col,
                        annotation_row = ann_row,
                        annotation_colors = ann_clr,
                        cluster_rows = F, cluster_cols = F,
                        clustering_distance_rows = 'euclidean',
                        clustering_distance_cols = 'euclidean',
                        clustering_method = 'ward.D2',
                        show_rownames = T, show_colnames = T
)
ggsave(filename = 'TOS_124markers_pathway_heatmapv3.pdf', plot_heatv3, width = 8, height = 6)


heat2_col_order <- data.frame(Patient = clust_col$labels[clust_col$order],
                              heat1_order = 1:length(clust_col$labels)) %>% 
  inner_join(ann_col %>% rownames_to_column('Patient')) %>% 
  arrange(TYPE, TOS)
plot_heatv4 <- pheatmap(heat_X[, heat2_col_order$Patient],
                        scale = "row",
                        color = mycolor,
                        fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
                        border_color = F,
                        annotation_col = ann_col,
                        annotation_row = ann_row,
                        annotation_colors = ann_clr,
                        cluster_rows = F, cluster_cols = F,
                        clustering_distance_rows = 'euclidean',
                        clustering_distance_cols = 'euclidean',
                        clustering_method = 'ward.D2',
                        show_rownames = T, show_colnames = T
)
ggsave(filename = 'TOS_124markers_pathway_heatmapv4.pdf', plot_heatv4, width = 8, height = 6)

lgd1 <- pheatmap(TOS %>% column_to_rownames('MS_file_name') %>% select(TOS),
         scale = "none",
         color = colorRampPalette(brewer.pal(11, 'PRGn'))(100),
         fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
         border_color = F,
         cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = T,
) %>% ggplotify::as.ggplot()
lgd2 <- pheatmap(ann_col %>% select(Age),
         scale = "none",
         color = colorRampPalette(brewer.pal(9, 'Greens'))(100),
         fontsize_col = 8, fontsize_row = 8, fontsize_number = 8,
         border_color = F,
         cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = T,
) %>% ggplotify::as.ggplot()
ggsave('TOS_124markers_pathway_heatmapv4_legend.pdf', lgd1 + lgd2, width = 2, height = 6)

