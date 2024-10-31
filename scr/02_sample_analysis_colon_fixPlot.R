# libraries ---------------------------------------------------------------
library(tidyverse)
library(MetaCycle)
library(ComplexHeatmap)
library(viridis)
library(cowplot)

# read data ---------------------------------------------------------------
# remove the genes that Cristina highlighted
df_model <- readRDS("../../out/object/df_model.rds") %>%
  separate(sample,into = c("gene","treat"),"_",remove = F) %>%
  filter(! gene %in% c("ROR alpha","Connexin 43","alpha Sintropin")) %>%
  select(-c(treat,gene))

# plotting the data -------------------------------------------------------

# test increase the resolution --------------------------------------------
# now use each model to produce a modelled range per gene
# x <- genes[1]
genes <- df_model %>%
  pull(sample) %>%
  unique()

df_predict <- lapply(genes,function(x){
  fit <- df_model %>%
    filter(sample == x) %>%
    pull(model) %>%
    .[[1]]
  
  params <- df_model %>%
    filter(sample == x) %>%
    pull(params) %>%
    .[[1]]
  
  newdata <- data.frame(tp = seq(from = 0,to=24,by = 1)) %>%
    mutate(phase = params$meta2d_phase,
           period = params$meta2d_period) %>%
    mutate(trend = tp-sum(tp)/nrow(.),
           amplitude = cos(2*pi*(tp-phase)/period)) %>%
    # add the prediction
    mutate(pred = predict(object = fit,newdata = .)) %>%
    mutate(sample = x)
  
  return(newdata)
  
}) %>%
  bind_rows()

# plot the heatmap usign the model value at higher resolution
mat3 <- df_predict %>%
  mutate(tp = as.factor(tp)) %>%
  select(sample,pred,tp) %>%
  pivot_wider(names_from = tp,values_from = pred) %>%
  column_to_rownames("sample")

# scale by sample
mat_scale3 <- t(scale(t(mat3)))
pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_sample_HigRes2.pdf",width = 4,height = 5)
# Heatmap(mat_scale,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale3,cluster_columns = F,cluster_rows = F,show_column_names = F)
dev.off()

# try to split the treatments side by side
mat_scale3 <- t(scale(t(mat3)))

list_mat_scale3 <- mat_scale3 %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(treat = str_extract(sample,"CTR|DSS"),
         sample2 = str_remove_all(sample,"_CTR|_DSS")) %>%
  split(.$treat) %>%
  lapply(function(x){
    mat <- x %>%
      rownames_to_column() %>%
      select(-c(rowname,treat,sample)) %>%
      column_to_rownames("sample2") %>%
      as.matrix()
    
    colnames(mat) <- str_remove(colnames(mat),"X")
    return(mat)
  })

# generate the list of heatmaps
heatmap_list <- pmap(list(list_mat_scale3,names(list_mat_scale3)), function(mat_x,name_mat) {
  Heatmap(mat_x, name = name_mat,cluster_columns = F,cluster_rows = F,column_title = name_mat,show_column_names = F)
})

pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_sample_HigRes2_split.pdf",width = 8,height = 5)
# Heatmap(mat_scale,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
draw(heatmap_list[[1]]+heatmap_list[[2]])
dev.off()

# reorder genes and plot the data
heatmap_list2 <- pmap(list(list_mat_scale3,names(list_mat_scale3)), function(mat_x,name_mat) {
  Heatmap(mat_x[c("NR1D1","BMAL1","CLOCK","PER1","PER2","CRY1","CRY2"),], name = name_mat,cluster_columns = F,cluster_rows = F,column_title = name_mat,show_column_names = F)
})

pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_sample_HigRes2_split2.pdf",width = 6,height = 3)
# Heatmap(mat_scale,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
draw(heatmap_list2[[1]]+heatmap_list2[[2]])
dev.off()

# try to plot the data as a scatter plot using the prediction values to smoothe the curve
df_model %>%
  unnest(data) %>%
  ggplot(aes(x=tp,y=value,col=treat)) +
  geom_point(shape=1) +
  geom_line(data = df_predict %>% separate(sample,into = c("geneName","treat"),"_"),aes(x=tp,y=pred,col=treat),linetype = "dashed") +
  facet_wrap(~geneName,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())+
  scale_color_manual(values = c("black","red"))
ggsave("../../out/plot/01_lanni_qpcr_gene_HigRes2.pdf",height = 9,width = 9)

# reorder the genes and change some graphical parameters
df_model %>%
  unnest(data) %>%
  mutate(geneName = factor(geneName,levels = c("PER1","PER2","CRY1","CRY2","NR1D1","BMAL1","CLOCK"))) %>%
  ggplot(aes(x=tp,y=value,col=treat)) +
  geom_point(shape=1) +
  geom_line(data = df_predict %>%
              separate(sample,into = c("geneName","treat"),"_") %>%
              mutate(geneName = factor(geneName,levels = c("PER1","PER2","CRY1","CRY2","NR1D1","BMAL1","CLOCK"))),
            aes(x=tp,y=pred,col=treat),linetype = "dashed") +
  facet_wrap(~geneName,scales = "free",as.table = F,nrow=2) +
  theme_cowplot() +
  theme(strip.background = element_blank())+
  scale_color_manual(values = c("black","red"))
ggsave("../../out/plot/01_lanni_qpcr_gene_HigRes3.pdf",height = 6,width = 13)
