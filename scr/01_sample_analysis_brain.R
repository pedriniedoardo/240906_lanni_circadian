# libraries ---------------------------------------------------------------
library(tidyverse)
library(MetaCycle)
library(ComplexHeatmap)
library(viridis)

# read data ---------------------------------------------------------------
df <- read_csv("../../data/241028_lanni_qpcr_brain.csv")

# wrangling ---------------------------------------------------------------
# extract the info from the sample column
df_full <- df %>%
  mutate(treat = str_extract(sample_full,pattern = "CTR|DSS"),
         replicate = str_extract(sample_full,pattern = "^\\d+ \\w"),
         replicate2 = str_extract(replicate,"\\d+"))

# shape the data to be accepted by MetaCycle
df_full2 <- df_full %>%
  mutate(sample = paste0(gene,"_",treat,"_",dataset)) %>%
  mutate(tp_rep = paste0("T_",time_point,"_",replicate2)) %>%
  select(sample,tp_rep,exp) %>%
  pivot_wider(names_from = tp_rep,values_from = exp)

# save it as a table
# write.csv(df_full2, file="../../data/241028_lanni_qpcr_brain_clean.csv", row.names=FALSE)

# EDA ---------------------------------------------------------------------
# explore the dataset
df_long <- df_full2 %>%
  pivot_longer(names_to = "tp_cat",values_to = "value",-sample) %>%
  separate(sample,into = c("geneName","treat","dataset"),sep="_",remove = F) %>%
  separate(tp_cat,into = c("X","tp","rep_id"),sep="_",remove = F) %>%
  dplyr::select(-X) %>%
  mutate(tp = as.numeric(tp))

# plot the data split by gene
df_long %>%
  ggplot(aes(x=tp,y=value,col=treat))+geom_point(shape=1) +
  facet_grid(dataset~geneName,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_color_manual(values = c("black","red"))

# estimate the period and phase with MetaCycle ----------------------------
# define the parameters
cyc <- meta2d(infile = "../../data/241028_lanni_qpcr_brain_clean.csv",
              filestyle = "csv",
              timepoints = rep(seq(0, 24, by=4), each=3),
              outputFile = FALSE,
              outRawData = TRUE)

# use the prediction to estimate the trend and amplitude paramters --------
# pull all the genes
genes <- cyc$meta$CycID

# define model, paramters and prediction for each gene
# x <- genes[1]
df_model <- lapply(genes, function(x){
  # pull the gene specific parameters
  params <- cyc$meta %>%
    filter(CycID %in% c(x)) %>%
    select(CycID,contains("meta"))
  
  # pull the data for the model
  df_test <- df_long %>%
    filter(sample %in% x) %>%
    mutate(phase = params$meta2d_phase,
           period = params$meta2d_period) %>%
    mutate(trend = tp-sum(tp)/nrow(.),
           amplitude = cos(2*pi*(tp-phase)/period))
  
  # model.matrix(~ trend+amplitude, data=df_test)
  # build the model
  fit <- lm(value ~ trend+amplitude, data=df_test)
  # summary(fit)
  
  # make the prediction using fit
  df_final <- df_test %>%
    mutate(pred = predict(object = fit,newdata = df_test))
  
  # add the data, the model and the gene name in a single object
  df_out <- as.tibble(df_final) %>%
    group_by(sample) %>%
    nest() %>%
    mutate(params = list(params)) %>%
    mutate(model = list(fit))
  
  
  return(df_out)
}) %>%
  bind_rows()

# save the object
# saveRDS(df_model, file="../../out/object/df_model_brain.rds")

# plotting the data -------------------------------------------------------
# plot all the genes with the relative model
df_model %>%
  unnest(data) %>%
  ggplot(aes(x=tp,y=value)) +
  geom_point(shape=1) +
  geom_line(aes(x=tp,y=pred),col="red",linetype = "dashed") +
  facet_grid(~sample,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())
# ggsave("../../out/plot/01_lanni_qpcr_sample.pdf",height = 12,width = 15)

# plot all the genes colored by treat
df_model %>%
  unnest(data) %>%
  ggplot(aes(x=tp,y=value,col=treat)) +
  geom_point(shape=1) +
  geom_line(aes(x=tp,y=pred,col=treat),linetype = "dashed") +
  facet_grid(dataset~geneName,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())+
  scale_color_manual(values = c("black","red"))
ggsave("../../out/plot/01_lanni_qpcr_gene_brain.pdf",height = 9,width = 22)

# add some jittering to the points
df_model %>%
  unnest(data) %>%
  ggplot(aes(x=tp,y=value,col=treat)) +
  geom_point(shape=1,position = position_jitter(width = 1)) +
  geom_line(aes(x=tp,y=pred,col=treat),linetype = "dashed") +
  facet_grid(dataset~geneName,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())+
  scale_color_manual(values = c("black","red"))
ggsave("../../out/plot/01_lanni_qpcr_gene_jitter_brain.pdf",height = 9,width = 22)

# plot the raw values as a heatmap
mat <- df_model %>%
  unnest(data) %>%
  mutate(tp = as.factor(tp)) %>%
  # select(sample,value,tp) %>%
  group_by(sample,tp) %>%
  summarise(avg_value = mean(value,na.rm = T)) %>%
  pivot_wider(names_from = tp,values_from = avg_value) %>%
  column_to_rownames("sample")

# scale by sample
mat_scale <- t(scale(t(mat)))

list_mat_scale <- mat_scale %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(dataset = str_extract(sample,"corteccia|ipotalamo|ippocampo"),
         sample2 = str_remove_all(sample,"_corteccia|_ipotalamo|_ippocampo")) %>%
  split(.$dataset) %>%
  lapply(function(x){
    mat <- x %>%
      rownames_to_column() %>%
      select(-c(rowname,dataset,sample)) %>%
      column_to_rownames("sample2") %>%
      as.matrix()
    
    colnames(mat) <- str_remove(colnames(mat),"X")
    return(mat)
  })

# generate the list of heatmaps
heatmap_list <- pmap(list(list_mat_scale,names(list_mat_scale)), function(mat_x,name_mat) {
  Heatmap(mat_x, name = name_mat,cluster_columns = F,cluster_rows = F,column_title = name_mat)
})

draw(heatmap_list[[1]]+heatmap_list[[2]]+heatmap_list[[3]])

list_mat_scale2 <- mat_scale %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(dataset = str_extract(sample,"corteccia|ipotalamo|ippocampo"),
         sample2 = str_remove_all(sample,"_corteccia|_ipotalamo|_ippocampo"),
         treat = str_extract(sample,"CTR|DSS"),
         sample3 = str_remove_all(sample,"CTR_|DSS_")) %>%
  split(.$treat) %>%
  lapply(function(x){
    mat <- x %>%
      rownames_to_column() %>%
      select(-c(rowname,dataset,sample2,treat,sample)) %>%
      column_to_rownames("sample3") %>%
      as.matrix()
    
    colnames(mat) <- str_remove(colnames(mat),"X")
    return(mat)
  })

# generate the list of heatmaps
heatmap_list2 <- pmap(list(list_mat_scale2,names(list_mat_scale2)), function(mat_x,name_mat) {
  Heatmap(mat_x, name = name_mat,cluster_columns = F,cluster_rows = F,column_title = name_mat)
})

draw(heatmap_list2[[1]]+heatmap_list2[[2]])

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

list_mat_scale3 <- mat_scale3 %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(dataset = str_extract(sample,"corteccia|ipotalamo|ippocampo"),
         sample2 = str_remove_all(sample,"_corteccia|_ipotalamo|_ippocampo")) %>%
  split(.$dataset) %>%
  lapply(function(x){
    mat <- x %>%
      rownames_to_column() %>%
      select(-c(rowname,dataset,sample)) %>%
      column_to_rownames("sample2") %>%
      as.matrix()
    
    colnames(mat) <- str_remove(colnames(mat),"X")
    return(mat)
  })

# generate the list of heatmaps
heatmap_list3 <- pmap(list(list_mat_scale3,names(list_mat_scale3)), function(mat_x,name_mat) {
  Heatmap(mat_x, name = name_mat,cluster_columns = F,cluster_rows = F,column_title = name_mat,show_column_names = F)
})

pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_sample_HigRes2_brain.pdf",width = 10,height = 5)
# Heatmap(mat_scale,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
draw(heatmap_list3[[1]]+heatmap_list3[[2]]+heatmap_list3[[3]])
dev.off()

# try to split the treatments side by side
list_mat_scale4 <- mat_scale3 %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(dataset = str_extract(sample,"corteccia|ipotalamo|ippocampo"),
         sample2 = str_remove_all(sample,"_corteccia|_ipotalamo|_ippocampo"),
         treat = str_extract(sample,"CTR|DSS"),
         sample3 = str_remove_all(sample,"CTR_|DSS_")) %>%
  split(.$treat) %>%
  lapply(function(x){
    mat <- x %>%
      rownames_to_column() %>%
      select(-c(rowname,dataset,sample2,treat,sample)) %>%
      column_to_rownames("sample3") %>%
      as.matrix()
    
    colnames(mat) <- str_remove(colnames(mat),"X")
    return(mat)
  })

# generate the list of heatmaps
heatmap_list4 <- pmap(list(list_mat_scale4,names(list_mat_scale4)), function(mat_x,name_mat) {
  Heatmap(mat_x, name = name_mat,cluster_columns = F,cluster_rows = F,column_title = name_mat,show_column_names = F)
})

pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_sample_HigRes2_brain_split.pdf",width = 8,height = 5)
# Heatmap(mat_scale,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
draw(heatmap_list4[[1]]+heatmap_list4[[2]])
dev.off()

# try to plot the data as a scatter plot using the prediction values to smoothe the curve
df_model %>%
  unnest(data) %>%
  ggplot(aes(x=tp,y=value,col=treat)) +
  geom_point(shape=1) +
  geom_line(data = df_predict %>% separate(sample,into = c("geneName","treat","dataset"),"_"),aes(x=tp,y=pred,col=treat),linetype = "dashed") +
  facet_grid(dataset~geneName,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())+
  scale_color_manual(values = c("black","red"))
ggsave("../../out/plot/01_lanni_qpcr_gene_HigRes2_brain.pdf",height = 9,width = 12)

df_model %>%
  unnest(data) %>%
  ggplot(aes(x=tp,y=value,col=treat)) +
  geom_point(shape=1) +
  geom_line(data = df_predict %>% separate(sample,into = c("geneName","treat","dataset"),"_"),aes(x=tp,y=pred,col=treat),linetype = "dashed") +
  facet_wrap(dataset~geneName,scales = "free",ncol = 7) +
  theme_bw() +
  theme(strip.background = element_blank())+
  scale_color_manual(values = c("black","red"))
ggsave("../../out/plot/01_lanni_qpcr_gene_HigRes3_brain.pdf",height = 9,width = 12)

df_model %>%
  unnest(data) %>%
  ggplot(aes(x=tp,y=value,col=treat)) +
  geom_point(shape=1) +
  geom_line(data = df_predict %>% separate(sample,into = c("geneName","treat","dataset"),"_"),aes(x=tp,y=pred,col=treat),linetype = "dashed") +
  facet_wrap(dataset~geneName,scales = "free",ncol = 7) +
  theme_cowplot() +
  theme(strip.background = element_blank())+
  scale_color_manual(values = c("black","red"))
ggsave("../../out/plot/01_lanni_qpcr_gene_HigRes3_brain2.pdf",height = 9,width = 15)

df_model %>%
  unnest(data) %>%
  mutate(geneName = factor(geneName,levels = c("NR1D1","BMAL1","CLOCK","PER1","PER2","CRY1","CRY2"))) %>%
  ggplot(aes(x=tp,y=value,col=treat)) +
  geom_point(shape=1) +
  geom_line(data = df_predict %>% separate(sample,into = c("geneName","treat","dataset"),"_") %>% mutate(geneName = factor(geneName,levels = c("NR1D1","BMAL1","CLOCK","PER1","PER2","CRY1","CRY2"))) ,aes(x=tp,y=pred,col=treat),linetype = "dashed") +
  facet_wrap(dataset~geneName,scales = "free",ncol = 7) +
  theme_cowplot() +
  theme(strip.background = element_blank())+
  scale_color_manual(values = c("black","red"))
ggsave("../../out/plot/01_lanni_qpcr_gene_HigRes3_brain3.pdf",height = 9,width = 15)
