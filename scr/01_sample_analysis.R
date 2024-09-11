# libraries ---------------------------------------------------------------
library(tidyverse)
library(MetaCycle)
library(ComplexHeatmap)
library(viridis)

# read data ---------------------------------------------------------------
df <- read_csv("../../data/240906_lanni_qpcr.csv")

# wrangling ---------------------------------------------------------------
# extract the info from the sample column
df_full <- df %>%
  mutate(treat = str_extract(sample_full,pattern = "CTR|DSS"),
         replicate = str_extract(sample_full,pattern = "^\\d+ \\w"),
         replicate2 = str_extract(replicate,"\\d+"))

# shape the data to be accepted by MetaCycle
df_full2 <- df_full %>%
  mutate(sample = paste0(gene,"_",treat)) %>%
  mutate(tp_rep = paste0("T_",time_point,"_",replicate2)) %>%
  select(sample,tp_rep,exp) %>%
  pivot_wider(names_from = tp_rep,values_from = exp)

# save it as a table
write.csv(df_full2, file="../../data/240906_lanni_qpcr_clean.csv", row.names=FALSE)

# EDA ---------------------------------------------------------------------
# explore the dataset
df_long <- df_full2 %>%
  pivot_longer(names_to = "tp_cat",values_to = "value",-sample) %>%
  separate(sample,into = c("geneName","treat"),sep="_",remove = F) %>%
  separate(tp_cat,into = c("X","tp","rep_id"),sep="_",remove = F) %>%
  dplyr::select(-X) %>%
  mutate(tp = as.numeric(tp))

# plot the data split by gene
df_long %>%
  ggplot(aes(x=tp,y=value,col=treat))+geom_point(shape=1)+facet_wrap(~geneName,scales = "free")+theme_bw()+theme(strip.background = element_blank())+scale_color_manual(values = c("black","red"))

# estimate the period and phase with MetaCycle ----------------------------
# define the parameters
cyc <- meta2d(infile = "../../data/240906_lanni_qpcr_clean.csv",
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

# plotting the data -------------------------------------------------------
# plot all the genes with the relative model
df_model %>%
  unnest(data) %>%
  ggplot(aes(x=tp,y=value)) +
  geom_point(shape=1) +
  geom_line(aes(x=tp,y=pred),col="red",linetype = "dashed") +
  facet_wrap(~sample,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())
ggsave("../../out/plot/01_lanni_qpcr_sample.pdf",height = 12,width = 15)

# plot all the genes colored by treat
df_model %>%
  unnest(data) %>%
  ggplot(aes(x=tp,y=value,col=treat)) +
  geom_point(shape=1) +
  geom_line(aes(x=tp,y=pred,col=treat),linetype = "dashed") +
  facet_wrap(~geneName,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())+
  scale_color_manual(values = c("black","red"))
ggsave("../../out/plot/01_lanni_qpcr_gene.pdf",height = 9,width = 9)

# add some jittering to the points
df_model %>%
  unnest(data) %>%
  ggplot(aes(x=tp,y=value,col=treat)) +
  geom_point(shape=1,position = position_jitter(width = 1)) +
  geom_line(aes(x=tp,y=pred,col=treat),linetype = "dashed") +
  facet_wrap(~geneName,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())+
  scale_color_manual(values = c("black","red"))
ggsave("../../out/plot/01_lanni_qpcr_gene_jitter.pdf",height = 9,width = 9)

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
pdf("../../out/plot/01_lanni_qpcr_heatmap_raw_sample.pdf",width = 4,height = 5)
# Heatmap(mat_scale,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale,cluster_columns = F,cluster_rows = F)
dev.off()

pdf("../../out/plot/01_lanni_qpcr_heatmap_raw_sample2.pdf",width = 4,height = 5)
# Heatmap(mat_scale,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale,cluster_columns = F,cluster_rows = T)
dev.off()

# scale the value by gene. keep the different samples alone (gene)
# 
mat_scale_gene <- mat %>%
  rownames_to_column("sample") %>%
  pivot_longer(names_to = "tp",values_to = "exp",-sample) %>%
  separate(sample,into = c("geneName","treat"),sep="_",remove = F) %>%
  group_by(geneName) %>%
  mutate(exp_scaled = (exp-mean(exp))/sd(exp)) %>%
  ungroup() %>%
  select(sample,tp,exp_scaled) %>%
  pivot_wider(names_from = tp,values_from = exp_scaled) %>%
  column_to_rownames("sample")%>%
  as.matrix()

pdf("../../out/plot/01_lanni_qpcr_heatmap_raw_gene.pdf",width = 4,height = 5)
# Heatmap(mat_scale_gene,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale_gene,cluster_columns = F,cluster_rows = F)
dev.off()

pdf("../../out/plot/01_lanni_qpcr_heatmap_raw_gene2.pdf",width = 4,height = 5)
# Heatmap(mat_scale_gene,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale_gene,cluster_columns = F,cluster_rows = T)
dev.off()

# plot the heatmap usign the predicted values from the model
mat2 <- df_model %>%
  unnest(data) %>%
  mutate(tp = as.factor(tp)) %>%
  # select(sample,value,tp) %>%
  group_by(sample,tp) %>%
  summarise(avg_pred = mean(pred,na.rm = T)) %>%
  pivot_wider(names_from = tp,values_from = avg_pred) %>%
  column_to_rownames("sample")

# scale by sample
mat_scale2 <- t(scale(t(mat2)))
pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_sample.pdf",width = 4,height = 5)
# Heatmap(mat_scale,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale2,cluster_columns = F,cluster_rows = F)
dev.off()

pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_sample2.pdf",width = 4,height = 5)
# Heatmap(mat_scale,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale2,cluster_columns = F,cluster_rows = T)
dev.off()

# scale the value by gene. keep the different samples alone (gene)
# 
mat_scale_gene2 <- mat2 %>%
  rownames_to_column("sample") %>%
  pivot_longer(names_to = "tp",values_to = "exp",-sample) %>%
  separate(sample,into = c("geneName","treat"),sep="_",remove = F) %>%
  group_by(geneName) %>%
  mutate(exp_scaled = (exp-mean(exp))/sd(exp)) %>%
  ungroup() %>%
  select(sample,tp,exp_scaled) %>%
  pivot_wider(names_from = tp,values_from = exp_scaled) %>%
  column_to_rownames("sample")%>%
  as.matrix()

pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_gene.pdf",width = 4,height = 5)
# Heatmap(mat_scale_gene,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale_gene2,cluster_columns = F,cluster_rows = F)
dev.off()

pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_gene2.pdf",width = 4,height = 5)
# Heatmap(mat_scale_gene,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale_gene2,cluster_columns = F,cluster_rows = T)
dev.off()

# test increase the resolution --------------------------------------------
# now use each model to produce a modelled range per gene
# x <- genes[1]
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
pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_sample_HigRes.pdf",width = 4,height = 5)
# Heatmap(mat_scale,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale3,cluster_columns = F,cluster_rows = F,show_column_names = F)
dev.off()

pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_sample_HigRes2.pdf",width = 4,height = 5)
# Heatmap(mat_scale,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale3,cluster_columns = F,cluster_rows = T,show_column_names = F)
dev.off()

# scale the value by gene. keep the different samples alone (gene)
# 
mat_scale_gene3 <- mat3 %>%
  rownames_to_column("sample") %>%
  pivot_longer(names_to = "tp",values_to = "exp",-sample) %>%
  separate(sample,into = c("geneName","treat"),sep="_",remove = F) %>%
  group_by(geneName) %>%
  mutate(exp_scaled = (exp-mean(exp))/sd(exp)) %>%
  ungroup() %>%
  select(sample,tp,exp_scaled) %>%
  pivot_wider(names_from = tp,values_from = exp_scaled) %>%
  column_to_rownames("sample")%>%
  as.matrix()

pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_gene_HigRes.pdf",width = 4,height = 5)
# Heatmap(mat_scale_gene,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale_gene3,cluster_columns = F,cluster_rows = F,show_column_names = F)
dev.off()

pdf("../../out/plot/01_lanni_qpcr_heatmap_prediction_gene_HigRes2.pdf",width = 4,height = 5)
# Heatmap(mat_scale_gene,cluster_columns = F,cluster_rows = F,col = viridis::viridis(option = "turbo",n = 10))
Heatmap(mat_scale_gene3,cluster_columns = F,cluster_rows = T,show_column_names = F)
dev.off()