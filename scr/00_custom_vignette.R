# AIM ---------------------------------------------------------------------
# load a sample dataset
# identify the main parameters, phase and perions using cellcycle
# use linear model to define the best trend and amplitude parameters
# try to plot the data with raw an prediction based on the model.
# use the model to increase the resolution of the data

# libraries ---------------------------------------------------------------
library(tidyverse)
library(MetaCycle)
library(ComplexHeatmap)
library(viridis)

# read in the data --------------------------------------------------------
# sample dataset
head(cycMouseLiverRNA[,1:5])

# save it as a table
write.csv(cycMouseLiverRNA, file="../../data/cycMouseLiverRNA.csv", row.names=FALSE)

# EDA ---------------------------------------------------------------------
# explore the dataset
df_long <- cycMouseLiverRNA %>%
  pivot_longer(names_to = "tp_cat",values_to = "value",-geneName) %>%
  mutate(tp = str_remove_all(tp_cat,"CT") %>%
           as.numeric())

# plot the data split by gene
df_long %>%
  ggplot(aes(x=tp,y=value))+geom_point(shape=1)+facet_wrap(~geneName,scales = "free")+theme_bw()+theme(strip.background = element_blank())

# estimate the period and phase with MetaCycle ----------------------------
# define the parameters
cyc <- meta2d(infile = "../../data/cycMouseLiverRNA.csv",
              filestyle = "csv",
              timepoints = 18:65,
              outputFile = FALSE,
              outRawData = TRUE)

# use the prediction to estimate the trend and amplitude paramters --------
# pull all the genes
genes <- cyc$meta$CycID

# define model, paramters and prediction for each gene
df_model <- lapply(genes, function(x){
  # pull the gene specific parameters
  params <- cyc$meta %>%
    filter(CycID %in% c(x)) %>%
    select(CycID,contains("meta"))
  
  # pull the data for the model
  df_test <- df_long %>%
    filter(geneName %in% x) %>%
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
    group_by(geneName) %>%
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
  facet_wrap(~geneName,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())

# plot the raw values as a heatmap
mat <- df_model %>%
  unnest(data) %>%
  mutate(tp = as.factor(tp)) %>%
  select(geneName,value,tp) %>%
  pivot_wider(names_from = tp,values_from = value) %>%
  column_to_rownames("geneName")

# scale the value by row (gene)
mat_scale <- t(scale(t(mat)))
Heatmap(mat_scale,cluster_columns = F,col = viridis::viridis(option = "turbo",n = 10))

# plot the heatmap usign the predicted values from the model
mat2 <- df_model %>%
  unnest(data) %>%
  mutate(tp = as.factor(tp)) %>%
  select(geneName,pred,tp) %>%
  pivot_wider(names_from = tp,values_from = pred) %>%
  column_to_rownames("geneName")

# scale the value by row
mat_scale2 <- t(scale(t(mat2)))
Heatmap(mat_scale2,cluster_columns = F,col = viridis::viridis(option = "turbo",n = 10))

# test increase the resolution --------------------------------------------
# now use each model to produce a modelled range per gene
df_predict <- lapply(genes,function(x){
  fit <- df_model %>%
    filter(geneName == x) %>%
    pull(model) %>%
    .[[1]]
  
  params <- df_model %>%
    filter(geneName == x) %>%
    pull(params) %>%
    .[[1]]
  
  newdata <- data.frame(tp = seq(from = 18,to=65,by = 0.1)) %>%
    mutate(phase = params$meta2d_phase,
           period = params$meta2d_period) %>%
    mutate(trend = tp-sum(tp)/nrow(.),
           amplitude = cos(2*pi*(tp-phase)/period)) %>%
    # add the prediction
    mutate(pred = predict(object = fit,newdata = .)) %>%
    mutate(geneName = x)
  
  return(newdata)
  
}) %>%
  bind_rows()

# plot the heatmap usign the model value at higher resolution
mat3 <- df_predict %>%
  mutate(tp = as.factor(tp)) %>%
  select(geneName,pred,tp) %>%
  pivot_wider(names_from = tp,values_from = pred) %>%
  column_to_rownames("geneName")

# scale the value by row
mat_scale3 <- t(scale(t(mat3)))
Heatmap(mat_scale3,cluster_columns = F,show_column_names = F,col = viridis::viridis(option = "turbo",n = 10))
