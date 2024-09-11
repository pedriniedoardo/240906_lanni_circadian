# libraries ---------------------------------------------------------------
library(tidyverse)
library(MetaCycle)

# read in the data --------------------------------------------------------
# sample dataset
head(cycMouseLiverRNA[,1:5])

# save it as a table
write.csv(cycMouseLiverRNA, file="../../data/cycMouseLiverRNA.csv", row.names=FALSE)

# explore the dataset
df_long <- cycMouseLiverRNA %>%
  pivot_longer(names_to = "tp_cat",values_to = "value",-geneName) %>%
  mutate(tp = str_remove_all(tp_cat,"CT") %>%
           as.numeric())

df_long %>%
  ggplot(aes(x=tp,y=value))+geom_point(shape=1)+facet_wrap(~geneName,scales = "free")+theme_bw()+theme(strip.background = element_blank())

# define the parameters
cyc <- meta2d(infile = "../../data/cycMouseLiverRNA.csv",
              filestyle = "csv",
              timepoints = 18:65,
              outputFile = FALSE,
              outRawData = TRUE)

# extract one sample variables
# cyc$meta %>%
#   filter(CycID %in% c("Rorc_1425792_a_at")) %>%
#   select(CycID,contains("CT")) %>%
#   pivot_longer(names_to = "tp_cat",values_to = "value",-CycID) %>%
#   mutate(tp = str_remove_all(tp_cat,"CT") %>%
#            as.numeric()) %>%
#   ggplot(aes(x=tp,y=value))+geom_point(shape=1)+facet_wrap(~CycID,scales = "free")+theme_bw()+theme(strip.background = element_blank())

params <- cyc$meta %>%
  filter(CycID %in% c("Rorc_1425792_a_at")) %>%
  select(CycID,contains("meta"))

# alternative model suggested in the metacycle package
# yi = B + trend * (ti-sum(y)/n) + amplitude *(cos(2*pi*(ti-phase)/period))
df_long %>%
  filter(geneName %in% c("Rorc_1425792_a_at")) %>%
  mutate(pred = params$meta2d_Base + params$meta2d_AMP*(cos(2*pi*(tp-params$meta2d_phase)/params$meta2d_period))) %>%
  ggplot(aes(x=tp,y=value)) +
  geom_point(shape=1) +
  geom_line(aes(x=tp,y=pred)) +
  facet_wrap(~geneName,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())

# define the model
# yi = B + trend * (ti-sum(t)/n) + amplitude *(cos(2*pi*(ti-phase)/period))
df_test <- df_long %>%
  filter(geneName %in% c("Rorc_1425792_a_at")) %>%
  mutate(phase = params$meta2d_phase,
         period = params$meta2d_period) %>%
  mutate(trend = tp-sum(tp)/nrow(.),
         amplitude = cos(2*pi*(tp-phase)/period))

model.matrix(~ trend+amplitude, data=df_test)
fit <- lm(value ~ trend+amplitude, data=df_test)
summary(fit)

# make the prediction using fit
df_test %>%
  mutate(pred = predict(object = fit,newdata = df_test)) %>%
  ggplot(aes(x=tp,y=value)) +
  geom_point(shape=1) +
  geom_line(aes(x=tp,y=pred)) +
  facet_wrap(~geneName,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())

# -------------------------------------------------------------------------
# try to include more genes
genes <- cyc$meta$CycID
x<- genes[[1]]

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

# plot all the genes with the relative model
df_model %>%
  unnest(data) %>%
  ggplot(aes(x=tp,y=value)) +
  geom_point(shape=1) +
  geom_line(aes(x=tp,y=pred),col="red",linetype = "dashed") +
  facet_wrap(~geneName,scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank())

# plot the heatmap with row values
mat <- df_model %>%
  unnest(data) %>%
  mutate(tp = as.factor(tp)) %>%
  select(geneName,value,tp) %>%
  pivot_wider(names_from = tp,values_from = value) %>%
  column_to_rownames("geneName")

# scale the value by row
mat_scale <- t(scale(t(mat)))

library(ComplexHeatmap)
Heatmap(mat_scale,cluster_columns = F)

# plot the heatmap usign the predicted values from the model
mat2 <- df_model %>%
  unnest(data) %>%
  mutate(tp = as.factor(tp)) %>%
  select(geneName,pred,tp) %>%
  pivot_wider(names_from = tp,values_from = pred) %>%
  column_to_rownames("geneName")

# scale the value by row
mat_scale2 <- t(scale(t(mat2)))

library(ComplexHeatmap)
Heatmap(mat_scale2,cluster_columns = F)

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

library(ComplexHeatmap)
Heatmap(mat_scale3,cluster_columns = F,show_column_names = F,col = viridis::viridis(option = "turbo",n = 10))

# help page ---------------------------------------------------------------

# save the objects to tables ----------------------------------------------
# sample dataset available in the package
head(cycSimu4h2d)
head(cycMouseLiverRNA)
head(cycYeastCycle)

# save the objects as a csv
write.csv(cycSimu4h2d, file="../../data/cycSimu4h2d.csv", row.names=FALSE)
write.csv(cycMouseLiverRNA, file="../../data/cycMouseLiverRNA.csv", row.names=FALSE)
write.csv(cycYeastCycle, file="../../data/cycYeastCycle.csv", row.names=FALSE)

# sample code to save as a tsv
write.table(cycMouseLiverProtein, file="../../data/cycMouseLiverProtein.txt",sep="\t", quote=FALSE, row.names=FALSE)

# explore the objects -----------------------------------------------------

# cycMouseLiverProtein
# analyze 'cycMouseLiverRNA.csv' with JTK_CYCLE
# check the structure of the dataset
dim(cycMouseLiverProtein)
head(cycMouseLiverProtein)

# notice that the dataset has some replicates per time points coded in the table
df_long_cycMouseLiverProtein <- cycMouseLiverProtein %>%
  pivot_longer(names_to = "sample",values_to = "value",-geneSymbol) %>%
  mutate(tp_cat = str_extract(sample,"CT\\d+")) %>%
  mutate(rep = str_extract(sample,"Rep\\d+")) %>%
  mutate(tp = str_remove_all(tp_cat,"CT") %>%
           as.numeric())

# summarise the trend per replicate
df_long_summary_cycMouseLiverProtein <- df_long_cycMouseLiverProtein %>%
  group_by(geneSymbol,tp) %>%
  summarise(avg_value = mean(value,na.rm = T))

# plot the data
df_long_cycMouseLiverProtein %>%
  ggplot(aes(x=tp,y=value))+
  geom_point(shape=1)+
  geom_line(data = df_long_summary_cycMouseLiverProtein,aes(x=tp,y=avg_value),col="red")+
  facet_wrap(~geneSymbol,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank())
ggsave("../../out/plot/01_cycMouseLiverProtein.pdf",width = 9,height = 6)

# cycSimu4h2d
# check the structure of the dataset
dim(cycSimu4h2d)
head(cycSimu4h2d)

# this dataset does not seems to have replicated
df_long_cycSimu4h2d <- cycSimu4h2d %>%
  pivot_longer(names_to = "tp_cat",values_to = "value",-curveID) %>%
  mutate(tp = str_remove_all(tp_cat,"CT") %>%
           as.numeric())

# summarise the trend per replicate
df_long_summary_cycSimu4h2d <- df_long_cycSimu4h2d %>%
  group_by(curveID,tp) %>%
  summarise(avg_value = mean(value,na.rm = T))

# plot the data
df_long_cycSimu4h2d %>%
  ggplot(aes(x=tp,y=value))+
  geom_point(shape=1)+
  geom_line(data = df_long_summary_cycSimu4h2d,aes(x=tp,y=avg_value),col="red")+
  facet_wrap(~curveID,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank())
ggsave("../../out/plot/01_cycSimu4h2d.pdf",width = 15,height = 12)

# cycYeastCycle
# check the structure of the dataset
dim(cycYeastCycle)
head(cycYeastCycle)

# this dataset does not seems to have replicated
df_long_cycYeastCycle <- cycYeastCycle %>%
  pivot_longer(names_to = "sample",values_to = "value",-ID_REF) %>%
  mutate(tp = str_remove_all(sample,"recov_|min") %>%
           as.numeric())

# summarise the trend per replicate
df_long_summary_cycYeastCycle <- df_long_cycYeastCycle %>%
  group_by(ID_REF,tp) %>%
  summarise(avg_value = mean(value,na.rm = T))

# plot the data
df_long_cycYeastCycle %>%
  ggplot(aes(x=tp,y=value))+
  geom_point(shape=1)+
  geom_line(data = df_long_summary_cycYeastCycle,aes(x=tp,y=avg_value),col="red")+
  facet_wrap(~ID_REF,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank())
ggsave("../../out/plot/01_cycYeastCycle.pdf",width = 12,height = 9)

# analysis ----------------------------------------------------------------
# meta2d(infile="cycMouseLiverRNA.csv",
#        filestyle="csv",
#        outdir="example",
#        timepoints=18:65,
#        cycMethod="JTK",
#        outIntegration="noIntegration")

# analyze 'cycMouseLiverProtein.txt' with JTK_CYCLE and Lomb-Scargle. the function is not returning anything, therefore there is no need to save the output.
meta2d(infile = "../../data/cycMouseLiverProtein.txt",
       filestyle = "txt",
       outdir = "../../out/table/",
       timepoints = rep(seq(0, 45, by=3), each=3),
       cycMethod = c("JTK","LS"),
       outIntegration = "noIntegration")

# meta2d(infile = "../../data/cycMouseLiverProtein.txt",
#        filestyle = "txt",
#        outdir = "../../out/table/",
#        timepoints = rep(seq(0, 45, by=3), each=3))

# analyze 'cycSimu4h2d.csv' with ARSER, JTK_CYCLE and Lomb-Scargle and output integration file with analysis results from each method
meta2d(infile = "../../data/cycSimu4h2d.csv",
       filestyle = "csv",
       outdir = "../../out/table",
       timepoints = "Line1")

# analyze 'cycYeastCycle.csv' with ARSER, JTK_CYCLE and Lomb-Scargle to detect transcripts associated with cell cycle, and only output integration file
meta2d(infile = "../../data/cycYeastCycle.csv",
       filestyle = "csv",
       outdir = "../../out/table",
       minper = 80,
       maxper = 96,
       timepoints = seq(2, 162, by=16),
       outIntegration = "onlyIntegration",
       ARSdefaultPer = 85,
       outRawData = TRUE)

# return analysis results instead of output them into files
cyc <- meta2d(infile = "../../data/cycYeastCycle.csv",
              filestyle = "csv",
              minper = 80,
              maxper = 96,
              timepoints = seq(2, 162, by = 16),
              outputFile = FALSE,
              ARSdefaultPer = 85,
              outRawData = TRUE)
head(cyc$ARS)
head(cyc$JTK)
head(cyc$LS)
head(cyc$meta)

# how to leverage the paremeter to plot a model



# sample from chinese blog ------------------------------------------------
# https://abego.cn/2019/05/31/the-rule-of-gene-expression-in-the-day-and-nigth/
cyc <- meta2d(infile = "../../data/cycYeastCycle.csv",filestyle="csv",
              minper = 80, maxper = 96,
              timepoints = seq(2, 162, by=16),
              outputFile = FALSE,
              ARSdefaultPer = 85,
              outRawData = TRUE)
head(cyc$ARS)
head(cyc$JTK)
head(cyc$LS)
head(cyc$meta)

head(cyc$JTK)
dat_raw <- cyc$JTK
dat <- dat_raw %>% filter(ADJ.P < 0.05 )
rownames(dat) <- dat[,1]
name <- colnames(dat)
name_t <- str_detect(name, '\\d')
name_f <- !str_detect(name, '\\d')
name_exp <- name[name_t]
name_wave <- name[name_f]
dat_exp <- dat[,name_exp]
dat_wave <- dat[,name_f]
dat_wave <- dat_wave[,-1]

require(pheatmap)
require(RColorBrewer)
color <- colorRampPalette(rev(brewer.pal(n = 3, name = "BuPu")))(2)
p3 <- pheatmap(dat_exp,
               color = color, scale = "row", cluster_cols = F,
               border_color = NA,
               cellheight =10,
               cellwidth = 10 ,
               legend =F, 
               show_rownames = T,
               treeheight_row = 0)

# require(reshape2)
dat_exp <- tibble::rownames_to_column(dat_exp)
dat_plot <- dat_exp %>%
  pivot_longer(names_to = "variable",values_to = "value",-rowname)

p1 <- ggplot(dat_plot, aes(x = variable, y = value, color = variable)) + 
  geom_point() +
  facet_grid(rowname ~. ,scales="free_y")

p2 <- p1 +theme_bw()+theme(legend.position="none") +
  labs(title = "the circadian rhythm gene expression",x = "Time", y = "FPKM") +
  theme(axis.text.x=element_text(face="bold", size=10, angle = 45, vjust = 0.5),
        axis.text.y=element_text(face="bold", size=6),
        axis.title.y=element_text(size=8, face="bold"),
        strip.text.y = element_text(size = 6, face ="bold", angle = 0 ))
