# libraries ---------------------------------------------------------------
library(tidyverse)
library(viridis)
library(cowplot)

# read data ---------------------------------------------------------------
df <- read_csv("../../data/240906_lanni_qpcr.csv")

# wrangling ---------------------------------------------------------------
# extract the info from the sample column
df_full <- df %>%
  mutate(treat = str_extract(sample_full,pattern = "CTR|DSS"),
         replicate = str_extract(sample_full,pattern = "^\\d+ \\w"),
         replicate2 = str_extract(replicate,"\\d+"))

# follow the method suggested by prof. Lanni
df_full_lanni <- df_full %>%
  group_by(time_point,gene,dataset,treat) %>%
  summarise(avg_exp = mean(exp,na.rm = T),
            sd_exp = sd(exp),
            n = n()) %>%
  ungroup()

# wrangling ---------------------------------------------------------------
# split the dataset per gene
list_df_full <- df_full %>%
  split(.$gene)

list_df_full_lanni <- df_full_lanni %>%
  split(.$gene)

# x <- list_df_full$CLOCK
# calculate the model parameters
df_lm <- lapply(list_df_full,function(x){
  model <- x %>%
    lm(exp~time_point+treat,data = .)
  
  test <- summary(model)
  broom::tidy(test) %>%
    filter(term == "treatDSS") %>%
    mutate(comp = "DSS vs CTR") %>%
    mutate(tissue = "colon")

}) %>%
  bind_rows(.id = "gene") %>%
  mutate(p.value.adj = p.adjust(p.value,method = "BH"))

df_lm %>%
  write_csv("../../out/table/03_df_lm_additive_colon.csv")  

# x <- list_df_full_lanni$CLOCK
df_lm_lanni <- lapply(list_df_full_lanni,function(x){
  model <- x %>%
    lm(avg_exp~time_point+treat,data = .)
  
  test <- summary(model)
  broom::tidy(test) %>%
    filter(term == "treatDSS") %>%
    mutate(comp = "DSS vs CTR") %>%
    mutate(tissue = "colon")
  
}) %>%
  bind_rows(.id = "gene") %>%
  mutate(p.value.adj = p.adjust(p.value,method = "BH"))

df_lm_lanni %>%
  write_csv("../../out/table/03_df_lm_additive_colon_lanni.csv")  

# lapply(list_df_full,function(x){
#   model <- x %>%
#     lm(exp~time_point*treat,data = .)
#   
#   test <- summary(model)
#   broom::tidy(test) %>%
#     filter(term == "treatDSS") %>%
#     mutate(comp = "DSS vs CTR")
#   
# }) %>%
#   bind_rows(.id = "gene") %>%
#   mutate(p.value.adj = p.adjust(p.value,method = "BH"))

# plot --------------------------------------------------------------------
df_full %>%
  ggplot(aes(x=treat,y=exp,col=time_point)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape=1,position = position_jitter(width = 0.2)) +
  facet_wrap(~gene,scales = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank()) +
  scale_color_viridis(option = "turbo")
ggsave("../../out/plot/03_sample_stat_colon.pdf",width = 10,height = 10)

df_full_lanni %>%
  ggplot(aes(x=treat,y=avg_exp,col=time_point)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(shape=1,position = position_jitter(width = 0.2)) +
  facet_wrap(~gene,scales = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank()) +
  scale_color_viridis(option = "turbo")
ggsave("../../out/plot/03_sample_stat_colon_lanni.pdf",width = 10,height = 10)
