source("F:/professdev/r_projects/gsr_analysis/code/m_biome_packages.R")

conflicts_prefer(dplyr::filter) # this line will prefer filter from dplyr
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)
conflicts_prefer(ggpubr::mutate)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(dplyr::desc)

#################### 16S Data for alpha diversity analysis ##################

# A. 16S data general

#metadata

cfvb_lsu_meta_p3 <- read_tsv("data/cfvb_16s_gsr_metadata.txt") |>
  select('sample_id'='sample-id','bagging'='treat_one', 
         'insecticide'='treat_two', variety,period) |> 
  filter(sample_id != "ctrl") |>   #drop the control sample
  filter(period == "3") |> 
  mutate(variety = case_when(variety == "cabfranc" ~ "CF",
                             variety == "vidal" ~ "VB"),
         bagging = case_when(bagging == "NNet" ~ "Nnet",
                             bagging == "Net" ~ "Net"),
         insecticide = case_when(insecticide == "MMax" ~ "Mmax",
                                 insecticide == "NMax" ~ "Nmax"))

cfvb_lsu_meta_p3$variety <- as.factor(cfvb_lsu_meta_p3$variety)
cfvb_lsu_meta_p3$bagging <- as.factor(cfvb_lsu_meta_p3$bagging)
cfvb_lsu_meta_p3$insecticide <- as.factor(cfvb_lsu_meta_p3$insecticide)

glimpse(cfvb_lsu_meta_p3)

#evenness
cfvb_lsu_even <- read_tsv("data/cfvb_16s_evenness.tsv") |>  
  rename("sample_id" = "...1", 'evenness'='pielou_evenness') #renames col 1

#faith
cfvb_lsu_faith <- read_tsv("data/cfvb_16s_faith_pd.tsv") |>  
  rename("sample_id" = "#SampleID") #renames col 1

#observed_otus
cfvb_lsu_obs_otus <- read_tsv("data/cfvb_16s_observed_otus.tsv") |>  
  rename("sample_id" = "...1", "obs_otus"="observed_features") #rename columns

#shannon
cfvb_lsu_shan <- read_tsv("data/cfvb_16s_shannon.tsv") |>  
  rename("sample_id" = "...1", "shannon"="shannon_entropy")

#otu
cfvb_lsu_otus <- read_tsv("data/cfvb_16s_rarefied_table_transp.txt", skip = 1) |> 
  rename("sample_id"="#OTU ID")

#ace
cfvb_lsu_ace <- read_tsv("data/cfvb_16s_ace.tsv") |>  
  rename("sample_id" = "...1") #rename columns

#chao1
cfvb_lsu_chao1 <- read_tsv("data/cfvb_16s_chao1.tsv") |>  
  rename("sample_id" = "...1")

#merge metadata, ace, chao1, evenness, faith, observed_otus, shannon and otus

cfvb_lsu_merged_p3 <- inner_join(cfvb_lsu_meta_p3, 
                                 cfvb_lsu_chao1, by='sample_id') %>%
  inner_join(., cfvb_lsu_faith, by='sample_id') %>%
  inner_join(., cfvb_lsu_obs_otus, by='sample_id') %>%
  inner_join(., cfvb_lsu_shan, by='sample_id') %>%
  inner_join(., cfvb_lsu_otus, by='sample_id') 
  

#cfvb_lsu_merged_p3 <- replace(cfvb_lsu_merged_p3,is.na(cfvb_lsu_merged_p3),0)

######################## Overall 16s Analysis ##################################
#Create 2x2 plot environment 
par(mfrow = c(2, 2))
hist(cfvb_lsu_merged_p3$chao1, main="Chao1 diversity", xlab="", breaks=10, 
     border="#000000", col="#99DDFF", las=1)
# hist(cfvb_lsu_merged_p3$ace, main="Ace diversity", xlab="", breaks=10, 
#      border="#000000", col="#BBCC33", las=1)
hist(cfvb_lsu_merged_p3$shannon, main="Shannon diversity", xlab="", breaks=10, 
     border="#000000", col="#FFAABB", las=1)
# hist(cfvb_lsu_merged_p3$evenness, main="Pielou's evenness", xlab="", breaks=10, 
#      border="#000000", col="#44BB99", las=1)
hist(cfvb_lsu_merged_p3$faith_pd, main="Phylegenetic diversity", xlab="", 
     breaks=10, border="#000000", col="#DDDDDD", las=1)
hist(cfvb_lsu_merged_p3$obs_otus, main="Community richness", xlab="", 
     breaks=10, border="#000000", col="#EEDD88", las=1)


# To test for normality

shapiro.test(cfvb_lsu_merged_p3$chao1)      # normally dist
shapiro.test(cfvb_lsu_merged_p3$shannon)
shapiro.test(cfvb_lsu_merged_p3$faith_pd)   # normally dist
shapiro.test(cfvb_lsu_merged_p3$obs_otus)   # normally dist


# a. Anova test for normally distributed data

# i-chao1
cfvb.lsu.chao1.ins.aov <- aov(chao1 ~ insecticide, data=cfvb_lsu_merged_p3)
summary(cfvb.lsu.chao1.ins.aov)
cfvb.lsu.chao1.bag.aov <- aov(chao1 ~ bagging, data=cfvb_lsu_merged_p3)
summary(cfvb.lsu.chao1.bag.aov)
cfvb.lsu.chao1.var.aov <- aov(chao1 ~ variety, data=cfvb_lsu_merged_p3)
summary(cfvb.lsu.chao1.var.aov)

# v-shannon
kruskal.test(shannon ~ insecticide, data=cfvb_lsu_merged_p3)
kruskal.test(shannon ~ bagging, data=cfvb_lsu_merged_p3)
kruskal.test(shannon ~ variety, data=cfvb_lsu_merged_p3)

# iii-Faith
cfvb.lsu.fd.ins.aov <- aov(faith_pd ~ insecticide, data=cfvb_lsu_merged_p3)
summary(cfvb.lsu.fd.ins.aov)
cfvb.lsu.fd.bag.aov <- aov(faith_pd ~ bagging, data=cfvb_lsu_merged_p3)
summary(cfvb.lsu.fd.bag.aov)
cfvb.lsu.fd.var.aov <- aov(faith_pd ~ variety, data=cfvb_lsu_merged_p3)
summary(cfvb.lsu.fd.var.aov)

# iv-Observed_otus
cfvb.lsu.obs.ins.aov <- aov(obs_otus ~ insecticide, data=cfvb_lsu_merged_p3)
summary(cfvb.lsu.obs.ins.aov)
cfvb.lsu.obs.bag.aov <- aov(obs_otus ~ bagging, data=cfvb_lsu_merged_p3)
summary(cfvb.lsu.obs.bag.aov)
cfvb.lsu.obs.var.aov <- aov(obs_otus ~ variety, data=cfvb_lsu_merged_p3)
summary(cfvb.lsu.obs.var.aov)


# Plotting relationships 16S General

# boxplots

my_labs <- c("shannon" = "shannon",
             "chao1" = "chao1",
             "faith_pd" = "faith",
             "obs_otus" = "observed otus") # labels for y-axis of boxplots

# a. shannon, chao1, otus and faith by insecticide
cfvb_lsu_alpha_div_ins_plot_p3 <- cfvb_lsu_merged_p3 |> 
  select(-period) |> 
  pivot_longer(cols = 5:8, names_to = 'indice', values_to = 'value') |>
  mutate(bagging=case_when(bagging=="Net"~"+Net",
                           bagging=="Nnet"~"-Net"),
         insecticide=case_when(insecticide=="Nmax"~"-Max",
                               insecticide=="Mmax"~"+Max")) |> 
  mutate(indice=factor(indice, levels = c("shannon", "chao1", 
                                          "obs_otus", "faith_pd"))) |>  
  select(insecticide, indice, value) |>  
  group_by(insecticide, indice) %>% 
  ggplot(., aes(x=insecticide, y=value)) + 
  geom_boxplot(aes(fill=insecticide),size=.5, outlier.shape = NA) +
  scale_fill_manual(values = c("#FFFFFF","#BFBFBF")) +
  facet_wrap(~indice, scales = "free",strip.position = "left",
             labeller = labeller(indice = my_labs)) +
  stat_compare_means(label = "p.signif", size=4) +
  labs(x = NULL, y = NULL, #removes the x- and y-axis labels
       title = " "
  ) +
  theme_bw() + theme(legend.position = "None") +
  theme(
    plot.title = element_text(hjust=0.5),
    strip.placement = "outside", #sends the panel labels from left to outside 
    #of y-axis
    strip.background = element_blank(), #removes background of strip
    strip.text = element_text(size=9,colour="#000000",family="sans"),
    panel.background = element_blank(),
    axis.title.y = element_text(size=9,family="sans",color="#000000"),
    axis.text = element_text(size=8,family="sans",color="#000000"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(hjust=.5,vjust=1,family="sans",color="#000000"))

ggsave("figures/p3_figures/cfvb_lsu_alpha_div_ins_plot_p3.pdf", 
       cfvb_lsu_alpha_div_ins_plot_p3, height = 4.5, width = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_alpha_div_ins_plot_p3.tif", 
       cfvb_lsu_alpha_div_ins_plot_p3, height = 4.5, width = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_alpha_div_ins_plot_p3.svg", 
       cfvb_lsu_alpha_div_ins_plot_p3, height = 4.5, width = 4, dpi = 300)

# b. shannon, chao1, otus and faith by bagging
cfvb_lsu_alpha_div_bag_plot_p3 <- cfvb_lsu_merged_p3 |> 
  select(-period) |> 
  pivot_longer(cols = 5:8, names_to = 'indice', values_to = 'value') |>
  mutate(bagging=case_when(bagging=="Net"~"+Net",
                           bagging=="Nnet"~"-Net"),
         insecticide=case_when(insecticide=="Mmax"~"+Max",
                               insecticide=="Nmax"~"-Max")) |> 
  mutate(indice=factor(indice, levels = c("shannon", "chao1", 
                                          "obs_otus", "faith_pd"))) |>  
  select(bagging, indice, value) |>  
  group_by(bagging, indice) %>% 
  ggplot(., aes(x=bagging, y=value)) + 
  geom_boxplot(aes(fill=bagging),size=.5, outlier.shape = NA) +
  scale_fill_manual(values = c("#FFFFFF","#BFBFBF")) +
  facet_wrap(~indice, scales = "free",strip.position = "left",
             labeller = labeller(indice = my_labs)) +
  stat_compare_means(label = "p.signif", size=4) +
  labs(x = NULL, y = NULL, #removes the x- and y-axis labels
       title = " "
  ) +
  theme_bw() + theme(legend.position = "None") +
  theme(
    plot.title = element_text(hjust=0.5),
    strip.placement = "outside", #sends the panel labels from left to outside 
    #of y-axis
    strip.background = element_blank(), #removes background of strip
    strip.text = element_text(size=9,colour="#000000",family="sans"),
    panel.background = element_blank(),
    axis.title.y = element_text(size=9,family="sans",color="#000000"),
    axis.text = element_text(size=8,family="sans",color="#000000"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(hjust=.5,vjust=1,family="sans",color="#000000"))

ggsave("figures/p3_figures/cfvb_lsu_alpha_div_bag_plot_p3.pdf", 
       cfvb_lsu_alpha_div_bag_plot_p3, height = 4.5, width = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_alpha_div_bag_plot_p3.tif", 
       cfvb_lsu_alpha_div_bag_plot_p3, height = 4.5, width = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_alpha_div_bag_plot_p3.svg", 
       cfvb_lsu_alpha_div_bag_plot_p3, height = 4.5, width = 4, dpi = 300)


# c. shannon, chao1, otus and faith by variety
cfvb_lsu_alpha_div_var_plot_p3 <- cfvb_lsu_merged_p3 |> 
  select(-period) |> 
  pivot_longer(cols = 5:8, names_to = 'indice', values_to = 'value') |>
  mutate(indice=factor(indice, levels = c("shannon", "chao1", 
                                          "obs_otus", "faith_pd"))) |>  
  select(variety, indice, value) |>  
  group_by(variety, indice) %>% 
  ggplot(., aes(x=variety, y=value)) + 
  geom_boxplot(aes(fill=variety),size=.5, outlier.shape = NA) +
  scale_fill_manual(values = c("#FFFFFF","#BFBFBF")) +
  facet_wrap(~indice, scales = "free",strip.position = "left",
             labeller = labeller(indice = my_labs)) +
  stat_compare_means(label = "p.signif", size=4) +
  labs(x = NULL, y = NULL, #removes the x- and y-axis labels
       title = " "
  ) +
  theme_bw() + theme(legend.position = "None") +
  theme(
    plot.title = element_text(hjust=0.5),
    strip.placement = "outside", #sends the panel labels from left to outside 
    #of y-axis
    strip.background = element_blank(), #removes background of strip
    strip.text = element_text(size=9,colour="black",family="sans"),
    panel.background = element_blank(),
    axis.title.y = element_text(size=9,family="sans",color="#000000"),
    axis.text = element_text(size=8,family="sans",color="black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(hjust=.5,vjust=1,family="sans",color="#000000"))

ggsave("figures/p3_figures/cfvb_lsu_alpha_div_var_plot_p3.pdf", 
       cfvb_lsu_alpha_div_var_plot_p3, height = 4.5, width = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_alpha_div_var_plot_p3.tif", 
       cfvb_lsu_alpha_div_var_plot_p3, height = 4.5, width = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_alpha_div_var_plot_p3.svg", 
       cfvb_lsu_alpha_div_var_plot_p3, height = 4.5, width = 4, dpi = 300)


########################  16s Analysis for CF #################################

cf_lsu_merged_p3 <- cfvb_lsu_merged_p3  |> 
  filter(variety=="CF")

#Create 2x2 plot environment 
par(mfrow = c(2, 2))
hist(cf_lsu_merged_p3$chao1, main="Chao1 diversity", xlab="", breaks=10, 
     border="#000000", col="#99DDFF", las=1)
# hist(cfvb_lsu_merged_p3$ace, main="Ace diversity", xlab="", breaks=10, 
#      border="#000000", col="#BBCC33", las=1)
hist(cf_lsu_merged_p3$shannon, main="Shannon diversity", xlab="", breaks=10, 
     border="#000000", col="#FFAABB", las=1)
# hist(cfvb_lsu_merged_p3$evenness, main="Pielou's evenness", xlab="", breaks=10, 
#      border="#000000", col="#44BB99", las=1)
hist(cf_lsu_merged_p3$faith_pd, main="Phylegenetic diversity", xlab="", 
     breaks=10, border="#000000", col="#DDDDDD", las=1)
hist(cf_lsu_merged_p3$obs_otus, main="Community richness", xlab="", 
     breaks=10, border="#000000", col="#EEDD88", las=1)


# To test for normality

shapiro.test(cf_lsu_merged_p3$chao1)      # normally dist
shapiro.test(cf_lsu_merged_p3$shannon)    # normally dist
shapiro.test(cf_lsu_merged_p3$faith_pd)   # normally dist
shapiro.test(cf_lsu_merged_p3$obs_otus)   # normally dist


# i-chao1
cf.lsu.chao1.ins.aov <- aov(chao1 ~ insecticide, data=cf_lsu_merged_p3)
summary(cf.lsu.chao1.ins.aov)
cf.lsu.chao1.bag.aov <- aov(chao1 ~ bagging, data=cf_lsu_merged_p3)
summary(cf.lsu.chao1.bag.aov)

# v-shannon
cf.lsu.shan.ins.aov <- aov(shannon ~ insecticide, data=cf_lsu_merged_p3)
summary(cf.lsu.shan.ins.aov)
cf.lsu.shan.bag.aov <- aov(shannon ~ bagging, data=cf_lsu_merged_p3)
summary(cf.lsu.shan.bag.aov)

# iii-Faith
cf.lsu.fd.ins.aov <- aov(faith_pd ~ insecticide, data=cf_lsu_merged_p3)
summary(cf.lsu.fd.ins.aov)
cf.lsu.fd.bag.aov <- aov(faith_pd ~ bagging, data=cf_lsu_merged_p3)
summary(cf.lsu.fd.bag.aov)

# iv-Observed_otus
cf.lsu.obs.ins.aov <- aov(obs_otus ~ insecticide, data=cf_lsu_merged_p3)
summary(cf.lsu.obs.ins.aov)
cf.lsu.obs.bag.aov <- aov(obs_otus ~ bagging, data=cf_lsu_merged_p3)
summary(cf.lsu.obs.bag.aov)



########################  16s Analysis for VB #################################

vb_lsu_merged_p3 <- cfvb_lsu_merged_p3  |> 
  filter(variety=="VB")

#Create 2x2 plot environment 
par(mfrow = c(2, 2))
hist(vb_lsu_merged_p3$chao1, main="Chao1 diversity", xlab="", breaks=10, 
     border="#000000", col="#99DDFF", las=1)
# hist(cfvb_lsu_merged_p3$ace, main="Ace diversity", xlab="", breaks=10, 
#      border="#000000", col="#BBCC33", las=1)
hist(vb_lsu_merged_p3$shannon, main="Shannon diversity", xlab="", breaks=10, 
     border="#000000", col="#FFAABB", las=1)
# hist(cfvb_lsu_merged_p3$evenness, main="Pielou's evenness", xlab="", breaks=10, 
#      border="#000000", col="#44BB99", las=1)
hist(vb_lsu_merged_p3$faith_pd, main="Phylegenetic diversity", xlab="", 
     breaks=10, border="#000000", col="#DDDDDD", las=1)
hist(vb_lsu_merged_p3$obs_otus, main="Community richness", xlab="", 
     breaks=10, border="#000000", col="#EEDD88", las=1)


# To test for normality

shapiro.test(vb_lsu_merged_p3$chao1)      
shapiro.test(vb_lsu_merged_p3$shannon)    # normally dist
shapiro.test(vb_lsu_merged_p3$faith_pd)   # normally dist
shapiro.test(vb_lsu_merged_p3$obs_otus)   


# i-chao1
kruskal.test(chao1 ~ insecticide, data=vb_lsu_merged_p3)
kruskal.test(chao1 ~ bagging, data=vb_lsu_merged_p3)

# v-shannon
vb.lsu.shan.ins.aov <- aov(shannon ~ insecticide, data=vb_lsu_merged_p3)
summary(vb.lsu.shan.ins.aov)
vb.lsu.shan.bag.aov <- aov(shannon ~ bagging, data=vb_lsu_merged_p3)
summary(vb.lsu.shan.bag.aov)

# iii-Faith
vb.lsu.fd.ins.aov <- aov(faith_pd ~ insecticide, data=vb_lsu_merged_p3)
summary(vb.lsu.fd.ins.aov)
vb.lsu.fd.bag.aov <- aov(faith_pd ~ bagging, data=vb_lsu_merged_p3)
summary(vb.lsu.fd.bag.aov)

# iv-Observed_otus
kruskal.test(obs_otus ~ insecticide, data=vb_lsu_merged_p3)
kruskal.test(obs_otus ~ bagging, data=vb_lsu_merged_p3)


##################### import cfvb_taxonomy file for 16S ######################

cfvb_lsu_tax_p3 <- read.table("data/cfvb_16s_taxonomy.tsv", 
                           header = T, sep = "\t") %>%
  select(Feature.ID, Taxon) %>% 
  mutate(Taxon = str_replace_all(Taxon, "d__", ""),#replace k_ with ,
         Taxon = str_replace_all(Taxon, "; p__", ","), #replace p_ with ,
         Taxon = str_replace_all(Taxon, "; c__", ","), #replace c_ with ,
         Taxon = str_replace_all(Taxon, "; o__", ","), #replace o_ with ,
         Taxon = str_replace_all(Taxon, "; f__", ","), #replace f_ with ,
         Taxon = str_replace_all(Taxon, "; g__", ","), #replace g_ with ,
         Taxon = str_replace_all(Taxon, "; s__", ","))%>%  #replace s_ with ,
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species"), ",")
#because of the warning message on NAs, we need to address this
#convert any NA values to empty spaces
cfvb_lsu_tax_p3[is.na(cfvb_lsu_tax_p3)] <- ""

#convert anything identified as "unidentified" into an empty space
cfvb_lsu_tax_p3[cfvb_lsu_tax_p3=="Unclassified"] <- "" 

#following code replaces those empty spaces with the Unclassified highest taxon
for (i in 1:nrow(cfvb_lsu_tax_p3)) {
  if (cfvb_lsu_tax_p3[i,7] != ""){
    cfvb_lsu_tax_p3$Species[i] <- paste(cfvb_lsu_tax_p3$Genus[i], sep = " ")
  }else if (cfvb_lsu_tax_p3[i,2] == ""){
    Kingdom <- paste("Unclassified", cfvb_lsu_tax_p3[i,1], sep = " ")
    cfvb_lsu_tax_p3[i, 2:7] <- Kingdom
  }else if (cfvb_lsu_tax_p3[i,3] == ""){
    Phylum <- paste("Unclassified", cfvb_lsu_tax_p3[i,2], sep = " ")
    cfvb_lsu_tax_p3[i, 3:7] <- Phylum
  }else if (cfvb_lsu_tax_p3[i,4] == ""){
    Class <- paste("Unclassified", cfvb_lsu_tax_p3[i,3], sep = " ")
    cfvb_lsu_tax_p3[i, 4:7] <- Class
  }else if (cfvb_lsu_tax_p3[i,5] == ""){
    Order <- paste("Unclassified", cfvb_lsu_tax_p3[i,4], sep = " ")
    cfvb_lsu_tax_p3[i, 5:7] <- Order 
  }else if (cfvb_lsu_tax_p3[i,6] == ""){
    Family <- paste("Unclassified", cfvb_lsu_tax_p3[i,5], sep = " ")
    cfvb_lsu_tax_p3[i, 6:7] <- Family
  }else if (cfvb_lsu_tax_p3[i,7] == ""){
    cfvb_lsu_tax_p3$Species[i] <- paste("Unclassified", 
                                        cfvb_lsu_tax_p3$Genus[i], sep = " ")
  }
}

#convert lsu_merged to long format and join to the cleaned taxonomy table

cfvb_lsu_merged_p3 %>% 
  select(-c(6:9)) %>% 
  pivot_longer(-c(sample_id, bagging, insecticide, 
                  variety, period), names_to = "OTU", 
               values_to = "Counts") %>% 
  inner_join(., cfvb_lsu_tax_p3, by=c('OTU'='Feature.ID')) %>% 
  write.table(., "results/cfvb_lsu_merged_long_p3.tsv", sep="\t", 
              row.names=FALSE, quote = FALSE)

# preparing bacteria phyla plot by variety

cfvb_lsu_merged_long_p3_ngs_cul <- read.table(
  "results/cfvb_lsu_merged_long_p3_ngs_cul_edited.txt",head=T,sep="\t"
) |> rename_all(tolower)

# rel abun for phylum by variety
cfvb_lsu_phylum_rel_abun_p3 <- cfvb_lsu_merged_long_p3_ngs_cul|> 
  group_by(variety) |> 
  mutate(rel_abun=counts/sum(counts)) |> 
  ungroup() |> 
  filter(method != "cul") |> 
  select(-c(otu, counts, method)) |> 
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", 
                        "family", "genus","species", "asv"), 
               names_to = "level", values_to = "taxon")


# bacterial pie chart phyla

cfvb_lsu_phylum_p3 <- cfvb_lsu_phylum_rel_abun_p3 |> 
  filter(level == "phylum") |>
  filter(rel_abun > 0) |>
  group_by(variety,taxon) |> 
  tally() |>
  #drop_na() |> 
  mutate(pct = n/sum(n),
         pct_label = paste0(round(pct * 100), "%")) |> #creates a label column
  ungroup()

cfvb_lsu_phylum_p3_pie <- cfvb_lsu_phylum_p3 %>%
  ggplot(., 
  aes(factor(x="",level=c('Firmicutes','Deinococcota','Pseudomonadota',
                          'Actinobacteriota','Bacteroidota',
                          'Unclassified Bacteria','Actinomycetota',
                          'Chloroflexota','Armatimonadota','Acidobacteriota')), 
                      y=pct, fill=taxon, label=pct_label)) +
  geom_col(color="#000000") + 
  geom_text(color="#000000", position=position_stack(vjust=.5)) +
  coord_polar(theta = 'y', start = 0) +
  facet_wrap(~variety) +
  scale_fill_manual(
    name = NULL,
    breaks = c('Firmicutes','Deinococcota','Pseudomonadota','Actinobacteriota',
               'Bacteroidota','Unclassified Bacteria','Actinomycetota',
               'Chloroflexota','Armatimonadota','Acidobacteriota'),
    labels=c('*Firmicutes*','*Deinococcota*','*Pseudomonadota*',
             '*Actinobacteriota*','*Bacteroidota*','Unclassified *Bacteria*',
             '*Actinomycetota*','*Chloroflexota*','*Armatimonadota*',
             '*Acidobacteriota*'),
    values = c("#708090","#99DDFF","#dddddd","#B452CD","#ffaabb","#eedd88",
               "#6e8b3d","#9ACD32","#ff0000","#0000cd","#77aadd","#696969",
               "#000000","#ee8866","#f4a460")) +
  theme_classic() +
  theme_void() + # this is shortcut to have all themes to element_blank()
  theme(
    panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF"),
    strip.text = element_text(size = 10, family = "sans"),
    strip.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF"),
    legend.title = element_blank(),
    legend.position = "bottom", #c(.5, .13),
    legend.text = element_markdown(size=8, family = "sans"),
    legend.key.size = unit(8, "pt"),
    legend.spacing.y = unit(.4, "cm")
    #legend.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")
  ) + guides(fill = guide_legend(byrow = TRUE))

ggsave("figures/p3_figures/cfvb_lsu_phylum_p3_pie.pdf", 
       cfvb_lsu_phylum_p3_pie, width = 5.5, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_phylum_p3_pie.svg", 
       cfvb_lsu_phylum_p3_pie, width = 5.5, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_phylum_p3_pie.tiff", 
       cfvb_lsu_phylum_p3_pie, width = 5.5, height = 4, dpi = 300)


# a. bacterial genus by insecticide b/w CF and VB stacked bar chats

cfvb_lsu_rel_abun_genus_var_ins_p3 <- cfvb_lsu_merged_long_p3_ngs_cul |> 
  select(-c(bagging,insecticide)) |> 
  group_by(spray, variety) |> 
  mutate(rel_abun=counts/sum(counts)) |> 
  ungroup() |>
  filter(method != "cul") |> 
  select(-c(otu, counts, method)) |> 
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", 
                        "family", "genus","species", "asv"), 
               names_to = "level", values_to = "taxon")

cfvb_lsu_genus_var_ins_p3 <- cfvb_lsu_rel_abun_genus_var_ins_p3 %>% 
  filter(level=="genus") %>% 
  filter(rel_abun > 0) |> 
  group_by(spray, taxon, variety) |> 
  summarise(mean_rel_abun = mean(100*sum(rel_abun)), .groups="drop")

#arrange rel_abun in increasing order
cfvb_lsu_genus_pool_var_ins_p3 <- cfvb_lsu_genus_var_ins_p3 |>  
  group_by(taxon) |>  
  summarise(pool=max(mean_rel_abun) < 6.81, .groups = "drop") #|>  
  #arrange(desc(max))

cfvb_lsu_genus_var_ins_plot_p3 <- inner_join(
  cfvb_lsu_genus_var_ins_p3, cfvb_lsu_genus_pool_var_ins_p3,
  by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |>  
  group_by(spray, taxon) |>  
  ggplot(aes(x=spray, y=mean_rel_abun, fill=taxon)) + 
  geom_col(show.legend = T) +
  facet_wrap(~variety, scales = "free") +
  scale_y_continuous(limits=c(0,100), expand = c(0,0)) +
  scale_fill_manual(
    name=NULL,
    breaks = c("Aeribacillus","Deinococcus","Escherichia","Geobacillus",
               "JC017","Other",'Tetragenococcus',"Thermoanaerobacterium",
               "Thermus","Unclassified Anoxybacillaceae","Unclassified Bacteria"),
    labels = c("*Aeribacillus*","*Deinococcus*","*Escherichia*",
               "*Geobacillus*","*JC017*","Other",'*Tetragenococcus*',
               "*Thermoanaerobacterium*","*Thermus*",
               "Unclassified *Anoxybacillaceae*","Unclassified *Bacteria*"),
    values = c("#708090","#99DDFF","#f3f3f3","#b6d7a8","#6aa84f","#8200ff",
               "#545454","#000000","#bcbcbc","#f0e1d5","#a384b4","#eedd88",
               "#6e8b3d","#9ACD32","#F4A460","#ffffff","#77aadd","#696969",
               "#0000cd","#ffaabb","#B452CD","#000000","#ee8866","#bbbbbb")
  ) +
  theme_bw() +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=9,family="sans",color="#000000"),
        plot.title = element_text(hjust=0.5),
        panel.border = element_blank(),
        axis.text = element_text(size=8,family="sans",color="black"),
        axis.title.y = element_text(size=9,family="sans",color="#000000"),
        legend.text = element_markdown(size=8,family="sans",color="#000000"),
        legend.key.size = unit(8, "pt"),
        legend.spacing.y = unit(.5, "cm"),
        axis.title.x = element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x=NULL, y="Mean Relative Abundance (%)",
       title = " ")
ggsave("figures/p3_figures/cfvb_lsu_genus_var_ins_plot_p3.pdf", 
       cfvb_lsu_genus_var_ins_plot_p3, width = 5, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_genus_var_ins_plot_p3.svg", 
       cfvb_lsu_genus_var_ins_plot_p3, width = 5, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_genus_var_ins_plot_p3.tiff", 
       cfvb_lsu_genus_var_ins_plot_p3, width = 5, height = 4, dpi = 300)

# b. bacterial genus by insecticide b/w CF and VB unstacked bar chats

cfvb_lsu_genus_var_ins_plot_p3b<- inner_join(
  cfvb_lsu_genus_var_ins_p3, cfvb_lsu_genus_pool_var_ins_p3,
  by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |>  
  mutate(taxon=factor(taxon),
         taxon=fct_reorder(taxon, mean_rel_abun, .desc=F)) |> 
  group_by(spray, taxon) |>  
  ggplot(aes(x=taxon,y=mean_rel_abun, fill=spray)) + 
  geom_col(show.legend = T, position=position_dodge()) +
  facet_wrap(~variety, scales = "free") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(breaks=c("+Max","-Max"),
                    labels=c("+Max", "-Max"),
                    values=c("#708090","#0000cd"))+
  scale_x_discrete(
      name=NULL,
      breaks = c("Aeribacillus","Deinococcus","Escherichia","Geobacillus",
                 "JC017","Other",'Tetragenococcus',"Thermoanaerobacterium",
                 "Thermus","Unclassified Anoxybacillaceae","Unclassified Bacteria"),
      labels = c("*Aeribacillus*","*Deinococcus*","*Escherichia*",
                 "*Geobacillus*","*JC017*","Other",'*Tetragenococcus*',
                 "*Thermoanaerobacterium*","*Thermus*",
                 "Unclassified *Anoxybacillaceae*","Unclassified *Bacteria*")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=9,family="sans",color="#000000"),
        plot.title = element_text(hjust=0.5),
        panel.border = element_blank(),
        panel.grid.major.x = element_line(),
        axis.text.y = element_markdown(size=8,family="sans",color="black"),
        axis.title.y = element_text(size=8,family="sans",color="#000000"),
        legend.text = element_text(size=8,family="sans",color="#000000"),
        legend.key.size = unit(9, "pt"),
        legend.spacing.y = unit(.5, "cm"),
        axis.ticks.y = element_line()) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x=NULL, y="Mean Relative Abundance (%)",
       title = " ")
ggsave("figures/p3_figures/cfvb_lsu_genus_var_ins_plot_p3b.pdf", 
       cfvb_lsu_genus_var_ins_plot_p3b, width = 7, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_genus_var_ins_plot_p3b.svg", 
       cfvb_lsu_genus_var_ins_plot_p3b, width = 7, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_genus_var_ins_plot_p3b.tiff", 
       cfvb_lsu_genus_var_ins_plot_p3b, width = 7, height = 4, dpi = 300)


# c. bacterial genus by netting b/w CF and VB stacked bar chats

cfvb_lsu_rel_abun_genus_var_net_p3 <- cfvb_lsu_merged_long_p3_ngs_cul |> 
  select(-c(bagging,insecticide)) |> 
  group_by(bag, variety) |> 
  mutate(rel_abun=counts/sum(counts)) |> 
  ungroup() |>
  filter(method != "cul") |> 
  select(-c(otu, counts, method)) |> 
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", 
                        "family", "genus","species", "asv"), 
               names_to = "level", values_to = "taxon")

cfvb_lsu_genus_var_net_p3 <- cfvb_lsu_rel_abun_genus_var_net_p3 %>% 
  filter(level=="genus") %>% 
  filter(rel_abun > 0) |> 
  group_by(bag, taxon, variety) |> 
  summarise(mean_rel_abun = mean(100*sum(rel_abun)), .groups="drop")

#arrange rel_abun in increasing order
cfvb_lsu_genus_pool_var_net_p3 <- cfvb_lsu_genus_var_net_p3 |>  
  group_by(taxon) |>  
  summarise(pool=max(mean_rel_abun) < 6.10, .groups = "drop") #|>  
  #arrange(desc(max))

cfvb_lsu_genus_var_net_plot_p3 <- inner_join(
  cfvb_lsu_genus_var_net_p3, cfvb_lsu_genus_pool_var_net_p3,
  by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |>  
  group_by(bag, taxon) |>  
  ggplot(aes(x=bag, y=mean_rel_abun, fill=taxon)) + 
  geom_col(show.legend = T) +
  facet_wrap(~variety, scales = "free") +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  scale_fill_manual(
    name=NULL,
    breaks = c("Aeribacillus","Deinococcus","Escherichia","Geobacillus",
               "JC017","Other",'Tetragenococcus',"Thermoanaerobacterium",
               "Thermus","Unclassified Anoxybacillaceae","Unclassified Bacteria"),
    labels = c("*Aeribacillus*","*Deinococcus*","*Escherichia*",
               "*Geobacillus*","*JC017*","Other",'*Tetragenococcus*',
               "*Thermoanaerobacterium*","*Thermus*",
               "Unclassified *Anoxybacillaceae*","Unclassified *Bacteria*"),
    values = c("#708090","#99DDFF","#f3f3f3","#b6d7a8","#6aa84f","#8200ff",
               "#545454","#000000","#bcbcbc","#f0e1d5","#a384b4","#eedd88",
               "#6e8b3d","#9ACD32","#F4A460","#ffffff","#77aadd","#696969",
               "#0000cd","#ffaabb","#B452CD","#000000","#ee8866","#bbbbbb")
  ) +
  theme_bw() +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=9, family="sans", color="#000000"),
        plot.title = element_text(hjust=0.5),
        panel.border = element_blank(),
        axis.text = element_text(size=8,family="sans",color="#000000"),
        axis.title.y = element_text(size=9,family="sans",color="#000000"),
        legend.text = element_markdown(size=8,family="sans",color="#000000"),
        legend.key.size = unit(8, "pt"),
        legend.spacing.y = unit(.5, "cm"),
        axis.title.x = element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x=NULL, y="Mean Relative Abundance (%)",
       title = " ")
ggsave("figures/p3_figures/cfvb_lsu_genus_var_net_plot_p3.pdf", 
       cfvb_lsu_genus_var_net_plot_p3, width = 5, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_genus_var_net_plot_p3.svg", 
       cfvb_lsu_genus_var_net_plot_p3, width = 5, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_genus_var_net_plot_p3.tiff", 
       cfvb_lsu_genus_var_net_plot_p3, width = 5, height = 4, dpi = 300)

# d. bacterial genus by netting b/w CF and VB unstacked bar chats

cfvb_lsu_genus_var_net_plot_p3b <- inner_join(
  cfvb_lsu_genus_var_net_p3, cfvb_lsu_genus_pool_var_net_p3,
  by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |>  
  mutate(taxon=factor(taxon),
         taxon=fct_reorder(taxon, mean_rel_abun, .desc=F)) |>
  group_by(bag, taxon) |>  
  ggplot(aes(x=taxon, y=mean_rel_abun, fill=bag)) + 
  geom_col(show.legend = T, position=position_dodge()) +
  facet_wrap(~variety, scales = "free") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(breaks=c("+Net","-Net"),
                    labels=c("+Net", "-Net"),
                    values=c("#708090","#0000cd"))+
  scale_x_discrete(
    name=NULL,
    breaks = c("Aeribacillus","Deinococcus","Escherichia","Geobacillus",
               "JC017","Other",'Tetragenococcus',"Thermoanaerobacterium",
               "Thermus","Unclassified Anoxybacillaceae","Unclassified Bacteria"),
    labels = c("*Aeribacillus*","*Deinococcus*","*Escherichia*",
               "*Geobacillus*","*JC017*","Other",'*Tetragenococcus*',
               "*Thermoanaerobacterium*","*Thermus*",
               "Unclassified *Anoxybacillaceae*","Unclassified *Bacteria*"))+
  theme_classic() +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=9,family="sans",color="#000000"),
        plot.title = element_text(hjust=0.5),
        panel.border = element_blank(),
        panel.grid.major.x = element_line(),
        axis.text.y = element_markdown(size=8,family="sans",color="black"),
        axis.title.y = element_text(size=8,family="sans",color="#000000"),
        legend.text = element_text(size=8,family="sans",color="#000000"),
        legend.key.size = unit(9, "pt"),
        legend.spacing.y = unit(.5, "cm"),
        axis.ticks = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x=NULL, y="Mean Relative Abundance (%)",
       title = " ")
ggsave("figures/p3_figures/cfvb_lsu_genus_var_net_plot_p3b.pdf", 
       cfvb_lsu_genus_var_net_plot_p3b, width = 7, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_genus_var_net_plot_p3b.svg", 
       cfvb_lsu_genus_var_net_plot_p3b, width = 7, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_genus_var_net_plot_p3b.tiff", 
       cfvb_lsu_genus_var_net_plot_p3b, width = 7, height = 4, dpi = 300)


####################### Venn Diagrams and heat maps ############################

# a. 16S net core biomes data

read_excel("data/cfvb_16s_taxonomy_table_p3_level-7.xlsx",3) |>
  pivot_longer(cols = c(
    "NNet...2","NNet...3","Net...4","Net...5","NNet...6","NNet...7","Net...8",
    "Net...9","NNet...10","NNet...11","Net...12","Net...13","NNet...14",
    "NNet...15","Net...16","Net...17","NNet...18","NNet...19","Net...20",
    "Net...21","NNet...22","NNet...23","Net...24","Net...25"),
    names_to='bagging', values_to='counts') %>%
  write.table(., file = "results/cfvb_16s_tab_4_cor_biom_bag_long_p3.csv",
              col.names = NA, row.names = T, sep = ",", quote = F)


# b. 16S insecticide core biomes data

read_excel("data/cfvb_16s_taxonomy_table_p3_level-7.xlsx", 4) |>
  pivot_longer(cols = c(
    "NMax...2","MMax...3","NMax...4","MMax...5","NMax...6","MMax...7",
    "NMax...8","MMax...9","NMax...10","MMax...11","NMax...12","MMax...13",
    "NMax...14","MMax...15","NMax...16","MMax...17","NMax...18","MMax...19",
    "NMax...20","MMax...21","NMax...22","MMax...23","NMax...24","MMax...25"),
    names_to='insecticide', values_to='counts') %>%
  write.table(., file = "results/cfvb_16s_tab_4_cor_biom_ins_long_p3.csv",
              col.names = NA, row.names = T, sep = ",", quote = F)

# c. 16S variety core biomes data

read_excel("data/cfvb_16s_taxonomy_table_p3_level-7.xlsx", 5) |>
  pivot_longer(cols = c(
    "CF...2","CF...3","CF...4","CF...5","CF...6","CF...7","CF...8","CF...9",
    "CF...10","CF...11","CF...12","CF...13","VB...14","VB...15","VB...16",
    "VB...17","VB...18","VB...19","VB...20","VB...21","VB...22","VB...23",
    "VB...24","VB...25"),
    names_to='variety', values_to='counts') %>%
  write.table(., file = "results/cfvb_16s_tab_4_cor_biom_var_long_p3.csv",
              col.names = NA, row.names = T, sep = ",", quote = F)


# venn diagrams from the combined ngs and cul data without otus

read_excel("results/cfvb_16s_tab_4_cor_biom_var_long_p3_edited.xlsx", 1) |>
  rename_all(tolower) |> 
  filter(counts > 0) |> 
  filter(variable=="CF") |> 
  arrange(genus) |>  
  select(genus,variable) |> 
  distinct(genus,variable) %>%
  write.table(., file="results/cfvb_16s_cf_distinct_genus_p3.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)

read_excel("results/cfvb_16s_tab_4_cor_biom_var_long_p3_edited.xlsx", 1) |>
  rename_all(tolower) |> 
  filter(counts > 0) |> 
  filter(variable=="VB") |> 
  arrange(genus) |>  
  select(genus,variable) |> 
  distinct(genus,variable) %>%
  write.table(., file="results/cfvb_16s_vb_distinct_genus_p3.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)

read_excel("results/cfvb_16s_tab_4_cor_biom_var_long_p3_edited.xlsx", 1) |>
  rename_all(tolower) |> 
  filter(counts > 0) |> 
  filter(variable=="+Max") |> 
  arrange(genus) |>  
  select(genus,variable) |> 
  distinct(genus,variable) %>%
  write.table(., file="results/cfvb_16s_mmax_distinct_genus_p3.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)

read_excel("results/cfvb_16s_tab_4_cor_biom_var_long_p3_edited.xlsx", 1) |>
  rename_all(tolower) |> 
  filter(counts > 0) |> 
  filter(variable=="-Max") |> 
  arrange(genus) |>  
  select(genus,variable) |> 
  distinct(genus,variable) %>%
  write.table(., file="results/cfvb_16s_nmax_distinct_genus_p3.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)

read_excel("results/cfvb_16s_tab_4_cor_biom_var_long_p3_edited.xlsx", 1) |>
  rename_all(tolower) |> 
  filter(counts > 0) |> 
  filter(variable=="+Net") |> 
  arrange(genus) |>  
  select(genus,variable) |> 
  distinct(genus,variable) %>%
  write.table(., file="results/cfvb_16s_net_distinct_genus_p3.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)

read_excel("results/cfvb_16s_tab_4_cor_biom_var_long_p3_edited.xlsx", 1) |>
  rename_all(tolower) |> 
  filter(counts > 0) |> 
  filter(variable=="-Net") |> 
  arrange(genus) |>  
  select(genus,variable) |> 
  distinct(genus,variable) %>%
  write.table(., file="results/cfvb_16s_nnet_distinct_genus_p3.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)


# import edited cfvb_16s_cf_vb_max_nmax_net_nnet_genus_merged data 4 venn diag.

cfvb_16s_cf_vb_max_nmax_net_nnet_species_merged <- read.table(
  "results/cfvb_16s_cf_vb_max_nmax_net_nnet_species_merged.txt", 
  header=T, sep="\t") |> view()

# creating list vector of venn diagram for cf_vb_mmax_nmax
cfvb_16s_cf_vb_max_nmax_species_merged_venn_list = list(
  CF = cfvb_16s_cf_vb_max_nmax_net_nnet_species_merged$species[1:73],
  VB = cfvb_16s_cf_vb_max_nmax_net_nnet_species_merged$species[74:138],
  Mmax = cfvb_16s_cf_vb_max_nmax_net_nnet_species_merged$species[139:214],
  Nmax = cfvb_16s_cf_vb_max_nmax_net_nnet_species_merged$species[215:276]
)
# Ploting the venn diagram for shared mycobiomes between -Max and +Max
cfvb_16s_cf_vb_max_nmax_genus_merged_venn_plot <- ggVennDiagram(
  cfvb_16s_cf_vb_max_nmax_genus_merged_venn_list, label_alpha = 0,
  category.names = c("CF", "VB", "+Max", "-Max")
) +
  scale_fill_gradient(low="#ffffff",high = "#99ddff") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.ticks.y.left = element_blank(),
    legend.position = "None"
  )
ggsave("figures/p3_figures/cfvb_16s_cf_vb_max_nmax_genus_venn_plot.pdf",
       cfvb_16s_cf_vb_max_nmax_genus_merged_venn_plot, 
       height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_16s_cf_vb_max_nmax_genus_venn_plot.tiff",
       cfvb_16s_cf_vb_max_nmax_genus_merged_venn_plot, 
       height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_16s_cf_vb_max_nmax_genus_venn_plot.svg",
       cfvb_16s_cf_vb_max_nmax_genus_merged_venn_plot, 
       height = 4, width = 5, dpi = 300)


# creating list vector of venn diagram for cf_vb_net_nnet
cfvb_16s_cf_vb_net_nnet_species_merged_venn_list = list(
  CF = cfvb_16s_cf_vb_max_nmax_net_nnet_species_merged$species[1:73],
  VB = cfvb_16s_cf_vb_max_nmax_net_nnet_species_merged$species[74:138],
  Net = cfvb_16s_cf_vb_max_nmax_net_nnet_species_merged$species[277:341],
  NNet = cfvb_16s_cf_vb_max_nmax_net_nnet_species_merged$species[342:416]
)
# Ploting the venn diagram for shared mycobiomes between -Net and +Net
cfvb_16s_cf_vb_net_nnet_species_merged_venn_plot <- ggVennDiagram(
  cfvb_16s_cf_vb_net_nnet_species_merged_venn_list, label_alpha = 0,
  category.names = c("CF", "VB", "+Net", "-Net")
) +
  scale_fill_gradient(low="#ffffff",high = "#b452cd") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.ticks.y.left = element_blank(),
    legend.position = "None"
  )
ggsave("figures/p3_figures/cfvb_16s_cf_vb_net_nnet_species_venn_plot.pdf",
       cfvb_16s_cf_vb_net_nnet_species_merged_venn_plot, 
       height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_16s_cf_vb_net_nnet_species_venn_plot.tiff",
       cfvb_16s_cf_vb_net_nnet_species_merged_venn_plot, 
       height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_16s_cf_vb_net_nnet_species_venn_plot.svg",
       cfvb_16s_cf_vb_net_nnet_species_merged_venn_plot, 
       height = 4, width = 5, dpi = 300)


# import edited bs_healthy_infected_species_merged data for venn diagram

cfvb_16s_cf_vb_max_nmax_net_nnet_genus_merged <- read.table(
  "results/cfvb_16s_cf_vb_max_nmax_net_nnet_genus_merged.txt",
  header=T, sep="\t") #|> view()


# creating list vector of venn diagram for cf_vb_mmax_nmax
cfvb_16s_cf_vb_max_nmax_genus_merged_venn_list = list(
  CF = cfvb_16s_cf_vb_max_nmax_net_nnet_genus_merged$genus[1:50],
  VB = cfvb_16s_cf_vb_max_nmax_net_nnet_genus_merged$genus[51:93],
  Mmax = cfvb_16s_cf_vb_max_nmax_net_nnet_genus_merged$genus[94:142],
  Nmax = cfvb_16s_cf_vb_max_nmax_net_nnet_genus_merged$genus[143:186]
)
# Ploting the venn diagram for shared mycobiomes between -Max and +Max
cfvb_16s_cf_vb_max_nmax_genus_merged_venn_plot <- ggVennDiagram(
  cfvb_16s_cf_vb_max_nmax_genus_merged_venn_list, label_alpha = 0,
  category.names = c("CF", "VB", "+Max", "-Max")
) +
  scale_fill_gradient(low="#ffffff",high = "#99ddff") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.ticks.y.left = element_blank(),
    legend.position = "None"
  )
ggsave("figures/p3_figures/cfvb_16s_cf_vb_max_nmax_genus_venn_plot.pdf",
       cfvb_16s_cf_vb_max_nmax_genus_merged_venn_plot,
       height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_16s_cf_vb_max_nmax_genus_venn_plot.tiff",
       cfvb_16s_cf_vb_max_nmax_genus_merged_venn_plot,
       height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_16s_cf_vb_max_nmax_genus_venn_plot.svg",
       cfvb_16s_cf_vb_max_nmax_genus_merged_venn_plot,
       height = 4, width = 5, dpi = 300)


# creating list vector of venn diagram for cf_vb_net_nnet
cfvb_16s_cf_vb_net_nnet_genus_merged_venn_list = list(
  CF = cfvb_16s_cf_vb_max_nmax_net_nnet_genus_merged$genus[1:50],
  VB = cfvb_16s_cf_vb_max_nmax_net_nnet_genus_merged$genus[51:93],
  Net = cfvb_16s_cf_vb_max_nmax_net_nnet_genus_merged$genus[187:231],
  NNet = cfvb_16s_cf_vb_max_nmax_net_nnet_genus_merged$genus[232:282]
)
# Ploting the venn diagram for shared mycobiomes between -Net and +Net
cfvb_16s_cf_vb_net_nnet_genus_merged_venn_plot <- ggVennDiagram(
  cfvb_16s_cf_vb_net_nnet_genus_merged_venn_list, label_alpha = 0,
  category.names = c("CF", "VB", "+Net", "-Net")
) +
  scale_fill_gradient(low="#ffffff",high = "#b452cd") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.ticks.y.left = element_blank(),
    legend.position = "None"
  )
ggsave("figures/p3_figures/cfvb_16s_cf_vb_net_nnet_genuss_venn_plot.pdf",
       cfvb_16s_cf_vb_net_nnet_genus_merged_venn_plot,
       height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_16s_cf_vb_net_nnet_genus_venn_plot.tiff",
       cfvb_16s_cf_vb_net_nnet_genus_merged_venn_plot,
       height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_16s_cf_vb_net_nnet_genus_venn_plot.svg",
       cfvb_16s_cf_vb_net_nnet_genus_merged_venn_plot,
       height = 4, width = 5, dpi = 300)


################################# Heatmap ####################################

cfvb_16s_cf_vb_max_nmax_net_nnet_top50_core_genus <- read_excel(
  "results/cfvb_16s_tab_4_cor_biom_var_long_p3_edited.xlsx") |> 
  select(-`Grand Total`) |>
  select(genus=Genus, everything()) |> 
  mutate_at(c('CF','VB','+Max','-Max','+Net','-Net'), 
            ~replace_na(.,0)) # replaces all NAs with zeros

# melt converts data into long format but then we loose our rownames
cfvb_16s_cf_vb_max_nmax_net_nnet_top50_core_genus_melt <- melt(
  cfvb_16s_cf_vb_max_nmax_net_nnet_top50_core_genus)

# #add column for fungal OTU name
# bs_16s_top50_core_biomes_melt$bact_otu <- rep(
#   row.names(bs_16s_top50_core_biomes), 2)

#view first six rows
head(cfvb_16s_cf_vb_max_nmax_net_nnet_top50_core_genus_melt)

# rescale values for all variables in melted data frame
cfvb_16s_cf_vb_max_nmax_net_nnet_top50_core_genus_melt <- ddply(
  cfvb_16s_cf_vb_max_nmax_net_nnet_top50_core_genus_melt, .(variable), 
  transform, rescale = rescale(value)) # rescale converts variables to a new range

# create heatmap using rescaled values
cfvb_16s_cf_vb_max_nmax_net_nnet_top50_core_genus_hm <- ggplot(
  cfvb_16s_cf_vb_max_nmax_net_nnet_top50_core_genus_melt, 
  aes(x=variable, y=genus)) +
  geom_tile(aes(fill = rescale), colour = "white") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradient(low = "#FFFFFF", high = "#0000cd", "#9400D3", "#CD950C",
                      expand = c(0,0)) +
  labs(x = "", y = "") +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=8, color="#000000"),
        axis.text.x = element_markdown(angle=0,hjust=.5,vjust=1,color="black")) +
  coord_fixed(ratio = 0.1) # resizes the boxes for each variable on the x-axis

ggsave("figures/p3_figures/cfvb_16s_cf_vb_ins_bag_top50_core_genus_hm.pdf", 
       cfvb_16s_cf_vb_max_nmax_net_nnet_top50_core_genus_hm, 
       height = 7, width = 8, dpi = 400)
ggsave("figures/p3_figures/cfvb_16s_cf_vb_ins_bag_top50_core_genus_hm.tiff", 
       cfvb_16s_cf_vb_max_nmax_net_nnet_top50_core_genus_hm, 
       height = 7, width = 8, dpi = 400)
ggsave("figures/p3_figures/cfvb_16s_cf_vb_ins_bag_top50_core_genus_hm.svg", 
       cfvb_16s_cf_vb_max_nmax_net_nnet_top50_core_genus_hm, 
       height = 7, width = 8, dpi = 400)

############################### Ordinations ####################################

# NMDS from cfvb_lsu_obs_otus table

# cfvb_lsu_otus_p3 <- cfvb_lsu_otus |>
#   mutate(sample_id2=sample_id) |>
#   select(sample_id,sample_id2, everything()) |> 
#   separate(sample_id2,
#            into = c("cul_treat", "perep"),
#            sep = "(?<=[A-Za-z])(?=[0-9])") |>
#   separate(perep,into=c("period","replicate"), sep="(?<=\\d)") |>
#   select(-c(cul_treat,replicate)) |>
#   filter(period==3) |>
#   select(-period) |> #this is a tibble and will not take row names
#   as.data.frame()
#   
# rownames(cfvb_lsu_otus_p3) <- cfvb_lsu_otus_p3$sample_id
# cfvb_lsu_otus_p3 <- cfvb_lsu_otus_p3[,-1]
# 
# cfvb_lsu_otus_p3 <- as.matrix(cfvb_lsu_otus_p3) #converts to matrix
# 
# #calculating distances in vegan
# 
# set.seed(04291982)
# 
# cfvb_lsu_dist_p3 <- vegdist(cfvb_lsu_otus_p3, method = "bray")
# 
# cfvb_lsu_nmds_p3 <- metaMDS(cfvb_lsu_dist_p3) |> 
#   scores() |> 
#   as_tibble(rownames="sample_id")
# 
# cfvb_lsu_nmds_meta_merged_p3 <- inner_join(cfvb_lsu_nmds_p3, cfvb_lsu_meta_p3,
#                                         by="sample_id") |> 
#   select(-period)
# 
# #plotting  
# cfvb_lsu_nmds_plot_p3 <- cfvb_lsu_nmds_meta_merged_p3 |> 
#   ggplot(aes(x=NMDS1, y=NMDS2, color=variety, shape=variety)) + 
#   geom_point(size=1.25) +
#   scale_color_manual(breaks=c("CF", "VB"),
#                      labels=c("CF", "VB"),
#                      values=c("#ff0000", "#0000cd"))+
#   theme_bw() +
#   theme(legend.background = element_blank(),
#         legend.text= element_text(size = 10, family = 'sans'),
#         legend.title = element_blank(),
#         panel.border = element_rect(color="#dddddd"),
#         legend.key.size = unit(.5, units = "cm"),
#         axis.text = element_text(colour = '#000000',size=9),
#         axis.title = element_text(color="#000000",size=10))
# 
# ggsave("figures/p3_figures/cfvb_lsu_nmds_plot_p3.tiff", 
#        cfvb_lsu_nmds_plot_p3, height = 4, width = 5, dpi = 300)



# OPTION 2

# a. get our metadata
cfvb_lsu_nmds_meta_p3 <- cfvb_lsu_merged_p3 |>
  select(sample_id,bag=bagging,spray=insecticide,variety) |> 
  mutate(bag=case_when(bag=="Net"~"+Net",
                       bag=="Nnet"~"-Net"),
         spray=case_when(spray=="Mmax"~"+Max",
                         spray=="Nmax"~"-Max"))

# b. get our otus
cfvb_lsu_nmds_otu_p3 <- cfvb_lsu_merged_p3 |>
  select(-c(2:9)) |>
  pivot_longer(-sample_id) |> 
  group_by(sample_id) |> 
  # summarise(N=sum(value)) |>  #this gives the total number of seqs in each
  # arrange(N) |> print(n=25) #sample and we arrange to get a value to rarefy
  mutate(N=sum(value)) |> #this shows the rarefied value
  filter(N >= 950) |> #rarefying to keep all samples having 1400 seqs
  select(-N) |> 
  pivot_wider(names_from="name", values_from="value", values_fill=0)

# c. get our distance matrix
cfvb_lsu_nmds_dist_p3 <- cfvb_lsu_nmds_otu_p3 |> 
  column_to_rownames("sample_id") |> 
  avgdist(sample=950) #calculates the bray curtis distances

# .d make ordination
set.seed(3) #required for the nmds to converge
cfvb_lsu_nmds_tbl_p3 <- metaMDS(cfvb_lsu_nmds_dist_p3) |> 
  scores() |> 
  as_tibble(rownames="sample_id") #converts to a tibble that we can use in ggplot

# e. merge nmds_meta and nmds_tbl
cfvb_lsu_nmds_4_plot_p3 <- inner_join(cfvb_lsu_nmds_meta_p3, 
                                      cfvb_lsu_nmds_tbl_p3, by="sample_id")

# test of significance before plotting
# ADONIS
# bagging
adonis2(cfvb_lsu_nmds_dist_p3~cfvb_lsu_nmds_meta_p3$bag,permutations=9999)
# insecticide
adonis2(cfvb_lsu_nmds_dist_p3~cfvb_lsu_nmds_meta_p3$spray,permutations=9999)

#interaction of bagging and insecticide
adonis2(cfvb_lsu_nmds_dist_p3~cfvb_lsu_nmds_meta_p3$spray*
          cfvb_lsu_nmds_meta_p3$bag,permutations=9999)
# variety
adonis2(cfvb_lsu_nmds_dist_p3~cfvb_lsu_nmds_meta_p3$variety,permutations=9999)
# get p-value

# str(adon_var) #check where the p-value is
# adon_var$`Pr(>F)`[1] #since it's a list, we put [1] to extract the first item

# Another statistical test of dispersion is using betadisper() in vegan

# betadisp_var <- anova(betadisper(cfvb_its_nmds_dist_p3,
#                                  cfvb_its_nmds_meta_p3$variety))
# 
# str(betadisp_var)
# betadisp_var$`Pr(>F)`[1]
# 
# #OR
# betadisp_var_permutest <- permutest(betadisper(cfvb_its_nmds_dist_p3,
#                                                cfvb_its_nmds_meta_p3$variety))
# str(betadisp_var_permutest)
# betadisp_var_permutest$tab$`Pr(>F)`[1]

# we prepare a centroid for our ellipses
cfvb_lsu_centroid_p3<- cfvb_lsu_nmds_4_plot_p3 |> 
  group_by(variety) |> 
  summarise(NMDS1=mean(NMDS1),
            NMDS2=mean(NMDS2))

# plotting
cfvb_lsu_nmds_plot_p3 <- cfvb_lsu_nmds_4_plot_p3|>
  ggplot(aes(x=NMDS1, y=NMDS2, color=variety)) +
  geom_point() +
  stat_ellipse(show.legend=F) +
  geom_point(data=cfvb_lsu_centroid_p3, size=3,shape=21,
             color="#000000",aes(fill=variety),show.legend=F) +
  scale_color_manual(values=c("#0000cd","#ff0000")) +
  scale_fill_manual(values=c("#0000cd","#ff0000")) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.background = element_blank(),
        legend.text= element_text(size = 10, family = 'sans'),
        legend.title = element_blank(),
        legend.key.size = unit(.5, units = "cm"),
        axis.text = element_text(colour = 'black'))

ggsave("figures/p3_figures/cfvb_lsu_nmds_plot_p3.pdf", 
       cfvb_lsu_nmds_plot_p3, height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_nmds_plot_p3.svg", 
       cfvb_lsu_nmds_plot_p3, height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_lsu_nmds_plot_p3.tiff", 
       cfvb_lsu_nmds_plot_p3, height = 4, width = 5, dpi = 300)


#################### ITS Data for alpha diversity analysis ##################

# A. ITS data general

#metadata

cfvb_its_meta_p3 <- read_tsv("data/cfvb_its_gsr_metadata.tsv") |> 
  select('sample_id'='sample-id','bagging'='treat_one', 
         'insecticide'='treat_two', variety,period) |> 
  filter(sample_id != "fctrl") |>  #drop the control sample
  filter(period == "3") |> 
  mutate(variety = case_when(variety == "carbfranc" ~ "CF",
                             variety == "vidalblanc" ~ "VB"),
         bagging = case_when(bagging == "NNet" ~ "Nnet",
                             bagging == "Net" ~ "Net"),
         insecticide = case_when(insecticide == "MMax" ~ "Mmax",
                                 insecticide == "NMax" ~ "Nmax"))

cfvb_its_meta_p3$variety <- as.factor(cfvb_its_meta_p3$variety)
cfvb_its_meta_p3$bagging <- as.factor(cfvb_its_meta_p3$bagging)
cfvb_its_meta_p3$insecticide <- as.factor(cfvb_its_meta_p3$insecticide)

glimpse(cfvb_its_meta_p3)

#evenness
cfvb_its_even <- read_tsv("data/cfvb_its_evenness.tsv") |>  
  rename("sample_id" = "...1", 'evenness'='pielou_evenness') #renames col 1

#faith
cfvb_its_faith <- read_tsv("data/cfvb_its_faith_pd.tsv") |>  
  rename("sample_id" = "#SampleID") #renames col 1

#observed_otus
cfvb_its_obs_otus <- read_tsv("data/cfvb_its_observed_otus.tsv") |>  
  rename("sample_id" = "...1", "obs_otus"="observed_features") #rename columns

#shannon
cfvb_its_shan <- read_tsv("data/cfvb_its_shannon.tsv") |>  
  rename("sample_id" = "...1", "shannon"="shannon_entropy")

#otu
cfvb_its_otus <- read_tsv("data/cfvb_its_rarefied_table_transp.txt", skip = 1) |> 
  rename("sample_id"="#OTU ID")

#ace
cfvb_its_ace <- read_tsv("data/cfvb_its_ace.tsv") |>  
  rename("sample_id" = "...1") #rename columns

#chao1
cfvb_its_chao1 <- read_tsv("data/cfvb_its_chao1.tsv") |>  
  rename("sample_id" = "...1")

#merge metadata, ace, chao1, evenness, faith, observed_otus, shannon and otus

cfvb_its_merged_p3 <- inner_join(cfvb_its_meta_p3, 
                                 cfvb_its_chao1, by='sample_id') %>%
  inner_join(., cfvb_its_faith, by='sample_id') %>%
  inner_join(., cfvb_its_obs_otus, by='sample_id') %>%
  inner_join(., cfvb_its_shan, by='sample_id') %>%
  inner_join(., cfvb_its_otus, by='sample_id')


######################## Overall ITs Analysis ##################################
#Create 2x2 plot environment 
par(mfrow = c(2, 2))
hist(cfvb_its_merged_p3$chao1, main="Chao1 diversity", xlab="", breaks=10, 
     border="#000000", col="#99DDFF", las=1)
# hist(cfvb_lsu_merged_p3$ace, main="Ace diversity", xlab="", breaks=10, 
#      border="#000000", col="#BBCC33", las=1)
hist(cfvb_its_merged_p3$shannon, main="Shannon diversity", xlab="", breaks=10, 
     border="#000000", col="#FFAABB", las=1)
# hist(cfvb_lsu_merged_p3$evenness, main="Pielou's evenness", xlab="", breaks=10, 
#      border="#000000", col="#44BB99", las=1)
hist(cfvb_its_merged_p3$faith_pd, main="Phylegenetic diversity", xlab="", 
     breaks=10, border="#000000", col="#DDDDDD", las=1)
hist(cfvb_its_merged_p3$obs_otus, main="Community richness", xlab="", 
     breaks=10, border="#000000", col="#EEDD88", las=1)


# To test for normality

shapiro.test(cfvb_its_merged_p3$chao1)      # normally dist
shapiro.test(cfvb_its_merged_p3$shannon)    # normally dist
shapiro.test(cfvb_its_merged_p3$faith_pd)   # normally dist
shapiro.test(cfvb_its_merged_p3$obs_otus)   # normally dist


# a. Anova test for normally distributed data

# i-chao1
cfvb.its.chao1.ins.aov <- aov(chao1 ~ insecticide, data=cfvb_its_merged_p3)
summary(cfvb.its.chao1.ins.aov)
cfvb.its.chao1.bag.aov <- aov(chao1 ~ bagging, data=cfvb_its_merged_p3)
summary(cfvb.its.chao1.bag.aov)
cfvb.its.chao1.var.aov <- aov(chao1 ~ variety, data=cfvb_its_merged_p3)
summary(cfvb.its.chao1.var.aov)

# v-shannon
cfvb.its.shan.ins.aov <- aov(shannon ~ insecticide, data=cfvb_its_merged_p3)
summary(cfvb.its.shan.ins.aov)
cfvb.its.shan.bag.aov <- aov(shannon ~ bagging, data=cfvb_its_merged_p3)
summary(cfvb.its.shan.bag.aov)
cfvb.its.shan.var.aov <- aov(shannon ~ variety, data=cfvb_its_merged_p3)
summary(cfvb.its.shan.var.aov)

# iii-Faith
cfvb.its.fd.ins.aov <- aov(faith_pd ~ insecticide, data=cfvb_its_merged_p3)
summary(cfvb.its.fd.ins.aov)
cfvb.its.fd.bag.aov <- aov(faith_pd ~ bagging, data=cfvb_its_merged_p3)
summary(cfvb.its.fd.bag.aov)
cfvb.its.fd.var.aov <- aov(faith_pd ~ variety, data=cfvb_its_merged_p3)
summary(cfvb.its.fd.var.aov)

# iv-Observed_otus
cfvb.its.obs.ins.aov <- aov(obs_otus ~ insecticide, data=cfvb_its_merged_p3)
summary(cfvb.its.obs.ins.aov)
cfvb.its.obs.bag.aov <- aov(obs_otus ~ bagging, data=cfvb_its_merged_p3)
summary(cfvb.its.obs.bag.aov)
cfvb.its.obs.var.aov <- aov(obs_otus ~ variety, data=cfvb_its_merged_p3)
summary(cfvb.its.obs.var.aov)


# Plotting relationships ITS General

# boxplots

my_labs <- c("shannon" = "shannon",
             "chao1" = "chao1",
             "faith_pd" = "faith",
             "obs_otus" = "observed otus") # labels for y-axis of boxplots

# a. shannon, chao1, otus and faith by insecticide
cfvb_its_alpha_div_ins_plot_p3 <- cfvb_its_merged_p3 |> 
  select(-period) |> 
  pivot_longer(cols = 5:8, names_to = 'indice', values_to = 'value') |>
  mutate(bagging=case_when(bagging=="Net"~"+Net",
                           bagging=="Nnet"~"-Net"),
         insecticide=case_when(insecticide=="Nmax"~"-Max",
                               insecticide=="Mmax"~"+Max")) |>
  mutate(indice=factor(indice, levels = c("shannon", "chao1", 
                                          "obs_otus", "faith_pd"))) |>  
  select(insecticide, indice, value) |>  
  group_by(insecticide, indice) %>% 
  ggplot(., aes(x=insecticide, y=value)) + 
  geom_boxplot(aes(fill=insecticide),size=.5, outlier.shape = NA) +
  scale_fill_manual(values = c("#FFFFFF","#BFBFBF")) +
  facet_wrap(~indice, scales = "free",strip.position = "left",
             labeller = labeller(indice = my_labs)) +
  stat_compare_means(label = "p.signif", size=4) +
  labs(x = NULL, y = NULL, #removes the x- and y-axis labels
       title = " "
  ) +
  theme_bw() + theme(legend.position = "None") +
  theme(
    plot.title = element_text(hjust=0.5),
    strip.placement = "outside", #sends the panel labels from left to outside 
    #of y-axis
    strip.background = element_blank(), #removes background of strip
    strip.text = element_text(size=9,colour="#000000",family="sans"),
    panel.background = element_blank(),
    axis.title.y = element_text(size=9,family="sans",color="#000000"),
    axis.text = element_text(size=8,family="sans",color="#000000"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(hjust=.5,vjust=1,family="sans",color="#000000"))

ggsave("figures/p3_figures/cfvb_its_alpha_div_ins_plot_p3_new.pdf", 
       cfvb_its_alpha_div_ins_plot_p3, height = 4.5, width = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_alpha_div_ins_plot_p3_new.tif", 
       cfvb_its_alpha_div_ins_plot_p3, height = 4.5, width = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_alpha_div_ins_plot_p3_new.svg", 
       cfvb_its_alpha_div_ins_plot_p3, height = 4.5, width = 4, dpi = 300)

# b. shannon, chao1, otus and faith by bagging
cfvb_its_alpha_div_bag_plot_p3 <- cfvb_its_merged_p3 |> 
  select(-period) |> 
  pivot_longer(cols = 5:8, names_to = 'indice', values_to = 'value') |>
  mutate(bagging=case_when(bagging=="Net"~"+Net",
                           bagging=="Nnet"~"-Net"),
         insecticide=case_when(insecticide=="Nmax"~"-Max",
                               insecticide=="Mmax"~"+Max")) |>
  mutate(indice=factor(indice, levels = c("shannon", "chao1", 
                                          "obs_otus", "faith_pd"))) |>  
  select(bagging, indice, value) |>  
  group_by(bagging, indice) %>% 
  ggplot(., aes(x=bagging, y=value)) + 
  geom_boxplot(aes(fill=bagging),size=.5, outlier.shape = NA) +
  scale_fill_manual(values = c("#FFFFFF","#BFBFBF")) +
  facet_wrap(~indice, scales = "free",strip.position = "left",
             labeller = labeller(indice = my_labs)) +
  stat_compare_means(label = "p.signif", size=4) +
  labs(x = NULL, y = NULL, #removes the x- and y-axis labels
       title = " "
  ) +
  theme_bw() + theme(legend.position = "None") +
  theme(
    plot.title = element_text(hjust=0.5),
    strip.placement = "outside", #sends the panel labels from left to outside 
    #of y-axis
    strip.background = element_blank(), #removes background of strip
    strip.text = element_text(size=9,colour="#000000",family="sans"),
    panel.background = element_blank(),
    axis.title.y = element_text(size=9,family="sans",color="#000000"),
    axis.text = element_text(size=8,family="sans",color="#000000"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(hjust=.5,vjust=1,family="sans",color="#000000"))

ggsave("figures/p3_figures/cfvb_its_alpha_div_bag_plot_p3_new.pdf", 
       cfvb_its_alpha_div_bag_plot_p3, height = 4.5, width = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_alpha_div_bag_plot_p3_new.tif", 
       cfvb_its_alpha_div_bag_plot_p3, height = 4.5, width = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_alpha_div_bag_plot_p3_new.svg", 
       cfvb_its_alpha_div_bag_plot_p3, height = 4.5, width = 4, dpi = 300)


# c. shannon, chao1, otus and faith by variety
cfvb_its_alpha_div_var_plot_p3 <- cfvb_its_merged_p3 |> 
  select(-period) |> 
  pivot_longer(cols = 5:8, names_to = 'indice', values_to = 'value') |>
  mutate(indice=factor(indice, levels = c("shannon", "chao1", 
                                          "obs_otus", "faith_pd"))) |>  
  select(variety, indice, value) |>  
  group_by(variety, indice) %>% 
  ggplot(., aes(x=variety, y=value)) + 
  geom_boxplot(aes(fill=variety),size=.5, outlier.shape = NA) +
  scale_fill_manual(values = c("#FFFFFF","#BFBFBF")) +
  facet_wrap(~indice, scales = "free",strip.position = "left",
             labeller = labeller(indice = my_labs)) +
  stat_compare_means(label = "p.signif", size=4) +
  labs(x = NULL, y = NULL, #removes the x- and y-axis labels
       title = " "
  ) +
  theme_bw() + theme(legend.position = "None") +
  theme(
    plot.title = element_text(hjust=0.5),
    strip.placement = "outside", #sends the panel labels from left to outside 
    #of y-axis
    strip.background = element_blank(), #removes background of strip
    strip.text = element_text(size=9,colour="#000000",family="sans"),
    panel.background = element_blank(),
    axis.title.y = element_text(size=9,family="sans",color="#000000"),
    axis.text = element_text(size=8,family="sans",color="#000000"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(hjust=.5,vjust=1,family="sans",color="#000000"))

ggsave("figures/p3_figures/cfvb_its_alpha_div_var_plot_p3_new.pdf", 
       cfvb_its_alpha_div_var_plot_p3, height = 4.5, width = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_alpha_div_var_plot_p3_new.tif", 
       cfvb_its_alpha_div_var_plot_p3, height = 4.5, width = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_alpha_div_var_plot_p3_new.svg", 
       cfvb_its_alpha_div_var_plot_p3, height = 4.5, width = 4, dpi = 300)



########################  ITs Analysis for CF #################################

cf_its_merged_p3 <- cfvb_its_merged_p3  |> 
  filter(variety=="CF")

#Create 2x2 plot environment 
par(mfrow = c(2, 2))
hist(cf_its_merged_p3$chao1, main="Chao1 diversity", xlab="", breaks=10, 
     border="#000000", col="#99DDFF", las=1)
# hist(cfvb_lsu_merged_p3$ace, main="Ace diversity", xlab="", breaks=10, 
#      border="#000000", col="#BBCC33", las=1)
hist(cf_its_merged_p3$shannon, main="Shannon diversity", xlab="", breaks=10, 
     border="#000000", col="#FFAABB", las=1)
# hist(cfvb_lsu_merged_p3$evenness, main="Pielou's evenness", xlab="", breaks=10, 
#      border="#000000", col="#44BB99", las=1)
hist(cf_its_merged_p3$faith_pd, main="Phylegenetic diversity", xlab="", 
     breaks=10, border="#000000", col="#DDDDDD", las=1)
hist(cf_its_merged_p3$obs_otus, main="Community richness", xlab="", 
     breaks=10, border="#000000", col="#EEDD88", las=1)


# To test for normality

shapiro.test(cf_its_merged_p3$chao1)      # normally dist
shapiro.test(cf_its_merged_p3$shannon)    # normally dist
shapiro.test(cf_its_merged_p3$faith_pd)   # normally dist
shapiro.test(cf_its_merged_p3$obs_otus)   # normally dist


# i-chao1
cf.its.chao1.ins.aov <- aov(chao1 ~ insecticide, data=cf_its_merged_p3)
summary(cf.its.chao1.ins.aov)
cf.its.chao1.bag.aov <- aov(chao1 ~ bagging, data=cf_its_merged_p3)
summary(cf.its.chao1.bag.aov)

# v-shannon
cf.its.shan.ins.aov <- aov(shannon ~ insecticide, data=cf_its_merged_p3)
summary(cf.its.shan.ins.aov)
cf.its.shan.bag.aov <- aov(shannon ~ bagging, data=cf_its_merged_p3)
summary(cf.its.shan.bag.aov)

# iii-Faith
cf.its.fd.ins.aov <- aov(faith_pd ~ insecticide, data=cf_its_merged_p3)
summary(cf.its.fd.ins.aov)
cf.its.fd.bag.aov <- aov(faith_pd ~ bagging, data=cf_its_merged_p3)
summary(cf.its.fd.bag.aov)

# iv-Observed_otus
cf.its.obs.ins.aov <- aov(obs_otus ~ insecticide, data=cf_its_merged_p3)
summary(cf.its.obs.ins.aov)
cf.its.obs.bag.aov <- aov(obs_otus ~ bagging, data=cf_its_merged_p3)
summary(cf.its.obs.bag.aov)


########################  ITs Analysis for VB #################################

vb_its_merged_p3 <- cfvb_its_merged_p3  |> 
  filter(variety=="VB")

#Create 2x2 plot environment 
par(mfrow = c(2, 2))
hist(vb_its_merged_p3$chao1, main="Chao1 diversity", xlab="", breaks=10, 
     border="#000000", col="#99DDFF", las=1)
# hist(cfvb_lsu_merged_p3$ace, main="Ace diversity", xlab="", breaks=10, 
#      border="#000000", col="#BBCC33", las=1)
hist(vb_its_merged_p3$shannon, main="Shannon diversity", xlab="", breaks=10, 
     border="#000000", col="#FFAABB", las=1)
# hist(cfvb_lsu_merged_p3$evenness, main="Pielou's evenness", xlab="", breaks=10, 
#      border="#000000", col="#44BB99", las=1)
hist(vb_its_merged_p3$faith_pd, main="Phylegenetic diversity", xlab="", 
     breaks=10, border="#000000", col="#DDDDDD", las=1)
hist(vb_its_merged_p3$obs_otus, main="Community richness", xlab="", 
     breaks=10, border="#000000", col="#EEDD88", las=1)


# To test for normality

shapiro.test(vb_its_merged_p3$chao1)      # normally dist
shapiro.test(vb_its_merged_p3$shannon)    # normally dist
shapiro.test(vb_its_merged_p3$faith_pd)   # normally dist
shapiro.test(vb_its_merged_p3$obs_otus)   # normally dist


# i-chao1
vb.its.chao1.ins.aov <- aov(chao1 ~ insecticide, data=vb_its_merged_p3)
summary(vb.its.chao1.ins.aov)
vb.its.chao1.bag.aov <- aov(chao1 ~ bagging, data=vb_its_merged_p3)
summary(vb.its.chao1.bag.aov)

# v-shannon
vb.its.shan.ins.aov <- aov(shannon ~ insecticide, data=vb_its_merged_p3)
summary(vb.its.shan.ins.aov)
vb.its.shan.bag.aov <- aov(shannon ~ bagging, data=vb_its_merged_p3)
summary(vb.its.shan.bag.aov)

# iii-Faith
vb.its.fd.ins.aov <- aov(faith_pd ~ insecticide, data=vb_its_merged_p3)
summary(vb.its.fd.ins.aov)
vb.its.fd.bag.aov <- aov(faith_pd ~ bagging, data=vb_its_merged_p3)
summary(vb.its.fd.bag.aov)

# iv-Observed_otus
vb.its.obs.ins.aov <- aov(obs_otus ~ insecticide, data=vb_its_merged_p3)
summary(vb.its.obs.ins.aov)
vb.its.obs.bag.aov <- aov(obs_otus ~ bagging, data=vb_its_merged_p3)
summary(vb.its.obs.bag.aov)



##################### import cfvb_taxonomy file for ITS ######################

cfvb_its_tax_p3 <- read.table("data/cfvb_its_taxonomy.tsv", 
                              header = T, sep = "\t") %>%
  select(Feature.ID, Taxon) %>% 
  mutate(Taxon = str_replace_all(Taxon, "k__", ""),#replace k_ with ,
         Taxon = str_replace_all(Taxon, ";p__", ","), #replace p_ with ,
         Taxon = str_replace_all(Taxon, ";c__", ","), #replace c_ with ,
         Taxon = str_replace_all(Taxon, ";o__", ","), #replace o_ with ,
         Taxon = str_replace_all(Taxon, ";f__", ","), #replace f_ with ,
         Taxon = str_replace_all(Taxon, ";g__", ","), #replace g_ with ,
         Taxon = str_replace_all(Taxon, ";s__", ","))%>%  #replace s_ with ,
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species"), ",")
#because of the warning message on NAs, we need to address this
#convert any NA values to empty spaces
cfvb_its_tax_p3[is.na(cfvb_its_tax_p3)] <- ""

#convert anything identified as "unidentified" into an empty space
cfvb_its_tax_p3[cfvb_its_tax_p3=="Unclassified"] <- "" 

#following code replaces those empty spaces with the Unclassified highest taxon
for (i in 1:nrow(cfvb_its_tax_p3)) {
  if (cfvb_its_tax_p3[i,7] != ""){
    cfvb_its_tax_p3$Species[i] <- paste(cfvb_its_tax_p3$Genus[i], sep = " ")
  }else if (cfvb_its_tax_p3[i,2] == ""){
    Kingdom <- paste("Unclassified", cfvb_its_tax_p3[i,1], sep = " ")
    cfvb_its_tax_p3[i, 2:7] <- Kingdom
  }else if (cfvb_its_tax_p3[i,3] == ""){
    Phylum <- paste("Unclassified", cfvb_its_tax_p3[i,2], sep = " ")
    cfvb_its_tax_p3[i, 3:7] <- Phylum
  }else if (cfvb_its_tax_p3[i,4] == ""){
    Class <- paste("Unclassified", cfvb_its_tax_p3[i,3], sep = " ")
    cfvb_its_tax_p3[i, 4:7] <- Class
  }else if (cfvb_its_tax_p3[i,5] == ""){
    Order <- paste("Unclassified", cfvb_its_tax_p3[i,4], sep = " ")
    cfvb_its_tax_p3[i, 5:7] <- Order 
  }else if (cfvb_its_tax_p3[i,6] == ""){
    Family <- paste("Unclassified", cfvb_its_tax_p3[i,5], sep = " ")
    cfvb_its_tax_p3[i, 6:7] <- Family
  }else if (cfvb_its_tax_p3[i,7] == ""){
    cfvb_its_tax_p3$Species[i] <- paste("Unclassified", 
                                        cfvb_its_tax_p3$Genus[i], sep = " ")
  }
}

#convert its_merged to long format and join to the cleaned taxonomy table

cfvb_its_merged_p3 %>% 
  select(-c(6:9)) %>% 
  pivot_longer(-c(sample_id, bagging, insecticide, 
                  variety, period), names_to = "OTU", 
               values_to = "Counts") %>% 
  inner_join(., cfvb_its_tax_p3, by=c('OTU'='Feature.ID')) %>% 
  write.table(., "results/cfvb_its_merged_long_p3_new.tsv", sep="\t", 
              row.names=FALSE, quote = FALSE)

# preparing bacteria phyla plot by variety

cfvb_its_merged_long_p3_ngs_cul <- read.table(
  "results/cfvb_its_merged_long_p3_ngs_cul_edited_new.txt",head=T,sep="\t"
) |> 
  mutate(bag = case_when(bagging == "Nnet" ~ "-Net",
                         bagging == "Net" ~ "+Net"),
         spray = case_when(insecticide == "Mmax" ~ "+Max",
                           insecticide == "Nmax" ~ "-Max")) |> 
  rename_all(tolower)

# rel abun for phylum by variety
cfvb_its_phylum_rel_abun_p3 <- cfvb_its_merged_long_p3_ngs_cul|> 
  group_by(variety) |> 
  mutate(rel_abun=counts/sum(counts)) |> 
  ungroup() |> 
  filter(method != "cul") |> 
  select(-c(otu, counts, method)) |> 
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", 
                        "family", "genus","species", "asv"), 
               names_to = "level", values_to = "taxon") 


# fungal pie chart phyla

cfvb_its_phylum_p3 <- cfvb_its_phylum_rel_abun_p3 |> 
  filter(level == "phylum") |>
  filter(rel_abun > 0) |>
  group_by(variety,taxon) |> 
  tally() |>
  #drop_na() |> 
  mutate(pct = n/sum(n),
         pct_label = paste0(round(pct * 100), "%")) |> #creates a label column
  ungroup()

cfvb_its_phylum_p3_pie <- cfvb_its_phylum_p3 %>%
  ggplot(., 
         aes(factor(x="",level=c('Ascomycota','Basidiomycota',
                                 'Unclassified Fungi')), 
             y=pct, fill=taxon, label=pct_label)) +
  geom_col(color="#000000") + 
  geom_text(color="#000000", position=position_stack(vjust=.5)) +
  coord_polar(theta = 'y', start = 0) +
  facet_wrap(~variety) +
  scale_fill_manual(
    name = NULL,
    breaks = c('Ascomycota','Basidiomycota','Unclassified Fungi'),
    labels=c('*Ascomycota*','*Basidiomycota*','Unclassified *Fungi*'),
    values = c("#708090","#99DDFF","#d3d3d3")) +
  theme_classic() +
  theme_void() + # this is shortcut to have all themes to element_blank()
  theme(
    panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF"),
    strip.text = element_text(size = 10, family = "sans"),
    strip.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF"),
    legend.title = element_blank(),
    legend.position = "bottom", #c(.5, .13),
    legend.text = element_markdown(size=8, family = "sans"),
    legend.key.size = unit(8, "pt"),
    legend.spacing.y = unit(.4, "cm")
    #legend.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")
  ) + guides(fill = guide_legend(byrow = TRUE))

ggsave("figures/p3_figures/cfvb_its_phylum_p3_pie_new.pdf", 
       cfvb_its_phylum_p3_pie, width = 5.5, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_phylum_p3_pie_new.svg", 
       cfvb_its_phylum_p3_pie, width = 5.5, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_phylum_p3_pie_new.tiff", 
       cfvb_its_phylum_p3_pie, width = 5.5, height = 4, dpi = 300)


# a. fungal genus by insecticide b/w CF and VB stacked bar chats

cfvb_its_rel_abun_genus_var_ins_p3 <- cfvb_its_merged_long_p3_ngs_cul |> 
  select(-c(bagging,insecticide,period)) |> 
  select(sample_id,variety,spray,bag,everything()) |> 
  group_by(spray, variety) |> 
  mutate(rel_abun=counts/sum(counts)) |> 
  ungroup() |> 
  filter(method != "cul") |> 
  select(-c(otu, counts, method)) |> 
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", 
                        "family", "genus","species", "asv"), 
               names_to = "level", values_to = "taxon")

cfvb_its_genus_var_ins_p3 <- cfvb_its_rel_abun_genus_var_ins_p3 %>% 
  filter(level=="genus") %>% 
  filter(rel_abun > 0) |> 
  group_by(spray, taxon, variety) |> 
  summarise(mean_rel_abun = mean(100*sum(rel_abun)), .groups="drop")

#arrange rel_abun in increasing order
cfvb_its_genus_pool_var_ins_p3 <- cfvb_its_genus_var_ins_p3 |>  
  group_by(taxon) |>  
  summarise(pool=max(mean_rel_abun) < 2.85, .groups = "drop") #|>  
  #arrange(desc(max))

cfvb_its_genus_var_ins_plot_p3 <- inner_join(
  cfvb_its_genus_var_ins_p3, cfvb_its_genus_pool_var_ins_p3,
  by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |>  
  group_by(spray, taxon) |>  
  ggplot(aes(x=spray, y=mean_rel_abun, fill=taxon)) + 
  geom_col(show.legend = T) +
  facet_wrap(~variety, scales = "free") +
  scale_y_continuous(expand = c(0,0),limits=c(0,100)) +
  scale_fill_manual(
    name=NULL,
    breaks = c("Alternaria","Cladosporium","Filobasidium","Gibberella",
               "Neopestalotiopsis","Nigrospora","Other",'Pseudopeyronellaea',
               "Sporobolomyces","Unclassified Ascomycota",
               "Unclassified Nectriaceae"),
    labels = c("*Alternaria*","*Cladosporium*","*Filobasidium*","*Gibberella*",
               "*Neopestalotiopsis*","*Nigrospora*","Other",
               '*Pseudopeyronellaea*',"*Sporobolomyces*",
               "Unclassified *Ascomycota*","Unclassified *Nectriaceae*"),
    values = c("#708090","#99DDFF","#f3f3f3","#b6d7a8","#6aa84f","#8200ff",
               "#545454","#000000","#bcbcbc","#f0e1d5","#a384b4","#eedd88",
               "#6e8b3d","#9ACD32","#F4A460","#ffffff","#77aadd","#696969",
               "#0000cd","#ffaabb","#B452CD","#000000","#ee8866","#bbbbbb")
  ) +
  theme_bw() +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=9,family="sans",color="#000000"),
        plot.title = element_text(hjust=0.5),
        panel.border = element_blank(),
        axis.text = element_text(size=8,family="sans",color="black"),
        axis.title.y = element_text(size=9,family="sans",color="#000000"),
        legend.text = element_markdown(size=8,family="sans",color="#000000"),
        legend.key.size = unit(8, "pt"),
        legend.spacing.y = unit(.5, "cm"),
        axis.title.x = element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x=NULL, y="Mean Relative Abundance (%)",
       title = " ")
ggsave("figures/p3_figures/cfvb_its_genus_var_ins_plot_p3_new.pdf", 
       cfvb_its_genus_var_ins_plot_p3, width = 5, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_genus_var_ins_plot_p3_new.svg", 
       cfvb_its_genus_var_ins_plot_p3, width = 5, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_genus_var_ins_plot_p3_new.tiff", 
       cfvb_its_genus_var_ins_plot_p3, width = 5, height = 4, dpi = 300)

# b. fungal genus by insecticide b/w CF and VB unstacked bar chats

cfvb_its_genus_var_ins_plot_p3b<- inner_join(
  cfvb_its_genus_var_ins_p3, cfvb_its_genus_pool_var_ins_p3,
  by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |>  
  mutate(taxon=factor(taxon),
         taxon=fct_reorder(taxon, mean_rel_abun, .desc=F)) |> 
  group_by(spray, taxon) |>  
  ggplot(aes(x=taxon,y=mean_rel_abun, fill=spray)) + 
  geom_col(show.legend = T, position=position_dodge()) +
  facet_wrap(~variety, scales = "free") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(breaks=c("+Max","-Max"),
                    labels=c("+Max", "-Max"),
                    values=c("#708090","#0000cd"))+
  scale_x_discrete(
    name=NULL,
    breaks = c("Alternaria","Cladosporium","Filobasidium","Gibberella",
               "Neopestalotiopsis","Nigrospora","Other",'Pseudopeyronellaea',
               "Sporobolomyces","Unclassified Ascomycota",
               "Unclassified Nectriaceae"),
    labels = c("*Alternaria*","*Cladosporium*","*Filobasidium*","*Gibberella*",
               "*Neopestalotiopsis*","*Nigrospora*","Other",
               '*Pseudopeyronellaea*',"*Sporobolomyces*",
               "Unclassified *Ascomycota*","Unclassified *Nectriaceae*")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=9,family="sans",color="#000000"),
        plot.title = element_text(hjust=0.5),
        panel.border = element_blank(),
        panel.grid.major.x = element_line(),
        axis.text.y = element_markdown(size=8,family="sans",color="black"),
        axis.title.y = element_text(size=8,family="sans",color="#000000"),
        legend.text = element_text(size=8,family="sans",color="#000000"),
        legend.key.size = unit(9, "pt"),
        legend.spacing.y = unit(.5, "cm"),
        axis.ticks.y = element_line()) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x=NULL, y="Mean Relative Abundance (%)",
       title = " ")
ggsave("figures/p3_figures/cfvb_its_genus_var_ins_plot_p3b_new.pdf", 
       cfvb_its_genus_var_ins_plot_p3b, width = 7, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_genus_var_ins_plot_p3b_new.svg", 
       cfvb_its_genus_var_ins_plot_p3b, width = 7, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_genus_var_ins_plot_p3b_new.tiff", 
       cfvb_its_genus_var_ins_plot_p3b, width = 7, height = 4, dpi = 300)


# c. fungal genus by netting b/w CF and VB stacked bar chats

cfvb_its_rel_abun_genus_var_net_p3 <- cfvb_its_merged_long_p3_ngs_cul |> 
  select(-c(bagging,insecticide)) |> 
  group_by(bag, variety) |> 
  mutate(rel_abun=counts/sum(counts)) |> 
  ungroup() |> 
  filter(method != "cul") |> 
  select(-c(otu, counts, method)) |> 
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", 
                        "family", "genus","species", "asv"), 
               names_to = "level", values_to = "taxon")

cfvb_its_genus_var_net_p3 <- cfvb_its_rel_abun_genus_var_net_p3 %>% 
  filter(level=="genus") %>% 
  filter(rel_abun > 0) |> 
  group_by(bag, taxon, variety) |> 
  summarise(mean_rel_abun = mean(100*sum(rel_abun)), .groups="drop")

#arrange rel_abun in increasing order
cfvb_its_genus_pool_var_net_p3 <- cfvb_its_genus_var_net_p3 |>  
  group_by(taxon) |>  
  summarise(pool=max(mean_rel_abun) < 2.9, .groups = "drop") #|>  
  #arrange(desc(max))

cfvb_its_genus_var_net_plot_p3 <- inner_join(
  cfvb_its_genus_var_net_p3, cfvb_its_genus_pool_var_net_p3,
  by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |>  
  group_by(bag, taxon) |>  
  ggplot(aes(x=bag, y=mean_rel_abun, fill=taxon)) + 
  geom_col(show.legend = T) +
  facet_wrap(~variety, scales = "free") +
  scale_y_continuous(expand = c(0,0), limits=c(0,100)) +
  scale_fill_manual(
    name=NULL,
    breaks = c("Alternaria","Cladosporium","Filobasidium","Gibberella",
               "Neopestalotiopsis","Nigrospora","Other",'Pseudopeyronellaea',
               "Sporobolomyces","Unclassified Ascomycota",
               "Unclassified Nectriaceae"),
    labels = c("*Alternaria*","*Cladosporium*","*Filobasidium*","*Gibberella*",
               "*Neopestalotiopsis*","*Nigrospora*","Other",
               '*Pseudopeyronellaea*',"*Sporobolomyces*",
               "Unclassified *Ascomycota*","Unclassified *Nectriaceae*"),
    values = c("#708090","#99DDFF","#f3f3f3","#b6d7a8","#6aa84f","#8200ff",
               "#545454","#000000","#bcbcbc","#f0e1d5","#a384b4","#eedd88",
               "#6e8b3d","#9ACD32","#F4A460","#ffffff","#77aadd","#696969",
               "#0000cd","#ffaabb","#B452CD","#000000","#ee8866","#bbbbbb")
  ) +
  theme_bw() +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=9, family="sans", color="#000000"),
        plot.title = element_text(hjust=0.5),
        panel.border = element_blank(),
        axis.text = element_text(size=8,family="sans",color="#000000"),
        axis.title.y = element_text(size=9,family="sans",color="#000000"),
        legend.text = element_markdown(size=8,family="sans",color="#000000"),
        legend.key.size = unit(8, "pt"),
        legend.spacing.y = unit(.5, "cm"),
        axis.title.x = element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x=NULL, y="Mean Relative Abundance (%)",
       title = " ")
ggsave("figures/p3_figures/cfvb_its_genus_var_net_plot_p3_new.pdf", 
       cfvb_its_genus_var_net_plot_p3, width = 5, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_genus_var_net_plot_p3_new.svg", 
       cfvb_its_genus_var_net_plot_p3, width = 5, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_genus_var_net_plot_p3_new.tiff", 
       cfvb_its_genus_var_net_plot_p3, width = 5, height = 4, dpi = 300)

# d. fungal genus by netting b/w CF and VB unstacked bar chats

cfvb_its_genus_var_net_plot_p3b <- inner_join(
  cfvb_its_genus_var_net_p3, cfvb_its_genus_pool_var_net_p3,
  by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |>  
  mutate(taxon=factor(taxon),
         taxon=fct_reorder(taxon, mean_rel_abun, .desc=F)) |>
  group_by(bag, taxon) |>  
  ggplot(aes(x=taxon, y=mean_rel_abun, fill=bag)) + 
  geom_col(show.legend = T, position=position_dodge()) +
  facet_wrap(~variety, scales = "free") +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(breaks=c("+Net","-Net"),
                    labels=c("+Net", "-Net"),
                    values=c("#708090","#0000cd"))+
  scale_x_discrete(
    name=NULL,
    breaks = c("Alternaria","Cladosporium","Filobasidium",
               "Neopestalotiopsis","Nigrospora","Other",'Pseudopeyronellaea',
               "Sporobolomyces","Unclassified Ascomycota",
               "Unclassified Saccharomycetales"),
    labels = c("*Alternaria*","*Cladosporium*","*Filobasidium*",
               "*Neopestalotiopsis*","*Nigrospora*","Other",
               '*Pseudopeyronellaea*',"*Sporobolomyces*",
               "Unclassified *Ascomycota*","Unclassified *Saccharomycetales*"))+
  theme_classic() +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=9,family="sans",color="#000000"),
        plot.title = element_text(hjust=0.5),
        panel.border = element_blank(),
        panel.grid.major.x = element_line(),
        axis.text.y = element_markdown(size=8,family="sans",color="black"),
        axis.title.y = element_text(size=8,family="sans",color="#000000"),
        legend.text = element_text(size=8,family="sans",color="#000000"),
        legend.key.size = unit(9, "pt"),
        legend.spacing.y = unit(.5, "cm"),
        axis.ticks = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x=NULL, y="Mean Relative Abundance (%)",
       title = " ")
ggsave("figures/p3_figures/cfvb_its_genus_var_net_plot_p3b_new.pdf", 
       cfvb_its_genus_var_net_plot_p3b, width = 7, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_genus_var_net_plot_p3b_new.svg", 
       cfvb_its_genus_var_net_plot_p3b, width = 7, height = 4, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_genus_var_net_plot_p3b_new.tiff", 
       cfvb_its_genus_var_net_plot_p3b, width = 7, height = 4, dpi = 300)


####################### Venn Diagrams and heat maps ############################

# a. ITS net core biomes data

read_excel("data/cfvb_its_taxonomy_table_p3_level-7_2.xlsx",3) |>
  pivot_longer(cols = c(
    "NNet...2","NNet...3","Net...4","Net...5","NNet...6","NNet...7","Net...8",
    "Net...9","NNet...10","NNet...11","Net...12","Net...13","NNet...14",
    "NNet...15","Net...16","Net...17","NNet...18","NNet...19","Net...20",
    "Net...21","NNet...22","NNet...23","Net...24","Net...25"),
    names_to='bagging', values_to='counts') %>%
  write.table(., file = "results/cfvb_its_tab_4_cor_biom_bag_long_p3_new.csv",
              col.names = NA, row.names = T, sep = ",", quote = F)


# b. 16S insecticide core biomes data

read_excel("data/cfvb_its_taxonomy_table_p3_level-7.xlsx", 4) |>
  pivot_longer(cols = c(
    "NMax...2","MMax...3","NMax...4","MMax...5","NMax...6","MMax...7",
    "NMax...8","MMax...9","NMax...10","MMax...11","NMax...12","MMax...13",
    "NMax...14","MMax...15","NMax...16","MMax...17","NMax...18","MMax...19",
    "NMax...20","MMax...21","NMax...22","MMax...23","NMax...24","MMax...25"),
    names_to='insecticide', values_to='counts') %>%
  write.table(., file = "results/cfvb_its_tab_4_cor_biom_ins_long_p3_new.csv",
              col.names = NA, row.names = T, sep = ",", quote = F)

# c. 16S variety core biomes data

read_excel("data/cfvb_its_taxonomy_table_p3_level-7.xlsx", 5) |>
  pivot_longer(cols = c(
    "CF...2","CF...3","CF...4","CF...5","CF...6","CF...7","CF...8","CF...9",
    "CF...10","CF...11","CF...12","CF...13","VB...14","VB...15","VB...16",
    "VB...17","VB...18","VB...19","VB...20","VB...21","VB...22","VB...23",
    "VB...24","VB...25"),
    names_to='variety', values_to='counts') %>%
  write.table(., file = "results/cfvb_its_tab_4_cor_biom_var_long_p3_new.csv",
              col.names = NA, row.names = T, sep = ",", quote = F)


# venn diagrams from the combined ngs and cul data without otus

read_excel("results/cfvb_its_tab_4_cor_biom_var_long_p3_edited_new.xlsx", 1) |>
  rename_all(tolower) |> 
  filter(counts > 0) |> 
  filter(variable=="CF") |> 
  arrange(genus) |>  
  select(genus,variable) |> 
  distinct(genus,variable) %>%
  write.table(., file="results/cfvb_its_cf_distinct_genus_p3_new.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)

read_excel("results/cfvb_its_tab_4_cor_biom_var_long_p3_edited_new.xlsx", 1) |>
  rename_all(tolower) |> 
  filter(counts > 0) |> 
  filter(variable=="VB") |> 
  arrange(genus) |>  
  select(genus,variable) |> 
  distinct(genus,variable) %>%
  write.table(., file="results/cfvb_its_vb_distinct_genus_p3_new.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)

read_excel("results/cfvb_its_tab_4_cor_biom_var_long_p3_edited_new.xlsx", 1) |>
  rename_all(tolower) |>
  filter(counts > 0) |>
  filter(variable=="+Max") |>
  arrange(genus) |>
  select(genus,variable) |>
  distinct(genus,variable) %>%
  write.table(., file="results/cfvb_its_mmax_distinct_genus_p3_new.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)

read_excel("results/cfvb_its_tab_4_cor_biom_var_long_p3_edited_new.xlsx", 1) |>
  rename_all(tolower) |>
  filter(counts > 0) |>
  filter(variable=="-Max") |>
  arrange(genus) |>
  select(genus,variable) |>
  distinct(genus,variable) %>%
  write.table(., file="results/cfvb_its_nmax_distinct_genus_p3_new.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)

read_excel("results/cfvb_its_tab_4_cor_biom_var_long_p3_edited_new.xlsx", 1) |>
  rename_all(tolower) |>
  filter(counts > 0) |>
  filter(variable=="+Net") |>
  arrange(genus) |>
  select(genus,variable) |>
  distinct(genus,variable) %>%
  write.table(., file="results/cfvb_its_net_distinct_genus_p3_new.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)

read_excel("results/cfvb_its_tab_4_cor_biom_var_long_p3_edited_new.xlsx", 1) |>
  rename_all(tolower) |>
  filter(counts > 0) |>
  filter(variable=="-Net") |>
  arrange(genus) |>
  select(genus,variable) |>
  distinct(genus,variable) %>%
  write.table(., file="results/cfvb_its_nnet_distinct_genus_p3_new.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)


# import edited cfvb_its_cf_vb_max_nmax_net_nnet_genus_merged data 4 venn diag.

cfvb_its_cf_vb_max_nmax_net_nnet_genus_merged <- read.table(
  "results/cfvb_its_cf_vb_max_nmax_net_nnet_genus_merged_new.txt", 
  header=T, sep="\t") #|> view()

# creating list vector of venn diagram for cf_vb_mmax_nmax
cfvb_its_cf_vb_max_nmax_genus_merged_venn_list = list(
  CF = cfvb_its_cf_vb_max_nmax_net_nnet_genus_merged$genus[1:54],
  VB = cfvb_its_cf_vb_max_nmax_net_nnet_genus_merged$genus[55:95]
)
# Ploting the venn diagram for shared mycobiomes between -Max and +Max
cfvb_its_cf_vb_max_nmax_genus_merged_venn_plot <- ggVennDiagram(
  cfvb_its_cf_vb_max_nmax_genus_merged_venn_list, label_alpha = 0,
  category.names = c("CF", "VB")
) +
  scale_fill_gradient(low="#ffffff",high = "#99ddff") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.ticks.y.left = element_blank(),
    legend.position = "None"
  )
#upon repeat of analysis, venn diag remained the same as previous

ggsave("figures/p3_figures/cfvb_its_cf_vb_max_nmax_genus_venn_plot.pdf",
       cfvb_its_cf_vb_max_nmax_genus_merged_venn_plot, 
       height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_cf_vb_max_nmax_genus_venn_plot.tiff",
       cfvb_its_cf_vb_max_nmax_genus_merged_venn_plot, 
       height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_cf_vb_max_nmax_genus_venn_plot.svg",
       cfvb_its_cf_vb_max_nmax_genus_merged_venn_plot, 
       height = 4, width = 5, dpi = 300)


# creating list vector of venn diagram for cf_vb_net_nnet
cfvb_its_cf_vb_net_nnet_genus_merged_venn_list = list(
  CF = cfvb_its_cf_vb_max_nmax_net_nnet_genus_merged$genus[1:54],
  VB = cfvb_its_cf_vb_max_nmax_net_nnet_genus_merged$genus[55:95],
  Net = cfvb_its_cf_vb_max_nmax_net_nnet_genus_merged$genus[184:227],
  NNet = cfvb_its_cf_vb_max_nmax_net_nnet_genus_merged$genus[228:274]
)
# Ploting the venn diagram for shared mycobiomes between -Net and +Net
cfvb_its_cf_vb_net_nnet_genus_merged_venn_plot <- ggVennDiagram(
  cfvb_its_cf_vb_net_nnet_genus_merged_venn_list, label_alpha = 0,
  category.names = c("CF", "VB", "+Net", "-Net")
) +
  scale_fill_gradient(low="#ffffff",high = "#b452cd") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.ticks.y.left = element_blank(),
    legend.position = "None"
  )
ggsave("figures/p3_figures/cfvb_its_cf_vb_net_nnet_genus_venn_plot.pdf",
       cfvb_its_cf_vb_net_nnet_genus_merged_venn_plot, 
       height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_cf_vb_net_nnet_genus_venn_plot.tiff",
       cfvb_its_cf_vb_net_nnet_genus_merged_venn_plot, 
       height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_cf_vb_net_nnet_genus_venn_plot.svg",
       cfvb_its_cf_vb_net_nnet_genus_merged_venn_plot, 
       height = 4, width = 5, dpi = 300)


################################# Heatmap ####################################

cfvb_its_cf_vb_max_nmax_net_nnet_top50_core_genus <- read_excel(
  "results/cfvb_its_tab_4_cor_biom_var_long_p3_edited.xlsx", 4) |> 
  select(-`Grand Total`) |>
  select(genus=Genus, everything()) |> 
  mutate_at(c('CF','VB','+Max','-Max','+Net','-Net'), 
            ~replace_na(.,0)) # replaces all NAs with zeros

# melt converts data into long format but then we loose our rownames
cfvb_its_cf_vb_max_nmax_net_nnet_top50_core_genus_melt <- melt(
  cfvb_its_cf_vb_max_nmax_net_nnet_top50_core_genus)

# #add column for fungal OTU name
# bs_16s_top50_core_biomes_melt$bact_otu <- rep(
#   row.names(bs_16s_top50_core_biomes), 2)

#view first six rows
head(cfvb_its_cf_vb_max_nmax_net_nnet_top50_core_genus_melt)

# rescale values for all variables in melted data frame
cfvb_its_cf_vb_max_nmax_net_nnet_top50_core_genus_melt <- ddply(
  cfvb_its_cf_vb_max_nmax_net_nnet_top50_core_genus_melt, .(variable), 
  transform, rescale = rescale(value)) # rescale converts variables to a new range

# create heatmap using rescaled values
cfvb_its_cf_vb_max_nmax_net_nnet_top50_core_genus_hm <- ggplot(
  cfvb_its_cf_vb_max_nmax_net_nnet_top50_core_genus_melt, 
  aes(x=variable, y=genus)) +
  geom_tile(aes(fill = rescale), colour = "white") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradient(low = "#FFFFFF", high = "#0000cd", "#9400D3", "#CD950C",
                      expand = c(0,0)) +
  labs(x = "", y = "") +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=8, color="#000000"),
        axis.text.x = element_markdown(angle=0,hjust=.5,vjust=1,color="black")) +
  coord_fixed(ratio = 0.1) # resizes the boxes for each variable on the x-axis

ggsave("figures/p3_figures/cfvb_its_cf_vb_ins_bag_top50_core_genus_hm.pdf", 
       cfvb_its_cf_vb_max_nmax_net_nnet_top50_core_genus_hm, 
       height = 7, width = 8, dpi = 400)
ggsave("figures/p3_figures/cfvb_its_cf_vb_ins_bag_top50_core_genus_hm.tiff", 
       cfvb_its_cf_vb_max_nmax_net_nnet_top50_core_genus_hm, 
       height = 7, width = 8, dpi = 400)
ggsave("figures/p3_figures/cfvb_its_cf_vb_ins_bag_top50_core_genus_hm.svg", 
       cfvb_its_cf_vb_max_nmax_net_nnet_top50_core_genus_hm, 
       height = 7, width = 8, dpi = 400)


############################### Ordinations ####################################

# NMDS from cfvb_its_otus table

# cfvb_its_otus_p3 <- cfvb_its_otus |>
#   mutate(sample_id2=sample_id) |>
#   select(sample_id,sample_id2, everything()) |> 
#   separate(sample_id2,
#            into = c("cul_treat", "perep"),
#            sep = "(?<=[A-Za-z])(?=[0-9])") |>
#   separate(perep,into=c("period","replicate"), sep="(?<=\\d)") |>
#   select(-c(cul_treat,replicate)) |>
#   filter(period==3) |>
#   select(-period) |> #this is a tibble and will not take row names
#   as.data.frame()
# 
# rownames(cfvb_its_otus_p3) <- cfvb_its_otus_p3$sample_id
# cfvb_its_otus_p3 <- cfvb_its_otus_p3[,-1]
# 
# cfvb_its_otus_p3 <- as.matrix(cfvb_its_otus_p3) #converts to matrix
# 
# #calculating distances in vegan
# 
# set.seed(0429)
# 
# cfvb_its_dist_p3 <- vegdist(cfvb_its_otus_p3, method = "bray")
# 
# cfvb_its_nmds_p3 <- metaMDS(cfvb_its_dist_p3) |> 
#   scores() |> 
#   as_tibble(rownames="sample_id")
# 
# cfvb_its_nmds_meta_merged_p3 <- inner_join(cfvb_its_nmds_p3, cfvb_its_meta_p3,
#                                            by="sample_id") |> 
#   select(-period)
# 
# #plotting  
# cfvb_its_nmds_plot_p3 <- cfvb_its_nmds_meta_merged_p3 |> 
#   ggplot(aes(x=NMDS1, y=NMDS2, color=variety, shape=variety)) + 
#   geom_point(size=1.25) +
#   scale_color_manual(breaks=c("CF", "VB"),
#                      labels=c("CF", "VB"),
#                      values=c("#ff0000", "#0000cd"))+
#   theme_bw() +
#   theme(legend.background = element_blank(),
#         legend.text= element_text(size = 10, family = 'sans'),
#         legend.title = element_blank(),
#         panel.border = element_rect(color="#dddddd"),
#         legend.key.size = unit(.5, units = "cm"),
#         axis.text = element_text(colour = '#000000',size=9),
#         axis.title = element_text(color="#000000",size=10))
# 
# ggsave("figures/p3_figures/cfvb_its_nmds_plot_p3.tiff", 
#        cfvb_its_nmds_plot_p3, height = 4, width = 5, dpi = 300)


# OPTION 2

# a. get our metadata
cfvb_its_nmds_meta_p3 <- cfvb_its_merged_p3 |>
  select(sample_id,bag=bagging,spray=insecticide,variety) |> 
  mutate(bag=case_when(bag=="Net"~"+Net",
                       bag=="Nnet"~"-Net"),
         spray=case_when(spray=="Mmax"~"+Max",
                         spray=="Nmax"~"-Max"))

# b. get our otus
cfvb_its_nmds_otu_p3 <- cfvb_its_merged_p3 |> 
  select(-c(2:9)) |>
  pivot_longer(-sample_id) |> 
  group_by(sample_id) |> 
  # summarise(N=sum(value)) |>  #this gives the total number of seqs in each
  # arrange(N) |> print(n=25) #sample and we arrange to get a value to rarefy
  mutate(N=sum(value)) |> #this shows the rarefied value
  filter(N >= 1400) |> #rarefying to keep all samples having 1400 seqs
  select(-N) |> 
  pivot_wider(names_from="name", values_from="value", values_fill=0)

# c. get our distance matrix
cfvb_its_nmds_dist_p3 <- cfvb_its_nmds_otu_p3 |> 
  column_to_rownames("sample_id") |> 
  avgdist(sample=1400) #calculates the bray curtis distances

# .d make ordination
set.seed(2) #required for the nmds to converge
cfvb_its_nmds_tbl_p3 <- metaMDS(cfvb_its_nmds_dist_p3) |> 
  scores() |> 
  as_tibble(rownames="sample_id") #converts to a tibble that we can use in ggplot

# e. merge nmds_meta and nmds_tbl
cfvb_its_nmds_4_plot_p3 <- inner_join(cfvb_its_nmds_meta_p3, 
                                      cfvb_its_nmds_tbl_p3, by="sample_id")

# test of significance before plotting
# ADONIS
# bagging
adon_bag <- adonis2(cfvb_its_nmds_dist_p3~cfvb_its_nmds_meta_p3$bag,permutations=9999)
str(adon_bag) #check where the p-value is
adon_bag$`Pr(>F)`[1] #since it's a list, we put [1] to extract the first item

# insecticide
adon_ins <- adonis2(cfvb_its_nmds_dist_p3~cfvb_its_nmds_meta_p3$spray,permutations=9999)
str(adon_ins) #check where the p-value is
adon_ins$`Pr(>F)`[1] #since it's a list, we put [1] to extract the first item

#interaction of bagging and insecticide
adonis2(cfvb_its_nmds_dist_p3~cfvb_its_nmds_meta_p3$spray*
          cfvb_its_nmds_meta_p3$bag,permutations=9999)
# variety
adon_var <- adonis2(cfvb_its_nmds_dist_p3~cfvb_its_nmds_meta_p3$variety,
                    permutations=9999)
# get p-value

str(adon_var) #check where the p-value is
adon_var$`Pr(>F)`[1] #since it's a list, we put [1] to extract the first item

# Another statistical test of dispersion is using betadisper() in vegan

betadisp_var <- anova(betadisper(cfvb_its_nmds_dist_p3,
                                 cfvb_its_nmds_meta_p3$bag))

str(betadisp_var)
betadisp_var$`Pr(>F)`[1]

#OR
betadisp_var_permutest <- permutest(betadisper(cfvb_its_nmds_dist_p3,
                                               cfvb_its_nmds_meta_p3$bag))
str(betadisp_var_permutest)
betadisp_var_permutest$tab$`Pr(>F)`[1]

# we prepare a centroid for our ellipses
cfvb_its_centroid_p3<- cfvb_its_nmds_4_plot_p3 |> 
  group_by(variety) |> 
  summarise(NMDS1=mean(NMDS1),
            NMDS2=mean(NMDS2))

# plotting
cfvb_its_nmds_plot_p3 <- cfvb_its_nmds_4_plot_p3|>
  ggplot(aes(x=NMDS1, y=NMDS2, color=variety)) +
  geom_point() +
  stat_ellipse(show.legend=F) +
  geom_point(data=cfvb_its_centroid_p3, size=3,shape=21,
             color="#000000",aes(fill=variety),show.legend=F) +
  scale_color_manual(values=c("#0000cd","#ff0000")) +
  scale_fill_manual(values=c("#0000cd","#ff0000")) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.background = element_blank(),
        legend.text= element_text(size = 10, family = 'sans'),
        legend.title = element_blank(),
        legend.key.size = unit(.5, units = "cm"),
        axis.text = element_text(colour = 'black'))

ggsave("figures/p3_figures/cfvb_its_nmds_plot_p3_new.pdf", 
       cfvb_its_nmds_plot_p3, height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_nmds_plot_p3_new.svg", 
       cfvb_its_nmds_plot_p3, height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_nmds_plot_p3_new.tiff", 
       cfvb_its_nmds_plot_p3, height = 4, width = 5, dpi = 300)


# plotting
# bagging effect
# we prepare a centroid for our ellipses
cfvb_its_centroid_p3_bag<- cfvb_its_nmds_4_plot_p3 |> 
  group_by(bag) |> 
  summarise(NMDS1=mean(NMDS1),
            NMDS2=mean(NMDS2))

cfvb_its_nmds_plot_p3_bag <- cfvb_its_nmds_4_plot_p3|>
  ggplot(aes(x=NMDS1, y=NMDS2, color=bag)) +
  geom_point() +
  stat_ellipse(show.legend=F) +
  geom_point(data=cfvb_its_centroid_p3_bag, size=3,shape=21,
             color="#000000",aes(fill=bag),show.legend=F) +
  scale_color_manual(values=c("#0000cd","#ff0000")) +
  scale_fill_manual(values=c("#0000cd","#ff0000")) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.background = element_blank(),
        legend.text= element_text(size = 10, family = 'sans'),
        legend.title = element_blank(),
        legend.key.size = unit(.5, units = "cm"),
        axis.text = element_text(colour = 'black'))

ggsave("figures/p3_figures/cfvb_its_nmds_plot_p3_bag.pdf", 
       cfvb_its_nmds_plot_p3_bag, height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_nmds_plot_p3_bag.svg", 
       cfvb_its_nmds_plot_p3_bag, height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_nmds_plot_p3_bag.tiff", 
       cfvb_its_nmds_plot_p3_bag, height = 4, width = 5, dpi = 300)


# insecticide effect
# we prepare a centroid for our ellipses
cfvb_its_centroid_p3_ins<- cfvb_its_nmds_4_plot_p3 |> 
  group_by(spray) |> 
  summarise(NMDS1=mean(NMDS1),
            NMDS2=mean(NMDS2))

cfvb_its_nmds_plot_p3_ins <- cfvb_its_nmds_4_plot_p3|>
  ggplot(aes(x=NMDS1, y=NMDS2, color=spray)) +
  geom_point() +
  stat_ellipse(show.legend=F) +
  geom_point(data=cfvb_its_centroid_p3_ins, size=3,shape=21,
             color="#000000",aes(fill=spray),show.legend=F) +
  scale_color_manual(values=c("#0000cd","#ff0000")) +
  scale_fill_manual(values=c("#0000cd","#ff0000")) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.background = element_blank(),
        legend.text= element_text(size = 10, family = 'sans'),
        legend.title = element_blank(),
        legend.key.size = unit(.5, units = "cm"),
        axis.text = element_text(colour = 'black'))

ggsave("figures/p3_figures/cfvb_its_nmds_plot_p3_ins.pdf", 
       cfvb_its_nmds_plot_p3_ins, height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_nmds_plot_p3_ins.svg", 
       cfvb_its_nmds_plot_p3_ins, height = 4, width = 5, dpi = 300)
ggsave("figures/p3_figures/cfvb_its_nmds_plot_p3_ins.tiff", 
       cfvb_its_nmds_plot_p3_ins, height = 4, width = 5, dpi = 300)




values = c("#708090","#99DDFF","#f3f3f3","#b6d7a8","#6aa84f","#8200ff",
           "#545454","#000000","#bcbcbc","#f0e1d5","#a384b4","#eedd88",
           "#6e8b3d","#9ACD32","#F4A460","#ffffff","#77aadd","#696969",
           "#0000cd","#ffaabb","#B452CD","#ee8866","#bbbbbb","#ff0000",
           "#6495ed","#00fa9a","#00ced1","#aa4499","#d3d3d3","#8b008b",
           "#44bb99","#00ff00")
