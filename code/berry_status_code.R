source("F:/professdev/r_projects/gsr_analysis/code/m_biome_packages.R")

conflicts_prefer(dplyr::filter) # this line will prefer filter from dplyr
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)
conflicts_prefer(ggpubr::mutate)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(dplyr::desc)

########################## Berry Status Analysis ###############################

######################## 16s Analysis ##################################

# rename readme and license files to lowercase ----------------------------
file.rename(from='README.md', to = 'readme.md')


#metadata

bs_lsu_meta <- read_tsv("data/bs_16s_gsr_metadata.tsv") |>
  select('sample_id'='sample-id','berry_status'='berry-status') |>
  filter(sample_id != "ctrl") |>   #drop the control sample
  mutate(berry_status = case_when(berry_status == "healthy" ~ "H",
                                  berry_status == "Infected" ~ "I"))

bs_lsu_meta$berry_status <- as.factor(bs_lsu_meta$berry_status)

glimpse(bs_lsu_meta)

#evenness
bs_lsu_even <- read_tsv("data/bs_16s_evenness.tsv") |>  
  rename("sample_id" = "...1", 'evenness'='pielou_evenness') #renames col 1

#faith
bs_lsu_faith <- read_tsv("data/bs_16s_faith_pd.tsv") |>  
  rename("sample_id" = "#SampleID") #renames col 1

#observed_otus
bs_lsu_obs_otus <- read_tsv("data/bs_16s_observed_otus.tsv") |>  
  rename("sample_id" = "...1", "obs_otus"="observed_features") #rename columns

#shannon
bs_lsu_shan <- read_tsv("data/bs_16s_shannon.tsv") |>  
  rename("sample_id" = "...1", "shannon"="shannon_entropy")

#otu
bs_lsu_otus <- read_tsv("data/bs_16s_rarefied_table_transp.txt", skip = 1) |> 
  rename("sample_id"="#OTU ID")

#ace
bs_lsu_ace <- read_tsv("data/bs_16s_ace.tsv") |>  
  rename("sample_id" = "...1") #rename columns

#chao1
bs_lsu_chao1 <- read_tsv("data/bs_16s_chao1.tsv") |>  
  rename("sample_id" = "...1")

#merge metadata, ace, chao1, evenness, faith, observed_otus, shannon and otus

bs_lsu_merged <- inner_join(bs_lsu_meta, bs_lsu_chao1, by='sample_id') %>%
  inner_join(., bs_lsu_ace, by='sample_id') %>%
  inner_join(., bs_lsu_even, by='sample_id') %>%
  inner_join(., bs_lsu_faith, by='sample_id') %>%
  inner_join(., bs_lsu_obs_otus, by='sample_id') %>% 
  inner_join(., bs_lsu_shan, by='sample_id') %>%
  inner_join(., bs_lsu_otus, by='sample_id')

# Normality Check

#Create 2x2 plot environment 
par(mfrow = c(2, 3))

#Plots
hist(bs_lsu_merged$chao1, main="Chao1 diversity", xlab="", breaks=10, 
     border="#000000", col="#99DDFF", las=1)
hist(bs_lsu_merged$ace, main="Ace diversity", xlab="", breaks=10, 
     border="#000000", col="#BBCC33", las=1)
hist(bs_lsu_merged$shannon, main="Shannon diversity", xlab="", breaks=10, 
     border="#000000", col="#FFAABB", las=1)
hist(bs_lsu_merged$evenness, main="Pielou's evenness", xlab="", breaks=10, 
     border="#000000", col="#44BB99", las=1)
hist(bs_lsu_merged$faith_pd, main="Phylegenetic diversity", xlab="", 
     breaks=10, border="#000000", col="#DDDDDD", las=1)
hist(bs_lsu_merged$obs_otus, main="Community richness", xlab="", 
     breaks=10, border="#000000", col="#EEDD88", las=1)


# Test for normalcy

shapiro.test(bs_lsu_merged$ace)        # normally dist
shapiro.test(bs_lsu_merged$chao1)      
shapiro.test(bs_lsu_merged$shannon)    # normally dist
shapiro.test(bs_lsu_merged$evenness)
shapiro.test(bs_lsu_merged$faith_pd)   # normally dist
shapiro.test(bs_lsu_merged$obs_otus)   # normally dist


# a. Anova test for normally distributed data

# i-ace
lsu.ace.bs.aov <- aov(ace ~ berry_status, data=bs_lsu_merged)
summary(lsu.ace.bs.aov)

# ii-chao1
kruskal.test(chao1 ~ berry_status, data=bs_lsu_merged)

# iii-shannon
lsu.shan.bs.aov <- aov(shannon ~ berry_status, data=bs_lsu_merged)
summary(lsu.shan.bs.aov)

# iv-evenness
kruskal.test(evenness ~ berry_status, data=bs_lsu_merged)

# v-Faith
lsu.fd.bs.aov <- aov(faith_pd ~ berry_status, data=bs_lsu_merged)
summary(lsu.fd.bs.aov)

# vi-Observed_otus
lsu.obs.bs.aov <- aov(obs_otus ~ berry_status, data=bs_lsu_merged)
summary(lsu.obs.bs.aov)


# a. chao1, shannon, otus and faith by berry status

my_labs <- c("shannon" = "shannon",
             "chao1" = "chao1",
             "faith_pd" = "faith",
             "obs_otus" = "observed otus") # labels for y-axis of boxplots

bs_lsu_alpha_div_plot <- bs_lsu_merged |> 
  select(-c(ace,evenness)) |> 
  pivot_longer(cols = 3:6, names_to = 'indice', values_to = 'value') |> 
  mutate(indice=factor(indice, levels = c("chao1","shannon", 
                                          "obs_otus", "faith_pd"))) |>  
  select(berry_status, indice, value) |>  
  group_by(berry_status, indice) %>% 
  ggplot(., aes(x=berry_status, y=value)) + 
  geom_boxplot(aes(fill=berry_status),size=.5, outlier.shape = NA) +
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
    strip.text = element_text(size = 12, colour = "black", family = "sans"),
    panel.background = element_blank(),
    axis.title.y = element_text(size = 12, family = "sans"),
    axis.text = element_text(size = 11, color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(hjust = .5, vjust = 1, family = "sans"))

ggsave("figures/bs_lsu_alpha_div_plot.tiff", bs_lsu_alpha_div_plot, 
       height = 4.5, width = 5, dpi = 300)
ggsave("figures/bs_lsu_alpha_div_plot.pdf", bs_lsu_alpha_div_plot, 
       height = 4.5, width = 5, dpi = 300)
ggsave("figures/bs_lsu_alpha_div_plot.svg", bs_lsu_alpha_div_plot, 
       height = 4.5, width = 5, dpi = 300)


####################### import bs_taxonomy file for 16S ########################

bs_lsu_tax <- read.table("data/bs_16s_taxonomy.tsv", 
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
bs_lsu_tax[is.na(bs_lsu_tax)] <- ""
#convert anything identified as "unidentified" into an empty space
bs_lsu_tax[bs_lsu_tax=="Unclassified"] <- "" 

#following code replaces those empty spaces with the Unclassified highest taxon
for (i in 1:nrow(bs_lsu_tax)) {
  if (bs_lsu_tax[i,7] != ""){
    bs_lsu_tax$Species[i] <- paste(bs_lsu_tax$Genus[i], sep = " ")
  }else if (bs_lsu_tax[i,2] == ""){
    Kingdom <- paste("Unclassified", bs_lsu_tax[i,1], sep = " ")
    bs_lsu_tax[i, 2:7] <- Kingdom
  }else if (bs_lsu_tax[i,3] == ""){
    Phylum <- paste("Unclassified", bs_lsu_tax[i,2], sep = " ")
    bs_lsu_tax[i, 3:7] <- Phylum
  }else if (bs_lsu_tax[i,4] == ""){
    Class <- paste("Unclassified", bs_lsu_tax[i,3], sep = " ")
    bs_lsu_tax[i, 4:7] <- Class
  }else if (bs_lsu_tax[i,5] == ""){
    Order <- paste("Unclassified", bs_lsu_tax[i,4], sep = " ")
    bs_lsu_tax[i, 5:7] <- Order 
  }else if (bs_lsu_tax[i,6] == ""){
    Family <- paste("Unclassified", bs_lsu_tax[i,5], sep = " ")
    bs_lsu_tax[i, 6:7] <- Family
  }else if (bs_lsu_tax[i,7] == ""){
    bs_lsu_tax$Species[i] <- paste("Unclassified", 
                                   bs_lsu_tax$Genus[i], sep = " ")
  }
}

#convert lsu_merged to long format and join to the cleaned taxonomy table

bs_lsu_merged %>%
  select(-c(5:10)) %>% 
  pivot_longer(-c(sample_id, berry_status), names_to = "OTU", 
               values_to = "Counts") %>% 
  inner_join(., bs_lsu_tax, by=c('OTU'='Feature.ID')) %>% 
  write.table(., "results/bs_lsu_merged_long.tsv", sep="\t", 
              row.names=FALSE, quote = FALSE)


# preparing bacteria phyla plot by berry_status

bs_lsu_merged_long_ngs_cul <- read.table(
  "results/bs_lsu_merged_long_ngs_cul_edited.txt",head=T,sep="\t"
) |> 
  rename_all(tolower) |> 
  group_by(berry_status) |> 
  mutate(rel_abun=counts/sum(counts)) |> 
  ungroup() |>
  filter(method != "cul") |> #drop culturable data
  select(-c(otu,counts,method,asv,kingdom,class,order,family,species)) |> 
  pivot_longer(cols = c("phylum","genus"), 
               names_to = "level", values_to = "taxon")


#Plotting

# a. bacterial bar chart phyla berry_status

bs_lsu_phylum_rel_abun <- bs_lsu_merged_long_ngs_cul %>% 
  filter(level=="phylum") %>% 
  group_by(berry_status, taxon) |> 
  summarise(mean_rel_abun = mean(100*sum(rel_abun)),.groups="drop")

#arrange rel_abun in increasing order
bs_lsu_phylum_pool <- bs_lsu_phylum_rel_abun %>% 
  group_by(taxon) %>% 
  summarise(pool=max(mean_rel_abun) < 4, .groups = "drop") #%>% 
  #arrange(desc(max))

bs_lsu_phylum_plot <- inner_join(
  bs_lsu_phylum_rel_abun, bs_lsu_phylum_pool, by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |>  
  group_by(berry_status, taxon) %>% 
  ggplot(., aes(x=berry_status, y=mean_rel_abun, fill=taxon)) + 
  geom_col(show.legend = T) +
  scale_y_continuous(limits = c(0,100),expand = c(0,0)) +
  scale_fill_manual(
    name=NULL,
    breaks = c("Actinobacteriota","Bacteroidota","Deinococcota",
               "Firmicutes","Other","Pseudomonadota","Unclassified Bacteria"),
    labels = c("*Actinobacteriota*","*Bacteroidota*","*Deinococcota*",
               "*Firmicutes*","Other","*Pseudomonadota*",
               "Unclassified *Bacteria*"),
    values = c("#708090","#99DDFF","#f3f3f3","#b6d7a8","#8200ff","#6aa84f",
               "#545454","#000000","#bcbcbc","#f0e1d5","#a384b4","#eedd88",
               "#6e8b3d","#9ACD32","#F4A460","#ffffff","#77aadd","#696969",
               "#0000cd","#ffaabb","#B452CD","#ee8866","#bbbbbb","#ff0000",
               "#6495ed","#00fa9a","#00ced1","#aa4499","#d3d3d3","#8b008b",
               "#44bb99","#00ff00")
  ) +
  theme_bw() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5),
        #panel.border = element_blank(),
        axis.text = element_text(size = 8, family = "sans", color="#000000"),
        axis.title.y = element_text(size = 9, family = "sans"),
        legend.text = element_markdown(size=8),
        legend.key.size = unit(8, "pt"),
        legend.spacing.y = unit(.5, "cm"),
        axis.title.x = element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x=NULL, y="Mean Relative Abundance (%)",
       title = " ")
ggsave("figures/bs_lsu_phylum_plot.pdf", 
       bs_lsu_phylum_plot, width = 4, height = 3, dpi = 300)
ggsave("figures/bs_lsu_phylum_plot.svg", 
       bs_lsu_phylum_plot, width = 4, height = 3, dpi = 300)
ggsave("figures/bs_lsu_phylum_plot.tiff", 
       bs_lsu_phylum_plot, width = 4, height = 3, dpi = 300)


# b. bacterial bar chart genera berry_status

bs_lsu_genus_rel_abun <- bs_lsu_merged_long_ngs_cul %>% 
  filter(level=="genus") %>% 
  group_by(berry_status, taxon) |> 
  summarise(mean_rel_abun = mean(100*sum(rel_abun)),.groups="drop")

#arrange rel_abun in increasing order
bs_lsu_genus_pool <- bs_lsu_genus_rel_abun %>% 
  group_by(taxon) %>% 
  summarise(pool=max(mean_rel_abun) < 7, .groups = "drop") #%>% 
  #arrange(desc(max))

bs_lsu_genus_plot <- inner_join(
  bs_lsu_genus_rel_abun, bs_lsu_genus_pool, by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |>  
  group_by(berry_status, taxon) %>% 
  ggplot(., aes(x=berry_status, y=mean_rel_abun, fill=taxon)) + 
  geom_col(show.legend = T) +
  scale_y_continuous(limits = c(0,100),expand = c(0,0)) +
  scale_fill_manual(
    name=NULL,
    breaks = c("Aeribacillus","Deinococcus","Escherichia","Gluconobacter",
               "JC017","Komagataeibacter","Other","Thermoanaerobacterium",
               "Thermus","Unclassified Anoxybacillaceae",
               "Unclassified Enterobacteriaceae"),
    labels = c("*Aeribacillus*","*Deinococcus*","*Escherichia*",
               "*Gluconobacter*","*JC017*","*Komagataeibacter*","Other",
               "*Thermoanaerobacterium*","*Thermus*",
               "Unclassified *Anoxybacillaceae*","Unclassified *Enterobacteriaceae*"),
    values = c("#708090","#99DDFF","#f3f3f3","#b6d7a8","#6aa84f","#000000",
               "#8200ff","#545454","#bcbcbc","#f0e1d5","#a384b4","#eedd88",
               "#6e8b3d","#9ACD32","#F4A460","#ffffff","#77aadd","#696969",
               "#0000cd","#ffaabb","#B452CD","#ee8866","#bbbbbb","#ff0000",
               "#6495ed","#00fa9a","#00ced1","#aa4499","#d3d3d3","#8b008b",
               "#44bb99","#00ff00")
  ) +
  theme_bw() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5),
        #panel.border = element_blank(),
        axis.text = element_text(size = 8, family = "sans", color="#000000"),
        axis.title.y = element_text(size = 9, family = "sans"),
        legend.text = element_markdown(size=8),
        legend.key.size = unit(8, "pt"),
        legend.spacing.y = unit(.5, "cm"),
        axis.title.x = element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x=NULL, y="Mean Relative Abundance (%)",
       title = " ")
ggsave("figures/bs_lsu_genus_plot.pdf", 
       bs_lsu_genus_plot, width = 4, height = 3, dpi = 300)
ggsave("figures/bs_lsu_genus_plot.svg", 
       bs_lsu_genus_plot, width = 4, height = 3, dpi = 300)
ggsave("figures/bs_lsu_genus_plot.tiff", 
       bs_lsu_genus_plot, width = 4, height = 3, dpi = 300)


# # venn diagrams from the combined ngs and cul data with otus
# bs_healthy <- bs_lsu_merged_long_ngs_cul %>% 
#   select(berry_status, level, taxon, rel_abun) %>%
#   group_by(taxon, berry_status) |> 
#   summarise(mean_rel_abun = mean(100*sum(rel_abun)), .groups = "drop") |> 
#   filter(mean_rel_abun > 0) |> 
#   filter(berry_status=="H") %>% 
#   arrange(taxon) %>% 
#   select(berry_status, taxon) |> 
#   distinct(berry_status, taxon) %>%
#   write.table(., file="results/bs_healthy_distinct_genus.tsv",
#               col.names = NA, row.names = T, sep = "\t", quote = F)
# 
# 
# bs_infected <- bs_lsu_merged_long_ngs_cul %>% 
#   select(berry_status, level, taxon, rel_abun) %>%
#   group_by(taxon, berry_status) |> 
#   summarise(mean_rel_abun = mean(100*sum(rel_abun)), .groups = "drop") |> 
#   filter(mean_rel_abun > 0) |> 
#   filter(berry_status=="I") %>% 
#   arrange(taxon) %>% 
#   select(berry_status, taxon) |> 
#   distinct(berry_status, taxon) %>%
#   write.table(., file="results/bs_infected_distinct_genus.tsv",
#               col.names = NA, row.names = T, sep = "\t", quote = F)
# 
# 
# # import edited merged_distinct_cf_vb_siteA_siteB data for venn diagram
# 
# bs_healthy_infected_genus_merged <- read.table(
#   "results/bs_healthy_infected_genus_merged.txt", 
#   header=T, sep="\t") #|> view()
# 
# # creating list vector for venn diagram
# bs_healthy_infected_genus_merged_venn_list = list(
#   H = bs_healthy_infected_genus_merged$taxon[1:52],
#   I = bs_healthy_infected_genus_merged$taxon[53:116]
# )
# # Ploting the venn diagram for shared mycobiomes across periods
# bs_healthy_infected_genus_merged_venn_plot <- ggVennDiagram(
#   bs_healthy_infected_genus_merged_venn_list, label_alpha = 0,
#   category.names = c("H", "I")
# ) +
#   scale_fill_gradient(low="#dddddd",high = "#b452cd") +
#   theme_bw() +
#   theme(
#     axis.title = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks.x.bottom = element_blank(),
#     axis.ticks.y.left = element_blank(),
#     legend.position = "None"
#   )
# 
# ggsave("figures/bs_healthy_infected_genus_merged_venn_plot.pdf",
#        bs_healthy_infected_genus_merged_venn_plot, 
#        height = 3, width = 3, dpi = 300)
# ggsave("figures/bs_healthy_infected_genus_merged_venn_plot.tiff",
#        bs_healthy_infected_genus_merged_venn_plot, 
#        height = 3, width = 3, dpi = 300)
# ggsave("figures/bs_healthy_infected_genus_merged_venn_plot.svg",
#        bs_healthy_infected_genus_merged_venn_plot, 
#        height = 3, width = 3, dpi = 300)


############## import bs_16S_taxonomy table for core biomes ###################

read_excel("data/bs_16s_taxonomy_table_level-7.xlsx", 3) %>%
  pivot_longer(cols = c(
    "Infected...2", "Infected...3","Infected...4","healthy...5","healthy...6",
    "healthy...7","healthy...8","Infected...9","Infected...10","Infected...11",
    "healthy...12","healthy...13","healthy...14","Infected...15","healthy...16",
    "healthy...17","Infected...18","Infected...19","Infected...20",
    "healthy...21","Infected...22","healthy...23","Infected...24","healthy...25",
    "Infected...26","healthy...27","Infected...28","healthy...29","Infected...30",
    "healthy...31","Infected...32","Infected...33","healthy...34","healthy...35"),
    names_to='berry_status', values_to='counts') %>%
  write.table(.,file="results/bs_16s_tab_4_cor_biome_long.csv",
              col.names = NA, row.names = T, sep = ",", quote = F)


# venn diagrams from the combined ngs and cul data without otus

read_excel("results/bs_16s_tab_4_cor_biome_long_edited.xlsx", 1) |>
  separate(Species, into=c("Genus", NA), sep=" ") |> #get just genus
  rename_all(tolower) |> 
  filter(counts > 0) |> 
  filter(berry_status=="I") |> 
  filter(method != "cul") |> 
  arrange(genus) |>  
  select(genus,berry_status) |> 
  distinct(genus,berry_status) %>%
  write.table(., file="results/bs_16s_infected_distinct_genera.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)


read_excel("results/bs_16s_tab_4_cor_biome_long_edited.xlsx", 1) |> 
  separate(Species, into=c("Genus", NA), sep=" ") |>
  rename_all(tolower) |> 
  filter(counts > 0) |> 
  filter(berry_status=="H") |> 
  filter(method != "cul") |>
  arrange(genus) |>  
  select(genus,berry_status) |> 
  distinct(genus,berry_status) %>%
  write.table(., file="results/bs_16s_helthy_distinct_genera.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)


# import edited bs_healthy_infected_species_merged data for venn diagram

bs_16s_healthy_infected_genera_merged <- read.table(
  "results/bs_16s_healthy_infected_genera_merged.txt", 
  header=T, sep="\t") #|> view()

# creating list vector for venn diagram
bs_16s_healthy_infected_genera_merged_venn_list = list(
  H = bs_16s_healthy_infected_genera_merged$genus[1:39],
  I = bs_16s_healthy_infected_genera_merged$genus[40:90]
)
# Ploting the venn diagram for shared mycobiomes across periods
bs_16s_healthy_infected_genera_merged_venn_plot <- ggVennDiagram(
  bs_16s_healthy_infected_genera_merged_venn_list, label_alpha = 0,
  category.names = c("H", "I")
) +
  scale_fill_gradient(low="#f2f2f2",high = "#6aa84f") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.ticks.y.left = element_blank(),
    legend.position = "None"
  )
ggsave("figures/bs_16s_healthy_infected_genera_merged_venn_plot.pdf",
       bs_16s_healthy_infected_genera_merged_venn_plot, 
       height = 3, width = 3, dpi = 300)
ggsave("figures/bs_16s_healthy_infected_genera_merged_venn_plot.tiff",
       bs_16s_healthy_infected_genera_merged_venn_plot, 
       height = 3, width = 3, dpi = 300)
ggsave("figures/bs_16s_healthy_infected_genera_merged_venn_plot.svg",
       bs_16s_healthy_infected_genera_merged_venn_plot, 
       height = 3, width = 3, dpi = 300)


################################# Heatmap ####################################
# my_cols <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
# 
# bs_16s_core_genera <- read_excel(
#   "results/bs_16s_tab_4_cor_biome_long_edited.xlsx", 5) |> 
#   select(-`Grand Total`)
# 
# bs_16s_df_rownames <- column_to_rownames(bs_16s_core_genera, "genus")
# 
# bs_16s_core_genera_matrix <- as.matrix(bs_16s_df_rownames[1:2])
# 
# heatmap(bs_16s_core_genera_matrix, 
#         scale = "column", # scale the values
#         Colv = NA,
#         col = my_cols) # drop the top phylogeny
# 
# library(RColorBrewer)
# coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

bs_16s_top50_core_genus <- read_excel(
  "results/bs_16s_tab_4_cor_biome_long_edited.xlsx", 5) |> 
  select(-`Grand Total`) #|> 
  #mutate_at(c('H','I'), ~replace_na(.,0)) # replaces all NAs with zeros

# melt converts data into long format but then we loose our rownames
bs_16s_top50_core_genus_melt <- melt(bs_16s_top50_core_genus)

#view first six rows
head(bs_16s_top50_core_genus_melt)

# rescale values for all variables in melted data frame
bs_16s_top50_core_genus_melt <- ddply(
  bs_16s_top50_core_genus_melt, .(variable), 
  transform, dist = scale(value)) # rescale converts variables to a new range

# create heatmap using rescaled values
bs_16s_top50_core_genus_hm <- ggplot(
  bs_16s_top50_core_genus_melt, aes(x=variable, y=genus)) +
  geom_tile(aes(fill = dist), colour = "white") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradient(low = "#FFFFFF", high = "#000000","#ff0000",
                      expand = c(0,0)) +
  labs(x = "", y = "") +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=8, color="#000000"),
        axis.text.x = element_markdown(angle=0,hjust=.5,vjust=1,color="black")) +
  coord_fixed(ratio = 0.1) # resizes the boxes for each variable on the x-axis

ggsave("figures/bs_16s_top50_core_genus_hm.pdf", 
       bs_16s_top50_core_genus_hm, height = 5, width = 4, dpi = 400)
ggsave("figures/bs_16s_top50_core_genus_hm.tiff", 
       bs_16s_top50_core_genus_hm, height = 5, width = 4, dpi = 400)
ggsave("figures/bs_16s_top50_core_genus_hm.svg", 
       bs_16s_top50_core_genus_hm, height = 5, width = 4, dpi = 400)



############################## Ordinations ####################################

# a. get our metadata
bs_lsu_nmds_meta <- bs_lsu_merged |>
  select(sample_id,berry_status)

# b. get our otus
bs_lsu_nmds_otu <- bs_lsu_merged |> 
  select(-c(2:8)) |>
  pivot_longer(-sample_id) |> 
  group_by(sample_id) |> 
  # summarise(N=sum(value)) |>  #this gives the total number of seqs in each
  # arrange(N) |> print(n=30) #sample and we arrange to get a value to rarefy
  mutate(N=sum(value)) |> #this shows the rarefied value
  filter(N >= 2500) |> #rarefying to keep all samples having 1400 seqs
  select(-N) |> 
  pivot_wider(names_from="name", values_from="value", values_fill=0)


# c. get our distance matrix
bs_lsu_nmds_dist <- bs_lsu_nmds_otu |> 
  column_to_rownames("sample_id") |> 
  avgdist(sample=2500)

# .d make ordination
set.seed(21) #required for the nmds to converge
bs_lsu_nmds_tbl <- metaMDS(bs_lsu_nmds_dist) |> 
  scores() |> 
  as_tibble(rownames="sample_id")

# e. merge nmds_meta and nmds_tbl
bs_lsu_nmds_4_plot <- inner_join(bs_lsu_nmds_meta, 
                                 bs_lsu_nmds_tbl, by="sample_id")

# test of significance before plotting
# ADONIS
# berry_status
adon_bs_lsu <- adonis2(bs_lsu_nmds_dist~bs_lsu_nmds_meta$berry_status,
                   permutations=9999)

str(adon_bs_lsu) #check where the p-value is
adon_bs_lsu$`Pr(>F)`[1]

# we prepare a centroid for our ellipses
bs_lsu_centroid <- bs_lsu_nmds_4_plot |> 
  group_by(berry_status) |> 
  summarise(NMDS1=mean(NMDS1),
            NMDS2=mean(NMDS2))

# plotting
bs_lsu_nmds_plot <- bs_lsu_nmds_4_plot|>
  ggplot(aes(x=NMDS1, y=NMDS2, color=berry_status)) +
  geom_point() +
  stat_ellipse(show.legend=F) +
  geom_point(data=bs_lsu_centroid, size=3,shape=21,
             color="#000000",aes(fill=berry_status),show.legend=F) +
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

ggsave("figures/bs_lsu_nmds_plot.pdf", 
       bs_lsu_nmds_plot, height = 4, width = 5, dpi = 300)
ggsave("figures/bs_lsu_nmds_plot.svg", 
       bs_lsu_nmds_plot, height = 4, width = 5, dpi = 300)
ggsave("figures/bs_lsu_nmds_plot.tiff", 
       bs_lsu_nmds_plot, height = 4, width = 5, dpi = 300)



################################ ITs Analysis ##################################

#metadata

bs_its_meta <- read_tsv("data/bs_its_gsr_metadata.tsv") |> 
  select('sample_id'='sample-id',berry_status) |>
  filter(sample_id != "fctrl") |>   #drop the control sample
  mutate(berry_status = case_when(berry_status == "healthy" ~ "H",
                                  berry_status == "Infected" ~ "I"))

bs_its_meta$berry_status <- as.factor(bs_its_meta$berry_status)

glimpse(bs_its_meta)

#evenness
bs_its_even <- read_tsv("data/bs_its_evenness.tsv") |>  
  rename("sample_id" = "...1", 'evenness'='pielou_evenness') #renames col 1

#faith
bs_its_faith <- read_tsv("data/bs_its_faith_pd.tsv") |>  
  rename("sample_id" = "#SampleID") #renames col 1

#observed_otus
bs_its_obs_otus <- read_tsv("data/bs_its_observed_otus.tsv") |>  
  rename("sample_id" = "...1", "obs_otus"="observed_features") #rename columns

#shannon
bs_its_shan <- read_tsv("data/bs_its_shannon.tsv") |>  
  rename("sample_id" = "...1", "shannon"="shannon_entropy")

#otu
bs_its_otus <- read_tsv("data/bs_its_rarefied_table_transp.txt", skip = 1) |> 
  rename("sample_id"="#OTU ID")

#ace
bs_its_ace <- read_tsv("data/bs_its_ace.tsv") |>  
  rename("sample_id" = "...1") #rename columns

#chao1
bs_its_chao1 <- read_tsv("data/bs_its_chao1.tsv") |>  
  rename("sample_id" = "...1")

#merge metadata, ace, chao1, evenness, faith, observed_otus, shannon and otus

bs_its_merged <- inner_join(bs_its_meta, bs_its_chao1, by='sample_id') %>%
  inner_join(., bs_its_ace, by='sample_id') %>%
  inner_join(., bs_its_even, by='sample_id') %>%
  inner_join(., bs_its_faith, by='sample_id') %>%
  inner_join(., bs_its_obs_otus, by='sample_id') %>% 
  inner_join(., bs_its_shan, by='sample_id') %>%
  inner_join(., bs_its_otus, by='sample_id')


# Normality Check

#Create 2x2 plot environment 
par(mfrow = c(2, 3))
hist(bs_its_merged$chao1, main="Chao1 diversity", xlab="", breaks=10, 
     border="#000000", col="#99DDFF", las=1)
hist(bs_its_merged$ace, main="Ace diversity", xlab="", breaks=10, 
     border="#000000", col="#BBCC33", las=1)
hist(bs_its_merged$shannon, main="Shannon diversity", xlab="", breaks=10, 
     border="#000000", col="#FFAABB", las=1)
hist(bs_its_merged$evenness, main="Pielou's evenness", xlab="", breaks=10, 
     border="#000000", col="#44BB99", las=1)
hist(bs_its_merged$faith_pd, main="Phylegenetic diversity", xlab="", 
     breaks=10, border="#000000", col="#DDDDDD", las=1)
hist(bs_its_merged$obs_otus, main="Community richness", xlab="", 
     breaks=10, border="#000000", col="#EEDD88", las=1)

# Test for normalcy

shapiro.test(bs_its_merged$ace)        
shapiro.test(bs_its_merged$chao1)      
shapiro.test(bs_its_merged$shannon)    # normally dist
shapiro.test(bs_its_merged$evenness)   # normally dist
shapiro.test(bs_its_merged$faith_pd)   # normally dist
shapiro.test(bs_its_merged$obs_otus)   # normally dist


# a. Anova test for normally distributed data

# i-ace
its.ace.bs.aov <- aov(ace ~ berry_status, data=bs_its_merged)
summary(its.ace.bs.aov)

# ii-chao1
its.chao1.bs.aov <- aov(chao1 ~ berry_status, data=bs_its_merged)
summary(its.chao1.bs.aov)

# iii-shannon
its.shan.bs.aov <- aov(shannon ~ berry_status, data=bs_its_merged)
summary(its.shan.bs.aov)

# iv-evenness
its.even.bs.aov <- aov(evenness ~ berry_status, data=bs_its_merged)
summary(its.even.bs.aov)

# v-Faith
its.fpd.bs.aov <- aov(faith_pd ~ berry_status, data=bs_its_merged)
summary(its.fpd.bs.aov)

# vi-Observed_otus
its.obs.bs.aov <- aov(obs_otus ~ berry_status, data=bs_its_merged)
summary(its.obs.bs.aov)


# plot alpha diversity
bs_its_alpha_div_plot <- bs_its_merged |> 
  select(-c(ace,evenness)) |>
  pivot_longer(cols = 3:6, names_to = 'indice', values_to = 'value') |> 
  mutate(indice=factor(indice, levels = c("chao1","shannon", 
                                          "obs_otus", "faith_pd"))) |>  
  select(berry_status, indice, value) |>  
  group_by(berry_status, indice) %>% 
  ggplot(., aes(x=berry_status, y=value)) + 
  geom_boxplot(aes(fill=berry_status),size=.5, outlier.shape = NA) +
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
    strip.text = element_text(size = 12, colour = "black", family = "sans"),
    panel.background = element_blank(),
    axis.title.y = element_text(size = 12, family = "sans"),
    axis.text = element_text(size = 11, color = "black"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(hjust = .5, vjust = 1, family = "sans"))

ggsave("figures/bs_its_alpha_div_plot_new.tiff", bs_its_alpha_div_plot, 
       height = 4.5, width = 5, dpi = 300)
ggsave("figures/bs_its_alpha_div_plot_new.pdf", bs_its_alpha_div_plot, 
       height = 4.5, width = 5, dpi = 300)
ggsave("figures/bs_its_alpha_div_plot_new.svg", bs_its_alpha_div_plot, 
       height = 4.5, width = 5, dpi = 300)

####################### import bs_taxonomy file for ITS ########################

bs_its_tax <- read.table("data/bs_its_taxonomy.tsv", 
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

#convert any NA values to empty spaces

bs_its_tax[is.na(bs_its_tax)] <- ""

#convert anything identified as "unidentified" into an empty space

bs_its_tax[bs_its_tax=="Unclassified"] <- "" 

#following code replaces those empty spaces with the Unclassified highest taxon
for (i in 1:nrow(bs_its_tax)) {
  if (bs_its_tax[i,7] != ""){
    bs_its_tax$Species[i] <- paste(bs_its_tax$Genus[i], sep = " ")
  }else if (bs_its_tax[i,2] == ""){
    Kingdom <- paste("Unclassified", bs_its_tax[i,1], sep = " ")
    bs_its_tax[i, 2:7] <- Kingdom
  }else if (bs_its_tax[i,3] == ""){
    Phylum <- paste("Unclassified", bs_its_tax[i,2], sep = " ")
    bs_its_tax[i, 3:7] <- Phylum
  }else if (bs_its_tax[i,4] == ""){
    Class <- paste("Unclassified", bs_its_tax[i,3], sep = " ")
    bs_its_tax[i, 4:7] <- Class
  }else if (bs_its_tax[i,5] == ""){
    Order <- paste("Unclassified", bs_its_tax[i,4], sep = " ")
    bs_its_tax[i, 5:7] <- Order 
  }else if (bs_its_tax[i,6] == ""){
    Family <- paste("Unclassified", bs_its_tax[i,5], sep = " ")
    bs_its_tax[i, 6:7] <- Family
  }else if (bs_its_tax[i,7] == ""){
    bs_its_tax$Species[i] <- paste("Unclassified", 
                                   bs_its_tax$Genus[i], sep = " ")
  }
}

#convert its_merged to long format and join to the cleaned taxonomy table

bs_its_merged %>% 
  select(-c(3:8)) %>% 
  pivot_longer(-c(sample_id, berry_status), names_to = "OTU", 
               values_to = "Counts") %>% 
  inner_join(., bs_its_tax, by=c('OTU'='Feature.ID')) %>% 
  write.table(., "results/bs_its_merged_long_new.tsv", sep="\t", 
              row.names=FALSE, quote = FALSE)


# preparing fungal phyla and genera plot by berry_status

bs_its_merged_long_ngs_cul <- read.table(
  "results/bs_its_merged_long_ngs_cul_edited_new.txt",head=T,sep="\t"
) |> 
  rename_all(tolower) |> 
  group_by(berry_status) |> 
  mutate(rel_abun=counts/sum(counts)) |> 
  ungroup() |> 
  filter(method != "cul") |> #drop culturable data
  select(-c(otu,counts,method,asv,kingdom,class,order,family,species)) |> 
  pivot_longer(cols = c("phylum","genus"), 
               names_to = "level", values_to = "taxon")

#Plotting

# a. fungal bar chart phyla berry_status

bs_its_phylum_rel_abun <- bs_its_merged_long_ngs_cul %>% 
  filter(level=="phylum") %>% 
  group_by(berry_status, taxon) |> 
  summarise(mean_rel_abun = mean(100*sum(rel_abun)),.groups="drop")

bs_its_phylum_plot <- bs_its_phylum_rel_abun %>% 
  ggplot(., aes(x=berry_status, y=mean_rel_abun, fill=taxon)) + 
  geom_col(show.legend = T) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(
    name=NULL,
    breaks = c("Ascomycota","Basidiomycota","Unclassified Fungi"),
    labels = c("*Ascomycota*","*Basidiomycota*","Unclassified *Fungi*"),
    values = c("#708090","#99DDFF","#000000","#b6d7a8","#8200ff","#6aa84f",
               "#545454","#f3f3f3","#bcbcbc","#f0e1d5","#a384b4","#eedd88",
               "#6e8b3d","#9ACD32","#F4A460","#ffffff","#77aadd","#696969",
               "#0000cd","#ffaabb","#B452CD","#ee8866","#bbbbbb","#ff0000",
               "#6495ed","#00fa9a","#00ced1","#aa4499","#d3d3d3","#8b008b",
               "#44bb99","#00ff00")
  ) +
  theme_bw() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5),
        #panel.border = element_blank(),
        axis.text = element_text(size = 8, family = "sans", color="#000000"),
        axis.title.y = element_text(size = 9, family = "sans"),
        legend.text = element_markdown(size=8),
        legend.key.size = unit(8, "pt"),
        legend.spacing.y = unit(.5, "cm"),
        axis.title.x = element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x=NULL, y="Mean Relative Abundance (%)",
       title = " ")
ggsave("figures/bs_its_phylum_plot_new.pdf", 
       bs_its_phylum_plot, width = 4, height = 3, dpi = 300)
ggsave("figures/bs_its_phylum_plot_new.svg", 
       bs_its_phylum_plot, width = 4, height = 3, dpi = 300)
ggsave("figures/bs_its_phylum_plot_new.tiff", 
       bs_its_phylum_plot, width = 4, height = 3, dpi = 300)


# b. fungi bar chart genera berry_status

bs_its_genus_rel_abun <- bs_its_merged_long_ngs_cul %>% 
  filter(level=="genus") %>% 
  group_by(berry_status, taxon) |> 
  summarise(mean_rel_abun = mean(100*sum(rel_abun)),.groups="drop")

#arrange rel_abun in increasing order
bs_its_genus_pool <- bs_its_genus_rel_abun %>% 
  group_by(taxon) %>% 
  summarise(pool=max(mean_rel_abun) < 4.31, .groups = "drop") #%>% 
  #arrange(desc(max))

bs_its_genus_plot <- inner_join(
  bs_its_genus_rel_abun, bs_its_genus_pool, by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |>  
  group_by(berry_status, taxon) %>% 
  ggplot(., aes(x=berry_status, y=mean_rel_abun, fill=taxon)) + 
  geom_col(show.legend = T) +
  scale_y_continuous(limits = c(0,100),expand = c(0,0)) +
  scale_fill_manual(
    name=NULL,
    breaks = c("Alternaria","Candida","Cladosporium","Filobasidium",
               "Neopestalotiopsis","Other","Pichia","Pseudopeyronellaea",
               "Sporobolomyces","Unclassified Ascomycota"),
    labels = c("*Alternaria*","*Candida*","*Cladosporium*","*Filobasidium*",
               "*Neopestalotiopsis*","Other","*Pichia*","*Pseudopeyronellaea*",
               "*Sporobolomyces*","Unclassified *Ascomycota*"),
    values = c("#708090","#99DDFF","#f3f3f3","#b6d7a8","#6aa84f","#8200ff",
               "#000000","#545454","#bcbcbc","#f0e1d5","#a384b4","#eedd88",
               "#6e8b3d","#9ACD32","#F4A460","#ffffff","#77aadd","#696969",
               "#0000cd","#ffaabb","#B452CD","#ee8866","#bbbbbb","#ff0000",
               "#6495ed","#00fa9a","#00ced1","#aa4499","#d3d3d3","#8b008b",
               "#44bb99","#00ff00")
  ) +
  theme_bw() +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust=0.5),
        #panel.border = element_blank(),
        axis.text = element_text(size = 8, family = "sans", color="#000000"),
        axis.title.y = element_text(size = 9, family = "sans"),
        legend.text = element_markdown(size=8),
        legend.key.size = unit(8, "pt"),
        legend.spacing.y = unit(.5, "cm"),
        axis.title.x = element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE)) +
  labs(x=NULL, y="Mean Relative Abundance (%)",
       title = " ")
ggsave("figures/bs_its_genus_plot_new.pdf", 
       bs_its_genus_plot, width = 4, height = 3, dpi = 300)
ggsave("figures/bs_its_genus_plot_new.svg", 
       bs_its_genus_plot, width = 4, height = 3, dpi = 300)
ggsave("figures/bs_its_genus_plot_new.tiff", 
       bs_its_genus_plot, width = 4, height = 3, dpi = 300)


# unstacked bar plots

bs_its_genus_plot_unstacked <- inner_join(
  bs_its_genus_rel_abun, bs_its_genus_pool, by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |> 
  mutate(taxon=factor(taxon),
         taxon=fct_reorder(taxon, mean_rel_abun, .desc=F)) |>
  group_by(berry_status, taxon) %>% 
  ggplot(., aes(x=taxon, y=mean_rel_abun, fill=berry_status)) + 
  geom_col(show.legend = T, position=position_dodge()) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(breaks=c("I","H"),
                    labels=c("Infected", "Healthy"),
                    values=c("#708090","#0000cd"))+
  scale_x_discrete(
    name=NULL,
    breaks = c("Alternaria","Candida","Cladosporium","Filobasidium",
               "Neopestalotiopsis","Other","Pichia","Pseudopeyronellaea",
               "Sporobolomyces","Unclassified Ascomycota",
               "Unclassified Dipodascaceae"),
    labels = c("*Alternaria*","*Candida*","*Cladosporium*","*Filobasidium*",
               "*Neopestalotiopsis*","Other","*Pichia*","*Pseudopeyronellaea*",
               "*Sporobolomyces*","Unclassified *Ascomycota*",
               "Unclassified *Dipodascaceae*")) +
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
ggsave("figures/bs_its_genus_plot_unstacked_new.pdf", 
       bs_its_genus_plot_unstacked, width = 4, height = 3, dpi = 300)
ggsave("figures/bs_its_genus_plot_unstacked_new.svg", 
       bs_its_genus_plot_unstacked, width = 4, height = 3, dpi = 300)
ggsave("figures/bs_its_genus_plot_unstacked_new.tiff", 
       bs_its_genus_plot_unstacked, width = 4, height = 3, dpi = 300)


############## import bs_ITS_taxonomy table for core biomes ###################

read_excel("data/bs_its_taxonomy_table_level-7_2.xlsx", 3) %>%
  pivot_longer(cols = c(
    "healthy...2","healthy...3","Infected...4","Infected...5","Infected...6",
    "Infected...7","Infected...8","Infected...9","Infected...10","Infected...11",
    "Infected...12","Infected...13","Infected...14","Infected...15",
    "Infected...16","Infected...17","Infected...18","Infected...19",
    "Infected...20","healthy...21","healthy...22","healthy...23","healthy...24",
    "healthy...25","healthy...26","healthy...27","healthy...28","healthy...29",
    "healthy...30","healthy...31","healthy...32","healthy...33","healthy...34",
    "healthy...35" ),
    names_to='berry_status', values_to='counts') %>%
  write.table(.,file="results/bs_its_tab_4_cor_biome_long_new.csv",
              col.names = NA, row.names = T, sep = ",", quote = F)


# venn diagrams from the combined ngs and cul data without otus

read_excel("results/bs_its_tab_4_cor_biome_long_edited_new.xlsx", 1) |> 
  separate(Species, into=c("Genus", NA), sep=" ") |> # get just the genus col
  rename_all(tolower) |> 
  filter(counts > 0) |> 
  filter(berry_status=="I") |>
  filter(method != "cul") |> 
  arrange(genus) |>  
  select(genus,berry_status) |> 
  distinct(genus,berry_status) %>%
  write.table(., file="results/bs_its_infected_distinct_genera_new.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)


read_excel("results/bs_its_tab_4_cor_biome_long_edited_new.xlsx", 1) |> 
  separate(Species, into=c("Genus", NA), sep=" ") |> # get just the genus col
  rename_all(tolower) |> 
  filter(counts > 0) |> 
  filter(berry_status=="H") |> 
  filter(method != "cul") |>
  arrange(genus) |>  
  select(genus,berry_status) |> 
  distinct(genus,berry_status) %>%
  write.table(., file="results/bs_its_helthy_distinct_genera_new.tsv",
              col.names = NA, row.names = T, sep = "\t", quote = F)


# import edited bs_its_healthy_infected_species_merged data for venn diagram

bs_its_healthy_infected_genus_merged <- read.table(
  "results/bs_its_healthy_infected_genera_merged_new.txt", 
  header=T, sep="\t") #|> view()

# creating list vector for venn diagram
bs_its_healthy_infected_genus_merged_venn_list = list(
  H = bs_its_healthy_infected_genus_merged$genus[1:38],
  I = bs_its_healthy_infected_genus_merged$genus[39:78]
)
# Ploting the venn diagram for shared mycobiomes across periods
bs_its_healthy_infected_genus_merged_venn_plot <- ggVennDiagram(
  bs_its_healthy_infected_genus_merged_venn_list, label_alpha = 0,
  category.names = c("H", "I")
) +
  scale_fill_gradient(low="#f2f2f2",high = "#ff0000") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.ticks.y.left = element_blank(),
    legend.position = "None"
  )
ggsave("figures/bs_its_healthy_infected_genus_merged_venn_plot_new.pdf",
       bs_its_healthy_infected_genus_merged_venn_plot, 
       height = 3, width = 3, dpi = 300)
ggsave("figures/bs_its_healthy_infected_genus_merged_venn_plot_new.tiff",
       bs_its_healthy_infected_genus_merged_venn_plot, 
       height = 3, width = 3, dpi = 300)
ggsave("figures/bs_its_healthy_infected_genus_merged_venn_plot_new.svg",
       bs_its_healthy_infected_genus_merged_venn_plot, 
       height = 3, width = 3, dpi = 300)


################################# Heatmap ####################################

bs_its_top50_core_genus <- read_excel(
  "results/bs_its_tab_4_cor_biome_long_edited_new.xlsx", 3) |> 
  select(-`Grand Total`) #|> 
#mutate_at(c('H','I'), ~replace_na(.,0)) # replaces all NAs with zeros

# melt converts data into long format but then we loose our rownames
bs_its_top50_core_genus_melt <- melt(bs_its_top50_core_genus)

#view first six rows
head(bs_its_top50_core_genus_melt)

# rescale values for all variables in melted data frame
bs_its_top50_core_genus_melt <- ddply(
  bs_its_top50_core_genus_melt, .(variable), 
  transform, dist = scale(value)) # rescale converts variables to a new range

# create heatmap using rescaled values
bs_its_top50_core_genus_hm <- ggplot(
  bs_its_top50_core_genus_melt, aes(x=variable, y=genus)) +
  geom_tile(aes(fill = dist), colour = "white") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradient(low = "#FFFFFF", high = "#000000","#ff0000",
                      expand = c(0,0)) +
  labs(x = "", y = "") +
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size=8, color="#000000"),
        axis.text.x = element_markdown(angle=0,hjust=.5,vjust=1,color="black")) +
  coord_fixed(ratio = 0.1) # resizes the boxes for each variable on the x-axis

ggsave("figures/bs_its_top50_core_genus_hm_new.pdf", 
       bs_its_top50_core_genus_hm, height = 5, width = 4, dpi = 400)
ggsave("figures/bs_its_top50_core_genus_hm_new.tiff", 
       bs_its_top50_core_genus_hm, height = 5, width = 4, dpi = 400)
ggsave("figures/bs_its_top50_core_genus_hm_new.svg", 
       bs_its_top50_core_genus_hm, height = 5, width = 4, dpi = 400)


############################## Ordinations ####################################

# a. get our metadata
bs_its_nmds_meta <- bs_its_merged |>
  select(sample_id,berry_status)

# b. get our otus
bs_its_nmds_otu <- bs_its_merged |> 
  select(-c(2:8)) |>
  pivot_longer(-sample_id) |> 
  group_by(sample_id) |> 
  # summarise(N=sum(value)) |>  #this gives the total number of seqs in each
  # arrange(N) |> print(n=30) #sample and we arrange to get a value to rarefy
  mutate(N=sum(value)) |> #this shows the rarefied value
  filter(N >= 25000) |> #rarefying to keep all samples having 1400 seqs
  select(-N) |> 
  pivot_wider(names_from="name", values_from="value", values_fill=0)


# c. get our distance matrix
bs_its_nmds_dist <- bs_its_nmds_otu |> 
  column_to_rownames("sample_id") |> 
  avgdist(sample=25000)

# .d make ordination
set.seed(21) #required for the nmds to converge
bs_its_nmds_tbl <- metaMDS(bs_its_nmds_dist) |> 
  scores() |> 
  as_tibble(rownames="sample_id")

# e. merge nmds_meta and nmds_tbl
bs_its_nmds_4_plot <- inner_join(bs_its_nmds_meta, 
                                 bs_its_nmds_tbl, by="sample_id")

# test of significance before plotting
# ADONIS
# berry_status
adon_bs <- adonis2(bs_its_nmds_dist~bs_its_nmds_meta$berry_status,
                    permutations=9999)

str(adon_bs) #check where the p-value is
adon_bs$`Pr(>F)`[1]

# we prepare a centroid for our ellipses
bs_its_centroid <- bs_its_nmds_4_plot |> 
  group_by(berry_status) |> 
  summarise(NMDS1=mean(NMDS1),
            NMDS2=mean(NMDS2))

# plotting
bs_its_nmds_plot <- bs_its_nmds_4_plot|>
  ggplot(aes(x=NMDS1, y=NMDS2, color=berry_status)) +
  geom_point() +
  stat_ellipse(show.legend=F) +
  geom_point(data=bs_its_centroid, size=3,shape=21,
             color="#000000",aes(fill=berry_status),show.legend=F) +
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

ggsave("figures/bs_its_nmds_plot_new.pdf", 
       bs_its_nmds_plot, height = 4, width = 5, dpi = 300)
ggsave("figures/bs_its_nmds_plot_new.svg", 
       bs_its_nmds_plot, height = 4, width = 5, dpi = 300)
ggsave("figures/bs_its_nmds_plot_new.tiff", 
       bs_its_nmds_plot, height = 4, width = 5, dpi = 300)

