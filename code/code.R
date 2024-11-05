# load required packages
source("F:/professdev/r_projects/m_biome_packages.R")
#choose packages to prefer
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)
conflicts_prefer(ggpubr::mutate)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(dplyr::desc)

# rename readme and license files to lowercase ----------------------------
file.rename(from='README.md', to = 'readme.md')

# import data -------------------------------------------------------------
#metadata
bs_meta <- read_tsv("data/metadata.tsv") |>
  select('sample_id'='sample-id','berry_status'='berry-status') |>
  filter(sample_id != "bctrl") |>
  mutate(berry_status = case_when(berry_status == "healthy" ~ "H",
                                  berry_status == "Infected" ~ "I"))
#convert to factor
bs_meta$berry_status <- as.factor(bs_meta$berry_status)
glimpse(bs_meta) #view

#evenness
bs_even <- read_tsv("data/evenness.tsv") |>  
  rename("sample_id" = "...1", 'evenness'='pielou_evenness') #renames col 1

#faith
bs_faith <- read_tsv("data/faith_pd.tsv") |>  
  rename("sample_id" = "#SampleID") #renames col 1

#observed_otus
bs_obs_otus <- read_tsv("data/observed_otus.tsv") |>  
  rename("sample_id" = "...1", "obs_otus"="observed_features") #rename columns

#shannon
bs_shan <- read_tsv("data/shannon.tsv") |>  
  rename("sample_id" = "...1", "shannon"="shannon_entropy")

#otu
#remove all OTUs that blasted as plant material from the taxonomy table
bs_otus <- read_tsv("data/rarefied_table.tsv", skip = 1) |> 
  rename("Group"="#OTU ID") |> #renamed Group to match with Lefse analysis
  filter(Group !='xyz') |> # xyz=plant, viral OTUs from blasting

# merge metadata and diversity indexes

bs_merged <- inner_join(bs_meta, bs_even, by='sample_id') %>%
  inner_join(., bs_faith, by='sample_id') %>%
  inner_join(., bs_obs_otus, by='sample_id') %>% 
  inner_join(., bs_shan, by='sample_id') %>%
  inner_join(., bs_otus, by=c('sample_id'='Group'))

############################### Alpha Diversity ################################

# Normality Check
# Shapiro test for normality ----------------------------------------------

shapiro.test(bs_merged$evenness)
shapiro.test(bs_merged$shannon)    # normally dist
shapiro.test(bs_merged$faith_pd)   # normally dist
shapiro.test(bs_merged$obs_otus)   # normally dist

# a. Anova test -----------------------------------------------------------

# i-evenness
kruskal.test(evenness ~ berry_status, data=bs_merged)

# ii-shannon
shan.bs.aov <- aov(shannon ~ berry_status, data=bs_merged)
summary(shan.bs.aov)

# iii-Faith
fd.bs.aov <- aov(faith_pd ~ berry_status, data=bs_merged)
summary(fd.bs.aov)

# iv-Observed_otus
obs.bs.aov <- aov(obs_otus ~ berry_status, data=bs_merged)
summary(obs.bs.aov)

# Alpha diversity plot ----------------------------------------------------

my_labs <- c("shannon" = "shannon",
             "evenness" = "evenness",
             "faith_pd" = "faith",
             "obs_otus" = "observed otus") # labels for y-axis of boxplots

# a. shannon, evenness, otus and faith by berry status
bs_alpha_div_plot <- bs_merged |> 
  pivot_longer(cols = 3:6, names_to = 'indice', values_to = 'value') |> 
  select(1:2,indice,value, everything()) |> 
  mutate(indice=factor(indice, levels = c("shannon","evenness", 
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
  scale_x_discrete(breaks=c('H','I'),
                   labels=c('Healthy','Infected')) +
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

ggsave("bs_alpha_div_plot.tiff",bs_alpha_div_plot,height=4.5,width=5,dpi=300)


####################### import bs_taxonomy file for 16S ########################

bs_tax <- read.table("data/taxonomy.tsv", header = T, sep = "\t") |> 
  select('OTU'='Feature.ID', Taxon) |> 
  mutate(Taxon = str_replace_all(Taxon, "d__", ""),#replace k_ with ,
         Taxon = str_replace_all(Taxon, "; p__", ","), #replace p_ with ,
         Taxon = str_replace_all(Taxon, "; c__", ","), #replace c_ with ,
         Taxon = str_replace_all(Taxon, "; o__", ","), #replace o_ with ,
         Taxon = str_replace_all(Taxon, "; f__", ","), #replace f_ with ,
         Taxon = str_replace_all(Taxon, "; g__", ","), #replace g_ with ,
         Taxon = str_replace_all(Taxon, "; s__", ",")) |>   #replace s_ with ,
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", 
                    "Family", "Genus", "Species"), ",")|> 
  mutate_all(na_if,"") #replace all empty strings with NA values

#because of the warning message on NAs, we need to address this
#convert any NA values to empty spaces
bs_tax[is.na(bs_tax)] <- ""
#convert anything identified as "unidentified" into an empty space
bs_tax[bs_tax=="Unclassified"] <- "" 

#following code replaces those empty spaces with the Unclassified highest taxon
for (i in 1:nrow(bs_tax)) {
  if (bs_tax[i,7] != ""){
    bs_tax$Species[i] <- paste(bs_tax$Genus[i], sep = " ")
  }else if (bs_tax[i,2] == ""){
    Kingdom <- paste("Unclassified", bs_tax[i,1], sep = " ")
    bs_tax[i, 2:7] <- Kingdom
  }else if (bs_tax[i,3] == ""){
    Phylum <- paste("Unclassified", bs_tax[i,2], sep = " ")
    bs_tax[i, 3:7] <- Phylum
  }else if (bs_tax[i,4] == ""){
    Class <- paste("Unclassified", bs_tax[i,3], sep = " ")
    bs_tax[i, 4:7] <- Class
  }else if (bs_tax[i,5] == ""){
    Order <- paste("Unclassified", bs_tax[i,4], sep = " ")
    bs_tax[i, 5:7] <- Order 
  }else if (bs_tax[i,6] == ""){
    Family <- paste("Unclassified", bs_tax[i,5], sep = " ")
    bs_tax[i, 6:7] <- Family
  }else if (bs_tax[i,7] == ""){
    bs_tax$Species[i] <- paste("Unclassified", 
                                   bs_tax$Genus[i], sep = " ")
  }
}

# edit to remove underscores and fill taxa names
# import edited taxonomy table
bs_tax_edited <- read_excel("bs_tax.xlsx",sheet = 2) |> 
  rename_all(tolower)

#create a vector of otus to match the rows in bs_lsu_tax_edited
bs_name_tmp <-sprintf("otu%s",seq(1:nrow(bs_tax_edited)))
bs_tax_edited$asv <- paste(bs_name_tmp,bs_tax_edited$asv, sep="")
bs_tax_edited$asv[1:5] #look at first 5 rows

bs_tax_edited <- bs_tax_edited |>
  mutate(method='ngs') |> 
  select(otu,asv,method,everything())

#convert lsu_merged to long format and join to the cleaned taxonomy table

bs_meta_otus_tax_long <- bs_merged |> 
  select(-c(3:6)) |>  
  pivot_longer(-c(sample_id, berry_status),names_to="otu",values_to="count") |>  
  inner_join(bs_tax_edited, by='otu')


########################### Beta diversity ####################################

# Test for significance using Adonis2 -------------------------------------

# a. prepare the metadata -------------------------------------------------
bs_nmds_meta <- bs_merged |>
  filter(sample_id !='cbt33') |> #lost during rarefaction
  select(sample_id,berry_status)

# b. prepare the otus -----------------------------------------------------
bs_nmds_otu <- read_tsv("bs_meta_otus_tax_long.tsv") |> 
  select(sample_id,otu,count) |>
  group_by(sample_id) |> 
  mutate(N=sum(count)) |> # the last 4 samples have 764,763,762 & 686 seqs
  filter(N >= n) |> # rarefy at n to keep max samples
  select(-N) |> 
  pivot_wider(names_from="otu", values_from="count", values_fill=0)

# c. get our distance matrix
bs_nmds_dist <- bs_nmds_otu |> 
  column_to_rownames("sample_id") |> 
  avgdist(sample=n) # n=numeric

# d. Adonis (PERMANOVA) for Healthy vs Infected ---------------------------
adon_bs <- adonis2(bs_nmds_dist ~ berry_status, permutations=9999, 
                       data=bs_nmds_meta)

# e. make ordination ------------------------------------------------------
set.seed(21) #required for the nmds to converge
bs_nmds_tbl <- metaMDS(bs_nmds_dist) |> 
  scores() |> 
  as_tibble(rownames="sample_id")

# f. merge nmds_meta and nmds_tbl -----------------------------------------
bs_nmds_4_plot <- inner_join(bs_nmds_meta,
                             bs_nmds_tbl, by="sample_id") |> 
  mutate(berry_status=case_when(berry_status=='H'~'Healthy',
                                berry_status=='I'~'Infected'))

# g. prepare centroids for the ellipses -----------------------------------
bs_centroid <- bs_nmds_4_plot |> 
  group_by(berry_status) |> 
  summarise(NMDS1=mean(NMDS1),
            NMDS2=mean(NMDS2))

# segment for berry_status ------------------------------------------------
bs_segment <- data.frame(
  bs_nmds_4_plot |> 
    select(berry_status,NMDS1,NMDS2) |>
    mutate(xend=if_else(berry_status=='Healthy',-0.195,0.160),
           yend=if_else(berry_status=='Healthy',0.0355,-0.0293))
)

# plotting for berry_status -----------------------------------------------
bs_nmds_plot <- bs_nmds_4_plot|>
  ggplot(aes(x=NMDS1, y=NMDS2, color=berry_status)) +
  geom_point(aes(colour=berry_status),shape=8,) +
  geom_point(data=bs_centroid, size=2,shape='.',
             color="#000",aes(fill=berry_status),show.legend=F) +
  stat_ellipse(show.legend=F) +
  scale_color_manual(values=c("#0000cd","#ff0000")) +
  geom_segment(bs_segment,mapping=aes(x=NMDS1,y=NMDS2,
                                          xend=xend,yend=yend)) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.background = element_blank(),
        legend.text= element_text(size = 10, family = 'sans'),
        legend.title = element_blank(),
        legend.key.size = unit(.5, units = "cm"),
        axis.text = element_text(colour = 'black'))

ggsave("bs_nmds_plot.pdf", bs_nmds_plot, height = 4, width = 5, dpi = 300)

# preparing bacteria phyla plot by berry_status

bs_merged_long_ngs_cul <- read_tsv("bs_meta_otus_tax_long_ngs_cul3.txt") |> 
  rename_all(tolower) |> 
  filter(method != "cul") |> #drop culturable data
  group_by(berry_status) |> 
  mutate(rel_abun=count/sum(count),
         .before='asv') |> #ads before asv col. (.after will ad after)
  ungroup() |>
  select(-c(otu,count,method,asv,kingdom,class,family,order,species)) |> 
  pivot_longer(cols = c("phylum","genus"), 
               names_to = "level", values_to = "taxon")

# a. bacteria genera by berry_status bar chart  ---------------------------
# edit for other taxonomic ranks

bs_genus_rel_abun <- bs_merged_long_ngs_cul %>% 
  filter(level=="genus") %>% 
  group_by(berry_status, taxon) |> 
  summarise(mean_rel_abun = mean(100*sum(rel_abun)),.groups="drop")

#arrange rel_abun in increasing order
bs_genus_pool <- bs_genus_rel_abun %>% 
  group_by(taxon) %>% 
  summarise(pool=max(mean_rel_abun) < n, .groups = "drop") # n=numeric

bs_genus_plot <- inner_join(
  bs_genus_rel_abun, bs_genus_pool, by="taxon") |>  
  mutate(taxon=if_else(pool, "Other", taxon)) |>  
  group_by(berry_status, taxon) %>% 
  ggplot(., aes(x=berry_status, y=mean_rel_abun, fill=taxon)) + 
  geom_col(show.legend = T) +
  scale_y_continuous(limits = c(0,100),expand = c(0,0)) +
  scale_x_discrete(breaks=c('H','I'),labels=c('Healthy','Infected')) +
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
    values = c("#708090","#77aadd","#f3f3f3","#0000cd","#d3d3d3","#000000",
               "#bbbbbb","#545454","#bcbcbc","#f0e1d5","#a384b4","#eedd88",
               "#6e8b3d","#9ACD32","#F4A460","#ffffff","#99DDFF","#696969",
               "#b6d7a8","#ffaabb","#B452CD","#ee8866","#8200ff","#ff0000",
               "#6495ed","#00fa9a","#00ced1","#aa4499","#6aa84f","#8b008b",
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
ggsave("bs_genus_plot.pdf",bs_genus_plot, width = 4, height = 3, dpi = 300)

# Venn Diagram for berry_status -------------------------------------------
# Healthy berries genera
read_tsv(
  "results/bs_meta_otus_tax_long_ngs_cul.txt") |> 
  rename_all(tolower) |> 
  filter(count > 0) |> 
  filter(method!='cul') |> 
  filter(berry_status=="I") |>
  arrange(genus) |>  
  distinct(genus, .keep_all=T) |>  #.keep_all keeps all other columns
  write_tsv(file="results/bs_infected_distinct_genera.tsv")

# Infected berries genera
read_tsv(
  "results/bs_meta_otus_tax_long_ngs_cul.txt") |> 
  rename_all(tolower) |> 
  filter(count > 0) |> 
  filter(method!='cul') |> 
  filter(berry_status=="H") |>
  arrange(genus) |>  
  distinct(genus, .keep_all=T) |>  
  write_tsv(file="results/bs_helthy_distinct_genera.tsv")

# import edited bs_healthy_infected_species_merged data -------------------

bs_genera_merged <- read_tsv("results/bs_genera_merged.txt") #|> view()

# creating list vector for venn diagram
bs_genera_venn_list = list(
  H = bs_genera_merged$genus[a:b], # a,b,c,d are numerics
  I = bs_genera_merged$genus[c:d]
)
# Ploting the venn diagram for shared mycobiomes across periods
bs_genera_venn_plot <- ggVennDiagram(bs_genera_venn_list, label_alpha = 0,
                                     category.names = c("Healthy", "Infected")
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
ggsave("bs_venn_plot.pdf",bs_genera_venn_plot,height=3,width=3,dpi=300)

# LefSe


# Filter data to use for LefSe --------------------------------------------
bs_otu_meta <- bs_merged |> 
  select(-c(3:6)) 

# Filter and save otu or shared file for lefse in mothur ------------------
bs_otu_meta |>
  select(-berry_status) |> 
  select('Group'='sample_id',everything()) |> #view() # 337 total cols
  mutate(label=0.04, #create a label column as mothur needs this for lefse
         numOtus=336) |> #mothur also needs this column for lefse
  select(label,Group,numOtus,everything()) |> 
  write_tsv('bs.shared')

# Get meta data -----------------------------------------------------------
bs_otu_meta |> 
  select('Group'='sample_id',berry_status) |> 
  write_tsv('bs.design')

# wrap the mothur code to execute from rstudio ----------------------------
command <- 'mothur/mothur.exe "#lefse(shared=bs.shared,design=bs.design)"'

# execute the above command to see mothur results -------------------------
system(command)

# read in the summary of results to build plot ----------------------------
asv_lefse_plot <- read_tsv('results/bs.0.04.lefse_summary') |> 
  drop_na(LDA) %>%
  inner_join(.,bs_lsu_tax_edited, by=c('OTU'='otu')) |>
  mutate(LDA=ifelse(Class=='I', -1*LDA, LDA), #makes Nnet to be negative
         OTU=fct_reorder(OTU, LDA))%>% #reorder OTUs by LDA
  ggplot(aes(x=LDA, y=asv, fill=Class)) +
  geom_col(width=.7) +
  geom_vline(xintercept=0) +
  labs(x='LDA Score (log10)',
       y=NULL) +
  scale_fill_manual(name=NULL, breaks=c('H','I'),
                    labels=c('Healthy','Infected'),
                    values=c('#0000FF','#000','#228b22','#ff0000')) +
  theme_bw() +
  theme(
    legend.position=c(.9,.94),
    legend.background=element_blank(),
    axis.title=element_text(size=12, family="sans", color="#000"),
    axis.text=element_text(size=10, family="sans", color="#000")
  )
ggsave("asv_lefse_plot.pdf", asv_lefse_plot, height = 4, width = 5, dpi = 300)


################################# Heatmap ####################################

# Culture-independent data ------------------------------------------------
core_genus <- read_excel("cor_biome_long_edited.xlsx", 3) |> 
  select(-`Grand Total`) #|> 
#mutate_at(c('H','I'), ~replace_na(.,0)) # replaces all NAs with zeros

# melt converts data into long format but then we loose our rownames
core_genus_long <- core_genus |> 
  pivot_longer(-genus, names_to='status', values_to='count')

#view first six rows
head(core_genus_long)

# rescale values for all variables in melted data frame
core_genus_long <- ddply(
  core_genus_long, .(status), 
  transform, dist = scale(count)) # rescale converts variables to a new range

# create heatmap using rescaled values
core_genus_hm <- ggplot(
  core_genus_long, aes(x=status, y=genus)) +
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

ggsave("core_genus_hm.pdf", core_genus_hm, height = 6, width = 4, dpi = 400)

# Culture-dependent data --------------------------------------------------

# Read in culturable data of abundances  ----------------------------------
core_genus_cul <- read_excel("cor_biome_long_edited3.xlsx", 5) |> 
  select(-`Grand Total`)

# Concert to long format --------------------------------------------------
core_genus_cul_long <- core_genus_cul |> 
  pivot_longer(-genus, names_to='status', values_to='count')

# Rescale data ------------------------------------------------------------
core_genus_cul_long <- ddply(
  core_genus_cul_long, .(status), 
  transform, dist = scale(count))

# Convert back to wide format and export to remove NAs --------------------
core_genus_cul_long |> 
  pivot_wider(names_from='status', values_from='dist') |> 
  write.csv("core_genus_cul_long.csv")

# Import data edited by removing introduced NAs ---------------------------
core_genus_cul_edited <- read_csv(
  "core_genus_cul_long.csv") |> 
  select(-count) %>% 
  remove_rownames %>% 
  column_to_rownames(var="genus")

# Create a matrix of distances --------------------------------------------
mat <- as.matrix(core_genus_cul_edited)

# Build heatmap  ----------------------------------------------------------
pheatmap(mat) # default in complete method

# More methods single,complete,median,centroid,mcquitty,ward.D2 -----------
hm_plot <- pheatmap(mat,
                    cluster_rows=T, cluster_cols=T,
                    clustering_method='average',
                    display_numbers=T,
                    fontsize_number=5,
                    number_color='#dcdcdc',
                    col=brewer.pal(8,'GnBu'), # others BuGn,YlGn etc
                    border_color='#dedede',
                    legend_breaks=c(-4,0,4),
                    legend_labels=c('low','medium','high'),
                    cutree_cols=2) 

# Save the plot -----------------------------------------------------------
pdf(paste0('hm_plot.pdf'),height=5, width=5)
hm_plot
dev.off()

############################## Rarefaction ####################################

# import unrarefied otus --------------------------------------------------
unrarefied_otus <- read_tsv(
  "data/unrarefied_table.tsv",skip = 1) |> # make sure to transpose table
  rename("sample_id"="#OTU ID")

# merge unrarefied otus and metadata --------------------------------------
unrarefied_otu_meta <- inner_join(meta,
                                  unrarefied_otus,
                                  by='sample_id')

# rarefy  -----------------------------------------------------------------
rarefied <- unrarefied_otu_meta |> 
  select(-berry_status) |> 
  as.data.frame() # convert because tibbles don't allow row names

# make sample_ids rownames ------------------------------------------------
rownames(rarefied) <- rarefied$sample_id
rarefied <- rarefied[,-1]# drop first column

str(rarefied)

# rarefy using the same threshold as for ordinations ----------------------
rarefy(rarefied, n) # n is a numeric 

# # bind elements of the list above together ------------------------------
rarecurve_data <- rarecurve(rarefied, step = 100)

# plot the rarefaction curve ----------------------------------------------
#map_dfr is a purr function
rarecurve_plot <- map_dfr(rarecurve_data, bind_rows) |> 
  bind_cols(sample_id = rownames(rarefied)) |> 
  pivot_longer(-sample_id) |>
  drop_na() |> 
  mutate(n_seqs=as.numeric(str_replace(name,"N",""))) |> 
  select(-name) |> 
  ggplot(aes(x=n_seqs, y=value, group=sample_id)) +
  geom_vline(xintercept=760, color="#cfcfcf") +
  geom_line() +
  scale_x_continuous(breaks=seq(0,1010,100), 
                     limits=c(0,1010),
                     expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0,100,20),
                     limits=c(0,100),
                     expand=c(0,0)) +
  labs(x = 'Sample size', y = 'Number of species') + 
  theme_classic() +
  theme(
    axis.title = element_text(size=12,family="sans",color="#000000"),
    axis.text = element_text(size=11,family="sans",color="#000000")
  )
ggsave("rarecurve_plot.pdf",rarecurve_plot, width = 5, height = 4, dpi = 300)