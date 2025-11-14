rm(list = ls())
setwd("E:/sunyingying/1111LS-GNSHF-project/1111data/01.IGC_humanswiss/01.GNHSF/09.figure/FigureS2_3_count_20250908")
library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyverse)
library(forcats)
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
#install.packages("ggsci")
library(ggsci)
#install.packages("paletteer")
library(paletteer)
#display.brewer.all()
mycolors<-c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
            "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"
)


# species_phylum_pepcount
sp_pepcount <- as.data.frame(read_tsv("pepcount_species_phylum_shortname.txt"))
sp_pepcount_head15 <- sp_pepcount[1:15,]

#mycolors <- brewer.pal(5, "Set1")
#ggplot(sp_pepcount_head40, aes(species, count, fill=phylum))+geom_bar(stat="identity", position="dodge")+coord_flip()
p1 <- sp_pepcount_head15 %>%
  mutate(species = fct_reorder(species, `Peptide count`)) %>%
  ggplot(aes(x=species, y=`Peptide count`, fill=phylum)) +
  geom_bar(stat="identity", width=.4) +
  coord_flip() +
  xlab("") +
  theme_classic() + scale_fill_manual(values=mycolors) + 
  theme(axis.ticks = element_blank(), axis.title.x = element_text(size = 25,face = "bold"), axis.text.x  = element_text(size=40), axis.text.y  = element_text(size=20))+
  theme(legend.title = element_text(size=30, face = "bold"), legend.text = element_text(size=30, face = "bold"))+
  geom_col(width = 0.8, position = position_dodge(width=0.8))

p1

ggsave("micro_species_pep_count_top15_mycolors_big.png", p1, width = 20, height = 10, dpi=600)
ggsave("micro_species_pep_count_top15_mycolors_big.pdf", p1, width = 20, height = 10)


# genus_phylum_pepcount
sp_pepcount <- as.data.frame(read_tsv("pepcount_genus_phylum_shortname.txt"))
sp_pepcount_head15 <- sp_pepcount[1:15,]

#mycolors <- brewer.pal(5, "Set1")
#ggplot(sp_pepcount_head40, aes(genus, count, fill=phylum))+geom_bar(stat="identity", position="dodge")+coord_flip()
p1 <- sp_pepcount_head15 %>%
  mutate(genus = fct_reorder(genus, `Peptide count`)) %>%
  ggplot(aes(x=genus, y=`Peptide count`, fill=phylum)) +
  geom_bar(stat="identity", width=.4) +
  coord_flip() +
  xlab("") +
  theme_classic() + scale_fill_manual(values=mycolors) + 
  theme(axis.ticks = element_blank(), axis.title.x = element_text(size = 25,face = "bold"), axis.text.x  = element_text(size=40), axis.text.y  = element_text(size=20))+
  theme(legend.title = element_text(size=30, face = "bold"), legend.text = element_text(size=30, face = "bold"))+
  geom_col(width = 0.8, position = position_dodge(width=0.8))

p1

ggsave("micro_genus_pep_count_top15_mycolors_big.jpeg", p1, width = 20, height = 10, dpi=600)
ggsave("micro_genus_pep_count_top15_mycolors_big.pdf", p1, width = 20, height = 10)


## phylum pep count
pp_pepcount <- as.data.frame(read_tsv("pepcount_phylum_shortname.txt"))
pp_pepcount_head4 <- pp_pepcount[1:4,]
pp <- pp_pepcount_head4 %>%
  mutate(taxon = fct_reorder(taxon, `count`)) %>%
  ggplot(aes(x=taxon, y=`count`)) +
  geom_bar(stat="identity", width=.4) +
  coord_flip() +
  xlab("") +
  theme_classic() + scale_fill_manual(values=mycolors) + 
  theme(axis.ticks = element_blank(), axis.title.x = element_text(size = 25,face = "bold"), axis.text.x  = element_text(size=40), axis.text.y  = element_text(size=20))+
geom_col(width = 0.8, position = position_dodge(width=0.8))

pp
ggsave("micro_phylum_pep_count_top4_big.jpeg", pp, width = 20, height = 10, dpi=600)
ggsave("micro_phylum_pep_count_top4_big.pdf", pp, width = 20, height = 10)


# cog_procount
## colornum = 11
mycolors<-c(brewer.pal(8, "Set1"), brewer.pal(9, "Paired")[1], brewer.pal(9, "Paired")[9], brewer.pal(8, "Pastel2")[8])

cog_procount <- as.data.frame(read_tsv("micro_procount_cog_des.txt"))
cog_procount <- as.data.frame(read_tsv("human_procount_cog_des.txt"))
cog_procount_head15 <- cog_procount[1:15,]


p2 <- cog_procount_head15 %>%
  mutate(cog = fct_reorder(COG, `Protein group count`)) %>%
  ggplot(aes(x=cog, y=`Protein group count`, fill=COGcat)) +
  geom_bar(stat="identity", width=.4) +
  coord_flip() +
  xlab("") +
  theme_classic() + scale_fill_manual(values=mycolors)+
  theme(axis.ticks = element_blank(), axis.title.x = element_text(size = 25,face = "bold"), axis.text.x  = element_text(size=20), axis.text.y  = element_text(size=15))+
  theme(legend.title = element_text(size=20, face = "bold"), legend.text = element_text(size=20, face = "bold"))+
  geom_col(width = 0.8, position = position_dodge(width=0.8))


p2

ggsave("micro_cog_procount_top15_mycolors.png", p2, width = 20, height = 10, dpi=600)
ggsave("micro_cog_procount_top15_mycolors.pdf", p2, width = 20, height = 10)


ggsave("human_cog_procount_top15_mycolors.png", p2, width = 20, height = 10, dpi=600)
ggsave("human_cog_procount_top15_mycolors.pdf", p2, width = 20, height = 10)

# cogcat_procount(cat2pro)
cat_procount <- as.data.frame(read_tsv("micro_procount_cogcat_des.txt"))
cat_procount <- as.data.frame(read_tsv("human_procount_cogcat_des.txt"))
cat_procount <- cat_procount[1:15, ]

p3 <- cat_procount %>%
  mutate(cogcat = fct_reorder(COGcat_des, `Protein group count`)) %>%
  ggplot(aes(x=cogcat, y=`Protein group count`)) +
  geom_bar(stat="identity", width=.4, alpha=.6) +
  coord_flip() +
  xlab("") +
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.title.x = element_text(size = 30,face = "bold"), axis.text.x  = element_text(size=30), axis.text.y  = element_text(size=25))+
  geom_col(width = 0.8, position = position_dodge(width=0.8))


p3
ggsave("micro_cogcat_procount_top15.jpeg", p3, width = 20, height = 10, dpi=600)
ggsave("micro_cogcat_procount_top15.pdf", p3, width = 20, height = 10)

ggsave("human_cogcat_procount_top15.jpeg", p3, width = 20, height = 10, dpi=600)
ggsave("human_cogcat_procount_top15.pdf", p3, width = 20, height = 10)

# kegg_pro
kegg_procount <- as.data.frame(read_tsv("micro_procount_kegg_des_modify.txt"))
kegg_procount <- as.data.frame(read_tsv("human_procount_kegg_des.txt"))
kegg_procount_head15 <- kegg_procount[1:15,]

p4 <- kegg_procount_head15 %>%
  mutate(`KEGG num` = fct_reorder(`KEGG`, `Protein group count`)) %>%
  ggplot(aes(x=`KEGG num`, y=`Protein group count`, fill=`KEGGcat`)) +
  geom_bar(stat="identity", width=.4) +
  coord_flip() +
  xlab("") +
  theme_classic() +scale_fill_manual(values=mycolors)+
  theme(axis.ticks = element_blank(), axis.title.x = element_text(size = 15,face = "bold"), axis.text.x  = element_text(size=20), axis.text.y  = element_text(size=20))+
  theme(legend.title = element_text(size=20, face = "bold"), legend.text = element_text(size=20, face = "bold"))+
  geom_col(width = 0.8, position = position_dodge(width=0.8))

p4

ggsave("micro_KEGG_procount_top15_mycolors.png", p4, width = 20, height = 10, dpi=600)
ggsave("micro_KEGG_procount_top15_mycolors.pdf", p4, width = 20, height = 10)

ggsave("human_KEGG_procount_top15_mycolors.png", p4, width = 20, height = 10, dpi=600)
ggsave("human_KEGG_procount_top15_mycolors.pdf", p4, width = 20, height = 10)

# kegg_cat_pro
keggcat_procount <- as.data.frame(read_tsv("micro_procount_keggcat_des.txt"))
keggcat_procount <- as.data.frame(read_tsv("human_procount_keggcat_des.txt"))
keggcat_procount <- keggcat_procount[1:15,]

p5 <- keggcat_procount %>%
  mutate(`KEGGcat` = fct_reorder(`KEGGcat`, `Protein group count`)) %>%
  ggplot(aes(x=`KEGGcat`, y=`Protein group count`)) +
  geom_bar(stat="identity", width=.4, alpha=.6) +
  coord_flip() +
  xlab("") +
  theme_classic()+
  theme(axis.ticks = element_blank(), axis.title.x = element_text(size = 30,face = "bold"), axis.text.x  = element_text(size=30), axis.text.y  = element_text(size=30))+
  geom_col(width = 0.8, position = position_dodge(width=0.8))

p5
ggsave("micro_KEGGcat_procount_top15.jpeg", p5, width = 20, height = 10, dpi=600)
ggsave("micro_KEGGcat_procount_top15.pdf", p5, width = 20, height = 10)

ggsave("human_KEGGcat_procount_top15.jpeg", p5, width = 20, height = 10, dpi=600)
ggsave("human_KEGGcat_procount_top15.pdf", p5, width = 20, height = 10)


