rm(list = ls())
library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyverse)
library(forcats)
library(RColorBrewer)
#install.packages("ggsci")
library(ggsci)
#install.packages("paletteer")
library(paletteer)

############################################################################################
anno <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microbiome_all.txt",sep="\t")
anno[anno==""] <- NA
anno1 <- unique(anno[,c(2,5,11)])
micro <- data.frame(data.table::fread("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_microprotein_sample.tsv"))
micro$Protein.Group <- unlist(lapply(micro$Protein.Group,function(x){
  strsplit(x," ")[[1]][1]
} ))
micro <- merge(anno1,micro,by.x = "Protein",by.y="Protein.Group")
map <- apply(micro[,4:ncol(micro)],2,function(x) { s.matrix <- micro[which(!is.na(x)),1:3]
  cog.map <- length(which(!is.na(s.matrix$COG)))/nrow(s.matrix)
  kegg.map <- length(which(!is.na(s.matrix$KEGG)))/nrow(s.matrix)
  return(c(length(which(!is.na(s.matrix$COG))),cog.map,length(which(!is.na(s.matrix$KEGG))),kegg.map,nrow(s.matrix)))
  })
rownames(map) <- c("cog.map","cog.prop","kegg.map","kegg.prop","total")

library(ggplot2)
library(ggthemes)
## map_cog
colors12<-c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
            "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"
)
mycolors<-c("#CC6677", brewer.pal(11,"Spectral")[10])
map_cog <- data.frame("sample"=colnames(map),"mapped"=map[2,],"unmapped"=1-map[2,],"arrange"=map[2,])
map_cog2 <- reshape2::melt(map_cog,id.vars = c("sample","arrange"))
map_cog2$condition = factor(map_cog2$variable, levels = c("mapped","unmapped"))
map.median <- round(median(map_cog$mapped),4)
unmap.median <- round(median(map_cog$unmapped),4)

p2 <- map_cog2 %>%
  mutate(sample = fct_reorder(sample, arrange, .desc = TRUE)) %>%
  ggplot(aes(x=sample, y=value, fill=condition)) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() + scale_fill_manual(values=mycolors) +xlab("Samples")+ylab("Proportion (%)")+scale_x_discrete(labels=NULL)+ggtitle("Function - COG")+
  theme_base()+
  theme(axis.ticks = element_blank())+
  geom_text(aes(label = paste0("Mapped Median ",100*map.median,"%"),x=1000,y=0.5),color="white",size = 5)+
  geom_text(aes(label =paste0("Unmapped Median ",100*unmap.median,"%"),x=1000,y=0.03),color="white",size=5)

p2
ggsave("map_cog.pdf", p2, width = 20, height = 10)

## map_kegg
map_kegg <- data.frame("sample"=colnames(map),"mapped"=map[4,],"unmapped"=1-map[4,],"arrange"=map[4,])
map_kegg2 <- reshape2::melt(map_kegg,id.vars = c("sample","arrange"))
map_kegg2$condition = factor(map_kegg2$variable, levels = c("mapped","unmapped"))
map.median <- round(median(map_kegg$mapped),4)
unmap.median <- round(median(map_kegg$unmapped),4)

p2 <- map_kegg2 %>%
  mutate(sample = fct_reorder(sample, arrange, .desc = TRUE)) %>%
  ggplot(aes(x=sample, y=value, fill=condition)) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() + scale_fill_manual(values=mycolors) +xlab("Samples")+ylab("Proportion (%)")+scale_x_discrete(labels=NULL)+ggtitle("Function - KO")+
  theme_base()+
  theme(axis.ticks = element_blank())+
  geom_text(aes(label = paste0("Mapped Median ",100*map.median,"%"),x=1000,y=0.5),color="white",size = 5)+
  geom_text(aes(label =paste0("Unmapped Median ",100*unmap.median,"%"),x=1000,y=0.1),color="white",size=5)

p2
ggsave("map_kegg.pdf", p2, width = 20, height = 10)



###########taxon######
library(openxlsx)
anno <- read.xlsx("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microprot_all_arrange20250916.xlsx")
anno[anno==""] <- NA
anno2 <- unique(anno[,c(1,14,15)])
micro <- data.frame(data.table::fread("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_microprotein_sample.tsv"))
micro$Protein.Group <- unlist(lapply(micro$Protein.Group,function(x){
  strsplit(x," ")[[1]][1]
} ))
micro <- merge(anno2,micro,by.x = "Protein",by.y="Protein.Group")
map <- apply(micro[,4:ncol(micro)],2,function(x) { s.matrix <- micro[which(!is.na(x)),1:3]
superkingdom <- length(which(!is.na(s.matrix$Taxon_rank_unipept)))
phylum <- length(which(s.matrix$Taxon_rank_unipept %in% c("genus","species","strain","class","order","family","phylum")))
genus <- length(which(s.matrix$Taxon_rank_unipept %in% c("genus","species","strain")))
species <- length(which(s.matrix$Taxon_rank_unipept %in% c("species","strain")))
return(c(superkingdom,phylum,genus,species,nrow(s.matrix),(nrow(s.matrix)-superkingdom)))
})
rownames(map) <- c("superkingdom","phylum","genus","species","total","unmapped")

map_taxon <- data.frame(t(map))
map_taxon_prop <- map_taxon/map_taxon$total
map_taxon_prop <- map_taxon_prop[,-5]
map_taxon_prop$arrange <- map_taxon_prop$unmapped
map_taxon_prop$sample <- rownames(map_taxon_prop)
map_taxon2 <- reshape2::melt(map_taxon_prop,id.vars = c("sample","arrange"))
map_taxon2$condition = factor(map_taxon2$variable, levels = c("superkingdom","phylum","genus","species","unmapped"))
super.map.median <- round(median(map_taxon_prop$superkingdom),4)
phylum.map.median <- round(median(map_taxon_prop$phylum),4)
genus.map.median <- round(median(map_taxon_prop$genus),4)
species.map.median <- round(median(map_taxon_prop$species),4)
unmap.median <- round(median(map_taxon_prop$unmapped),4)



library(ggplot2)
library(ggthemes)
unique(map_taxon2$condition)
display.brewer.all()
mycolors<-c(rev(brewer.pal(11,"Spectral")[2:5]), brewer.pal(11,"Spectral")[10])

p1 <- map_taxon2 %>%
  mutate(sample = fct_reorder(sample, arrange, .desc = F)) %>%
  ggplot(aes(x=sample, y=value, fill=condition)) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() +xlab("Samples")+ylab("Proportion (%)")+scale_x_discrete(labels=NULL)+
  theme(axis.ticks = element_blank(), axis.title.x = element_text(size = 25,face = "bold"), axis.text.x  = element_text(size=40), axis.text.y  = element_text(size=40))+
  theme(legend.title = element_text(size=30, face = "bold"), legend.text = element_text(size=30, face = "bold"))+
  scale_fill_manual(values=mycolors)

p1
ggsave("map_taxon.png", p1, width = 20, height = 10, dpi=600)
ggsave("map_taxon.pdf", p1, width = 20, height = 10)

