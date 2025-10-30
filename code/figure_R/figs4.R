##figs4a####

library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
library(stringr)

gene <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_humanprotein_sample.tsv")
gene$sum <-  rowSums(gene[,2:ncol(gene)],na.rm=T)
gene$order <- order(gene$sum,decreasing = T)
gene <- gene[,c(1,ncol(gene)-1,ncol(gene))]
symbol <- gene$Protein.Group
symbol <- data.frame(str_split_fixed(symbol,"[|]|_HUMAN",10))$X2
eg = bitr(symbol, fromType="UNIPROT", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id <- as.character(eg[,2])
id_uniprot <- gene[,1]

ego.mf <- enrichGO(gene=id,OrgDb = org.Hs.eg.db,ont="MF",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
gomf.map <- summary(ego.mf)
dotplot(ego.mf)
write.csv(gomf.map,"GNHSF_humanprot_to_GOmf_20250406.csv")
pdf("GNHSF_human_ego_mf.pdf")
dotplot(ego.mf,showCategory=10)
dev.off()

##abundant human proteins####
human <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_humanprotein_sample.tsv",sep="\t")
sum <- data.frame(rowSums(human[,2:ncol(human)],na.rm = T))
rownames(sum) <- human$Protein.Group


####figs4b####
human <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_humanprotein_cogcat_sample.txt",sep="\t")
catsum <- rowSums(human[,2:ncol(human)],na.rm = T)
catpercent <- data.frame("Category"=human$COGcat,"Sum"=catsum)
library (dplyr)
catpercent <- catpercent %>%  mutate(percent = (Sum/sum(Sum))*100)
catpercent <- arrange(catpercent,desc(percent))
catpercent$Category[11:nrow(catpercent)] <- "Others"
catpercent <- aggregate(catpercent[2:3],by=list(catpercent$Category),sum)
catpercent$pert= paste0( round(catpercent$percent,1),"%")
library(ggpubr)
ggpie(catpercent,"percent",label="Group.1",lab.pos = "out",fill="Group.1")+scale_fill_brewer(palette ="Set3")+
geom_text(aes(x = 1.2, label = pert), 
          position = position_stack(vjust = 0.5),
          size = 3, color = "black", fontface = "bold")
ggsave("figs4b.pdf")

micro <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_microprotein_cogcat_sample_1385.txt",sep="\t")
catsum <- rowSums(micro[,2:ncol(micro)],na.rm = T)
catpercent <- data.frame("Category"=micro$COGcat,"Sum"=catsum)
library (dplyr)
catpercent <- catpercent %>%  mutate(percent = (Sum/sum(Sum))*100)
catpercent <- arrange(catpercent,desc(percent))
catpercent$Category[11:nrow(catpercent)] <- "Others"
catpercent <- aggregate(catpercent[2:3],by=list(catpercent$Category),sum)
catpercent$pert= paste0( round(catpercent$percent,1),"%")
library(ggpubr)
ggpie(catpercent,"percent",label="Group.1",lab.pos = "out",fill="Group.1")+scale_fill_brewer(palette ="Set3")+
  geom_text(aes(x = 1.2, label = pert), 
            position = position_stack(vjust = 0.5),
            size = 3, color = "black", fontface = "bold")
ggsave("figs3e.pdf")




####figs4c#####
rm(list=ls())
setwd("~/Documents/gnhsf/data202509/code/figs4")
pro.genus <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1039/GNHSF_diann_IGC_humanswiss_filter5_genus_sample_1039.txt",sep="\t")
gen.genus <- read.delim("~/Documents/gnhsf/data202509/metagenomics/metagenome_genus_sample_1039_read_count.tsv",sep="\t")

##change old gen name
name.change <- read.delim("~/Documents/gnhsf/data202509/metadata/ncbi_name_change2023.txt",sep="\t")
for(i in 1:nrow(name.change)){
    gen.genus$X <- gsub(name.change$old.name[i],name.change$new.name[i],gen.genus$X)
}


pro.genus[is.na(pro.genus)] <- 0
gen.genus[is.na(gen.genus)] <- 0


rownames(pro.genus) <- pro.genus[,1]
rownames(gen.genus) <- gen.genus[,1]
pro.genus <- pro.genus[-1]
gen.genus <- gen.genus[-1]

x <- pro.genus
y <- gen.genus

x.mean <- data.frame("fea"=rownames(x),"mean.pro"=rowSums(x)/1039)
y.mean <- data.frame("fea"=rownames(y),"mean.gen"=rowSums(y)/1039)


library(stringr)
x.mean$fea<- data.frame(str_split_fixed(x.mean$fea,"__|[|]",20))$X14
y.mean$fea<- data.frame(str_split_fixed(y.mean$fea,"__|[|]",20))$X12

library(dplyr)

x.mean <- arrange(x.mean,desc(mean.pro))
x.mean1 <- x.mean[1:20,]

y.mean <- arrange(y.mean,desc(mean.gen))

y.mean1 <- y.mean[which(y.mean$fea %in% x.mean1$fea),]


mean1 <- merge(x.mean,y.mean1,by="fea")

library(ggplot2)
library(ggpubr)
library(ggrepel)
ggplot(mean1,mapping = aes(x=mean.pro,y=mean.gen))+geom_point()+geom_smooth(method="lm")+scale_x_log10()+scale_y_log10()+
  theme_bw()+stat_cor(data=mean1,method="spearman")+
  ggtitle("shared genus comparation")
ggsave(paste0("gnhsf_shared_genus_correlation20250402.pdf"),width = 10,height = 10)


####figs4d#####
pro.species <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1039/GNHSF_diann_IGC_humanswiss_filter5_species_sample_1039.txt",sep="\t")
gen.species <- read.delim("~/Documents/gnhsf/data202509/metagenomics/metagenome_species_sample_1039_read_count.tsv",sep="\t")


##change old gen name
name.change <- read.delim("~/Documents/gnhsf/data202509/metadata/ncbi_name_change2023.txt",sep="\t")
for(i in 1:nrow(name.change)){
  gen.species$X <- gsub(name.change$old.name[i],name.change$new.name[i],gen.species$X)
}


pro.species[is.na(pro.species)] <- 0
gen.species[is.na(gen.species)] <- 0


rownames(pro.species) <- pro.species[,1]
rownames(gen.species) <- gen.species[,1]
pro.species <- pro.species[-1]
gen.species <- gen.species[-1]

x <- pro.species
y <- gen.species


x.mean <- data.frame("fea"=rownames(x),"mean.pro"=rowSums(x)/1039)
y.mean <- data.frame("fea"=rownames(y),"mean.gen"=rowSums(y)/1039)

library(stringr)
x.mean$fea<- data.frame(str_split_fixed(x.mean$fea,"__|[|]",20))$X16
y.mean$fea<- data.frame(str_split_fixed(y.mean$fea,"__|[|]",20))$X14

library(dplyr)

x.mean <- arrange(x.mean,desc(mean.pro))
x.mean1 <- x.mean[1:50,]

y.mean <- arrange(y.mean,desc(mean.gen))

y.mean1 <- y.mean[which(y.mean$fea %in% x.mean$fea),]

mean1 <- merge(x.mean,y.mean1,by="fea")

library(ggplot2)
library(ggpubr)
ggplot(mean1,mapping = aes(x=mean.pro,y=mean.gen))+geom_point()+scale_x_log10()+scale_y_log10()+
  geom_smooth(method="lm",se=T)+
  theme_bw()+stat_cor(data=mean1,method="spearman")+
  ggtitle("shared species comparation")
ggsave(paste0("gnhsf_all_shared_compare_species.pdf"),width = 10,height = 10)
