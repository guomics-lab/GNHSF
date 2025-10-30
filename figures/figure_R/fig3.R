
#####age.glm.glmm#######
age <- read.delim("age.glm.glmm_20250417.csv",sep=",")
age1 <- age[which(age$Pr...t...y<0.05),]
kegg.list2 <- age1$feature[age1$level=="kegg"]
cog.list2 <- age1$feature[age1$level=="cog"]
library(MicrobiomeProfiler)
#run_MicrobiomeProfiler()
library(openxlsx)
anno <- read.xlsx("kegg_ko_database20230327.xlsx")
ko.anno <- subset(anno,ko %in% age1$feature)
ko.anno <- merge(age1,ko.anno,by.x="feature",by.y="ko")
plot.ko <- unique(ko.anno[,c(21,5)])
library(dplyr)
plot.ko <- arrange(plot.ko,Estimate.x)
barplot(plot.ko$Estimate,names.arg = plot.ko$des,las=2,cex.names=0.2,horiz = TRUE)

core.prokegg.enrich <- enrichKO(kegg.list2,qvalueCutoff = 0.05,minGSSize = 10,maxGSSize = 500,pvalueCutoff = 0.1)
prokegg <- core.prokegg.enrich@result
dotplot(core.prokegg.enrich)
library(stringr)
prokegg$bg <- as.numeric(data.frame(str_split_fixed(prokegg$BgRatio,"/",2))$X1) / as.numeric(data.frame(str_split_fixed(prokegg$BgRatio,"/",2))$X2)
prokegg$gene <- as.numeric(data.frame(str_split_fixed(prokegg$GeneRatio,"/",2))$X1) / as.numeric(data.frame(str_split_fixed(prokegg$GeneRatio,"/",2))$X2)
prokegg$foldenrich <- prokegg$gene/prokegg$bg
prokegg1 <- prokegg[which(prokegg$qvalue<0.05),]
library(tidyr)
prokegg2 <- prokegg1 %>%as_tibble() %>%    separate_rows(geneID, sep = "/")

prokegg2 <- merge(prokegg2,age1,by.x="geneID",by.y = "feature")
prokegg2$Description <- factor(prokegg2$Description,levels=rev(prokegg1$Description))
library(ggplot2)
ggplot(prokegg2,aes(x=Description,y=Estimate.x))+geom_point()+coord_flip()


####age micro####
annotate1 <- read.xlsx("microprot_annotation.xlsx")
age.micro <- subset(age1,level=="micro")
age.micro1 <- merge(age.micro,annotate1,by.x="feature",by.y="glm.match")
#write.csv(age.micro1,"age_glm_glmm_micro.csv")
age.micro1 <- read.csv("age_glm_glmm_micro.csv")
library(tidyr)
age.micro2 <-age.micro1 %>% as_tibble() %>%    separate_rows(Taxon_rank_unipept,Taxon_name_unipept, sep = ";")
age.micro2 <- age.micro2[which(age.micro2$Taxon_rank_unipept %in% c("phylum","genus","species")),]
library(reshape2)
age.micro3 <- dcast(age.micro2,feature+Protein+Eggpro_des+Eggpro_prefer_name+COG_NCBIcat+COG_NCBIcat_des+Estimate.x+Estimate.y~Taxon_rank_unipept)
age.micro3$name <- age.micro3$Eggpro_prefer_name
age.micro3$value <- abs(age.micro3$Estimate.x)
library(ggraph)
library(tidygraph)
library(ccgraph)
nodes <- gather_graph_node(age.micro3,index=colnames(age.micro3)[c(10,9,11,12)],value="value")
edge <- gather_graph_edge(age.micro3,index=colnames(age.micro3)[c(10,9,11,12)])
graph<- tbl_graph(nodes,edge)

get.inform <- create_layout(graph , 'circlepack')

ggraph(graph, 'circlepack',weight=node.size) + 
  geom_node_circle(aes(fill = factor(node.level)),alpha=0.5) + 
  coord_fixed() +
  theme_graph()+ 
  geom_node_text(
    aes(label=node.short_name,color=leaf),repel=T,
    fontface="bold",
    size=3
  )+ scale_fill_manual(values=c("#66bc9870","#aad09d70","#e3ea96","#fcdc89"))+scale_color_manual(values = c("black","darkblue")) +theme_map()
ggsave("identification_category_count.pdf",width=10,height=10)

library(sankeyD3plus)
sangi <- age.micro3[,c(6,12,11,9,10,13)]
sangi$color <- as.factor(sangi$genus)
sangi$color <- as.numeric(sangi$color)
#sangi$color[is.na(sangi$color)] <- "17"
sangi <- sangi[-which(is.na(sangi$genus)),]
c1.1 <- na.omit(sangi[,c(1,2,7)])
colnames(c1.1) <- c("source","target","color")
c1.1.2 <- na.omit(sangi[,c(2,3,7)])
colnames(c1.1.2) <- c("source","target","color")
c1.2 <- na.omit(sangi[,c(3,4,7)])
colnames(c1.2) <- c("source","target","color")
c1.2.2 <- na.omit(sangi[which(is.na(sangi$species)),c(2,4,7)])
colnames(c1.2.2) <- c("source","target","color")
c1.3 <- na.omit(sangi[,c(4,5,7)])
colnames(c1.3) <- c("source","target","color")
c1.3.2 <- na.omit(sangi[which(is.na(sangi$genus)),c(2,5,7)])
colnames(c1.3.2) <- c("source","target","color")

c1 <- rbind(c1.1,c1.1.2,c1.2,c1.3,c1.2.2,c1.3.2)
c1.5 <- data.frame(table(c1))
c1.5 <- subset(c1.5,Freq !=0)

id1.1 <- unique(data.frame("cat"="Category","value"=sangi$COG_NCBIcat_des))
id1 <- unique(data.frame("cat"="Protein","value"=sangi$name))
id2 <- na.omit(unique(data.frame("cat"="Genus","value"=sangi$genus)))
id3 <- unique(data.frame("cat"="Phylum","value"=sangi$phylum))
id4 <-  na.omit(unique(data.frame("cat"="Species","value"=sangi$species)))


color_map <- c(
  "0" = "#E6F598",  "1" = "#F39A9B", "2" = "#a7c6e0", "3" = "#b1dcb0",
  "4" = "#FFc78f", "5" = "#FFFFa5", "6" = "#d8b5a1", "7" = "#Fbc8e3",
  "8" = "#cc8fcc", "9" = "#b8dfd3", "10" = "#FC8D62", "11" = "#8DA0CB",
  "12" = "#E78AC3", "13" = "#A6D854", "14" = "#FFB6B0", "15" = "#8C510A",
  "16" = "#D8B365", "17" = "#5AB4AC", "18" = "#D9D900", "19" = "#FF1493",
  "20" = "#E6F598"
)



id <- rbind(id1.1,id1,id2,id3,id4)
id$id <- rownames(id)
id$xpos <- factor(id$cat,levels=c("Category","Protein","Species","Genus","Phylum"))
id$xpos <- as.numeric(id$xpos)-1
id$color <- id$xpos+5
id$color <- color_map[id$color]

c1.5$IDsource <- match(c1.5$source, id$value) - 1 
c1.5$IDtarget <- match(c1.5$target, id$value) - 1
c1.5$colour <- as.character(c1.5$color)
c1.5$colour <- color_map[c1.5$colour]

library(sankeyD3plus)
p <- sankeyNetwork(Links = c1.5, Nodes = id,units = "normalized",nodePadding = 2,nodeWidth = 10,width=1200,height=500, LinkGroup = "colour",linkColor = "colour",
                   Source = 'IDsource', Target = 'IDtarget', Value = 'Freq',dragY = T, dragX = F,showNodeValues = F,scaleNodeBreadthsByString=T,
                   NodeColor = "color",NodeID = 'value',NodePosX = "xpos", #xAxisDomain = c("Category","Protein","Species","Genus","Phylum"),
                   fontSize = 12)

p


saveNetwork(p,file_name =  "age_micro", selfcontained = TRUE,html = "create")

library(webshot2)
webshot2::webshot("age_micro.html", "age_micro.pdf")



####venn#####
age1 <- age[which(age$Pr...t...y<0.05),]
age.venn <- age[,c(2,13,10,15,16)]
age.venn1 <- age.venn

aaa <- split(age.venn1,age.venn1$level)
library(Vennerable)  
lapply(aaa,function(x){
  between.subject <- x$feature[which(x$q.x<0.1)]
  within.subject <- x$feature[which(x$Pr...t...y<0.05)]
  if(length(within.subject)==0){
    within.subject=NA
  }
  data<-compute.Venn(Venn(list("Between-subject correlation"=between.subject,"Within-subject correlation"=within.subject))) 
  gp <- VennThemes(data)
  gp[["Face"]][["11"]]$fill <-"#9392be"
  gp[["Face"]][["01"]]$fill <-"#d0e7ed"
  gp[["Face"]][["10"]]$fill <-"#d5e4a8"
  gp[["Set"]][["Set2"]]$col <- NA
  gp[["Set"]][["Set1"]]$col <- NA
  gp[["SetText"]][["Set1"]]$col <-"black"
  gp[["SetText"]][["Set2"]]$col <-"black"
  pdf(paste0(unique(x$level),".venn.pdf"))
  plot(data,gp=gp)  
  dev.off()  
}  )

####age species 1385###

genus <- read.delim("GNHSF_diann1.8_IGC_humanswiss_genus_sample_1385_2_rmallNA.tsv",sep="\t")
species <- read.delim("GNHSF_diann1.8_IGC_humanswiss_species_sample_1385_2_rmallNA.tsv",sep="\t")
meta <- read.delim("GNHSF_metadata_1385_20220515.csv",sep=",")

genus$Taxon <- gsub("[|]",".",genus$Taxon)
genus1 <- genus[which(genus$Taxon %in% age1$feature),]
colnames(genus1)[1] <- "feature"
species$Taxon <- gsub("[|]",".",species$Taxon)
species1 <- species[which(species$Taxon %in% age1$feature),]
colnames(species1)[1] <- "feature"

feature <- rbind(genus1,species1)
feature <- data.frame(t(feature))
colnames(feature) <- feature[1,]
feature <- feature[-1,]
feature$sample <- rownames(feature)
meta1 <- meta[,1:3]

feature1 <- merge(meta1,feature,by="sample")

feature1[,2:ncol(feature1)] <- apply(feature1[,2:ncol(feature1)],2,function(x) as.numeric(x))
feature1[is.na(feature1)] <- 0
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(ggpointdensity)

ggplot(feature1,aes(x=Age,y=k__Bacteria.p__Firmicutes.c__Negativicutes.o__Acidaminococcales.f__Acidaminococcaceae.g__Phascolarctobacterium))+geom_pointdensity()+
  scale_color_continuous(low = "#d9bdd8", high = "#83247f")+
  geom_smooth(method="lm",color="#3c2260",fill="#d9bdd8",size=1.2)+stat_cor()+theme_classic()
#stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=5)
ggsave("Phascolarctobacterium_age.pdf")


ggplot(feature1,aes(x=Age,y=k__Bacteria.p__Firmicutes.c__Negativicutes.o__Acidaminococcales.f__Acidaminococcaceae.g__Phascolarctobacterium.s__Phascolarctobacterium_faecium))+geom_pointdensity()+
  scale_color_continuous(low = "#d9bdd8", high = "#83247f")+
  geom_smooth(method="lm",color="#3c2260",fill="#d9bdd8",size=1.2)+stat_cor()+theme_classic()
#stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=5)
ggsave("Phascolarctobacterium_faecium_age.pdf")


ggplot(feature1,aes(x=Age,y=k__Bacteria.p__Firmicutes.c__Negativicutes.o__Veillonellales.f__Veillonellaceae.g__Dialister))+geom_pointdensity()+
  scale_color_continuous(low = "#f6b6af", high = "#ee3432")+
  geom_smooth(method="lm",color="#3c2260",fill="#f6b6af",size=1.2)+stat_cor()+theme_classic()
#stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=5)
ggsave("Dialister_age.pdf")

