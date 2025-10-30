rm(list=ls())

##glm related T2D species#######
all <- data.frame(data.table::fread("~/Documents/gnhsf/data202509/code/fig2/all_GLM.csv",sep=","))
all.t2d <- subset(all,metadata=="dm_cl")
all.t2d.sig <- subset(all.t2d,q.sig !="nsig")
table(all.t2d.sig$level)
#cog        genus humanprotein         kegg microprotein      species 
#269           23           20          297         1586           24 
write.xlsx(all.t2d.sig,"all_t2d_sig.xlsx")

t2d.species <- subset(all.t2d.sig,level=="species")

t2d.micro <- subset(all.t2d.sig,level=="microprotein")

micro.anno <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/08.link/microbial_prot_to_unique_species_link.txt",sep="\t")
micro.anno$glm.match <- gsub("-",".",micro.anno$Protein)
micro.anno$glm.match <- paste0("X",micro.anno$glm.match)
micro.anno$Taxon_name_unipept <- gsub("_"," ",micro.anno$Taxon_name_unipept)

t2d.micro.species <- merge(t2d.micro,micro.anno,by.x = "feature",by.y="glm.match")
micro.species <- data.frame(table(t2d.micro.species$Taxon_name_unipept))
#total 87
micro.species <- micro.species[-which(micro.species$Freq %in% c(1,2,3,4)),]
#leave at least 2 proteins, 35
#leave at least 3 proteins, 27
#leave at least 5 proteins, 16

t2d.species2 <- data.frame("Var1"=unlist(lapply(strsplit(t2d.species$feature,"[.]s_"), '[',2)),"Freq"=t2d.species$Estimate)
library(openxlsx)
t2d.micro.species1 <- subset(t2d.micro.species,Taxon_name_unipept %in% micro.species$Var1)
write.xlsx(t2d.micro.species1,"t2d_micro_species_anno202509.xlsx")

species <- rbind(t2d.species2,micro.species)
species$Var1 <- gsub("_"," ",species$Var1)

length(unique(species$Var1))       
##total 33 species,pep5####       


#gengerate table s9#
all.t2d.taxon <- subset(all.t2d,level %in% c("genus","species"))
feature.name<- strsplit(all.t2d.taxon$feature,"[.]|__")
feature.name1 <- do.call(c,lapply(feature.name, function(x) {
  aa <- x[length(x)]
  return(aa)
}))
feature.name1 <- gsub("^s_","",feature.name1)
feature.name1 <- gsub("_"," ",feature.name1)
all.t2d.taxon$feature2 <- feature.name1

glm.species <- all.t2d.taxon[which(all.t2d.taxon$feature2 %in% species$Var1),]
species.pmid <- read.xlsx("33species.xlsx")
micro <- read.xlsx("t2d_micro_species_anno_1_20250923.xlsx")
micro <- subset(micro,!is.na(Cluster))

glm.species <- merge(glm.species,species.pmid,by.x="feature2",by.y="Var1")
glm.species.micro <- merge(glm.species,micro,by.x="feature2",by.y="Taxon_name_unipept",all = T)
write.xlsx(glm.species.micro,"tables9.xlsx")

##figs8b####
genus.list <- all.t2d.sig$feature[which(all.t2d.sig$level %in% c("genus"))]
species.list <- unique(species$Var1)

all.t2d.taxon <- subset(all.t2d,level %in% c("genus","species"))
feature.name<- strsplit(all.t2d.taxon$feature,"[.]|__")
feature.name1 <- do.call(c,lapply(feature.name, function(x) {
  aa <- x[length(x)]
  return(aa)
}))
feature.name1 <- gsub("^s_","",feature.name1)
feature.name1 <- gsub("_"," ",feature.name1)
all.t2d.taxon$feature2 <- feature.name1

taxon.pro1 <- all.t2d.taxon[which(all.t2d.taxon$feature %in% c(genus.list,species.list)),]
taxon.pro2 <- all.t2d.taxon[which(all.t2d.taxon$feature2 %in% c(genus.list,species.list)),]
taxon.pro <- rbind(taxon.pro1,taxon.pro2)

genus.gen <- read.delim("~/Documents/gnhsf/data202509/code/glm/glm_1039/genus_NA90_all_meta_glm_res_adj_ex4.csv",sep=",")
species.gen <- read.delim("~/Documents/gnhsf/data202509/code/glm/glm_1039/species_NA90_all_meta_glm_res_adj_ex4.csv",sep=",")
genus.gen$level <- "genus"
species.gen$level <- "species"
taxa.gen <- rbind(genus.gen,species.gen)
feature.name<- strsplit(taxa.gen$feature,"[.]|__")
feature.name1 <- do.call(c,lapply(feature.name, function(x) {
  aa <- x[length(x)]
  return(aa)
}))
feature.name1 <- gsub("^s_","",feature.name1)
feature.name1 <- gsub("_"," ",feature.name1)
taxa.gen$feature2 <- feature.name1
taxa.gen.t2d <- subset(taxa.gen,metadata=="dm_cl")


plot.list <- c(taxon.pro$feature2,taxa.gen.t2d$feature2[which(taxa.gen.t2d$q.sig!="nsig")])
taxon.pro1 <- all.t2d.taxon[which(all.t2d.taxon$feature2 %in% plot.list),]
taxon.pro1$cat <- "metaproteomics"
taxon.pro1 <- taxon.pro1[,c(2,13,15,16,17)]
taxon.gene1 <- taxa.gen.t2d[which(taxa.gen.t2d$feature2 %in% plot.list),]
taxon.gene1$cat <- "metagenomics"
taxon.gene1 <- taxon.gene1[,c(1,9,10,11,12)]
taxon.plot <- rbind(taxon.gene1,taxon.pro1)

taxon.plot2 <- reshape2::dcast(taxon.plot,feature2+level~cat,value.var = c("Estimate"))
taxon.plot2.sig <- reshape2::dcast(taxon.plot,feature2+level~cat,value.var = c("q.sig"))
taxon.plot2.sig[taxon.plot2.sig=="nsig"] <- ""
taxon.plot2.sig[is.na(taxon.plot2.sig)] <- ""
rownames(taxon.plot2) <- taxon.plot2$feature2


library(ComplexHeatmap)
library(dplyr)
breaks <- c(seq(-1, 0, length.out = 8), seq(0, 1.5, length.out = 8)) %>% unique()
colors <- c(
  colorRampPalette(c("cyan3", "lightblue"))(7), 
  colorRampPalette(c("thistle2","lightcoral"))(7)  
)

pdf("t2d_taxa_metapro_metagen_match.pdf")
pheatmap(taxon.plot2[,3:4],display_numbers = as.matrix(taxon.plot2.sig[,3:4]),cluster_rows = F,cluster_cols = F,na_col = "white",breaks = breaks,color=colors,row_split=taxon.plot2$level,annotation_row = taxon.plot2[,2,drop=F],show_rownames = TRUE)
dev.off()


aa <- pheatmap(taxon.plot2[,3:4],display_numbers = as.matrix(taxon.plot2.sig[,3:4]),cluster_rows = F,cluster_cols = F,na_col = "white",breaks = breaks,color=colors,row_split=taxon.plot2$level,annotation_row = taxon.plot2[,2,drop=F],show_rownames = TRUE)
row_order <- row_order(aa)
taxon.ep <- taxon.plot2[unlist(row_order),]
write.csv(taxon.ep,"figs8b_rebuttal_validata.csv")



library(stringr)
#read species 1385, for cluster network####
species1385 <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_filter5_species_sample_1385.txt",sep="\t")
name <- data.frame(str_split_fixed(species1385$Taxon,"[|]|__",20))$X16
name <- gsub("_"," ",name)
species1385$Taxon <- name
species1385 <- species1385[which(species1385$Taxon %in% species$Var1),]
species1385[is.na(species1385)] <- 0

library(Hmisc)
sp.cor <- rcorr(t(species1385[,-1]),type="spearman")
r.cor <- sp.cor$r
p.cor <- sp.cor$P
colnames(r.cor) <- species1385$Taxon
rownames(r.cor) <- species1385$Taxon

r.cor2 <- r.cor
r.cor2[which(abs(r.cor2)<0.3)] <- 0
r.cor2[which(r.cor2==1)]=0
p.cor <- p.adjust(p.cor,method="BH")
r.cor3 <- r.cor2
r.cor3[which(p.cor>0.05)] <- 0
write.csv(r.cor,"t2d_species_cor_origin1.csv")
##fig6a to gephi for network####
write.csv(r.cor3,"t2d_species_cor_0.3_5pep.csv")




##fig6a generate pie### 
t2d.micro.species <- read.xlsx("t2d_micro_species_anno_1.xlsx")
t2d.micro.species[t2d.micro.species==""] <- NA
t2d.micro.species$COG_NCBIcat <- ifelse(t2d.micro.species$COG_NCBIcat %in% c("C","G","E","F","I",'M',"H","O","J"),t2d.micro.species$COG_NCBIcat,"Others")

micro.anno.pie <- t2d.micro.species[which(t2d.micro.species$Taxon_name_unipept %in% species$Var1),]

micro.pie <- split(micro.anno.pie,micro.anno.pie$Taxon_name_unipept)

library(ggpubr)
library (dplyr)

for (i in 1:29) {
  x=micro.pie[[i]]
  name =names(micro.pie)[i]
  x.up <- subset(x,Estimate>0)
  if(nrow(x.up)>0){
    cogcat <- data.frame(table(x.up$COG_NCBIcat))
    catpercent <- cogcat %>%  mutate(percent = (Freq/sum(Freq))*100)
    catpercent <- arrange(catpercent,desc(percent))
    catpercent$pert= paste0( round(catpercent$percent,1),"%")
    ggpie(catpercent,"Freq",label = catpercent$Var1,fill="Var1",lab.pos = "in")+scale_fill_manual(values = c("C"="#9ac2e8","G"="#da7a83","E"="#d2df4e","F"="#c8a063","I"="#b4a6cf"
                                                                                                             ,"J"="#90ddcc","M"="#e75541","O"="#fee030","H"="#119789","Others"="#f0e4dd"))
    ggsave(paste0(name,"_pie_up.pdf"))                                                                                                         
  }
  x.down <- subset(x,Estimate<0)
  if(nrow(x.down)>0){
    cogcat <- data.frame(table(x.down$COG_NCBIcat))
    catpercent <- cogcat %>%  mutate(percent = (Freq/sum(Freq))*100)
    catpercent <- arrange(catpercent,desc(percent))
    catpercent$pert= paste0( round(catpercent$percent,1),"%")
    ggpie(catpercent,"Freq",label = catpercent$Var1,fill="Var1",lab.pos = "in")+scale_fill_manual(values = c("C"="#9ac2e8","G"="#da7a83","E"="#d2df4e","F"="#c8a063","I"="#b4a6cf"   ,"J"="#90ddcc","M"="#e75541","O"="#fee030","H"="#119789","Others"="#f0e4dd"))
      ggsave(paste0(name,"_pie_down.pdf"))                                                                                                         
  }
}

## all signif species and their signif proteins, and unsignif species with at least 5 signif proteins


for (i in 1:28) {
  x=micro.pie[[i]]
  name =names(micro.pie)[i]
    cogcat <- data.frame(table(x$COG_NCBIcat))
    catpercent <- cogcat %>%  mutate(percent = (Freq/sum(Freq))*100)
    catpercent <- arrange(catpercent,desc(percent))
    catpercent$pert= paste0( round(catpercent$percent,1),"%")
    ggpie(catpercent,"Freq",label = catpercent$Var1,fill="Var1",lab.pos = "in")+scale_fill_manual(values = c("C"="#9ac2e8","G"="#da7a83","E"="#d2df4e","F"="#c8a063","I"="#b4a6cf"
                                                                                                             ,"J"="#90ddcc","M"="#e75541","O"="#fee030","H"="#119789","Others"="#f0e4dd"))
    ggsave(paste0(name,"_pie.pdf"))                                                                                                         
}


##figs9a#####
t2d.micro.species <- read.xlsx("t2d_micro_species_anno_1.xlsx")
t2d.micro.species[t2d.micro.species==""] <- NA

data <- t2d.micro.species
cluster <-  data[which(!is.na(data$Cluster)),]
cluster.1 <- na.omit(cluster[,c(34,31,22)])
colnames(cluster.1) <-c("Cluster","Species","Protein category")

c.1 <- cluster.1[,c(1,2)]
colnames(c.1) <- c("source","target")
c.2 <- cluster.1[,c(2:3)]
colnames(c.2) <- c("source","target")
c <- rbind(c.1,c.2)
c.3 <- data.frame(table(c))
c.3 <- subset(c.3,Freq !=0)

id1 <- unique(data.frame("cat"="Cluster","value"=cluster.1$Cluster))
library(dplyr)
id1 <- arrange(id1,value)
id2 <- unique(data.frame("cat"="Species","value"=cluster.1$Species))
id3 <- unique(data.frame("cat"="Protein","value"=cluster.1$'Protein category'))
id <- rbind(id1,id2,id3)

c.3$IDsource <- match(c.3$source, id$value) - 1 
c.3$IDtarget <- match(c.3$target, id$value) - 1

c1.id <- c(unique(cluster$Taxon_name_unipept[ cluster$Cluster=="Cluster 1"]),"Cluster 1")
c2.id <- c(unique(cluster$Taxon_name_unipept[ cluster$Cluster=="Cluster 2"]),"Cluster 2")
c3.id <- c(unique(cluster$Taxon_name_unipept[ cluster$Cluster=="Cluster 3"]),"Cluster 3")

c.3$color <- ifelse(c.3$source %in% c1.id,"Cluster1",ifelse(c.3$source %in% c2.id,"Cluster2","Cluster3"))
id$color <- ifelse(id$value %in% c1.id,"Cluster1",ifelse(id$value %in% c2.id,"Cluster2","Cluster3"))


color <- 'd3.scaleOrdinal() .domain([ "Cluster1", "Cluster2", "Cluster3"]) 
.range(["#ffb4ac","#c7dff3", "#cbe8a3"])'


library(networkD3)
p <- sankeyNetwork(Links = c.3, Nodes = id,colourScale = color,
                   Source = 'IDsource', Target = 'IDtarget', Value = 'Freq', 
                   NodeID = 'value', NodeGroup = 'color', LinkGroup = "color",
                   fontSize = 12, sinksRight = FALSE,iterations = 4000)
p

saveNetwork(p,"figs9a_sangi.html")


####figs9b#### both gnhs and fh cohort, t2d####
glm1385 <- data.frame(data.table::fread("~/Documents/gnhsf/data202509/code/fig2/all_GLM.csv",sep=","))
glm.t2d1385 <- subset(glm1385,q.sig !="nsig"&metadata=="dm_cl")
glm.t2d1385$feature <- gsub("__","_",glm.t2d1385$feature )

table(glm.t2d1385$q.sig,glm.t2d1385$level)
#     cog genus humanprotein kegg microprotein species
#*   139    19           16  149          784      18
#**   88     3            4   90          456       4
#***  42     1            0   58          346       2

fh <- read.csv("FH_glm_all.csv")
glm.t2d104 <- subset(fh,q.sig !="nsig"&metadata=="dm_clyes")
length(which(glm.t2d104$feature %in% glm.t2d1385$feature))
##159 overlap###

overlap <- glm.t2d104[ which(glm.t2d104$feature %in% glm.t2d1385$feature),]
write.xlsx(overlap,"tableS10.xlsx")



table(overlap$level)

#cog        genus humanprotein         kegg microprotein      species 
#10            2            1           16          129            1 
library(openxlsx)
t2d.micro.species <- read.xlsx("t2d_micro_species_anno_1_20250923.xlsx")
overlap.anno <- t2d.micro.species[which(t2d.micro.species$feature %in% overlap$feature),]
#overlap 82 microprotein anno into species###
table(overlap.anno$Taxon_name_unipept)
write.xlsx(overlap.anno,"fh_validate_prot_to_species.xlsx")

tables10 <- read.xlsx("tableS10.xlsx")
tables10 <- merge(tables10,overlap.anno[,c(16,31)],by.y="Protein",by.x="Feature",all.x=T)
write.xlsx(tables10,"tableS10.xlsx")

#####t2d fh 82 prot####
fh <- overlap.anno
fh[is.na(fh)] <- "unknown"
fh$group <- ifelse(fh$Taxon_name_unipept=="Megasphaera elsdenii","G1",ifelse(fh$Taxon_name_unipept=="Megasphaera","G2","G3"))
fh$Protein <- paste0(fh$feature,"(",fh$Eggpro_prefer_name,")")
fh.1 <- fh[,c(16,22,31,34)]
fh.2 <- fh.1[,c(1,2,4)]
colnames(fh.2)[1:2] <- c("source","target")
fh.3 <- fh.1[,c(2,3,4)]
colnames(fh.3)[1:2] <- c("source","target")
fh.4 <- rbind(fh.2,fh.3)

fh.5 <- data.frame(table(fh.4))
fh.5 <- subset(fh.5,Freq !=0)
fh.5$Cluster <- gsub(" ","",fh.5$Cluster)

id1 <- unique(data.frame("cat"="Metaprotein","value"=fh.1$Protein))
id2 <- unique(data.frame("cat"="COGcat","value"=fh.1$COG_NCBIcat_des))
id3 <- unique(data.frame("cat"="Species","value"=fh.1$Taxon_name_unipept))
id <- rbind(id1,id2,id3)

fh.5$IDsource <- match(fh.5$source, id$value) - 1 
fh.5$IDtarget <- match(fh.5$target, id$value) - 1
fh.5$color=fh.5$Cluster
color <- 'd3.scaleOrdinal() .domain([ "Cluster1", "Cluster2", "Cluster3","unknown"]) 
.range(["#ffb4ac","#c7dff3", "#cbe8a3","#d6b9d0"])'


library(networkD3)
p <- sankeyNetwork(Links = fh.5, Nodes = id,nodePadding = 5,colourScale = color, 
                   Source = 'IDsource', Target = 'IDtarget', Value = 'Freq', 
                   NodeID = 'value', NodeGroup = 'cat', LinkGroup = "color",
                   fontSize = 12, sinksRight = FALSE,iterations = 5000)
p
saveNetwork(p,"sangi_fh_validate_micro.html")

#fig s11 a####
library(openxlsx)
micro.anno <- read.xlsx("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microprot_2phylum2genus2species_one20250916.xlsx")
glm <- data.frame(data.table::fread("~/Documents/gnhsf/data202509/code/fig2/all_GLM.csv"))
micro.glm <- subset(glm,level=="microprotein")

micro.anno$glm.match <- gsub("-",".",micro.anno$Protein)
micro.anno$glm.match <- paste0("X",micro.anno$glm.match)
micro.anno$glm.match <- gsub("[.]$","",micro.anno$glm.match)
mega.anno <- micro.anno[grep("Megasphaera_elsdenii",micro.anno$Taxon_all_unipept),]

micro.glm <- merge(micro.glm,mega.anno,by.x="feature",by.y="glm.match")
micro.glm.dm <- subset(micro.glm,metadata %in% c("dm_cl","dm_med"))
write.xlsx(micro.glm.dm,"mega_microprot_glm_dm.xlsx")
#manual impute gene name

micro.glm.dm <- read.xlsx("mega_microprot_glm_dm_1.xlsx")
library(reshape2)
mega.esti <- dcast(micro.glm.dm,Protein+Eggpro_des+Eggpro_prefer_name+COG+COG_NCBIdes~metadata,value.var=c("Estimate"))
mega.esti.sig <- dcast(micro.glm.dm,Protein+Eggpro_des+Eggpro_prefer_name+COG+COG_NCBIdes~metadata,value.var=c("q.sig"))
del=which(mega.esti.sig$dm_med=="nsig")
mega.esti1 <- mega.esti[-del,]
rownames(mega.esti1) <- paste0(mega.esti1$Protein,"(",mega.esti1$Eggpro_prefer_name,")")
mega.esti.sig1 <- mega.esti.sig[-del,]
mega.esti.sig1[mega.esti.sig1=="nsig"] <- ""
rownames(mega.esti.sig1) <- paste0(mega.esti1$Protein,"(",mega.esti1$Eggpro_prefer_name,")")
library(ComplexHeatmap)
library(dplyr)
breaks <- c(seq(-0.4, 0, length.out = 6), seq(0, 0.6, length.out = 6)) %>% unique()
colors <- c(
  colorRampPalette(c("cyan3", "lightblue"))(5),  # 负数区间6个
  colorRampPalette(c("thistle2","lightcoral"))(5)   # 正数区间6个
)
pdf("mega_prot_dm_glm_est.pdf")
pheatmap(mega.esti1[,6:7],display_numbers = as.matrix(mega.esti.sig1[,6:7]),cluster_cols = F,breaks = breaks,color=colors)
dev.off()

aa <- pheatmap(mega.esti1[,6:7],display_numbers = as.matrix(mega.esti.sig1[,6:7]),cluster_cols = F,breaks = breaks,color=colors)
row_order <- row_order(aa)
megalist.ep <- mega.esti1[unlist(row_order),]
write.csv(megalist.ep,"figs11a_rebuttal_validata.csv")



####fig6B#### glm and ml features, network
rm(list=ls())
library(openxlsx)



glm1385 <- data.frame(data.table::fread("~/Documents/gnhsf/data202509/code/fig2/all_GLM.csv",sep=","))
glm.t2d1385 <- subset(glm1385,q.sig !="nsig"&metadata=="dm_cl")
glm.t2d1385$feature <- gsub("__","_",glm.t2d1385$feature )

table(glm.t2d1385$q.sig,glm.t2d1385$level)
#     cog genus humanprotein kegg microprotein species
#*   139    19           16  149          784      18
#**   88     3            4   90          456       4
#***  42     1            0   58          346       2

fh <- read.csv("FH_glm_all.csv")
glm.t2d104 <- subset(fh,q.sig !="nsig"&metadata=="dm_clyes")
length(which(glm.t2d104$feature %in% glm.t2d1385$feature))
##159 overlap###

overlap <- glm.t2d104[ which(glm.t2d104$feature %in% glm.t2d1385$feature),]
table(overlap$level)

#cog        genus humanprotein         kegg microprotein      species 
#10            2            1           16          129            1 



level <- c("species","kegg","humanprotein","microprotein")
#level <- c("species","kegg","humanprotein","microprotein","genus","cog")
fold <- c("fold_0","fold_1","fold_2","fold_3","fold_4")
list <- c()
for (x in level) {
  for (y in fold) {
    list1 <- read.delim(paste0("MLoutput/",x,"/selected_features/",y,"/xgb_high_freq_features.csv"),sep=',')
    list1$level <- x
    list1$fold <- y
    list <- rbind(list,list1)
  }
}
write.xlsx(list,"Tables11_ml_t2d_related.xlsx")
list.ml <- unique(list[,c(1,4)])
table(list.ml$level)

overlap$feature <- gsub("^X","",overlap$feature)
glm.list <- subset(overlap,level %in% c("kegg","humanprotein","microprotein"))
glm.list <- glm.list[,c(6,14)]
#del glm.list- m.elsdenii , cause list.ml already conclude it

total.list <- unique(rbind(list.ml,glm.list))



#micro-prot link
#prot-ko link
#micro-ko link
#humanp-ko, corr
#humanp-species, corr
#between each level, corr

total.df <- data.frame(matrix(nrow=length(total.list$feature),ncol=length(total.list$feature)))
#humanprotein    kegg microprotein      species 
#31           96          163           55 

rownames(total.df) <- total.list$feature
colnames(total.df) <- total.list$feature

#link1. micro-prot.link
prot.link <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/08.link/microbial_prot_to_unique_species_link.txt",sep="\t")
prot.link$Protein <- gsub("-",'.',prot.link$Protein)
prot.link <- prot.link[,c(1,15,16)]
library(stringr)
prot.link2 <- data.frame(str_split_fixed(prot.link$Taxon_all_unipept,"[|]",10))
prot.link2$species <- prot.link$Taxon_name_unipept
prot.link2[prot.link2==""] <- NA
prot.link2$X8 <- ifelse(is.na(prot.link2$X8),paste0("s__",prot.link2$species),prot.link2$X8)
prot.link$Taxon_all_unipept <- apply(prot.link2[,1:8],1,function(x)  paste(x,collapse  = "."))

prot.link.m <- prot.link[which(prot.link$Protein %in% total.list$feature), ]
length(unique(prot.link.m$Taxon_all_unipept))
#95 prot, to 23 species
length(which(unique(prot.link.m$Taxon_all_unipept) %in% total.list$feature))
#10 species match

for(x in 1:nrow(prot.link.m)){
  aa <- which(colnames(total.df)==prot.link.m$Protein[x])
  bb <- which(rownames(total.df)==prot.link.m$Taxon_all_unipept[x])
  if(length(aa)>0& length(bb)>0){
    total.df[aa,bb] <- 1
    total.df[bb,aa] <- 1
    print(x)
  }
}

#link2. prot-ko link
prot.ko <- data.frame(data.table::fread("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microbiome_all.txt",sep="\t"))
prot.ko <- unique(prot.ko[,c(2,11)])
prot.ko$Protein <- gsub("-",'.',prot.ko$Protein)
prot.ko1 <- prot.ko[which(prot.ko$Protein %in% total.list$feature),]
#total 163, match
library(tidyr)
prot.ko1 <-prot.ko1 %>% as_tibble() %>%    separate_rows(KEGG, sep = "/")
prot.ko1[prot.ko1==""] <- NA
prot.ko1 <- na.omit(prot.ko1)
length(which(unique(prot.ko1$KEGG) %in% total.list$feature))
#17 kegg match
for(x in 1:nrow(prot.ko1)){
  aa <- which(colnames(total.df)==prot.ko1$Protein[x])
  bb <- which(rownames(total.df)==prot.ko1$KEGG[x])
  if(length(aa)>0& length(bb)>0){
    total.df[aa,bb] <- 1
    total.df[bb,aa] <- 1
    print(x)
  }
}

#link3, micro-ko link
micro.ko <- data.frame(data.table::fread("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microbiome_all.txt",sep="\t"))
micro.ko <- unique(micro.ko[,c(11,16)])
micro.ko[micro.ko==""] <- NA
micro.ko <- na.omit(micro.ko)
micro.ko <- micro.ko[grep("s__",micro.ko$Taxon_all_unipept),]
library(stringr)
micro.ko2 <- data.frame(str_split_fixed(micro.ko$Taxon_all_unipept,"[|]",10))
micro.ko$Taxon_all_unipept <- apply(micro.ko2[,1:8],1,function(x)  paste(x,collapse  = "."))
library(tidyr)
micro.ko1 <-micro.ko %>% as_tibble() %>%    separate_rows(KEGG, sep = "/")
micro.ko2 <- micro.ko1[which(micro.ko1$Taxon_all_unipept %in% total.list$feature),]
length(unique(micro.ko2$Taxon_all_unipept))
#all 55 match
length(which(unique(micro.ko2$KEGG) %in% total.list$feature))
#70 kegg match

for(x in 1:nrow(micro.ko2)){
  aa <- which(colnames(total.df)==micro.ko2$Taxon_all_unipept[x])
  bb <- which(rownames(total.df)==micro.ko2$KEGG[x])
  if(length(aa)>0& length(bb)>0){
    total.df[aa,bb] <- 1
    total.df[bb,aa] <- 1
    print(x)
  }
}

# corr1, human, ko, species, microprot

human <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_humanprotein_sample_1385.txt",sep="\t")
human$Protein.Group <- data.frame(str_split_fixed(human$Protein.Group," ",2))$X1
human$Protein.Group <- gsub("[|]",".",human$Protein.Group)
human.match <- human[which(human$Protein.Group %in% total.list$feature),]
#all 31 match

ko <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_microprotein_kegg_sample_1385.txt",sep="\t")
ko.match <- ko[which(ko$KEGGnum %in% total.list$feature),]
#all 96 match

micro <- data.frame(data.table::fread("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_microprotein_sample_1385.txt",sep="\t"))
micro$Protein.Group <- data.frame(str_split_fixed(micro$Protein.Group," ",2))$X1
micro$Protein.Group <- gsub("-",".",micro$Protein.Group)
micro.match <- micro[which(micro$Protein.Group %in% total.list$feature),]
#all 163 match

species <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_filter5_species_sample_1385.txt",sep="\t")
species$Taxon <- gsub("[|]",".",species$Taxon)
species.match <- species[which(species$Taxon %in% total.list$feature),]
#all 55 match

add.corr <- function(x){
  rownames(x) <- x[,1]
  x <- x[,-1]
  x[is.na(x)] <- 0
  corr <- cor(t(x),use="everything",method='spearman')
  corr[corr==1] <- 0
  corr1 <- reshape2::melt(corr)
  corr1$value <- ifelse(abs(corr1$value)>0.8,corr1$value,NA)
  corr1 <- na.omit(corr1)
  return(corr1)
}
corr.human <- add.corr(human.match)
corr.ko <- add.corr(ko.match)
corr.species <- add.corr(species.match)
corr.micro <- add.corr(micro.match)

# corr2, human-ko, huamn-species

add.corr.inter <- function(x,y){
  rownames(x) <- x[,1]
  x <- x[,-1]
  x[is.na(x)] <- 0
  rownames(y) <- y[,1]
  y <- y[,-1]
  y[is.na(y)] <- 0
  corr <- cor(x=t(x),y=t(y),use="everything",method='spearman')
  corr[corr==1] <- 0
  corr1 <- reshape2::melt(corr)
  corr1$value <- ifelse(abs(corr1$value)>0.3,corr1$value,NA)
  corr1 <- na.omit(corr1)
  return(corr1)
}

corr.human.ko <- add.corr.inter(human.match,ko.match)
corr.human.species <- add.corr.inter(human.match,species.match)

totalcor <- rbind(corr.human,corr.ko,corr.species,corr.micro,corr.human.ko,corr.human.species)
totalcor <- na.omit(totalcor)

accur <- data.frame(table(totalcor$Var1))

for(x in 1:nrow(totalcor)){
  aa <- which(colnames(total.df)==totalcor$Var1[x])
  bb <- which(rownames(total.df)==totalcor$Var2[x])
  if(length(aa)>0& length(bb)>0){
    total.df[aa,bb] <- totalcor$value[x]
    total.df[bb,aa] <- totalcor$value[x]
    print(x)
  }
}
##fig6b to gephi for network####
write.csv(total.df,"t2d_network_cor0.3.csv")



#exprot t2d prot. ml + glm# 
prot.link <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/08.link/microbial_prot_to_unique_species_link.txt",sep="\t")
prot.link$Protein <- gsub("-",'.',prot.link$Protein)
prot.link.m <- prot.link[which(prot.link$Protein %in% total.list$feature), ]
write.csv(prot.link.m,"network_prot_to_species.csv",row.names = F)
