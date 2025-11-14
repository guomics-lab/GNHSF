#glm medicine 202509
rm(list=ls())

setwd("/Users/shuangliang/Documents/Guomics/gnhsf/glm_medicine202509")
######################## GLM for med ###########################
sampledf <- read.csv("GNHSF_sample_inform_normalized_1385_2_int.csv")
sampledf$Bristol.scale <- sampledf$Bristol.scale-1
sampledf$Bristol.scale = factor(sampledf$Bristol.scale, levels = c("0","1","2"))

##changemed
sampledf$DLP.med.combine <- ifelse(sampledf$DLP==0,"0",ifelse(sampledf$DLP==1&sampledf$DLP.med==1,"2","1"))
table(sampledf$DLP.med.combine)

sampledf$DM.med.combine <- ifelse(sampledf$DM==0,"0",ifelse(sampledf$DM==1&sampledf$DM.med==1,"2","1"))
table(sampledf$DM.med.combine)

sampledf$HBP.med.combine <- ifelse(sampledf$HBP==0,"0",ifelse(sampledf$HBP==1&sampledf$HBP.med==1,"2","1"))
table(sampledf$HBP.med.combine)

sampledf_tmp <- sampledf[,c(1,2,3,12,33:35)]
sampledf_tmp$Sex <- factor(sampledf_tmp$Sex)
sampledf_tmp$DLP.med.combine <- factor(sampledf_tmp$DLP.med.combine)
sampledf_tmp$DM.med.combine <- factor(sampledf_tmp$DM.med.combine)
sampledf_tmp$HBP.med.combine <- factor(sampledf_tmp$HBP.med.combine)
sampledf_tmp$DLP.weights <- ifelse(sampledf_tmp$DLP.med.combine == 0, 1 / 673,
                                   ifelse(sampledf_tmp$DLP.med.combine == 1, 1 / 156,
                                          1 / 197))
sampledf_tmp$DLP.weights <- sampledf_tmp$DLP.weights / sum(sampledf_tmp$DLP.weights,na.rm = T)

sampledf_tmp$HBP.weights <- ifelse(sampledf_tmp$HBP.med.combine == 0, 1 / 833,
                                   ifelse(sampledf_tmp$HBP.med.combine == 1, 1 / 18,
                                          1 / 434))
sampledf_tmp$HBP.weights <- sampledf_tmp$HBP.weights / sum(sampledf_tmp$HBP.weights,na.rm = T)

sampledf_tmp$DM.weights <- ifelse(sampledf_tmp$DM.med.combine == 0, 1 / 1153,
                                  ifelse(sampledf_tmp$DM.med.combine == 1, 1 / 28,
                                         1 / 126))
sampledf_tmp$DM.weights <- sampledf_tmp$DM.weights / sum(sampledf_tmp$DM.weights,na.rm = T)



func.gml.common <- function(type){
  df <- as.data.frame(read.delim(file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_1385_NA90_INT.tsv"), header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE))
  dft <- as.data.frame(t(df))
  colnames(dft) <- as.character(lapply(strsplit(as.character(colnames(dft)), split=" "), head, n=1))
  if (type == "microprotein"){
    colnames(dft) <- paste0("X", colnames(dft))
  }
  colnames(dft) <- gsub("\\|", ".", colnames(dft))
  colnames(dft) <- gsub("\\-", ".", colnames(dft))
  colnames(dft) <- gsub("\\.$", "", colnames(dft))
  dft$sample <- rownames(dft)
  
  gml.res.final.DLP <- data.frame()
  gml.res.final.HBP <- data.frame()
  gml.res.final.DM <- data.frame()
  for (i in 1:(ncol(dft)-1)){
    cat("\r",type,i)
    df.gml <- merge(dft[,c(i, ncol(dft))], sampledf_tmp, by.x = "sample", by.y = "sample")
    fomular.DLP <- paste0(colnames(df.gml)[2], " ~ Sex + Age + Bristol.scale + DLP.med.combine")
    gml.DLP <- glm(fomular.DLP, data = df.gml, family = gaussian,weights=DLP.weights)
    gml.sum.DLP <- summary(gml.DLP)
    gml.res.DLP <- as.data.frame(gml.sum.DLP$coefficients)
    gml.res.DLP$feature <- rep(colnames(df.gml)[2], nrow(gml.res.DLP))
    gml.res.DLP$metadata <- rownames(gml.res.DLP)
    gml.res.row <- c()
    gml.res.add.DLP <- data.frame()
    if (length(c(grep("DLP", rownames(gml.res.DLP)))) > 0){
      gml.res.row <- c(grep("DLP", rownames(gml.res.DLP)))
      gml.res.add.DLP <- gml.res.DLP[gml.res.row,]
      gml.res.add.DLP$level <- type
      gml.res.add.DLP$clino <- "DLP"
      
      gml.res.final.DLP <- rbind(gml.res.final.DLP, gml.res.add.DLP)
    }
    
    fomular.HBP <- paste0(colnames(df.gml)[2], " ~ Sex + Age + Bristol.scale + HBP.med.combine")
    gml.HBP <- glm(fomular.HBP, data = df.gml, family = gaussian,weights=HBP.weights)
    gml.sum.HBP <- summary(gml.HBP)
    gml.res.HBP <- as.data.frame(gml.sum.HBP$coefficients)
    gml.res.HBP$feature <- rep(colnames(df.gml)[2], nrow(gml.res.HBP))
    gml.res.HBP$metadata <- rownames(gml.res.HBP)
    gml.res.row <- c()
    gml.res.add.HBP <- data.frame()
    if (length(c(grep("HBP", rownames(gml.res.HBP)))) > 0){
      gml.res.row <- c(grep("HBP", rownames(gml.res.HBP)))
      gml.res.add.HBP <- gml.res.HBP[gml.res.row,]
      gml.res.add.HBP$level <- type
      gml.res.add.HBP$clino <- "HBP"
      
      gml.res.final.HBP <- rbind(gml.res.final.HBP, gml.res.add.HBP)
      
    }
    
    fomular.DM <- paste0(colnames(df.gml)[2], " ~ Sex + Age + Bristol.scale + DM.med.combine")
    gml.DM <- glm(fomular.DM, data = df.gml, family = gaussian,weights=DM.weights)
    gml.sum.DM <- summary(gml.DM)
    gml.res.DM <- as.data.frame(gml.sum.DM$coefficients)
    gml.res.DM$feature <- rep(colnames(df.gml)[2], nrow(gml.res.DM))
    gml.res.DM$metadata <- rownames(gml.res.DM)
    gml.res.row <- c()
    gml.res.add.DM <- data.frame()
    if (length(c(grep("DM", rownames(gml.res.DM)))) > 0){
      gml.res.row <- c(grep("DM", rownames(gml.res.DM)))
      gml.res.add.DM <- gml.res.DM[gml.res.row,]
      gml.res.add.DM$clino <- "T2D"
      gml.res.add.DM$level <- type
      gml.res.final.DM <- rbind(gml.res.final.DM, gml.res.add.DM)
    }
  }
  gml.res.final <- rbind(gml.res.final.DM,gml.res.final.HBP,gml.res.final.DLP)
  
  
  ##need change#####
  gml.res.final1 <- do.call(rbind,lapply(split(gml.res.final,gml.res.final$level),function(x){
    x$q <- p.adjust(x$`Pr(>|t|)`,method="BH")
    q.sig <- c()
    for (i in 1:nrow(x)) {
      if (x$q[i] < 0.001){
        q.sig[i] <- "***"
      }else if (x$q[i] < 0.01){
        q.sig[i] <- "**"
      }else if (x$q[i] < 0.05){
        q.sig[i] <- "*"
      }else  {
        q.sig[i] <- "nsig"
      }
    }
    x$q.sig <- q.sig
    return(x)
  }
  ))
  write.csv(gml.res.final1, file = paste0("med.glm.", type, "20250423.csv"))
}





func.gml.common("cog")
func.gml.common("kegg")
func.gml.common("genus")
func.gml.common("species")
func.gml.common("humanprotein")
func.gml.common("microprotein")

gml.human <- read.delim("med.glm.humanprotein20250423.csv",sep=",")
human.name <- data.frame(stringr::str_split_fixed(gml.human$feature,"_|[.]",10))
gml.human$feature <- human.name$X4
write.csv(gml.human,"med.glm.humanprotein20250423.csv",row.names = F)


all <- do.call(rbind, lapply(c("genus","species","cog","kegg","humanprotein","microprotein"), function(x) {
  data=read.csv(paste0("med.glm.", x, "20250423.csv"))
  return(data)
}))


write.csv(all,"med.glm.all.csv")


###plot####
rm(list=ls())
all <- read.csv("med.glm.all.csv")
all <- all[,-1]

all.sig <- subset(all,q.sig!="nsig")

stats <- data.frame(table(all.sig$metadata,all.sig$level))

library(openxlsx)
anno.human <- read.xlsx("~/Documents/gnhsf/data202509/metadata/uniprotkb_human.xlsx")
colnames(anno.human)[1] <- "protein"
anno.kegg <- read.xlsx("~/Documents/gnhsf/data202509/metadata/kegg_database20230328.xlsx")
anno.kegg0 <- anno.kegg
anno.kegg <- unique(anno.kegg[,c(5,4)])
anno.cog <- read.xlsx("~/Documents/gnhsf/data202509/metadata/01.NCBI.COG.list.xlsx")
anno.micro <- read.xlsx("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microprot_2phylum2genus2species_all20250916.xlsx")
anno.micro$glm.match <- paste0("X",anno.micro$Protein)
anno.micro$glm.match <- gsub("-",".",anno.micro$glm.match)




library(reshape2)
all.sig.1 <- dcast(all.sig,feature+level~metadata,value.var = "Estimate")



#####dlp#####

dlp <- all.sig.1[,1:4]
dlp <- dlp[-which(is.na(dlp$DLP.med.combine1)& is.na(dlp$DLP.med.combine2)),]

library(ggvenn)
list1 <- list("disease"=dlp$feature[which(!is.na(dlp$DLP.med.combine1))],"med"=dlp$feature[which(!is.na(dlp$DLP.med.combine2))])
ggvenn(list1)

dlp2 <- dlp[-c(which(dlp$DLP.med.combine1>0 & dlp$DLP.med.combine2>0),which(dlp$DLP.med.combine1<0 & dlp$DLP.med.combine2<0)),]
dlp3 <- dlp[which(!is.na(dlp$DLP.med.combine1)),]
dlp3$change <- ifelse(dlp3$DLP.med.combine1<0,"down.up","up.down")

dlp.human <- subset(dlp3,level=="humanprotein")
library(dplyr)
dlp.human <- merge(dlp.human,anno.human,by.x="feature",by.y="protein",all.x=T)
dlp.human <- arrange(dlp.human,DLP.med.combine1)

dlp.kegg <- subset(dlp3,level=="kegg")
library(dplyr)
dlp.kegg <- merge(dlp.kegg,anno.kegg,by.x="feature",by.y="kegg",all.x=T)
dlp.kegg <- arrange(dlp.kegg,DLP.med.combine1)

dlp.cog <- subset(dlp3,level=="cog")
library(dplyr)
dlp.cog <- merge(dlp.cog,anno.cog,by.x="feature",by.y="COGnum")
dlp.cog <- arrange(dlp.cog,DLP.med.combine1)

dlp.micro <- subset(dlp3,level=="microprotein")
library(dplyr)
dlp.micro <- merge(dlp.micro,anno.micro,by.x="feature",by.y="glm.match")
dlp.micro <- arrange(dlp.micro,DLP.med.combine1)

dlp.genus <- subset(dlp3,level=="genus")
dlp.species <- subset(dlp3,level=="species")


dlp4 <- dlp2[-which(!is.na(dlp2$DLP.med.combine1)),]
dlp4$change <- ifelse(dlp4$DLP.med.combine2<0,"nsig.down","nsig.up")

dlp2.human <- subset(dlp4,level=="humanprotein")
library(dplyr)
dlp2.human <- merge(dlp2.human,anno.human,by.x="feature",by.y="protein")
dlp2.human <- arrange(dlp2.human,DLP.med.combine1)

dlp2.kegg <- subset(dlp4,level=="kegg")
library(dplyr)
dlp2.kegg <- merge(dlp2.kegg,anno.kegg,by.x="feature",by.y="kegg")
dlp2.kegg <- arrange(dlp2.kegg,DLP.med.combine1)

dlp2.cog <- subset(dlp4,level=="cog")
library(dplyr)
dlp2.cog <- merge(dlp2.cog,anno.cog,by.x="feature",by.y="COGnum")
dlp2.cog <- arrange(dlp2.cog,DLP.med.combine1)

dlp2.micro <- subset(dlp4,level=="microprotein")
library(dplyr)
dlp2.micro <- merge(dlp2.micro,anno.micro,by.x="feature",by.y="glm.match")
dlp2.micro <- arrange(dlp2.micro,DLP.med.combine1)

dlp2.genus <- subset(dlp4,level=="genus")
dlp2.species <- subset(dlp4,level=="species")




####dm###. IG, myrosin,ras ,acyl-CoA dehydrogenase,collagen，Cytochrome oxidase subunit ,Cytokeratin
##	Large ribosomal subunit protein up, small up
DM <- all.sig.1[,c(1:2,5,6)]
DM <- DM[-which(is.na(DM$DM.med.combine1)& is.na(DM$DM.med.combine2)),]

library(ggvenn)
list1 <- list("disease"=DM$feature[which(!is.na(DM$DM.med.combine1))],"med"=DM$feature[which(!is.na(DM$DM.med.combine2))])
ggvenn(list1)

DM2 <- DM[-c(which(DM$DM.med.combine1>0 & DM$DM.med.combine2>0),which(DM$DM.med.combine1<0 & DM$DM.med.combine2<0)),]
DM3 <- DM2[which(!is.na(DM2$DM.med.combine1)),]
DM3$change <- ifelse(DM3$DM.med.combine1<0,"down.up","up.down")

DM.human <- subset(DM3,level=="humanprotein")
library(dplyr)
DM.human <- merge(DM.human,anno.human,by.x="feature",by.y="protein",all.x=T)
DM.human <- arrange(DM.human,DM.med.combine1)

DM.kegg <- subset(DM3,level=="kegg")
library(dplyr)
DM.kegg <- merge(DM.kegg,anno.kegg,by.x="feature",by.y="kegg")
DM.kegg <- arrange(DM.kegg,DM.med.combine1)

DM.cog <- subset(DM3,level=="cog")
library(dplyr)
DM.cog <- merge(DM.cog,anno.cog,by.x="feature",by.y="COGnum")
DM.cog <- arrange(DM.cog,DM.med.combine1)

DM.micro <- subset(DM3,level=="microprotein")
library(dplyr)
DM.micro <- merge(DM.micro,anno.micro,by.x="feature",by.y="glm.match")
DM.micro <- arrange(DM.micro,DM.med.combine1)

DM.genus <- subset(DM3,level=="genus")
DM.species <- subset(DM3,level=="species")

DM4 <- DM2[-which(!is.na(DM2$DM.med.combine1)),]
DM4$change <- ifelse(DM4$DM.med.combine2<0,"nsig.down","nsig.up")

DM2.human <- subset(DM4,level=="humanprotein")
library(dplyr)
DM2.human <- merge(DM2.human,anno.human,by.x="feature",by.y="protein",all.x=T)
DM2.human <- arrange(DM2.human,DM.med.combine1)

DM2.kegg <- subset(DM4,level=="kegg")
library(dplyr)
DM2.kegg <- merge(DM2.kegg,anno.kegg,by.x="feature",by.y="kegg")
DM2.kegg <- arrange(DM2.kegg,DM.med.combine1)

DM2.cog <- subset(DM4,level=="cog")
library(dplyr)
DM2.cog <- merge(DM2.cog,anno.cog,by.x="feature",by.y="COGnum")
DM2.cog <- arrange(DM2.cog,DM.med.combine1)

DM2.micro <- subset(DM4,level=="microprotein")
library(dplyr)
DM2.micro <- merge(DM2.micro,anno.micro,by.x="feature",by.y="glm.match")
DM2.micro <- arrange(DM2.micro,DM.med.combine1)

DM2.genus <- subset(DM4,level=="genus")
DM2.species <- subset(DM4,level=="species")




####HBP### Cytokeratin,ras,myosin
## ribosomal subunit protein  down
HBP <- all.sig.1[,c(1:2,7,8)]
HBP <- HBP[-which(is.na(HBP$HBP.med.combine1)& is.na(HBP$HBP.med.combine2)),]

library(ggvenn)
list1 <- list("disease"=HBP$feature[which(!is.na(HBP$HBP.med.combine1))],"med"=HBP$feature[which(!is.na(HBP$HBP.med.combine2))])
ggvenn(list1)

HBP2 <- HBP[-c(which(HBP$HBP.med.combine1>0 & HBP$HBP.med.combine2>0),which(HBP$HBP.med.combine1<0 & HBP$HBP.med.combine2<0)),]
HBP3 <- HBP2[which(!is.na(HBP2$HBP.med.combine1)),]
HBP3$change <- ifelse(HBP3$HBP.med.combine1<0,"down.up","up.down")

HBP.human <- subset(HBP3,level=="humanprotein")
library(dplyr)
HBP.human <- merge(HBP.human,anno.human,by.x="feature",by.y="protein",all.x=T)
HBP.human <- arrange(HBP.human,HBP.med.combine1)

HBP.kegg <- subset(HBP3,level=="kegg")
library(dplyr)
HBP.kegg <- merge(HBP.kegg,anno.kegg,by.x="feature",by.y="kegg")
HBP.kegg <- arrange(HBP.kegg,HBP.med.combine1)

HBP.cog <- subset(HBP3,level=="cog")
library(dplyr)
HBP.cog <- merge(HBP.cog,anno.cog,by.x="feature",by.y="COGnum")
HBP.cog <- arrange(HBP.cog,HBP.med.combine1)

HBP.micro <- subset(HBP3,level=="microprotein")
library(dplyr)
HBP.micro <- merge(HBP.micro,anno.micro,by.x="feature",by.y="glm.match")
HBP.micro <- arrange(HBP.micro,HBP.med.combine1)

HBP.genus <- subset(HBP3,level=="genus")
HBP.species <- subset(HBP3,level=="species")

HBP4 <- HBP2[-which(!is.na(HBP2$HBP.med.combine1)),]
HBP4$change <- ifelse(HBP4$HBP.med.combine2<0,"nsig.down","nsig.up")

HBP2.human <- subset(HBP4,level=="humanprotein")
library(dplyr)
HBP2.human <- merge(HBP2.human,anno.human,by.x="feature",by.y="protein",all.x =T)
HBP2.human <- arrange(HBP2.human,HBP.med.combine1)

HBP2.kegg <- subset(HBP4,level=="kegg")
library(dplyr)
HBP2.kegg <- merge(HBP2.kegg,anno.kegg,by.x="feature",by.y="kegg")
HBP2.kegg <- arrange(HBP2.kegg,HBP.med.combine1)

HBP2.cog <- subset(HBP4,level=="cog")
library(dplyr)
HBP2.cog <- merge(HBP2.cog,anno.cog,by.x="feature",by.y="COGnum")
HBP2.cog <- arrange(HBP2.cog,HBP.med.combine1)

HBP2.micro <- subset(HBP4,level=="microprotein")
library(dplyr)
HBP2.micro <- merge(HBP2.micro,anno.micro,by.x="feature",by.y="glm.match")
HBP2.micro <- arrange(HBP2.micro,HBP.med.combine1)

HBP2.genus <- subset(HBP4,level=="genus")
HBP2.species <- subset(HBP4,level=="species")


###human sig fig5 D-G####
human.sig <- subset(all.sig.1,level=="humanprotein")
human.sig <- merge(human.sig,anno.human,by.x="feature",by.y="protein")


micro.sig <- subset(all.sig.1,level=="microprotein")
micro.sig <- merge(micro.sig,anno.micro,by.x="feature",by.y="glm.match")

species.sig <- subset(all.sig.1,level=="species")

genus.sig <- subset(all.sig.1,level=="genus")


human.sig.list <- unique(c(dlp.human$feature,DM.human$feature,HBP.human$feature))
all.human <- subset(all, feature %in% human.sig.list)
all.human <- merge(all.human,anno.human,by.x="feature",by.y="protein")
all.human$Entry.Name <- gsub("_HUMAN","",all.human$Entry.Name)
all.human$des <- paste0(all.human$Entry.Name,": ",all.human$Protein.names)
library(reshape2)

####myosin 
myolist <- unique(c(dlp.human$feature[grep("myosin|Myosin",dlp.human$Protein.names)],DM.human$feature[grep("myosin|Myosin",DM.human$Protein.names)],
                    HBP.human$feature[grep("myosin|Myosin",HBP.human$Protein.names)]))
myosin.list <- all.human[which(all.human$feature %in% myolist),]
group <- dcast(myosin.list,feature+des~metadata,value.var ="Estimate")
group.sig <- dcast(myosin.list,feature+des~metadata,value.var ="q.sig")
group.sig[group.sig=="nsig"] <- ""

rownames(group) <- group$des

group <- group[,c(3:8)]
rownames(group.sig) <- group.sig$des
group.sig <- group.sig[,c(3:8)]

library(ComplexHeatmap)
breaks <- c(seq(-0.4, 0, length.out = 6), seq(0, 0.4, length.out = 6)) %>% unique()
colors <- c(
  colorRampPalette(c("cyan3", "lightblue"))(5),  
  colorRampPalette(c("thistle2","lightcoral"))(5)  
)
pdf("myosin.pdf")
pheatmap(group,display_numbers = as.matrix(group.sig),cluster_cols = F,breaks = breaks,color=colors)
dev.off()


####immuno
library(reshape2)
immunolist <- unique(c(DM.human$feature[grep("Immun|immun",DM.human$Protein.names)]))
immuno.list <- all.human[which(all.human$feature %in% immunolist),]
group <- dcast(immuno.list,feature+des~metadata,value.var ="Estimate")
group.sig <- dcast(immuno.list,feature+des~metadata,value.var ="q.sig")
group.sig[group.sig=="nsig"] <- ""

rownames(group) <- group$des
group <- group[,c(5:6)]
rownames(group.sig) <- group$des
group.sig <- group.sig[,5:6]

library(ComplexHeatmap)
breaks <- c(seq(-0.4, 0, length.out = 6), seq(0, 0.4, length.out = 6)) %>% unique()
colors <- c(
  colorRampPalette(c("cyan3", "lightblue"))(5),  
  colorRampPalette(c("thistle2","lightcoral"))(5)   
)
pdf("immuno.pdf")
pheatmap(group,display_numbers = as.matrix(group.sig),cluster_cols = F,breaks = breaks,color=colors)
dev.off()


##Cytochrome 
library(reshape2)
cytolist <- unique(c(DM.human$feature[grep("Cytochrome|cytochrome",DM.human$Protein.names)]))
cyto.list <- all.human[which(all.human$feature %in% cytolist),]
group <- dcast(cyto.list,feature+des~metadata,value.var ="Estimate")
group.sig <- dcast(cyto.list,feature+des~metadata,value.var ="q.sig")
group.sig[group.sig=="nsig"] <- ""

rownames(group) <- group$des
group <- group[,c(5:6)]
rownames(group.sig) <- group$des
group.sig <- group.sig[,5:6]

library(ComplexHeatmap)
breaks <- c(seq(-0.4, 0, length.out = 6), seq(0, 0.4, length.out = 6)) %>% unique()
colors <- c(
  colorRampPalette(c("cyan3", "lightblue"))(5),  
  colorRampPalette(c("thistle2","lightcoral"))(5)  
)
pdf("cytochrome.pdf")
pheatmap(group,display_numbers = as.matrix(group.sig),cluster_cols = F,breaks = breaks,color=colors)
dev.off()

####ribo
ribolist <- unique(c(DM.human$feature[grep("riboso|Riboso",DM.human$Protein.names)],DM2.human$feature[grep("riboso|Riboso",DM2.human$Protein.names)],
                     HBP.human$feature[grep("riboso|Riboso",HBP.human$Protein.names)],HBP2.human$feature[grep("riboso|Riboso",HBP2.human$Protein.names)]))
ribo.list <- all.human[which(all.human$feature %in% ribolist),]
group <- dcast(ribo.list,feature+des~metadata,value.var ="Estimate")
group.sig <- dcast(ribo.list,feature+des~metadata,value.var ="q.sig")
group.sig[group.sig=="nsig"] <- ""
rownames(group) <- group$des
group <- group[,c(5:8)]
rownames(group.sig) <- group$des
group.sig <- group.sig[,5:8]

library(ComplexHeatmap)
breaks <- c(seq(-0.4, 0, length.out = 6), seq(0, 0.4, length.out = 6)) %>% unique()
colors <- c(
  colorRampPalette(c("cyan3", "lightblue"))(5), 
  colorRampPalette(c("thistle2","lightcoral"))(5)  
)
pdf("riboso.pdf")
pheatmap(group,display_numbers = as.matrix(group.sig),cluster_cols = F,breaks = breaks,color=colors)
dev.off()


#####microbe##########

####for microbial protein analysis， medication########
library(stringr)
dlp.micro <- subset(dlp3,level=="microprotein")
library(dplyr)
dlp.micro <- merge(dlp.micro,anno.micro,by.x="feature",by.y="glm.match")
dlp.micro <- arrange(dlp.micro,DLP.med.combine1)

library(tidyr)
dlp.micro2 <-dlp.micro %>% as_tibble() %>%    separate_rows(Taxon_rank_unipept,Taxon_name_unipept, sep = ";")


dlp.micro2 <- dlp.micro2[which(dlp.micro2$Taxon_rank_unipept %in% c("phylum","genus","species")),]
library(reshape2)
dlp.micro3 <- dcast(dlp.micro2,feature+Protein+Eggpro_des+Eggpro_prefer_name+COG+COG_NCBIdes+KEGG+KEGG_NCBIcat_des+DLP.med.combine1+DLP.med.combine2+change~Taxon_rank_unipept,value.var = c("Taxon_name_unipept"))
dlp.micro3[dlp.micro3=="NA"] <- NA
dlp.pro.table <- data.frame(table(dlp.micro3$species))

DM.micro <- subset(DM3,level=="microprotein")
library(dplyr)
DM.micro <- merge(DM.micro,anno.micro,by.x="feature",by.y="glm.match")
DM.micro <- arrange(DM.micro,DM.med.combine1)
library(tidyr)
DM.micro2 <-DM.micro %>% as_tibble() %>%    separate_rows(Taxon_rank_unipept,Taxon_name_unipept, sep = ";")
DM.micro2 <- DM.micro2[which(DM.micro2$Taxon_rank_unipept %in% c("phylum","genus","species")),]
library(reshape2)
DM.micro3 <- dcast(DM.micro2,feature+Protein+Eggpro_des+Eggpro_prefer_name+COG+COG_NCBIdes+KEGG+KEGG_NCBIcat_des+DM.med.combine1+DM.med.combine2+change~Taxon_rank_unipept,value.var = c("Taxon_name_unipept"))
DM.micro3[DM.micro3=="NA"] <- NA
dm.pro.table <- data.frame(table(DM.micro3$species,DM.micro3$change))




HBP.micro <- subset(HBP3,level=="microprotein")
library(dplyr)
HBP.micro <- merge(HBP.micro,anno.micro,by.x="feature",by.y="glm.match")
HBP.micro <- arrange(HBP.micro,HBP.med.combine1)
library(tidyr)
HBP.micro2 <-HBP.micro %>% as_tibble() %>%    separate_rows(Taxon_rank_unipept,Taxon_name_unipept, sep = ";")
HBP.micro2 <- HBP.micro2[which(HBP.micro2$Taxon_rank_unipept %in% c("phylum","genus","species")),]
library(reshape2)
HBP.micro3 <- dcast(HBP.micro2,feature+Protein+Eggpro_des+Eggpro_prefer_name+COG+COG_NCBIdes+KEGG+KEGG_NCBIcat_des+HBP.med.combine1+HBP.med.combine2+change~Taxon_rank_unipept,value.var = c("Taxon_name_unipept"))
HBP.micro3[HBP.micro3=="NA"] <- NA
HBP.pro.table <- data.frame(table(HBP.micro3$species,HBP.micro3$change))



####microbe protein fig5b####
library(ggplot2)
library(ggthemes)
dm.pro.table <- data.frame(table(DM.micro3$species))
dm.pro.table <- arrange(dm.pro.table,desc(Freq))
dm.pro.table1 <- data.frame(table(DM.micro3$species,DM.micro3$change))
dm.pro.table1 <- dm.pro.table1[which(dm.pro.table1$Var1 %in% dm.pro.table$Var1[1:5]),]
dm.pro.table1$Var1 <- factor(dm.pro.table1$Var1,levels=rev(as.character(dm.pro.table$Var1[1:5])))
ggplot(dm.pro.table1,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat="identity")+coord_flip()+scale_fill_manual(values = c("#c3dfe7","#f1a7a1"))+ggtitle("T2D")+theme_base()
ggsave("dm_med_micro_prot_stats.pdf",width=7,height=7)

HBP.pro.table <- data.frame(table(HBP.micro3$species))
HBP.pro.table <- arrange(HBP.pro.table,desc(Freq))
HBP.pro.table1 <- data.frame(table(HBP.micro3$species,HBP.micro3$change))
HBP.pro.table1 <- HBP.pro.table1[which(HBP.pro.table1$Var1 %in% HBP.pro.table$Var1[1:5]),]
HBP.pro.table1$Var1 <- factor(HBP.pro.table1$Var1,levels=rev(as.character(HBP.pro.table$Var1[1:5])))
ggplot(HBP.pro.table1,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat="identity")+coord_flip()+scale_fill_manual(values = c("#c3dfe7","#f1a7a1"))+ggtitle("HBP")+theme_base()
ggsave("hbp_med_micro_prot_stats.pdf",width=7,height=7)


dlp.pro.table <- data.frame(table(dlp.micro3$species))
dlp.pro.table <- arrange(dlp.pro.table,desc(Freq))
dlp.pro.table1 <- data.frame(table(dlp.micro3$species,dlp.micro3$change))
dlp.pro.table1 <- dlp.pro.table1[which(dlp.pro.table1$Var1 %in% dlp.pro.table$Var1[1:5]),]
dlp.pro.table1$Var1 <- factor(dlp.pro.table1$Var1,levels=rev(as.character(dlp.pro.table$Var1[1:5])))
ggplot(dlp.pro.table1,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat="identity")+coord_flip()+scale_fill_manual(values = c("#c3dfe7","#f1a7a1"))+ggtitle("DLP")+theme_base()
ggsave("dlp_med_micro_prot_stats.pdf",width=7,height=7)



###fig5c#####
library(dplyr)
library(MicrobiomeProfiler)
library(clusterProfiler)

DM.micro3$metadata <- "DM"
colnames(DM.micro3)[9:10] <-c( "med0","med1")
dlp.micro3$metadata <- "dlp"
colnames(dlp.micro3)[9:10] <-c( "med0","med1")
HBP.micro3$metadata <- "HBP"
colnames(HBP.micro3)[9:10] <-c( "med0","med1")
micro <- rbind(DM.micro3,dlp.micro3,HBP.micro3)
micro <- micro[-which(is.na(micro$species)),]
micro$metadata <- factor(micro$metadata,levels=c("DM","HBP","dlp"))

######add anno kegg########
anno.kegg1 <- anno.kegg0
anno.kegg1 <- unique(anno.kegg1[,1:3])
anno.kegg1$C <- substr(anno.kegg1$C,7,nchar(anno.kegg1$C))
anno.kegg1$C <- data.frame(str_split_fixed(anno.kegg1$C," \\[",2))$X1


target="s__Faecalibacterium_prausnitzii"
micro.target <- subset(micro,species==target)
micro.target <- micro.target[-which(is.na(micro.target$KEGG)),]
target.enrich <- compareCluster(KEGG~metadata+change,data=micro.target,fun="enrichKO",source_from = "MicrobiomeProfiler",qvalueCutoff = 0.1,pvalueCutoff=0.05,maxGSSize=230)
compare <- target.enrich@compareClusterResult
compare <- merge(compare,anno.kegg1,by.x="Description",by.y="C",all.x=T)
write.csv(compare,"compare.test.fae.csv")
# mannual impute KEGG A and B

compare <- read.xlsx("compare.test.fae_1.xlsx")
compare <- compare[,-1]
compare <- arrange(compare,B)
compare$Description <- factor(compare$Description,levels=rev(unique(compare$Description)))
target.enrich@compareClusterResult <- compare

library(RColorBrewer)
color1 <- palette("Paired")
color2 <- palette("Pastel1")
color <- c(color2,color1[-1])
p <- dotplot(target.enrich, includeAll=TRUE, size="RichFactor",showCategory=50) 
p_colorbar <- ggplot(unique(compare[,c(1,17)]), aes(x = 1, y = Description, fill = B)) +
  geom_tile(color = "white", width = 0.1) +
  scale_fill_manual(values = color ) +
  theme_void() + theme(legend.position = "left")   

library(patchwork)

final_plot <- p_colorbar + p + plot_layout(widths = c(0.1, 1))
final_plot
ggsave("Faecalibacterium_med.pdf",width=10,height=10)




#####mega 2######
library(dplyr)
DM2.micro <- subset(DM4,level=="microprotein")
library(dplyr)
DM2.micro <- merge(DM2.micro,anno.micro,by.x="feature",by.y="glm.match")
DM2.micro <- arrange(DM2.micro,DM.med.combine1)

library(tidyr)
DM2.micro2 <-DM2.micro %>% as_tibble() %>%    separate_rows(Taxon_rank_unipept,Taxon_name_unipept, sep = ";")

DM2.micro2 <- DM2.micro2[which(DM2.micro2$Taxon_rank_unipept %in% c("phylum","genus","species")),]
library(reshape2)
DM2.micro3 <- dcast(DM2.micro2,feature+Protein+Eggpro_des+Eggpro_prefer_name+COG+COG_NCBIdes+KEGG+KEGG_NCBIcat_des+DM.med.combine1+DM.med.combine2+change~Taxon_rank_unipept,value.var = "Taxon_name_unipept")

DM2.micro3$species <- gsub("s__","",DM2.micro3$species)
DM2.pro.table <- data.frame(table(DM2.micro3$species))
DM2.pro.table <- DM2.pro.table[-which(DM2.pro.table$Var1=="NA"),]
DM2.pro.table <- arrange(DM2.pro.table,desc(Freq))
DM2.pro.table$species <- paste0(DM2.pro.table$Var1,"(",DM2.pro.table$Freq,")")
DM2.pro.table1 <- data.frame(table(DM2.micro3$species,DM2.micro3$change))
DM2.pro.table1 <- DM2.pro.table1[-which(DM2.pro.table1$Var1=="NA"),]
DM2.pro.table1 <- DM2.pro.table1[which(DM2.pro.table1$Var1 %in% DM2.pro.table$Var1[1:40]),]
DM2.pro.table1$Var1 <- factor(DM2.pro.table1$Var1,levels=rev(as.character(DM2.pro.table$Var1[1:40])))
ggplot(DM2.pro.table1,aes(x=Var1,y=Freq,fill=Var2))+geom_bar(stat="identity")+coord_flip()+scale_fill_manual(values = c("#c3dfe7","#f1a7a1"))+theme_base()
ggsave("DM2_med_micro_prot_stats.pdf",width=7,height=7)


dm2.micro4 <- DM2.micro3[-which(DM2.micro3$species=="NA"),]
dm2.micro4 <- merge(dm2.micro4,DM2.pro.table,by.x="species",by.y="Var1",all.x=T)
dm2.micro4$species <- factor(dm2.micro4$species,levels=rev(as.character(DM2.pro.table$Var1[1:40])))
dm2.micro4$species.y <- factor(dm2.micro4$species.y,levels=rev(as.character(DM2.pro.table$species[1:40])))

dm2.micro4 <- subset(dm2.micro4,species %in% DM2.pro.table$Var1[1:40])

table(dm2.micro4$species)

##figs11 b####
library(ggplot2)
ggplot(dm2.micro4,aes(x=species.y,y=DM.med.combine2))+geom_boxplot()+coord_flip()+theme_base()
ggsave("DM2_med_micro_prot_estimate_number20250512.pdf",width=7,height=7)


###stats total#####
dlp2$meta <- "Dyslipidemia"
dlp2$group <- ifelse(is.na(dlp2$DLP.med.combine1),"medication specific features","medication associated features")
colnames(dlp2)[3:4] <- c("med.combine1","med.combine2")

DM2$meta <- "Diabetes"
DM2$group <- ifelse(is.na(DM2$DM.med.combine1),"medication specific features","medication associated features")
colnames(DM2)[3:4] <- c("med.combine1","med.combine2")

HBP2$meta <- "Hypertension"
HBP2$group <- ifelse(is.na(HBP2$HBP.med.combine1),"medication specific features","medication associated features")
colnames(HBP2)[3:4] <- c("med.combine1","med.combine2")

stats.total <- rbind(dlp2,DM2,HBP2)  
write.table(stats.total,"med_mets_total.txt",sep="\t")  

table(stats.total$group)





tables8 <- read.delim("med_mets_total.txt",sep="\t")
tables8$meta <- ifelse(tables8$meta=="Diabetes","Type 2 Diabetes",tables8$meta)
table(tables8$group)
library(reshape2)
all2 <- all
all2$metadata2 <- sapply(all2$metadata, function(x) ifelse(length(grep("combine1",x))>0,"non-disease_vs_unmedicated disease","non-disease_vs_medicated disease"))
 
all2.estimate <-  dcast(all2,feature+clino+level~metadata2,value.var="Estimate")
colnames(all2.estimate)[4:5] <- paste0("Beta.coefficient.",colnames(all2.estimate)[4:5])
all2.q <-  dcast(all2,feature+clino+level~metadata2,value.var="q")
colnames(all2.q)[4:5] <- paste0("Q.",colnames(all2.q)[4:5])
all2.qsig <-  dcast(all2,feature+clino+level~metadata2,value.var="q.sig")
colnames(all2.qsig)[4:5] <- paste0("Q.sig.",colnames(all2.qsig)[4:5])

all2.all <- Reduce(function(x,y) merge(x,y,by=c("feature","clino","level"),all.x=TRUE),list(all2.estimate,all2.q,all2.qsig),accumulate =FALSE)
all2.all$clino <- ifelse(all2.all$clino=="T2D","Type 2 Diabetes",ifelse(all2.all$clino=="DLP","Dyslipidemia","Hypertension"))
all2.all <- merge(all2.all,tables8[,c(2,11,12)],by.x=c("feature","clino"),by.y=c("feature_name","meta"),all.y=T)
all2.all <- all2.all[,c(3,1,2,5,7,9,4,6,8,10)]
all2.all$annotation <- NA
all2.all$annotation <- ifelse(all2.all$group=="medication specific features" & all2.all$`Beta.coefficient.non-disease_vs_medicated disease`>0,
                              "Disease-Null, Treatment-Positive Associations", 
                              ifelse(all2.all$group=="medication specific features" & all2.all$`Beta.coefficient.non-disease_vs_medicated disease`<0,
                                     "Disease-Null, Treatment-Negative Associations",NA))

all2.all$annotation <- ifelse(all2.all$group=="medication associated features" &all2.all$`Beta.coefficient.non-disease_vs_unmedicated disease`>0 &all2.all$`Q.sig.non-disease_vs_medicated disease`=="nsig",
                              "Disease-Positive, Treatment-Negative/Null Associations", 
                              ifelse(all2.all$group=="medication associated features" &all2.all$`Beta.coefficient.non-disease_vs_unmedicated disease`<0 &all2.all$`Q.sig.non-disease_vs_medicated disease`=="nsig",
                                     "Disease-Negative, Treatment-Positive/Null Associations",all2.all$annotation))

                              
all2.all$annotation <- ifelse(all2.all$group=="medication associated features" &all2.all$`Beta.coefficient.non-disease_vs_unmedicated disease`>0 &all2.all$`Q.sig.non-disease_vs_medicated disease`!="nsig" & all2.all$`Beta.coefficient.non-disease_vs_medicated disease`<0,
                              "Disease-Positive, Treatment-Negative/Null Associations", 
                              ifelse(all2.all$group=="medication associated features" &all2.all$`Beta.coefficient.non-disease_vs_unmedicated disease`<0 &all2.all$`Q.sig.non-disease_vs_medicated disease`!="nsig" & all2.all$`Beta.coefficient.non-disease_vs_medicated disease`>0,
                                     "Disease-Negative, Treatment-Positive/Null Associations", all2.all$annotation))
table(all2.all$annotation)
library(openxlsx)

all2.all$feature <- gsub("^X","",all2.all$feature)

micro.link <- read.xlsx("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microprot_all_arrange20250916.xlsx")
micro.link$Protein <- gsub("-",".",micro.link$Protein)

all2.all <- merge(all2.all,micro.link,by.x="feature",by.y="Protein",all.x=T)

write.xlsx(all2.all,"tableS8202510.xlsx")                                     
                                     
                                     

                                     