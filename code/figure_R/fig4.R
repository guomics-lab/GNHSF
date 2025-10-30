rm(list=ls())
gc()

setwd("~/Documents/gnhsf/data202509/code/fig4_shared_met")

library(openxlsx)
###annotate
annotate <- read.xlsx("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microprot_2phylum2genus2species_all20250916.xlsx")
annotate$glm.match <- paste0("X",annotate$Protein)
annotate$glm.match <- gsub("-",".",annotate$glm.match)
annotate$glm.match <- gsub("[.]$","",annotate$glm.match)

cog <- read.xlsx("~/Documents/gnhsf/data202509/metadata/01.NCBI.COG.list.xlsx")
kegg <- read.xlsx("~/Documents/gnhsf/data202509/metadata/kegg_database20230328.xlsx")
human <- read.delim("~/Documents/gnhsf/data202509/metadata/gnhsf_human20250827.txt",sep="\t")
human <- unique(human[,c(1,5,6,7,2,3)])
colnames(human) <- c("prot","name","des","ko","a","b")
human2 <- aggregate(human,by=list(human$prot),function(x) paste0(na.omit(unique(x)),collapse =";"))
human2 <- human2[,-1]
human2$name <- paste(human2$prot,human2$name,sep="_")

kegg <- kegg[,c(6,5,4,3,2)]
kegg2 <- aggregate(kegg,by=list(kegg$kegg),function(x) paste0(na.omit(unique(x)),collapse =";"))
kegg2 <- kegg2[,-1]



#input
all <- data.frame(data.table::fread("~/Documents/gnhsf/data202509/code/fig2/all_GLM.csv",sep=","))
library(dplyr)
table(all$q.sig)
###
###old  *      **     ***       +    nsig 
###old 5662    2432    2620    6678 1931588 

##new  *      **     ***       +    nsig 
##new  15199    6296    4754   18385 2641731 


glm.cog <- subset(all,level=="cog")
glm.cog <- merge(glm.cog,cog,by.x="feature",by.y="COGnum",all.x=T)
#write.xlsx(glm.cog,"COG_GLM202504.xlsx")

glm.kegg <- subset(all,level=="kegg")
glm.kegg <- merge(glm.kegg,kegg2,by.x="feature",by.y="kegg",all.x=T)
#write.xlsx(glm.kegg,"KEGG_GLM202504.xlsx")

glm.genus <- subset(all,level=="genus")
#write.xlsx(glm.genus,"Genus_GLM202504.xlsx")

glm.species <- subset(all,level=="species")
#write.xlsx(glm.species,"species_GLM202504.xlsx")

glm.human <- subset(all,level=="humanprotein")
glm.human$protein <- data.frame(stringr::str_split_fixed(glm.human$feature,"[.]",3))$X2
glm.human <- merge(glm.human,human2,by.x="protein",by.y="prot",all.x=T)
glm.human$feature <- glm.human$protein
#write.xlsx(glm.human,"human_GLM202504.xlsx")

glm.micro <- subset(all,level=="microprotein")
glm.micro <- merge(glm.micro,annotate,by.x="feature",by.y="glm.match",all.x=T)
#write.xlsx(glm.micro,"micro_GLM202504.xlsx")




#####met shared and specific#######
met.human <- glm.human[glm.human$Cinfo_type=="Metabolic disorders",]
met.human.sig <- met.human[met.human$q.sig !="nsig",]
met.human.number <- data.frame(table(met.human.sig$feature))
met.human.aggre <-aggregate(met.human.sig[,c(1,2,3,14,17:21)],by=list(met.human.sig$feature),function(x) paste0(na.omit(unique(x)),collapse =";"))
met.human.stats <- merge(met.human.number,met.human.aggre,by.x="Var1",by.y="Group.1")
write.xlsx(met.human.stats,"met.human.sig.stats.xlsx")
met.human.list <- unique(met.human.stats$Var1)
met.human.plot <- met.human[met.human$feature %in% met.human.list,]
library(reshape2)
met.human.est <- dcast(met.human.plot,feature+name~metadata,value.var = "Estimate")
met.human.est$name <- ifelse(is.na(met.human.est$name),met.human.est$feature,met.human.est$name)
rownames(met.human.est) <-met.human.est$name
met.human.est.sig <- dcast(met.human.plot,feature+name~metadata,value.var = "q.sig")
met.human.est.sig$name <- ifelse(is.na(met.human.est.sig$name),met.human.est.sig$feature,met.human.est.sig$name)
rownames(met.human.est.sig) <- met.human.est.sig$name
met.human.est.sig[met.human.est.sig=="nsig"] <- ""
library(pheatmap)
breaks <- c(seq(-0.6, 0, length.out = 5), seq(0, 0.4, length.out = 5)) %>% unique()
colors <- c(
  colorRampPalette(c("darkblue", "lightblue"))(4),  # 负数区间4个
  colorRampPalette(c("lightcoral", "darkred"))(4)   # 正数区间3个
  )
pheatmap(met.human.est[,-1:-2],scale="none",show_rownames =T,breaks=breaks,color=colors,display_numbers = as.matrix(met.human.est.sig[,-1:-2]))

met.human.est2 <- met.human.est[,-1:-2]
met.human.est2$annotate <- NA
met.human.est.sig2 <- met.human.est.sig[,-1:-2]
for (i in 1:nrow(met.human.est2)) {
  # 获取当前行数据和p值
  row_values <- as.numeric(met.human.est2[i, 1:(ncol(met.human.est2)-1)])
  row_p <- as.character(met.human.est.sig2[i, ])
  # 条件1: 至少2个显著
  if (sum(row_p != "") >= 2) {
    met.human.est2$annotate[i] <- "share"
  }
  
  # 条件1.1: 全正/全负且至少2个显著
  if ((all(row_values > 0) || all(row_values < 0)) && 
      sum(row_p != "") >= 2) {
    met.human.est2$annotate[i] <- "share.all.same.dir"
  }
  
  # 条件2: 3正1负且负显著
  if (sum(row_values > 0) == 3 && sum(row_values < 0) == 1) {
    neg_col <- which(row_values < 0)
    if (row_p[neg_col] != "") {
      met.human.est2$annotate[i] <- "specific"
    }
  }
  
  # 条件3: 3负1正且正显著
  if (sum(row_values < 0) == 3 && sum(row_values > 0) == 1) {
    pos_col <- which(row_values > 0)
    if (row_p[pos_col] != "") {
      met.human.est2$annotate[i] <- "specific"
    }
  }
}
met.human.est2$feature <- rownames(met.human.est2)
met.human.est3 <- merge(met.human.est2,human2,by.x="feature",by.y="name",all.x=T)

write.xlsx(met.human.est3,"Met_human_est_annotate.xlsx")


#micro
met.micro <- glm.micro[glm.micro$Cinfo_type=="Metabolic disorders",]
met.micro.sig <- met.micro[met.micro$q.sig !="nsig",]
met.micro.number <- data.frame(table(met.micro.sig$feature))
met.micro.aggre <-aggregate(met.micro.sig[,c(1,2,3,13,17:34)],by=list(met.micro.sig$feature),function(x) paste0(na.omit(unique(x)),collapse =";"))
met.micro.stats <- merge(met.micro.number,met.micro.aggre,by.x="Var1",by.y="Group.1")
write.xlsx(met.micro.stats,"met.micro.sig.stats.xlsx")

met.micro.list <- unique(met.micro.stats$Var1)
met.micro.plot <- met.micro[met.micro$feature %in% met.micro.list,]
library(reshape2)
met.micro.est <- dcast(met.micro.plot,feature+Eggpro_prefer_name~metadata,value.var = "Estimate")
met.micro.est$Eggpro_prefer_name <- paste(met.micro.est$feature,met.micro.est$Eggpro_prefer_name,sep="_")
rownames(met.micro.est) <- met.micro.est$Eggpro_prefer_name
met.micro.est.sig <- dcast(met.micro.plot,feature+Eggpro_prefer_name~metadata,value.var = "q.sig")
met.micro.est.sig$Eggpro_prefer_name <- paste(met.micro.est.sig$feature,met.micro.est.sig$Eggpro_prefer_name,sep="_")
rownames(met.micro.est.sig) <- met.micro.est.sig$Eggpro_prefer_name
met.micro.est.sig[met.micro.est.sig=="nsig"] <- ""
library(pheatmap)
breaks <- c(seq(-0.6, 0, length.out = 5), seq(0, 0.4, length.out = 5)) %>% unique()
colors <- c(
  colorRampPalette(c("darkblue", "lightblue"))(4),  # 负数区间4个
  colorRampPalette(c("lightcoral", "darkred"))(4)   # 正数区间3个
)
pheatmap(met.micro.est[,-1:-2],scale="none",show_rownames =T,breaks=breaks,color=colors,display_numbers = as.matrix(met.micro.est.sig[,-1:-2]))

met.micro.est2 <- met.micro.est[,-1:-2]
met.micro.est2$annotate <- NA
met.micro.est.sig2 <- met.micro.est.sig[,-1:-2]
for (i in 1:nrow(met.micro.est2)) {
  # 获取当前行数据和p值
  row_values <- as.numeric(met.micro.est2[i, 1:(ncol(met.micro.est2)-1)])
  row_p <- as.character(met.micro.est.sig2[i, ])
  
  # 条件1: 至少2个显著
  if (sum(row_p != "") >= 2) {
    met.micro.est2$annotate[i] <- "share"
  }
  
  # 条件1.1: 全正/全负且至少2个显著
  if ((all(row_values > 0) || all(row_values < 0)) && 
      sum(row_p != "") >= 2) {
    met.micro.est2$annotate[i] <- "share.all.same.dir"
  }
  
  # 条件2: 3正1负且负显著
  if (sum(row_values > 0) == 3 && sum(row_values < 0) == 1) {
    neg_col <- which(row_values < 0)
    if (row_p[neg_col] != "") {
      met.micro.est2$annotate[i] <- "specific"
    }
  }
  
  # 条件3: 3负1正且正显著
  if (sum(row_values < 0) == 3 && sum(row_values > 0) == 1) {
    pos_col <- which(row_values > 0)
    if (row_p[pos_col] != "") {
      met.micro.est2$annotate[i] <- "specific"
    }
  }
}
met.micro.est2$feature <- met.micro.est$feature
met.micro.est2$prot <- rownames(met.micro.est2)
met.micro.est3 <- merge(met.micro.est2,annotate,by.x="feature",by.y="glm.match")

write.xlsx(met.micro.est3,"Met_micro_est_annotate.xlsx")








##genus,for discovery research, describe "+"####
met.genus <- glm.genus[glm.genus$Cinfo_type=="Metabolic disorders",]
met.genus.sig <- met.genus[met.genus$q.sig !="nsig",]
met.genus.number <- data.frame(table(met.genus.sig$feature))
met.genus.aggre <-aggregate(met.genus.sig[,c(1,2,13)],by=list(met.genus.sig$feature),function(x) paste0(na.omit(unique(x)),collapse =";"))
met.genus.stats <- merge(met.genus.number,met.genus.aggre,by.x="Var1",by.y="Group.1")
write.xlsx(met.genus.stats,"met.genus.sig.stats.xlsx")

met.genus.list <- unique(met.genus.stats$Var1)
met.genus.plot <- met.genus[met.genus$feature %in% met.genus.list,]
met.genus.plot$q.sig <- ifelse(met.genus.plot$q>0.05 & met.genus.plot$q<0.1,"+",met.genus.plot$q.sig)
library(reshape2)
met.genus.est <- dcast(met.genus.plot,feature~metadata,value.var = "Estimate")
rownames(met.genus.est) <- met.genus.est[,1]
met.genus.est.sig <- dcast(met.genus.plot,feature~metadata,value.var = "q.sig")
rownames(met.genus.est.sig) <- met.genus.est.sig[,1]
met.genus.est.sig[met.genus.est.sig=="nsig"] <- ""
library(pheatmap)
breaks <- c(seq(-0.6, 0, length.out = 5), seq(0, 0.4, length.out = 5)) %>% unique()
colors <- c(
  colorRampPalette(c("darkblue", "lightblue"))(4),  # 负数区间4个
  colorRampPalette(c("lightcoral", "darkred"))(4)   # 正数区间3个
)
pheatmap(met.genus.est[,-1],scale="none",show_rownames =T,breaks=breaks,color=colors,display_numbers = as.matrix(met.genus.est.sig[,-1]))

met.genus.est2 <- met.genus.est[,-1]
met.genus.est2$annotate <- NA
met.genus.est.sig2 <- met.genus.est.sig[,-1]
for (i in 1:nrow(met.genus.est2)) {
  # 获取当前行数据和p值
  row_values <- as.numeric(met.genus.est2[i, 1:(ncol(met.genus.est2)-1)])
  row_p <- as.character(met.genus.est.sig2[i, ])
  
  # 条件1: 至少2个显著
  if (sum(row_p != "") >= 2) {
    met.genus.est2$annotate[i] <- "share"
  }
  
  # 条件1.1: 全正/全负且至少2个显著
  if ((all(row_values > 0) || all(row_values < 0)) && 
      sum(row_p != "") >= 2) {
    met.genus.est2$annotate[i] <- "share.all.same.dir"
  }
  
  
  # 条件2: 3正1负且负显著
  if (sum(row_values > 0) == 3 && sum(row_values < 0) == 1) {
    neg_col <- which(row_values < 0)
    if (row_p[neg_col] != "") {
      met.genus.est2$annotate[i] <- "specific"
      next
    }
  }
  
  # 条件3: 3负1正且正显著
  if (sum(row_values < 0) == 3 && sum(row_values > 0) == 1) {
    pos_col <- which(row_values > 0)
    if (row_p[pos_col] != "") {
      met.genus.est2$annotate[i] <- "specific"
    }
  }
}

met.genus.est2$feature <- rownames(met.genus.est2)
write.xlsx(met.genus.est2,"Met_genus_est_annotate.xlsx")






##species,for discovery research, describe "+"####
met.species <- glm.species[glm.species$Cinfo_type=="Metabolic disorders",]
met.species.sig <- met.species[met.species$q.sig !="nsig",]
met.species.number <- data.frame(table(met.species.sig$feature))
met.species.aggre <-aggregate(met.species.sig[,c(1,2,10)],by=list(met.species.sig$feature),function(x) paste0(na.omit(unique(x)),collapse =";"))
met.species.stats <- merge(met.species.number,met.species.aggre,by.x="Var1",by.y="Group.1")
write.xlsx(met.species.stats,"met.species.sig.stats.xlsx")

met.species.list <- unique(met.species.stats$Var1)
met.species.plot <- met.species[met.species$feature %in% met.species.list,]
met.species.plot$q.sig <- ifelse(met.species.plot$q>0.05 & met.species.plot$q<0.1,"+",met.species.plot$q.sig)
library(reshape2)
met.species.est <- dcast(met.species.plot,feature~metadata,value.var = "Estimate")
rownames(met.species.est) <- met.species.est[,1]
met.species.est.sig <- dcast(met.species.plot,feature~metadata,value.var = "q.sig")
rownames(met.species.est.sig) <- met.species.est.sig[,1]
met.species.est.sig[met.species.est.sig=="nsig"] <- ""
library(pheatmap)
breaks <- c(seq(-0.6, 0, length.out = 5), seq(0, 0.4, length.out = 5)) %>% unique()
colors <- c(
  colorRampPalette(c("darkblue", "lightblue"))(4),  # 负数区间4个
  colorRampPalette(c("lightcoral", "darkred"))(4)   # 正数区间3个
)
pheatmap(met.species.est[,-1],scale="none",show_rownames =T,breaks=breaks,color=colors,display_numbers = as.matrix(met.species.est.sig[,-1]))

met.species.est2 <- met.species.est[,-1]
met.species.est2$annotate <- NA
met.species.est.sig2 <- met.species.est.sig[,-1]
for (i in 1:nrow(met.species.est2)) {
  # 获取当前行数据和p值
  row_values <- as.numeric(met.species.est2[i, 1:(ncol(met.species.est2)-1)])
  row_p <- as.character(met.species.est.sig2[i, ])
  
  # 条件1: 至少2个显著
  if (sum(row_p != "") >= 2) {
    met.species.est2$annotate[i] <- "share"
  }
  
  # 条件1.1: 全正/全负且至少2个显著
  if ((all(row_values > 0) || all(row_values < 0)) && 
      sum(row_p != "") >= 2) {
    met.species.est2$annotate[i] <- "share.all.same.dir"
  }
  
  # 条件2: 3正1负且负显著
  if (sum(row_values > 0) == 3 && sum(row_values < 0) == 1) {
    neg_col <- which(row_values < 0)
    if (row_p[neg_col] != "") {
      met.species.est2$annotate[i] <- "specific"
      next
    }
  }
  
  # 条件3: 3负1正且正显著
  if (sum(row_values < 0) == 3 && sum(row_values > 0) == 1) {
    pos_col <- which(row_values > 0)
    if (row_p[pos_col] != "") {
      met.species.est2$annotate[i] <- "specific"
    }
  }
}
met.species.est2$feature <- rownames(met.species.est2)
write.xlsx(met.species.est2,"Met_species_est_annotate.xlsx")



met.cog <- glm.cog[glm.cog$Cinfo_type=="Metabolic disorders",]
met.cog.sig <- met.cog[met.cog$q.sig !="nsig",]
met.cog.number <- data.frame(table(met.cog.sig$feature))
met.cog.aggre <-aggregate(met.cog.sig[,c(1,2,13,16:18)],by=list(met.cog.sig$feature),function(x) paste0(na.omit(unique(x)),collapse =";"))
met.cog.stats <- merge(met.cog.number,met.cog.aggre,by.x="Var1",by.y="Group.1")
write.xlsx(met.cog.stats,"met.cog.sig.stats.xlsx")

met.cog.list <- unique(met.cog.stats$Var1)
met.cog.plot <- met.cog[met.cog$feature %in% met.cog.list,]
library(reshape2)
met.cog.est <- dcast(met.cog.plot,feature~metadata,value.var = "Estimate")
rownames(met.cog.est) <- met.cog.est[,1]
met.cog.est.sig <- dcast(met.cog.plot,feature~metadata,value.var = "q.sig")
rownames(met.cog.est.sig) <- met.cog.est.sig[,1]
met.cog.est.sig[met.cog.est.sig=="nsig"] <- ""
library(pheatmap)
breaks <- c(seq(-0.6, 0, length.out = 5), seq(0, 0.4, length.out = 5)) %>% unique()
colors <- c(
  colorRampPalette(c("darkblue", "lightblue"))(4),  # 负数区间4个
  colorRampPalette(c("lightcoral", "darkred"))(4)   # 正数区间3个
)
pheatmap(met.cog.est[,-1],scale="none",show_rownames =T,breaks=breaks,color=colors,display_numbers = as.matrix(met.cog.est.sig[,-1]))

met.cog.est2 <- met.cog.est[,-1]
met.cog.est2$annotate <- NA
met.cog.est.sig2 <- met.cog.est.sig[,-1]
for (i in 1:nrow(met.cog.est2)) {
  # 获取当前行数据和p值
  row_values <- as.numeric(met.cog.est2[i, 1:(ncol(met.cog.est2)-1)])
  row_p <- as.character(met.cog.est.sig2[i, ])
  
  # 条件1: 至少2个显著
  if (sum(row_p != "") >= 2) {
    met.cog.est2$annotate[i] <- "share"
  }
  
  # 条件1.1: 全正/全负且至少2个显著
  if ((all(row_values > 0) || all(row_values < 0)) && 
      sum(row_p != "") >= 2) {
    met.cog.est2$annotate[i] <- "share.all.same.dir"
  }
  
  # 条件2: 3正1负且负显著
  if (sum(row_values > 0) == 3 && sum(row_values < 0) == 1) {
    neg_col <- which(row_values < 0)
    if (row_p[neg_col] != "") {
      met.cog.est2$annotate[i] <- "specific"
      next
    }
  }
  
  # 条件3: 3负1正且正显著
  if (sum(row_values < 0) == 3 && sum(row_values > 0) == 1) {
    pos_col <- which(row_values > 0)
    if (row_p[pos_col] != "") {
      met.cog.est2$annotate[i] <- "specific"
    }
  }
}
met.cog.est2$feature <- rownames(met.cog.est2)
met.cog.est3 <- merge(met.cog.est2,cog,by.x="feature",by.y="COGnum")

write.xlsx(met.cog.est3,"Met_cog_est_annotate.xlsx")


met.kegg <- glm.kegg[glm.kegg$Cinfo_type=="Metabolic disorders",]
met.kegg.sig <- met.kegg[met.kegg$q.sig !="nsig",]
met.kegg.number <- data.frame(table(met.kegg.sig$feature))
met.kegg.aggre <-aggregate(met.kegg.sig[,c(1,2,3,13,16:19)],by=list(met.kegg.sig$feature),function(x) paste0(na.omit(unique(x)),collapse =";"))
met.kegg.stats <- merge(met.kegg.number,met.kegg.aggre,by.x="Var1",by.y="Group.1")
write.xlsx(met.kegg.stats,"met.kegg.sig.stats.xlsx")

met.kegg.list <- unique(met.kegg.stats$Var1)
met.kegg.plot <- met.kegg[met.kegg$feature %in% met.kegg.list,]
met.kegg.plot$KO <- ifelse(is.na(met.kegg.plot$KO),met.kegg.plot$feature,met.kegg.plot$KO)
library(reshape2)
met.kegg.est <- dcast(met.kegg.plot,KO~metadata,value.var = "Estimate")
rownames(met.kegg.est) <- met.kegg.est[,1]
met.kegg.est.sig <- dcast(met.kegg.plot,KO~metadata,value.var = "q.sig")
rownames(met.kegg.est.sig) <- met.kegg.est.sig[,1]
met.kegg.est.sig[met.kegg.est.sig=="nsig"] <- ""
library(pheatmap)
breaks <- c(seq(-0.6, 0, length.out = 5), seq(0, 0.4, length.out = 5)) %>% unique()
colors <- c(
  colorRampPalette(c("darkblue", "lightblue"))(4),  # 负数区间4个
  colorRampPalette(c("lightcoral", "darkred"))(4)   # 正数区间3个
)
pheatmap(met.kegg.est[,-1],scale="none",show_rownames =T,breaks=breaks,color=colors,display_numbers = met.kegg.est.sig[,-1])

met.kegg.est2 <- met.kegg.est[,-1]
met.kegg.est2$annotate <- NA
met.kegg.est.sig2 <- met.kegg.est.sig[,-1]
for (i in 1:nrow(met.kegg.est2)) {
  # 获取当前行数据和p值
  row_values <- as.numeric(met.kegg.est2[i, 1:(ncol(met.kegg.est2)-1)])
  row_p <- as.character(met.kegg.est.sig2[i, ])
  
  # 条件1: 至少2个显著
  if (sum(row_p != "") >= 2) {
    met.kegg.est2$annotate[i] <- "share"
  }
  
  # 条件1.1: 全正/全负且至少2个显著
  if ((all(row_values > 0) || all(row_values < 0)) && 
      sum(row_p != "") >= 2) {
    met.kegg.est2$annotate[i] <- "share.all.same.dir"
  }
  
  # 条件2: 3正1负且负显著
  if (sum(row_values > 0) == 3 && sum(row_values < 0) == 1) {
    neg_col <- which(row_values < 0)
    if (row_p[neg_col] != "") {
      met.kegg.est2$annotate[i] <- "specific"
      next
    }
  }
  
  # 条件3: 3负1正且正显著
  if (sum(row_values < 0) == 3 && sum(row_values > 0) == 1) {
    pos_col <- which(row_values > 0)
    if (row_p[pos_col] != "") {
      met.kegg.est2$annotate[i] <- "specific"
    }
  }
}
met.kegg.est2$feature <- rownames(met.kegg.est2)
met.kegg.est3 <- merge(met.kegg.est2,kegg2,by.x="feature",by.y="KO",all.x=T)
write.xlsx(met.kegg.est3,"Met_kegg_est_annotate.xlsx")


###met shared heatmap#######
met.human.est2$level <- "human prot"
met.cog.est2$level <- "cog"
met.kegg.est2$level <- "kegg"
met.genus.est2$level <- "genus"
met.species.est2$level <- "species"
met.micro.est2$level <- "microbial prot"

met.micro.est2 <- met.micro.est2[,-6]
colnames(met.micro.est2)[6] <- "feature"

met.est <- rbind(met.human.est2,met.cog.est2,met.kegg.est2,met.genus.est2,met.species.est2,met.micro.est2)
met.est.sig <- rbind(met.human.est.sig2,met.cog.est.sig2,met.kegg.est.sig2,met.genus.est.sig2,met.species.est.sig2,met.micro.est.sig2)

met.est2 <- subset(met.est,level!="microbial prot")
met.est2 <- subset(met.est2,level!="cog")

met.est2 <- subset(met.est2,annotate=="share.all.same.dir")
met.est.sig2 <- met.est.sig[which(rownames(met.est.sig) %in% rownames(met.est2)),]

met.est3 <- merge(met.est2,kegg2,by.x="feature",by.y="kegg",all.x=T)
met.est3$des <- ifelse(is.na(met.est3$des),met.est3$feature,met.est3$des)
met.est.sig3 <- met.est.sig2[match(met.est3$feature,rownames(met.est.sig2)),]
rownames(met.est3) <- met.est3$des
library(dplyr)
anno_row <- met.est3[,7,drop=F]
library(ComplexHeatmap)
breaks <- c(seq(-0.4, 0, length.out = 6), seq(0, 0.4, length.out = 6)) %>% unique()
colors <- c(
  colorRampPalette(c("cyan3", "lightblue"))(5),  # 负数区间4个
  colorRampPalette(c("thistle2","lightcoral"))(5)   # 正数区间3个
)


##fig4a####
aa <- pheatmap(met.est3[,2:5],annotation_row = anno_row,show_rownames = T,row_split=anno_row$level,breaks=breaks,color=colors,display_numbers=as.matrix(met.est.sig3),fontsize_row=6)
row_order <- row_order(aa)

met.est3.ep <- met.est3[unlist(row_order),]
write.csv(met.est3.ep,"fig4a_rebuttal_validate.csv")

####fig4c met shared microbial #####
micro.plot <- subset(met.micro.est3,annotate=="share.all.same.dir")
micro.up.plot <- micro.plot[which(micro.plot$dm_cl>0),]
library(tidyr)
micro.up.plot2 <-micro.up.plot %>% as_tibble() %>%    separate_rows(Taxon_rank_unipept,Taxon_name_unipept, sep = ";")
micro.up.plot2[micro.up.plot2=="NA"] <- NA
library(reshape2)
micro.up.plot3 <- dcast(micro.up.plot2,feature+Protein+Eggpro_des+Eggpro_prefer_name+COG+COG_NCBIdes+COG_NCBIcat+COG_NCBIcat_des+dm_cl+dys_cl+hyper_cl+mets_cl~Taxon_rank_unipept,value.var = "Taxon_name_unipept")
micro.up.plot3$name <- paste0(micro.up.plot3$Protein,"(",micro.up.plot3$Eggpro_prefer_name,")")
write.xlsx(micro.up.plot3,"met_micro.up.plot3.xlsx")
## mannual annotate microbial proteins 


micro.up.plot3 <- read.xlsx("met_micro.up.plot3_2.xlsx")
sangi <- micro.up.plot3[,c(16,8,15,13,14)]
sangi$color <- as.factor(sangi$phylum)
sangi$color <- as.numeric(sangi$color)

c1.1 <- na.omit(sangi[,c(2,1,6)])
colnames(c1.1) <- c("source","target","color")
c1.1.2 <- na.omit(sangi[,c(1,3,6)])
colnames(c1.1.2) <- c("source","target","color")
c1.2 <- na.omit(sangi[,c(3,4,6)])
colnames(c1.2) <- c("source","target","color")
c1.2.2 <- na.omit(sangi[which(is.na(sangi$species)),c(1,4,6)])
colnames(c1.2.2) <- c("source","target","color")
c1.3 <- na.omit(sangi[,c(4,5,6)])
colnames(c1.3) <- c("source","target","color")
c1.3.2 <- na.omit(sangi[which(is.na(sangi$genus)),c(1,5,6)])
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
  "4" = "#E6F598",  "0" = "#F39A9b", "1" = "#a7c6e0", "2" = "#b1dcb0",
  "3" = "#FFc78f", "5" = "#FFFFa5", "6" = "#d8b5a1", "7" = "#Fbc8e3",
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
p <- sankeyNetwork(Links = c1.5, Nodes = id,units = "normalized",nodePadding = 0,nodeWidth = 10,width=1200,linkColor = "colour",
                   Source = 'IDsource', Target = 'IDtarget', Value = 'Freq',dragY = T, dragX = F,showNodeValues = F,scaleNodeBreadthsByString=T,
                  NodeColor = "color",NodeID = 'value', NodeGroup = 'cat', LinkGroup = "colour",NodePosX = "xpos", #xAxisDomain = c("Category","Protein","Species","Genus","Phylum"),
             fontSize = 12,iterations =0)
p

saveNetwork(p,file_name =  "met_micro_up", selfcontained = TRUE,html = "create")

library(webshot2)
webshot2::webshot("met_micro_up.html", "met_micro_up.pdf")
pdf_data <- pdftools::pdf_data("met_micro_up.pdf")[[1]]  # 获取第一页数据
order <- pdf_data$text[pdf_data$x==241]

micro.up.plot3 <- micro.up.plot3[match(order,micro.up.plot3$name),]
rownames(met.micro.est.sig2) <- met.micro.est.sig$feature
micro.up.sig <- met.micro.est.sig2[match(micro.up.plot3$feature,rownames(met.micro.est.sig2)),]
heatmap.plot <- as.matrix(micro.up.plot3[,9:12])
rownames(heatmap.plot) <- micro.up.plot3$name
library(ComplexHeatmap)
pdf("micro_met_up.pdf")
pheatmap(heatmap.plot,show_rownames = T,breaks=breaks,color=colors,display_numbers=as.matrix(micro.up.sig),fontsize_row=6,cluster_rows = F,cluster_cols = F)
dev.off()

micro.down.plot <- micro.plot[which(micro.plot$dm_cl<0),]
micro.down.plot <- micro.down.plot[-which(is.na(micro.down.plot$Taxon_all_unipept)),]
library(tidyr)
micro.down.plot2 <-micro.down.plot %>% as_tibble() %>%    separate_rows(Taxon_rank_unipept,Taxon_name_unipept, sep = ";")
micro.down.plot2[micro.down.plot2=="NA"] <- NA
library(reshape2)
micro.down.plot3 <- dcast(micro.down.plot2,feature+Protein+Eggpro_des+Eggpro_prefer_name+COG+COG_NCBIdes+COG_NCBIcat+COG_NCBIcat_des+dm_cl+dys_cl+hyper_cl+mets_cl~Taxon_rank_unipept,value.var = "Taxon_name_unipept")
micro.down.plot3$name <- paste0(micro.down.plot3$Protein,"(",micro.down.plot3$Eggpro_prefer_name,")")
write.xlsx(micro.down.plot3,"met_micro.down.plot3.xlsx")
#manual impute gene names

library(openxlsx)
micro.down.plot3 <- read.xlsx("met_micro.down.plot3_2.xlsx")
micro.down.plot4 <- subset(micro.down.plot3,genus=="g__Faecalibacterium")
micro.down.plot5 <- subset(micro.down.plot3,genus!="g__Faecalibacterium")


sangi <- micro.down.plot4[,c(16,8,15,13,14)]
#sangi <- micro.down.plot5[,c(16,8,15,13,14)]
sangi$color <- as.factor(ifelse(is.na(sangi$genus),"unknown",sangi$genus))
sangi$color <- as.numeric(sangi$color)

c1.1 <- na.omit(sangi[,c(2,1,6)])
colnames(c1.1) <- c("source","target","color")
c1.1.2 <- na.omit(sangi[,c(1,3,6)])
colnames(c1.1.2) <- c("source","target","color")
c1.2 <- na.omit(sangi[,c(3,4,6)])
colnames(c1.2) <- c("source","target","color")
c1.2.2 <- na.omit(sangi[which(is.na(sangi$species)),c(1,4,6)])
colnames(c1.2.2) <- c("source","target","color")
c1.3 <- na.omit(sangi[,c(4,5,6)])
colnames(c1.3) <- c("source","target","color")
c1.3.2 <- na.omit(sangi[which(is.na(sangi$genus)),c(1,5,6)])
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
  "4" = "#E6F598",  "0" = "#F39A9b", "1" = "#a7c6e0", "2" = "#b1dcb0",
  "3" = "#FFc78f", "5" = "#FFFFa5", "6" = "#d8b5a1", "7" = "#Fbc8e3",
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
id[is.na(id)] <- "unknown"

c1.5$IDsource <- match(c1.5$source, id$value) - 1 
c1.5$IDtarget <- match(c1.5$target, id$value) - 1
c1.5$colour <- as.character(c1.5$color)
c1.5$colour <- color_map[c1.5$colour]


library(sankeyD3plus)
p <- sankeyNetwork(Links = c1.5, Nodes = id,units = "normalized",nodePadding = 1,nodeWidth = 10,width=1200,linkColor = "colour",
                   Source = 'IDsource', Target = 'IDtarget', Value = 'Freq',dragY = T, dragX = F,showNodeValues = F,scaleNodeBreadthsByString=T,
                   NodeColor = "color",NodeID = 'value', NodeGroup = 'cat', LinkGroup = "colour",NodePosX = "xpos", #xAxisDomain = c("Category","Protein","Species","Genus","Phylum"),
                   fontSize = 12,orderByPath=F,curvature=0.5,iterations = 5000)
p

saveNetwork(p,file_name =  "met_micro_down_1", selfcontained = TRUE,html = "create")

library(webshot2)　
webshot2::webshot("met_micro_down.html", "met_micro_down.pdf")
pdf_data <- pdftools::pdf_data("met_micro_down_1.pdf")[[1]]  # 获取第一页数据
order <- pdf_data$text[pdf_data$x==206]
#position will change


micro.down.plot3 <- micro.down.plot3[match(order,micro.down.plot3$name),]
rownames(met.micro.est.sig2) <- met.micro.est.sig$feature
micro.down.sig <- met.micro.est.sig2[match(micro.down.plot3$feature,rownames(met.micro.est.sig2)),]
heatmap.plot <- as.matrix(micro.down.plot3[,9:12])
rownames(heatmap.plot) <- micro.down.plot3$name
library(ComplexHeatmap)
breaks <- c(seq(-0.4, 0, length.out = 6), seq(0, 0.4, length.out = 6)) %>% unique()
colors <- c(
  colorRampPalette(c("cyan3", "lightblue"))(5),  # 负数区间4个
  colorRampPalette(c("thistle2","lightcoral"))(5)   # 正数区间3个
)
pdf("micro_met_down1.pdf")
pheatmap(heatmap.plot,show_rownames = T,breaks=breaks,color=colors,display_numbers=as.matrix(micro.down.sig),fontsize_row=6,cluster_rows = F,cluster_cols = F)
dev.off()





###met.same.dir. kegg. pathway.enrich#####
met.est.enrich <- subset(met.est,annotate=="share.all.same.dir")
kegg.list <- met.est.enrich$feature[met.est.enrich$level=="kegg"]
kegg.list <- data.frame(str_split_fixed(kegg.list," ",2))$X1


library(MicrobiomeProfiler)
#run_MicrobiomeProfiler()

core.prokegg.enrich <- enrichKO(kegg.list,qvalueCutoff = 0.2,minGSSize = 10,maxGSSize = 500,pvalueCutoff = 0.2)
prokegg <- core.prokegg.enrich@result
library(stringr)
prokegg$bg <- as.numeric(data.frame(str_split_fixed(prokegg$BgRatio,"/",2))$X1) / as.numeric(data.frame(str_split_fixed(prokegg$BgRatio,"/",2))$X2)
prokegg$gene <- as.numeric(data.frame(str_split_fixed(prokegg$GeneRatio,"/",2))$X1) / as.numeric(data.frame(str_split_fixed(prokegg$GeneRatio,"/",2))$X2)
prokegg$foldenrich <- prokegg$gene/prokegg$bg



##fig4b#####
enrichplot::dotplot(core.prokegg.enrich,showCategory=6,color="pvalue",size="FoldEnrichment")



