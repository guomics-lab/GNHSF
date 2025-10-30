rm(list=ls())

##arrange glm results####
cate <- c("Anthropometrics" ,"Biochemical measures" ,"Diet" ,"Lifestyle"  , "Medicine" ,"Metabolic disorders" , "Self-report diseases" ,"Sociodemographics"   )
level <- c("cog","kegg","genus","species","humanprotein","microprotein")

file <- expand.grid(level,cate)

all.glm <- data.frame()
for (x in 1:nrow(file)) {
level.x <-  as.character(file$Var1[x] )
cate.x <- as.character(file$Var2[x])
test <- read.delim(paste0("~/Documents/gnhsf/data202509/code/glm/glm/",level.x,"_NA90_",cate.x,"_glm_res_adj_category.csv"),sep=",")
test$level <- level.x
all.glm <- rbind(all.glm,test)
cat("\r",x)
}
write.csv(all.glm,"all_GLM.csv",row.names = F)

##figs7b####
all.glm <- data.frame(data.table::fread("~/Documents/gnhsf/data202509/code/fig2/all_GLM.csv",sep=","))
all.sig <- subset(all.glm,q.sig!="nsig")
write.csv(all.sig,"ALL_GLM_SIG.csv",row.names = F)
stats1 <- data.frame(table(all.sig$metadata)) 
stats <-data.frame(table(all.sig[,14:15])) 
stats$level <- factor(stats$level,levels=(c("humanprotein",'microprotein',"cog","kegg","genus","species")))
library(ggplot2)
library(ggthemes)

color <- c("Anthropometrics"="#BE3B29","Biochemical measures"="#0472b7","Sociodemographics"="#e28727","Diet"="#0f864e","Lifestyle"="#7876b3","Medicine"="#6f99ad","Metabolic disorders"="#fedc91"
,"Self report diseases"="#ee4c98")

ggplot(stats,aes(x=Cinfo_type,y=Freq,fill=Cinfo_type))+geom_bar(stat="identity")+coord_flip()+scale_fill_manual(values=color)+theme_base()+facet_grid(.~level,scales = "free_x")
ggsave("figs7b_glm_number_by_categroy.pdf",width=12,height = 5)



#all.self <- subset(all.glm,Cinfo_type=="Self-report diseases")
#table(all.self$q.sig)
#nsig

table(all.sig$q.sig)
#   *    **   *** 
#  15199  6296  4754 



##fig2a table####
stats1 <- data.frame(table(all.sig[,c(1)]))
#dmmed,bbs,tg,dm,tea,hypermed,age
stats2 <-data.frame(table(all.sig[, c(1,15)])) 
stats3 <- data.frame(table(all.sig$Cinfo_type))

##fig2b heatmap####
for (x in level) {
text <- subset(all.glm,level==x)
write.csv(text,paste0("~/Documents/gnhsf/data202509/code/glm/glm/",x, "_NA90_all_meta_glm_res_adj_categroy.csv"),row.names = F)   
}

library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(reshape2)
library(stringr)
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
pie(rep(1, 25), col = c25)

mycolors <- c25

##annotate####
library(openxlsx)
anno.kegg <- read.xlsx("~/Documents/gnhsf/data202509/metadata/kegg_database20230328.xlsx")
anno.cog <- read.xlsx("~/Documents/gnhsf/data202509/metadata/01.NCBI.COG.list.xlsx")
anno.cog <- anno.cog[, c(1, 2)]
anno.human <- read.csv("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_human_all.txt",sep="\t")
anno.human <- unique(anno.human[,-1])
anno.human$Protein <- as.character(lapply(strsplit(as.character(anno.human$Protein), split="\\|"), tail, n=1))
anno.human$Protein <- gsub("_HUMAN", "", anno.human$Protein)
anno.human <- anno.human[, c(1,6,7)]


metatype.vect <- c("Sociodemographics", "Metabolic disorders", "Lifestyle", "Biochemical measures", "Anthropometrics","Medicine")
meta <- c("bristol_scale","dm_cl","TG","tea_f","dm_med" ,"age")
top <- 10

#####genus####


level="genus"
df0 <- read.csv(paste0("~/Documents/gnhsf/data202509/code/glm/glm/",level, "_NA90_all_meta_glm_res_adj_categroy.csv"))
df <- subset(df0,metadata %in% meta)
df$metadata <- paste0(df$metadata, "_", df$metadata_value)
dfsig <- df[df$q.sig != "nsig", ]
dfsig <- dfsig %>% group_by(metadata)
dfsig.arr <- dfsig %>% arrange(desc(abs(Estimate)), .by_group = T)
metau <- unique(dfsig.arr$metadata)
feature <- c()
for (i in 1:length(metau)){
  dfnow <- dfsig.arr[dfsig.arr$metadata == metau[i], ]
  if (nrow(dfnow) < top){
    top.new <- nrow(dfnow)
  }else{
    top.new <- top
  }
  feature.add <- dfnow[c(1:top.new), ]$feature
  feature <- append(feature, feature.add)
}
coef.all <- data.frame()
qsig.all <- data.frame()
for (i in 1:length(metau)){
  coef <- df[((df$feature %in% feature) & (df$metadata == metau[i])), c("feature", "Estimate")]
  qsig <- df[((df$feature %in% feature) & (df$metadata == metau[i])), c("feature", "q.sig")]
  qsig$q.sig <- gsub("nsig", "", qsig$q.sig)
  colnames(coef)[2] <- metau[i]
  colnames(qsig)[2] <- metau[i]
  if (i == 1){
    coef.all <- coef
    qsig.all <- qsig
  }else{
    coef.all <- merge(coef.all, coef, by.x = "feature", by.y = "feature")
    qsig.all <- merge(qsig.all, qsig, by.x = "feature", by.y = "feature")
  }
}
rownames(coef.all) <- coef.all$feature
rownames(qsig.all) <- qsig.all$feature
rownames(coef.all) <- as.character(lapply(strsplit(as.character(rownames(coef.all)), split="\\."), tail, n=1))
rownames(qsig.all) <- as.character(lapply(strsplit(as.character(rownames(qsig.all)), split="\\."), tail, n=1))
coef.all.1 <- as.data.frame(coef.all[, -1])
qsig.all.1 <- as.data.frame(qsig.all[, -1])
colnames(coef.all.1) <- colnames(coef.all)[-1]
rownames(coef.all.1) <- rownames(coef.all)
colnames(qsig.all.1) <- colnames(qsig.all)[-1]
rownames(qsig.all.1) <- rownames(qsig.all)
coef.all <- as.matrix(coef.all.1)
qsig.all <- as.matrix(qsig.all.1)
tax <- as.data.frame(str_split(feature, "\\.", simplify = T))
tax <- tax %>% distinct()
anno.new <- as.data.frame(tax[, 3])
rownames(anno.new) <- tax[, ncol(tax)]
colnames(anno.new) <- "feature"

anno.new.u <- unique(anno.new$feature)
mycolors.now <- mycolors[1:length(anno.new.u)]
anno.colors.list <- list(setNames(mycolors.now, anno.new.u))
names(anno.colors.list) <- c("feature")


pwidth = 10
pdf(file = paste0("GNHSF_1385_", level, "_NA90_covariate_glm_res_adj_ex4_heatmap.pdf"), width = pwidth, height = 10)
print(pheatmap(coef.all, cluster_rows = T, cluster_cols = F, display_numbers = qsig.all,
               annotation_row = anno.new, fontsize_row = 10, fontsize_col = 10,
               annotation_colors = anno.colors.list,
               fontsize_number = 8, fontsize = 10, border_color = "#E5E5E5",
               angle_col = "45", na_col = "white", number_color = "black",
               cellwidth = 12, cellheight = 10))
dev.off()





####kegg####
level <- "kegg"
df0 <- read.csv(paste0("~/Documents/gnhsf/data202509/code/glm/glm/",level, "_NA90_all_meta_glm_res_adj_categroy.csv"))
df <- subset(df0,metadata %in% meta)
dfname <- anno.kegg[, c(6,5)]
dfname$des <- as.character(lapply(strsplit(as.character(dfname$des), split="\\;"), tail, n=1))
dfname$des <- as.character(lapply(strsplit(as.character(dfname$des), split="\\[EC"), head, n=1))
dfmerge <- unique(merge(df, dfname, by.x = "feature", by.y = "kegg",all.x=T))
df$metadata <- paste0(df$metadata, "_", df$metadata_value)
dfsig <- df[df$q.sig != "nsig", ]
dfsig <- dfsig %>% group_by(metadata)
dfsig.arr <- dfsig %>% arrange(desc(abs(Estimate)), .by_group = T)
metau <- unique(dfsig.arr$metadata)
feature <- c()
for (i in 1:length(metau)){
  dfnow <- dfsig.arr[dfsig.arr$metadata == metau[i], ]
  if (nrow(dfnow) < top){
    top.new <- nrow(dfnow)
  }else{
    top.new <- top
  }
  feature.add <- dfnow[c(1:top.new), ]$feature
  feature <- append(feature, feature.add)
}
coef.all <- data.frame()
qsig.all <- data.frame()
for (i in 1:length(metau)){
  coef <- df[((df$feature %in% feature) & (df$metadata == metau[i])), c("feature", "Estimate")]
  qsig <- df[((df$feature %in% feature) & (df$metadata == metau[i])), c("feature", "q.sig")]
  qsig$q.sig <- gsub("nsig", "", qsig$q.sig)
  colnames(coef)[2] <- metau[i]
  colnames(qsig)[2] <- metau[i]
  if (i == 1){
    coef.all <- coef
    qsig.all <- qsig
  }else{
    coef.all <- merge(coef.all, coef, by.x = "feature", by.y = "feature")
    qsig.all <- merge(qsig.all, qsig, by.x = "feature", by.y = "feature")
  }
}
rownames(coef.all) <- coef.all$feature
rownames(qsig.all) <- coef.all$feature
coef.all <- coef.all[,-1]
qsig.all <- qsig.all[,-1]

anno<- anno.kegg[anno.kegg$kegg %in% feature,]

anno <- anno[,c(6,5,2)] %>% distinct()
anno.new <- aggregate(anno,by=list(anno$kegg),FUN=function(x) paste(unique(x),collapse=";"))
anno.new$des <- as.character(lapply(strsplit(as.character(anno.new$des), split="\\;"), tail, n=1))
anno.new$des <- as.character(lapply(strsplit(as.character(anno.new$des), split="\\[EC"), head, n=1))
anno.new <- anno.new[match(rownames(coef.all),anno.new$Group.1),]    
anno.new$anno <-  paste(anno.new$kegg,anno.new$des,sep =":")
rownames(anno.new) <- anno.new$anno
anno.new1 <- as.data.frame(anno.new[, 4,drop=F])
colnames(anno.new1) <-"feature"

anno.new.u1 <- unique(anno.new1$feature)
anno.new.u <- c(anno.new.u1[-grep(";",anno.new.u1)],anno.new.u1[grep(";",anno.new.u1)])
mycolors.now <- mycolors[1:length(anno.new.u)]
anno.colors.list <- list(setNames(mycolors.now, anno.new.u))
names(anno.colors.list) <- c("feature")


rownames(coef.all) <- anno.new$anno
write.table(anno.new,"heatmap_ko_anno.txt",sep="\t")

pdf(file = paste0("GNHSF_1385_", level, "_NA90_covariate_glm_res_adj_ex4_heatmap.pdf"), width = 20, height = 10)
print(pheatmap(coef.all, cluster_rows = T, cluster_cols = F, display_numbers = as.matrix(qsig.all),
               annotation_row = anno.new1, fontsize_row = 10, fontsize_col = 10,
               annotation_colors = anno.colors.list,
               fontsize_number = 8, fontsize = 10, border_color = "#E5E5E5",
               angle_col = "45", na_col = "white", number_color = "black",
               cellwidth = 12, cellheight = 10))
dev.off()




###human protein####
level="humanprotein"
df0 <- read.csv(paste0("~/Documents/gnhsf/data202509/code/glm/glm/",level, "_NA90_all_meta_glm_res_adj_categroy.csv"))
df <- subset(df0,metadata %in% meta)
df$feature <- as.character(lapply(strsplit(as.character(df$feature), split="\\."), tail, n=1))
df$feature <- gsub("_HUMAN", "", df$feature)
df$metadata <- paste0(df$metadata, "_", df$metadata_value)
dfsig <- df[df$q.sig != "nsig", ]
dfsig <- dfsig %>% group_by(metadata)
dfsig.arr <- dfsig %>% arrange(desc(abs(Estimate)), .by_group = T)
metau <- unique(dfsig.arr$metadata)
feature <- c()
for (i in 1:length(metau)){
  dfnow <- dfsig.arr[dfsig.arr$metadata == metau[i], ]
  if (nrow(dfnow) < top){
    top.new <- nrow(dfnow)
  }else{
    top.new <- top
  }
  feature.add <- dfnow[c(1:top.new), ]$feature
  feature <- append(feature, feature.add)
}
coef.all <- data.frame()
qsig.all <- data.frame()
for (i in 1:length(metau)){
  coef <- df[((df$feature %in% feature) & (df$metadata == metau[i])), c("feature", "Estimate")]
  qsig <- df[((df$feature %in% feature) & (df$metadata == metau[i])), c("feature", "q.sig")]
  qsig$q.sig <- gsub("nsig", "", qsig$q.sig)
  colnames(coef)[2] <- metau[i]
  colnames(qsig)[2] <- metau[i]
  if (i == 1){
    coef.all <- coef
    qsig.all <- qsig
  }else{
    coef.all <- merge(coef.all, coef, by.x = "feature", by.y = "feature")
    qsig.all <- merge(qsig.all, qsig, by.x = "feature", by.y = "feature")
  }
}
rownames(coef.all) <- coef.all$feature
rownames(qsig.all) <- qsig.all$feature
coef.all.1 <- as.data.frame(coef.all[, -1])
qsig.all.1 <- as.data.frame(qsig.all[, -1])
colnames(coef.all.1) <- colnames(coef.all)[-1]
rownames(coef.all.1) <- rownames(coef.all)
colnames(qsig.all.1) <- colnames(qsig.all)[-1]
rownames(qsig.all.1) <- rownames(qsig.all)
coef.all <- as.matrix(coef.all.1)
qsig.all <- as.matrix(qsig.all.1)



anno <- anno.human[which(anno.human$Protein %in% feature),]
write.csv(anno,"human_heatmap_to_manual_anno_cat.csv",row.names = F)

anno <- read.xlsx("human_heatmap_to_manual_anno_cat_2.xlsx")
anno <- anno[match(rownames(coef.all),anno$Protein),]
rownames(anno) <- anno$Protein
anno.new <- data.frame(   anno[,c(-1,-2)] )
anno.new[anno.new==""] <- "unknown"
colnames(anno.new) <- "feature"
rownames(anno.new) <- rownames(anno)
anno.new.u <- unique(anno.new$feature)
anno.new.u[anno.new.u==""] <- "unknown"
mycolors.now <- mycolors[1:length(anno.new.u)]
anno.colors.list <-list(setNames(mycolors.now, anno.new.u))


rownames(coef.all) <- anno$Protein

pwidth = 10
pdf(file = paste0("GNHSF_1385_", level, "_NA90_covariate_glm_res_adj_ex4_heatmap.pdf"), width = pwidth, height = 10)
print(pheatmap(coef.all, cluster_rows = T, cluster_cols = F, display_numbers = as.matrix(qsig.all),
               annotation_row = anno.new, fontsize_row = 10, fontsize_col = 10,
               annotation_colors = anno.colors.list,
               fontsize_number = 8, fontsize = 10, border_color = "#E5E5E5",
               angle_col = "45", na_col = "white", number_color = "black",
               cellwidth = 12, cellheight = 10))
dev.off()

