rm(list=ls())


setwd("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_954")
#abundance NA90，impute0 and int###
## NA cutoff 0.9 and impute 0

func.prev <- function(x){
  sum(is.na(x))
}

genus <- read.delim("GNHSF_diann_IGC_humanswiss_filter5_genus_sample_954.txt",sep="\t")
species <- read.delim("GNHSF_diann_IGC_humanswiss_filter5_species_sample_954.txt",sep="\t")
cog <- read.delim("GNHSF_diann_IGC_humanswiss_microprotein_cog_sample_954.txt",sep="\t")
kegg <- read.delim("GNHSF_diann_IGC_humanswiss_microprotein_kegg_sample_954.txt",sep="\t")
humanprot <- read.delim("GNHSF_diann_IGC_humanswiss_humanprotein_sample_954.txt",sep="\t")
microprot <- read.delim("GNHSF_diann_IGC_humanswiss_microprotein_sample_954.txt",sep="\t")



func.NA <- function(df){
  type=deparse(substitute(df))
  df$prevalence <- ((ncol(df)-1) - as.numeric(apply(df[, 2:ncol(df)], 1, func.prev))) / (ncol(df)-1)
  df_NA90 <- df[df$prevalence > 0.1,]
  df_NA90 <- df_NA90[,-ncol(df_NA90)]
  df_NA90[is.na(df_NA90)] <- 0
  write.table(df_NA90, file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_954_NA90.tsv"), sep = "\t", row.names = F)
}
func.NA(genus)
func.NA(species)
func.NA(cog)
func.NA(kegg)
func.NA(humanprot)
func.NA(microprot)

## INT transform
func.INT <- function(x){
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}
func.INT.abun <- function(type){
  df <- as.data.frame(read.delim(file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_954_NA90.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE))
  df.int <- as.data.frame(t(apply(df[,2:ncol(df)], 1, func.INT)), check.name = F)
  df.int <- cbind(df[,1], df.int)
  colnames(df.int)[1] <- colnames(df)[1]
  write.table(df.int, file = paste0("GNHSF_diann_IGC_humanswiss_", type, "_sample_954_NA90_INT.tsv"), sep = "\t", row.names = F)
}
func.INT.abun("cog")
func.INT.abun("kegg")
func.INT.abun("species")
func.INT.abun("genus")
func.INT.abun("microprot")
func.INT.abun("humanprot")

##readmetadata####
metadata1 <- read.csv("~/Documents/gnhsf/data202509/metadata/GNHSF_sample_inform_normalized_954_int.csv")
metadata1$Bristol.scale <- metadata1$Bristol.scale-1
metadata1$Marriage <- metadata1$Marriage-1
metadata2 <- read.csv("~/Documents/gnhsf/data202509/metadata/17foodgroups_impute_1967.csv",sep=",")
colnames(metadata2)[1] <- "sample"
metadata3 <- read.csv("~/Documents/gnhsf/data202509/metadata/22disease_1967_all.csv",sep=",")
metadata3 <- metadata3[,c(1,6:25)]
metadata3 <- metadata3[,c(1,which(colnames(metadata3) %in% c("insomnia","gout","gastric_duodenal_ulcer","urolithiasis","cataract")))]
metadata3[metadata3=="yes"]=1
metadata3[metadata3=="no"]=0
metadata3[,2:ncol(metadata3)] <- lapply(metadata3[,2:ncol(metadata3),drop=F],as.factor)

metadata <- Reduce(function(x,y) merge(x,y,by="sample"),list(metadata1,metadata2,metadata3),accumulate=F)

metadata$id <- substr(metadata$sample,3,nchar(metadata$sample))
metadata$time <- metadata$Age
metadata$stage <- substr(metadata$sample,1,2)
metadata <- arrange(metadata,sample)

clino <- colnames(metadata)[c(-1,-55:-ncol(metadata))]
clino.del.covari <- clino[c(-1,-2,-11)]
factor.list <- c( "Sex" ,"Income","Edu","Marriage"  , "Tea"  ,"Smoke"  ,"Alc" ,"Bristol.scale","Hyper.TGmia" , "HBP" ,  "HBP.med"  , "Hyper.LDLemia","Hypo.HDLemia" , "Hypercholesterolemia" ,"MetS"   , "DM"  ,"DM.med"  , 
                  "DLP" ,  "DLP.med","id","stage","gout","gastric_duodenal_ulcer","urolithiasis","cataract","insomnia"    )
metadata[,which(colnames(metadata) %in% factor.list)] <-  lapply(metadata[,which(colnames(metadata) %in% factor.list), drop = FALSE], as.factor)




###glmm######
library(lme4)
library(lmerTest)
library(MuMIn)



get.glmm <- function(m){
  feature1 <- read.delim(paste0("GNHSF_diann_IGC_humanswiss_",m,"_sample_954_NA90_INT.tsv"),sep="\t")
  rownames(feature1) <- feature1[,1]
  feature1 <- feature1[,-1]
  feature2 <- data.frame(t(feature1))
  library(dplyr)
  feature2 <- feature2[match(metadata$sample,rownames(feature2)),]
  
  feature.data1 <- data.frame()
  for (x in clino.del.covari) {
    feature.data <- data.frame()
    for(i in 1:ncol(feature2)){
      metadata$feature <- feature2[,i]
      formula <- formula(paste0("feature~ Age + Sex + Bristol.scale +",x,"+ (1|id)"))
      
      glmm_model <- lmer(formula,data = metadata)
      glm.sum <- summary(glmm_model)
      feature.data0 <- data.frame(glm.sum$coefficients[-1:-5,,drop=F])
      feature.data0$feature <- colnames(feature2)[i]
      feature.data0$metadata <- rownames(feature.data0)
      feature.data <- rbind(feature.data,feature.data0)
      cat("\r",m,x,i)
    }
    feature.data1 <- rbind(feature.data1,feature.data)
  }
  
  
  for(i in 1:ncol(feature2)){
    metadata$feature <- feature2[,i]
    formula <- formula(paste0("feature~ Age + Sex + Bristol.scale + (1|id)"))
    
    glmm_model <- lmer(formula,data = metadata)
    glm.sum <- summary(glmm_model)
    feature.data0 <- data.frame(glm.sum$coefficients[-1,])
    feature.data0$feature <- colnames(feature2)[i]
    feature.data0$metadata <- rownames(feature.data0)
    feature.data1 <- rbind(feature.data1,feature.data0)
    cat("\r",i)
  }
  write.table(feature.data1,paste0("GNHSF_954_GLMM_",m,"_adjust_id.txt"),sep="\t",row.names = F)
  print("done")
}

get.glmm("genus")
get.glmm("species")
get.glmm("cog")
get.glmm("kegg")
get.glmm("humanprot")
get.glmm("microprot")




###use parallel#####
library(lme4)
library(lmerTest)
library(MuMIn)
library(parallel)
library(dplyr)
library(pbmcapply)  

get.glmm <- function(m){
  feature1 <- read.delim(paste0("GNHSF_diann_IGC_humanswiss_",m,"_sample_954_NA90_INT.tsv"),sep="\t")
  rownames(feature1) <- feature1[,1]
  feature1 <- feature1[,-1]
  feature2 <- data.frame(t(feature1))
  feature2 <- feature2[match(metadata$sample,rownames(feature2)),]
  
  # 并行处理函数
  process_feature <- function(i) {
    metadata$feature <- feature2[,i]
    results <- list()
    
    # 处理带协变量的模型
    for (x in clino.del.covari) {
      formula <- formula(paste0("feature~ Age + Sex + Bristol.scale +",x,"+ (1|id)"))
      glmm_model <- lmer(formula, data = metadata)
      glm.sum <- summary(glmm_model)
      feature.data0 <- data.frame(glm.sum$coefficients[-1:-5,,drop=F])
      feature.data0$feature <- colnames(feature2)[i]
      feature.data0$metadata <- rownames(feature.data0)
      results[[length(results)+1]] <- feature.data0
    }
    
    # 处理基础模型
    formula <- formula(paste0("feature~ Age + Sex + Bristol.scale + (1|id)"))
    glmm_model <- lmer(formula, data = metadata)
    glm.sum <- summary(glmm_model)
    feature.data0 <- data.frame(glm.sum$coefficients[-1,])
    feature.data0$feature <- colnames(feature2)[i]
    feature.data0$metadata <- rownames(feature.data0)
    results[[length(results)+1]] <- feature.data0
    
    return(do.call(rbind, results))
  }
  
  # 检测CPU核心数
  num_cores <- detectCores() - 1
  
  # 使用pbmclapply显示进度条
  cat("开始并行处理", ncol(feature2), "个特征，使用", num_cores, "个核心...\n")
  feature_results <- pbmclapply(1:ncol(feature2), process_feature, mc.cores = num_cores)
  
  # 合并结果
  feature.data1 <- do.call(rbind, feature_results)
  
  write.table(feature.data1, paste0("~/Documents/gnhsf/data202509/code/fig3_age/GNHSF_954_GLMM_",m,"_adjust_id.txt"), sep="\t", row.names = F)
  cat("完成处理:", m, "\n")
  return("done")
}

# 并行运行不同的数据类型
data_types <- c("genus", "species", "cog", "kegg", "humanprot", "microprot")
data_types <- c("microprot")

for(dtype in data_types) {
  cat("\n=== 处理数据类型:", dtype, "===\n")
  get.glmm(dtype)
}



##adjust p value####
library(stringr)
micro <- read.delim("~/Documents/gnhsf/data202509/code/fig3_age/GNHSF_954_GLMM_microprot_adjust_id.txt",sep="\t")
feature.name <- data.frame(str_split_fixed(micro$feature,"[.][.]",2))
feature.name$X3 <- data.frame(str_split_fixed(feature.name$X1,"[.]strand",2))$X1
micro$feature <- feature.name$X3
micro$feature <- gsub("^X","",micro$feature)
micro$feature <- paste0("X",micro$feature)
write.table(micro,"~/Documents/gnhsf/data202509/code/fig3_age/GNHSF_954_GLMM_microprot_adjust_id.txt",sep="\t",row.names = F)

human <- read.delim("~/Documents/gnhsf/data202509/code/fig3_age/GNHSF_954_GLMM_humanprot_adjust_id.txt",sep="\t")
feature.name <- data.frame(str_split_fixed(human$feature,"HUMAN[.]",2))
feature.name$X3 <- paste0(feature.name$X1,"HUMAN")
human$feature <- feature.name$X3
write.table(human,"~/Documents/gnhsf/data202509/code/fig3_age/GNHSF_954_GLMM_humanprot_adjust_id.txt",sep="\t",row.names = F)



metadata <- read.csv("~/Documents/gnhsf/data202509/metadata/metadata_classify (1).csv")
metadata <- metadata[,2:3]
adjp <- function(m){
  df <- read.delim(paste0("~/Documents/gnhsf/data202509/code/fig3_age/GNHSF_954_GLMM_",m,"_adjust_id.txt"),sep="\t")
  df$meta <- gsub("[0-9]$","",df$metadata)
  meta.list <- data.frame(table(df$meta))
  df <- merge(df,metadata,by="meta")
  df1 <- do.call(rbind,lapply(split(df,df$Cinfo_type),function(x){
    x$adj.p <- p.adjust(x$Pr...t..,method="BH")
    return(x)} ))
  df1$adj.p.annotate <- NA
  df1$adj.p.annotate <- ifelse(df1$adj.p<0.1,"*",df1$adj.p.annotate)
  df1$adj.p.annotate <- ifelse(df1$adj.p<0.05,"*",df1$adj.p.annotate)
  df1$adj.p.annotate <- ifelse(df1$adj.p<0.01,"**",df1$adj.p.annotate)
  df1$adj.p.annotate <- ifelse(df1$adj.p<0.001,"***",df1$adj.p.annotate)
  
  df1$level <- m
  write.table(df1,paste0("~/Documents/gnhsf/data202509/code/fig3_age/GNHSF_954_GLMM_",m,"_adjust_id_adjp.txt"),sep="\t",row.names=F)
}
adjp("genus")
adjp("species")
adjp("cog")
adjp("kegg")
adjp("humanprot")
adjp("microprot")

levels=c("cog","kegg","genus","species","humanprot","microprot")
all <- do.call(rbind, lapply(levels,function(x){
  df <- read.delim(paste0("~/Documents/gnhsf/data202509/code/fig3_age/GNHSF_954_GLMM_",x,"_adjust_id_adjp.txt"),sep="\t")
  return(df)}))
colnames(all)[10:11] <- c("q","q.sig")
write.table(all,"~/Documents/gnhsf/data202509/code/fig3_age/GNHSF_954_GLMM_all_adjust_id_adjp.txt",sep="\t")
########stats########
cate <- c("Anthropometrics" ,"Biochemical measures" ,"Diet" ,"Lifestyle"  , "Medicine" ,"Metabolic disorders" , "Self-report diseases" ,"Sociodemographics"   )

catee <- function(x){
  df1 <- subset(all,Cinfo_type==x)
  write.table(df1,paste0(x,"_glmm202504.txt"),sep="\t")
}

for (i in cate) {
  catee(i)
  cat(i)
}

all.sig <- subset(all,q.sig !="nsig")

#######age#######
age.glmm <- subset(all,metadata=="Age")
age.glmm$q.sig <- ifelse(age.glmm$q>0.05&age.glmm$q<0.1,"+",age.glmm$q.sig)
write.csv(age.glmm,"~/Documents/gnhsf/data202509/code/fig3_age/GLMM_age20250417.csv")
colnames(age.glmm) <- paste0(colnames(age.glmm),"_glmm954")
table(age.glmm$q.sig,age.glmm$level)


glm <- data.frame(data.table::fread("~/Documents/gnhsf/data202509/code/fig2/all_GLM.csv",sep=","))
age.glm <- subset(glm,metadata=="age")
age.glm2 <- age.glm[,c(2,5,6,12:15)]
colnames(age.glm2) <- paste0(colnames(age.glm2),"_glm1385")

table(age.glm2$q.sig_glm1385,age.glm2$level_glm1385)

age.all <- merge(age.glmm,age.glm2,by.x="feature_glmm954",by.y="feature_glm1385",all=T)
write.csv(age.all,"~/Documents/gnhsf/data202509/code/fig3_age/age.glm.glmm_20250417.csv")



rm(list=ls())
#####age.glm.glmm#######
age <- read.delim("~/Documents/gnhsf/data202509/code/fig3_age/age.glm.glmm_20250417.csv",sep=",")

table(age$q.sig_glmm954)

#glmm. p 0.05,  species 0, genus 1, cog50, kegg 72, human 14,  micro 270
#glmm  p 0.01,  species 0, genus 1,  cog 21, kegg 36, human 10, micro 93

####fig3a. venn#####
age1 <- age[which(age$q.sig_glm1385 !="nsig"),]
age.venn <- age1[,c(2,19,8,11,14,16)]
age.venn1 <- age.venn

age.level <- split(age.venn1,age.venn1$level_glm1385)
library(Vennerable)  
lapply(age.level,function(x){
  between.subject <- x$feature_glmm954[which(x$q_glm1385<0.05)]
  within.subject <- x$feature_glmm954[which(x$Pr...t.._glmm954<0.01)]
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
  pdf(paste0("~/Documents/gnhsf/data202509/code/fig3_age/",unique(x$level),".venn.pdf"))
  plot(data,gp=gp)  
  dev.off()  
}  )




####fig3b. age genus 1385####

genus <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_filter5_genus_sample_1385.txt",sep="\t")
meta <- read.delim("~/Documents/gnhsf/data202509/metadata/GNHSF_sample_inform_normalized_1385.csv",sep=",")

genus$Taxon <- gsub("[|]",".",genus$Taxon)
genus1 <- genus[which(genus$Taxon %in% age1$feature),]
colnames(genus1)[1] <- "feature"

feature <- genus1
feature <- data.frame(t(feature))
colnames(feature) <- feature[1,]
feature <- feature[-1,]
feature$sample <- rownames(feature)
meta1 <- meta[,c(1,3)]

feature1 <- merge(meta1,feature,by="sample")

feature1[,2:ncol(feature1)] <- apply(feature1[,2:ncol(feature1)],2,function(x) as.numeric(x))
feature1[is.na(feature1)] <- 0
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(ggpointdensity)

ggplot(feature1,aes(x=age,y=d__Bacteria.k__Bacillati.p__Bacillota.c__Negativicutes.o__Acidaminococcales.f__Acidaminococcaceae.g__Phascolarctobacterium))+geom_pointdensity()+
  scale_color_continuous(low = "#d9bdd8", high = "#83247f")+
  geom_smooth(method="lm",color="#3c2260",fill="#d9bdd8",size=1.2)+stat_cor()+theme_classic()
#stat_poly_eq(aes(label=paste(..eq.label..,..adj.rr.label..,..p.value.label..,sep = "~~~~")),formula = y~x,parse=T,size=5)
ggsave("~/Documents/gnhsf/data202509/code/fig3_age/Phascolarctobacterium_age.pdf")



setwd("~/Documents/gnhsf/data202509/code/fig3_age")

##fig3c####
kegg.list2 <- na.omit(age1$feature_glmm954[age1$level_glmm954=="kegg" & age1$Pr...t.._glmm954 <0.01])
library(MicrobiomeProfiler)
#run_MicrobiomeProfiler()
library(openxlsx)
anno <- read.xlsx("~/Documents/gnhsf/data202509/metadata/kegg_database20230328.xlsx")
ko.anno <- subset(anno,kegg %in% kegg.list2)
ko.anno <- unique(merge(age1,ko.anno[,4:6],by.x="feature_glmm954",by.y="kegg"))
plot.ko <- unique(ko.anno[,c(20,4)])
library(dplyr)
plot.ko <- arrange(plot.ko,Estimate_glmm954)
plot.ko <- plot.ko[c(1:5,29:36),]
pdf("~/Documents/gnhsf/data202509/code/fig3_age/kegg_estimate.pdf")
barplot(plot.ko$Estimate_glmm954,names.arg = plot.ko$KO,las=2,cex.names=0.2,horiz = TRUE)
dev.off()

##fig3d####
core.prokegg.enrich <- enrichKO(kegg.list2,qvalueCutoff = 0.05,minGSSize = 10,maxGSSize = 500,pvalueCutoff = 0.1)
prokegg <- core.prokegg.enrich@result
pdf("~/Documents/gnhsf/data202509/code/fig3_age/kegg_pathway.pdf")
enrichplot::dotplot(core.prokegg.enrich,showCategory=20)
dev.off()

####age micro####
###fig 3e####
annotate1 <- read.xlsx("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microprot_2phylum2genus2species_all20250916.xlsx")
annotate1$Protein  <- gsub("^X","",annotate1$Protein)
annotate1$Protein <- paste0("X",annotate1$Protein)
annotate1$Protein <-  gsub("[-]",".",annotate1$Protein)

age.micro <- subset(age1,level_glm1385=="microprotein" & Pr...t.._glmm954 <0.01)
age.micro1 <- merge(age.micro,annotate1,by.x="feature_glmm954",by.y="Protein")
write.csv(age.micro1,"~/Documents/gnhsf/data202509/code/fig3_age/age_glm_glmm_micro.csv",row.names = F)
## manual annotate COG cat ####

age.micro1 <- read.xlsx("~/Documents/gnhsf/data202509/code/fig3_age/age_glm_glmm_micro_1.xlsx")
library(tidyr)
age.micro2 <-age.micro1 %>% as_tibble() %>%    separate_rows(Taxon_rank_unipept,Taxon_name_unipept, sep = ";")
age.micro2 <- age.micro2[which(age.micro2$Taxon_rank_unipept %in% c("phylum","genus","species")),]
age.micro2[age.micro2=="NA"] <- NA
age.micro2 <- age.micro2[-which(is.na(age.micro2$Taxon_name_unipept)),]
library(reshape2)
age.micro3 <- dcast(age.micro2,feature_glmm954+Eggpro_des+Eggpro_prefer_name+COG_NCBIcat+COG_NCBIcat_des+Estimate_glm1385+Estimate_glmm954~Taxon_rank_unipept,value.var=c("Taxon_name_unipept"))
age.micro3$name <- age.micro3$Eggpro_prefer_name
age.micro3$value <- abs(age.micro3$Estimate_glm1385)

library(ggraph)
library(tidygraph)
library(ccgraph)
library(ggthemes)
nodes <- gather_graph_node(age.micro3,index=colnames(age.micro3)[c(9,8,10,11)],value="value")
edge <- gather_graph_edge(age.micro3,index=colnames(age.micro3)[c(9,8,10,11)])
graph<- tbl_graph(nodes,edge)
ggraph(graph, 'circlepack',weight=node.size) + 
  geom_node_circle(data = . %>% filter(node.short_name != "NA"),aes(fill = factor(node.level)),alpha=0.5) + 
  coord_fixed() +
  theme_graph()+ 
  geom_node_text(data = . %>% filter(node.short_name != "NA"),
    aes(label=node.short_name,color=leaf),repel=T,
    fontface="bold",
    size=3
  )+ scale_fill_manual(values=c("#e9d8e1","#91d4f0","#f1f4ca","#fcdc89"))+scale_color_manual(values = c("black","darkblue")) +theme_map()
ggsave("~/Documents/gnhsf/data202509/code/fig3_age/identification_category_count.pdf",width=10,height=10)



###fig 3f####
age.micro4 <- split(age.micro3,age.micro3$phylum)


get_sangi <- function(sangi.network){
library(sankeyD3plus)
sangi.network$COG_NCBIcat_des <- ifelse(is.na(sangi.network$COG_NCBIcat_des),"Unknown",sangi.network$COG_NCBIcat_des)
sangi <- sangi.network[,c(5,11,10,8,9,12)]
sangi$color <- sangi$genus
sangi$color <- as.factor(ifelse(is.na(sangi$color),"Unkonw",sangi$color) )
sangi$color <- as.numeric(sangi$color)
#sangi$color[is.na(sangi$color)] <- "17"
#sangi <- sangi[-which(is.na(sangi$genus)),]
c1.1 <- na.omit(sangi[,c(1,2,7)])
colnames(c1.1) <- c("source","target","color")
c1.1.2 <- na.omit(sangi[,c(2,3,7)])
colnames(c1.1.2) <- c("source","target","color")
c1.2 <- na.omit(sangi[,c(3,4,7)])
colnames(c1.2) <- c("source","target","color")
c1.2.2 <- sangi[which(is.na(sangi$species)& !is.na(sangi$genus)),c(2,4,7)]
colnames(c1.2.2) <- c("source","target","color")
c1.3 <- na.omit(sangi[,c(4,5,7)])
colnames(c1.3) <- c("source","target","color")
c1.3.2 <-sangi[which(is.na(sangi$genus)& !is.na(sangi$phylum)),c(2,5,7)]
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
p <- sankeyNetwork(Links = c1.5, Nodes = id,units = "normalized",nodePadding = 1,nodeWidth = 10,width=1200,height=500, LinkGroup = "colour",linkColor = "colour",
                   Source = 'IDsource', Target = 'IDtarget', Value = 'Freq',dragY = T, dragX = F,showNodeValues = F,scaleNodeBreadthsByString=T,
                   NodeColor = "color",NodeID = 'value',NodePosX = "xpos", #xAxisDomain = c("Category","Protein","Species","Genus","Phylum"),
                   fontSize = 12,iterations=4000)

p
file_name <- paste0("age_micro_",unique(sangi$phylum))

saveNetwork(p,file_name = file_name, selfcontained = TRUE,html = "create")
}


lapply(age.micro4,get_sangi)

dialister <- subset(age.micro3,genus %in% c("g__Dialister","g__Veillonella","g__Coprococcus"))
get_sangi(dialister)
nondialister <- subset(age.micro4[[2]],!(genus %in% c("g__Dialister","g__Veillonella","g__Coprococcus")))
get_sangi(nondialister)

