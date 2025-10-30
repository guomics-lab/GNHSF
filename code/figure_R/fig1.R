
setwd("~/Documents/gnhsf/data202509/code/fig1")
rm(list=ls())

library(data.table)
file <- list.files("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample")

data.test <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_filter2_genus_sample.txt",sep="\t")
samplename <- colnames(data.test)[-1]
ident.sample <- data.frame("sample"=samplename)
ident.sample.total <- c()
for (x in file) {
  data <- fread(paste0("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/",x),sep="\t")
  sums <- rowSums(data[,2:ncol(data)],na.rm=T)
  if(length(which(sums==0))>0){
    data <- data[-sums,]
  }
  data <- as.data.frame(data)
  data <- data[,match(samplename,colnames(data))]
  ident.data <- apply(data,2,function(x) length(which(!is.na(x))))
  ident.sample$x <- ident.data
  colnames(ident.sample)[ncol(ident.sample)] <- x
  ident.sample.total <- c(ident.sample.total,nrow(data))
  cat("\r",x)
}



meta <- read.csv("~/Documents/gnhsf/data202509/GNHSF_sample_ident_20220307.csv")

meta <- merge(meta[,2:9],ident.sample,by.y="sample",by.x="label3",all.x=T)

write.csv(meta,"~/Documents/gnhsf/data202509/GNHSF_sample_ident_202508.csv")


mean <- colMeans(meta[,9:ncol(meta)],na.rm=T)
sd <-  apply(meta[9:ncol(meta)],2,function(x) sd(x,na.rm=T))
options(scipen = 999)
mean
sd

names(ident.sample.total) <- names(mean)

colnames(ident.sample) <- gsub("GNHSF_diann_IGC_humanswiss_","",colnames(ident.sample))
colnames(ident.sample) <- gsub("_sample.txt|_sample.tsv","",colnames(ident.sample))

list <- c("humanprotein$","microprotein$","microprotein_kegg$","microprotein_cog$","genus","species")
ident.sample1 <- ident.sample[,grep(paste(list,collapse = "|"),colnames(ident.sample))] 
rownames(ident.sample1) <- ident.sample$sample
ident.sample1 <- ident.sample1[c(3,4,5,6,1,2)]

##generate pr matrix####
pr <- data.frame(fread("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/diann/GNHSF_metaExpertPro_lib_2514_DIA_diann1.8.1.pr_matrix.tsv",sep="\t"))
pr1 <- pr[,c(10:ncol(pr))]
library(stringr)
replace<-data.frame(str_split_fixed(colnames(pr1),"GNHSF_",3))
colnames(pr1) <- replace$X2

replace<-data.frame(str_split_fixed(meta$name,"GNHSF_",3))
replace$replace <- ifelse(replace$X3=="",replace$X2,replace$X3)
meta$replace.name <- replace$replace
sample.list <- subset(meta,label2=="sample")

pr2 <- pr1[,c(1,which(colnames(pr1) %in% sample.list$replace.name))]
rows <- rowSums(pr2[,2:ncol(pr2)],na.rm=T)
pr2 <- pr2[-which(rows==0),]
#total ident
nrow(pr2)
pr.ident <- c()
for (x in 2:ncol(pr2)) {
  pr.ident[x] <- length(which(!is.na(pr2[,x])))
}
median(pr.ident,na.rm = T)
sd(pr.ident,na.rm=T)

pr2.ident <- data.frame(pr.ident)
pr2.ident$sample <- colnames(pr2)

metadata <- merge(meta,pr2.ident,by.x="replace.name",by.y="sample",all.x=T)
write.csv(meta,"~/Documents/gnhsf/data202509/GNHSF_sample_ident_202508.csv")
write.csv(meta,"~/Documents/gnhsf/data202509/metadata/GNHSF_sample_ident_202508.csv")

##fig 1 d####
sample.metadata <- subset(metadata,label2=="sample")
ident <- sample.metadata[,c(2,91,78,84,66,70,80,82)]
ident <- reshape::melt(ident)
ident$variable <- factor(ident$variable,levels=c("pr.ident","GNHSF_diann_IGC_humanswiss_humanprotein_sample.tsv","GNHSF_diann_IGC_humanswiss_microprotein_sample.tsv",
                                                 "GNHSF_diann_IGC_humanswiss_filter5_genus_sample.txt","GNHSF_diann_IGC_humanswiss_filter5_species_sample.txt",
                                                 "GNHSF_diann_IGC_humanswiss_microprotein_kegg_sample.txt",
                                                 "GNHSF_diann_IGC_humanswiss_microprotein_cog_sample.txt"))
col <- c("#A8c5e6","#85bf53","#d696a4","#348083","#877ea7","#f89b31","#7895a3")

library(ggplot2)
ggplot(ident,aes(x=1,y=value,fill=variable))+geom_boxplot()+facet_wrap(.~variable,scales="free_y",nrow=1,axis = "all_y")+ggthemes::theme_base()+
  scale_fill_manual(values=col)
ggsave("fig1dident.pdf",width=15)


##fig s2c####
col <- c("#85bf53","#d696a4","#f89b31","#7895a3","#348083","#877ea7")
cv <- c()
pdf("sample_ient_hist.pdf",width=10,height=5)
par(mfrow=c(2,4))
for( i in 1:ncol(ident.sample1)){
  cv[i] <- sd(ident.sample1[,i])/mean(ident.sample1[,i])
  hist(ident.sample1[,i],breaks=30,freq=T,include.lowest=T,right=T,xlab=paste0(colnames(ident.sample1)[i]," richness"),main=paste0("CV = ",round(cv[i],2)),col=col[i])
}
dev.off()

##fig s2d####
human <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_humanprotein_sample.tsv",sep="\t")
micro <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_microprotein_sample.tsv",sep="\t"

##quant
human.sum <- as.data.frame(colSums(human[,2:ncol(human)],na.rm = T))
colnames(human.sum) <- "abundance"
micro.sum <- as.data.frame(colSums(micro[,2:ncol(micro)],na.rm = T))
colnames(micro.sum) <- "abundance"

human.sum <- human.sum[match(rownames(micro.sum),rownames(human.sum)),,drop=F]

ratio <- micro.sum$abundance/human.sum$abundance

human.sum$log <- log2(human.sum$abundance)
micro.sum$log <- log2(micro.sum$abundance)

median <- median(ratio,na.rm=T)
sd <-  sd(ratio,na.rm=T)

library(ggplot2)

human.sum$ratio <- ratio
p1=ggplot(human.sum,aes(x=ratio))+geom_histogram(fill="#eecfa1",color="black",bins=50)+xlab("Ratio of Human to Microbial Protein Abundance per Individual")+theme_bw()+ylab("Count")+
  geom_vline(xintercept = median,linetype="dashed",color="red",size=1)+geom_text(aes(x=median+2,y=150,label="Median = 5.899 , SD = 4.826"))
ggsave("ratio.pdf",width=6,height=5)


###ident
micro <- micro[,match(colnames(human),colnames(micro))]
human.ident <- apply(human[,2:ncol(human)],2,function(x) length(which(!is.na(x))))
micro.ident <- apply(micro[,2:ncol(micro)],2,function(x) length(which(!is.na(x))))
ident.ratio <- micro.ident/human.ident
human.sum$ident.ratio <- ident.ratio
median <- median(ident.ratio,na.rm=T)
sd <- sd(ident.ratio,na.rm=T)
p2=ggplot(human.sum,aes(x=ident.ratio))+geom_histogram(fill="#bce6ff",color="black",bins=50)+xlab("Ratio of Human to Microbial Protein identification per Individual")+theme_bw()+ylab("Count")+
  geom_vline(xintercept = median,linetype="dashed",color="red",size=1)+geom_text(aes(x=median+10,y=100,label="Median = 29.984 , SD = 12.219"))
ggsave("identratio.pdf",width=6,height=5)

library(patchwork)
p1+p2
ggsave("ratio_identratio.pdf",width=10,height=5)
