rm(list=ls())



setwd("~/Documents/gnhsf/data202509/code/figs1")
####BC matrix#####

######biorep#####
bio.bc.file <- list.files("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/biorep_BC",full.names = T)

standardize_ident <- function(x) {
  if(length(grep("^b",x)>0)){
    parts <- strsplit(x, "_")
    batch <- sapply(parts, `[`, 1)
    position <- sapply(parts, `[`, 2)
    batch <- sub("^b", "", batch)  
    batch <- sprintf("b%02d", as.integer(batch))  
    pos_num <- sub("^([0-9]+)([a-zA-Z]*)$", "\\1", position)  
    pos_text <- sub("^([0-9]+)([a-zA-Z]*)$", "\\2", position) 
    pos_num <- sprintf("%02d", as.integer(pos_num)) 
    position <- paste0(pos_num, pos_text)
    rename <- paste0(batch, "_", position)
  }else{
    rename <- x
  }
}
biobcmatrix <- data.frame(matrix(nrow=169))
for (x in bio.bc.file) {
  bc.test <- read.table(x,sep="\t")
  colnames(bc.test) <- sapply(colnames(bc.test) ,standardize_ident)
  rownames(bc.test) <- sapply(rownames(bc.test) ,standardize_ident)
  bio1 <- bc.test
  bio2 <- gsub("_22","_02",colnames(bio1))
  bio2 <- gsub("_23","_21",bio2)
  bio2 <- gsub("biorep","",bio2)
  
  bioname <- data.frame(colnames(bio1),bio2)
  bioname$cat <- ifelse(bioname$colnames.bio1.==bioname$bio2,0,1)
  bionames <- subset(bioname,bioname$cat==1)
  
  biobcx <- c()
  for (i in 1:nrow(bionames)) {
    a <- bionames$colnames.bio1.[i]
    b <- bionames$bio2[i]
    biobcx[i] <- as.numeric(subset(bio1,rownames(bio1)==a,colnames(bio1)==b))
  }
  biobcmatrix$biobc <- biobcx
  level <- strsplit(x,"humanswiss_",10)[[1]]
  level1 <- level[length(level)]
  level1 <- gsub("[.]txt|[.]tsv","",level1)
  colnames(biobcmatrix)[ncol(biobcmatrix)] <- level1
  cat("\r", x)
}
biobcmatrix <- biobcmatrix[,-1]

####techrep########
tech.bc.file <- list.files("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/techrep_BC",full.names = T)

standardize_ident <- function(x) {
  if(length(grep("^b",x)>0)){
    parts <- strsplit(x, "_")
    batch <- sapply(parts, `[`, 1)
    position <- sapply(parts, `[`, 2)
    batch <- sub("^b", "", batch)  
    batch <- sprintf("b%02d", as.integer(batch))  
    pos_num <- sub("^([0-9]+)([a-zA-Z]*)$", "\\1", position) 
    pos_text <- sub("^([0-9]+)([a-zA-Z]*)$", "\\2", position) 
    pos_num <- sprintf("%02d", as.integer(pos_num))  
    position <- paste0(pos_num, pos_text)
    rename <- paste0(batch, "_", position)
  }else{
    rename <- x
  }
}
techbcmatrix <- data.frame(matrix(nrow=185))
for (x in tech.bc.file) {
  bc.test <- read.table(x,sep="\t")
  colnames(bc.test) <- sapply(colnames(bc.test) ,standardize_ident)
  rownames(bc.test) <- sapply(rownames(bc.test) ,standardize_ident)
  tech1 <- bc.test
  tech2 <- gsub("_25","_10",colnames(tech1))
  tech2 <- gsub("_26","_20",tech2)
  tech2 <- gsub("techrep","",tech2)
  
  techname <- data.frame(colnames(tech1),tech2)
  techname$cat <- ifelse(techname$colnames.tech1.==techname$tech2,0,1)
  technames <- subset(techname,techname$cat==1)
  
  techbcx <- c()
  for (i in 1:nrow(technames)) {
    a <- technames$colnames.tech1.[i]
    b <- technames$tech2[i]
    techbcx[i] <- as.numeric(subset(tech1,rownames(tech1)==a,colnames(tech1)==b))
  }
  techbcmatrix$techbc <- techbcx
  level <- strsplit(x,"humanswiss_",10)[[1]]
  level1 <- level[length(level)]
  level1 <- gsub("[.]txt|[.]tsv","",level1)
  colnames(techbcmatrix)[ncol(techbcmatrix)] <- level1
  cat("\r", x)
}
techbcmatrix <- techbcmatrix[,-1]





###pool#############

info <- read.csv("~/Documents/gnhsf/data202509/metadata/GNHSF_pool_QC_20210922.csv")
colnames(info)[3] <- "qc"
pool.bc.file <- list.files("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/pool_BC",full.names = T)

poolbcmatrix <- data.frame(matrix(nrow=3648))
for (x in pool.bc.file) {

pooltest<- read.delim(x,sep="\t",header=T)
colnames(pooltest) <- sapply(colnames(pooltest) ,standardize_ident)
rownames(pooltest) <- sapply(rownames(pooltest) ,standardize_ident)

colnames(pooltest) <- gsub("_27","",colnames(pooltest))
rownames(pooltest) <-  gsub("_27","",rownames(pooltest))
colnames(pooltest) <- gsub("pool_2120","b04add",colnames(pooltest))
colnames(pooltest) <- gsub("pool_2146","b05add",colnames(pooltest))
rownames(pooltest) <- gsub("pool_2120","b04add",rownames(pooltest))
rownames(pooltest) <- gsub("pool_2146","b05add",rownames(pooltest))

pooltest$X <- rownames(pooltest)

testlabel <- merge(info,pooltest,by.x="batch",by.y="X")

pooltest1 <- subset(testlabel,testlabel$pool=="pool1")
pooltest1 <- subset(pooltest1,select = pooltest1$batch)
pooltest2 <- subset(testlabel,testlabel$pool=="pool2")
pooltest2 <- subset(pooltest2,select = pooltest2$batch)

pooltestbc <- c(pooltest1[upper.tri(pooltest1,diag=F)],pooltest2[upper.tri(pooltest2,diag=F)])
poolbcmatrix$x <- pooltestbc
level <- strsplit(x,"humanswiss_",10)[[1]]
level1 <- level[length(level)]
level1 <- gsub("[.]txt|[.]tsv","",level1)
colnames(poolbcmatrix)[ncol(poolbcmatrix)] <- level1
cat("\r", x)

}



######qc###############

##qc humanprot###
qc.bc.file <- list.files("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/qc_BC",full.names = T)


qcbcmatrix <- data.frame(matrix(nrow=1274))

for (x in qc.bc.file) {
  
qctest<- read.delim(x,sep="\t",header=T)
colnames(qctest) <- sapply(colnames(qctest) ,standardize_ident)
rownames(qctest) <- sapply(rownames(qctest) ,standardize_ident)

colnames(qctest) <- gsub("_24","",colnames(qctest))
colnames(qctest) <- gsub("_10","",colnames(qctest))

qctest$X <- rownames(qctest)
  qctest$X <- gsub("_24","",qctest$X)
  qctest$X <- gsub("_10","",qctest$X)
  
rownames(qctest) <- qctest$X


testlabel <- merge(info,qctest,by.x="batch",by.y="X")

qctest1 <- subset(testlabel,testlabel$qc=="qc1")
qctest1 <- subset(qctest1,select = qctest1$batch)
qctest2 <- subset(testlabel,testlabel$qc=="qc2")
qctest2 <- subset(qctest2,select = qctest2$batch)
qctest3 <- subset(testlabel,testlabel$qc=="qc3")
qctest3 <- subset(qctest3,select = qctest3$batch)
qctest4 <- subset(testlabel,testlabel$qc=="qc4")
qctest4 <- subset(qctest4,select = qctest4$batch)
qctest5 <- subset(testlabel,testlabel$qc=="qc5")
qctest5 <- subset(qctest5,select = qctest5$batch)

qctestbc <- c(qctest1[upper.tri(qctest1,diag=F)],qctest2[upper.tri(qctest2,diag=F)],qctest3[upper.tri(qctest3,diag=F)],qctest4[upper.tri(qctest4,diag=F)],qctest5[upper.tri(qctest5,diag=F)])
qcbcmatrix$x <- qctestbc
level <- strsplit(x,"humanswiss_",10)[[1]]
level1 <- level[length(level)]
level1 <- gsub("[.]txt|[.]tsv","",level1)
colnames(qcbcmatrix)[ncol(qcbcmatrix)] <- level1
cat("\r", x)
}





##correlation#####
######biorep#####
bio.corr.file <- list.files("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/biorep",full.names = T)

standardize_ident <- function(x) {
  if(length(grep("^b",x)>0)){
    parts <- strsplit(x, "_")
    batch <- sapply(parts, `[`, 1)
    position <- sapply(parts, `[`, 2)
    batch <- sub("^b", "", batch)  
    batch <- sprintf("b%02d", as.integer(batch))  
    pos_num <- sub("^([0-9]+)([a-zA-Z]*)$", "\\1", position) 
    pos_text <- sub("^([0-9]+)([a-zA-Z]*)$", "\\2", position)  
    pos_num <- sprintf("%02d", as.integer(pos_num)) 
    position <- paste0(pos_num, pos_text)
    rename <- paste0(batch, "_", position)
  }else{
    rename <- x
  }
}
biocorrmatrix <- data.frame(matrix(nrow=169))
for (x in bio.corr.file) {
  corr.test <- read.table(x,sep="\t",header = T)
  colnames(corr.test) <- sapply(colnames(corr.test) ,standardize_ident)
  bio1 <- corr.test[,-1]
  colnames(bio1) <- gsub("_22","_02",colnames(bio1))
  colnames(bio1) <- gsub("_23","_21",colnames(bio1))
  colnames(bio1) <- gsub("biorep","",colnames(bio1))
  
  bioname <- unique(colnames(bio1))

  biocorrx <- c()
  for (i in 1:length(bioname)) {
    aa <- bio1[,grep(bioname[i],colnames(bio1))]
    aa <- as.data.frame(apply(aa, 2, as.numeric))
    biocorrx[i] <- cor(aa[,1],aa[,2],use ="pairwise.complete.obs",method = "spearman")
      }
  biocorrmatrix$biocorr <- biocorrx
  level <- strsplit(x,"humanswiss_",10)[[1]]
  level1 <- level[length(level)]
  level1 <- gsub("[.]txt|[.]tsv","",level1)
  colnames(biocorrmatrix)[ncol(biocorrmatrix)] <- level1
  cat("\r", x)
}
biocorrmatrix <- biocorrmatrix[,-1]



####techrep########
tech.corr.file <- list.files("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/techrep",full.names = T)

standardize_ident <- function(x) {
  if(length(grep("^b",x)>0)){
    parts <- strsplit(x, "_")
    batch <- sapply(parts, `[`, 1)
    position <- sapply(parts, `[`, 2)
    batch <- sub("^b", "", batch)  
    batch <- sprintf("b%02d", as.integer(batch))  
    pos_num <- sub("^([0-9]+)([a-zA-Z]*)$", "\\1", position)  
    pos_text <- sub("^([0-9]+)([a-zA-Z]*)$", "\\2", position)  
    pos_num <- sprintf("%02d", as.integer(pos_num)) 
    position <- paste0(pos_num, pos_text)
    rename <- paste0(batch, "_", position)
  }else{
    rename <- x
  }
}
techcorrmatrix <- data.frame(matrix(nrow=185))
for (x in tech.corr.file) {
  corr.test <- read.table(x,sep="\t",header = T)
  colnames(corr.test) <- sapply(colnames(corr.test) ,standardize_ident)
  tech1 <- corr.test[,-1]
  colnames(tech1) <- gsub("_25","_10",colnames(tech1))
  colnames(tech1) <- gsub("_26","_20",colnames(tech1))
  colnames(tech1) <- gsub("techrep","",colnames(tech1))
  
  techname <- unique(colnames(tech1))
  
  techcorrx <- c()
  for (i in 1:length(techname)) {
    aa <- tech1[,grep(techname[i],colnames(tech1))]
    aa <- as.data.frame(apply(aa, 2, as.numeric))
    techcorrx[i] <- cor(aa[,1],aa[,2],use ="pairwise.complete.obs",method = "spearman")
  }
  techcorrmatrix$techcorr <- techcorrx
  level <- strsplit(x,"humanswiss_",10)[[1]]
  level1 <- level[length(level)]
  level1 <- gsub("[.]txt|[.]tsv","",level1)
  colnames(techcorrmatrix)[ncol(techcorrmatrix)] <- level1
  cat("\r", x)
}
techcorrmatrix <- techcorrmatrix[,-1]





###pool#############

info <- read.csv("~/Documents/gnhsf/data202509/metadata/GNHSF_pool_QC_20210922.csv")
colnames(info)[3] <- "qc"
pool.corr.file <- list.files("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/pool",full.names = T)

poolcorrmatrix <- data.frame(matrix(nrow=3648))
for (x in pool.corr.file) {
  
  pooltest<- read.delim(x,sep="\t",header=T)
  colnames(pooltest) <- sapply(colnames(pooltest) ,standardize_ident)

  colnames(pooltest) <- gsub("_27","",colnames(pooltest))
  colnames(pooltest) <- gsub("pool_2120","b04add",colnames(pooltest))
  colnames(pooltest) <- gsub("pool_2146","b05add",colnames(pooltest))

  pooltest <- data.frame(t(pooltest))
  

  colnames(pooltest) <- pooltest[1,]
  pooltest <- pooltest[-1,]
  pooltest$X <- rownames(pooltest)
  
  
  testlabel <- merge(info,pooltest,by.x="batch",by.y="X")
  
  pooltest1 <- subset(testlabel,testlabel$pool=="pool1")
  pooltest1 <- data.frame(t(pooltest1))
  pooltest1 <- apply(pooltest1[-1:-6,],2,as.numeric)
  pooltest1cor <- cor(pooltest1,use ="pairwise.complete.obs",method = "spearman")
  pooltest2 <- subset(testlabel,testlabel$pool=="pool2")
  pooltest2 <- data.frame(t(pooltest2))
  pooltest2 <- apply(pooltest2[-1:-6,],2,as.numeric)
  pooltest2cor <- cor(pooltest2,use ="pairwise.complete.obs",method = "spearman")
  pooltestcor <- c(pooltest1cor[upper.tri(pooltest1cor,diag=F)],pooltest2cor[upper.tri(pooltest2cor,diag=F)])
  
  poolcorrmatrix$x <- pooltestcor
  level <- strsplit(x,"humanswiss_",10)[[1]]
  level1 <- level[length(level)]
  level1 <- gsub("[.]txt|[.]tsv","",level1)
  colnames(poolcorrmatrix)[ncol(poolcorrmatrix)] <- level1
  cat("\r", x)
  
}



######qc###############

qc.corr.file <- list.files("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/qc",full.names = T)


qccorrmatrix <- data.frame(matrix(nrow=1274))

for (x in qc.corr.file) {
  
  qctest<- read.delim(x,sep="\t",header=T)
  colnames(qctest) <- sapply(colnames(qctest) ,standardize_ident)
  colnames(qctest) <- gsub("_24","",colnames(qctest))
  colnames(qctest) <- gsub("_10","",colnames(qctest))
  
  qctest <- data.frame(t(qctest))
  
  colnames(qctest) <- qctest[1,]
  qctest <- qctest[-1,]
  qctest$X <- rownames(qctest)
  
  
  testlabel <- merge(info,qctest,by.x="batch",by.y="X")
  
  qctest1 <- subset(testlabel,testlabel$qc=="qc1")
  qctest1 <- data.frame(t(qctest1))
  qctest1 <- apply(qctest1[-1:-6,],2,as.numeric)
  qctest1cor <- cor(qctest1,use ="pairwise.complete.obs",method = "spearman")
  qctest2 <- subset(testlabel,testlabel$qc=="qc2")
  qctest2 <- data.frame(t(qctest2))
  qctest2 <- apply(qctest2[-1:-6,],2,as.numeric)
  qctest2cor <- cor(qctest2,use ="pairwise.complete.obs",method = "spearman")
  
  qctest3 <- subset(testlabel,testlabel$qc=="qc3")
  qctest3 <- data.frame(t(qctest3))
  qctest3 <- apply(qctest3[-1:-6,],2,as.numeric)
  qctest3cor <- cor(qctest3,use ="pairwise.complete.obs",method = "spearman")
  qctest4 <- subset(testlabel,testlabel$qc=="qc4")
  qctest4 <- data.frame(t(qctest4))
  qctest4 <- apply(qctest4[-1:-6,],2,as.numeric)
  qctest4cor <- cor(qctest4,use ="pairwise.complete.obs",method = "spearman")
  qctest4 <- subset(testlabel,testlabel$qc=="qc4")
  
  qctest5 <- subset(testlabel,testlabel$qc=="qc5")
  qctest5 <- data.frame(t(qctest5))
  qctest5 <- apply(qctest5[-1:-6,],2,as.numeric)
  qctest5cor <- cor(qctest5,use ="pairwise.complete.obs",method = "spearman")
  
  qctestcor <- c(qctest1cor[upper.tri(qctest1cor,diag=F)],qctest2cor[upper.tri(qctest2cor,diag=F)],qctest3cor[upper.tri(qctest3cor,diag=F)],qctest4cor[upper.tri(qctest4cor,diag=F)],qctest5cor[upper.tri(qctest5cor,diag=F)])

    qccorrmatrix$x <- qctestcor
  level <- strsplit(x,"humanswiss_",10)[[1]]
  level1 <- level[length(level)]
  level1 <- gsub("[.]txt|[.]tsv","",level1)
  colnames(qccorrmatrix)[ncol(qccorrmatrix)] <- level1
  cat("\r", x)
}




###fig s1 b####
pdf("Rep_BC.pdf",width=10,height=4)
par(mfrow=c(1,4))
BCbio <- data.frame("human_protein"=biobcmatrix$humanprotein_biorep,"microbial_protein"=biobcmatrix$microprotein_biorep,"genus"=biobcmatrix$filter5_genus_biorep,"species"=biobcmatrix$filter5_species_biorep,"cog"=biobcmatrix$microprotein_cog_biorep,"kegg"=biobcmatrix$microprotein_kegg_biorep)
library(vioplot)
vioplot(BCbio$human_protein,BCbio$microbial_protein,BCbio$genus,BCbio$species,BCbio$cog,BCbio$kegg,
        col=c("#8bc63e","#d696a4","#348083","#877ea7","#f89b31","#7895a3"),names = c("Human Protein","Microbial Protein","Genus","Species","COG","KEGG"),main = "Intra-batch biological replicates",ylim=c(0,1))
text(labels=c(round(median(biobcmatrix$humanprotein_biorep),2),round(median(biobcmatrix$microprotein_biorep),2),round(median(biobcmatrix$filter5_genus_biorep),2),round(median(biobcmatrix$filter5_species_biorep),2),round(median(biobcmatrix$microprotein_cog_biorep),2),round(median(biobcmatrix$microprotein_kegg_biorep),2)),x=c(1,2,3,4,5,6),y=0.8)
text(labels="Median",x=3.5,y=0.95)


BCtech <- data.frame("human_protein"=techbcmatrix$humanprotein_techrep,"microbial_protein"=techbcmatrix$microprotein_techrep,"genus"=techbcmatrix$filter5_genus_techrep,"species"=techbcmatrix$filter5_species_techrep,"cog"=techbcmatrix$microprotein_cog_techrep,"kegg"=techbcmatrix$microprotein_kegg_techrep)
vioplot(BCtech$human_protein,BCtech$microbial_protein,BCtech$genus,BCtech$species,BCtech$cog,BCtech$kegg,
        col=c("#8bc63e","#d696a4","#348083","#877ea7","#f89b31","#7895a3"),names = c("Human Protein","Microbial Protein","Genus","Species","COG","KEGG"),main = "Intra-batch technical replicates",ylim=c(0,1))
text(labels=c(round(median(techbcmatrix$humanprotein_techrep),2),round(median(techbcmatrix$microprotein_techrep),2),round(median(techbcmatrix$filter5_genus_techrep),2),round(median(techbcmatrix$filter5_species_techrep),2),round(median(techbcmatrix$microprotein_cog_techrep),2),round(median(techbcmatrix$microprotein_kegg_techrep),2)),x=c(1,2,3,4,5,6),y=0.8)
text(labels="Median",x=3.5,y=0.95)


BCqc <- data.frame("human_protein"=qcbcmatrix$humanprotein_qc,"microbial_protein"=qcbcmatrix$microprotein_qc,"genus"=qcbcmatrix$filter5_genus_qc,"species"=qcbcmatrix$filter5_species_qc,"cog"=qcbcmatrix$microprotein_cog_qc,"kegg"=qcbcmatrix$microprotein_kegg_qc)
vioplot(BCqc$human_protein,BCqc$microbial_protein,BCqc$genus,BCqc$species,BCqc$cog,BCqc$kegg,
        col=c("#8bc63e","#d696a4","#348083","#877ea7","#f89b31","#7895a3"),names = c("Human Protein","Microbial Protein","Genus","Species","COG","KEGG"),main = "Inter-batch biological replicates",ylim=c(0,1))
text(labels=c(round(median(qcbcmatrix$humanprotein_qc),2),round(median(qcbcmatrix$microprotein_qc),2),round(median(qcbcmatrix$filter5_genus_qc),2),round(median(qcbcmatrix$filter5_species_qc),2),round(median(qcbcmatrix$microprotein_cog_qc),2),round(median(qcbcmatrix$microprotein_kegg_qc),2)),x=c(1,2,3,4,5,6),y=0.8)
text(labels="Median",x=3.5,y=0.95)

BCpool <- data.frame("human_protein"=poolbcmatrix$humanprotein_pool,"microbial_protein"=poolbcmatrix$microprotein_pool,"genus"=poolbcmatrix$filter5_genus_pool,"species"=poolbcmatrix$filter5_species_pool,"cog"=poolbcmatrix$microprotein_cog_pool,"kegg"=poolbcmatrix$microprotein_kegg_pool)
vioplot(BCpool$human_protein,BCpool$microbial_protein,BCpool$genus,BCpool$species,BCpool$cog,BCpool$kegg,
        col=c("#8bc63e","#d696a4","#348083","#877ea7","#f89b31","#7895a3"),names = c("Human Protein","Microbial Protein","Genus","Species","COG","KEGG"),main = "Inter-batch technical replicates",ylim=c(0,1))
text(labels=c(round(median(poolbcmatrix$humanprotein_pool),2),round(median(poolbcmatrix$microprotein_pool),2),round(median(poolbcmatrix$filter5_genus_pool),2),round(median(poolbcmatrix$filter5_species_pool),2),round(median(poolbcmatrix$microprotein_cog_pool),2),round(median(poolbcmatrix$microprotein_kegg_pool),2)),x=c(1,2,3,4,5,6),y=0.8)
text(labels="Median",x=3.5,y=0.95)

dev.off()

##need to add permanova####


###fig s1a####
pdf("Rep_Cor.pdf",width=10,height=4)
par(mfrow=c(1,4))

corrbio <- data.frame("human_protein"=biocorrmatrix$humanprotein_biorep,"microbial_protein"=biocorrmatrix$microprotein_biorep,"genus"=biocorrmatrix$filter5_genus_biorep,"species"=biocorrmatrix$filter5_species_biorep,"cog"=biocorrmatrix$microprotein_cog_biorep,"kegg"=biocorrmatrix$microprotein_kegg_biorep)
vioplot(corrbio$human_protein,corrbio$microbial_protein,corrbio$genus,corrbio$species,corrbio$cog,corrbio$kegg,
        col=c("#8bc63e","#d696a4","#348083","#877ea7","#f89b31","#7895a3"),names = c("Human Protein","Microbial Protein","Genus","Species","COG","KEGG"),main = "Intra-batch biological replicates",ylim=c(0,1))
text(labels=c(round(median(biocorrmatrix$humanprotein_biorep,na.rm=T),2),round(median(biocorrmatrix$microprotein_biorep,na.rm=T),2),round(median(biocorrmatrix$filter5_genus_biorep,na.rm=T),2),round(median(biocorrmatrix$filter5_species_biorep,na.rm=T),2),round(median(biocorrmatrix$microprotein_cog_biorep,na.rm=T),2),round(median(biocorrmatrix$microprotein_kegg_biorep,na.rm=T),2)),x=c(1,2,3,4,5,6),y=0.8)
text(labels="Median",x=3.5,y=0.95)

corrtech <- data.frame("human_protein"=techcorrmatrix$humanprotein_techrep,"microbial_protein"=techcorrmatrix$microprotein_techrep,"genus"=techcorrmatrix$filter5_genus_techrep,"species"=techcorrmatrix$filter5_species_techrep,"cog"=techcorrmatrix$microprotein_cog_techrep,"kegg"=techcorrmatrix$microprotein_kegg_techrep)
vioplot(corrtech$human_protein,corrtech$microbial_protein,corrtech$genus,corrtech$species,corrtech$cog,corrtech$kegg,
        col=c("#8bc63e","#d696a4","#348083","#877ea7","#f89b31","#7895a3"),names = c("Human Protein","Microbial Protein","Genus","Species","COG","KEGG"),main = "Intra-batch technical replicates",ylim=c(0,1))
text(labels=c(round(median(techcorrmatrix$humanprotein_techrep),2),round(median(techcorrmatrix$microprotein_techrep),2),round(median(techcorrmatrix$filter5_genus_techrep),2),round(median(techcorrmatrix$filter5_species_techrep),2),round(median(techcorrmatrix$microprotein_cog_techrep),2),round(median(techcorrmatrix$microprotein_kegg_techrep),2)),x=c(1,2,3,4,5,6),y=0.8)
text(labels="Median",x=3.5,y=0.95)

corrqc <- data.frame("human_protein"=qccorrmatrix$humanprotein_qc,"microbial_protein"=qccorrmatrix$microprotein_qc,"genus"=qccorrmatrix$filter5_genus_qc,"species"=qccorrmatrix$filter5_species_qc,"cog"=qccorrmatrix$microprotein_cog_qc,"kegg"=qccorrmatrix$microprotein_kegg_qc)
vioplot(corrqc$human_protein,corrqc$microbial_protein,corrqc$genus,corrqc$species,corrqc$cog,corrqc$kegg,
        col=c("#8bc63e","#d696a4","#348083","#877ea7","#f89b31","#7895a3"),names = c("Human Protein","Microbial Protein","Genus","Species","COG","KEGG"),main = "Inter-batch biological replicates",ylim=c(0,1))
text(labels=c(round(median(qccorrmatrix$humanprotein_qc),2),round(median(qccorrmatrix$microprotein_qc),2),round(median(qccorrmatrix$filter5_genus_qc),2),round(median(qccorrmatrix$filter5_species_qc),2),round(median(qccorrmatrix$microprotein_cog_qc),2),round(median(qccorrmatrix$microprotein_kegg_qc),2)),x=c(1,2,3,4,5,6),y=0.8)
text(labels="Median",x=3.5,y=0.95)

corrpool <- data.frame("human_protein"=poolcorrmatrix$humanprotein_pool,"microbial_protein"=poolcorrmatrix$microprotein_pool,"genus"=poolcorrmatrix$filter5_genus_pool,"species"=poolcorrmatrix$filter5_species_pool,"cog"=poolcorrmatrix$microprotein_cog_pool,"kegg"=poolcorrmatrix$microprotein_kegg_pool)
vioplot(corrpool$human_protein,corrpool$microbial_protein,corrpool$genus,corrpool$species,corrpool$cog,corrpool$kegg,
        col=c("#8bc63e","#d696a4","#348083","#877ea7","#f89b31","#7895a3"),names = c("Human Protein","Microbial Protein","Genus","Species","COG","KEGG"),main = "Inter-batch technical replicates",ylim=c(0,1))
text(labels=c(round(median(poolcorrmatrix$humanprotein_pool),2),round(median(poolcorrmatrix$microprotein_pool),2),round(median(poolcorrmatrix$filter5_genus_pool),2),round(median(poolcorrmatrix$filter5_species_pool),2),round(median(poolcorrmatrix$microprotein_cog_pool),2),round(median(poolcorrmatrix$microprotein_kegg_pool),2)),x=c(1,2,3,4,5,6),y=0.8)
text(labels="Median",x=3.5,y=0.95)

dev.off()



##fig 1e#####
library(vioplot)
corrtotal <- rbind(corrpool,corrqc,corrtech,corrbio)
corrtotal <- na.omit(corrtotal)
pdf("total_qc_rep.pdf",width=8,height=4)
vioplot(corrtotal$human_protein,corrtotal$microbial_protein,corrtotal$genus,corrtotal$species,corrtotal$kegg,corrtotal$cog,
        col=c("#8bc63e","#d696a4","#348083","#877ea7","#f89b31","#7895a3"),names = c("Human Protein","Microbial Protein","Genus","Species","KEGG","COG"),main ="QC replicates",ylim=c(0,1))
text(labels=c(round(median(corrtotal$human_protein),2),round(median(corrtotal$microbial_protein),2),round(median(corrtotal$genus),2),round(median(corrtotal$species),2),round(median(corrtotal$kegg),2),round(median(corrtotal$cog),2)),x=c(1,2,3,4,5,6),y=0.8)
text(labels="Median",x=3.5,y=0.95)
dev.off()

###fig1e 2,random sample###
human <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_humanprotein_sample.tsv",sep="\t")
micro <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_microprotein_sample.tsv")
genus <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_filter5_genus_sample.txt")
species <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_filter5_species_sample.txt")
cog <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_microprotein_cog_sample.txt")
kegg <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_microprotein_kegg_sample.txt")

sample.comb <- as.data.frame(combn(colnames(genus)[2:ncol(genus)],2))
sample.comb2 <- sample(sample.comb,500)

get.corr.bc <- function(x){
  biocorrx <- c()
  bcdisx <- c()
  for(i in 1:500){
  xx <- x[,c(which(colnames(x) %in% sample.comb2[,i]))]
  biocorrx[i] <- cor(xx[,1],xx[,2],use ="pairwise.complete.obs",method = "spearman")
  xx[is.na(xx)]=0
  bcdisx[i] <- vegan::vegdist(t(xx),method="bray")
  }
  return(list("corr"=biocorrx,"bc"=bcdisx))
}

human.list <- get.corr.bc(human)
micro.list <- get.corr.bc(micro)
genus.list <- get.corr.bc(genus)
species.list <- get.corr.bc(species)
cog.list <- get.corr.bc(cog)
kegg.list <- get.corr.bc(kegg)

library(vioplot)
pdf("random_rep.pdf",width=8,height=4)
vioplot(human.list$corr,micro.list$corr,genus.list$corr,species.list$corr,kegg.list$corr,cog.list$corr,
        col=c("#8bc63e","#d696a4","#348083","#877ea7","#f89b31","#7895a3"),names = c("Human Protein","Microbial Protein","Genus","Species","KEGG","COG"),main ="QC replicates",ylim=c(0,1))
text(labels=c(round(median(human.list$corr),2),round(median(micro.list$corr),2),round(median(genus.list$corr),2),round(median(species.list$corr),2),round(median(kegg.list$corr),2),round(median(cog.list$corr),2)),x=c(1,2,3,4,5,6),y=0.8)
text(labels="Median",x=3.5,y=0.95)
dev.off()

pdf("random_dis.pdf",width=8,height=4)
vioplot(human.list$bc,micro.list$bc,genus.list$bc,species.list$bc,kegg.list$bc,cog.list$bc,
        col=c("#8bc63e","#d696a4","#348083","#877ea7","#f89b31","#7895a3"),names = c("Human Protein","Microbial Protein","Genus","Species","KEGG","COG"),main ="QC replicates",ylim=c(0,1))
text(labels=c(round(median(human.list$bc),2),round(median(micro.list$bc),2),round(median(genus.list$bc),2),round(median(species.list$bc),2),round(median(kegg.list$bc),2),round(median(cog.list$bc),2)),x=c(1,2,3,4,5,6),y=0.8)
text(labels="Median",x=3.5,y=0.95)
dev.off()


####figs1c#########PCoA########

rm(list=ls())

standardize_ident <- function(x) {
  if(length(grep("^b",x)>0)){
    parts <- strsplit(x, "_")
    batch <- sapply(parts, `[`, 1)
    position <- sapply(parts, `[`, 2)
    batch <- sub("^b", "", batch)  
    batch <- sprintf("b%02d", as.integer(batch))  
    pos_num <- sub("^([0-9]+)([a-zA-Z]*)$", "\\1", position)  
    pos_text <- sub("^([0-9]+)([a-zA-Z]*)$", "\\2", position)  
    pos_num <- sprintf("%02d", as.integer(pos_num)) 
    position <- paste0(pos_num, pos_text)
    rename <- paste0(batch, "_", position)
  }else{
    rename <- x
  }
}

library(data.table)
library(Matrix)
library(vegan)
library(ape)
library(ggplot2)


#metadata
metadata <- read.delim("~/Documents/gnhsf/data202509/metadata/GNHSF_sample_ident_20220307.csv",sep=",")
metadata2 <- read.delim("~/Documents/gnhsf/data202509/metadata/GNHSF_pool_QC_20210922.csv",sep=",")
metadata3 <- merge(metadata2,metadata,by="batch",all.y=T)
metadata3$group <- ifelse(metadata3$label1=="pool",metadata3$pool,ifelse(metadata3$label1=="qc",metadata3$QC,""))
# input all micro and human proteins
pepall <- fread("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/all/GNHSF_diann_IGC_humanswiss_protein_all.tsv", sep = "\t")
rownames(pepall) <- pepall$Protein.Group

aa <- apply(pepall[, 2:ncol(pepall)], 1, function(x) length(which(is.na(x))) / 2514)
#pepall1 <- pepall[-which(aa > 0.8), ]
pepall1 <- pepall

pepall1[is.na(pepall1)] <- 0

pepall1_t <- data.frame(t(pepall1[, -1]))
pepall1_t$sample <- rownames(pepall1_t)

pepall1_t$sample <- sapply(pepall1_t$sample, standardize_ident)

pep_data <- merge(metadata3[, c(9,ncol(metadata3))], pepall1_t, by = "sample")

pep_data2 <- pep_data[, 3:(ncol(pep_data))]
rownames(pep_data2) <- pep_data$sample
empty_rows <- rowSums(pep_data2 != 0 & !is.na(pep_data2)) == 0
pep_data2_clean <- pep_data2[!empty_rows, ]
batch_info_clean <- pep_data$group [!empty_rows]
pep_data2_clean <- pep_data2_clean[, colSums(pep_data2_clean != 0 & !is.na(pep_data2_clean)) > 0]

pep_data2_clean <- as.matrix(pep_data2_clean)

cat("Total samples:", nrow(pep_data2_clean), "\n")
cat("Total features:", ncol(pep_data2_clean), "\n")

library(parallelDist)

  dist_matrix <- parallelDist(as.matrix(pep_data2_clean), 
                              method = "bray", 
                              threads = parallel::detectCores() - 1)

dist <- as.matrix(dist_matrix)
pcoa_result <- ape::pcoa(dist, correction = "none")

pcoa_coords <- pcoa_result$vectors[, 1:3]
explained_var <- pcoa_result$values$Relative_eig[1:3]
library(ggthemes)
pcoa_data <- data.frame(
  PCo1 = pcoa_coords[, 1],
  PCo2 = pcoa_coords[, 2],
  PCo3 = pcoa_coords[, 3],
  batch = as.factor(batch_info_clean)
)

p_pcoa_12 <- ggplot(pcoa_data, aes(x = PCo1, y = PCo2, color = batch)) +
  geom_point(size = 1.2) +
  labs(
    title = "PCoA - Batch Effect Assessment",
    x = sprintf("PCo1 (%.1f%%)", explained_var[1] * 100),
    y = sprintf("PCo2 (%.1f%%)", explained_var[2] * 100)
  ) + scale_color_brewer(palette = "Paired")+
  theme_base() +stat_ellipse(data=pcoa_data,
                             geom = "polygon",level=0.9,
                             linetype = 1,size=0.8,
                             aes(fill=batch),
                             alpha=0.1,
                             show.legend = T)+
  scale_fill_brewer(palette = "Paired")
print(p_pcoa_12)
