setwd("~/Documents/gnhsf/data202509/code")



##getBC


library(data.table)
get_bc <- function(x){
data.list <- list.files(x)
dir.create(paste0(x,"_BC"))
for (i in 1:length(data.list)) {
data.test <- fread(paste0(x,"/",data.list[i]),sep="\t")
data.test <- data.test[,-1]
data.test[is.na(data.test)] <- 0

library(vegan)
bc.data <- as.matrix(vegdist(t(data.test),method="bray",binary=F))
write.table(bc.data,paste0(x,"_BC/BC_",data.list[i]),sep="\t")
cat("\r",data.list[i],i,"/",length(data.list))
}
}

get_bc("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/biorep")
get_bc("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/techrep")
get_bc("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/pool")
get_bc("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/qc")
