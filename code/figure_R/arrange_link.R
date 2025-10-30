###arrange link####
rm(list=ls())


library(stringr)
##generate link for phylum###
link <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microbiome_all.txt",sep="\t")
link.to.phylum <- link[grep("p__",link$Taxon_all_unipept),]
link.to.phylum$Taxon_rank_unipept="phylum"
link.to.phylum$Taxon_name_unipept <- data.frame(str_split_fixed(link.to.phylum$Taxon_all_unipept,"[|]|__",20))$X6
link2 <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microbiome_all_species2genus.txt",sep="\t")
link2 <- rbind(link2,link.to.phylum)
 write.table(link2,"~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microbiome_all_species2genus.txt",sep="\t") 
  

 
 rm(list=ls())
 ####extract each level#####
link <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microbiome_all_species2genus.txt",sep="\t")

#species link
link.species <- subset(link,Taxon_rank_unipept=="species")
link.species2 <- unique(link.species[,c(2,17)])
protein.taxon.freq <- data.frame(table(link.species2$Protein))
length(which(protein.taxon.freq$Freq==1))
#total42826 protein, 41133 proteins with unique species
protein.species <- protein.taxon.freq$Var1[which(protein.taxon.freq$Freq==1)]
link.species3 <- unique(link.species[which(link.species$Protein %in% protein.species),c(-1,-18)])
#link.species3 <- link.species3[-which(duplicated(link.species3$Protein)),]
write.table(link.species3,"microbial_prot_to_unique_species_link.txt",sep="\t",row.names = F)



##genus link
link.genus <- subset(link,Taxon_rank_unipept=="genus")
link.genus2 <- unique(link.genus[,c(2,17)])
protein.taxon.freq <- data.frame(table(link.genus2$Protein))
length(which(protein.taxon.freq$Freq==1))
#total 56827 protein, 56135 proteins with unique genus
protein.genus <- protein.taxon.freq$Var1[which(protein.taxon.freq$Freq==1)]
link.genus3 <- unique(link.genus[which(link.genus$Protein %in% protein.genus),c(-1,-18)])
link.genus3 <- link.genus3[-which(duplicated(link.genus3$Protein)),]
write.table(link.genus3,"microbial_prot_to_unique_genus_link.txt",sep="\t",row.names = F)


##phylum link
link.phylum <- subset(link,Taxon_rank_unipept=="phylum")
link.phylum2 <- unique(link.phylum[,c(2,17)])
protein.taxon.freq <- data.frame(table(link.phylum2$Protein))
length(which(protein.taxon.freq$Freq==1))
#total 69519 protein, 69383  proteins with unique phylum
protein.phylum <- protein.taxon.freq$Var1[which(protein.taxon.freq$Freq==1)]
link.phylum3 <- unique(link.phylum[which(link.phylum$Protein %in% protein.phylum),c(-1,-18)])
link.phylum3 <- link.phylum3[-which(duplicated(link.phylum3$Protein)),]
write.table(link.phylum3,"microbial_prot_to_unique_phylum_link.txt",sep="\t",row.names = F)



##all to the depth #### not LCA, use the unqiue and depth####
library(stringr)
link <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microbiome_all.txt",sep="\t")
link1 <- unique(link[,-1])
freq <- data.frame(table(link1$Protein))
freq <- freq[-which(freq$Freq==1),]
prot.list <- as.character(freq$Var1)

unq.link <- link1[-which(link1$Protein %in% prot.list),]
for (x in prot.list) {
  data=link1[which(link1$Protein==x),]
    num <- sapply(data$Taxon_all_unipept,function(x) stringr::str_count(x,pattern = "[|]"))
    while (length(which(num == max(num, na.rm = TRUE))) > 1) {
      if (all(num == 0, na.rm = TRUE) || max(num, na.rm = TRUE) <= 0) {
        break
      }
            num[which(num == max(num, na.rm = TRUE))] <- 0
    }  
    if(all(num==0)){
      final.data <-data[1,]
      final.data$Taxon_rank_unipept <- NA
      final.data$Taxon_all_unipept <- NA
      final.data$Taxon_name_unipept <- NA
      final.data$Taxon_id_unipept <- NA
    }else{
    final.data <- data[which(num==max(num)),]}
unq.link <- rbind(unq.link,final.data)
cat("\r",which(prot.list==x))
}

link2 <- unq.link
write.xlsx(link2,"~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microprot_all_arrange20250916.xlsx")



##to phylum,to genus, tospecies$$$
library(stringr)
link2 <- readxl::read_xlsx("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microprot_all_arrange20250916.xlsx")
link2[link2==""] <- NA
data.total <- data.frame()
prot.list <- link2$Protein
for (x in na.omit(prot.list)) {
  data=link2[which(link2$Protein==x),]
  lca <- data.frame(str_split_fixed(data$Taxon_all_unipept,"[|]",9))
  lca[lca==""] <- NA
  colnames(lca) <- c("domain","kingdom", "phylum" , "class", "order", "family", "genus", "species","strain")
  if(!is.na(lca$species)[1]){
    data$Taxon_name_unipept<- strsplit(lca$species[1],"__",2)[[1]][2]
    data$Taxon_rank_unipept <- "species"
  }else if(!is.na(lca$genus)[1]){
    data$Taxon_name_unipept<- strsplit(lca$genus[1],"__",2)[[1]][2]
    data$Taxon_rank_unipept <- "genus"
    }else if(!is.na(lca$phylum)[1]){
    data$Taxon_name_unipept<- strsplit(lca$phylum[1],"__",2)[[1]][2]
    data$Taxon_rank_unipept <- "phylum"
    }else{
      data$Taxon_name_unipept<- NA
      data$Taxon_rank_unipept <- NA
    }
  data.total <- rbind(data.total,data)
  cat("\r",which(prot.list==x))
}
library(openxlsx)
write.xlsx(data.total,"~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microprot_2phylum2genus2species_one20250916.xlsx")


###link2####
##to phylum,to genus, tospecies$$$
link2 <- readxl::read_xlsx("~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microprot_all_arrange20250916.xlsx")
link2[link2==""] <- NA
data.total <- data.frame()
prot.list <- link2$Protein
for (x in na.omit(prot.list)) {
  data=link2[which(link2$Protein==x),]
  lca <- data.frame(str_split_fixed(data$Taxon_all_unipept,"[|]",9))
  lca[lca==""] <- NA
  colnames(lca) <- c("domain","kingdom", "phylum" , "class", "order", "family", "genus", "species","strain")
  data$Taxon_name_unipept <- paste(lca$phylum,lca$genus,lca$species,sep = ";") 
  data$Taxon_rank_unipept <- paste("phylum","genus","species",sep = ";") 
  data.total <- rbind(data.total,data)
  cat("\r",which(prot.list==x))
}
write.xlsx(data.total,"~/Documents/gnhsf/data202509/metaproteomics/08.link/GNHSF_link_microprot_2phylum2genus2species_all20250916.xlsx")
