rm(list=ls())

###########metaproteome#######
########cross 1385#######

cross.phylum <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_filter5_phylum_sample_1385.txt",sep="\t")
cross.genus <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_filter5_genus_sample_1385.txt",sep="\t")
cross.species <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_filter5_species_sample_1385.txt",sep="\t")
cross.cog <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_microprotein_cog_sample_1385.txt",sep="\t")
cross.kegg <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_microprotein_kegg_sample_1385.txt",sep="\t")
cross.metaprot <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_microprotein_sample_1385.txt",sep="\t")
cross.humanprot <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385/GNHSF_diann_IGC_humanswiss_humanprotein_sample_1385.txt",sep="\t")

rownames(cross.phylum) <- cross.phylum[,1]
cross.phylum <- cross.phylum[,-1]

rownames(cross.genus) <- cross.genus[,1]
cross.genus <- cross.genus[,-1]

rownames(cross.species) <- cross.species[,1]
cross.species <- cross.species[,-1]

rownames(cross.cog) <- cross.cog[,1]
cross.cog <- cross.cog[,-1]

rownames(cross.kegg) <- cross.kegg[,1]
cross.kegg <- cross.kegg[,-1]

rownames(cross.metaprot) <- cross.metaprot[,1]
cross.metaprot <- cross.metaprot[,-1]

rownames(cross.humanprot) <- cross.humanprot[,1]
cross.humanprot <- cross.humanprot[,-1]

cross.phylum$prevalence <- apply(cross.phylum,1,function(x) sum(!is.na(x))/1385)
crosscore.phylum80 <- rownames(cross.phylum)[which(cross.phylum$prevalence>=0.8)]

cross.genus$prevalence <- apply(cross.genus,1,function(x) sum(!is.na(x))/1385)
crosscore.genus80 <- rownames(cross.genus)[which(cross.genus$prevalence>=0.8)]

cross.species$prevalence <- apply(cross.species,1,function(x) sum(!is.na(x))/1385)
crosscore.species80 <- rownames(cross.species)[which(cross.species$prevalence>=0.8)]

cross.cog$prevalence <- apply(cross.cog,1,function(x) sum(!is.na(x))/1385)
crosscore.cog80 <- rownames(cross.cog)[which(cross.cog$prevalence>=0.8)]

cross.kegg$prevalence <- apply(cross.kegg,1,function(x) sum(!is.na(x))/1385)
crosscore.kegg80 <- rownames(cross.kegg)[which(cross.kegg$prevalence>=0.8)]
crossvariable.kegg80 <- rownames(cross.kegg)[which(cross.kegg$prevalence<0.8)]

cross.metaprot$prevalence <- apply(cross.metaprot,1,function(x) sum(!is.na(x))/1385)
crosscore.metaprot80 <- rownames(cross.metaprot)[which(cross.metaprot$prevalence>=0.8)]

cross.humanprot$prevalence <- apply(cross.humanprot,1,function(x) sum(!is.na(x))/1385)
crosscore.humanprot80 <- rownames(cross.humanprot)[which(cross.humanprot$prevalence>=0.8)]

###long 954######
long.phylum <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_954/GNHSF_diann_IGC_humanswiss_filter5_phylum_sample_954.txt",sep="\t")
long.genus <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_954/GNHSF_diann_IGC_humanswiss_filter5_genus_sample_954.txt",sep="\t")
long.species <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_954/GNHSF_diann_IGC_humanswiss_filter5_species_sample_954.txt",sep="\t")
long.cog <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_954/GNHSF_diann_IGC_humanswiss_microprotein_cog_sample_954.txt",sep="\t")
long.kegg <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_954/GNHSF_diann_IGC_humanswiss_microprotein_kegg_sample_954.txt",sep="\t")
long.metaprot <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_954/GNHSF_diann_IGC_humanswiss_microprotein_sample_954.txt",sep="\t")
long.humanprot <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_954/GNHSF_diann_IGC_humanswiss_humanprotein_sample_954.txt",sep="\t")

time.sample <- intersect(colnames(cross.cog),colnames(long.cog))


rownames(long.phylum) <- long.phylum[,1]
long.phylum <- long.phylum[,-1]
long.phylum <- long.phylum[,-which(colnames(long.phylum) %in% time.sample)]

rownames(long.genus) <- long.genus[,1]
long.genus <- long.genus[,-1]
long.genus <- long.genus[,-which(colnames(long.genus) %in% time.sample)]

rownames(long.species) <- long.species[,1]
long.species <- long.species[,-1]
long.species <- long.species[,-which(colnames(long.species) %in% time.sample)]

rownames(long.cog) <- long.cog[,1]
long.cog <- long.cog[,-1]
long.cog <- long.cog[,-which(colnames(long.cog) %in% time.sample)]

rownames(long.kegg) <- long.kegg[,1]
long.kegg <- long.kegg[,-1]
long.kegg <- long.kegg[,-which(colnames(long.kegg) %in% time.sample)]

rownames(long.metaprot) <- long.metaprot[,1]
long.metaprot <- long.metaprot[,-1]
long.metaprot <- long.metaprot[,-which(colnames(long.metaprot) %in% time.sample)]

rownames(long.humanprot) <- long.humanprot[,1]
long.humanprot <- long.humanprot[,-1]
long.humanprot <- long.humanprot[,-which(colnames(long.humanprot) %in% time.sample)]

long.phylum$prevalence <- apply(long.phylum,1,function(x) sum(!is.na(x))/477)
longcore.phylum80 <- rownames(long.phylum)[which(long.phylum$prevalence>=0.8)]
longvariable.phylum80 <- rownames(long.phylum)[which(long.phylum$prevalence<0.8)]

long.genus$prevalence <- apply(long.genus,1,function(x) sum(!is.na(x))/477)
longcore.genus80 <- rownames(long.genus)[which(long.genus$prevalence>=0.8)]

long.species$prevalence <- apply(long.species,1,function(x) sum(!is.na(x))/477)
longcore.species80 <- rownames(long.species)[which(long.species$prevalence>=0.8)]

long.cog$prevalence <- apply(long.cog,1,function(x) sum(!is.na(x))/477)
longcore.cog80 <- rownames(long.cog)[which(long.cog$prevalence>=0.8)]

long.kegg$prevalence <- apply(long.kegg,1,function(x) sum(!is.na(x))/477)
longcore.kegg80 <- rownames(long.kegg)[which(long.kegg$prevalence>=0.8)]

long.metaprot$prevalence <- apply(long.metaprot,1,function(x) sum(!is.na(x))/477)
longcore.metaprot80 <- rownames(long.metaprot)[which(long.metaprot$prevalence>=0.8)]


long.humanprot$prevalence <- apply(long.humanprot,1,function(x) sum(!is.na(x))/477)
longcore.humanprot80 <- rownames(long.humanprot)[which(long.humanprot$prevalence>=0.8)]



####overlap#####
library(qpcR)
core.fea80 <- data.frame( qpcR:::cbind.na("core_phylum"=intersect(longcore.phylum80,crosscore.phylum80),
                                          "core_genus"=intersect(longcore.genus80,crosscore.genus80),
                                          "core_species"=intersect(longcore.species80,crosscore.species80),
                                          "core_cog"=intersect(longcore.cog80,crosscore.cog80),
                                          "core_kegg"=intersect(longcore.kegg80,crosscore.kegg80),
                                          "core_microprot"=intersect(longcore.metaprot80,crosscore.metaprot80),
                                          "corehumanprot"=intersect(longcore.humanprot80,crosscore.humanprot80)  ))
write.table(core.fea80,"gnhsf_metaproteome_core_fea_prevalence80_20250909.txt",sep="\t",row.names = F)


##fig s5a####
pdf("corefea_bar.pdf",width=6,height=5)
barplot(height=c(length(na.omit(core.fea80$corehumanprot)),length(na.omit(core.fea80$core_microprot)),length(na.omit(core.fea80$core_phylum)),length(na.omit(core.fea80$core_genus)),length(na.omit(core.fea80$core_species)),length(na.omit(core.fea80$core_kegg)),length(na.omit(core.fea80$core_cog))),
names.arg = c("Human protein","Microbial protein","Phylum","Genus","Species","KO","COG")  ,
col = c("#A8c5e6","#85bf53","#d696a4","#348083","#877ea7","#f89b31","#7895a3"),ylab =" # of core features")
text(labels=c(length(na.omit(core.fea80$corehumanprot)),length(na.omit(core.fea80$core_microprot)),length(na.omit(core.fea80$core_phylum)),length(na.omit(core.fea80$core_genus)),length(na.omit(core.fea80$core_species)),length(na.omit(core.fea80$core_kegg)),length(na.omit(core.fea80$core_cog))),
     x=c(1:7),
     y=c(length(na.omit(core.fea80$corehumanprot)),length(na.omit(core.fea80$core_microprot)),length(na.omit(core.fea80$core_phylum)),length(na.omit(core.fea80$core_genus)),length(na.omit(core.fea80$core_species)),length(na.omit(core.fea80$core_kegg)),length(na.omit(core.fea80$core_cog))))

dev.off()


##fig s5b####
rm(list=ls())
core.fea80 <- read.delim("gnhsf_metaproteome_core_fea_prevalence80_20250909.txt",sep="\t")
library(stringr)
genussabu <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_filter5_genus_sample.txt",sep="\t")
genusabun <- rowSums(genussabu[,2:ncol(genussabu)],na.rm = T)
genusabunnamec <- data.frame(str_split_fixed(genussabu$Taxon,"[|]",10))$X7
genus <- data.frame("genus"=genusabunnamec,"logFC"=genusabun)
coretaxon <- data.frame(str_split_fixed(na.omit(core.fea80$core_genus),"[|]",10))

phylumcore <- unique(coretaxon$X3)
genuscore <- unique(coretaxon$X7)
taxon <- data.frame(matrix(ncol=length(genuscore),nrow=length(phylumcore),dimnames = list(phylumcore,genuscore)))

for (i in 1:length(genuscore)) {
  p <- coretaxon[i,3]
  g <- coretaxon[i,7]
  taxon[p,g] <- 1
}
taxon[is.na(taxon)] <- 0

taxon$total <- rowSums(taxon)
library(dplyr)
taxon <- arrange(taxon,total)
taxon2 <- data.frame(t(taxon))
library(tibble)
taxon3 <- add_column(taxon2,"genus"=rownames(taxon2),.before=1)
taxon4 <- merge(taxon3,genus,by="genus")
rownames(taxon4) <- taxon4[,1]
taxon4 <- taxon4[,-1]
taxon4$logFC <- scale(as.numeric(taxon4$logFC))

write.csv(taxon4,"GNHSF_coregenus_match20220323.csv")
library(GOplot)
##最后一列的列名必须是logFC才能识别

GOChord(taxon4,gene.order = "logFC",lfc.col= c("darkred","lightcoral","white"),gene.size = 3,ribbon.col = brewer.pal(8, "Set3"))
ggsave("GNHSF_coregenus_chord_20220323.pdf")


library(openxlsx)
#### core microbial proteins####
coremicro <- na.omit(core.fea80$core_microprot)
coremicro <- sapply(coremicro, function(x) strsplit(x," ")[[1]][1])
names(coremicro) <- NULL
link.genus<- read.delim("~/Documents/gnhsf/data202509/metaproteomics/08.link/microbial_prot_to_unique_genus_link.txt",sep="\t")
link.genus.core <- link.genus[which(link.genus$Protein %in% coremicro),]
write.xlsx(link.genus.core,"tableS3.xlsx")
link.species<- read.delim("~/Documents/gnhsf/data202509/metaproteomics/08.link/microbial_prot_to_unique_species_link.txt",sep="\t")
link.species.core <- link.species[which(link.species$Protein %in% coremicro),]

link.phylum <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/08.link/microbial_prot_to_unique_phylum_link.txt",sep="\t")
link.phylum.core <- link.phylum[which(link.phylum$Protein %in% coremicro),]
table(link.phylum.core$Taxon_name_unipept)
genus.talbe <- data.frame(table(link.genus.core$Taxon_name_unipept))

coregenus <- data.frame("core.genus"=na.omit(core.fea80$core_genus))
coregenus$core.genus <- sapply(coregenus$core.genus, function(x) strsplit(x,"[|]|__")[[1]][14])
coregenus$core<- "core"

genus.table <- merge(genus.talbe,coregenus,by.x="Var1",by.y="core.genus",all=T)
genus.table <- arrange(genus.table,desc(Freq))
top20.genus <- as.character(genus.table$Var1[1:20])
#figs5c####
library(ggpubr)
link.genus.core.plot <- subset(link.genus.core,Taxon_name_unipept %in% top20.genus)
pdf("core_prot_number_to_genus.pdf")
ggpie::ggpie(link.genus.core.plot,group_key = "Taxon_name_unipept",count_type = "full",label_pos = "out")
dev.off()



##figs5 d####
top8.genus <- as.character(genus.table$Var1[1:8])
library(tidyr)
genus8 <-  link.genus.core[which(link.genus.core$Taxon_name_unipept %in% top8.genus),]
genus8 <- genus8 %>% as_tibble() %>%    separate_rows(COG_NCBIcat_des, sep = ";")
genus8[genus8==""] <- NA
genus8 <- genus8[-which(is.na(genus8$COG_NCBIcat_des)),]
library(ggplot2)
library(ggthemes)
ggplot(genus8,aes(x=COG_NCBIcat_des,fill=Taxon_name_unipept))+geom_bar()+facet_wrap(.~Taxon_name_unipept,nrow = 2,scales = "free_x")+coord_flip()+theme_base()+theme(legend.position = "none")
ggsave("core_genus_to_protein_stats.pdf",width=15)


##fig s6a,b####
keggcore <- na.omit(core.fea80$core_kegg)
cogcore <-  na.omit(core.fea80$core_cog)
cogcat <- read.delim("~/Documents/gnhsf/data202509/metadata/01.NCBI.COG.list.txt",sep="\t",quote="")
cogcore <- cogcat[which(cogcat$COGnum %in% cogcore),]
library(openxlsx)
keggcat <- read.xlsx("~/Documents/gnhsf/data202509/metadata/kegg_database20230328.xlsx")
keggcore <- keggcat[which(keggcat$kegg %in% keggcore),]
keggcore <- unique(keggcore[,c(1,5)])

kegg <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_microprotein_kegg_sample.txt",sep="\t")
cog <- read.delim("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample/GNHSF_diann_IGC_humanswiss_microprotein_cog_sample.txt",sep="\t")
keggabun <- rowSums(kegg[,2:ncol(kegg)],na.rm=T)
cogabun <- rowSums(cog[,2:ncol(cog)],na.rm=T)
keggabun <- data.frame("kegg"=kegg$KEGGnum,"abun"=keggabun)
cogabun <- data.frame("cog"=cog$COGnum,"abun"=cogabun)

keggplot <- merge(keggcore,keggabun,by.x="kegg",by.y = "kegg")
cogplot <- merge(cogcore,cogabun,by.x="COGnum",by.y = "cog")
for (i in 1:nrow(cogplot)) {
  cogplot$NCBICOGcat1[i] <-   ifelse(nchar(cogplot$NCBICOGcat[i])==1,  cogplot$NCBICOGcat[i],
                                     paste(substring(cogplot$NCBICOGcat[i],c(1,2,3),c(1,2,3)),collapse = "_")
  )
}
library(tidyr)
cogplot1 <- separate_rows(cogplot,NCBICOGcat1,convert=T)
cogplot2 <- subset(cogplot1,NCBICOGcat1 != "")
library(ggplot2)

ggplot(keggplot,aes(x=A,y=abun,color=A))+geom_point(aes(size=scale(abun)),position = "jitter")+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(legend.position = "none")+labs(x="",y="Abundance")
ggsave("GNHSF_core_kegg_A_20250323.pdf",width=20,height=10)

ggplot(cogplot2,aes(x=NCBICOGcat1,y=abun,color=NCBICOGcat1))+geom_point(aes(size=scale(abun)),position="jitter")+theme_classic()+
  theme(legend.position = "none")+labs(x="",y="Abundance")
ggsave("GNHSF_core_cog20250323.pdf",width=20,height=10)

##figs6c####
library(MicrobiomeProfiler)
core.genekegg.enrich <- enrichKO(na.omit(core.fea80$core_kegg))
corekegg.enrich <- core.genekegg.enrich@result
corekegg.enrich$is.sig  <- ifelse(corekegg.enrich$FoldEnrichment>2 & corekegg.enrich$qvalue<0.05,"sig","nsig")
table(corekegg.enrich$is.sig)
write.csv(corekegg.enrich,"tableS4.csv")
# nsig  sig 
# 138   61 
pdf("corekegg_pathway.pdf",height=15,width=8)
dotplot(core.genekegg.enrich,showCategory=20)
dev.off()
