rm(list = ls())
setwd("E:/sunyingying/1111LS-GNSHF-project/1111data/01.IGC_humanswiss/01.GNHSF/09.figure/Figure4_permanova_maaslin/permanova_20250908/heatmap")
library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyverse)
library(forcats)
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
#install.packages("ggsci")
library(ggsci)
#install.packages("paletteer")
library(paletteer)
#display.brewer.all()
library("scales")
show_col(pal_nejm("default")(8))
type <- "cog"
func.psig <- function(type){
  df <- read.csv(file = paste0("PERMANOVA_", type, "_new_rmNA_extract_info_addadjp_sort.csv"))
  psig <- c()
  for (i in 1:nrow(df)){
    if (df$P.value[i] < 0.001){
      psig[i] <- "***"
    }else if (df$P.value[i] < 0.01){
      psig[i] <- "**"
    }else if (df$P.value[i] < 0.05){
      psig[i] <- "*"
    }else{
      psig[i] <- "nsig"
    }
  }
  df$psig <- psig
  write.csv(df, file = paste0("PERMANOVA_", type, "_new_rmNA_extract_info_addadjp_sort_psig.csv"), row.names = F)
}
func.psig("cog")
func.psig("kegg")
func.psig("humanprotein")
func.psig("microprotein")
func.psig("genus")
func.psig("species")

func.plot <- function(type){
  permanov<-read.csv(paste0("PERMANOVA_", type, "_new_rmNA_extract_info_addadjp_sort.csv"))
  cat<- rev(permanov$Cinfo)
  permanov$Cinfo = factor(permanov$Cinfo, levels = cat)
  head(permanov)
  p1 <- permanov %>%
    ggplot(aes(x=Cinfo, y=R2_adjusted, fill=Cinfo_type)) +
    geom_bar(stat="identity", width=.4, alpha=.6) +
    coord_flip() +
    xlab("") +
    theme_classic() + scale_fill_nejm() + 
    theme(axis.title.x = element_text(size = 25,face = "bold"), axis.text.x  = element_text(size=40), 
          axis.text.y  = element_text(size=30))+
    theme(legend.title = element_text(size=30, face = "bold"), legend.text = element_text(size=25, face = "bold"))+
    geom_col(width = 0.8, position = position_dodge(width=0.8))
  
  p1
  #pvalue <- permanov$`P-value`
  #adjustp <- p.adjust(pvalue, method = "BH", n=length(pvalue))
  #new<-cbind(permanov, adjustp)
  
  ggsave(paste0("PERMANOVA_", type, "_new_classify_rmNA.png"), p1, width = 20, height = 10, dpi=600)
  ggsave(paste0("PERMANOVA_", type, "_new_classify_rmNA.pdf"), p1, width = 20, height = 10)
  
}
func.plot("cog")
func.plot("kegg")
func.plot("humanprotein")
func.plot("microprotein")
func.plot("genus")
func.plot("species")

