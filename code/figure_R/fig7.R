setwd("~/Documents/gnhsf/data202509/metaproteomics/07.matrix/sample_1385")

rm(list=ls())


metadata <- read.csv("~/Documents/gnhsf/data202509/metadata/GNHSF_sample_inform_normalized_1385.csv")

#fig7a####

species <- read.delim("GNHSF_diann_IGC_humanswiss_filter5_species_sample_1385.txt",sep="\t")
library(dplyr)
species2 <- species %>%
  `rownames<-`(.[,1]) %>%
  select(-1) %>%
  t() %>%
  as.data.frame()

library(tibble)
species2 <- add_column(species2,"sample"=rownames(species2),.before=1)

species3 <- merge(species2,metadata,by="sample")

mega <- species3[, grep("dm|Megasphaera_elsdenii",colnames(species3))]

#dm and non-dm 
mega1 <- mega[-which(is.na(mega$dm_cl)),]
mega1$preva <- ifelse(is.na(mega1$`d__Bacteria|k__Bacillati|p__Bacillota|c__Negativicutes|o__Veillonellales|f__Veillonellaceae|g__Megasphaera|s__Megasphaera_elsdenii`),0,1)

preva <- data.frame(table(mega1[,c(2,4)]))
   aggregate(preva$Freq,by=list(preva$dm_cl),sum)
#yes 223, no 1153
   dm.pre <- preva$Freq[which(preva$dm_cl=="yes" & preva$preva==1)]/223
   nondm.pre <- preva$Freq[which(preva$dm_cl=="no" & preva$preva==1)]/1153
   
mega2 <- mega1[-which(is.na(mega1$"d__Bacteria|k__Bacillati|p__Bacillota|c__Negativicutes|o__Veillonellales|f__Veillonellaceae|g__Megasphaera|s__Megasphaera_elsdenii")),]
mega2$m.elsdenii <- log10(mega2$`d__Bacteria|k__Bacillati|p__Bacillota|c__Negativicutes|o__Veillonellales|f__Veillonellaceae|g__Megasphaera|s__Megasphaera_elsdenii`)

#dmmed and dmnonmed 
mega3 <- mega[which(mega$dm_cl=="yes" & !is.na(mega$dm_med)),]
mega3$preva <- ifelse(is.na(mega3$`d__Bacteria|k__Bacillati|p__Bacillota|c__Negativicutes|o__Veillonellales|f__Veillonellaceae|g__Megasphaera|s__Megasphaera_elsdenii`),0,1)

preva <- data.frame(table(mega3[,c(3,4)]))
aggregate(preva$Freq,by=list(preva$dm_med),sum)
#yes 126, no 28
dmmed.pre <- preva$Freq[which(preva$dm_med=="yes" & preva$preva==1)]/126
nondmmed.pre <- preva$Freq[which(preva$dm_med=="no" & preva$preva==1)]/28

mega4 <- mega3[-which(is.na(mega3$`d__Bacteria|k__Bacillati|p__Bacillota|c__Negativicutes|o__Veillonellales|f__Veillonellaceae|g__Megasphaera|s__Megasphaera_elsdenii`)),]
mega4$m.elsdenii <- log10(mega4$`d__Bacteria|k__Bacillati|p__Bacillota|c__Negativicutes|o__Veillonellales|f__Veillonellaceae|g__Megasphaera|s__Megasphaera_elsdenii`)


###preva bar　
pdf("M.elsenii.prevalence.pdf",width=3,height=4)
 barplot(height=c(100*nondm.pre,100*dm.pre,100*nondmmed.pre,100*dmmed.pre),names.arg = c("Non-T2D","T2D","NMDG","MDG"),ylim=c(0,100),ylab = "Prevalence (%)",col="#C8C8AC")
 text(1:4,c(100*nondm.pre+3,100*dm.pre+3,100*nondmmed.pre+3,100*dmmed.pre+3),labels=c(round(100*nondm.pre,2),round(100*dm.pre,2),round(100*nondmmed.pre,2),round(100*dmmed.pre,2)))
dev.off()

##combine
mega2$group <- ifelse(mega2$dm_cl=="no","Non-T2D","T2D")
mega4$group <- ifelse(mega4$dm_med=="no","NMDG","MDG")

plot <- rbind(mega2,mega4)
plot$group <- factor(plot$group,levels=c("Non-T2D","T2D","NMDG","MDG"))
library(ggpubr)
ggboxplot(plot,y="m.elsdenii",x="group",fill="group",add = "jitter")+stat_compare_means(comparisons = list(c("Non-T2D","T2D"),c("NMDG","MDG"),c("Non-T2D","NMDG"),c("Non-T2D","MDG")),method = "t.test",aes(label = after_stat(p.signif)))+
  scale_fill_brewer(palette="Set3")+ylab("Log10(relatvie intensity of M. elsdenii)")+theme(legend.position = "none")
ggsave("M.elsdenii.abundance.pdf",width=4,height=4)




#fig7b####
library(openxlsx)
data <- read.xlsx("mega_growth_curve20240126.xlsx",sheet=2)
data1 <- data.frame(t(data))
colnames(data1) <- data1[1,]
data1 <- data1[-1,]
data1$sample <- rownames(data1)

library(reshape2)
data2 <- melt(data1,id.vars = c("group","sample"))
data2$Time <- as.numeric(as.character(data2$variable))

data2$value <- as.numeric(data2$value)
library(tidyr)
library(dplyr)
df_filtered <- data2 %>%
  filter(value >= 0.2 & value <= 0.8)

slopes <- df_filtered %>%
  group_by(group,sample) %>%
  do({
    fit <- lm(value ~ Time, data = .)
    data.frame(Slope = coef(fit)[["Time"]])
  })
slopes$group <- factor(slopes$group,levels=c("control","Acar 20mM","met 20mM","gli 10μM"))


library(ggplot2)
library(ggpubr)
ggplot(slopes, aes(x = group, y = Slope, fill = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +scale_fill_manual(values=c("#d4815d","#99b5d5","#cfb0b9","#a1c8b2"))+stat_compare_means(comparisons = list(c("control","Acar 20mM"),c("control","met 20mM"),c("control","gli 10μM")),method = "t.test",aes(label = after_stat(p.signif)))+
  labs(title = "Growth rate comparison", y = "Slope (growth rate)", x = "Group")+theme_classic()
ggsave("od_slope.pdf")


rm(list=ls())
library(openxlsx)
data <- read.xlsx("mega_growth_curve20240126.xlsx")
data1 <- data

library(reshape2)
data2 <- melt(data,id.vars = c("X1"))
data2$value <- as.numeric(data2$value)
data2$X1 <- factor(data2$X1,levels = c(0,1,3,5,7,8,9))
data2$variable <- factor(data2$variable)

data2$log.cell <- log10(data2$value)

data3 <- aggregate(data2,by=list(data2$variable,data2$X1),mean)
data4 <- aggregate(data2$value,by=list(data2$variable,data2$X1),sd)

colnames(data3) <- c("Group","Time","x1","X2","mean","log10.mean")
colnames(data4) <- c("Group","Time","sd")
data5 <- merge(data3,data4,by=c("Group","Time"))
ggplot(data=data5,aes(x=Time,y=mean,fill=Group))+geom_point(aes(color=Group),size=2)+geom_smooth(aes(group=Group,color=Group),se=F,method = "loess",span=0.6)+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,color=Group),width=0.1)+
  theme_classic()+ylab("Cell number")
ggsave("Mega_growth_cell_number.pdf")



data5=data2[which(data2$X1 != 0),]
data5$id=as.factor(rep(1:6,24))
data4=data5

library(tidyverse)
library(ggpubr)
library(rstatix)


data5=data4[which(data4$variable %in% c("control","acar20")),]
colnames(data5)[1]="group"

acar.res.aov <- anova_test( data = data5, dv = value, wid = id, within = c(variable,group) ) 
get_anova_table(acar.res.aov)


data6=data4[which(data4$variable %in% c("control","met20")),]
colnames(data6)[1]="group"

met.res.aov <- anova_test( data = data6, dv = value, wid = id, within = c(variable,group) ) 
get_anova_table(met.res.aov)

data7=data4[which(data4$variable %in% c("control","gli50")),]
colnames(data7)[1]="group"

gli.res.aov <- anova_test( data = data7, dv = value, wid = id, within = c(variable,group) ) 
get_anova_table(gli.res.aov)




#fig7d####
library(openxlsx)
glu <- read.xlsx("mouse_glucose.xlsx")
library(reshape2)
glu1 <- melt(glu)
glu1$variable <- as.factor(glu1$variable)

glu2 <- aggregate(value~group+variable,data=glu1, mean)
glu3 <- aggregate(value~group+variable,data=glu1, sd)

colnames(glu2) <- c("group","day","mean")
glu2$sd <- glu3$value

library(ggplot2)
ggplot(data=glu2,mapping = aes(x=day,y=mean,color=group))+geom_smooth(aes(group=group,fill=group))+geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(0.05))+
  theme_classic()
ggsave("model_glucose.pdf")

# two way ANOVA with repeated measurements
glu1$id <- as.factor(rep(1:7))
glu1$variable <- as.factor(glu1$variable)
glu1$group=as.factor(glu1$group)

library(tidyverse)
library(ggpubr)
library(rstatix)

res.aov <- anova_test(
  data = glu1, dv = value,wid = id,
  within = c(group, variable)
)
get_anova_table(res.aov)


#fig7e####
scfa <- read.xlsx("scfa_content_20250103.xlsx",sheet=2)
scfa1 <- subset(scfa,group %in% c("CD","HFD","LR","Mega"))
scfa2$group <- as.factor(scfa2$group)
compare <- list(c("CD","Mega"),c("HFD","Mega"),c("Mega","LR"))
colnames(scfa2)[colnames(scfa2) == "umol/g"] <- "umol_g"
ggboxplot(data=scfa2,x="group",y="umol_g",fill="group",group="group")+stat_compare_means(comparisons = compare)
ggsave("butyrate20250106.pdf")
