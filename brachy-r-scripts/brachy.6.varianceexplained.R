rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(reshape2)
library(scales)

setwd("/Users/eludwig/Desktop/brachy/")
DFs_132 <- read.csv("brachy_control_single_data_132_accessions.csv")

DFs_132$Treatment = NA
DFs_132$Treatment[grep("20", DFs_132$Treatment.2)] <- "Drought"
DFs_132$Treatment[grep("100", DFs_132$Treatment.2)] <- "Control"

heatDFs_132 <- read.csv("brachy_heat_single_data_132_accessions.csv")


heatDFs_132$Treatment = NA
heatDFs_132$Treatment[grep("20", heatDFs_132$Treatment.2)] <- "Heat-drought"
heatDFs_132$Treatment[grep("100", heatDFs_132$Treatment.2)] <- "Heat"


bothDF <- merge(DFs_132, heatDFs_132, by = c("Accession", "Treatment", "DAP", "Replicate", "percent_damage", "cor_area", "cor_height_above_reference",
                                     "hue_circular_mean", "hue_circular_std","cor_width", "cor_convex_hull_area", "solidity", "cor_perimeter",
                                     "cor_longest_path", "cor_ellipse_major_axis", "cor_ellipse_minor_axis"), all = TRUE)
bothDF <- dplyr::select(bothDF, "Accession", "Treatment", "DAP", "Replicate", "percent_damage", "cor_area", "cor_height_above_reference",
                        "hue_circular_mean", "hue_circular_std","cor_width", "cor_convex_hull_area", "solidity", "cor_perimeter",
                        "cor_longest_path", "cor_ellipse_major_axis", "cor_ellipse_minor_axis")

bothDF <- filter(bothDF, DAP != 47)
bothDF <- filter(bothDF, DAP != 49)
bothDF <- filter(bothDF, DAP != 46)

bothDF$group = NA
bothDF$group[grep("16", bothDF$DAP)] <- "16"
bothDF$group[grep("18", bothDF$DAP)] <- "18"
bothDF$group[grep("20", bothDF$DAP)] <- "20"
bothDF$group[grep("22", bothDF$DAP)] <- "22"
bothDF$group[grep("24", bothDF$DAP)] <- "24"
bothDF$group[grep("26", bothDF$DAP)] <- "26/27"
bothDF$group[grep("27", bothDF$DAP)] <- "26/27"
bothDF$group[grep("28", bothDF$DAP)] <- "28/29"
bothDF$group[grep("29", bothDF$DAP)] <- "28/29"
bothDF$group[grep("30", bothDF$DAP)] <- "30/31"
bothDF$group[grep("31", bothDF$DAP)] <- "30/31"
bothDF$group[grep("32", bothDF$DAP)] <- "32/33"
bothDF$group[grep("33", bothDF$DAP)] <- "32/33"
bothDF$group[grep("34", bothDF$DAP)] <- "34/35"
bothDF$group[grep("35", bothDF$DAP)] <- "34/35"
bothDF$group[grep("36", bothDF$DAP)] <- "36/37"
bothDF$group[grep("37", bothDF$DAP)] <- "36/37"
bothDF$group[grep("38", bothDF$DAP)] <- "38/39"
bothDF$group[grep("39", bothDF$DAP)] <- "38/39"
bothDF$group[grep("40", bothDF$DAP)] <- "40/41"
bothDF$group[grep("41", bothDF$DAP)] <- "40/41"
bothDF$group[grep("42", bothDF$DAP)] <- "42/43"
bothDF$group[grep("43", bothDF$DAP)] <- "42/43"
bothDF$group[grep("44", bothDF$DAP)] <- "44/45"
bothDF$group[grep("45", bothDF$DAP)] <- "44/45"




########################################
### Time Point 1 (22 DAP) ###
########################################
#keep only data from 22 DAP
tp1 <- filter(bothDF, group == "22")
tp1 <- drop_na(tp1)


#make sure all accessions have data for all treatments
table(tp1$Accession, tp1$Treatment)

#remove accessions with missing data
tp1 <- filter(tp1, Accession != "Bar2" | Accession != "BdTR1H" | Accession != "BdTR3S")

###############################################################
### TP1: Heritability ### control and drought 
###############################################################


ctl_drt_tp1 <- filter(tp1, Treatment == "Control" | Treatment == "Drought")
ctl_drt_tp1 <- dplyr::select(ctl_drt_tp1, -c("DAP", "group"))

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
colnames(ctl_drt_tp1)[5] <- "area"
colnames(ctl_drt_tp1)[6] <- "height_above_reference"
colnames(ctl_drt_tp1)[9] <- "width"
colnames(ctl_drt_tp1)[10] <- "convex_hull_area"
colnames(ctl_drt_tp1)[12] <- "perimeter"
colnames(ctl_drt_tp1)[13] <- "longest_path"
colnames(ctl_drt_tp1)[14] <- "ellipse_major_axis"
colnames(ctl_drt_tp1)[15] <- "ellipse_minor_axis"

tail(ctl_drt_tp1)
shapes <- c("percent_damage", "area", "height_above_reference",
            "hue_circular_mean", "hue_circular_std", "width", "convex_hull_area", "solidity", "perimeter",
            "longest_path", "ellipse_major_axis", "ellipse_minor_axis")
Htp1_drt <- c()
for(e in shapes){
  modeltp1_drt <- lmer(eval(parse(text=e))~(1|Accession)+(1|Treatment)+(1|Accession:Treatment),data = ctl_drt_tp1)
  retp1_drt<-as.numeric(VarCorr(modeltp1_drt))
  restp1_drt<-attr(VarCorr(modeltp1_drt), "sc")^2
  interaction.vartp1_drt <- retp1_drt[1]
  Accession.vartp1_drt<-retp1_drt[2]
  treatment.vartp1_drt<-retp1_drt[3]
  tot.vartp1_drt<-sum(retp1_drt,restp1_drt)
  unexptp1_drt <- 1-sum(retp1_drt)/sum(retp1_drt,restp1_drt)
  htp1_drt <- c((Accession.vartp1_drt/tot.vartp1_drt),
          (treatment.vartp1_drt/tot.vartp1_drt),
          (interaction.vartp1_drt/tot.vartp1_drt),
          unexptp1_drt)
  Htp1_drt <- rbind(Htp1_drt,htp1_drt)
}
Htp1_drt <- data.frame(Htp1_drt,row.names = shapes)
Htp1_drt$Shape <- rownames(Htp1_drt)
rownames(Htp1_drt) <- NULL
colnames(Htp1_drt) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
Htp1_drt$Shape <-  ordered(Htp1_drt$Shape,levels=Htp1_drt$Shape[order(Htp1_drt$Unexplained)])
write.csv(Htp1_drt, "/Users/eludwig/Downloads/tp1_drought_variance_explained.csv")
Htp1_drt_melt <- melt(Htp1_drt,id=c("Shape"))
Htp1_drt_melt$value <- Htp1_drt_melt$value*100
Htp1_drt_melt$variable <- ordered(Htp1_drt_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(Htp1_drt_melt)

#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#004D40", "#0072B2", "#F0E442", "#D55E00", "#636985")

TP1_drt <- ggplot(data=Htp1_drt_melt,aes(Shape,value))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Drought Time Point 1")+
  theme(plot.title = element_text(size = 12))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill = NA, linewidth = 1,linetype = 1))+
  theme(legend.position = "right")+
  theme(legend.text=element_text(size=14))+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
TP1_drt

ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_drought_TP1.png",width=5.5,height=5.5,plot = TP1_drt, dpi = 300)
ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_drought_TP1.pdf",width=5.5,height=5.5,plot = TP1_drt, dpi = 300)


###############################################################
### TP1: Heritability ### control and heat
###############################################################
ctl_heat_tp1 <- filter(tp1, Treatment == "Control" | Treatment == "Heat")
ctl_heat_tp1 <- dplyr::select(ctl_heat_tp1,-c("DAP", "group"))

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
colnames(ctl_heat_tp1)[5] <- "area"
colnames(ctl_heat_tp1)[6] <- "height_above_reference"
colnames(ctl_heat_tp1)[9] <- "width"
colnames(ctl_heat_tp1)[10] <- "convex_hull_area"
colnames(ctl_heat_tp1)[12] <- "perimeter"
colnames(ctl_heat_tp1)[13] <- "longest_path"
colnames(ctl_heat_tp1)[14] <- "ellipse_major_axis"
colnames(ctl_heat_tp1)[15] <- "ellipse_minor_axis"

tail(ctl_heat_tp1)
shapes <- c("percent_damage", "area", "height_above_reference",
            "hue_circular_mean", "hue_circular_std", "width", "convex_hull_area", "solidity", "perimeter",
            "longest_path", "ellipse_major_axis", "ellipse_minor_axis")
Htp1_heat <- c()
for(e in shapes){
  modeltp1_heat <- lmer(eval(parse(text=e))~(1|Accession)+(1|Treatment)+(1|Accession:Treatment),data = ctl_heat_tp1)
  retp1_heat<-as.numeric(VarCorr(modeltp1_heat))
  restp1_heat<-attr(VarCorr(modeltp1_heat), "sc")^2
  interaction.vartp1_heat <- retp1_heat[1]
  Accession.vartp1_heat<-retp1_heat[2]
  treatment.vartp1_heat<-retp1_heat[3]
  tot.vartp1_heat<-sum(retp1_heat,restp1_heat)
  unexptp1_heat <- 1-sum(retp1_heat)/sum(retp1_heat,restp1_heat)
  htp1_heat <- c((Accession.vartp1_heat/tot.vartp1_heat),
          (treatment.vartp1_heat/tot.vartp1_heat),
          (interaction.vartp1_heat/tot.vartp1_heat),
          unexptp1_heat)
  Htp1_heat <- rbind(Htp1_heat,htp1_heat)
}
Htp1_heat <- data.frame(Htp1_heat,row.names = shapes)
Htp1_heat$Shape <- rownames(Htp1_heat)
rownames(Htp1_heat) <- NULL
colnames(Htp1_heat) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
Htp1_heat$Shape <-  ordered(Htp1_heat$Shape,levels=Htp1_heat$Shape[order(Htp1_heat$Unexplained)])
write.csv(Htp1_heat, "/Users/eludwig/Downloads/tp1_heat_variance_explained.csv")
Htp1_heat_melt <- melt(Htp1_heat,id=c("Shape"))
Htp1_heat_melt$value <- Htp1_heat_melt$value*100
Htp1_heat_melt$variable <- ordered(Htp1_heat_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(Htp1_heat_melt)

TP1_heat <- ggplot(data=Htp1_heat_melt,aes(Shape,value))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Heat Time Point 1")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill = NA, linewidth = 1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
TP1_heat

ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heat_TP1.png",width=5.5,height=5.5,plot = TP1_heat, dpi = 300)
ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heat_TP1.pdf",width=5.5,height=5.5,plot = TP1_heat, dpi = 300)



###############################################################
### TP1: Heritability ### control and heat-drought
###############################################################
ctl_heatdrt_tp1 <- filter(tp1, Treatment == "Control" | Treatment == "Heat-drought")
ctl_heatdrt_tp1 <- dplyr::select(ctl_heatdrt_tp1, -c("DAP", "group"))

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
colnames(ctl_heatdrt_tp1)[5] <- "area"
colnames(ctl_heatdrt_tp1)[6] <- "height_above_reference"
colnames(ctl_heatdrt_tp1)[9] <- "width"
colnames(ctl_heatdrt_tp1)[10] <- "convex_hull_area"
colnames(ctl_heatdrt_tp1)[12] <- "perimeter"
colnames(ctl_heatdrt_tp1)[13] <- "longest_path"
colnames(ctl_heatdrt_tp1)[14] <- "ellipse_major_axis"
colnames(ctl_heatdrt_tp1)[15] <- "ellipse_minor_axis"

tail(ctl_heatdrt_tp1)
shapes <- c("percent_damage", "area", "height_above_reference",
            "hue_circular_mean", "hue_circular_std", "width", "convex_hull_area", "solidity", "perimeter",
            "longest_path", "ellipse_major_axis", "ellipse_minor_axis")

Htp1_heatdrt <- c()
for(e in shapes){
  modeltp1_heatdrt <- lmer(eval(parse(text=e))~(1|Accession)+(1|Treatment)+(1|Accession:Treatment),data = ctl_heatdrt_tp1)
  retp1_heatdrt<-as.numeric(VarCorr(modeltp1_heatdrt))
  restp1_heatdrt<-attr(VarCorr(modeltp1_heatdrt), "sc")^2
  interaction.vartp1_heatdrt <- retp1_heatdrt[1]
  Accession.vartp1_heatdrt<-retp1_heatdrt[2]
  treatment.vartp1_heatdrt<-retp1_heatdrt[3]
  tot.vartp1_heatdrt<-sum(retp1_heatdrt,restp1_heatdrt)
  unexptp1_heatdrt <- 1-sum(retp1_heatdrt)/sum(retp1_heatdrt,restp1_heatdrt)
  htp1_heatdrt <- c((Accession.vartp1_heatdrt/tot.vartp1_heatdrt),
          (treatment.vartp1_heatdrt/tot.vartp1_heatdrt),
          (interaction.vartp1_heatdrt/tot.vartp1_heatdrt),
          unexptp1_heatdrt)
  Htp1_heatdrt <- rbind(Htp1_heatdrt,htp1_heatdrt)
}
Htp1_heatdrt <- data.frame(Htp1_heatdrt,row.names = shapes)
Htp1_heatdrt$Shape <- rownames(Htp1_heatdrt)
rownames(Htp1_heatdrt) <- NULL
colnames(Htp1_heatdrt) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
Htp1_heatdrt$Shape <-  ordered(Htp1_heatdrt$Shape,levels=Htp1_heatdrt$Shape[order(Htp1_heatdrt$Unexplained)])
write.csv(Htp1_heatdrt, "/Users/eludwig/Downloads/tp1_heatdrought_variance_explained.csv")
Htp1_heatdrt_melt <- melt(Htp1_heatdrt,id=c("Shape"))
Htp1_heatdrt_melt$value <- Htp1_heatdrt_melt$value*100
Htp1_heatdrt_melt$variable <- ordered(Htp1_heatdrt_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(Htp1_heatdrt_melt)

TP1_heatdrt <- ggplot(data=Htp1_heatdrt_melt,aes(Shape,value))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Heat-Drought Time Point 1")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill = NA, linewidth = 1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
TP1_heatdrt

ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heatdrought_TP1.png",width=5.5,height=5.5,plot = TP1_heatdrt, dpi = 300)
ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heatdrought_TP1.pdf",width=5.5,height=5.5,plot = TP1_heatdrt, dpi = 300)





########################################
### Time Point 2 (28/29 DAP) ###
########################################
#keep only data from 22 DAP
tp2 <- filter(bothDF, group == "28/29")
tp2 <- drop_na(tp2)

#make sure all accessions have data for all treatments
table(tp2$Accession, tp2$Treatment)

#remove accessions with missing data
tp2 <- filter(tp2, Accession != "Bar2" | Accession != "BdTR1H" | Accession != "BdTR3S")


###############################################################
### TP2: Heritability ### control and drought 
###############################################################
ctl_drt_tp2 <- filter(tp2, Treatment == "Control" | Treatment == "Drought")
ctl_drt_tp2 <- dplyr::select(ctl_drt_tp2, -c("DAP", "group"))

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
colnames(ctl_drt_tp2)[5] <- "area"
colnames(ctl_drt_tp2)[6] <- "height_above_reference"
colnames(ctl_drt_tp2)[9] <- "width"
colnames(ctl_drt_tp2)[10] <- "convex_hull_area"
colnames(ctl_drt_tp2)[12] <- "perimeter"
colnames(ctl_drt_tp2)[13] <- "longest_path"
colnames(ctl_drt_tp2)[14] <- "ellipse_major_axis"
colnames(ctl_drt_tp2)[15] <- "ellipse_minor_axis"

tail(ctl_drt_tp2)
shapes <- c("percent_damage", "area", "height_above_reference",
            "hue_circular_mean", "hue_circular_std", "width", "convex_hull_area", "solidity", "perimeter",
            "longest_path", "ellipse_major_axis", "ellipse_minor_axis")
Htp2_drt <- c()
for(e in shapes){
  modeltp2_drt <- lmer(eval(parse(text=e))~(1|Accession)+(1|Treatment)+(1|Accession:Treatment),data = ctl_drt_tp2)
  retp2_drt<-as.numeric(VarCorr(modeltp2_drt))
  restp2_drt<-attr(VarCorr(modeltp2_drt), "sc")^2
  interaction.vartp2_drt <- retp2_drt[1]
  Accession.vartp2_drt<-retp2_drt[2]
  treatment.vartp2_drt<-retp2_drt[3]
  tot.vartp2_drt<-sum(retp2_drt,restp2_drt)
  unexptp2_drt <- 1-sum(retp2_drt)/sum(retp2_drt,restp2_drt)
  htp2_drt <- c((Accession.vartp2_drt/tot.vartp2_drt),
                (treatment.vartp2_drt/tot.vartp2_drt),
                (interaction.vartp2_drt/tot.vartp2_drt),
                unexptp2_drt)
  Htp2_drt <- rbind(Htp2_drt,htp2_drt)
}
Htp2_drt <- data.frame(Htp2_drt,row.names = shapes)
Htp2_drt$Shape <- rownames(Htp2_drt)
rownames(Htp2_drt) <- NULL
colnames(Htp2_drt) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
Htp2_drt$Shape <-  ordered(Htp2_drt$Shape,levels=Htp2_drt$Shape[order(Htp2_drt$Unexplained)])
write.csv(Htp2_drt, "/Users/eludwig/Downloads/tp2_drought_variance_explained.csv")
Htp2_drt_melt <- melt(Htp2_drt,id=c("Shape"))
Htp2_drt_melt$value <- Htp2_drt_melt$value*100
Htp2_drt_melt$variable <- ordered(Htp2_drt_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(Htp2_drt_melt)

#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#004D40", "#0072B2", "#F0E442", "#D55E00", "#636985")

TP2_drt <- ggplot(data=Htp2_drt_melt,aes(Shape,value))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Drought Time Point 2")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill = NA, linewidth = 1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
TP2_drt

ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_drought_TP2.png",width=5.5,height=5.5,plot = TP2_drt, dpi = 300)
ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_drought_TP2.pdf",width=5.5,height=5.5,plot = TP2_drt, dpi = 300)


###############################################################
### TP2: Heritability ### control and heat
###############################################################
ctl_heat_tp2 <- filter(tp2, Treatment == "Control" | Treatment == "Heat")
ctl_heat_tp2 <- dplyr::select(ctl_heat_tp2, -c("DAP", "group"))

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
colnames(ctl_heat_tp2)[5] <- "area"
colnames(ctl_heat_tp2)[6] <- "height_above_reference"
colnames(ctl_heat_tp2)[9] <- "width"
colnames(ctl_heat_tp2)[10] <- "convex_hull_area"
colnames(ctl_heat_tp2)[12] <- "perimeter"
colnames(ctl_heat_tp2)[13] <- "longest_path"
colnames(ctl_heat_tp2)[14] <- "ellipse_major_axis"
colnames(ctl_heat_tp2)[15] <- "ellipse_minor_axis"

tail(ctl_heat_tp2)
shapes <- c("percent_damage", "area", "height_above_reference",
            "hue_circular_mean", "hue_circular_std", "width", "convex_hull_area", "solidity", "perimeter",
            "longest_path", "ellipse_major_axis", "ellipse_minor_axis")
Htp2_heat <- c()
for(e in shapes){
  modeltp2_heat <- lmer(eval(parse(text=e))~(1|Accession)+(1|Treatment)+(1|Accession:Treatment),data = ctl_heat_tp2)
  retp2_heat<-as.numeric(VarCorr(modeltp2_heat))
  restp2_heat<-attr(VarCorr(modeltp2_heat), "sc")^2
  interaction.vartp2_heat <- retp2_heat[1]
  Accession.vartp2_heat<-retp2_heat[2]
  treatment.vartp2_heat<-retp2_heat[3]
  tot.vartp2_heat<-sum(retp2_heat,restp2_heat)
  unexptp2_heat <- 1-sum(retp2_heat)/sum(retp2_heat,restp2_heat)
  htp2_heat <- c((Accession.vartp2_heat/tot.vartp2_heat),
                 (treatment.vartp2_heat/tot.vartp2_heat),
                 (interaction.vartp2_heat/tot.vartp2_heat),
                 unexptp2_heat)
  Htp2_heat <- rbind(Htp2_heat,htp2_heat)
}
Htp2_heat <- data.frame(Htp2_heat,row.names = shapes)
Htp2_heat$Shape <- rownames(Htp2_heat)
rownames(Htp2_heat) <- NULL
colnames(Htp2_heat) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
Htp2_heat$Shape <-  ordered(Htp2_heat$Shape,levels=Htp2_heat$Shape[order(Htp2_heat$Unexplained)])
write.csv(Htp2_heat, "/Users/eludwig/Downloads/tp2_heat_variance_explained.csv")
Htp2_heat_melt <- melt(Htp2_heat,id=c("Shape"))
Htp2_heat_melt$value <- Htp2_heat_melt$value*100
Htp2_heat_melt$variable <- ordered(Htp2_heat_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(Htp2_heat_melt)

TP2_heat <- ggplot(data=Htp2_heat_melt,aes(Shape,value))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Heat Time Point 2")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill = NA, linewidth = 1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
TP2_heat

ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heat_TP2.png",width=5.5,height=5.5,plot = TP2_heat, dpi = 300)
ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heat_TP2.pdf",width=5.5,height=5.5,plot = TP2_heat, dpi = 300)



###############################################################
### TP2: Heritability ### control and heat-drought
###############################################################
ctl_heatdrt_tp2 <- filter(tp2, Treatment == "Control" | Treatment == "Heat-drought")
ctl_heatdrt_tp2 <- dplyr::select(ctl_heatdrt_tp2, -c("DAP", "group"))

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
colnames(ctl_heatdrt_tp2)[5] <- "area"
colnames(ctl_heatdrt_tp2)[6] <- "height_above_reference"
colnames(ctl_heatdrt_tp2)[9] <- "width"
colnames(ctl_heatdrt_tp2)[10] <- "convex_hull_area"
colnames(ctl_heatdrt_tp2)[12] <- "perimeter"
colnames(ctl_heatdrt_tp2)[13] <- "longest_path"
colnames(ctl_heatdrt_tp2)[14] <- "ellipse_major_axis"
colnames(ctl_heatdrt_tp2)[15] <- "ellipse_minor_axis"

tail(ctl_heatdrt_tp2)
shapes <- c("percent_damage", "area", "height_above_reference",
            "hue_circular_mean", "hue_circular_std", "width", "convex_hull_area", "solidity", "perimeter",
            "longest_path", "ellipse_major_axis", "ellipse_minor_axis")

Htp2_heatdrt <- c()
for(e in shapes){
  modeltp2_heatdrt <- lmer(eval(parse(text=e))~(1|Accession)+(1|Treatment)+(1|Accession:Treatment),data = ctl_heatdrt_tp2)
  retp2_heatdrt<-as.numeric(VarCorr(modeltp2_heatdrt))
  restp2_heatdrt<-attr(VarCorr(modeltp2_heatdrt), "sc")^2
  interaction.vartp2_heatdrt <- retp2_heatdrt[1]
  Accession.vartp2_heatdrt<-retp2_heatdrt[2]
  treatment.vartp2_heatdrt<-retp2_heatdrt[3]
  tot.vartp2_heatdrt<-sum(retp2_heatdrt,restp2_heatdrt)
  unexptp2_heatdrt <- 1-sum(retp2_heatdrt)/sum(retp2_heatdrt,restp2_heatdrt)
  htp2_heatdrt <- c((Accession.vartp2_heatdrt/tot.vartp2_heatdrt),
                    (treatment.vartp2_heatdrt/tot.vartp2_heatdrt),
                    (interaction.vartp2_heatdrt/tot.vartp2_heatdrt),
                    unexptp2_heatdrt)
  Htp2_heatdrt <- rbind(Htp2_heatdrt,htp2_heatdrt)
}
Htp2_heatdrt <- data.frame(Htp2_heatdrt,row.names = shapes)
Htp2_heatdrt$Shape <- rownames(Htp2_heatdrt)
rownames(Htp2_heatdrt) <- NULL
colnames(Htp2_heatdrt) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
Htp2_heatdrt$Shape <-  ordered(Htp2_heatdrt$Shape,levels=Htp2_heatdrt$Shape[order(Htp2_heatdrt$Unexplained)])
write.csv(Htp2_heatdrt, "/Users/eludwig/Downloads/tp2_heatdrought_variance_explained.csv")
Htp2_heatdrt_melt <- melt(Htp2_heatdrt,id=c("Shape"))
Htp2_heatdrt_melt$value <- Htp2_heatdrt_melt$value*100
Htp2_heatdrt_melt$variable <- ordered(Htp2_heatdrt_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(Htp2_heatdrt_melt)

TP2_heatdrt <- ggplot(data=Htp2_heatdrt_melt,aes(Shape,value))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Heat-Drought Time Point 2")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill = NA, linewidth = 1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
TP2_heatdrt

ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heatdrought_TP2.png",width=5.5,height=5.5,plot = TP2_heatdrt, dpi = 300)
ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heatdrought_TP2.pdf",width=5.5,height=5.5,plot = TP2_heatdrt, dpi = 300)


########################################
### Time Point 3 (36/37 DAP) ###
########################################
#keep only data from 36/37 DAP
tp3 <- filter(bothDF, group == "36/37")
tp3 <- drop_na(tp3)

#make sure all accessions have data for all treatments
table(tp3$Accession, tp3$Treatment)

#remove accessions with missing data
tp3 <- filter(tp3, Accession != "Bar2" | Accession != "BdTR1H" | Accession != "BdTR3S")

#tp3_real <- tp3
###############################################################
### TP3: Heritability ### control and drought 
###############################################################
ctl_drt_tp3 <- filter(tp3, Treatment == "Control" | Treatment == "Drought")
ctl_drt_tp3 <- dplyr::select(ctl_drt_tp3, -c("DAP", "group"))

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
colnames(ctl_drt_tp3)[5] <- "area"
colnames(ctl_drt_tp3)[6] <- "height_above_reference"
colnames(ctl_drt_tp3)[9] <- "width"
colnames(ctl_drt_tp3)[10] <- "convex_hull_area"
colnames(ctl_drt_tp3)[12] <- "perimeter"
colnames(ctl_drt_tp3)[13] <- "longest_path"
colnames(ctl_drt_tp3)[14] <- "ellipse_major_axis"
colnames(ctl_drt_tp3)[15] <- "ellipse_minor_axis"

tail(ctl_drt_tp3)
shapes <- c("percent_damage", "area", "height_above_reference",
            "hue_circular_mean", "hue_circular_std", "width", "convex_hull_area", "solidity", "perimeter",
            "longest_path", "ellipse_major_axis", "ellipse_minor_axis")
Htp3_drt <- c()
e="longest_path"
for(e in shapes){
  modeltp3_drt <- lmer(eval(parse(text=e))~(1|Accession)+(1|Treatment)+(1|Accession:Treatment),data = ctl_drt_tp3)
  retp3_drt<-as.numeric(VarCorr(modeltp3_drt))
  restp3_drt<-attr(VarCorr(modeltp3_drt), "sc")^2
  interaction.vartp3_drt <- retp3_drt[1]
  Accession.vartp3_drt<-retp3_drt[2]
  treatment.vartp3_drt<-retp3_drt[3]
  tot.vartp3_drt<-sum(retp3_drt,restp3_drt)
  unexptp3_drt <- 1-sum(retp3_drt)/sum(retp3_drt,restp3_drt)
  htp3_drt <- c((Accession.vartp3_drt/tot.vartp3_drt),
                (treatment.vartp3_drt/tot.vartp3_drt),
                (interaction.vartp3_drt/tot.vartp3_drt),
                unexptp3_drt)
  Htp3_drt <- rbind(Htp3_drt,htp3_drt)
}
Htp3_drt <- data.frame(Htp3_drt,row.names = shapes)
Htp3_drt$Shape <- rownames(Htp3_drt)
rownames(Htp3_drt) <- NULL
colnames(Htp3_drt) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
Htp3_drt$Shape <-  ordered(Htp3_drt$Shape,levels=Htp3_drt$Shape[order(Htp3_drt$Unexplained)])
write.csv(Htp3_drt, "/Users/eludwig/Downloads/tp3_drought_variance_explained.csv")
Htp3_drt_melt <- melt(Htp3_drt,id=c("Shape"))
Htp3_drt_melt$value <- Htp3_drt_melt$value*100
Htp3_drt_melt$variable <- ordered(Htp3_drt_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(Htp3_drt_melt)

#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#004D40", "#0072B2", "#F0E442", "#D55E00", "#636985")

TP3_drt <- ggplot(data=Htp3_drt_melt,aes(Shape,value))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Drought Time Point 3")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill = NA, linewidth = 1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
TP3_drt

ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_drought_TP3.png",width=5.5,height=5.5,plot = TP3_drt, dpi = 300)
ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_drought_TP3.pdf",width=5.5,height=5.5,plot = TP3_drt, dpi = 300)


###############################################################
### TP3: Heritability ### control and heat
###############################################################
ctl_heat_tp3 <- filter(tp3, Treatment == "Control" | Treatment == "Heat")
ctl_heat_tp3 <- dplyr::select(ctl_heat_tp3, -c("DAP", "group"))

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
colnames(ctl_heat_tp3)[5] <- "area"
colnames(ctl_heat_tp3)[6] <- "height_above_reference"
colnames(ctl_heat_tp3)[9] <- "width"
colnames(ctl_heat_tp3)[10] <- "convex_hull_area"
colnames(ctl_heat_tp3)[12] <- "perimeter"
colnames(ctl_heat_tp3)[13] <- "longest_path"
colnames(ctl_heat_tp3)[14] <- "ellipse_major_axis"
colnames(ctl_heat_tp3)[15] <- "ellipse_minor_axis"

tail(ctl_heat_tp3)
shapes <- c("percent_damage", "area", "height_above_reference",
            "hue_circular_mean", "hue_circular_std", "width", "convex_hull_area", "solidity", "perimeter",
            "longest_path", "ellipse_major_axis", "ellipse_minor_axis")
Htp3_heat <- c()
for(e in shapes){
  modeltp3_heat <- lmer(eval(parse(text=e))~(1|Accession)+(1|Treatment)+(1|Accession:Treatment),data = ctl_heat_tp3)
  retp3_heat<-as.numeric(VarCorr(modeltp3_heat))
  restp3_heat<-attr(VarCorr(modeltp3_heat), "sc")^2
  interaction.vartp3_heat <- retp3_heat[1]
  Accession.vartp3_heat<-retp3_heat[2]
  treatment.vartp3_heat<-retp3_heat[3]
  tot.vartp3_heat<-sum(retp3_heat,restp3_heat)
  unexptp3_heat <- 1-sum(retp3_heat)/sum(retp3_heat,restp3_heat)
  htp3_heat <- c((Accession.vartp3_heat/tot.vartp3_heat),
                 (treatment.vartp3_heat/tot.vartp3_heat),
                 (interaction.vartp3_heat/tot.vartp3_heat),
                 unexptp3_heat)
  Htp3_heat <- rbind(Htp3_heat,htp3_heat)
}
Htp3_heat <- data.frame(Htp3_heat,row.names = shapes)
Htp3_heat$Shape <- rownames(Htp3_heat)
rownames(Htp3_heat) <- NULL
colnames(Htp3_heat) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
Htp3_heat$Shape <-  ordered(Htp3_heat$Shape,levels=Htp3_heat$Shape[order(Htp3_heat$Unexplained)])
write.csv(Htp3_heat, "/Users/eludwig/Downloads/tp3_heat_variance_explained.csv")
Htp3_heat_melt <- melt(Htp3_heat,id=c("Shape"))
Htp3_heat_melt$value <- Htp3_heat_melt$value*100
Htp3_heat_melt$variable <- ordered(Htp3_heat_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(Htp3_heat_melt)

TP3_heat <- ggplot(data=Htp3_heat_melt,aes(Shape,value))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Heat Time Point 3")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill = NA, linewidth = 1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
TP3_heat

ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heat_TP3.png",width=5.5,height=5.5,plot = TP3_heat, dpi = 300)
ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heat_TP3.pdf",width=5.5,height=5.5,plot = TP3_heat, dpi = 300)



###############################################################
### TP3: Heritability ### control and heat-drought
###############################################################
ctl_heatdrt_tp3 <- filter(tp3, Treatment == "Control" | Treatment == "Heat-drought")
ctl_heatdrt_tp3 <- dplyr::select(ctl_heatdrt_tp3, -c("DAP", "group"))

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
colnames(ctl_heatdrt_tp3)[5] <- "area"
colnames(ctl_heatdrt_tp3)[6] <- "height_above_reference"
colnames(ctl_heatdrt_tp3)[9] <- "width"
colnames(ctl_heatdrt_tp3)[10] <- "convex_hull_area"
colnames(ctl_heatdrt_tp3)[12] <- "perimeter"
colnames(ctl_heatdrt_tp3)[13] <- "longest_path"
colnames(ctl_heatdrt_tp3)[14] <- "ellipse_major_axis"
colnames(ctl_heatdrt_tp3)[15] <- "ellipse_minor_axis"

tail(ctl_heatdrt_tp3)
shapes <- c("percent_damage", "area", "height_above_reference",
            "hue_circular_mean", "hue_circular_std", "width", "convex_hull_area", "solidity", "perimeter",
            "longest_path", "ellipse_major_axis", "ellipse_minor_axis")

Htp3_heatdrt <- c()
for(e in shapes){
  modeltp3_heatdrt <- lmer(eval(parse(text=e))~(1|Accession)+(1|Treatment)+(1|Accession:Treatment),data = ctl_heatdrt_tp3)
  retp3_heatdrt<-as.numeric(VarCorr(modeltp3_heatdrt))
  restp3_heatdrt<-attr(VarCorr(modeltp3_heatdrt), "sc")^2
  interaction.vartp3_heatdrt <- retp3_heatdrt[1]
  Accession.vartp3_heatdrt<-retp3_heatdrt[2]
  treatment.vartp3_heatdrt<-retp3_heatdrt[3]
  tot.vartp3_heatdrt<-sum(retp3_heatdrt,restp3_heatdrt)
  unexptp3_heatdrt <- 1-sum(retp3_heatdrt)/sum(retp3_heatdrt,restp3_heatdrt)
  htp3_heatdrt <- c((Accession.vartp3_heatdrt/tot.vartp3_heatdrt),
                    (treatment.vartp3_heatdrt/tot.vartp3_heatdrt),
                    (interaction.vartp3_heatdrt/tot.vartp3_heatdrt),
                    unexptp3_heatdrt)
  Htp3_heatdrt <- rbind(Htp3_heatdrt,htp3_heatdrt)
}
Htp3_heatdrt <- data.frame(Htp3_heatdrt,row.names = shapes)
Htp3_heatdrt$Shape <- rownames(Htp3_heatdrt)
rownames(Htp3_heatdrt) <- NULL
colnames(Htp3_heatdrt) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
Htp3_heatdrt$Shape <-  ordered(Htp3_heatdrt$Shape,levels=Htp3_heatdrt$Shape[order(Htp3_heatdrt$Unexplained)])
write.csv(Htp3_heatdrt, "/Users/eludwig/Downloads/tp3_heatdrought_variance_explained.csv")
Htp3_heatdrt_melt <- melt(Htp3_heatdrt,id=c("Shape"))
Htp3_heatdrt_melt$value <- Htp3_heatdrt_melt$value*100
Htp3_heatdrt_melt$variable <- ordered(Htp3_heatdrt_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(Htp3_heatdrt_melt)

TP3_heatdrt <- ggplot(data=Htp3_heatdrt_melt,aes(Shape,value))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Heat-Drought Time Point 3")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill = NA, linewidth = 1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
TP3_heatdrt

ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heatdrought_TP3.png",width=5.5,height=5.5,plot = TP3_heatdrt, dpi = 300)
ggsave("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heatdrought_TP3.pdf",width=5.5,height=5.5,plot = TP3_heatdrt, dpi = 300)


###############################################################
### Time Point 4 (44/45 DAP) ###
###############################################################
#keep only data from 44/45 DAP
tp4 <- filter(bothDF, group == "44/45")
tp4 <- drop_na(tp4)

#make sure all accessions have data for all treatments
table(tp4$Accession, tp4$Treatment)

#remove accessions with missing data
tp4 <- filter(tp4, Accession != "Bar2" | Accession != "BdTR1H" | Accession != "BdTR3S")

###############################################################
### TP4: Heritability ### control and drought 
###############################################################
ctl_drt_tp4 <- filter(tp4, Treatment == "Control" | Treatment == "Drought")
ctl_drt_tp4 <- dplyr::select(ctl_drt_tp4, -c("DAP", "group"))

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
colnames(ctl_drt_tp4)[5] <- "area"
colnames(ctl_drt_tp4)[6] <- "height_above_reference"
colnames(ctl_drt_tp4)[9] <- "width"
colnames(ctl_drt_tp4)[10] <- "convex_hull_area"
colnames(ctl_drt_tp4)[12] <- "perimeter"
colnames(ctl_drt_tp4)[13] <- "longest_path"
colnames(ctl_drt_tp4)[14] <- "ellipse_major_axis"
colnames(ctl_drt_tp4)[15] <- "ellipse_minor_axis"

tail(ctl_drt_tp4)
shapes <- c("percent_damage", "area", "height_above_reference",
            "hue_circular_mean", "hue_circular_std", "width", "convex_hull_area", "solidity", "perimeter",
            "longest_path", "ellipse_major_axis", "ellipse_minor_axis")
Htp4_drt <- c()
for(e in shapes){
  modeltp4_drt <- lmer(eval(parse(text=e))~(1|Accession)+(1|Treatment)+(1|Accession:Treatment),data = ctl_drt_tp4)
  retp4_drt<-as.numeric(VarCorr(modeltp4_drt))
  restp4_drt<-attr(VarCorr(modeltp4_drt), "sc")^2
  interaction.vartp4_drt <- retp4_drt[1]
  Accession.vartp4_drt<-retp4_drt[2]
  treatment.vartp4_drt<-retp4_drt[3]
  tot.vartp4_drt<-sum(retp4_drt,restp4_drt)
  unexptp4_drt <- 1-sum(retp4_drt)/sum(retp4_drt,restp4_drt)
  htp4_drt <- c((Accession.vartp4_drt/tot.vartp4_drt),
          (treatment.vartp4_drt/tot.vartp4_drt),
          (interaction.vartp4_drt/tot.vartp4_drt),
          unexptp4_drt)
  Htp4_drt <- rbind(Htp4_drt,htp4_drt)
}
Htp4_drt <- data.frame(Htp4_drt,row.names = shapes)
Htp4_drt$Shape <- rownames(Htp4_drt)
rownames(Htp4_drt) <- NULL
colnames(Htp4_drt) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
Htp4_drt$Shape <-  ordered(Htp4_drt$Shape,levels=Htp4_drt$Shape[order(Htp4_drt$Unexplained)])
write.csv(Htp4_drt, "/Users/eludwig/Downloads/tp4_drought_variance_explained.csv")
Htp4_drt_melt <- melt(Htp4_drt,id=c("Shape"))
Htp4_drt_melt$value <- Htp4_drt_melt$value*100
Htp4_drt_melt$variable <- ordered(Htp4_drt_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(Htp4_drt_melt)

#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#004D40", "#0072B2", "#F0E442", "#D55E00", "#636985")


TP4_drt <- ggplot(data=Htp4_drt_melt,aes(Shape,value))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Drought Time Point 4")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill = NA, linewidth = 1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
TP4_drt

ggsave(plot = TP4_drt,"/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_drought_tp4.png",width=5.5,height=5.5, dpi = 300)
ggsave(plot = TP4_drt,"/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_drought_tp4.pdf",width=5.5,height=5.5, dpi = 300)


###############################################################
### TP4: Heritability ### control and heat
###############################################################
ctl_heat_tp4 <- filter(tp4, Treatment == "Control" | Treatment == "Heat")
ctl_heat_tp4 <- dplyr::select(ctl_heat_tp4, -c("DAP", "group"))

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
colnames(ctl_heat_tp4)[5] <- "area"
colnames(ctl_heat_tp4)[6] <- "height_above_reference"
colnames(ctl_heat_tp4)[9] <- "width"
colnames(ctl_heat_tp4)[10] <- "convex_hull_area"
colnames(ctl_heat_tp4)[12] <- "perimeter"
colnames(ctl_heat_tp4)[13] <- "longest_path"
colnames(ctl_heat_tp4)[14] <- "ellipse_major_axis"
colnames(ctl_heat_tp4)[15] <- "ellipse_minor_axis"

tail(ctl_heat_tp4)
shapes <- c("percent_damage", "area", "height_above_reference",
            "hue_circular_mean", "hue_circular_std", "width", "convex_hull_area", "solidity", "perimeter",
            "longest_path", "ellipse_major_axis", "ellipse_minor_axis")
Htp4_heat <- c()
for(e in shapes){
  modeltp4_heatdrt <- lmer(eval(parse(text=e))~(1|Accession)+(1|Treatment)+(1|Accession:Treatment),data = ctl_heat_tp4)
  retp4_heatdrt<-as.numeric(VarCorr(modeltp4_heatdrt))
  restp4_heatdrt<-attr(VarCorr(modeltp4_heatdrt), "sc")^2
  interaction.vartp4_heatdrt <- retp4_heatdrt[1]
  Accession.vartp4_heatdrt<-retp4_heatdrt[2]
  treatment.vartp4_heatdrt<-retp4_heatdrt[3]
  tot.vartp4_heatdrt<-sum(retp4_heatdrt,restp4_heatdrt)
  unexptp4_heatdrt <- 1-sum(retp4_heatdrt)/sum(retp4_heatdrt,restp4_heatdrt)
  htp4_heat <- c((Accession.vartp4_heatdrt/tot.vartp4_heatdrt),
          (treatment.vartp4_heatdrt/tot.vartp4_heatdrt),
          (interaction.vartp4_heatdrt/tot.vartp4_heatdrt),
          unexptp4_heatdrt)
  Htp4_heat <- rbind(Htp4_heat,htp4_heat)
}
Htp4_heat <- data.frame(Htp4_heat,row.names = shapes)
Htp4_heat$Shape <- rownames(Htp4_heat)
rownames(Htp4_heat) <- NULL
colnames(Htp4_heat) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
Htp4_heat$Shape <-  ordered(Htp4_heat$Shape,levels=Htp4_heat$Shape[order(Htp4_heat$Unexplained)])
write.csv(Htp4_heat, "/Users/eludwig/Downloads/tp4_heat_variance_explained.csv")
Htp4_heat_melt <- melt(Htp4_heat,id=c("Shape"))
Htp4_heat_melt$value <- Htp4_heat_melt$value*100
Htp4_heat_melt$variable <- ordered(Htp4_heat_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(Htp4_heat_melt)

TP4_heat <- ggplot(data=Htp4_heat_melt,aes(Shape,value))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Heat Time Point 4")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill = NA, linewidth = 1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
TP4_heat

ggsave(plot = TP4_heat, "/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heat_tp4.png",width=5.5,height=5.5, dpi = 300)
ggsave(plot = TP4_heat, "/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heat_tp4.pdf",width=5.5,height=5.5, dpi = 300)



###############################################################
### TP4: Heritability ### control and heat-drought
###############################################################
ctl_heatdrt_tp4 <- filter(tp4, Treatment == "Control" | Treatment == "Heat-drought")
ctl_heatdrt_tp4 <- dplyr::select(ctl_heatdrt_tp4, -c("DAP", "group"))

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
colnames(ctl_heatdrt_tp4)[5] <- "area"
colnames(ctl_heatdrt_tp4)[6] <- "height_above_reference"
colnames(ctl_heatdrt_tp4)[9] <- "width"
colnames(ctl_heatdrt_tp4)[10] <- "convex_hull_area"
colnames(ctl_heatdrt_tp4)[12] <- "perimeter"
colnames(ctl_heatdrt_tp4)[13] <- "longest_path"
colnames(ctl_heatdrt_tp4)[14] <- "ellipse_major_axis"
colnames(ctl_heatdrt_tp4)[15] <- "ellipse_minor_axis"

tail(ctl_heatdrt_tp4)
shapes <- c("percent_damage", "area", "height_above_reference",
            "hue_circular_mean", "hue_circular_std", "width", "convex_hull_area", "solidity", "perimeter",
            "longest_path", "ellipse_major_axis", "ellipse_minor_axis")

Htp4_heatdrt <- c()
for(e in shapes){
  modeltp4_heatdrt <- lmer(eval(parse(text=e))~(1|Accession)+(1|Treatment)+(1|Accession:Treatment),data = ctl_heatdrt_tp4)
  retp4_heatdrt<-as.numeric(VarCorr(modeltp4_heatdrt))
  restp4_heatdrt<-attr(VarCorr(modeltp4_heatdrt), "sc")^2
  interaction.vartp4_heatdrt <- retp4_heatdrt[1]
  Accession.vartp4_heatdrt<-retp4_heatdrt[2]
  treatment.vartp4_heatdrt<-retp4_heatdrt[3]
  tot.vartp4_heatdrt<-sum(retp4_heatdrt,restp4_heatdrt)
  unexptp4_heatdrt <- 1-sum(retp4_heatdrt)/sum(retp4_heatdrt,restp4_heatdrt)
  htp4_heatdrt <- c((Accession.vartp4_heatdrt/tot.vartp4_heatdrt),
          (treatment.vartp4_heatdrt/tot.vartp4_heatdrt),
          (interaction.vartp4_heatdrt/tot.vartp4_heatdrt),
          unexptp4_heatdrt)
  Htp4_heatdrt <- rbind(Htp4_heatdrt,htp4_heatdrt)
}
Htp4_heatdrt <- data.frame(Htp4_heatdrt,row.names = shapes)
Htp4_heatdrt$Shape <- rownames(Htp4_heatdrt)
rownames(Htp4_heatdrt) <- NULL
colnames(Htp4_heatdrt) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
Htp4_heatdrt$Shape <-  ordered(Htp4_heatdrt$Shape,levels=Htp4_heatdrt$Shape[order(Htp4_heatdrt$Unexplained)])
write.csv(Htp4_heatdrt, "/Users/eludwig/Downloads/tp4_heatdrought_variance_explained.csv")
Htp4_heatdrt_melt <- melt(Htp4_heatdrt,id=c("Shape"))
Htp4_heatdrt_melt$value <- Htp4_heatdrt_melt$value*100
Htp4_heatdrt_melt$variable <- ordered(Htp4_heatdrt_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(Htp4_heatdrt_melt)

TP4_heatdrt <- ggplot(data=Htp4_heatdrt_melt,aes(Shape,value))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Heat-Drought Time Point 4")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill = NA, linewidth = 1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
TP4_heatdrt

ggsave(plot = TP4_heatdrt, "/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heatdrought_tp4.png",width=5.5,height=5.5, dpi = 300)
ggsave(plot = TP4_heatdrt, "/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained_heatdrought_tp4.pdf",width=5.5,height=5.5, dpi = 300)




###############################################################
### Combine plots ###
###############################################################

# need to set legend.position="right" (or anything other than none) for Day 1, then change for others
legend <- get_legend(TP1_drt)
#pheno_pca_1 <- pheno_pca_1 + theme(legend.position = "none")


### All plots organized by treatment

# combine PCA plots from all days and add legend
var_expl_combined <- ggarrange(TP1_drt, TP2_drt, TP3_drt, TP4_drt,
                               TP1_heat, TP2_heat, TP3_heat, TP4_heat,
                               TP1_heatdrt, TP2_heatdrt, TP3_heatdrt, TP4_heatdrt, 
                                ncol = 4, nrow = 3, common.legend = TRUE, legend = "top") +
  theme(plot.margin = unit(c(0.5, 0, 0, 0), "cm")) 
# Add overall title to combined PCAs
var_expl_combined <- annotate_figure(var_expl_combined,
                                     top = text_grob("Variance Explained between Stress Treatments and Control", color = "black", face = "bold", size = 28),
                                     bottom = text_grob("Trait", color = "black", size = 24),
                                     left = text_grob("Variance Explained (%)", color = "black", rot = 90, size = 24)) +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
var_expl_combined

# save combined PCA plot as pdf and png
#ggsave(plot = pheno_pca_combined,"/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/phenotype_pca.png", width = 7.5, height = 10, dpi = 300)
ggsave(plot = var_expl_combined,"/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/variance_explained.pdf", width = 20, height = 15, dpi = 300)
ggsave(plot = var_expl_combined,"/Users/eludwig/Downloads/variance_explained.pdf", width = 20, height = 15, dpi = 300)





