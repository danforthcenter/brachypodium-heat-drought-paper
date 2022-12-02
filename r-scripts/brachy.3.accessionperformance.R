rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
#library(vioplot)
library(FactoMineR)
library(factoextra)
library(missMDA)
library(ggpubr)
library(plotly)

setwd("/Users/eludwig/Desktop/brachy/")
DFs <- read.csv("brachy_control_VIS_SV_single_data_compiled_currated.csv")

DFs$Treatment = NA
DFs$Treatment[grep("20", DFs$Treatment.2)] <- "Drought"
DFs$Treatment[grep("100", DFs$Treatment.2)] <- "Control"

control <- filter(DFs, Treatment == "Control", DAP == 44)
drought <- filter(DFs, Treatment == "Drought", DAP == 44)

heatDFs <- read.csv("brachy_heat_VIS_SV_single_data_compiled_currated.csv")

heatDFs$Treatment = NA
heatDFs$Treatment[grep("20", heatDFs$Treatment.2)] <- "Heat-drought"
heatDFs$Treatment[grep("100", heatDFs$Treatment.2)] <- "Heat"

heat <- filter(heatDFs, Treatment == "Heat", DAP == 45)
heatdrought <- filter(heatDFs, Treatment == "Heat-drought", DAP == 45)

###########################################
### Phenotype PCA ###
###########################################


bothDF <- merge(DFs, heatDFs, by = c("Accession", "Treatment", "DAP", "Replicate", "percent_damage", "cor_area", "cor_height_above_reference",
                                     "hue_circular_mean", "cor_width", "cor_height", "cor_convex_hull_area", "solidity", "cor_perimeter",
                                     "cor_longest_path", "cor_ellipse_major_axis", "cor_ellipse_minor_axis", "ellipse_angle", "ellipse_eccentricity"), all = TRUE)
bothDF <- dplyr::select(bothDF, "Accession", "Treatment", "DAP", "Replicate", "percent_damage", "cor_area", "cor_height_above_reference",
                 "hue_circular_mean", "cor_width", "cor_height", "cor_convex_hull_area", "solidity", "cor_perimeter",
                 "cor_longest_path", "cor_ellipse_major_axis", "cor_ellipse_minor_axis", "ellipse_angle", "ellipse_eccentricity")

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


avg_damage <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_percent_damage = mean(percent_damage, na.rm = TRUE))
avg_area <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_area = mean(cor_area))
avg_height <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_height = mean(cor_height_above_reference))
avg_hue <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_hue = mean(hue_circular_mean))
avg_width <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_width = mean(cor_width))
avg_height2 <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_height2 = mean(cor_height))
avg_convex <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_convex_hull_area = mean(cor_convex_hull_area))
avg_solidity <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_solidity = mean(solidity))
avg_perimeter <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_perimeter = mean(cor_perimeter))
avg_longest_path <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_longest_path = mean(cor_longest_path))
avg_ellipse_major_axis <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_ellipse_major_axis = mean(cor_ellipse_major_axis))
avg_ellipse_minor_axis <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_ellipse_minor_axis = mean(cor_ellipse_minor_axis))
avg_ellipse_angle <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_ellipse_angle = mean(ellipse_angle))
avg_ellipse_eccentricity <- group_by(bothDF, Accession, Treatment, DAP, group) %>% summarise(avg_ellipse_eccentricity = mean(ellipse_eccentricity))



# combine all the phenotype data frames together for PCA
pheno_pca_data <- merge(avg_damage, avg_area, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)
pheno_pca_data <- merge(pheno_pca_data, avg_height, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)
pheno_pca_data <- merge(pheno_pca_data, avg_hue, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)
pheno_pca_data <- merge(pheno_pca_data, avg_width, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)
pheno_pca_data <- merge(pheno_pca_data, avg_height2, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)
pheno_pca_data <- merge(pheno_pca_data, avg_convex, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)
pheno_pca_data <- merge(pheno_pca_data, avg_solidity, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)
pheno_pca_data <- merge(pheno_pca_data, avg_perimeter, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)
pheno_pca_data <- merge(pheno_pca_data, avg_longest_path, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)
pheno_pca_data <- merge(pheno_pca_data, avg_ellipse_major_axis, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)
pheno_pca_data <- merge(pheno_pca_data, avg_ellipse_minor_axis, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)
pheno_pca_data <- merge(pheno_pca_data, avg_ellipse_angle, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)
pheno_pca_data <- merge(pheno_pca_data, avg_ellipse_eccentricity, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)


# Do PCA with phenotype data
#pheno_pca_data_imputed <- as.data.frame(imputePCA(pheno_pca_data[,4:17], ncp = 5)[[1]]) #impute missing values
pheno.pca = PCA(pheno_pca_data[,5:18], scale.unit=TRUE, ncp=5, graph=F)

# make data frame with Accession and treatment info and PC values
pheno_pca_df = data.frame(
  "Accession" = pheno_pca_data$Accession,
  "Treatment" = pheno_pca_data$Treatment,
  "DAP" = pheno_pca_data$DAP,
  "group" = pheno_pca_data$group,
  "PC1" = pheno.pca$ind$coord[, 1],
  "PC2" = pheno.pca$ind$coord[, 2],
  "PC3" = pheno.pca$ind$coord[, 3],
  "PC4" = pheno.pca$ind$coord[, 4])

# find the percent variance explained by each PC
percentofvariance = data.frame(pheno.pca$eig) 
percentofvariance$PC <- row.names(percentofvariance)
row.names(percentofvariance) <- NULL
percentofvariance$PC <- row.names(percentofvariance)
percentofvariance$PC <- as.numeric(percentofvariance$PC)
percentofvariance$percentage.of.variance <- as.numeric(percentofvariance$percentage.of.variance)
str(percentofvariance)

pheno.loadings <- sweep(pheno.pca$var$coord,2,sqrt(pheno.pca$eig[1:ncol(pheno.pca$var$coord),1]),FUN="/")

pheno_pca_df <- merge(pheno_pca_df, pheno_pca_data, by = c("Accession", "Treatment", "DAP", "group"), all = FALSE)


# make colorblind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#004D40", "#0072B2", "#F0E442", "#D55E00", "#636985")

pheno_pca_plot <- pheno_pca_df%>%
  mutate(group = factor(paste0(group, " DAP")))%>%
  ggplot(aes(x = PC1, y = PC2, color = Treatment))+
  facet_wrap(~group, nrow=5)+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_point(size=0.5)+
  theme_minimal()+
  theme(panel.grid = element_blank(), panel.border = element_rect(fill="transparent"))+
  scale_colour_manual(values=cbPalette)+
  stat_ellipse()+
  xlab("PC 1 (59.90%)")+
  ylab("PC 2 (13.41%)")+
  theme(legend.position=c(.50,-0.10))+
  guides(color=guide_legend(nrow=1))+
  ggtitle("PCA of Phenotypic Measurements over Time") +
  theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm"))
pheno_pca_plot
ggsave(plot = pheno_pca_plot,"/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/phenotype_pca.png", width = 7.5, height = 10, dpi = 300)
ggsave(plot = pheno_pca_plot,"/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/phenotype_pca.pdf", width = 7.5, height = 10, dpi = 300)


###########################################
### Setting up comparison between control and stresses for overlap ###
###########################################


#damage_lines <- filter(avg_damage, Accession == "BdTR10A" | Accession == "BdTR9K" | Accession == "BdTR11A" | Accession == "ABR5")
damage_lines15 <- filter(avg_damage, group == "44/45")
damage_lines15 <- ungroup(damage_lines15)
damage_lines15 <- dplyr::select(damage_lines15, -c(DAP, group))
damage_lines_dif <- spread(damage_lines15, key = Treatment, value = avg_percent_damage)
damage_lines_dif <- drop_na(damage_lines_dif)
damage_lines_dif$drought_dif = NA
damage_lines_dif$heat_dif = NA
damage_lines_dif$heatdrought_dif = NA
damage_lines_dif$drought_dif <- (damage_lines_dif$Control-damage_lines_dif$Drought)
damage_lines_dif$heat_dif <- (damage_lines_dif$Control-damage_lines_dif$Heat)
damage_lines_dif$heatdrought_dif <- (damage_lines_dif$Control-damage_lines_dif$`Heat-drought`)
damage_lines_dif <- dplyr::select(damage_lines_dif, Accession, drought_dif, heat_dif, heatdrought_dif)
damage_lines_dif <- damage_lines_dif %>% rename(Drought = drought_dif,Heat = heat_dif, `Heat-drought`= heatdrought_dif)
damage_lines_dif <- gather(damage_lines_dif, Treatment, damage_dif, Drought:`Heat-drought`)
damage_lines_dif <- join(damage_lines_dif, damage_lines15, by= c("Accession", "Treatment"))

damage_drought <- filter(damage_lines_dif, Treatment == "Drought")
damage_heat <- filter(damage_lines_dif, Treatment == "Heat")
damage_heatdrought <- filter(damage_lines_dif, Treatment == "Heat-drought")

damage_plot<-ggplot(damage_lines_dif, aes(x=Treatment, y=damage_dif, group=Accession)) +
  geom_line(aes(color=Accession))+
  geom_point(aes(color=Accession,size=(avg_percent_damage)))
ggplotly(damage_plot)

damage_plot_drought<-ggplot(damage_drought, aes(x=Accession, y=damage_dif)) +
  geom_point(aes(color=Accession,size=(avg_percent_damage)))+
  ylim(-90, 25)+
  ggtitle("Drought Compared to Control")
ggplotly(damage_plot_drought)

damage_plot_heat<-ggplot(damage_heat, aes(x=Accession, y=damage_dif)) +
  geom_point(aes(color=Accession,size=(avg_percent_damage)))+
  ylim(-90, 25)+
  ggtitle("Heat Compared to Control")
ggplotly(damage_plot_heat)

damage_plot_heatdrought<-ggplot(damage_heatdrought, aes(x=Accession, y=damage_dif)) +
  geom_point(aes(color=Accession,size=(avg_percent_damage)))+
  ylim(-90, 25)+
  ggtitle("Heat-drought Compared to Control")
ggplotly(damage_plot_heatdrought)


# combine PCA plots from all days and add legend
damage_combined <- ggarrange(damage_plot_drought, damage_plot_heat, damage_plot_heatdrought, legend = "none",
                                ncol = 3, nrow = 1)
damage_combined




#area_lines <- filter(avg_area, Accession == "BdTR10A" | Accession == "BdTR9K" | Accession == "BdTR11A" | Accession == "ABR5")
area_lines15 <- filter(avg_area, group == "44/45")
area_lines15 <- ungroup(area_lines15)
area_lines15 <- dplyr::select(area_lines15, -c(DAP, group))
area_lines_dif <- spread(area_lines15, key = Treatment, value = avg_area)
area_lines_dif <- drop_na(area_lines_dif)
area_lines_dif$drought_dif = NA
area_lines_dif$heat_dif = NA
area_lines_dif$heatdrought_dif = NA
area_lines_dif$drought_dif <- (area_lines_dif$Control-area_lines_dif$Drought)
area_lines_dif$heat_dif <- (area_lines_dif$Control-area_lines_dif$Heat)
area_lines_dif$heatdrought_dif <- (area_lines_dif$Control-area_lines_dif$`Heat-drought`)
area_lines_dif <- dplyr::select(area_lines_dif, Accession, drought_dif, heat_dif, heatdrought_dif)
area_lines_dif <- area_lines_dif %>% rename(Drought = drought_dif,Heat = heat_dif, `Heat-drought`= heatdrought_dif)
area_lines_dif <- gather(area_lines_dif, Treatment, area_dif, Drought:`Heat-drought`)
area_lines_dif <- join(area_lines_dif, area_lines15, by= c("Accession", "Treatment"))

area_drought <- filter(area_lines_dif, Treatment == "Drought")
area_heat <- filter(area_lines_dif, Treatment == "Heat")
area_heatdrought <- filter(area_lines_dif, Treatment == "Heat-drought")



area_plot<-ggplot(area_lines_dif, aes(x=Treatment, y=area_dif, group=Accession)) +
  geom_line(aes(color=Accession))+
  geom_point(aes(color=Accession))+
  scale_y_reverse()
ggplotly(area_plot)


area_plot_drought<-ggplot(area_drought, aes(x=Accession, y=area_dif)) +
  geom_point(aes(color=Accession,size=(avg_area)))+
  ylim(55, 0)+ #range 8 to 39
  ggtitle("Drought Compared to Control")
ggplotly(area_plot_drought)

area_plot_heat<-ggplot(area_heat, aes(x=Accession, y=area_dif)) +
  geom_point(aes(color=Accession,size=(avg_area)))+
  ylim(55, 0)+ #range 14 to 47
  ggtitle("Heat Compared to Control")
ggplotly(area_plot_heat)

area_plot_heatdrought<-ggplot(area_heatdrought, aes(x=Accession, y=area_dif)) +
  geom_point(aes(color=Accession,size=(avg_area)))+
  ylim(55, 0)+ #range 15 to 49
  ggtitle("Heat-drought Compared to Control")
ggplotly(area_plot_heatdrought)


# combine PCA plots from all days and add legend
area_combined <- ggarrange(area_plot_drought, area_plot_heat, area_plot_heatdrought, legend = "none",
                                ncol = 3, nrow = 1)
area_combined





#height_lines <- filter(avg_height, Accession == "BdTR10A" | Accession == "BdTR9K" | Accession == "BdTR11A" | Accession == "ABR5")
height_lines15 <- filter(avg_height, group == "44/45")
height_lines15 <- ungroup(height_lines15)
height_lines15 <- dplyr::select(height_lines15, -c(DAP, group))
height_lines_dif <- spread(height_lines15, key = Treatment, value = avg_height)
height_lines_dif <- drop_na(height_lines_dif)
height_lines_dif$drought_dif = NA
height_lines_dif$heat_dif = NA
height_lines_dif$heatdrought_dif = NA
height_lines_dif$drought_dif <- (height_lines_dif$Control-height_lines_dif$Drought)
height_lines_dif$heat_dif <- (height_lines_dif$Control-height_lines_dif$Heat)
height_lines_dif$heatdrought_dif <- (height_lines_dif$Control-height_lines_dif$"Heat-drought")
height_lines_dif <- dplyr::select(height_lines_dif, Accession, drought_dif, heat_dif, heatdrought_dif)
height_lines_dif <- height_lines_dif %>% rename(Drought = drought_dif,Heat = heat_dif, `Heat-drought`= heatdrought_dif)
height_lines_dif <- gather(height_lines_dif, Treatment, height_dif, Drought:`Heat-drought`)
height_lines_dif <- join(height_lines_dif, height_lines15, by= c("Accession", "Treatment"))

height_drought <- filter(height_lines_dif, Treatment == "Drought")
height_heat <- filter(height_lines_dif, Treatment == "Heat")
height_heatdrought <- filter(height_lines_dif, Treatment == "Heat-drought")


height_plot<-ggplot(height_lines_dif, aes(x=Treatment, y=height_dif, group=Accession)) +
  geom_line(aes(color=Accession))+
  geom_point(aes(color=Accession))+
  scale_y_reverse()
ggplotly(height_plot)


height_plot_drought<-ggplot(height_drought, aes(x=Accession, y=height_dif)) +
  geom_point(aes(color=Accession,size=(avg_height)))+
  ylim(20, -5)+ #range -.6 to 9
  ggtitle("Drought Compared to Control")
ggplotly(height_plot_drought)

height_plot_heat<-ggplot(height_heat, aes(x=Accession, y=height_dif)) +
  geom_point(aes(color=Accession,size=(avg_height)))+
  ylim(20, -5)+ #range 2 to 14
  ggtitle("Heat Compared to Control")
ggplotly(height_plot_heat)

height_plot_heatdrought<-ggplot(height_heatdrought, aes(x=Accession, y=height_dif)) +
  geom_point(aes(color=Accession,size=(avg_height)))+
  ylim(20, -5)+ #range 4 to 14
  ggtitle("Heat-drought Compared to Control")
ggplotly(height_plot_heatdrought)


# combine PCA plots from all days and add legend
height_combined <- ggarrange(height_plot_drought, height_plot_heat, height_plot_heatdrought, legend = "none",
                           ncol = 3, nrow = 1)
height_combined


#calculate means of columns and make df
means_of_dif <- data.frame(matrix(ncol = 9, nrow = 1))
x <- c("height_drought", "height_heat", "height_heatdrought", "area_drought", "area_heat", "area_heatdrought",
       "damage_drought", "damage_heat", "damage_heatdrought")
colnames(means_of_dif) <- x
means_of_dif$height_drought <- mean(height_drought$height_dif)
means_of_dif$height_heat <- mean(height_heat$height_dif)
means_of_dif$height_heatdrought <- mean(height_heatdrought$height_dif)
means_of_dif$area_drought <- mean(area_drought$area_dif)
means_of_dif$area_heat <- mean(area_heat$area_dif)
means_of_dif$area_heatdrought <- mean(area_heatdrought$area_dif)
means_of_dif$damage_drought <- mean(damage_drought$damage_dif)
means_of_dif$damage_heat <- mean(damage_heat$damage_dif)
means_of_dif$damage_heatdrought <- mean(damage_heatdrought$damage_dif)


#calculate means of columns and make df
means <- data.frame(matrix(ncol = 9, nrow = 1))
x <- c("height_drought", "height_heat", "height_heatdrought", "area_drought", "area_heat", "area_heatdrought",
       "damage_drought", "damage_heat", "damage_heatdrought")
colnames(means) <- x
means$height_drought <- mean(height_drought$avg_height)
means$height_heat <- mean(height_heat$avg_height)
means$height_heatdrought <- mean(height_heatdrought$avg_height)
means$area_drought <- mean(area_drought$avg_area)
means$area_heat <- mean(area_heat$avg_area)
means$area_heatdrought <- mean(area_heatdrought$avg_area)
means$damage_drought <- mean(damage_drought$avg_percent_damage)
means$damage_heat <- mean(damage_heat$avg_percent_damage)
means$damage_heatdrought <- mean(damage_heatdrought$avg_percent_damage)


qqnorm(height_drought$avg_height)
qqline(height_drought$avg_height)

qqnorm(height_heat$avg_height)
qqline(height_heat$avg_height)

qqnorm(height_heatdrought$avg_height)
qqline(height_heatdrought$avg_height)


qqnorm(area_drought$avg_area)
qqline(area_drought$avg_area)

qqnorm(area_heat$avg_area)
qqline(area_heat$avg_area)

qqnorm(area_heatdrought$avg_area)
qqline(area_heatdrought$avg_area)


qqnorm(damage_drought$avg_percent_damage)
qqline(damage_drought$avg_percent_damage)

qqnorm(damage_heat$avg_percent_damage)
qqline(damage_heat$avg_percent_damage)

qqnorm(damage_heatdrought$avg_percent_damage)
qqline(damage_heatdrought$avg_percent_damage)



t.test(height_heat$avg_height, height_heatdrought$avg_height, alternative = "two.sided", var.equal = TRUE)
t.test(area_heat$avg_area, area_heatdrought$avg_area, alternative = "two.sided", var.equal = TRUE)
t.test(damage_drought$avg_percent_damage, damage_heatdrought$avg_percent_damage, alternative = "two.sided", var.equal = TRUE)

t.test(height_heat$height_dif, height_heatdrought$height_dif, alternative = "two.sided", var.equal = TRUE)
t.test(area_heat$area_dif, area_heatdrought$area_dif, alternative = "two.sided", var.equal = TRUE)
t.test(damage_drought$damage_dif, damage_heatdrought$damage_dif, alternative = "two.sided", var.equal = TRUE)


###########################################
### messing around and testing plots ###
###########################################
test <- filter(DFs, DAP > 10)
test <- filter(test, Treatment == "WL")




# Notched Boxplot of Tooth Growth Against 2 Crossed Factors
# boxes colored for ease of interpretation 
box <- boxplot(cor_area~Accession, data=DFs, notch=TRUE)
        #main="Tooth Growth", xlab="Suppliment and Dose")
box

# Notched Boxplot of Tooth Growth Against 2 Crossed Factors
# boxes colored for ease of interpretation 
box <- boxplot(cor_area~Accession, data=DFs, notch=TRUE)
#main="Tooth Growth", xlab="Suppliment and Dose")
box

ggsave(data = box, "/Users/ellaludwig/Desktop/geno_pattern/boxplot.pdf", width = 7.25, height = 6, dpi = 300)




vioplot(test$cor_area ~ test$Accession)

x1 <- test$cor_area[test$Accession=="Abi1"]
x2 <- test$cor_area[test$Accession=="Adi-3"]
x3 <- test$cor_area[test$Accession=="BdTR11C"]
x4 <- test$cor_area[test$Accession=="BdTR1G"]
x5 <- test$cor_area[test$Accession=="BdTR11A"]
x6 <- test$cor_area[test$Accession=="BdTR12B"]
x7 <- test$cor_area[test$Accession=="BdTR1C"]
x8 <- test$cor_area[test$Accession=="BdTR10J"]
x9 <- test$cor_area[test$Accession=="BdTR13M"]
x10 <- test$cor_area[test$Accession=="BdTR12A"]
x11 <- test$cor_area[test$Accession=="BdTR5G"]
x12 <- test$cor_area[test$Accession=="ABR5"]
x13 <- test$cor_area[test$Accession=="Mon3"]
x14 <- test$cor_area[test$Accession=="Gaz-2"]
x15 <- test$cor_area[test$Accession=="Bar2"]


pdf(file="/Users/ellaludwig/Desktop/geno_pattern/drought_boxplot_area.pdf",width = 8,height = 12)

boxplot(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, notch = TRUE,
        names=c("Abi1", "Adi-3", "BdTR11C", "BdTR1G", "BdTR11A", "BdTR12B", "BdTR1C", "BdTR10J", "BdTR13M", "BdTR12A", "BdTR5G", "ABR5", "Mon3", "Gaz-2", "Bar2"), 
        col="grey")
dev.off()




###########################################
# percentage of damaged plant tissue setup #
###########################################

# control and drought #
avg_damage <- ungroup(avg_damage)

ctl_damage <- filter(avg_damage, Treatment == "Control") #make a DF with all images in the control treatment
drt_damage <- filter (avg_damage, Treatment == "Drought") #make a DF with all images in the drought treamtment

ctl_damage <- filter(ctl_damage, group == "44/45")
ctl_damage <- dplyr::select(ctl_damage, Accession, avg_percent_damage)
ctl_damage <- ctl_damage %>% arrange(avg_percent_damage)
write.csv(ctl_damage, "/Users/eludwig/Google Drive/My Drive/brachy_project//brachy-figures/percent_damage_control.csv")

drt_damage <- filter(drt_damage, group == "44/45")
drt_damage <- dplyr::select(drt_damage, Accession, avg_percent_damage)
drt_damage <- drt_damage %>% arrange(avg_percent_damage)
write.csv(drt_damage, "/Users/eludwig/Google Drive/My Drive/brachy_project//brachy-figures/percent_damage_drought.csv")


# heat and heat-drought #
heat_damage <- filter(avg_damage, Treatment == "Heat") #make a DF with all images in the control treatment
heatdrt_damage <- filter (avg_damage, Treatment == "Heat-drought") #make a DF with all images in the drought treamtment

heat_damage <- filter(heat_damage, group == "44/45")
heat_damage <- dplyr::select(heat_damage, Accession, avg_percent_damage)
heat_damage <- heat_damage %>% arrange(avg_percent_damage)
write.csv(heat_damage, "/Users/eludwig/Google Drive/My Drive/brachy_project//brachy-figures/percent_damage_heat.csv")


heatdrt_damage <- filter(heatdrt_damage, group == "44/45")
heatdrt_damage <- dplyr::select(heatdrt_damage, Accession, avg_percent_damage)
heatdrt_damage <- heatdrt_damage %>% arrange(avg_percent_damage)
write.csv(heatdrt_damage, "/Users/eludwig/Google Drive/My Drive/brachy_project//brachy-figures/percent_damage_heatdrought.csv")




###########################################
# biomass/area setup#
###########################################

# control and drought #

#calculate average area and create avg df
avg_area <- ungroup(avg_area)

ctl_area <- filter(avg_area, Treatment == "Control") #make a DF with all images in the control treatment
drt_area <- filter(avg_area, Treatment == "Drought") #make a DF with all images in the control treatment

ctl_area <- filter(ctl_area, group == "44/45")
ctl_area <- dplyr::select(ctl_area, Accession, avg_area)
ctl_area <- ctl_area %>% arrange(desc(avg_area))
write.csv(ctl_area, "/Users/eludwig/Google Drive/My Drive/brachy_project//brachy-figures/area_control.csv")

drt_area <- filter(drt_area, group == "44/45")
drt_area <- dplyr::select(drt_area, Accession, avg_area)
drt_area <- drt_area %>% arrange(desc(avg_area))
write.csv(drt_area, "/Users/eludwig/Google Drive/My Drive/brachy_project//brachy-figures/area_drought.csv")

# heat and heat-drought #
heat_area <- filter(avg_area, Treatment == "Heat") #make a DF with all images in the control treatment
heatdrt_area <- filter (avg_area, Treatment == "Heat-drought") #make a DF with all images in the drought treamtment

heat_area <- filter(heat_area, group == "44/45")
heat_area <- dplyr::select(heat_area, Accession, avg_area)
heat_area <- heat_area %>% arrange(desc(avg_area))
write.csv(heat_area, "/Users/eludwig/Google Drive/My Drive/brachy_project//brachy-figures/area_heat.csv")

heatdrt_area <- filter(heatdrt_area, group == "44/45")
heatdrt_area <- dplyr::select(heatdrt_area, Accession, avg_area)
heatdrt_area <- heatdrt_area %>% arrange(desc(avg_area))
write.csv(heatdrt_area, "/Users/eludwig/Google Drive/My Drive/brachy_project//brachy-figures/area_heatdrought.csv")




###########################################
# height setup #
###########################################

# control and drought #

#calculate average area and create avg df
avg_height <- ungroup(avg_height)

ctl_height <- filter(avg_height, Treatment == "Control") #make a DF with all images in the control treatment
drt_height <- filter(avg_height, Treatment == "Drought") #make a DF with all images in the control treatment

ctl_height <- filter(ctl_height, group == "44/45")
ctl_height <- dplyr::select(ctl_height, Accession, avg_height)
ctl_height <- ctl_height %>% arrange(desc(avg_height))
write.csv(ctl_height, "/Users/eludwig/Google Drive/My Drive/brachy_project//brachy-figures/height_control.csv")

drt_height <- filter(drt_height, group == "44/45")
drt_height <- dplyr::select(drt_height, Accession, avg_height)
drt_height <- drt_height %>% arrange(desc(avg_height))
write.csv(drt_height, "/Users/eludwig/Google Drive/My Drive/brachy_project//brachy-figures/height_drought.csv")

# heat and heat-drought #
heat_height <- filter(avg_height, Treatment == "Heat") #make a DF with all images in the control treatment
heatdrt_height <- filter (avg_height, Treatment == "Heat-drought") #make a DF with all images in the drought treamtment

heat_height <- filter(heat_height, group == "44/45")
heat_height <- dplyr::select(heat_height, Accession, avg_height)
heat_height <- heat_height %>% arrange(desc(avg_height))
write.csv(heat_height, "/Users/eludwig/Google Drive/My Drive/brachy_project//brachy-figures/height_heat.csv")

heatdrt_height <- filter(heatdrt_height, group == "44/45")
heatdrt_height <- dplyr::select(heatdrt_height, Accession, avg_height)
heatdrt_height <- heatdrt_height %>% arrange(desc(avg_height))
write.csv(heatdrt_height, "/Users/eludwig/Google Drive/My Drive/brachy_project//brachy-figures/height_heatdrought.csv")


#####
#**********************************************
### Lists of good and bad performing lines ###
#**********************************************
#####
############################################
### % damage Comparisons ###
############################################

### Percent damage ###

# look at overlap of difference between stresses and control

#damage_drought (129)
#damage_heat (129)
#damage_heatdrought (129)

damage_drought <- as.data.frame(damage_drought)
damage_heat <- as.data.frame(damage_heat)
damage_heatdrought <- as.data.frame(damage_heatdrought)

damage_drought <- drop_na(dplyr::select(damage_drought, Accession, damage_dif))
damage_drought$damage_dif <- as.numeric(damage_drought$damage_dif)
top_drt_damage <- damage_drought %>% top_n(32) 
bottom_drt_damage <- damage_drought %>% top_n(-32) 

damage_heat <- drop_na(dplyr::select(damage_heat, Accession, damage_dif))
top_heat_damage <- damage_heat %>% top_n(32) 
bottom_heat_damage <- damage_heat %>% top_n(-32) 

damage_heatdrought <- drop_na(dplyr::select(damage_heatdrought, Accession, damage_dif))
top_heatdrt_damage <- damage_heatdrought %>% top_n(32) 
bottom_heatdrt_damage <- damage_heatdrought %>% top_n(-32) 

###see which accessions overlap that perform well across stress treatments
# compare heat to heat-drought
list1 <- unique(top_heatdrt_damage$Accession)
top_heat_damage$compare.geno <- top_heat_damage$Accession %in% list1
all_overlap_stressed_top <- filter(top_heat_damage, compare.geno == TRUE)
length(unique(all_overlap_stressed_top$Accession))
all_overlap_stressed_top = dplyr::select(all_overlap_stressed_top, -compare.geno)
list2 <- unique(all_overlap_stressed_top$Accession)

#compare heat and drought
list1 <- unique(top_heat_damage$Accession)
top_drt_damage$compare.geno <- top_drt_damage$Accession %in% list1
all_overlap_stressed_top <- filter(top_drt_damage, compare.geno == TRUE)
length(unique(all_overlap_stressed_top$Accession))
all_overlap_stressed_top = dplyr::select(all_overlap_stressed_top, -compare.geno)
list3 <- unique(all_overlap_stressed_top$Accession)

#compare heat-drought and drought
list1 <- unique(top_heatdrt_damage$Accession)
top_drt_damage$compare.geno <- top_drt_damage$Accession %in% list1
all_overlap_stressed_top <- filter(top_drt_damage, compare.geno == TRUE)
length(unique(all_overlap_stressed_top$Accession))
all_overlap_stressed_top = dplyr::select(all_overlap_stressed_top, -compare.geno)
list4 <- unique(all_overlap_stressed_top$Accession)

#compare heat/heat-drought to drought to find overall overlap
top_drt_damage$compare.geno <- top_drt_damage$Accession %in% list2
all_overlap_stressed_top <- filter(top_drt_damage, compare.geno == TRUE)




### See which accessions overlap that perform poorly ###

# compare heat to heat-drought
list21 <- unique(bottom_heatdrt_damage$Accession)
bottom_heat_damage$compare.geno <- bottom_heat_damage$Accession %in% list21
all_overlap_stressed_bottom <- filter(bottom_heat_damage, compare.geno == TRUE)
length(unique(all_overlap_stressed_bottom$Accession))
all_overlap_stressed_bottom = dplyr::select(all_overlap_stressed_bottom, -compare.geno)
list22 <- unique(all_overlap_stressed_bottom$Accession)

#compare heat and drought
list21 <- unique(bottom_heat_damage$Accession)
bottom_drt_damage$compare.geno <- bottom_drt_damage$Accession %in% list21
all_overlap_stressed_bottom <- filter(bottom_drt_damage, compare.geno == TRUE)
length(unique(all_overlap_stressed_bottom$Accession))
all_overlap_stressed_bottom = dplyr::select(all_overlap_stressed_bottom, -compare.geno)
list23 <- unique(all_overlap_stressed_bottom$Accession)

#compare heat-drought and drought
list21 <- unique(bottom_heatdrt_damage$Accession)
bottom_drt_damage$compare.geno <- bottom_drt_damage$Accession %in% list21
all_overlap_stressed_bottom <- filter(bottom_drt_damage, compare.geno == TRUE)
length(unique(all_overlap_stressed_bottom$Accession))
all_overlap_stressed_bottom = dplyr::select(all_overlap_stressed_bottom, -compare.geno)
list24 <- unique(all_overlap_stressed_bottom$Accession)

#compare heat/heat-drought to drought to find overall overlap
bottom_drt_damage$compare.geno <- bottom_drt_damage$Accession %in% list22
all_overlap_stressed_bottom <- filter(bottom_drt_damage, compare.geno == TRUE)


#####

#control damage = ctl_damage (137, 34)
#drought damage = drt_damage (137, 34)
#heat damage = heat_damage (143, 35)
#heat-drought damage = heatdrt_damage (139, 34)


# look at overlap using actual data values
top_ctl_damage <- ctl_damage %>% top_n(-34)
bottom_ctl_damage <- ctl_damage %>% top_n(34)

top_drt_damage <- drt_damage %>% top_n(-34)
bottom_drt_damage <- drt_damage %>% top_n(34)

top_heat_damage <- heat_damage %>% top_n(-35)
bottom_heat_damage <- heat_damage %>% top_n(35)

top_heatdrt_damage <- heatdrt_damage %>% top_n(-34)
bottom_heatdrt_damage <- heatdrt_damage %>% top_n(34)

### See which accessions overlap that perform poorly ###
#compare control and drought
list1 <- unique(top_ctl_damage$Accession)
top_drt_damage$compare.geno <- top_drt_damage$Accession %in% list1
all_overlap_stressed_top2 <- filter(top_drt_damage, compare.geno == TRUE)
length(unique(all_overlap_stressed_top2$Accession))
all_overlap_stressed_top2 = dplyr::select(all_overlap_stressed_top2, -compare.geno)
list5 <- unique(all_overlap_stressed_top2$Accession)

# compare control to heat
list1 <- unique(top_ctl_damage$Accession)
top_heat_damage$compare.geno <- top_heat_damage$Accession %in% list1
all_overlap_stressed_top2 <- filter(top_heat_damage, compare.geno == TRUE)
length(unique(all_overlap_stressed_top2$Accession))
all_overlap_stressed_top2 = dplyr::select(all_overlap_stressed_top2, -compare.geno)
list6 <- unique(all_overlap_stressed_top2$Accession)

# compare control to heat-drought
list1 <- unique(top_ctl_damage$Accession)
top_heatdrt_damage$compare.geno <- top_heatdrt_damage$Accession %in% list1
all_overlap_stressed_top2 <- filter(top_heatdrt_damage, compare.geno == TRUE)
length(unique(all_overlap_stressed_top2$Accession))
all_overlap_stressed_top2 = dplyr::select(all_overlap_stressed_top2, -compare.geno)
list7 <- unique(all_overlap_stressed_top2$Accession)

#compare control/drought to heat to find overlap
top_heat_damage$compare.geno <- top_heat_damage$Accession %in% list5
all_overlap_stressed_top2 <- filter(top_heat_damage, compare.geno == TRUE)
list8 <- unique(all_overlap_stressed_top2$Accession)

#compare control/drought to heat-drought to find overlap
top_heatdrt_damage$compare.geno <- top_heatdrt_damage$Accession %in% list5
all_overlap_stressed_top2 <- filter(top_heatdrt_damage, compare.geno == TRUE)
list9 <- unique(all_overlap_stressed_top2$Accession)

#compare control/heat to heat-drought to find overlap
top_heatdrt_damage$compare.geno <- top_heatdrt_damage$Accession %in% list6
all_overlap_stressed_top2 <- filter(top_heatdrt_damage, compare.geno == TRUE)
list10 <- unique(all_overlap_stressed_top2$Accession)

#compare control/heat & heat-drought with drought to find overlap
top_drt_damage$compare.geno <- top_drt_damage$Accession %in% list10
all_overlap_stressed_top2 <- filter(top_drt_damage, compare.geno == TRUE)



### See which accessions overlap that perform well ###
#compare control and drought
list21 <- unique(bottom_ctl_damage$Accession)
bottom_drt_damage$compare.geno <- bottom_drt_damage$Accession %in% list21
all_overlap_stressed_bottom2 <- filter(bottom_drt_damage, compare.geno == TRUE)
length(unique(all_overlap_stressed_bottom2$Accession))
all_overlap_stressed_bottom2 = dplyr::select(all_overlap_stressed_bottom2, -compare.geno)
list25 <- unique(all_overlap_stressed_bottom2$Accession)

# compare control to heat
list21 <- unique(bottom_ctl_damage$Accession)
bottom_heat_damage$compare.geno <- bottom_heat_damage$Accession %in% list21
all_overlap_stressed_bottom2 <- filter(bottom_heat_damage, compare.geno == TRUE)
length(unique(all_overlap_stressed_bottom2$Accession))
all_overlap_stressed_bottom2 = dplyr::select(all_overlap_stressed_bottom2, -compare.geno)
list26 <- unique(all_overlap_stressed_bottom2$Accession)

# compare control to heat-drought
list21 <- unique(bottom_ctl_damage$Accession)
bottom_heatdrt_damage$compare.geno <- bottom_heatdrt_damage$Accession %in% list21
all_overlap_stressed_bottom2 <- filter(bottom_heatdrt_damage, compare.geno == TRUE)
length(unique(all_overlap_stressed_bottom2$Accession))
all_overlap_stressed_bottom2 = dplyr::select(all_overlap_stressed_bottom2, -compare.geno)
list27 <- unique(all_overlap_stressed_bottom2$Accession)

#compare control/drought to heat to find overlap
bottom_heat_damage$compare.geno <- bottom_heat_damage$Accession %in% list25
all_overlap_stressed_bottom2 <- filter(bottom_heat_damage, compare.geno == TRUE)
list28 <- unique(all_overlap_stressed_bottom2$Accession)

#compare control/drought to heat-drought to find overlap
bottom_heatdrt_damage$compare.geno <- bottom_heatdrt_damage$Accession %in% list25
all_overlap_stressed_bottom2 <- filter(bottom_heatdrt_damage, compare.geno == TRUE)
list29 <- unique(all_overlap_stressed_bottom2$Accession)

#compare control/heat to heat-drought to find overlap
bottom_heatdrt_damage$compare.geno <- bottom_heatdrt_damage$Accession %in% list26
all_overlap_stressed_bottom2 <- filter(bottom_heatdrt_damage, compare.geno == TRUE)
list30 <- unique(all_overlap_stressed_bottom2$Accession)

#compare control/heat & heat-drought with drought to find overlap
bottom_drt_damage$compare.geno <- bottom_drt_damage$Accession %in% list30
all_overlap_stressed_bottom2 <- filter(bottom_drt_damage, compare.geno == TRUE)

#******************************************
############################################
### Biomass/Area Comparisons ###
############################################



###*************************
# look at overlap of difference between stresses and control
#area_drought (129)
#area_heat (129)
#area_heatdrought (129)

area_drought <- as.data.frame(area_drought)
area_heat <- as.data.frame(area_heat)
area_heatdrought <- as.data.frame(area_heatdrought)

area_drought <- drop_na(dplyr::select(area_drought, Accession, area_dif))
top_drt_area <- area_drought %>% top_n(-32)
bottom_drt_area <- area_drought %>% top_n(32)

area_heat <- drop_na(dplyr::select(area_heat, Accession, area_dif))
top_heat_area <- area_heat %>% top_n(-32)
bottom_heat_area <- area_heat %>% top_n(32)

area_heatdrought <- drop_na(dplyr::select(area_heatdrought, Accession, area_dif))
top_heatdrt_area <- area_heatdrought %>% top_n(-32)
bottom_heatdrt_area <- area_heatdrought %>% top_n(32)


### See which accessions overlap that perform well ###

# compare heat to heat-drought
list31 <- unique(top_heatdrt_area$Accession)
top_heat_area$compare.geno <- top_heat_area$Accession %in% list31
all_overlap_area_top <- filter(top_heat_area, compare.geno == TRUE)
length(unique(all_overlap_area_top$Accession))
all_overlap_area_top = dplyr::select(all_overlap_area_top, -compare.geno)
list32 <- unique(all_overlap_area_top$Accession)

#compare heat and drought
list31 <- unique(top_heat_area$Accession)
top_drt_area$compare.geno <- top_drt_area$Accession %in% list31
all_overlap_area_top <- filter(top_drt_area, compare.geno == TRUE)
length(unique(all_overlap_area_top$Accession))
all_overlap_area_top = dplyr::select(all_overlap_area_top, -compare.geno)
list33 <- unique(all_overlap_area_top$Accession)

#compare heat-drought and drought
list31 <- unique(top_heatdrt_area$Accession)
top_drt_area$compare.geno <- top_drt_area$Accession %in% list31
all_overlap_area_top <- filter(top_drt_area, compare.geno == TRUE)
length(unique(all_overlap_area_top$Accession))
all_overlap_area_top = dplyr::select(all_overlap_area_top, -compare.geno)
list34 <- unique(all_overlap_area_top$Accession)

#compare heat/heat-drought to drought to find overall overlap
top_drt_area$compare.geno <- top_drt_area$Accession %in% list32
all_overlap_area_top <- filter(top_drt_area, compare.geno == TRUE)


### See which accessions overlap that perform poorly ###
# compare heat to heat-drought
list41 <- unique(bottom_heatdrt_area$Accession)
bottom_heat_area$compare.geno <- bottom_heat_area$Accession %in% list41
all_overlap_area_bottom <- filter(bottom_heat_area, compare.geno == TRUE)
length(unique(all_overlap_area_bottom$Accession))
all_overlap_area_bottom = dplyr::select(all_overlap_area_bottom, -compare.geno)
list42 <- unique(all_overlap_area_bottom$Accession)

#compare heat and drought
list41 <- unique(bottom_heat_area$Accession)
bottom_drt_area$compare.geno <- bottom_drt_area$Accession %in% list41
all_overlap_area_bottom <- filter(bottom_drt_area, compare.geno == TRUE)
length(unique(all_overlap_area_bottom$Accession))
all_overlap_area_bottom = dplyr::select(all_overlap_area_bottom, -compare.geno)
list43 <- unique(all_overlap_area_bottom$Accession)

#compare heat-drought and drought
list41 <- unique(bottom_heatdrt_area$Accession)
bottom_drt_area$compare.geno <- bottom_drt_area$Accession %in% list41
all_overlap_area_bottom <- filter(bottom_drt_area, compare.geno == TRUE)
length(unique(all_overlap_area_bottom$Accession))
all_overlap_area_bottom = dplyr::select(all_overlap_area_bottom, -compare.geno)
list44 <- unique(all_overlap_area_bottom$Accession)

#compare heat/heat-drought to drought to find overall overlap
bottom_drt_area$compare.geno <- bottom_drt_area$Accession %in% list42
all_overlap_area_bottom <- filter(bottom_drt_area, compare.geno == TRUE)




#####

#control area = ctl_area (137, 34)
#drought area = drt_area (137, 34)
#heat area = heat_area (143, 35)
#heat-drought area = heatdrt_area (139, 34)

# look at overlap using actual data values
top_ctl_area <- ctl_area %>% top_n(34)
bottom_ctl_area <- ctl_area %>% top_n(-34)

top_drt_area <- drt_area %>% top_n(34)
bottom_drt_area <- drt_area %>% top_n(-34)

top_heat_area <- heat_area %>% top_n(35)
bottom_heat_area <- heat_area %>% top_n(-35)

top_heatdrt_area <- heatdrt_area %>% top_n(34)
bottom_heatdrt_area <- heatdrt_area %>% top_n(-34)

### See which accessions overlap that perform well ###
#compare control and drought
list31 <- unique(top_ctl_area$Accession)
top_drt_area$compare.geno <- top_drt_area$Accession %in% list31
all_overlap_area_top2 <- filter(top_drt_area, compare.geno == TRUE)
length(unique(all_overlap_area_top2$Accession))
all_overlap_area_top2 = dplyr::select(all_overlap_area_top2, -compare.geno)
list35 <- unique(all_overlap_area_top2$Accession)

# compare control to heat
list31 <- unique(top_ctl_area$Accession)
top_heat_area$compare.geno <- top_heat_area$Accession %in% list31
all_overlap_area_top2 <- filter(top_heat_area, compare.geno == TRUE)
length(unique(all_overlap_area_top2$Accession))
all_overlap_area_top2 = dplyr::select(all_overlap_area_top2, -compare.geno)
list36 <- unique(all_overlap_area_top2$Accession)

# compare control to heat-drought
list31 <- unique(top_ctl_area$Accession)
top_heatdrt_area$compare.geno <- top_heatdrt_area$Accession %in% list31
all_overlap_area_top2 <- filter(top_heatdrt_area, compare.geno == TRUE)
length(unique(all_overlap_area_top2$Accession))
all_overlap_area_top2 = dplyr::select(all_overlap_area_top2, -compare.geno)
list37 <- unique(all_overlap_area_top2$Accession)

#compare control/drought to heat to find overlap
top_heat_area$compare.geno <- top_heat_area$Accession %in% list35
all_overlap_area_top2 <- filter(top_heat_area, compare.geno == TRUE)
list38 <- unique(all_overlap_area_top2$Accession)

#compare control/drought to heat-drought to find overlap
top_heatdrt_area$compare.geno <- top_heatdrt_area$Accession %in% list35
all_overlap_area_top2 <- filter(top_heatdrt_area, compare.geno == TRUE)
list39 <- unique(all_overlap_area_top2$Accession)

#compare control/heat to heat-drought to find overlap
top_heatdrt_area$compare.geno <- top_heatdrt_area$Accession %in% list36
all_overlap_area_top2 <- filter(top_heatdrt_area, compare.geno == TRUE)
list40 <- unique(all_overlap_area_top2$Accession)

#compare control/heat & heat-drought with drought to find overlap
top_drt_area$compare.geno <- top_drt_area$Accession %in% list40
all_overlap_area_top2 <- filter(top_drt_area, compare.geno == TRUE)




### See which accessions overlap that perform poorly ###

#compare control and drought
list41 <- unique(bottom_ctl_area$Accession)
bottom_drt_area$compare.geno <- bottom_drt_area$Accession %in% list41
all_overlap_area_bottom2 <- filter(bottom_drt_area, compare.geno == TRUE)
length(unique(all_overlap_area_bottom2$Accession))
all_overlap_area_bottom2 = dplyr::select(all_overlap_area_bottom2, -compare.geno)
list45 <- unique(all_overlap_area_bottom2$Accession)

# compare control to heat
list41 <- unique(bottom_ctl_area$Accession)
bottom_heat_area$compare.geno <- bottom_heat_area$Accession %in% list41
all_overlap_area_bottom2 <- filter(bottom_heat_area, compare.geno == TRUE)
length(unique(all_overlap_area_bottom2$Accession))
all_overlap_area_bottom2 = dplyr::select(all_overlap_area_bottom2, -compare.geno)
list46 <- unique(all_overlap_area_bottom2$Accession)

# compare control to heat-drought
list41 <- unique(bottom_ctl_area$Accession)
bottom_heatdrt_area$compare.geno <- bottom_heatdrt_area$Accession %in% list41
all_overlap_area_bottom2 <- filter(bottom_heatdrt_area, compare.geno == TRUE)
length(unique(all_overlap_area_bottom2$Accession))
all_overlap_area_bottom2 = dplyr::select(all_overlap_area_bottom2, -compare.geno)
list47 <- unique(all_overlap_area_bottom2$Accession)

#compare control/drought to heat to find overlap
bottom_heat_area$compare.geno <- bottom_heat_area$Accession %in% list45
all_overlap_area_bottom2 <- filter(bottom_heat_area, compare.geno == TRUE)
list48 <- unique(all_overlap_area_bottom2$Accession)

#compare control/drought to heat-drought to find overlap
bottom_heatdrt_area$compare.geno <- bottom_heatdrt_area$Accession %in% list45
all_overlap_area_bottom2 <- filter(bottom_heatdrt_area, compare.geno == TRUE)
list49 <- unique(all_overlap_area_bottom2$Accession)

#compare control/heat to heat-drought to find overlap
bottom_heatdrt_area$compare.geno <- bottom_heatdrt_area$Accession %in% list46
all_overlap_area_bottom2 <- filter(bottom_heatdrt_area, compare.geno == TRUE)
list50 <- unique(all_overlap_area_bottom2$Accession)

#compare control/heat & heat-drought with drought to find overlap
bottom_drt_area$compare.geno <- bottom_drt_area$Accession %in% list50
all_overlap_area_bottom2 <- filter(bottom_drt_area, compare.geno == TRUE)




############################################
### Height Comparisons ###
############################################
### Height Data ###

# look at overlap of difference between stresses and control
#height_drought 
#height_heat
#height_heatdrought

height_drought <- as.data.frame(height_drought)
height_heat <- as.data.frame(height_heat)
height_heatdrought <- as.data.frame(height_heatdrought)

height_drought <- drop_na(dplyr::select(height_drought, Accession, height_dif))
top_drt_height <- height_drought %>% top_n(-32)
bottom_drt_height <- height_drought %>% top_n(32)

height_heat <- drop_na(dplyr::select(height_heat, Accession, height_dif))
top_heat_height <- height_heat %>% top_n(-32)
bottom_heat_height <- height_heat %>% top_n(32)

height_heatdrought <- drop_na(dplyr::select(height_heatdrought, Accession, height_dif))
top_heatdrt_height <- height_heatdrought %>% top_n(-32)
bottom_heatdrt_height <- height_heatdrought %>% top_n(32)


### See which accessions overlap that perform well ###

# compare heat to heat-drought
list51 <- unique(top_heatdrt_height$Accession)
top_heat_height$compare.geno <- top_heat_height$Accession %in% list51
all_overlap_height_top <- filter(top_heat_height, compare.geno == TRUE)
length(unique(all_overlap_height_top$Accession))
all_overlap_height_top = dplyr::select(all_overlap_height_top, -compare.geno)
list52 <- unique(all_overlap_height_top$Accession)

#compare heat and drought
list51 <- unique(top_heat_height$Accession)
top_drt_height$compare.geno <- top_drt_height$Accession %in% list51
all_overlap_height_top <- filter(top_drt_height, compare.geno == TRUE)
length(unique(all_overlap_height_top$Accession))
all_overlap_height_top = dplyr::select(all_overlap_height_top, -compare.geno)
list53 <- unique(all_overlap_height_top$Accession)

#compare heat-drought and drought
list51 <- unique(top_heatdrt_height$Accession)
top_drt_height$compare.geno <- top_drt_height$Accession %in% list51
all_overlap_height_top <- filter(top_drt_height, compare.geno == TRUE)
length(unique(all_overlap_height_top$Accession))
all_overlap_height_top = dplyr::select(all_overlap_height_top, -compare.geno)
list54 <- unique(all_overlap_height_top$Accession)

#compare heat/heat-drought to drought to find overall overlap
top_drt_height$compare.geno <- top_drt_height$Accession %in% list52
all_overlap_height_top <- filter(top_drt_height, compare.geno == TRUE)



### See which accessions overlap that perform poorly ###

# compare heat to heat-drought
list61 <- unique(bottom_heatdrt_height$Accession)
bottom_heat_height$compare.geno <- bottom_heat_height$Accession %in% list61
all_overlap_height_bottom <- filter(bottom_heat_height, compare.geno == TRUE)
length(unique(all_overlap_height_bottom$Accession))
all_overlap_height_bottom = dplyr::select(all_overlap_height_bottom, -compare.geno)
list62 <- unique(all_overlap_height_bottom$Accession)

#compare heat and drought
list61 <- unique(bottom_heat_height$Accession)
bottom_drt_height$compare.geno <- bottom_drt_height$Accession %in% list61
all_overlap_height_bottom <- filter(bottom_drt_height, compare.geno == TRUE)
length(unique(all_overlap_height_bottom$Accession))
all_overlap_height_bottom = dplyr::select(all_overlap_height_bottom, -compare.geno)
list63 <- unique(all_overlap_height_bottom$Accession)

#compare heat-drought and drought
list61 <- unique(bottom_heatdrt_height$Accession)
bottom_drt_height$compare.geno <- bottom_drt_height$Accession %in% list61
all_overlap_height_bottom <- filter(bottom_drt_height, compare.geno == TRUE)
length(unique(all_overlap_height_bottom$Accession))
all_overlap_height_bottom = dplyr::select(all_overlap_height_bottom, -compare.geno)
list64 <- unique(all_overlap_height_bottom$Accession)

#compare heat/heat-drought to drought to find overall overlap
bottom_drt_height$compare.geno <- bottom_drt_height$Accession %in% list62
all_overlap_height_bottom <- filter(bottom_drt_height, compare.geno == TRUE)

#####


#control height = ctl_height (137, 34)
#drought height = drt_height (137, 34)
#heat height = heat_height (143, 35)
#heat-drought height = heatdrt_height (139, 34)

# look at overlap using actual data values
top_ctl_height <- ctl_height %>% top_n(34)
bottom_ctl_height <- ctl_height %>% top_n(-34)

top_drt_height <- drt_height %>% top_n(34)
bottom_drt_height <- drt_height %>% top_n(-34)

top_heat_height <- heat_height %>% top_n(35)
bottom_heat_height <- heat_height %>% top_n(-35)

top_heatdrt_height <- heatdrt_height %>% top_n(34)
bottom_heatdrt_height <- heatdrt_height %>% top_n(-34)

### See which accessions overlap that perform well ###

#compare control and drought
list51 <- unique(top_ctl_height$Accession)
top_drt_height$compare.geno <- top_drt_height$Accession %in% list51
all_overlap_height_top2 <- filter(top_drt_height, compare.geno == TRUE)
length(unique(all_overlap_height_top2$Accession))
all_overlap_height_top2 = dplyr::select(all_overlap_height_top2, -compare.geno)
list55 <- unique(all_overlap_height_top2$Accession)

# compare control to heat
list51 <- unique(top_ctl_height$Accession)
top_heat_height$compare.geno <- top_heat_height$Accession %in% list51
all_overlap_height_top2 <- filter(top_heat_height, compare.geno == TRUE)
length(unique(all_overlap_height_top2$Accession))
all_overlap_height_top2 = dplyr::select(all_overlap_height_top2, -compare.geno)
list56 <- unique(all_overlap_height_top2$Accession)

# compare control to heat-drought
list51 <- unique(top_ctl_height$Accession)
top_heatdrt_height$compare.geno <- top_heatdrt_height$Accession %in% list51
all_overlap_height_top2 <- filter(top_heatdrt_height, compare.geno == TRUE)
length(unique(all_overlap_height_top2$Accession))
all_overlap_height_top2 = dplyr::select(all_overlap_height_top2, -compare.geno)
list57 <- unique(all_overlap_height_top2$Accession)

#compare control/drought to heat to find overlap
top_heat_height$compare.geno <- top_heat_height$Accession %in% list55
all_overlap_height_top2 <- filter(top_heat_height, compare.geno == TRUE)
list58 <- unique(all_overlap_height_top2$Accession)

#compare control/drought to heat-drought to find overlap
top_heatdrt_height$compare.geno <- top_heatdrt_height$Accession %in% list55
all_overlap_height_top2 <- filter(top_heatdrt_height, compare.geno == TRUE)
list59 <- unique(all_overlap_height_top2$Accession)

#compare control/heat to heat-drought to find overlap
top_heatdrt_height$compare.geno <- top_heatdrt_height$Accession %in% list56
all_overlap_height_top2 <- filter(top_heatdrt_height, compare.geno == TRUE)
list60 <- unique(all_overlap_height_top2$Accession)

#compare control/heat & heat-drought with drought to find overlap
top_drt_height$compare.geno <- top_drt_height$Accession %in% list60
all_overlap_height_top2 <- filter(top_drt_height, compare.geno == TRUE)

### See which accessions overlap that perform poorly ###

#compare control and drought
list61 <- unique(bottom_ctl_height$Accession)
bottom_drt_height$compare.geno <- bottom_drt_height$Accession %in% list61
all_overlap_height_bottom2 <- filter(bottom_drt_height, compare.geno == TRUE)
length(unique(all_overlap_height_bottom2$Accession))
all_overlap_height_bottom2 = dplyr::select(all_overlap_height_bottom2, -compare.geno)
list65 <- unique(all_overlap_height_bottom2$Accession)

# compare control to heat
list61 <- unique(bottom_ctl_height$Accession)
bottom_heat_height$compare.geno <- bottom_heat_height$Accession %in% list61
all_overlap_height_bottom2 <- filter(bottom_heat_height, compare.geno == TRUE)
length(unique(all_overlap_height_bottom2$Accession))
all_overlap_height_bottom2 = dplyr::select(all_overlap_height_bottom2, -compare.geno)
list66 <- unique(all_overlap_height_bottom2$Accession)

# compare control to heat-drought
list61 <- unique(bottom_ctl_height$Accession)
bottom_heatdrt_height$compare.geno <- bottom_heatdrt_height$Accession %in% list61
all_overlap_height_bottom2 <- filter(bottom_heatdrt_height, compare.geno == TRUE)
length(unique(all_overlap_height_bottom2$Accession))
all_overlap_height_bottom2 = dplyr::select(all_overlap_height_bottom2, -compare.geno)
list67 <- unique(all_overlap_height_bottom2$Accession)

#compare control/drought to heat to find overlap
bottom_heat_height$compare.geno <- bottom_heat_height$Accession %in% list65
all_overlap_height_bottom2 <- filter(bottom_heat_height, compare.geno == TRUE)
list68 <- unique(all_overlap_height_bottom2$Accession)

#compare control/drought to heat-drought to find overlap
bottom_heatdrt_height$compare.geno <- bottom_heatdrt_height$Accession %in% list65
all_overlap_height_bottom2 <- filter(bottom_heatdrt_height, compare.geno == TRUE)
list69 <- unique(all_overlap_height_bottom2$Accession)

#compare control/heat to heat-drought to find overlap
bottom_heatdrt_height$compare.geno <- bottom_heatdrt_height$Accession %in% list66
all_overlap_height_bottom2 <- filter(bottom_heatdrt_height, compare.geno == TRUE)
list70 <- unique(all_overlap_height_bottom2$Accession)

#compare control/heat & heat-drought with drought to find overlap
bottom_drt_height$compare.geno <- bottom_drt_height$Accession %in% list70
all_overlap_height_bottom2 <- filter(bottom_drt_height, compare.geno == TRUE)


#############################################
### Overall Comparisons ###
#############################################
# Only Stress Treatments
# Get overall list of accessions that perform well for ANY phenotype
all_overlap_any_top <- merge(all_overlap_area_top, all_overlap_height_top, by = "Accession", all = TRUE)
all_overlap_any_top <- merge(all_overlap_any_top, all_overlap_stressed_top, by = "Accession", all = TRUE)
all_overlap_any_top <- dplyr::select(all_overlap_any_top, Accession)

# Get overall list of accessions that perform poorly for ANY phenotype
all_overlap_any_bottom <- merge(all_overlap_area_bottom, all_overlap_height_bottom, by = "Accession", all = TRUE)
all_overlap_any_bottom <- merge(all_overlap_any_bottom, all_overlap_stressed_bottom, by = "Accession", all = TRUE)
all_overlap_any_bottom <- dplyr::select(all_overlap_any_bottom, Accession)

# Get overall list of accessions that perform well for ALL phenotypes
all_overlap_all_top <- merge(all_overlap_area_top, all_overlap_height_top, by = "Accession", all = FALSE)
all_overlap_all_top <- merge(all_overlap_all_top, all_overlap_stressed_top, by = "Accession", all = FALSE)
all_overlap_all_top <- dplyr::select(all_overlap_all_top, Accession)

# Get overall list of accessions that perform poorly for ALL phenotypes
all_overlap_all_bottom <- merge(all_overlap_area_bottom, all_overlap_height_bottom, by = "Accession", all = FALSE)
all_overlap_all_bottom <- merge(all_overlap_all_bottom, all_overlap_stressed_bottom, by = "Accession", all = FALSE)
all_overlap_all_bottom <- dplyr::select(all_overlap_all_bottom, Accession)

#####
# Including Control
# Get overall list of accessions that perform well for ANY phenotype
all_overlap_any_top2 <- merge(all_overlap_area_top2, all_overlap_height_top2, by = "Accession", all = TRUE)
all_overlap_any_top2 <- merge(all_overlap_any_top2, all_overlap_stressed_top2, by = "Accession", all = TRUE)
all_overlap_any_top2 <- dplyr::select(all_overlap_any_top2, Accession)

# Get overall list of accessions that perform poorly for ANY phenotype
all_overlap_any_bottom2 <- merge(all_overlap_area_bottom2, all_overlap_height_bottom2, by = "Accession", all = TRUE)
all_overlap_any_bottom2 <- merge(all_overlap_any_bottom2, all_overlap_stressed_bottom2, by = "Accession", all = TRUE)
all_overlap_any_bottom2 <- dplyr::select(all_overlap_any_bottom2, Accession)

# Get overall list of accessions that perform well for ALL phenotypes
all_overlap_all_top2 <- merge(all_overlap_area_top2, all_overlap_height_top2, by = "Accession", all = FALSE)
all_overlap_all_top2 <- merge(all_overlap_all_top2, all_overlap_stressed_top2, by = "Accession", all = FALSE)
all_overlap_all_top2 <- dplyr::select(all_overlap_all_top2, Accession)

# Get overall list of accessions that perform poorly for ALL phenotypes
all_overlap_all_bottom2 <- merge(all_overlap_area_bottom2, all_overlap_height_bottom2, by = "Accession", all = FALSE)
all_overlap_all_bottom2 <- merge(all_overlap_all_bottom2, all_overlap_stressed_bottom2, by = "Accession", all = FALSE)
all_overlap_all_bottom2 <- dplyr::select(all_overlap_all_bottom2, Accession)







##### 
#********************************************
### Find overlap between different phenotypes within treatments ###
#********************************************

#############################################
### Overlap in Control ###
#############################################

#top_ctl_damage <- ctl_damage %>% top_n(34) BAD
#bottom_ctl_damage <- ctl_damage %>% top_n(-34) GOOD

#top_ctl_area <- ctl_area %>% top_n(34) GOOD
#bottom_ctl_area <- ctl_area %>% top_n(-34) BAD

#top_ctl_height <- ctl_height %>% top_n(34) GOOD
#bottom_ctl_height <- ctl_height %>% top_n(-34) BAD


### See which accessions overlap that perform well ###

# compare damage to area 
list71 <- unique(top_ctl_area$Accession)
top_ctl_damage$compare.geno <- top_ctl_damage$Accession %in% list71
all_overlap_ctl_top <- filter(top_ctl_damage, compare.geno == TRUE)
length(unique(all_overlap_ctl_top$Accession))
all_overlap_ctl_top = dplyr::select(all_overlap_ctl_top, -compare.geno)
list72 <- unique(all_overlap_ctl_top$Accession)

#compare height to damage
list71 <- unique(top_ctl_damage$Accession)
top_ctl_height$compare.geno <- top_ctl_height$Accession %in% list71
all_overlap_ctl_top <- filter(top_ctl_height, compare.geno == TRUE)
length(unique(all_overlap_ctl_top$Accession))
all_overlap_ctl_top = dplyr::select(all_overlap_ctl_top, -compare.geno)
list73 <- unique(all_overlap_ctl_top$Accession)

#compare height to area
list71 <- unique(top_ctl_area$Accession)
top_ctl_height$compare.geno <- top_ctl_height$Accession %in% list71
all_overlap_ctl_top <- filter(top_ctl_height, compare.geno == TRUE)
length(unique(all_overlap_ctl_top$Accession))
all_overlap_ctl_top = dplyr::select(all_overlap_ctl_top, -compare.geno)
list74 <- unique(all_overlap_ctl_top$Accession)

#compare damage/area to height to find overall overlap
top_drt_height$compare.geno <- top_drt_height$Accession %in% list72
all_overlap_ctl_top <- filter(top_drt_height, compare.geno == TRUE)


### See which accessions overlap that perform poorly ###

# compare damage to area
list81 <- unique(bottom_ctl_area$Accession)
bottom_ctl_damage$compare.geno <- bottom_ctl_damage$Accession %in% list81
all_overlap_ctl_bottom <- filter(bottom_ctl_damage, compare.geno == TRUE)
length(unique(all_overlap_ctl_bottom$Accession))
all_overlap_ctl_bottom = dplyr::select(all_overlap_ctl_bottom, -compare.geno)
list82 <- unique(all_overlap_ctl_bottom$Accession)

#compare damage and height
list81 <- unique(bottom_ctl_damage$Accession)
bottom_ctl_height$compare.geno <- bottom_ctl_height$Accession %in% list81
all_overlap_ctl_bottom <- filter(bottom_ctl_height, compare.geno == TRUE)
length(unique(all_overlap_ctl_bottom$Accession))
all_overlap_ctl_bottom = dplyr::select(all_overlap_ctl_bottom, -compare.geno)
list83 <- unique(all_overlap_ctl_bottom$Accession)

#compare area and height
list81 <- unique(bottom_ctl_area$Accession)
bottom_ctl_height$compare.geno <- bottom_ctl_height$Accession %in% list81
all_overlap_ctl_bottom <- filter(bottom_ctl_height, compare.geno == TRUE)
length(unique(all_overlap_ctl_bottom$Accession))
all_overlap_ctl_bottom = dplyr::select(all_overlap_ctl_bottom, -compare.geno)
list84 <- unique(all_overlap_ctl_bottom$Accession)

#compare damage/area to height to find overall overlap
bottom_ctl_height$compare.geno <- bottom_ctl_height$Accession %in% list82
all_overlap_ctl_bottom <- filter(bottom_ctl_height, compare.geno == TRUE)


#############################################
### Overlap in Drought ###
#############################################

#top_drt_damage <- drt_damage %>% top_n(34)
#bottom_drt_damage <- drt_damage %>% top_n(-34)

#top_drt_area <- drt_area %>% top_n(34)
#bottom_drt_area <- drt_area %>% top_n(-34)

#top_drt_height <- drt_height %>% top_n(34)
#bottom_drt_height <- drt_height %>% top_n(-34)


### See which accessions overlap that perform well ###

# compare damage to area 
list91 <- unique(top_drt_area$Accession)
top_drt_damage$compare.geno <- top_drt_damage$Accession %in% list91
all_overlap_drt_top <- filter(top_drt_damage, compare.geno == TRUE)
length(unique(all_overlap_drt_top$Accession))
all_overlap_drt_top = dplyr::select(all_overlap_drt_top, -compare.geno)
list92 <- unique(all_overlap_drt_top$Accession)

#compare height to damage
list91 <- unique(top_drt_damage$Accession)
top_drt_height$compare.geno <- top_drt_height$Accession %in% list91
all_overlap_drt_top <- filter(top_drt_height, compare.geno == TRUE)
length(unique(all_overlap_drt_top$Accession))
all_overlap_drt_top = dplyr::select(all_overlap_drt_top, -compare.geno)
list93 <- unique(all_overlap_drt_top$Accession)

#compare height to area
list91 <- unique(top_drt_area$Accession)
top_drt_height$compare.geno <- top_drt_height$Accession %in% list91
all_overlap_drt_top <- filter(top_drt_height, compare.geno == TRUE)
length(unique(all_overlap_drt_top$Accession))
all_overlap_drt_top = dplyr::select(all_overlap_drt_top, -compare.geno)
list94 <- unique(all_overlap_drt_top$Accession)

#compare damage/area to height to find overall overlap
top_drt_height$compare.geno <- top_drt_height$Accession %in% list92
all_overlap_drt_top <- filter(top_drt_height, compare.geno == TRUE)


### See which accessions overlap that perform poorly ###

# compare damage to area
list101 <- unique(bottom_drt_area$Accession)
bottom_drt_damage$compare.geno <- bottom_drt_damage$Accession %in% list101
all_overlap_drt_bottom <- filter(bottom_drt_damage, compare.geno == TRUE)
length(unique(all_overlap_drt_bottom$Accession))
all_overlap_drt_bottom = dplyr::select(all_overlap_drt_bottom, -compare.geno)
list102 <- unique(all_overlap_drt_bottom$Accession)

#compare damage and height
list101 <- unique(bottom_drt_damage$Accession)
bottom_drt_height$compare.geno <- bottom_drt_height$Accession %in% list101
all_overlap_drt_bottom <- filter(bottom_drt_height, compare.geno == TRUE)
length(unique(all_overlap_drt_bottom$Accession))
all_overlap_drt_bottom = dplyr::select(all_overlap_drt_bottom, -compare.geno)
list103 <- unique(all_overlap_drt_bottom$Accession)

#compare area and height
list101 <- unique(bottom_drt_area$Accession)
bottom_drt_height$compare.geno <- bottom_drt_height$Accession %in% list101
all_overlap_drt_bottom <- filter(bottom_drt_height, compare.geno == TRUE)
length(unique(all_overlap_drt_bottom$Accession))
all_overlap_drt_bottom = dplyr::select(all_overlap_drt_bottom, -compare.geno)
list104 <- unique(all_overlap_drt_bottom$Accession)

#compare damage/area to height to find overall overlap
bottom_drt_height$compare.geno <- bottom_drt_height$Accession %in% list102
all_overlap_drt_bottom <- filter(bottom_drt_height, compare.geno == TRUE)



#############################################
### Overlap in Heat ###
#############################################

#top_heat_damage <- heat_damage %>% top_n(35)
#bottom_heat_damage <- heat_damage %>% top_n(-35)

#top_heat_area <- heat_area %>% top_n(35)
#bottom_heat_area <- heat_area %>% top_n(-35)

#top_heat_height <- heat_height %>% top_n(35)
#bottom_heat_height <- heat_height %>% top_n(-35)


### See which accessions overlap that perform well ###

# compare damage to area 
list111 <- unique(top_heat_area$Accession)
top_heat_damage$compare.geno <- top_heat_damage$Accession %in% list111
all_overlap_heat_top <- filter(top_heat_damage, compare.geno == TRUE)
length(unique(all_overlap_heat_top$Accession))
all_overlap_heat_top = dplyr::select(all_overlap_heat_top, -compare.geno)
list112 <- unique(all_overlap_heat_top$Accession)

#compare height to damage
list111 <- unique(top_heat_damage$Accession)
top_heat_height$compare.geno <- top_heat_height$Accession %in% list111
all_overlap_heat_top <- filter(top_heat_height, compare.geno == TRUE)
length(unique(all_overlap_heat_top$Accession))
all_overlap_heat_top = dplyr::select(all_overlap_heat_top, -compare.geno)
list113 <- unique(all_overlap_heat_top$Accession)

#compare height to area
list111 <- unique(top_heat_area$Accession)
top_heat_height$compare.geno <- top_heat_height$Accession %in% list111
all_overlap_heat_top <- filter(top_heat_height, compare.geno == TRUE)
length(unique(all_overlap_heat_top$Accession))
all_overlap_heat_top = dplyr::select(all_overlap_heat_top, -compare.geno)
list114 <- unique(all_overlap_heat_top$Accession)

#compare damage/area to height to find overall overlap
top_heat_height$compare.geno <- top_heat_height$Accession %in% list112
all_overlap_heat_top <- filter(top_heat_height, compare.geno == TRUE)


### See which accessions overlap that perform poorly ###

# compare damage to area
list121 <- unique(bottom_heat_area$Accession)
bottom_heat_damage$compare.geno <- bottom_heat_damage$Accession %in% list121
all_overlap_heat_bottom <- filter(bottom_heat_damage, compare.geno == TRUE)
length(unique(all_overlap_heat_bottom$Accession))
all_overlap_heat_bottom = dplyr::select(all_overlap_heat_bottom, -compare.geno)
list122 <- unique(all_overlap_heat_bottom$Accession)

#compare damage and height
list121 <- unique(bottom_heat_damage$Accession)
bottom_heat_height$compare.geno <- bottom_heat_height$Accession %in% list121
all_overlap_heat_bottom <- filter(bottom_heat_height, compare.geno == TRUE)
length(unique(all_overlap_heat_bottom$Accession))
all_overlap_heat_bottom = dplyr::select(all_overlap_heat_bottom, -compare.geno)
list123 <- unique(all_overlap_heat_bottom$Accession)

#compare area and height
list121 <- unique(bottom_heat_area$Accession)
bottom_heat_height$compare.geno <- bottom_heat_height$Accession %in% list121
all_overlap_heat_bottom <- filter(bottom_heat_height, compare.geno == TRUE)
length(unique(all_overlap_heat_bottom$Accession))
all_overlap_heat_bottom = dplyr::select(all_overlap_heat_bottom, -compare.geno)
list124 <- unique(all_overlap_heat_bottom$Accession)

#compare damage/area to height to find overall overlap
bottom_heat_height$compare.geno <- bottom_heat_height$Accession %in% list122
all_overlap_heat_bottom <- filter(bottom_heat_height, compare.geno == TRUE)




#############################################
### Overlap in Heat-drought ###
#############################################

#top_heatdrt_damage <- heatdrt_damage %>% top_n(34)
#bottom_heatdrt_damage <- heatdrt_damage %>% top_n(-34)

#top_heatdrt_area <- heatdrt_area %>% top_n(34)
#bottom_heatdrt_area <- heatdrt_area %>% top_n(-34)

#top_heatdrt_height <- heatdrt_height %>% top_n(34)
#bottom_heatdrt_height <- heatdrt_height %>% top_n(-34)


### See which accessions overlap that perform well ###

# compare damage to area 
list131 <- unique(top_heatdrt_area$Accession)
top_heatdrt_damage$compare.geno <- top_heatdrt_damage$Accession %in% list131
all_overlap_heatdrt_top <- filter(top_heatdrt_damage, compare.geno == TRUE)
length(unique(all_overlap_heatdrt_top$Accession))
all_overlap_heatdrt_top = dplyr::select(all_overlap_heatdrt_top, -compare.geno)
list132 <- unique(all_overlap_heatdrt_top$Accession)

#compare height to damage
list131 <- unique(top_heatdrt_damage$Accession)
top_heatdrt_height$compare.geno <- top_heatdrt_height$Accession %in% list131
all_overlap_heatdrt_top <- filter(top_heatdrt_height, compare.geno == TRUE)
length(unique(all_overlap_heatdrt_top$Accession))
all_overlap_heatdrt_top = dplyr::select(all_overlap_heatdrt_top, -compare.geno)
list133 <- unique(all_overlap_heatdrt_top$Accession)

#compare height to area
list131 <- unique(top_heatdrt_area$Accession)
top_heatdrt_height$compare.geno <- top_heatdrt_height$Accession %in% list131
all_overlap_heatdrt_top <- filter(top_heatdrt_height, compare.geno == TRUE)
length(unique(all_overlap_heatdrt_top$Accession))
all_overlap_heatdrt_top = dplyr::select(all_overlap_heatdrt_top, -compare.geno)
list134 <- unique(all_overlap_heatdrt_top$Accession)

#compare damage/area to height to find overall overlap
top_heatdrt_height$compare.geno <- top_heatdrt_height$Accession %in% list132
all_overlap_heatdrt_top <- filter(top_heatdrt_height, compare.geno == TRUE)


### See which accessions overlap that perform poorly ###

# compare damage to area
list141 <- unique(bottom_heatdrt_area$Accession)
bottom_heatdrt_damage$compare.geno <- bottom_heatdrt_damage$Accession %in% list141
all_overlap_heatdrt_bottom <- filter(bottom_heatdrt_damage, compare.geno == TRUE)
length(unique(all_overlap_heatdrt_bottom$Accession))
all_overlap_heatdrt_bottom = dplyr::select(all_overlap_heatdrt_bottom, -compare.geno)
list142 <- unique(all_overlap_heatdrt_bottom$Accession)

#compare damage and height
list141 <- unique(bottom_heatdrt_damage$Accession)
bottom_heatdrt_height$compare.geno <- bottom_heatdrt_height$Accession %in% list141
all_overlap_heatdrt_bottom <- filter(bottom_heatdrt_height, compare.geno == TRUE)
length(unique(all_overlap_heatdrt_bottom$Accession))
all_overlap_heatdrt_bottom = dplyr::select(all_overlap_heatdrt_bottom, -compare.geno)
list143 <- unique(all_overlap_heatdrt_bottom$Accession)

#compare area and height
list141 <- unique(bottom_heatdrt_area$Accession)
bottom_heatdrt_height$compare.geno <- bottom_heatdrt_height$Accession %in% list141
all_overlap_heatdrt_bottom <- filter(bottom_heatdrt_height, compare.geno == TRUE)
length(unique(all_overlap_heatdrt_bottom$Accession))
all_overlap_heatdrt_bottom = dplyr::select(all_overlap_heatdrt_bottom, -compare.geno)
list144 <- unique(all_overlap_heatdrt_bottom$Accession)

#compare damage/area to height to find overall overlap
bottom_heatdrt_height$compare.geno <- bottom_heatdrt_height$Accession %in% list142
all_overlap_heatdrt_bottom <- filter(bottom_heatdrt_height, compare.geno == TRUE)


