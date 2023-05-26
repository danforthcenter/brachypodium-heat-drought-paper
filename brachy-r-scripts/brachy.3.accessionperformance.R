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
library(gplots)
library(RColorBrewer)

setwd("/Users/eludwig/Desktop/brachy/")
DFs <- read.csv("brachy_control_VIS_SV_single_data_compiled_currated.csv")

DFs$Treatment = NA
DFs$Treatment[grep("20", DFs$Treatment.2)] <- "Drought"
DFs$Treatment[grep("100", DFs$Treatment.2)] <- "Control"

heatDFs <- read.csv("brachy_heat_VIS_SV_single_data_compiled_currated.csv")

heatDFs$Treatment = NA
heatDFs$Treatment[grep("20", heatDFs$Treatment.2)] <- "Heat-drought"
heatDFs$Treatment[grep("100", heatDFs$Treatment.2)] <- "Heat"



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
  ylab("PC 2 (13.42%)")+
  theme(legend.position=c(.50,-0.10))+
  guides(color=guide_legend(nrow=1))+
  ggtitle("PCA of Phenotypic Measurements over Time") +
  theme(plot.margin = unit(c(1, 1, 2.5, 1), "cm"))
pheno_pca_plot

#ggsave(plot = pheno_pca_plot,"/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/phenotype_pca.png", width = 7.5, height = 10, dpi = 300)
#ggsave(plot = pheno_pca_plot,"/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/phenotype_pca.pdf", width = 7.5, height = 10, dpi = 300)
ggsave(plot = pheno_pca_plot,"/Users/eludwig/Downloads/phenotype_pca.pdf", width = 7.5, height = 10, dpi = 300)


###########################################
### heatmap of performance ###
###########################################

# make avg_area, avg_height, avg_damage dfs in "Phenotype PCA" section above

# keep only data on last day
avg_area_final <- filter(avg_area, group == "44/45")
avg_area_final <- ungroup(avg_area_final)
avg_height_final <- filter(avg_height, group == "44/45")
avg_height_final <- ungroup(avg_height_final)
avg_damage_final <- filter(avg_damage, group == "44/45")
avg_damage_final <- ungroup(avg_damage_final)

# find accessions with data for all treatments
with_all <- avg_area_final %>% group_by(Accession) %>% summarise(number_of_treatments=n_distinct(Treatment))
with_all <- filter(with_all, number_of_treatments == 4)
with_all <- select(with_all, Accession)
all_list <- unique(with_all$Accession)

# keep only accessions with data for all treatments
avg_area_final$compare.geno <- avg_area_final$Accession %in% all_list
avg_area_final <- filter(avg_area_final, compare.geno == TRUE)
avg_area_final <- select(avg_area_final, -c(compare.geno, DAP, group))

avg_height_final$compare.geno <- avg_height_final$Accession %in% all_list
avg_height_final <- filter(avg_height_final, compare.geno == TRUE)
avg_height_final <- select(avg_height_final, -c(compare.geno, DAP, group))

avg_damage_final$compare.geno <- avg_damage_final$Accession %in% all_list
avg_damage_final <- filter(avg_damage_final, compare.geno == TRUE)
avg_damage_final <- select(avg_damage_final, -c(compare.geno, DAP, group))


### Make df with difference between control and stresses ###
# change data from long to wide format for heatmap and to find differences more easily
avg_area_final <- avg_area_final %>% pivot_wider(names_from = Treatment, values_from = avg_area)
avg_height_final <- avg_height_final %>% pivot_wider(names_from = Treatment, values_from = avg_height)
avg_damage_final <- avg_damage_final %>% pivot_wider(names_from = Treatment, values_from = avg_percent_damage)

# difference between area in stresses and control
area_lines_dif <- avg_area_final
area_lines_dif$drought_dif = NA
area_lines_dif$heat_dif = NA
area_lines_dif$heatdrought_dif = NA
area_lines_dif$drought_dif <- ((area_lines_dif$Drought-area_lines_dif$Control)/area_lines_dif$Control)
area_lines_dif$heat_dif <- ((area_lines_dif$Heat-area_lines_dif$Control)/area_lines_dif$Control)
area_lines_dif$heatdrought_dif <- ((area_lines_dif$`Heat-drought`-area_lines_dif$Control)/area_lines_dif$Control)
area_lines_dif <- dplyr::select(area_lines_dif, Accession, Control, drought_dif, heat_dif, heatdrought_dif)
area_lines_dif <- area_lines_dif %>% rename(Drought = drought_dif,Heat = heat_dif, `Heat-drought`= heatdrought_dif)

# difference between height in stresses and control
height_lines_dif <- avg_height_final
height_lines_dif$drought_dif = NA
height_lines_dif$heat_dif = NA
height_lines_dif$heatdrought_dif = NA
height_lines_dif$drought_dif <- ((height_lines_dif$Drought-height_lines_dif$Control)/height_lines_dif$Control)
height_lines_dif$heat_dif <- ((height_lines_dif$Heat-height_lines_dif$Control)/height_lines_dif$Control)
height_lines_dif$heatdrought_dif <- ((height_lines_dif$`Heat-drought`-height_lines_dif$Control)/height_lines_dif$Control)
height_lines_dif <- dplyr::select(height_lines_dif, Accession, Control, drought_dif, heat_dif, heatdrought_dif)
height_lines_dif <- height_lines_dif %>% rename(Drought = drought_dif,Heat = heat_dif, `Heat-drought`= heatdrought_dif)

# difference between percent damage in stresses and control
damage_lines_dif <- avg_damage_final
damage_lines_dif$drought_dif = NA
damage_lines_dif$heat_dif = NA
damage_lines_dif$heatdrought_dif = NA
damage_lines_dif$drought_dif <- (damage_lines_dif$Control-damage_lines_dif$Drought)
damage_lines_dif$heat_dif <- (damage_lines_dif$Control-damage_lines_dif$Heat)
damage_lines_dif$heatdrought_dif <- (damage_lines_dif$Control-damage_lines_dif$`Heat-drought`)
damage_lines_dif <- dplyr::select(damage_lines_dif, Accession, Control, drought_dif, heat_dif, heatdrought_dif)
damage_lines_dif <- damage_lines_dif %>% rename(Drought = drought_dif,Heat = heat_dif, `Heat-drought`= heatdrought_dif)

# get z-score for values for heatmaps
z_area <- as.data.frame(area_lines_dif)
z_area[,2:5] <- scale(z_area[,2:5])
z_height <- as.data.frame(height_lines_dif)
z_height[,2:5] <- scale(z_height[,2:5])
z_damage <- as.data.frame(damage_lines_dif)
z_damage$Control <- (z_damage$Control *-1)
z_damage[,2:5] <- scale(z_damage[,2:5])

# rename columns to include trait name so when we combine all three dfs the columns will have different names
colnames(z_area)[2] <- "area_control"
colnames(z_area)[3] <- "area_drought"
colnames(z_area)[4] <- "area_heat"
colnames(z_area)[5] <- "area_heatdrought"

colnames(z_height)[2] <- "height_control"
colnames(z_height)[3] <- "height_drought"
colnames(z_height)[4] <- "height_heat"
colnames(z_height)[5] <- "height_heatdrought"

colnames(z_damage)[2] <- "damage_control"
colnames(z_damage)[3] <- "damage_drought"
colnames(z_damage)[4] <- "damage_heat"
colnames(z_damage)[5] <- "damage_heatdrought"

# combine dfs for all 3 traits to find means of z-scores across traits and treatments
#then arrange by descending order for heatmap order
combined <- merge(z_area, z_height, by = "Accession")
combined <- merge(combined, z_damage, by = "Accession")
combined$Means <-apply(combined[,2:13],1,mean)
combined <- combined %>% arrange(desc(Means))

# separate traits to make it easier to plot them separately in heatmaps (only have data in df for that exact heatmap)
z_area <- select(combined, "Accession", "area_control", "area_drought", "area_heat", "area_heatdrought")
z_height <- select(combined, "Accession", "height_control", "height_drought", "height_heat", "height_heatdrought")
z_damage <- select(combined, "Accession", "damage_control", "damage_drought", "damage_heat", "damage_heatdrought")

# rename columns so heatmap axis labeling is easier (can just use column names)
colnames(z_area)[2] <- "Control"
colnames(z_area)[3] <- "Drought"
colnames(z_area)[4] <- "Heat"
colnames(z_area)[5] <- "Heat-drought"

colnames(z_height)[2] <- "Control"
colnames(z_height)[3] <- "Drought"
colnames(z_height)[4] <- "Heat"
colnames(z_height)[5] <- "Heat-drought"

colnames(z_damage)[2] <- "Control"
colnames(z_damage)[3] <- "Drought"
colnames(z_damage)[4] <- "Heat"
colnames(z_damage)[5] <- "Heat-drought"

# make distance matrix based off of one trait and then do hierarchical clustering (didn't use for figure)
# to use this in heatmap: set Rowv = hclust_avg, and dendrogram = "row"
#dist_mat <- dist(z_damage, method = 'euclidean')
#hclust_avg <- as.dendrogram(hclust(dist_mat, method = 'complete'))
#plot(hclust_avg)

#create the area heatmap
row.names_area <- z_area[1:129,1] # save the row names
row.names_area <- t(row.names_area)
row.names_area <- as.character(row.names_area)
area_data <- z_area[,2:5] # save the DF without the first column
plot.mtx_area <- as.matrix(area_data) #make the DF a matrix
rownames(plot.mtx_area) <- c(row.names_area) #add back row names

# make custom color palette
color.palette = colorRampPalette(c("lightyellow", "lightyellow", "lightblue", "midnightblue", "midnightblue"),space="rgb")
# save heatmap
pdf(file="/Users/eludwig/Downloads/area_performance_heatmap.pdf",width = 12,height = 18)
# plot the heatmap
heatmaparea <- heatmap.2(plot.mtx_area, 
                              Colv=FALSE,
                              Rowv = FALSE,
                              dendrogram = "none",
                              xlab = "Treatment",
                              ylab = "Accession",
                              main = "Plant Performance: Area",
                              col = color.palette,
                              #margin= c(5,5),
                              #scale = c("column"),
                              symbreaks = TRUE,
                              srtCol =0,
                              adjCol = c(0.5,.5),
                              #adjRow = c(0.1,0.5),
                              offsetRow = 0,
                              key.ylab = NA,
                              key.xlab = "performance",
                              keysize=2,
                              lhei=c(2, 15),
                              lwid = c(2,6),
                              density.info = "none",
                              trace = "none") # plot the heatmap
dev.off()



#create the height heatmap
row.names_height <- z_height[1:129,1] # save the row names
row.names_height <- t(row.names_height)
row.names_height <- as.character(row.names_height)
height_data <- z_height[,2:5] # save the DF without the first column
plot.mtx_height <- as.matrix(height_data) #make the DF a matrix
rownames(plot.mtx_height) <- c(row.names_height) #add back row names

# make custom color palette
color.palette = colorRampPalette(c("lightyellow", "lightyellow", "lightblue", "midnightblue", "midnightblue"),space="rgb")
# save heatmap
pdf(file="/Users/eludwig/Downloads/height_performance_heatmap.pdf",width = 12,height = 18)
# plot the heatmap
heatmapheight <- heatmap.2(plot.mtx_height, 
                         Colv=FALSE,
                         Rowv = FALSE,
                         dendrogram = "none",
                         xlab = "Treatment",
                         ylab = "Accession",
                         main = "Plant Performance: height",
                         col = color.palette,
                         #margin= c(5,5),
                         #scale = ("column"),
                         srtCol =0,
                         adjCol = c(0.5,.5),
                         #adjRow = c(0.1,0.5),
                         offsetRow = 0,
                         key.ylab = NA,
                         key.xlab = "performance",
                         keysize=2,
                         lhei=c(2, 15),
                         lwid = c(2,6),
                         density.info = "none",
                         trace = "none") # plot the heatmap
dev.off()


#create the percent damage heatmap
row.names_damage <- z_damage[1:129,1] # save the row names
row.names_damage <- t(row.names_damage)
row.names_damage <- as.character(row.names_damage)
damage_data <- z_damage[,2:5] # save the DF without the first column
plot.mtx_damage <- as.matrix(damage_data) #make the DF a matrix
rownames(plot.mtx_damage) <- c(row.names_damage) #add back row names

#make custom color palette
color.palette2 = colorRampPalette(c("lightyellow","lightyellow","lightyellow","lightyellow","lightyellow","lightyellow","lightblue","midnightblue","midnightblue","midnightblue","midnightblue", "midnightblue", "midnightblue"),space="rgb")
# save heatmap
pdf(file="/Users/eludwig/Downloads/damage_performance_heatmap.pdf",width = 12,height = 18)
# plot the heatmap
heatmapdamage <- heatmap.2(plot.mtx_damage, 
                         Colv=FALSE,
                         Rowv = FALSE,
                         dendrogram = "none",
                         xlab = "Treatment",
                         ylab = "Accession",
                         main = "Plant Performance: Damage",
                         col = color.palette2,
                         #margin= c(5,5),
                         #scale = ("column"),
                         srtCol =0,
                         adjCol = c(0.5,.5),
                         #adjRow = c(0.1,0.5),
                         offsetRow = 0,
                         key.ylab = NA,
                         key.xlab = "performance",
                         keysize=2,
                         lhei=c(2, 15),
                         lwid = c(2,6),
                         density.info = "none",
                         trace = "none") # plot the heatmap
dev.off()


###########################################
### Find accessions that perform well or poorly in and across stresses ###
###########################################
drought <- select(combined, "Accession", "area_drought", "height_drought", "damage_drought")
drought$Mean <-apply(drought[,2:4],1,mean)
drought <- drought %>% arrange(desc(Mean))
heat <- select(combined, "Accession", "area_heat", "height_heat", "damage_heat")
heat$Mean <-apply(heat[,2:4],1,mean)
heat <- heat %>% arrange(desc(Mean))
heatdrought <- select(combined, "Accession", "area_heatdrought", "height_heatdrought", "damage_heatdrought")
heatdrought$Mean <-apply(heatdrought[,2:4],1,mean)
heatdrought <- heatdrought %>% arrange(desc(Mean))


# make dfs for accessions ranked by performance in each treatment
colnames(area_lines_dif)[3] <- "area_drought_dif"
colnames(height_lines_dif)[3] <- "height_drought_dif"
colnames(damage_lines_dif)[3] <- "percent_damage_drought_dif"

colnames(area_lines_dif)[4] <- "area_heat_dif"
colnames(height_lines_dif)[4] <- "height_heat_dif"
colnames(damage_lines_dif)[4] <- "percent_damage_heat_dif"

colnames(area_lines_dif)[5] <- "area_heatdrought_dif"
colnames(height_lines_dif)[5] <- "height_heatdrought_dif"
colnames(damage_lines_dif)[5] <- "percent_damage_heatdrought_dif"

# df for accessions in drought including dif between drought and control and z-scores and mean
ranked_drought <- left_join(drought, select(area_lines_dif, -c("area_heat_dif", "area_heatdrought_dif", "Control")), by = "Accession")
ranked_drought <- left_join(ranked_drought, select(height_lines_dif, -c("height_heat_dif", "height_heatdrought_dif", "Control")), by = "Accession")
ranked_drought <- left_join(ranked_drought, select(damage_lines_dif, -c("percent_damage_heat_dif", "percent_damage_heatdrought_dif", "Control")), by = "Accession")
ranked_drought <- ranked_drought %>% 
  rename("z_area_dif" = "area_drought",
         "z_height_dif" = "height_drought",
         "z_damage_dif" = "damage_drought",
         "area_dif" = "area_drought_dif",
         "height_dif" = "height_drought_dif",
         "damage_dif" = "percent_damage_drought_dif",
         "mean_z" = "Mean")
ranked_drought <- ranked_drought[, c(1, 6, 7, 8, 2, 3, 4, 5)]
write.csv(ranked_drought, "./accessions_ranked_drought.csv")

# df for accessions in heat including dif between heat and control and z-scores and mean
ranked_heat <- left_join(heat, select(area_lines_dif, -c("area_drought_dif", "area_heatdrought_dif", "Control")), by = "Accession")
ranked_heat <- left_join(ranked_heat, select(height_lines_dif, -c("height_drought_dif", "height_heatdrought_dif", "Control")), by = "Accession")
ranked_heat <- left_join(ranked_heat, select(damage_lines_dif, -c("percent_damage_drought_dif", "percent_damage_heatdrought_dif", "Control")), by = "Accession")
ranked_heat <- ranked_heat %>% 
  rename("z_area_dif" = "area_heat",
         "z_height_dif" = "height_heat",
         "z_damage_dif" = "damage_heat",
         "area_dif" = "area_heat_dif",
         "height_dif" = "height_heat_dif",
         "damage_dif" = "percent_damage_heat_dif",
         "mean_z" = "Mean")
ranked_heat <- ranked_heat[, c(1, 6, 7, 8, 2, 3, 4, 5)]
write.csv(ranked_heat, "./accessions_ranked_heat.csv")

# df for accessions in heat-drought including dif between heat-drought and control and z-scores and mean
ranked_heatdrought <- left_join(heatdrought, select(area_lines_dif, -c("area_heat_dif", "area_drought_dif", "Control")), by = "Accession")
ranked_heatdrought <- left_join(ranked_heatdrought, select(height_lines_dif, -c("height_heat_dif", "height_drought_dif", "Control")), by = "Accession")
ranked_heatdrought <- left_join(ranked_heatdrought, select(damage_lines_dif, -c("percent_damage_heat_dif", "percent_damage_drought_dif", "Control")), by = "Accession")
ranked_heatdrought <- ranked_heatdrought %>% 
  rename("z_area_dif" = "area_heatdrought",
         "z_height_dif" = "height_heatdrought",
         "z_damage_dif" = "damage_heatdrought",
         "area_dif" = "area_heatdrought_dif",
         "height_dif" = "height_heatdrought_dif",
         "damage_dif" = "percent_damage_heatdrought_dif",
         "mean_z" = "Mean")
ranked_heatdrought <- ranked_heatdrought[, c(1, 6, 7, 8, 2, 3, 4, 5)]
write.csv(ranked_heatdrought, "./accessions_ranked_heatdrought.csv")




# make dfs for top and bottom 25% in each treatment
top_drought <- drop_na(dplyr::select(drought, Accession, Means)) %>% top_n(32) 
top_heat <- drop_na(dplyr::select(heat, Accession, Means)) %>% top_n(32) 
top_heatdrought <- drop_na(dplyr::select(heatdrought, Accession, Means)) %>% top_n(32) 

bottom_drought <- drop_na(dplyr::select(drought, Accession, Means)) %>% top_n(-32) 
bottom_heat <- drop_na(dplyr::select(heat, Accession, Means)) %>% top_n(-32) 
bottom_heatdrought <- drop_na(dplyr::select(heatdrought, Accession, Means)) %>% top_n(-32) 


###see which accessions overlap that perform well across stress treatments for all 3 traits
# compare drought to heat
list1 <- unique(top_drought$Accession)
top_heat$compare.geno <- top_heat$Accession %in% list1
top_overlap_drought_heat <- filter(top_heat, compare.geno == TRUE)
length(unique(top_overlap_drought_heat$Accession))
list2 <- unique(top_overlap_drought_heat$Accession)

#compare drought and heat-drought
top_heatdrought$compare.geno <- top_heatdrought$Accession %in% list1
top_overlap_drought_heatdrought <- filter(top_heatdrought, compare.geno == TRUE)
length(unique(top_overlap_drought_heatdrought$Accession))
list3 <- unique(top_overlap_drought_heat$Accession)

#compare heat and heat-drought
list4 <- unique(top_heat$Accession)
top_heatdrought$compare.geno <- top_heatdrought$Accession %in% list4
top_overlap_heat_heatdrought <- filter(top_heatdrought, compare.geno == TRUE)
length(unique(top_overlap_heat_heatdrought$Accession))
list5 <- unique(top_overlap_heat_heatdrought$Accession)

#compare heat/heat-drought to drought to find overall overlap
top_drought$compare.geno <- top_drought$Accession %in% list5
top_overlap_all <- filter(top_drought, compare.geno == TRUE)



###see which accessions overlap that perform poorly across stress treatments for all 3 traits
# compare drought to heat
list10 <- unique(bottom_drought$Accession)
bottom_heat$compare.geno <- bottom_heat$Accession %in% list10
bottom_overlap_drought_heat <- filter(bottom_heat, compare.geno == TRUE)
length(unique(bottom_overlap_drought_heat$Accession))
list20 <- unique(bottom_overlap_drought_heat$Accession)

#compare drought and heat-drought
bottom_heatdrought$compare.geno <- bottom_heatdrought$Accession %in% list10
bottom_overlap_drought_heatdrought <- filter(bottom_heatdrought, compare.geno == TRUE)
length(unique(bottom_overlap_drought_heatdrought$Accession))
list30 <- unique(bottom_overlap_drought_heat$Accession)

#compare heat and heat-drought
list40 <- unique(bottom_heat$Accession)
bottom_heatdrought$compare.geno <- bottom_heatdrought$Accession %in% list40
bottom_overlap_heat_heatdrought <- filter(bottom_heatdrought, compare.geno == TRUE)
length(unique(bottom_overlap_heat_heatdrought$Accession))
list50 <- unique(bottom_overlap_heat_heatdrought$Accession)

#compare heat/heat-drought to drought to find overall overlap
bottom_drought$compare.geno <- bottom_drought$Accession %in% list50
bottom_overlap_all <- filter(bottom_drought, compare.geno == TRUE)








###########################################
### Find means of traits in different treatments ###
###########################################

# keep only data on last day
avg_area_final <- filter(avg_area, group == "44/45")
avg_area_final <- ungroup(avg_area_final)
avg_height_final <- filter(avg_height, group == "44/45")
avg_height_final <- ungroup(avg_height_final)
avg_damage_final <- filter(avg_damage, group == "44/45")
avg_damage_final <- ungroup(avg_damage_final)

# filter data
avg_area_final <- select(avg_area_final, -c(DAP, group))
avg_height_final <- select(avg_height_final, -c(DAP, group))
avg_damage_final <- select(avg_damage_final, -c(DAP, group))


# change data from long to wide format to calculate column means
avg_area_final <- avg_area_final %>% pivot_wider(names_from = Treatment, values_from = avg_area)
avg_height_final <- avg_height_final %>% pivot_wider(names_from = Treatment, values_from = avg_height)
avg_damage_final <- avg_damage_final %>% pivot_wider(names_from = Treatment, values_from = avg_percent_damage)

# save data for these traits on final day as csv
write.csv(avg_area_final, "/Users/eludwig/Google Drive/My Drive/brachy_project/brachy_paper/supplements/avg_area_final_day.csv")
write.csv(avg_height_final, "/Users/eludwig/Google Drive/My Drive/brachy_project/brachy_paper/supplements/avg_height_final_day.csv")
write.csv(avg_damage_final, "/Users/eludwig/Google Drive/My Drive/brachy_project/brachy_paper/supplements/avg_percent_damage_final_day.csv")


#calculate means of columns and make df
means <- data.frame(matrix(ncol = 12, nrow = 1))
x <- c("height_control", "height_drought", "height_heat", "height_heatdrought","area_control", "area_drought", "area_heat", 
       "area_heatdrought", "damage_control", "damage_drought", "damage_heat", "damage_heatdrought")
colnames(means) <- x
means$height_control <- mean(avg_height_final$Control, na.rm = TRUE)
means$height_drought <- mean(avg_height_final$Drought, na.rm = TRUE)
means$height_heat <- mean(avg_height_final$Heat, na.rm = TRUE)
means$height_heatdrought <- mean(avg_height_final$`Heat-drought`, na.rm = TRUE)

means$area_control <- mean(avg_area_final$Control, na.rm = TRUE)
means$area_drought <- mean(avg_area_final$Drought, na.rm = TRUE)
means$area_heat <- mean(avg_area_final$Heat, na.rm = TRUE)
means$area_heatdrought <- mean(avg_area_final$`Heat-drought`, na.rm = TRUE)

means$damage_control <- mean(avg_damage_final$Control, na.rm = TRUE)
means$damage_drought <- mean(avg_damage_final$Drought, na.rm = TRUE)
means$damage_heat <- mean(avg_damage_final$Heat, na.rm = TRUE)
means$damage_heatdrought <- mean(avg_damage_final$`Heat-drought`, na.rm = TRUE)




# difference between area in stresses and control
area_lines_dif <- avg_area_final
area_lines_dif$drought_dif = NA
area_lines_dif$heat_dif = NA
area_lines_dif$heatdrought_dif = NA
area_lines_dif$drought_dif <- (area_lines_dif$Drought-area_lines_dif$Control)
area_lines_dif$heat_dif <- (area_lines_dif$Heat-area_lines_dif$Control)
area_lines_dif$heatdrought_dif <- (area_lines_dif$`Heat-drought`-area_lines_dif$Control)
area_lines_dif <- dplyr::select(area_lines_dif, Accession, drought_dif, heat_dif, heatdrought_dif)

# difference between height in stresses and control
height_lines_dif <- avg_height_final
height_lines_dif$drought_dif = NA
height_lines_dif$heat_dif = NA
height_lines_dif$heatdrought_dif = NA
height_lines_dif$drought_dif <- (height_lines_dif$Drought-height_lines_dif$Control)
height_lines_dif$heat_dif <- (height_lines_dif$Heat-height_lines_dif$Control)
height_lines_dif$heatdrought_dif <- (height_lines_dif$`Heat-drought`-height_lines_dif$Control)
height_lines_dif <- dplyr::select(height_lines_dif, Accession, drought_dif, heat_dif, heatdrought_dif)

# difference between percent damage in stresses and control
damage_lines_dif <- avg_damage_final
damage_lines_dif$drought_dif = NA
damage_lines_dif$heat_dif = NA
damage_lines_dif$heatdrought_dif = NA
damage_lines_dif$drought_dif <- (damage_lines_dif$Drought-damage_lines_dif$Control)
damage_lines_dif$heat_dif <- (damage_lines_dif$Heat-damage_lines_dif$Control)
damage_lines_dif$heatdrought_dif <- (damage_lines_dif$`Heat-drought`-damage_lines_dif$Control)
damage_lines_dif <- dplyr::select(damage_lines_dif, Accession, drought_dif, heat_dif, heatdrought_dif)

# save data for these traits on final day as csv
write.csv(area_lines_dif, "/Users/eludwig/Google Drive/My Drive/brachy_project/brachy_paper/supplements/avg_area_difference_between_control_and_stresses_final_day.csv")
write.csv(height_lines_dif, "/Users/eludwig/Google Drive/My Drive/brachy_project/brachy_paper/supplements/avg_height_difference_between_control_and_stresses_final_day.csv")
write.csv(damage_lines_dif, "/Users/eludwig/Google Drive/My Drive/brachy_project/brachy_paper/supplements/avg_percent_difference_between_control_and_stresses_damage_final_day.csv")

#calculate means of columns and make df
means_of_dif <- data.frame(matrix(ncol = 9, nrow = 1))
x <- c("height_drought", "height_heat", "height_heatdrought", "area_drought", "area_heat", "area_heatdrought",
       "damage_drought", "damage_heat", "damage_heatdrought")
colnames(means_of_dif) <- x
means_of_dif$height_drought <- mean(height_lines_dif$drought_dif, na.rm = TRUE)
means_of_dif$height_heat <- mean(height_lines_dif$heat_dif, na.rm = TRUE)
means_of_dif$height_heatdrought <- mean(height_lines_dif$heatdrought_dif, na.rm = TRUE)

means_of_dif$area_drought <- mean(area_lines_dif$drought_dif, na.rm = TRUE)
means_of_dif$area_heat <- mean(area_lines_dif$heat_dif, na.rm = TRUE)
means_of_dif$area_heatdrought <- mean(area_lines_dif$heatdrought_dif, na.rm = TRUE)

means_of_dif$damage_drought <- mean(damage_lines_dif$drought_dif, na.rm = TRUE)
means_of_dif$damage_heat <- mean(damage_lines_dif$heat_dif, na.rm = TRUE)
means_of_dif$damage_heatdrought <- mean(damage_lines_dif$heatdrought_dif, na.rm = TRUE)



