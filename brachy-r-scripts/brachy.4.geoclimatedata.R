rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gplots)
library(plotly)
library(FactoMineR)
library(factoextra)
library(missMDA)
library(corrplot)
library("Hmisc")
library("data.table")


#***************************************************
### First Time Start Here ### (Only have to do this once)
#***************************************************

#read in trait data (plantcv outputs)
setwd("/Users/eludwig/Desktop/brachy/")
DFs <- read.csv("brachy_control_VIS_SV_single_data_compiled_currated.csv")

heatDFs <- read.csv("brachy_heat_VIS_SV_single_data_compiled_currated.csv")



#compare Accessions in control and heat treatments and get list of all Accessions that are in either experiment
geno_control <- as.data.frame(unique(DFs$Accession))
colnames(geno_control)[1] <- "Accession"
geno_heat <- as.data.frame(unique(heatDFs$Accession))
colnames(geno_heat)[1] <- "Accession"
bothDF <- merge(geno_control, geno_heat, by.x = "Accession", by.y = "Accession", all = TRUE)
length(unique(bothDF$Accession))
in_either <- unique(bothDF$Accession)

DFs$compare.geno <- DFs$Accession %in% in_either
heatDFs$compare.geno <- heatDFs$Accession %in% in_either
DFs_149 <- filter(DFs, compare.geno == TRUE)
heatDFs_149 <- filter(heatDFs, compare.geno == TRUE)

#**********************************************
### geographical & climate data ###
#**********************************************

###############################################
### Get WorldClim Data ### (Only have to do this once)
###############################################

#download package from github
devtools::install_github("kapitzas/WorldClimTiles")
#require package (also need raster package)
library("WorldClimTiles")
library(raster)

#use country GADM boundaries to see which tiles we need to download from WorldClim (Spain and Turkey)
#look up country GADM codes to figure out what the ones you need are called
#Level 0 indicates country level boundaries, higher number is more specific (state, county, etc)
esp <- getData("GADM", country = "ESP", level = 0)
tilenames <- tile_name(esp)
tilenames

tur <- getData("GADM", country = "TUR", level = 0)
tilenames2 <- tile_name(tur)
tilenames2

#create list of tiles that you need and want to download and then merge
#can look up which tiles you need 
tilenames_merged <- c("15","16","17")
tilenames_merged <- as.character(tilenames_merged)

#download tiles from WorldClim
tiles <- tile_get(tilenames_merged, "bio")

#merge tiles (can be a LONG STEP depending on your computer's memory/RAM, takes a few mins per layer per tile so about 13 hrs total)
merged <- tile_merge(tiles)

#copy merged tiles so don't have to repeat merging step in case something goes wrong
merged_copy <- merged

#Plot one raster layer of merged tiles to check it
plot(merged[[1]], main="Annual Mean Temperature")

#save merged raster to disk so don't have to repeat previous steps
writeRaster(merged, filename = "mergedtiles.tif", format="GTiff", overwrite=TRUE)


###############################################
### Organize and Filter WorldClim Data ### (Only have to do this once)
###############################################

#read in merged raster
merged <- stack("/Users/eludwig/Desktop/brachy/mergedtiles.tif")

# Define the extent of the study region we want
e <- extent(c(-7, 45, 35, 47))
# crop worldclim data to this extent
climate <- crop(merged, e)

#plot one raster layer to check cropping
plot(climate[[1]], main="Annual Mean Temperature")

#Convert climate from rasterstack format to dataframe format
climate <- as.data.frame(climate, xy=TRUE)
write.csv(climate, "./worldclim_data.csv")


#read in locations csv file and keep relevant info (file with my accessions collection locations)
locations <- read.csv("./brachy_collection_locations.csv")
locations = dplyr::select(locations, Accession, Collection_location, Latitude, Longitude, Elevation)
locations = drop_na((locations))

#Split collection_location column to have country alone in column (some included state and city info that we don't need)
locations$Collection_location <- as.character(locations$Collection_location)
locations <- separate(data= locations, col=Collection_location, into= c("A", "B", "C", "Country"), fill = "left", sep = "\\, ")
locations <- dplyr::select(locations, Accession, Country, Latitude, Longitude, Elevation)

locations$compare.geno <- locations$Accession %in% in_either
locations_149 <- filter(locations, compare.geno == TRUE)

#save locations_149
write.csv(locations_149, "./brachy_locations_149.csv", row.names = FALSE)

#read in climate data
climate <- read.csv("worldclim_data.csv")

#### Combine climate data with location data

#Calculate the euclidian distance between the accession collection locations and WorldClim data points
back <- apply(t(cbind(locations$Latitude,locations$Longitude)),MARGIN = 2,function(i){
  climate$diff_x <- climate$x - i[2]
  climate$diff_y <- climate$y - i[1]
  climate$euclid <- with(climate,sqrt(diff_x^2 + diff_y^2))
  climate[which.min(climate$euclid),]
})

back_df <- data.frame(do.call(rbind,back),stringsAsFactors = F)
#rename bioclimatic variable column names
back_df <- back_df %>% 
  rename(
    bio1 = mergedtiles.1,
    bio2 = mergedtiles.2,
    bio3 = mergedtiles.3,
    bio4 = mergedtiles.4,
    bio5 = mergedtiles.5,
    bio6 = mergedtiles.6,
    bio7 = mergedtiles.7,
    bio8 = mergedtiles.8,
    bio9 = mergedtiles.9,
    bio10 = mergedtiles.10,
    bio11 = mergedtiles.11,
    bio12 = mergedtiles.12,
    bio13 = mergedtiles.13,
    bio14 = mergedtiles.14,
    bio15 = mergedtiles.15,
    bio16 = mergedtiles.16,
    bio17 = mergedtiles.17,
    bio18 = mergedtiles.18,
    bio19 = mergedtiles.19,
  )

#bind locations df with back_df and keep only values with climate data within euclidean distance <sqrt(2)
locations_new <- cbind(locations,back_df)
head(locations_new)
filtered <- locations_new[locations_new$euclid < sqrt(2),]

#make a histogram of the euclidian distance between collection locations and climate data points to check distribution
hist(back_df$euclid)

#compare Accessions in control and heat treatments and keep only the ones in either experiment
geno_control <- as.data.frame(unique(DFs$Accession))
colnames(geno_control)[1] <- "Accession"
geno_heat <- as.data.frame(unique(heatDFs$Accession))
colnames(geno_heat)[1] <- "Accession"
bothDF <- merge(geno_control, geno_heat, by.x = "Accession", by.y = "Accession", all = TRUE)
length(unique(bothDF$Accession))
in_either <- unique(bothDF$Accession)

locations_new$compare.geno <- locations_new$Accession %in% in_either
locations_new <- filter(locations_new, compare.geno == TRUE)
length(unique(locations_new$Accession))
locations_new = dplyr::select(locations_new, -compare.geno)

#remove columns we don't need
locations_new <- dplyr::select(locations_new, -c(X, x, y, diff_x, diff_y))

#save locations dataframe
write.csv(locations_new, "./brachy_final_locations.csv", row.names = FALSE)
write.csv(locations_new, "/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/brachy_final_locations.csv", row.names = FALSE)


#####

#**********************************************
### After first time start here ###
#**********************************************

setwd("/Users/eludwig/Desktop/brachy/")
locations_new <- read.csv("brachy_final_locations.csv")
locations_new <- dplyr::select(locations_new, -c(euclid))
locations_new$Elevation <- as.numeric(locations_new$Elevation)


###############################################
### Assess where certain accessions are in relation to the rest ###
###############################################
# subset those accessions to be able to color them differently than the rest
Bd210 <- subset(locations_new, Accession == "Bd21-0")
#BdTR9K <- subset(locations_new, Accession == "BdTR9K")
#BdTR10A <- subset(locations_new, Accession == "BdTR10A")
#ABR5 <- subset(locations_new, Accession == "ABR5")
#ABR6 <- subset(locations_new, Accession == "ABR6")
#Arc1 <- subset(locations_new, Accession == "Arc1")
#BdTR11A <- subset(locations_new, Accession == "BdTR11A")
#BdTR12B <- subset(locations_new, Accession == "BdTR12B")
#Mon3 <- subset(locations_new, Accession == "Mon3")
BdTR1E <- subset(locations_new, Accession == "BdTR1E")
BdTR10E <- subset(locations_new, Accession == "BdTR10E")
Adi18 <- subset(locations_new, Accession == "Adi-18")
BdTR5J <- subset(locations_new, Accession == "BdTR5J")



bio8 <- select(locations_new, Accession, Country, Latitude, Longitude, Elevation, bio8)
bio8 <- subset(locations_new, Accession == "BdTR10E" | Accession == "BdTR5J")
bio8$bio8 <- (bio8$bio8/10)
BdTR10E$bio8 <- (BdTR10E$bio8/10)
BdTR5J$bio8 <- (BdTR5J$bio8/10)

# plot a layer of the climate data
g <- ggplot(bio8, aes(Accession,bio8)) +
  geom_col() +
  geom_col(data=BdTR10E, fill="blue") +
  geom_col(data=BdTR5J, fill="orange") +
  ylab("Temperature (°C)") +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 16))+
  theme(title = element_text(size = 14))+
  ggtitle("Bio8 = Mean Temperature of Wettest Quarter") 
g
ggsave(plot = g,"/Users/eludwig/Desktop/bio8.png", width = 6, height = 4, dpi = 300)


  
#geom_point(data=BdTR9K, colour="purple")   
  #geom_point(data=Mon3, colour="pink")   


  


test <- filter(avg_height, imgdays == "15")
ggplot(test, aes(Accession,avg_height, color = Treatment)) +
  geom_point() 


test2 <- filter(avg_area, imgdays == "15")
ggplot(test2, aes(Accession,avg_area, color = Treatment)) +
  geom_point() 

test3 <- merge(test, test2, by = c("Accession", "Treatment"), all = FALSE)
test3$compactness <- test3$avg_area / test3$avg_height 

ABR5 <- subset(test3, Accession == "ABR5")
BdTR9K <- subset(test3, Accession == "BdTR9K")

ggplot(test3, aes(Accession,compactness, color = Treatment)) +
  geom_point() +
  geom_point(data=ABR5, colour="red") +  
  geom_point(data=BdTR9K, colour="purple")   


###############################################
### Climate Correlation Plot and PCA ###
###############################################


#Make correlation plot with climate variables

bio <-select(locations_new, -Accession,-Country,-Latitude,-Longitude)
bio_cor <- cor(bio)
pdf(file= "./climatecorrelationplot.pdf", width = 6, height = 7)
png(file= "./climatecorrelationplot.png", width = 1200, height = 1400)
corrplot(bio_cor, method = 'circle', type = 'lower', order = 'FPC',
         tl.col = "black", tl.srt = 45, tl.cex = 3, cl.cex = 3,
         title = "Correlation Plot of Climate Variables", mar=c(0,0,2,0))
dev.off()


# Do PCA with climate data
locations_new$Elevation <- as.numeric(locations_new$Elevation)
bio.pca = PCA(locations_new[,5:24], scale.unit=TRUE, ncp=5, graph=T)

climate_pca_df = data.frame(
  "Accession" = locations_new$Accession,
  "Country" = locations_new$Country,
  "Latitude" = locations_new$Latitude,
  "Longitude" = locations_new$Longitude,
  "PC1" = bio.pca$ind$coord[, 1],
  "PC2" = bio.pca$ind$coord[, 2],
  "PC3" = bio.pca$ind$coord[, 3],
  "PC4" = bio.pca$ind$coord[, 4])

#save climate_pca
write.csv(climate_pca_df, "/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/brachy_climate_pca_df.csv", row.names = FALSE)

percentofvariance = data.frame(bio.pca$eig) 
percentofvariance$PC <- row.names(percentofvariance)
row.names(percentofvariance) <- NULL
percentofvariance$PC <- row.names(percentofvariance)
percentofvariance$PC <- as.numeric(percentofvariance$PC)
percentofvariance$percentage.of.variance <- as.numeric(percentofvariance$percentage.of.variance)
str(percentofvariance)

#plot scree plot
scree <- ggplot(percentofvariance, aes(PC, percentage.of.variance, group=1))+
  labs(title = "Climate PCA Scree Plot", x = "Principal Component", y = "Percentage of Variance Explained")+
  #scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"))+
  scale_x_discrete(limit = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"))+
  geom_line()+
  geom_point()
scree
ggsave(plot = scree,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/climatescreeplot.png", width = 7.25, height = 6, dpi = 300)
ggsave(plot = scree,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/climatescreeplot.pdf", width = 7.25, height = 6, dpi = 300)



#calculate loadings for climate PCA
# loadings are significant when > sqrt(1/19) which = 0.229415733870562 (19 is # of climate variables in PCA)
climate.loadings <- sweep(bio.pca$var$coord,2,sqrt(bio.pca$eig[1:ncol(bio.pca$var$coord),1]),FUN="/")
write.csv(climate.loadings, "/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/climate_loadings.csv")

# make colorblind friendly palette
cbPalette <- c("#999999", "#0072B2", "#56B4E9", "#004D40", "#E69F00", "#F0E442", "#D55E00", "#636985")
# To use for fills, add
#scale_fill_manual(values=cbPalette)
# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)


climate_pca <- ggplot(climate_pca_df, aes(x = PC1, y = PC2, color = Country))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  geom_point()+
  theme_minimal()+
  theme(panel.grid = element_blank(), panel.border = element_rect(fill="transparent"))+
  scale_colour_manual(values=cbPalette)+
  #scale_fill_manual(values=cbPalette)+
  stat_ellipse()+
  xlab("PC 1 (47.92%)")+
  ylab("PC 2 (21.82%)")+
  ggtitle("PCA of Accession Collection Location Climate")
climate_pca
ggsave(plot = climate_pca,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/climatepca.png", width = 7.25, height = 6, dpi = 300)
ggsave(plot = climate_pca,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/climatepca.pdf", width = 7.25, height = 6, dpi = 300)






############################################
### Identify any associations between climate and phenotype (need to do this after gwas script "brachy.7.gapitgwas.R")
###############################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp4")
height_phenotype_tp4 <- read.table("./height_phenotype_tp4_GAPIT.txt", header = TRUE)

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp4")
area_phenotype_tp4 <- read.table("./area_phenotype_tp4_GAPIT.txt", header = TRUE)

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp4")
damage_phenotype_tp4 <- read.table("./damage_phenotype_tp4_GAPIT.txt", header = TRUE)

climate_pca_df <- read.csv("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/brachy_climate_pca_df.csv")
colnames(climate_pca_df)[1] <- "Accession"

locations_new <- read.csv("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/brachy_final_locations.csv")
colnames(locations_new)[1] <- "Accession"

climate_asso_height <- join(locations_new, height_phenotype_tp4, by="Accession")
climate_asso_height <- climate_asso_height %>% 
  rename(
    height_drought = Drought,
    height_heat = Heat,
    height_heatdrought = Heat.Drought)

climate_asso_area <- join(locations_new, area_phenotype_tp4, by="Accession")
climate_asso_area <- climate_asso_area %>% 
  rename(
    area_drought = Drought,
    area_heat = Heat,
    area_heatdrought = Heat.Drought)

climate_asso_damage <- join(locations_new, damage_phenotype_tp4, by="Accession")
climate_asso_damage <- climate_asso_damage %>% 
  rename(
    damage_drought = Drought,
    damage_heat = Heat,
    damage_heatdrought = Heat.Drought)


climate_p_height <-do.call(rbind,lapply(c("Elevation",paste0("bio",1:19)), function(col){
  do.call(rbind, lapply(c("height_drought", "height_heat", "height_heatdrought"), function(condition){
    p<-cor.test(climate_asso_height[[col]], climate_asso_height[[condition]], method="spearman")$p.value
    data.frame(cond=condition, col=col, p.value=p)
  }))
}))


climate_p_area <-do.call(rbind,lapply(c("Elevation",paste0("bio",1:19)), function(col){
  do.call(rbind, lapply(c("area_drought", "area_heat", "area_heatdrought"), function(condition){
    p<-cor.test(climate_asso_area[[col]], climate_asso_area[[condition]], method="spearman")$p.value
    data.frame(cond=condition, col=col, p.value=p)
  }))
}))


climate_p_damage <-do.call(rbind,lapply(c("Elevation",paste0("bio",1:19)), function(col){
  do.call(rbind, lapply(c("damage_drought", "damage_heat", "damage_heatdrought"), function(condition){
    p<-cor.test(climate_asso_damage[[col]], climate_asso_damage[[condition]], method="spearman")$p.value
    data.frame(cond=condition, col=col, p.value=p)
  }))
}))


climate_p_height$sig<-ifelse(climate_p_height$p.value<=0.001, "***",
                             ifelse(climate_p_height$p.value <=0.01, "**",
                                    ifelse(climate_p_height$p.value<=0.05, "*", "")))
#climate_p_height$col<-factor(climate_p_height$col, levels=c("Elevation", paste0("bio",19:1)), ordered=T)
#ggplot(climate_p_height, aes(y=col,x=cond,label=sig,fill=p.value))+ geom_tile()+geom_text(color="white")+ theme_minimal()


climate_p_area$sig<-ifelse(climate_p_area$p.value<=0.001, "***",
                             ifelse(climate_p_area$p.value <=0.01, "**",
                                    ifelse(climate_p_area$p.value<=0.05, "*", "")))
#climate_p_area$col<-factor(climate_p_area$col, levels=c("Elevation", paste0("bio",19:1)), ordered=T)
#ggplot(climate_p_area, aes(y=col,x=cond,label=sig,fill=p.value))+ geom_tile()+geom_text(color="white")+ theme_minimal()


climate_p_damage$sig<-ifelse(climate_p_damage$p.value<=0.001, "***",
                             ifelse(climate_p_damage$p.value <=0.01, "**",
                                    ifelse(climate_p_damage$p.value<=0.05, "*", "")))
#climate_p_damage$col<-factor(climate_p_damage$col, levels=c("Elevation", paste0("bio",19:1)), ordered=T)
#ggplot(climate_p_damage, aes(y=col,x=cond,label=sig,fill=p.value))+ geom_tile()+geom_text(color="white")+ theme_minimal()


climate_p_height$cond[grep("height_heatdrought", climate_p_height$cond)] <- "Height Heat-drought"
climate_p_height$cond[grep("height_drought", climate_p_height$cond)] <- "Height Drought"
climate_p_height$cond[grep("height_heat", climate_p_height$cond)] <- "Height Heat"


climate_p_area$cond[grep("area_heatdrought", climate_p_area$cond)] <- "Area Heat-drought"
climate_p_area$cond[grep("area_drought", climate_p_area$cond)] <- "Area Drought"
climate_p_area$cond[grep("area_heat", climate_p_area$cond)] <- "Area Heat"

climate_p_damage$cond[grep("damage_heatdrought", climate_p_damage$cond)] <- "Percent Damage Heat-drought"
climate_p_damage$cond[grep("damage_drought", climate_p_damage$cond)] <- "Percent Damage Drought"
climate_p_damage$cond[grep("damage_heat", climate_p_damage$cond)] <- "Percent Damage Heat"

climate_p <- rbind(climate_p_height, climate_p_area)
climate_p <- rbind(climate_p, climate_p_damage)


climate_p$col<-factor(climate_p$col, levels=c("Elevation", paste0("bio",19:1)), ordered=T)

climate_heatmap <- ggplot(climate_p, aes(y=col,x=cond,label=sig,fill=p.value))+
  geom_tile()+geom_text(color="white")+
  theme_minimal()+
  xlab("Trait and stress treatment")+
  ylab("Climate variable")+
  labs(fill = "p-value")+
  theme(
    axis.ticks = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1))
climate_heatmap
ggsave(plot = climate_heatmap,"/Users/eludwig/Downloads/climate_association_heatmap.png", width = 5, height = 8, dpi = 600)
ggsave(plot = climate_heatmap,"/Users/eludwig/Downloads/climate_association_heatmap.pdf", width = 5, height = 8)


############################################
### Identify any associations between climate and phenotype (need to do this after gwas script "brachy.7.gapitgwas.R")
###############################################
setwd("/Users/eludwig/Desktop/brachy/")
DFs <- read.csv("brachy_control_VIS_SV_single_data_compiled_currated.csv")

DFs$Treatment = NA
DFs$Treatment[grep("20", DFs$Treatment.2)] <- "Drought"
DFs$Treatment[grep("100", DFs$Treatment.2)] <- "Control"

heatDFs <- read.csv("brachy_heat_VIS_SV_single_data_compiled_currated.csv")

heatDFs$Treatment = NA
heatDFs$Treatment[grep("20", heatDFs$Treatment.2)] <- "Heat-drought"
heatDFs$Treatment[grep("100", heatDFs$Treatment.2)] <- "Heat"

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

avg_area_final <- avg_area_final %>% pivot_wider(names_from = Treatment, values_from = avg_area)
avg_height_final <- avg_height_final %>% pivot_wider(names_from = Treatment, values_from = avg_height)
avg_damage_final <- avg_damage_final %>% pivot_wider(names_from = Treatment, values_from = avg_percent_damage)

locations_new <- read.csv("/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/brachy_final_locations.csv")
colnames(locations_new)[1] <- "Accession"



climate_asso_height2 <- join(locations_new, avg_height_final, by="Accession")
climate_asso_height2 <- climate_asso_height2 %>% 
  rename(
    height_control = Control,
    height_drought = Drought,
    height_heat = Heat,
    height_heatdrought = "Heat-drought")

climate_asso_area2 <- join(locations_new, avg_area_final, by="Accession")
climate_asso_area2 <- climate_asso_area2 %>% 
  rename(
    area_control = Control,
    area_drought = Drought,
    area_heat = Heat,
    area_heatdrought = "Heat-drought")

climate_asso_damage2 <- join(locations_new, avg_damage_final, by="Accession")
climate_asso_damage2 <- climate_asso_damage2 %>% 
  rename(
    damage_control = Control,
    damage_drought = Drought,
    damage_heat = Heat,
    damage_heatdrought = "Heat-drought")


climate_p_height2 <-do.call(rbind,lapply(c("Elevation",paste0("bio",1:19)), function(col){
  do.call(rbind, lapply(c("height_control", "height_drought", "height_heat", "height_heatdrought"), function(condition){
    p<-cor.test(climate_asso_height2[[col]], climate_asso_height2[[condition]], method="spearman")$p.value
    data.frame(cond=condition, col=col, p.value=p)
  }))
}))


climate_p_area2 <-do.call(rbind,lapply(c("Elevation",paste0("bio",1:19)), function(col){
  do.call(rbind, lapply(c("area_control", "area_drought", "area_heat", "area_heatdrought"), function(condition){
    p<-cor.test(climate_asso_area2[[col]], climate_asso_area2[[condition]], method="spearman")$p.value
    data.frame(cond=condition, col=col, p.value=p)
  }))
}))


climate_p_damage2 <-do.call(rbind,lapply(c("Elevation",paste0("bio",1:19)), function(col){
  do.call(rbind, lapply(c("damage_control", "damage_drought", "damage_heat", "damage_heatdrought"), function(condition){
    p<-cor.test(climate_asso_damage2[[col]], climate_asso_damage2[[condition]], method="spearman")$p.value
    data.frame(cond=condition, col=col, p.value=p)
  }))
}))


climate_p_height2$sig<-ifelse(climate_p_height2$p.value<=0.001, "***",
                             ifelse(climate_p_height2$p.value <=0.01, "**",
                                    ifelse(climate_p_height2$p.value<=0.05, "*", "")))
#climate_p_height$col<-factor(climate_p_height$col, levels=c("Elevation", paste0("bio",19:1)), ordered=T)
#ggplot(climate_p_height, aes(y=col,x=cond,label=sig,fill=p.value))+ geom_tile()+geom_text(color="white")+ theme_minimal()


climate_p_area2$sig<-ifelse(climate_p_area2$p.value<=0.001, "***",
                           ifelse(climate_p_area2$p.value <=0.01, "**",
                                  ifelse(climate_p_area2$p.value<=0.05, "*", "")))
#climate_p_area$col<-factor(climate_p_area$col, levels=c("Elevation", paste0("bio",19:1)), ordered=T)
#ggplot(climate_p_area, aes(y=col,x=cond,label=sig,fill=p.value))+ geom_tile()+geom_text(color="white")+ theme_minimal()


climate_p_damage2$sig<-ifelse(climate_p_damage2$p.value<=0.001, "***",
                             ifelse(climate_p_damage2$p.value <=0.01, "**",
                                    ifelse(climate_p_damage2$p.value<=0.05, "*", "")))
#climate_p_damage$col<-factor(climate_p_damage$col, levels=c("Elevation", paste0("bio",19:1)), ordered=T)
#ggplot(climate_p_damage, aes(y=col,x=cond,label=sig,fill=p.value))+ geom_tile()+geom_text(color="white")+ theme_minimal()


climate_p_height2$cond[grep("height_heatdrought", climate_p_height2$cond)] <- "Height Heat-drought"
climate_p_height2$cond[grep("height_control", climate_p_height2$cond)] <- "Height Control"
climate_p_height2$cond[grep("height_drought", climate_p_height2$cond)] <- "Height Drought"
climate_p_height2$cond[grep("height_heat", climate_p_height2$cond)] <- "Height Heat"


climate_p_area2$cond[grep("area_heatdrought", climate_p_area2$cond)] <- "Area Heat-drought"
climate_p_area2$cond[grep("area_control", climate_p_area2$cond)] <- "Area Control"
climate_p_area2$cond[grep("area_drought", climate_p_area2$cond)] <- "Area Drought"
climate_p_area2$cond[grep("area_heat", climate_p_area2$cond)] <- "Area Heat"

climate_p_damage2$cond[grep("damage_heatdrought", climate_p_damage2$cond)] <- "Percent Damage Heat-drought"
climate_p_damage2$cond[grep("damage_control", climate_p_damage2$cond)] <- "Percent Damage Control"
climate_p_damage2$cond[grep("damage_drought", climate_p_damage2$cond)] <- "Percent Damage Drought"
climate_p_damage2$cond[grep("damage_heat", climate_p_damage2$cond)] <- "Percent Damage Heat"

climate_p2 <- rbind(climate_p_height2, climate_p_area2)
climate_p2 <- rbind(climate_p2, climate_p_damage2)


climate_p2$col<-factor(climate_p2$col, levels=c("Elevation", paste0("bio",19:1)), ordered=T)

climate_heatmap2 <- ggplot(climate_p2, aes(y=col,x=cond,label=sig,fill=p.value))+
  geom_tile()+geom_text(color="white")+
  theme_minimal()+
  ggtitle("Associations between climate variables and traits in experimental conditions",)+
  xlab("Trait and treatment")+
  ylab("Climate variable")+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
climate_heatmap2
ggsave(plot = climate_heatmap2,"/Users/eludwig/Downloads/climate_association_heatmap2.png", width = 6, height = 8, dpi = 600)
ggsave(plot = climate_heatmap2,"/Users/eludwig/Downloads/climate_association_heatmap2.pdf", width = 6, height = 8)



###############################################
# make map of population collection locations
###############################################


library(tidyverse)
library(sf)
library(raster)
library(mapview)
library(ggspatial)



mapview(locations_new, xcol = "Longitude", ycol = "Latitude", crs = 4269, grid = FALSE,legend = TRUE, col)

setwd("/Users/eludwig/Desktop/brachy/")
world_shape <- st_read("./ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
plot(world_shape, col="khaki")


locations_sf <- st_as_sf(locations_new, coords = c("Longitude", "Latitude"), crs = 4326)
locations_sf <- st_crop(locations_sf, extent(-10, 47, 30, 47))


map <- ggplot(data = world_shape) +
  geom_sf(fill="#e5f5e0")+
  annotation_scale(location = "bl", width_hint = 0.25) +
  #annotation_north_arrow(location = "bl", which_north = "true", 
  #pad_x = unit(1.20, "in"), pad_y = unit(0, "in"),
  #style = north_arrow_fancy_orienteering) +
  geom_sf(data = locations_sf, col="#FF6600")+
  coord_sf(xlim = c(-10,47), ylim = c(30,47),expand = TRUE)+
  xlab("Longitude")+
  ylab("Latitude")+
  ggtitle("B. distachyon Population Collection Locations")+
  theme(panel.grid.major = element_line(color = gray(.15), linetype = "dashed", size = 0.15), 
        panel.background = element_rect(fill = "#deebf7"))
map
ggsave(plot = map,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/collectionlocationmap2.png", width = 6, height = 4, dpi = 300)
ggsave(plot = map,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/collectionlocationmap2.pdf", width = 6, height = 4, dpi = 300)



devtools::install_github("tylermorganwall/rayshader")
library(rayshader)
library(magick)

#Here, I load a map with the raster package.
loadzip = tempfile() 
download.file("https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip", loadzip)
localtif = raster::raster(unzip(loadzip, "wc2.1_30s_elev.tif"))
unlink(loadzip)

e <- extent(c(-7, 45, 35, 47))
# crop worldclim data to this extent
localtif_cropped <- crop(localtif, e)

#And convert it to a matrix:
elmat = raster_to_matrix(localtif_cropped)

#We use another one of rayshader's built-in textures:
elmat %>%
  sphere_shade(sunangle = 45, texture = "imhof2") %>%
  plot_map()

elmat %>%
  sphere_shade(texture = "imhof2") %>%
  add_water(detect_water(elmat), color = "imhof3") %>%
  plot_map()

#“imhof1”, “imhof2”, “imhof3”, imhof4“,”desert“,”bw“,”unicorn".


elmat %>%
  sphere_shade(texture = "imhof2") %>%
  add_water(detect_water(elmat), color = "imhof3") %>%
  plot_3d(elmat, zscale = 12, fov = 0, theta = 0, zoom = .77, phi = 65, 
          windowsize = c(1000, 800), water = TRUE, waterdepth = 0, wateralpha = 0.5, watercolor = "lightblue",
          waterlinecolor = "white", waterlinealpha = 0.5)
Sys.sleep(0.2)
render_snapshot()

render_points(extent = attr(montereybay,"extent"), 
              lat = unlist(bird_track_lat), long = unlist(bird_track_long), 
              altitude = z_out, zscale=50, color="red")
render_highquality(point_radius = 1, samples = 256)



moss_landing_coord = c(36.806807, -121.793332) 
x_vel_out = -0.001 + rnorm(1000)[1:500]/1000
y_vel_out = rnorm(1000)[1:500]/200
z_out = c(seq(0,2000,length.out = 180), seq(2000,0,length.out=10), 
          seq(0,2000,length.out = 100), seq(2000,0,length.out=10))

bird_track_lat = list()
bird_track_long = list()
bird_track_lat[[1]] = moss_landing_coord[1]
bird_track_long[[1]] = moss_landing_coord[2]

for(i in 2:500) {
  bird_track_lat[[i]] = bird_track_lat[[i-1]] + y_vel_out[i]
  bird_track_long[[i]] = bird_track_long[[i-1]] + x_vel_out[i]
}

montereybay %>%
  sphere_shade(zscale = 10, texture = "imhof1") %>%
  plot_3d(montereybay, zscale = 50, fov = 70, theta = 270, phi = 30, 
          windowsize = c(1000, 800), zoom = 0.6,  
          water = TRUE, waterdepth = 0, wateralpha = 0.5, watercolor = "#233aa1",
          waterlinecolor = "white", waterlinealpha = 0.5)
Sys.sleep(0.2)
render_highquality(lightdirection = c(-45,45), lightaltitude  = 30, clamp_value = 10, 
                   samples = 256, camera_lookat= c(0,-50,0),
                   ground_material = diffuse(color="grey50",checkercolor = "grey20", checkerperiod = 100),
                   clear = TRUE)


