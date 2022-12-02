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
library(raster)
library("WorldClimTiles") # first need to download from github on line 56
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
require("WorldClimTiles")

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
merged <- stack("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-heatdrought/mergedtiles.tif")

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

### After first time start here ###

###############################################
### Climate Correlation Plot and PCA ###
###############################################


setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-heatdrought/")
locations_new <- read.csv("brachy_final_locations.csv")
locations_new <- dplyr::select(locations_new, -c(euclid))
locations_new$Elevation <- as.numeric(locations_new$Elevation)

turkey <- filter(locations_new, Country == "Turkey" | Country == "Iraq")
spain <- filter(locations_new, Country == "Spain" | Country == "France" | Country == "Slovenia")


###############################################
### Assess where certain accessions are in relation to the rest ###
###############################################
# subset those accessions to be able to color them differently than the rest
#Bd210 <- subset(locations_new, Accession == "Bd21-0")
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





# plot a layer of the climate data
ggplot(locations_new, aes(Accession,bio19)) +
  geom_point() +
  geom_point(data=BdTR1E, colour="green") + 
  geom_point(data=BdTR10E, colour="blue") +
  geom_point(data=Adi18, colour="red") +  
  geom_point(data=BdTR5J, colour="orange")  
  
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


#####

#Make correlation plot with climate variables

bio <-select(locations_new, -Accession,-Country,-Latitude,-Longitude)
bio_cor <- cor(bio)
pdf(file= "./climatecorrelationplot.pdf", width = 6, height = 7)
png(file= "./climatecorrelationplot.png", width = 1200, height = 1400)
corrplot(bio_cor, method = 'circle', type = 'lower', order = 'FPC',
         tl.col = "black", tl.srt = 45, tl.cex = 3, cl.cex = 3,
         title = "Correlation Plot of Climate Variables", mar=c(0,0,2,0))

#climatecorr
dev.off()


# Do PCA with climate data
locations_new$Elevation <- as.numeric(locations_new$Elevation)
#locations_new_imputed <- as.data.frame(imputePCA(locations_new[,5:24], ncp = 5)[[1]]) #impute missing values
#bio.pca = PCA(locations_new_imputed, scale.unit=TRUE, ncp=5, graph=T)
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

climate_asso <- join(climate_pca_df, damage_phenotype_tp4, by="Accession")

climate_asso <- climate_asso[,c("Accession", "PC1", "PC2", "PC3", "PC4", "Drought", "Heat", "Heat.Drought")]
head(climate_asso)
climate_asso <- climate_asso[!duplicated(climate_asso),]
with(climate_asso, cor.test(PC1,Drought,method = "spearman")) 
with(climate_asso, cor.test(PC1,Heat,method = "spearman"))
with(climate_asso, cor.test(PC1,Heat.Drought,method = "spearman")) 
with(climate_asso, cor.test(PC2,Drought,method = "spearman")) 
with(climate_asso, cor.test(PC2,Heat,method = "spearman"))
with(climate_asso, cor.test(PC2,Heat.Drought,method = "spearman")) 
with(climate_asso, cor.test(PC3,Drought,method = "spearman")) 
with(climate_asso, cor.test(PC3,Heat,method = "spearman"))
with(climate_asso, cor.test(PC3,Heat.Drought,method = "spearman")) 

# with Elevation
#PC1,Drought
#PC1,Heat
#PC1,Heat.Drought

#PC2,Drought
#PC2,Heat
#PC2,Heat.Drought, damage = 0.002062

#PC3,Drought, height = 0.04116 
#PC3,Heat, height = 0.002478
#PC3,Heat.Drought, height = 0.009578




###############################################
# make map of population collection locations
###############################################


library(tidyverse)
library(sf)
library(mapview)
library(ggspatial)


mapview(locations_new, xcol = "Longitude", ycol = "Latitude", crs = 4269, grid = FALSE,legend = TRUE, col)

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-heatdrought/")
world_shape <- st_read("./ne_10m_admin_0_countries/ne_10m_admin_0_countries.shp")
#plot(world_shape, col="khaki")

#e <- extent(c(-7, 45, 35, 47))
## Crop to the desired extent, then plot
world_shape <- st_crop(world_shape, extent(-10, 47, 30, 47))
plot(world_shape[4])

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




