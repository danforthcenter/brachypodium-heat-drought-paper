rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
#library(readr)
#library(forcats)
#library(RColorBrewer)
library(gplots)
#library(plotly)
library(FactoMineR)
library(factoextra)
#library(missMDA)
#library(reshape2)
#library(lubridate)
#library(corrplot)
#library(car)


setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-drought/")
DFs <- read.csv("brachy_control_VIS_SV_single_data_compiled_currated.csv")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-heatdrought/")
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


############################################################
### PCA of Population Structure ###
############################################################
# Read in GBS structure data
gbs_pca <- read.table("brachypca.eigenvec", sep = " ", stringsAsFactors = FALSE)
# Change ARC1 to Arc1 to match locations Accession format
gbs_pca[23, 2] = "Arc1"
colnames(gbs_pca)[2] <- "Accession"

# See what PCs 1 and 2 for GBS look like
plot(gbs_pca[,3],gbs_pca[,4])

#keep only Accessions in gbs_pca that are in either experiment
gbs_pca$compare.geno <- gbs_pca$Accession %in% in_either
gbs_pca <- filter(gbs_pca, compare.geno == TRUE)
length(unique(gbs_pca$Accession))
gbs_pca = dplyr::select(gbs_pca, -compare.geno)

write.csv(gbs_pca, "./filtered_gbs_data.csv")


#Accessions_with_gbs <- as.data.frame(unique(gbs_pca$Accession))
#write.csv(Accessions_with_gbs, "./Accessions_with_gbs_data.csv")


# Check where in PCA the Accessions are that don't have all data to make sure not all clustered together
#ggplot(gbs_pca, aes(V3,V4))+
  #geom_point(aes(color=compare.geno))


# after first time start here
gbs_pca <- read.csv("./filtered_gbs_data.csv")
gbs_pca <- dplyr::select(gbs_pca, -c(X,V1))

# Read in climate pca data
climate_pca <- read.csv("brachy_climate_pca_df.csv")

# Join GBS PCA data with climate PCA data
joined_pca_df <- join(climate_pca, gbs_pca[1:5], by="Accession")

head(joined_pca_df)

# make colorblind friendly palette
cbPalette <- c("#999999", "#0072B2", "#56B4E9", "#004D40", "#E69F00", "#F0E442", "#D55E00", "#636985")
# To use for fills, add
#scale_fill_manual(values=cbPalette)
# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

# Plot population structure PCA
pop_pca <- ggplot(joined_pca_df, aes(x = V3, y = V4, color = Country))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  theme_minimal()+
  theme(panel.grid = element_blank(), panel.border = element_rect(fill="transparent"))+
  scale_colour_manual(values=cbPalette)+
  scale_fill_manual(values=cbPalette)+
  stat_ellipse()+
  xlab("PC 1 (18.76%)")+
  ylab("PC 2 (17.56%)")+
  ggtitle("PCA of Population Structure")+
  #geom_text(label=joined_pca_df$Accession, nudge_x = .01,check_overlap = T)+
  geom_point()
pop_pca
ggsave(plot = pop_pca,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/poppca.png", width = 7.25, height = 6, dpi = 300)
ggsave(plot = pop_pca,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/poppca.pdf", width = 7.25, height = 6, dpi = 300)

ggplotly(pop_pca)
