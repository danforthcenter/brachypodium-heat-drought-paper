rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(gplots)
library(reshape2)
#install.packages("devtools")
#devtools::install_github("jiabowang/GAPIT3",force=TRUE)
#library(GAPIT3)

setwd("/Users/eludwig/Desktop/brachy/")
DFs <- read.csv("brachy_control_VIS_SV_single_data_compiled_currated.csv")

heatDFs <- read.csv("brachy_heat_VIS_SV_single_data_compiled_currated.csv")

accessions_with_gbs <- read.csv("accessions_with_gbs_data.csv")
colnames(accessions_with_gbs)[2] <- "Accession"
with_gbs <- unique(accessions_with_gbs$Accession)

# keep only accessions with GBS data
DFs$compare.geno <- DFs$Accession %in% with_gbs
heatDFs$compare.geno <- heatDFs$Accession %in% with_gbs
DFs_144 <- filter(DFs, compare.geno == TRUE)
heatDFs_144 <- filter(heatDFs, compare.geno == TRUE)


DFs_144 <- filter(DFs_144, DAP != 46)
DFs_144 <- filter(DFs_144, DAP != 47)

DFs_144$group = NA
DFs_144$group[grep("16", DFs_144$DAP)] <- "16"
DFs_144$group[grep("18", DFs_144$DAP)] <- "18"
DFs_144$group[grep("20", DFs_144$DAP)] <- "20"
DFs_144$group[grep("22", DFs_144$DAP)] <- "22"
DFs_144$group[grep("24", DFs_144$DAP)] <- "24"
DFs_144$group[grep("26", DFs_144$DAP)] <- "26/27"
DFs_144$group[grep("28", DFs_144$DAP)] <- "28/29"
DFs_144$group[grep("30", DFs_144$DAP)] <- "30/31"
DFs_144$group[grep("32", DFs_144$DAP)] <- "32/33"
DFs_144$group[grep("34", DFs_144$DAP)] <- "36/37"
DFs_144$group[grep("36", DFs_144$DAP)] <- "36/37"
DFs_144$group[grep("38", DFs_144$DAP)] <- "38/39"
DFs_144$group[grep("40", DFs_144$DAP)] <- "40/41"
DFs_144$group[grep("42", DFs_144$DAP)] <- "42/43"
DFs_144$group[grep("44", DFs_144$DAP)] <- "44/45"


heatDFs_144 <- filter(heatDFs_144, DAP != 47)
heatDFs_144 <- filter(heatDFs_144, DAP != 49)

heatDFs_144$group = NA
heatDFs_144$group[grep("16", heatDFs_144$DAP)] <- "16"
heatDFs_144$group[grep("18", heatDFs_144$DAP)] <- "18"
heatDFs_144$group[grep("20", heatDFs_144$DAP)] <- "20"
heatDFs_144$group[grep("22", heatDFs_144$DAP)] <- "22"
heatDFs_144$group[grep("24", heatDFs_144$DAP)] <- "24"
heatDFs_144$group[grep("26", heatDFs_144$DAP)] <- "26/27"
heatDFs_144$group[grep("27", heatDFs_144$DAP)] <- "26/27"
heatDFs_144$group[grep("29", heatDFs_144$DAP)] <- "28/29"
heatDFs_144$group[grep("31", heatDFs_144$DAP)] <- "30/31"
heatDFs_144$group[grep("33", heatDFs_144$DAP)] <- "32/33"
heatDFs_144$group[grep("35", heatDFs_144$DAP)] <- "36/37"
heatDFs_144$group[grep("37", heatDFs_144$DAP)] <- "36/37"
heatDFs_144$group[grep("39", heatDFs_144$DAP)] <- "38/39"
heatDFs_144$group[grep("41", heatDFs_144$DAP)] <- "40/41"
heatDFs_144$group[grep("43", heatDFs_144$DAP)] <- "42/43"
heatDFs_144$group[grep("45", heatDFs_144$DAP)] <- "44/45"




################################################################
### Read in hapmap and get in correct format ### 
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit")
# read in hapmap file
hapmap_og <- read.csv("./brachy_hapmap.csv")
#rename Accession names to match other files
colnames(hapmap_og)[21] <- "Adi-13"
colnames(hapmap_og)[22] <- "Adi-14"
colnames(hapmap_og)[23] <- "Adi-15"
colnames(hapmap_og)[24] <- "Adi-17"
colnames(hapmap_og)[25] <- "Adi-18"
colnames(hapmap_og)[26] <- "Adi-1"
colnames(hapmap_og)[27] <- "Adi-2"
colnames(hapmap_og)[28] <- "Adi-3"
colnames(hapmap_og)[29] <- "Adi-4"
colnames(hapmap_og)[30] <- "Adi-5"
colnames(hapmap_og)[31] <- "Adi-7"
colnames(hapmap_og)[32] <- "Adi-8"
colnames(hapmap_og)[33] <- "Adi-9"
colnames(hapmap_og)[34] <- "Arc1"
colnames(hapmap_og)[37] <- "Bd1-1"
colnames(hapmap_og)[38] <- "Bd21-0"
colnames(hapmap_og)[39] <- "Bd21-3"
colnames(hapmap_og)[134] <- "Bis-1"
colnames(hapmap_og)[136] <- "Gaz-2"
colnames(hapmap_og)[137] <- "Gaz-3"
colnames(hapmap_og)[138] <- "Gaz-4"
colnames(hapmap_og)[139] <- "Gaz-5"
colnames(hapmap_og)[140] <- "Gaz-7"
colnames(hapmap_og)[141] <- "Gaz-8"
colnames(hapmap_og)[143] <- "Kah-2"
colnames(hapmap_og)[144] <- "Kah-3"
colnames(hapmap_og)[145] <- "Kah-6"
colnames(hapmap_og)[146] <- "Koz-2"
colnames(hapmap_og)[147] <- "Koz-3"
colnames(hapmap_og)[154] <- "Tek-4"
colnames(hapmap_og)[155] <- "Tek-5"

# rename columns to remove any random symbols
colnames(hapmap_og)[1] <- "rs"
colnames(hapmap_og)[6] <- "assembly"
#fix snp column format
hapmap_og <- separate(hapmap_og, rs, c("rs1","rs2"))
hapmap_og <- unite(hapmap_og, rs, c(rs1, rs2))

#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only 144 accessions in both experiments
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% with_gbs
hapmap_og_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_og_melt$Accession))
hapmap_og_melt = dplyr::select(hapmap_og_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap_og <- reshape2::dcast(hapmap_og_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap_og <- dplyr::select(hapmap_og, -rn)

#remove centromere snp
hapmap_og <- hapmap_og[!(hapmap_og$chrom=="Bd1_centromere_containing_Bradi1g41430"),]

# Rename chromosomes to numbers
hapmap_og[hapmap_og=="BD1"]<-"1"
hapmap_og[hapmap_og=="BD2"]<-"2"
hapmap_og[hapmap_og=="BD3"]<-"3"
hapmap_og[hapmap_og=="BD4"]<-"4"
hapmap_og[hapmap_og=="BD5"]<-"5"

# save hapmap_og as tsv for GAPIT input
#write_tsv(hapmap_og, "./hapmap_GAPIT.tsv")
#write.csv(hapmap_og, "./hapmap_GAPIT.csv")

################################################################
### make file with Population Structure PCs ###
################################################################
# Only necessary on first run to see if necessary to include as covariates
# Didn't end up using them in final run

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-heatdrought/")
gbs_pca <- read.csv("filtered_gbs_data.csv")
pop_structure <- dplyr::select(gbs_pca, Accession, V3, V4, V5, V6)
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit")
write.table(pop_structure, "./popstructurePC_GAPIT.txt")

#read in pop structure file
pop_structure <- read.table("./popstructurePC_GAPIT.txt", header = TRUE)
pop_structure$V3 <- as.numeric(pop_structure$V3)
pop_structure$V4 <- as.numeric(pop_structure$V4)
pop_structure$V5 <- as.numeric(pop_structure$V5)
pop_structure$V6 <- as.numeric(pop_structure$V6)

################################################################
### make file with climate PCs
################################################################
# Only necessary on first run to see if necessary to include as covariates
# Didn't end up using them in final run

# Read in climate pca data
climate_pca <- read.csv("brachy_climate_pca_df.csv")
climate_cov <- dplyr::select(climate_pca, Accession, PC1, PC2, PC3, PC4)
write.table(climate_cov, "./climatePC_GAPIT.txt")

str(climate_cov)

#read in climate covariate file
climate_cov <- read.table("./climatePC_GAPIT.txt", header = TRUE)
climate_cov$PC1 <- as.numeric(climate_cov$PC1)
climate_cov$PC2 <- as.numeric(climate_cov$PC2)
climate_cov$PC3 <- as.numeric(climate_cov$PC3)
climate_cov$PC4 <- as.numeric(climate_cov$PC4)



#####



# ***TP1*** Comparison GWAS with GAPIT (organized by trait)
################################################################
### GAPIT GWAS for Height (difference)
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp1")

# Height Early Day: control and drought
####
avgheight <- group_by(DFs_144, Accession, Treatment, group) %>% summarise(avg_height = mean(cor_height_above_reference, na.rm = TRUE))

avgheight_tp1 <- aggregate(data = avgheight[avgheight$group == "22",], avg_height~Accession+Treatment, FUN= function(i)mean(i))
head(avgheight_tp1)

# check to make sure each phenotype and treatment has data
table(avgheight_tp1$Accession, avgheight_tp1$Treatment)
avgheight_tp1[avgheight_tp1 == "WW"] <- "control"
avgheight_tp1[avgheight_tp1 == "WL"] <- "drought"

#make df with only WW data from drought exp
droughtWW_avgheight_tp1 <- avgheight_tp1[avgheight_tp1$Treatment == "control",]

#make df with only WL data from drought experiment (control temp)
droughtWL_avgheight_tp1 <- avgheight_tp1[avgheight_tp1$Treatment == "drought",]

height_dif4 <- setNames(data.frame(do.call("rbind",lapply(split(avgheight_tp1, avgheight_tp1$Accession), function(i)diff(i$avg_height)))), c("height_dif"))
height_dif4$Accession <- row.names(height_dif4)
head(height_dif4)
row.names(height_dif4) <- NULL


####
#control and heat
####
heatDFs_144$cor_height_above_reference <- as.numeric(heatDFs_144$cor_height_above_reference)
avgheightheat <- group_by(heatDFs_144, Accession, Treatment, group) %>% summarise(avg_height = mean(cor_height_above_reference, na.rm = TRUE))
 

heatavgheight_tp1 <- aggregate(data = avgheightheat[avgheightheat$group == "22",], avg_height~Accession+Treatment, FUN= function(i)mean(i))
head(heatavgheight_tp1)

table(heatavgheight_tp1$Accession, heatavgheight_tp1$Treatment)
heatavgheight_tp1[heatavgheight_tp1 == "WW"] <- "heat"
heatavgheight_tp1[heatavgheight_tp1 == "WL"] <- "heat.drought"

# make df with only WW data from heat and control temps
heatWW_avgheight_tp1 <- heatavgheight_tp1[heatavgheight_tp1$Treatment == "heat",]

WW_avgheight_tp1 <- rbind(heatWW_avgheight_tp1, droughtWW_avgheight_tp1)

#check to make sure all accessions have data for both treatments on last day
table(WW_avgheight_tp1$Accession, WW_avgheight_tp1$Treatment)

height_dif5 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgheight_tp1, WW_avgheight_tp1$Accession), function(i)diff(i$avg_height)))), c("height_dif2"))
height_dif5$Accession <- row.names(height_dif5)
head(height_dif5)
row.names(height_dif5) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgheight_tp1 <- heatavgheight_tp1[heatavgheight_tp1$Treatment == "heat.drought",]

control_and_heatWL_avgheight_tp1 <- rbind(heatWL_avgheight_tp1, droughtWW_avgheight_tp1)

#check to make sure all accessions have data for both treatments on last day
table(control_and_heatWL_avgheight_tp1$Accession, control_and_heatWL_avgheight_tp1$Treatment)

height_dif6 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgheight_tp1, control_and_heatWL_avgheight_tp1$Accession), function(i)diff(i$avg_height)))), c("height_dif3"))
height_dif6$Accession <- row.names(height_dif6)
head(height_dif6)
row.names(height_dif6) <- NULL

# create phenotype (Y) file for height
temporary <- join(height_dif4, height_dif5, by="Accession")
height_phenotype_tp1 <- join(temporary, height_dif6, by="Accession")
height_phenotype_tp1 <- drop_na(height_phenotype_tp1)

#organize it
height_phenotype_tp1 <- dplyr::select(height_phenotype_tp1, Accession, height_dif, height_dif2, height_dif3)
# rename columns
colnames(height_phenotype_tp1)[2] <- "Drought"
colnames(height_phenotype_tp1)[3] <- "Heat"
colnames(height_phenotype_tp1)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(height_phenotype_tp1, "./height_phenotype_tp1_GAPIT.txt")


#### After making phenotype file once, can start here:
#read in phenotype file
height_phenotype_tp1 <- read.table("./height_phenotype_tp1_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(height_phenotype_tp1)

#check phenotype file
hist(height_phenotype_tp1$Drought)
mean(height_phenotype_tp1$Drought)
range(height_phenotype_tp1$Drought)
sd(height_phenotype_tp1$Drought)
which(is.na(height_phenotype_tp1$Drought))

hist(height_phenotype_tp1$Heat)
mean(height_phenotype_tp1$Heat)
range(height_phenotype_tp1$Heat)
sd(height_phenotype_tp1$Heat)
which(is.na(height_phenotype_tp1$Heat))

hist(height_phenotype_tp1$Heat.Drought)
mean(height_phenotype_tp1$Heat.Drought)
range(height_phenotype_tp1$Heat.Drought)
sd(height_phenotype_tp1$Heat.Drought)
which(is.na(height_phenotype_tp1$Heat.Drought))

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(height_phenotype_tp1$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap_og as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp1")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp1")
#run GAPIT
height_gapit <- GAPIT(
  Y = height_phenotype_tp1,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU"),
  Multiple_analysis = TRUE
)



################################################################
### GAPIT GWAS for Area (difference)
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp1")

# Area Early Day: control and drought
####
avgarea <- group_by(DFs_144, Accession, Treatment, group) %>% summarise(avg_area = mean(cor_area, na.rm = TRUE))

avgarea_tp1 <- aggregate(data = avgarea[avgarea$group == "22",], avg_area~Accession+Treatment, FUN= function(i)mean(i))
head(avgarea_tp1)

# check to make sure each phenotype and treatment has data
table(avgarea_tp1$Accession, avgarea_tp1$Treatment)
avgarea_tp1[avgarea_tp1 == "WW"] <- "control"
avgarea_tp1[avgarea_tp1 == "WL"] <- "drought"

#make df with only WL data from drought experiment (control temp)
droughtWL_avgarea_tp1 <- avgarea_tp1[avgarea_tp1$Treatment == "drought",]

area_dif4 <- setNames(data.frame(do.call("rbind",lapply(split(avgarea_tp1, avgarea_tp1$Accession), function(i)diff(i$avg_area)))), c("area_dif"))
area_dif4$Accession <- row.names(area_dif4)
head(area_dif4)
row.names(area_dif4) <- NULL


####
#control and heat
####
heatDFs_144$cor_area <- as.numeric(heatDFs_144$cor_area)
avgareaheat <- group_by(heatDFs_144, Accession, Treatment, group) %>% summarise(avg_area = mean(cor_area, na.rm = TRUE))
 

heatavgarea_tp1 <- aggregate(data = avgareaheat[avgareaheat$group == "22",], avg_area~Accession+Treatment, FUN= function(i)mean(i))
head(heatavgarea_tp1)

table(heatavgarea_tp1$Accession, heatavgarea_tp1$Treatment)
heatavgarea_tp1[heatavgarea_tp1 == "WW"] <- "heat"
heatavgarea_tp1[heatavgarea_tp1 == "WL"] <- "heat.drought"


# make df with only WW data from heat and control temps
heatWW_avgarea_tp1 <- heatavgarea_tp1[heatavgarea_tp1$Treatment == "heat",]

droughtWW_avgarea_tp1 <- avgarea_tp1[avgarea_tp1$Treatment == "control",]

WW_avgarea_tp1 <- rbind(heatWW_avgarea_tp1, droughtWW_avgarea_tp1)

#check to make sure all accessions have data for both treatments on last day
table(WW_avgarea_tp1$Accession, WW_avgarea_tp1$Treatment)

area_dif5 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgarea_tp1, WW_avgarea_tp1$Accession), function(i)diff(i$avg_area)))), c("area_dif2"))
area_dif5$Accession <- row.names(area_dif5)
head(area_dif5)
row.names(area_dif5) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgarea_tp1 <- heatavgarea_tp1[heatavgarea_tp1$Treatment == "heat.drought",]

control_and_heatWL_avgarea_tp1 <- rbind(heatWL_avgarea_tp1, droughtWW_avgarea_tp1)

#check to make sure all accessions have data for both treatments on last day
table(control_and_heatWL_avgarea_tp1$Accession, control_and_heatWL_avgarea_tp1$Treatment)

area_dif6 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgarea_tp1, control_and_heatWL_avgarea_tp1$Accession), function(i)diff(i$avg_area)))), c("area_dif3"))
area_dif6$Accession <- row.names(area_dif6)
head(area_dif6)
row.names(area_dif6) <- NULL

# create phenotype (Y) file for area
temporary <- join(area_dif4, area_dif5, by="Accession")
area_phenotype_tp1 <- join(temporary, area_dif6, by="Accession")
area_phenotype_tp1 <- drop_na(area_phenotype_tp1)

#organize it
area_phenotype_tp1 <- dplyr::select(area_phenotype_tp1, Accession, area_dif, area_dif2, area_dif3)
# rename columns
colnames(area_phenotype_tp1)[2] <- "Drought"
colnames(area_phenotype_tp1)[3] <- "Heat"
colnames(area_phenotype_tp1)[4] <- "Heat.Drought"


#### After making phenotype file once, can start here:
# save file for GAPIT input
write.table(area_phenotype_tp1, "./area_phenotype_tp1_GAPIT.txt")
#read in phenotype file
area_phenotype_tp1 <- read.table("./area_phenotype_tp1_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(area_phenotype_tp1)

#check phenotype file
hist(area_phenotype_tp1$Drought)
mean(area_phenotype_tp1$Drought)
range(area_phenotype_tp1$Drought)
sd(area_phenotype_tp1$Drought)
which(is.na(area_phenotype_tp1$Drought))

hist(area_phenotype_tp1$Heat)
mean(area_phenotype_tp1$Heat)
range(area_phenotype_tp1$Heat)
sd(area_phenotype_tp1$Heat)
which(is.na(area_phenotype_tp1$Heat))

hist(area_phenotype_tp1$Heat.Drought)
mean(area_phenotype_tp1$Heat.Drought)
range(area_phenotype_tp1$Heat.Drought)
sd(area_phenotype_tp1$Heat.Drought)
which(is.na(area_phenotype_tp1$Heat.Drought))

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(area_phenotype_tp1$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp1")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp1")
#run GAPIT
area_gapit <- GAPIT(
  Y = area_phenotype_tp1,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU"),
  Multiple_analysis = TRUE
)


################################################################
### GAPIT GWAS for Percent Damage (difference)
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp1")

# Percent Damage Early Day: control and drought
####
avgdamage <- group_by(DFs_144, Accession, Treatment, group) %>% summarise(avg_damage = mean(percent_damage, na.rm = TRUE))

avgdamage_tp1 <- aggregate(data = avgdamage[avgdamage$group == "22",], avg_damage~Accession+Treatment, FUN= function(i)mean(i))
head(avgdamage_tp1)

# check to make sure each phenotype and treatment has data
table(avgdamage_tp1$Accession, avgdamage_tp1$Treatment)
avgdamage_tp1[avgdamage_tp1 == "WW"] <- "control"
avgdamage_tp1[avgdamage_tp1 == "WL"] <- "drought"

#make df with only WL data from drought experiment (control temp)
droughtWL_avgdamage_tp1 <- avgdamage_tp1[avgdamage_tp1$Treatment == "drought",]

damage_dif4 <- setNames(data.frame(do.call("rbind",lapply(split(avgdamage_tp1, avgdamage_tp1$Accession), function(i)diff(i$avg_damage)))), c("damage_dif"))
damage_dif4$Accession <- row.names(damage_dif4)
head(damage_dif4)
row.names(damage_dif4) <- NULL


####
#control and heat
####
heatDFs_144$percent_damage <- as.numeric(heatDFs_144$percent_damage)
avgdamageheat <- group_by(heatDFs_144, Accession, Treatment, group) %>% summarise(avg_damage = mean(percent_damage, na.rm = TRUE))

heatavgdamage_tp1 <- aggregate(data = avgdamageheat[avgdamageheat$group == "22",], avg_damage~Accession+Treatment, FUN= function(i)mean(i))
head(heatavgdamage_tp1)

table(heatavgdamage_tp1$Accession, heatavgdamage_tp1$Treatment)
heatavgdamage_tp1[heatavgdamage_tp1 == "WW"] <- "heat"
heatavgdamage_tp1[heatavgdamage_tp1 == "WL"] <- "heat.drought"


# make df with only WW data from heat and control temps
heatWW_avgdamage_tp1 <- heatavgdamage_tp1[heatavgdamage_tp1$Treatment == "heat",]

droughtWW_avgdamage_tp1 <- avgdamage_tp1[avgdamage_tp1$Treatment == "control",]

WW_avgdamage_tp1 <- rbind(heatWW_avgdamage_tp1, droughtWW_avgdamage_tp1)

#check to make sure all accessions have data for both treatments on last day
table(WW_avgdamage_tp1$Accession, WW_avgdamage_tp1$Treatment)

damage_dif5 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgdamage_tp1, WW_avgdamage_tp1$Accession), function(i)diff(i$avg_damage)))), c("damage_dif2"))
damage_dif5$Accession <- row.names(damage_dif5)
head(damage_dif5)
row.names(damage_dif5) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgdamage_tp1 <- heatavgdamage_tp1[heatavgdamage_tp1$Treatment == "heat.drought",]

control_and_heatWL_avgdamage_tp1 <- rbind(heatWL_avgdamage_tp1, droughtWW_avgdamage_tp1)

#check to make sure all accessions have data for both treatments on last day
table(control_and_heatWL_avgdamage_tp1$Accession, control_and_heatWL_avgdamage_tp1$Treatment)

damage_dif6 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgdamage_tp1, control_and_heatWL_avgdamage_tp1$Accession), function(i)diff(i$avg_damage)))), c("damage_dif3"))
damage_dif6$Accession <- row.names(damage_dif6)
head(damage_dif6)
row.names(damage_dif6) <- NULL

# create phenotype (Y) file for damage
temporary <- join(damage_dif4, damage_dif5, by="Accession")
damage_phenotype_tp1 <- join(temporary, damage_dif6, by="Accession")
damage_phenotype_tp1 <- drop_na(damage_phenotype_tp1)

#organize it
damage_phenotype_tp1 <- dplyr::select(damage_phenotype_tp1, Accession, damage_dif, damage_dif2, damage_dif3)
# rename columns
colnames(damage_phenotype_tp1)[2] <- "Drought"
colnames(damage_phenotype_tp1)[3] <- "Heat"
colnames(damage_phenotype_tp1)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(damage_phenotype_tp1, "./damage_phenotype_tp1_GAPIT.txt")


#### After making phenotype file once, can start here:
#read in phenotype file
damage_phenotype_tp1 <- read.table("./damage_phenotype_tp1_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(damage_phenotype_tp1)

#check phenotype file
hist(damage_phenotype_tp1$Drought)
mean(damage_phenotype_tp1$Drought)
range(damage_phenotype_tp1$Drought)
sd(damage_phenotype_tp1$Drought)
which(is.na(damage_phenotype_tp1$Drought))

hist(damage_phenotype_tp1$Heat)
mean(damage_phenotype_tp1$Heat)
range(damage_phenotype_tp1$Heat)
sd(damage_phenotype_tp1$Heat)
which(is.na(damage_phenotype_tp1$Heat))

hist(damage_phenotype_tp1$Heat.Drought)
mean(damage_phenotype_tp1$Heat.Drought)
range(damage_phenotype_tp1$Heat.Drought)
sd(damage_phenotype_tp1$Heat.Drought)
which(is.na(damage_phenotype_tp1$Heat.Drought))

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(damage_phenotype_tp1$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap_og as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp1")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp1")
#run GAPIT
damage_gapit <- GAPIT(
  Y = damage_phenotype_tp1,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU"),
  Multiple_analysis = TRUE
)

#####


#*******************************************************
# ***TP2*** Comparison GWAS with GAPIT (organized by trait)
################################################################
### GAPIT GWAS for Height (difference)
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp2")

# Height TP2: control and drought
####
avgheight <- group_by(DFs_144, Accession, Treatment, group) %>% summarise(avg_height = mean(cor_height_above_reference, na.rm = TRUE))
 
avgheight_tp2 <- aggregate(data = avgheight[avgheight$group == "28/29" ,], avg_height~Accession+Treatment, FUN= function(i)mean(i))
head(avgheight_tp2)

# check to make sure each phenotype and treatment has data
table(avgheight_tp2$Accession, avgheight_tp2$Treatment)
avgheight_tp2[avgheight_tp2 == "WW"] <- "control"
avgheight_tp2[avgheight_tp2 == "WL"] <- "drought"

#make df with only WW data from drought exp
droughtWW_avgheight_tp2 <- avgheight_tp2[avgheight_tp2$Treatment == "control",]

#make df with only WL data from drought experiment (control temp)
droughtWL_avgheight_tp2 <- avgheight_tp2[avgheight_tp2$Treatment == "drought",]

height_dif47 <- setNames(data.frame(do.call("rbind",lapply(split(avgheight_tp2, avgheight_tp2$Accession), function(i)diff(i$avg_height)))), c("height_dif"))
height_dif47$Accession <- row.names(height_dif47)
head(height_dif47)
row.names(height_dif47) <- NULL


####
#control and heat
####
heatDFs_144$cor_height_above_reference <- as.numeric(heatDFs_144$cor_height_above_reference)
avgheightheat <- group_by(heatDFs_144, Accession, Treatment, group) %>% summarise(avg_height = mean(cor_height_above_reference, na.rm = TRUE))
 
heatavgheight_tp2 <- aggregate(data = avgheightheat[avgheightheat$group == "28/29" ,], avg_height~Accession+Treatment, FUN= function(i)mean(i))
head(heatavgheight_tp2)

table(heatavgheight_tp2$Accession, heatavgheight_tp2$Treatment)
heatavgheight_tp2[heatavgheight_tp2 == "WW"] <- "heat"
heatavgheight_tp2[heatavgheight_tp2 == "WL"] <- "heat.drought"

# make df with only WW data from heat and control temps
heatWW_avgheight_tp2 <- heatavgheight_tp2[heatavgheight_tp2$Treatment == "heat",]

WW_avgheight_tp2 <- rbind(heatWW_avgheight_tp2, droughtWW_avgheight_tp2)

#check to make sure all accessions have data for both treatments on last day
table(WW_avgheight_tp2$Accession, WW_avgheight_tp2$Treatment)

height_dif57 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgheight_tp2, WW_avgheight_tp2$Accession), function(i)diff(i$avg_height)))), c("height_dif2"))
height_dif57$Accession <- row.names(height_dif57)
head(height_dif57)
row.names(height_dif57) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgheight_tp2 <- heatavgheight_tp2[heatavgheight_tp2$Treatment == "heat.drought",]

control_and_heatWL_avgheight_tp2 <- rbind(heatWL_avgheight_tp2, droughtWW_avgheight_tp2)

#check to make sure all accessions have data for both treatments on last day
table(control_and_heatWL_avgheight_tp2$Accession, control_and_heatWL_avgheight_tp2$Treatment)

height_dif67 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgheight_tp2, control_and_heatWL_avgheight_tp2$Accession), function(i)diff(i$avg_height)))), c("height_dif3"))
height_dif67$Accession <- row.names(height_dif67)
head(height_dif67)
row.names(height_dif67) <- NULL

# create phenotype (Y) file for height
temporary <- join(height_dif47, height_dif57, by="Accession")
height_phenotype_tp2 <- join(temporary, height_dif67, by="Accession")
height_phenotype_tp2 <- drop_na(height_phenotype_tp2)

#organize it
height_phenotype_tp2 <- dplyr::select(height_phenotype_tp2, Accession, height_dif, height_dif2, height_dif3)
# rename columns
colnames(height_phenotype_tp2)[2] <- "Drought"
colnames(height_phenotype_tp2)[3] <- "Heat"
colnames(height_phenotype_tp2)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(height_phenotype_tp2, "./height_phenotype_tp2_GAPIT.txt")


#### After making phenotype file once, can start here:
#read in phenotype file
height_phenotype_tp2 <- read.table("./height_phenotype_tp2_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(height_phenotype_tp2)

#check phenotype file
hist(height_phenotype_tp2$Drought)
mean(height_phenotype_tp2$Drought)
range(height_phenotype_tp2$Drought)
sd(height_phenotype_tp2$Drought)
which(is.na(height_phenotype_tp2$Drought))

hist(height_phenotype_tp2$Heat)
mean(height_phenotype_tp2$Heat)
range(height_phenotype_tp2$Heat)
sd(height_phenotype_tp2$Heat)
which(is.na(height_phenotype_tp2$Heat))

hist(height_phenotype_tp2$Heat.Drought)
mean(height_phenotype_tp2$Heat.Drought)
range(height_phenotype_tp2$Heat.Drought)
sd(height_phenotype_tp2$Heat.Drought)
which(is.na(height_phenotype_tp2$Heat.Drought))

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(height_phenotype_tp2$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap_og as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp2")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp2")
#run GAPIT
height_gapit <- GAPIT(
  Y = height_phenotype_tp2,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU"),
  Multiple_analysis = TRUE
)

################################################################
### GAPIT GWAS for Area (difference)
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp2")

# Area Early Day: control and drought
####
avgarea <- group_by(DFs_144, Accession, Treatment, group) %>% summarise(avg_area = mean(cor_area, na.rm = TRUE))

avgarea_tp2 <- aggregate(data = avgarea[avgarea$group == "28/29" ,], avg_area~Accession+Treatment, FUN= function(i)mean(i))
head(avgarea_tp2)

# check to make sure each phenotype and treatment has data
table(avgarea_tp2$Accession, avgarea_tp2$Treatment)
avgarea_tp2[avgarea_tp2 == "WW"] <- "control"
avgarea_tp2[avgarea_tp2 == "WL"] <- "drought"

#make df with only WL data from drought experiment (control temp)
droughtWL_avgarea_tp2 <- avgarea_tp2[avgarea_tp2$Treatment == "drought",]

area_dif47 <- setNames(data.frame(do.call("rbind",lapply(split(avgarea_tp2, avgarea_tp2$Accession), function(i)diff(i$avg_area)))), c("area_dif"))
area_dif47$Accession <- row.names(area_dif47)
head(area_dif47)
row.names(area_dif47) <- NULL


####
#control and heat
####
heatDFs_144$cor_area <- as.numeric(heatDFs_144$cor_area)
avgareaheat <- group_by(heatDFs_144, Accession, Treatment, group) %>% summarise(avg_area = mean(cor_area, na.rm = TRUE))
 

heatavgarea_tp2 <- aggregate(data = avgareaheat[avgareaheat$group == "28/29" ,], avg_area~Accession+Treatment, FUN= function(i)mean(i))
head(heatavgarea_tp2)

table(heatavgarea_tp2$Accession, heatavgarea_tp2$Treatment)
heatavgarea_tp2[heatavgarea_tp2 == "WW"] <- "heat"
heatavgarea_tp2[heatavgarea_tp2 == "WL"] <- "heat.drought"


# make df with only WW data from heat and control temps
heatWW_avgarea_tp2 <- heatavgarea_tp2[heatavgarea_tp2$Treatment == "heat",]

droughtWW_avgarea_tp2 <- avgarea_tp2[avgarea_tp2$Treatment == "control",]

WW_avgarea_tp2 <- rbind(heatWW_avgarea_tp2, droughtWW_avgarea_tp2)

#check to make sure all accessions have data for both treatments on last day
table(WW_avgarea_tp2$Accession, WW_avgarea_tp2$Treatment)

area_dif57 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgarea_tp2, WW_avgarea_tp2$Accession), function(i)diff(i$avg_area)))), c("area_dif2"))
area_dif57$Accession <- row.names(area_dif57)
head(area_dif57)
row.names(area_dif57) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgarea_tp2 <- heatavgarea_tp2[heatavgarea_tp2$Treatment == "heat.drought",]

control_and_heatWL_avgarea_tp2 <- rbind(heatWL_avgarea_tp2, droughtWW_avgarea_tp2)

#check to make sure all accessions have data for both treatments on last day
table(control_and_heatWL_avgarea_tp2$Accession, control_and_heatWL_avgarea_tp2$Treatment)

area_dif67 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgarea_tp2, control_and_heatWL_avgarea_tp2$Accession), function(i)diff(i$avg_area)))), c("area_dif3"))
area_dif67$Accession <- row.names(area_dif67)
head(area_dif67)
row.names(area_dif67) <- NULL

# create phenotype (Y) file for area
temporary <- join(area_dif47, area_dif57, by="Accession")
area_phenotype_tp2 <- join(temporary, area_dif67, by="Accession")
area_phenotype_tp2 <- drop_na(area_phenotype_tp2)

#organize it
area_phenotype_tp2 <- dplyr::select(area_phenotype_tp2, Accession, area_dif, area_dif2, area_dif3)
# rename columns
colnames(area_phenotype_tp2)[2] <- "Drought"
colnames(area_phenotype_tp2)[3] <- "Heat"
colnames(area_phenotype_tp2)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(area_phenotype_tp2, "./area_phenotype_tp2_GAPIT.txt")


#### After making phenotype file once, can start here:
#read in phenotype file
area_phenotype_tp2 <- read.table("./area_phenotype_tp2_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(area_phenotype_tp2)

#check phenotype file
hist(area_phenotype_tp2$Drought)
mean(area_phenotype_tp2$Drought)
range(area_phenotype_tp2$Drought)
sd(area_phenotype_tp2$Drought)
which(is.na(area_phenotype_tp2$Drought))

hist(area_phenotype_tp2$Heat)
mean(area_phenotype_tp2$Heat)
range(area_phenotype_tp2$Heat)
sd(area_phenotype_tp2$Heat)
which(is.na(area_phenotype_tp2$Heat))

hist(area_phenotype_tp2$Heat.Drought)
mean(area_phenotype_tp2$Heat.Drought)
range(area_phenotype_tp2$Heat.Drought)
sd(area_phenotype_tp2$Heat.Drought)
which(is.na(area_phenotype_tp2$Heat.Drought))

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(area_phenotype_tp2$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp2")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp2")
#run GAPIT
area_gapit <- GAPIT(
  Y = area_phenotype_tp2,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU"),
  Multiple_analysis = TRUE
)


################################################################
### GAPIT GWAS for Percent Damage (difference)
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp2")

# Percent Damage Early Day: control and drought
####
avgdamage <- group_by(DFs_144, Accession, Treatment, group) %>% summarise(avg_damage = mean(percent_damage, na.rm = TRUE))
 
avgdamage_tp2 <- aggregate(data = avgdamage[avgdamage$group == "28/29" ,], avg_damage~Accession+Treatment, FUN= function(i)mean(i))
head(avgdamage_tp2)

# check to make sure each phenotype and treatment has data
table(avgdamage_tp2$Accession, avgdamage_tp2$Treatment)
avgdamage_tp2[avgdamage_tp2 == "WW"] <- "control"
avgdamage_tp2[avgdamage_tp2 == "WL"] <- "drought"

#make df with only WL data from drought experiment (control temp)
droughtWL_avgdamage_tp2 <- avgdamage_tp2[avgdamage_tp2$Treatment == "drought",]

damage_dif47 <- setNames(data.frame(do.call("rbind",lapply(split(avgdamage_tp2, avgdamage_tp2$Accession), function(i)diff(i$avg_damage)))), c("damage_dif"))
damage_dif47$Accession <- row.names(damage_dif47)
head(damage_dif47)
row.names(damage_dif47) <- NULL


####
#control and heat
####
heatDFs_144$percent_damage <- as.numeric(heatDFs_144$percent_damage)
avgdamageheat <- group_by(heatDFs_144, Accession, Treatment, group) %>% summarise(avg_damage = mean(percent_damage, na.rm = TRUE))

heatavgdamage_tp2 <- aggregate(data = avgdamageheat[avgdamageheat$group == "28/29" ,], avg_damage~Accession+Treatment, FUN= function(i)mean(i))
head(heatavgdamage_tp2)

table(heatavgdamage_tp2$Accession, heatavgdamage_tp2$Treatment)
heatavgdamage_tp2[heatavgdamage_tp2 == "WW"] <- "heat"
heatavgdamage_tp2[heatavgdamage_tp2 == "WL"] <- "heat.drought"


# make df with only WW data from heat and control temps
heatWW_avgdamage_tp2 <- heatavgdamage_tp2[heatavgdamage_tp2$Treatment == "heat",]

droughtWW_avgdamage_tp2 <- avgdamage_tp2[avgdamage_tp2$Treatment == "control",]

WW_avgdamage_tp2 <- rbind(heatWW_avgdamage_tp2, droughtWW_avgdamage_tp2)

#check to make sure all accessions have data for both treatments on last day
table(WW_avgdamage_tp2$Accession, WW_avgdamage_tp2$Treatment)

damage_dif57 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgdamage_tp2, WW_avgdamage_tp2$Accession), function(i)diff(i$avg_damage)))), c("damage_dif2"))
damage_dif57$Accession <- row.names(damage_dif57)
head(damage_dif57)
row.names(damage_dif57) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgdamage_tp2 <- heatavgdamage_tp2[heatavgdamage_tp2$Treatment == "heat.drought",]

control_and_heatWL_avgdamage_tp2 <- rbind(heatWL_avgdamage_tp2, droughtWW_avgdamage_tp2)

#check to make sure all accessions have data for both treatments on last day
table(control_and_heatWL_avgdamage_tp2$Accession, control_and_heatWL_avgdamage_tp2$Treatment)

damage_dif67 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgdamage_tp2, control_and_heatWL_avgdamage_tp2$Accession), function(i)diff(i$avg_damage)))), c("damage_dif3"))
damage_dif67$Accession <- row.names(damage_dif67)
head(damage_dif67)
row.names(damage_dif67) <- NULL

# create phenotype (Y) file for damage
temporary <- join(damage_dif47, damage_dif57, by="Accession")
damage_phenotype_tp2 <- join(temporary, damage_dif67, by="Accession")
damage_phenotype_tp2 <- drop_na(damage_phenotype_tp2)

#organize it
damage_phenotype_tp2 <- dplyr::select(damage_phenotype_tp2, Accession, damage_dif, damage_dif2, damage_dif3)
# rename columns
colnames(damage_phenotype_tp2)[2] <- "Drought"
colnames(damage_phenotype_tp2)[3] <- "Heat"
colnames(damage_phenotype_tp2)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(damage_phenotype_tp2, "./damage_phenotype_tp2_GAPIT.txt")


#### After making phenotype file once, can start here:
#read in phenotype file
damage_phenotype_tp2 <- read.table("./damage_phenotype_tp2_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(damage_phenotype_tp2)

#check phenotype file
hist(damage_phenotype_tp2$Drought)
mean(damage_phenotype_tp2$Drought)
range(damage_phenotype_tp2$Drought)
sd(damage_phenotype_tp2$Drought)
which(is.na(damage_phenotype_tp2$Drought))

hist(damage_phenotype_tp2$Heat)
mean(damage_phenotype_tp2$Heat)
range(damage_phenotype_tp2$Heat)
sd(damage_phenotype_tp2$Heat)
which(is.na(damage_phenotype_tp2$Heat))

hist(damage_phenotype_tp2$Heat.Drought)
mean(damage_phenotype_tp2$Heat.Drought)
range(damage_phenotype_tp2$Heat.Drought)
sd(damage_phenotype_tp2$Heat.Drought)
which(is.na(damage_phenotype_tp2$Heat.Drought))

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(damage_phenotype_tp2$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap_og as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp2")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp2")
#run GAPIT
damage_gapit <- GAPIT(
  Y = damage_phenotype_tp2,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU"),
  Multiple_analysis = TRUE
)

#####

#*************************************************************


# ***TP3*** Comparison GWAS with GAPIT (organized by trait)
################################################################
### GAPIT GWAS for Height (difference)
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp3")

# Height Middle Day: control and drought
####
avgheight <- group_by(DFs_144, Accession, Treatment, group) %>% summarise(avg_height = mean(cor_height_above_reference, na.rm = TRUE))
 
avgheight_tp3 <- aggregate(data = avgheight[avgheight$group == "36/37",], avg_height~Accession+Treatment, FUN= function(i)mean(i))
head(avgheight_tp3)

# check to make sure each phenotype and treatment has data
table(avgheight_tp3$Accession, avgheight_tp3$Treatment)
avgheight_tp3[avgheight_tp3 == "WW"] <- "control"
avgheight_tp3[avgheight_tp3 == "WL"] <- "drought"

#make df with only WL data from drought experiment (control temp)
droughtWL_avgheight_tp3 <- avgheight_tp3[avgheight_tp3$Treatment == "drought",]

height_dif7 <- setNames(data.frame(do.call("rbind",lapply(split(avgheight_tp3, avgheight_tp3$Accession), function(i)diff(i$avg_height)))), c("height_dif"))
height_dif7$Accession <- row.names(height_dif7)
head(height_dif7)
row.names(height_dif7) <- NULL


####
#control and heat
####
heatDFs_144$cor_height_above_reference <- as.numeric(heatDFs_144$cor_height_above_reference)
avgheightheat <- group_by(heatDFs_144, Accession, Treatment, group) %>% summarise(avg_height = mean(cor_height_above_reference, na.rm = TRUE))
 

heatavgheight_tp3 <- aggregate(data = avgheightheat[avgheightheat$group == "36/37",], avg_height~Accession+Treatment, FUN= function(i)mean(i))
head(heatavgheight_tp3)

table(heatavgheight_tp3$Accession, heatavgheight_tp3$Treatment)
heatavgheight_tp3[heatavgheight_tp3 == "WW"] <- "heat"
heatavgheight_tp3[heatavgheight_tp3 == "WL"] <- "heat.drought"

# make df with only WW data from heat and control temps
heatWW_avgheight_tp3 <- heatavgheight_tp3[heatavgheight_tp3$Treatment == "heat",]

droughtWW_avgheight_tp3 <- avgheight_tp3[avgheight_tp3$Treatment == "control",]

WW_avgheight_tp3 <- rbind(heatWW_avgheight_tp3, droughtWW_avgheight_tp3)

#check to make sure all accessions have data for both treatments on last day
table(WW_avgheight_tp3$Accession, WW_avgheight_tp3$Treatment)

height_dif8 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgheight_tp3, WW_avgheight_tp3$Accession), function(i)diff(i$avg_height)))), c("height_dif2"))
height_dif8$Accession <- row.names(height_dif8)
head(height_dif8)
row.names(height_dif8) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgheight_tp3 <- heatavgheight_tp3[heatavgheight_tp3$Treatment == "heat.drought",]

control_and_heatWL_avgheight_tp3 <- rbind(heatWL_avgheight_tp3, droughtWW_avgheight_tp3)

#check to make sure all accessions have data for both treatments on last day
table(control_and_heatWL_avgheight_tp3$Accession, control_and_heatWL_avgheight_tp3$Treatment)

height_dif9 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgheight_tp3, control_and_heatWL_avgheight_tp3$Accession), function(i)diff(i$avg_height)))), c("height_dif3"))
height_dif9$Accession <- row.names(height_dif9)
head(height_dif9)
row.names(height_dif9) <- NULL

# create phenotype (Y) file for height
temporary <- join(height_dif7, height_dif8, by="Accession")
height_phenotype_tp3 <- join(temporary, height_dif9, by="Accession")
height_phenotype_tp3 <- drop_na(height_phenotype_tp3)

#organize it
height_phenotype_tp3 <- dplyr::select(height_phenotype_tp3, Accession, height_dif, height_dif2, height_dif3)
# rename columns
colnames(height_phenotype_tp3)[2] <- "Drought"
colnames(height_phenotype_tp3)[3] <- "Heat"
colnames(height_phenotype_tp3)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(height_phenotype_tp3, "./height_phenotype_tp3_GAPIT.txt")


#### After making phenotype file once, can start here:
#read in phenotype file
height_phenotype_tp3 <- read.table("./height_phenotype_tp3_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(height_phenotype_tp3)

#check phenotype file
hist(height_phenotype_tp3$Drought)
mean(height_phenotype_tp3$Drought)
range(height_phenotype_tp3$Drought)
sd(height_phenotype_tp3$Drought)
which(is.na(height_phenotype_tp3$Drought))

hist(height_phenotype_tp3$Heat)
mean(height_phenotype_tp3$Heat)
range(height_phenotype_tp3$Heat)
sd(height_phenotype_tp3$Heat)
which(is.na(height_phenotype_tp3$Heat))

hist(height_phenotype_tp3$Heat.Drought)
mean(height_phenotype_tp3$Heat.Drought)
range(height_phenotype_tp3$Heat.Drought)
sd(height_phenotype_tp3$Heat.Drought)
which(is.na(height_phenotype_tp3$Heat.Drought))

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(height_phenotype_tp3$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap_og as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp3")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp3")
#run GAPIT
height_gapit <- GAPIT(
  Y = height_phenotype_tp3,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU"),
  Multiple_analysis = TRUE
)

################################################################
### GAPIT GWAS for Area (difference)
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp3")

# Area Middle Day: control and drought
####
avgarea <- group_by(DFs_144, Accession, Treatment, group) %>% summarise(avg_area = mean(cor_area, na.rm = TRUE))
 
avgarea_tp3 <- aggregate(data = avgarea[avgarea$group == "36/37",], avg_area~Accession+Treatment, FUN= function(i)mean(i))
head(avgarea_tp3)

# check to make sure each phenotype and treatment has data
table(avgarea_tp3$Accession, avgarea_tp3$Treatment)
avgarea_tp3[avgarea_tp3 == "WW"] <- "control"
avgarea_tp3[avgarea_tp3 == "WL"] <- "drought"

#make df with only WL data from drought experiment (control temp)
droughtWL_avgarea_tp3 <- avgarea_tp3[avgarea_tp3$Treatment == "drought",]

area_dif7 <- setNames(data.frame(do.call("rbind",lapply(split(avgarea_tp3, avgarea_tp3$Accession), function(i)diff(i$avg_area)))), c("area_dif"))
area_dif7$Accession <- row.names(area_dif7)
head(area_dif7)
row.names(area_dif7) <- NULL


####
#control and heat
####
heatDFs_144$cor_area <- as.numeric(heatDFs_144$cor_area)
avgareaheat <- group_by(heatDFs_144, Accession, Treatment, group) %>% summarise(avg_area = mean(cor_area, na.rm = TRUE))
 

heatavgarea_tp3 <- aggregate(data = avgareaheat[avgareaheat$group == "36/37",], avg_area~Accession+Treatment, FUN= function(i)mean(i))
head(heatavgarea_tp3)

table(heatavgarea_tp3$Accession, heatavgarea_tp3$Treatment)
heatavgarea_tp3[heatavgarea_tp3 == "WW"] <- "heat"
heatavgarea_tp3[heatavgarea_tp3 == "WL"] <- "heat.drought"


# make df with only WW data from heat and control temps
heatWW_avgarea_tp3 <- heatavgarea_tp3[heatavgarea_tp3$Treatment == "heat",]

droughtWW_avgarea_tp3 <- avgarea_tp3[avgarea_tp3$Treatment == "control",]

WW_avgarea_tp3 <- rbind(heatWW_avgarea_tp3, droughtWW_avgarea_tp3)

#check to make sure all accessions have data for both treatments on last day
table(WW_avgarea_tp3$Accession, WW_avgarea_tp3$Treatment)

area_dif8 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgarea_tp3, WW_avgarea_tp3$Accession), function(i)diff(i$avg_area)))), c("area_dif2"))
area_dif8$Accession <- row.names(area_dif8)
head(area_dif8)
row.names(area_dif8) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgarea_tp3 <- heatavgarea_tp3[heatavgarea_tp3$Treatment == "heat.drought",]

control_and_heatWL_avgarea_tp3 <- rbind(heatWL_avgarea_tp3, droughtWW_avgarea_tp3)

#check to make sure all accessions have data for both treatments on last day
table(control_and_heatWL_avgarea_tp3$Accession, control_and_heatWL_avgarea_tp3$Treatment)

area_dif9 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgarea_tp3, control_and_heatWL_avgarea_tp3$Accession), function(i)diff(i$avg_area)))), c("area_dif3"))
area_dif9$Accession <- row.names(area_dif9)
head(area_dif9)
row.names(area_dif9) <- NULL

# create phenotype (Y) file for area
temporary <- join(area_dif7, area_dif8, by="Accession")
area_phenotype_tp3 <- join(temporary, area_dif9, by="Accession")
area_phenotype_tp3 <- drop_na(area_phenotype_tp3)

#organize it
area_phenotype_tp3 <- dplyr::select(area_phenotype_tp3, Accession, area_dif, area_dif2, area_dif3)
# rename columns
colnames(area_phenotype_tp3)[2] <- "Drought"
colnames(area_phenotype_tp3)[3] <- "Heat"
colnames(area_phenotype_tp3)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(area_phenotype_tp3, "./area_phenotype_tp3_GAPIT.txt")



#### After making phenotype file once, can start here:
#read in phenotype file
area_phenotype_tp3 <- read.table("./area_phenotype_tp3_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(area_phenotype_tp3)

#check phenotype file
hist(area_phenotype_tp3$Drought)
mean(area_phenotype_tp3$Drought)
range(area_phenotype_tp3$Drought)
sd(area_phenotype_tp3$Drought)
which(is.na(area_phenotype_tp3$Drought))

hist(area_phenotype_tp3$Heat)
mean(area_phenotype_tp3$Heat)
range(area_phenotype_tp3$Heat)
sd(area_phenotype_tp3$Heat)
which(is.na(area_phenotype_tp3$Heat))

hist(area_phenotype_tp3$Heat.Drought)
mean(area_phenotype_tp3$Heat.Drought)
range(area_phenotype_tp3$Heat.Drought)
sd(area_phenotype_tp3$Heat.Drought)
which(is.na(area_phenotype_tp3$Heat.Drought))

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(area_phenotype_tp3$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap_og as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp3")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp3")
#run GAPIT
area_gapit <- GAPIT(
  Y = area_phenotype_tp3,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU"),
  Multiple_analysis = TRUE
)


################################################################
### GAPIT GWAS for Percent Damage (difference)
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp3")

# Percent Damage Middle Day: control and drought
####
avgdamage <- group_by(DFs_144, Accession, Treatment, group) %>% summarise(avg_damage = mean(percent_damage, na.rm = TRUE))

avgdamage_tp3 <- aggregate(data = avgdamage[avgdamage$group == "36/37",], avg_damage~Accession+Treatment, FUN= function(i)mean(i))
head(avgdamage_tp3)

# check to make sure each phenotype and treatment has data
table(avgdamage_tp3$Accession, avgdamage_tp3$Treatment)
avgdamage_tp3[avgdamage_tp3 == "WW"] <- "control"
avgdamage_tp3[avgdamage_tp3 == "WL"] <- "drought"

#make df with only WL data from drought experiment (control temp)
droughtWL_avgdamage_tp3 <- avgdamage_tp3[avgdamage_tp3$Treatment == "drought",]

damage_dif7 <- setNames(data.frame(do.call("rbind",lapply(split(avgdamage_tp3, avgdamage_tp3$Accession), function(i)diff(i$avg_damage)))), c("damage_dif"))
damage_dif7$Accession <- row.names(damage_dif7)
head(damage_dif7)
row.names(damage_dif7) <- NULL


####
#control and heat
####
heatDFs_144$percent_damage <- as.numeric(heatDFs_144$percent_damage)
avgdamageheat <- group_by(heatDFs_144, Accession, Treatment, group) %>% summarise(avg_damage = mean(percent_damage, na.rm = TRUE))

heatavgdamage_tp3 <- aggregate(data = avgdamageheat[avgdamageheat$group == "36/37",], avg_damage~Accession+Treatment, FUN= function(i)mean(i))
head(heatavgdamage_tp3)

table(heatavgdamage_tp3$Accession, heatavgdamage_tp3$Treatment)
heatavgdamage_tp3[heatavgdamage_tp3 == "WW"] <- "heat"
heatavgdamage_tp3[heatavgdamage_tp3 == "WL"] <- "heat.drought"

# make df with only WW data from heat and control temps
heatWW_avgdamage_tp3 <- heatavgdamage_tp3[heatavgdamage_tp3$Treatment == "heat",]

droughtWW_avgdamage_tp3 <- avgdamage_tp3[avgdamage_tp3$Treatment == "control",]

WW_avgdamage_tp3 <- rbind(heatWW_avgdamage_tp3, droughtWW_avgdamage_tp3)

#check to make sure all accessions have data for both treatments on last day
table(WW_avgdamage_tp3$Accession, WW_avgdamage_tp3$Treatment)

damage_dif8 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgdamage_tp3, WW_avgdamage_tp3$Accession), function(i)diff(i$avg_damage)))), c("damage_dif2"))
damage_dif8$Accession <- row.names(damage_dif8)
head(damage_dif8)
row.names(damage_dif8) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgdamage_tp3 <- heatavgdamage_tp3[heatavgdamage_tp3$Treatment == "heat.drought",]

control_and_heatWL_avgdamage_tp3 <- rbind(heatWL_avgdamage_tp3, droughtWW_avgdamage_tp3)

#check to make sure all accessions have data for both treatments on last day
table(control_and_heatWL_avgdamage_tp3$Accession, control_and_heatWL_avgdamage_tp3$Treatment)

damage_dif9 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgdamage_tp3, control_and_heatWL_avgdamage_tp3$Accession), function(i)diff(i$avg_damage)))), c("damage_dif3"))
damage_dif9$Accession <- row.names(damage_dif9)
head(damage_dif9)
row.names(damage_dif9) <- NULL

# create phenotype (Y) file for damage
temporary <- join(damage_dif7, damage_dif8, by="Accession")
damage_phenotype_tp3 <- join(temporary, damage_dif9, by="Accession")
damage_phenotype_tp3 <- drop_na(damage_phenotype_tp3)

#organize it
damage_phenotype_tp3 <- dplyr::select(damage_phenotype_tp3, Accession, damage_dif, damage_dif2, damage_dif3)
# rename columns
colnames(damage_phenotype_tp3)[2] <- "Drought"
colnames(damage_phenotype_tp3)[3] <- "Heat"
colnames(damage_phenotype_tp3)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(damage_phenotype_tp3, "./damage_phenotype_tp3_GAPIT.txt")



#### After making phenotype file once, can start here:
#read in phenotype file
damage_phenotype_tp3 <- read.table("./damage_phenotype_tp3_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(damage_phenotype_tp3)

#check phenotype file
hist(damage_phenotype_tp3$Drought)
mean(damage_phenotype_tp3$Drought)
range(damage_phenotype_tp3$Drought)
sd(damage_phenotype_tp3$Drought)
which(is.na(damage_phenotype_tp3$Drought))

hist(damage_phenotype_tp3$Heat)
mean(damage_phenotype_tp3$Heat)
range(damage_phenotype_tp3$Heat)
sd(damage_phenotype_tp3$Heat)
which(is.na(damage_phenotype_tp3$Heat))

hist(damage_phenotype_tp3$Heat.Drought)
mean(damage_phenotype_tp3$Heat.Drought)
range(damage_phenotype_tp3$Heat.Drought)
sd(damage_phenotype_tp3$Heat.Drought)
which(is.na(damage_phenotype_tp3$Heat.Drought))

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(damage_phenotype_tp3$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap_og as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp3")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp3")
#run GAPIT
damage_gapit <- GAPIT(
  Y = damage_phenotype_tp3,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU"),
  Multiple_analysis = TRUE
)

#####


#*************************************************************
# ***TP4*** Comparison GWAS with GAPIT (organized by trait)
################################################################
### GAPIT GWAS for Height (difference)
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp4")

# Height Final Day: control and drought
####
avgheight <- group_by(DFs_144, Accession, Treatment, group) %>% summarise(avg_height = mean(cor_height_above_reference, na.rm = TRUE))

avgheight_tp4 <- aggregate(data = avgheight[avgheight$group == "44/45",], avg_height~Accession+Treatment, FUN= function(i)mean(i))
head(avgheight_tp4)

# check to make sure each phenotype and treatment has data
table(avgheight_tp4$Accession, avgheight_tp4$Treatment)
avgheight_tp4[avgheight_tp4 == "WW"] <- "control"
avgheight_tp4[avgheight_tp4 == "WL"] <- "drought"

#make df with only WL data from drought experiment (control temp)
droughtWL_avgheight_tp4 <- avgheight_tp4[avgheight_tp4$Treatment == "drought",]

height_dif <- setNames(data.frame(do.call("rbind",lapply(split(avgheight_tp4, avgheight_tp4$Accession), function(i)diff(i$avg_height)))), c("height_dif"))
height_dif$Accession <- row.names(height_dif)
head(height_dif)
row.names(height_dif) <- NULL


####
#control and heat
####
avgheightheat <- group_by(heatDFs_144, Accession, Treatment, group) %>% summarise(avg_height = mean(cor_height_above_reference, na.rm = TRUE))
 

heatavgheight_tp4 <- aggregate(data = avgheightheat[avgheightheat$group == "44/45",], avg_height~Accession+Treatment, FUN= function(i)mean(i))
head(heatavgheight_tp4)

table(heatavgheight_tp4$Accession, heatavgheight_tp4$Treatment)
heatavgheight_tp4[heatavgheight_tp4 == "WW"] <- "heat"
heatavgheight_tp4[heatavgheight_tp4 == "WL"] <- "heat.drought"

# make df with only WW data from heat and control temps
heatWW_avgheight_tp4 <- heatavgheight_tp4[heatavgheight_tp4$Treatment == "heat",]

droughtWW_avgheight_tp4 <- avgheight_tp4[avgheight_tp4$Treatment == "control",]

WW_avgheight_tp4 <- rbind(heatWW_avgheight_tp4, droughtWW_avgheight_tp4)

#check to make sure all accessions have data for both treatments on last day
table(WW_avgheight_tp4$Accession, WW_avgheight_tp4$Treatment)

height_dif2 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgheight_tp4, WW_avgheight_tp4$Accession), function(i)diff(i$avg_height)))), c("height_dif2"))
height_dif2$Accession <- row.names(height_dif2)
head(height_dif2)
row.names(height_dif2) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgheight_tp4 <- heatavgheight_tp4[heatavgheight_tp4$Treatment == "heat.drought",]

control_and_heatWL_avgheight_tp4 <- rbind(heatWL_avgheight_tp4, droughtWW_avgheight_tp4)

#check to make sure all accessions have data for both treatments on last day
table(control_and_heatWL_avgheight_tp4$Accession, control_and_heatWL_avgheight_tp4$Treatment)

height_dif3 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgheight_tp4, control_and_heatWL_avgheight_tp4$Accession), function(i)diff(i$avg_height)))), c("height_dif3"))
height_dif3$Accession <- row.names(height_dif3)
head(height_dif3)
row.names(height_dif3) <- NULL

# create phenotype (Y) file for height
temporary <- join(height_dif, height_dif2, by="Accession")
height_phenotype_tp4 <- join(temporary, height_dif3, by="Accession")
height_phenotype_tp4 <- drop_na(height_phenotype_tp4)

#organize it
height_phenotype_tp4 <- dplyr::select(height_phenotype_tp4, Accession, height_dif, height_dif2, height_dif3)
# rename columns
colnames(height_phenotype_tp4)[2] <- "Drought"
colnames(height_phenotype_tp4)[3] <- "Heat"
colnames(height_phenotype_tp4)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(height_phenotype_tp4, "./height_phenotype_tp4_GAPIT.txt")


#### After making phenotype file once, can start here:
#read in phenotype file
height_phenotype_tp4 <- read.table("./height_phenotype_tp4_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(height_phenotype_tp4)

#check phenotype file
hist(height_phenotype_tp4$Drought)
mean(height_phenotype_tp4$Drought)
range(height_phenotype_tp4$Drought)
sd(height_phenotype_tp4$Drought)
which(is.na(height_phenotype_tp4$Drought))

hist(height_phenotype_tp4$Heat)
mean(height_phenotype_tp4$Heat)
range(height_phenotype_tp4$Heat)
sd(height_phenotype_tp4$Heat)
which(is.na(height_phenotype_tp4$Heat))

hist(height_phenotype_tp4$Heat.Drought)
mean(height_phenotype_tp4$Heat.Drought)
range(height_phenotype_tp4$Heat.Drought)
sd(height_phenotype_tp4$Heat.Drought)
which(is.na(height_phenotype_tp4$Heat.Drought))

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(height_phenotype_tp4$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap_og as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp4")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp4")
#run GAPIT
height_gapit <- GAPIT(
  Y = height_phenotype_tp4,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU"),
  Multiple_analysis = TRUE
)

g

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/test2")

kinship_pheno <- select(accessions_with_gbs, Accession, X)
kinship_pheno <- cbind(kinship_pheno, new)

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(kinship_pheno$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap_og as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

#setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/height_dif_tp4")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

#run GAPIT
test_gapit <- GAPIT(
  Y = kinship_pheno,
  G = hapmap,
  model = c("MLM"),
  Multiple_analysis = TRUE,
  kinship.cluster = "average",
  kinship.group = "Mean",
  Inter.Plot=T
)



################################################################
### GAPIT GWAS for Area (difference)
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp4")

# Area Final Day: control and drought
####
avgarea <- group_by(DFs_144, Accession, Treatment, group) %>% summarise(avg_area = mean(cor_area, na.rm = TRUE))

avgarea_tp4 <- aggregate(data = avgarea[avgarea$group == "44/45",], avg_area~Accession+Treatment, FUN= function(i)mean(i))
head(avgarea_tp4)

# check to make sure each phenotype and treatment has data
table(avgarea_tp4$Accession, avgarea_tp4$Treatment)
avgarea_tp4[avgarea_tp4 == "WW"] <- "control"
avgarea_tp4[avgarea_tp4 == "WL"] <- "drought"

#make df with only WL data from drought experiment (control temp)
droughtWL_avgarea_tp4 <- avgarea_tp4[avgarea_tp4$Treatment == "drought",]

area_dif <- setNames(data.frame(do.call("rbind",lapply(split(avgarea_tp4, avgarea_tp4$Accession), function(i)diff(i$avg_area)))), c("area_dif"))
area_dif$Accession <- row.names(area_dif)
head(area_dif)
row.names(area_dif) <- NULL


####
#control and heat
####
avgareaheat <- group_by(heatDFs_144, Accession, Treatment, group) %>% summarise(avg_area = mean(cor_area, na.rm = TRUE))
 

heatavgarea_tp4 <- aggregate(data = avgareaheat[avgareaheat$group == "44/45",], avg_area~Accession+Treatment, FUN= function(i)mean(i))
head(heatavgarea_tp4)

table(heatavgarea_tp4$Accession, heatavgarea_tp4$Treatment)
heatavgarea_tp4[heatavgarea_tp4 == "WW"] <- "heat"
heatavgarea_tp4[heatavgarea_tp4 == "WL"] <- "heat.drought"

# make df with only WW data from heat and control temps
heatWW_avgarea_tp4 <- heatavgarea_tp4[heatavgarea_tp4$Treatment == "heat",]

droughtWW_avgarea_tp4 <- avgarea_tp4[avgarea_tp4$Treatment == "control",]

WW_avgarea_tp4 <- rbind(heatWW_avgarea_tp4, droughtWW_avgarea_tp4)

#check to make sure all accessions have data for both treatments on last day
table(WW_avgarea_tp4$Accession, WW_avgarea_tp4$Treatment)

area_dif2 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgarea_tp4, WW_avgarea_tp4$Accession), function(i)diff(i$avg_area)))), c("area_dif2"))
area_dif2$Accession <- row.names(area_dif2)
head(area_dif2)
row.names(area_dif2) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgarea_tp4 <- heatavgarea_tp4[heatavgarea_tp4$Treatment == "heat.drought",]

control_and_heatWL_avgarea_tp4 <- rbind(heatWL_avgarea_tp4, droughtWW_avgarea_tp4)

#check to make sure all accessions have data for both treatments on last day
table(control_and_heatWL_avgarea_tp4$Accession, control_and_heatWL_avgarea_tp4$Treatment)

area_dif3 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgarea_tp4, control_and_heatWL_avgarea_tp4$Accession), function(i)diff(i$avg_area)))), c("area_dif3"))
area_dif3$Accession <- row.names(area_dif3)
head(area_dif3)
row.names(area_dif3) <- NULL

# create phenotype (Y) file for area
temporary <- join(area_dif, area_dif2, by="Accession")
area_phenotype_tp4 <- join(temporary, area_dif3, by="Accession")
area_phenotype_tp4 <- drop_na(area_phenotype_tp4)

#organize it
area_phenotype_tp4 <- dplyr::select(area_phenotype_tp4, Accession, area_dif, area_dif2, area_dif3)
# rename columns
colnames(area_phenotype_tp4)[2] <- "Drought"
colnames(area_phenotype_tp4)[3] <- "Heat"
colnames(area_phenotype_tp4)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(area_phenotype_tp4, "./area_phenotype_tp4_GAPIT.txt")



#### After making phenotype file once, can start here:
#read in phenotype file
area_phenotype_tp4 <- read.table("./area_phenotype_tp4_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(area_phenotype_tp4)

#check phenotype file
hist(area_phenotype_tp4$Drought)
mean(area_phenotype_tp4$Drought)
range(area_phenotype_tp4$Drought)
sd(area_phenotype_tp4$Drought)
which(is.na(area_phenotype_tp4$Drought))

hist(area_phenotype_tp4$Heat)
mean(area_phenotype_tp4$Heat)
range(area_phenotype_tp4$Heat)
sd(area_phenotype_tp4$Heat)
which(is.na(area_phenotype_tp4$Heat))

hist(area_phenotype_tp4$Heat.Drought)
mean(area_phenotype_tp4$Heat.Drought)
range(area_phenotype_tp4$Heat.Drought)
sd(area_phenotype_tp4$Heat.Drought)
which(is.na(area_phenotype_tp4$Heat.Drought))

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(area_phenotype_tp4$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap_og as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp4")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/area_dif_tp4")
#run GAPIT
area_gapit <- GAPIT(
  Y = area_phenotype_tp4,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU"),
  Multiple_analysis = TRUE
)



################################################################
### GAPIT GWAS for Percent Damage (difference)
################################################################
setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp4")

# Percent Damage Final Day: control and drought
####
avgdamage <- group_by(DFs_144, Accession, Treatment, group) %>% summarise(avg_damage = mean(percent_damage, na.rm = TRUE))
 
avgdamage_tp4 <- aggregate(data = avgdamage[avgdamage$group == "44/45",], avg_damage~Accession+Treatment, FUN= function(i)mean(i))
head(avgdamage_tp4)

# check to make sure each phenotype and treatment has data
table(avgdamage_tp4$Accession, avgdamage_tp4$Treatment)
avgdamage_tp4[avgdamage_tp4 == "WW"] <- "control"
avgdamage_tp4[avgdamage_tp4 == "WL"] <- "drought"

#make df with only WL data from drought experiment (control temp)
droughtWL_avgdamage_tp4 <- avgdamage_tp4[avgdamage_tp4$Treatment == "drought",]

damage_dif <- setNames(data.frame(do.call("rbind",lapply(split(avgdamage_tp4, avgdamage_tp4$Accession), function(i)diff(i$avg_damage)))), c("damage_dif"))
damage_dif$Accession <- row.names(damage_dif)
head(damage_dif)
row.names(damage_dif) <- NULL


####
#control and heat
####
heatDFs_144$percent_damage <- as.numeric(heatDFs_144$percent_damage)
avgdamageheat <- group_by(heatDFs_144, Accession, Treatment, group) %>% summarise(avg_damage = mean(percent_damage, na.rm = TRUE))

heatavgdamage_tp4 <- aggregate(data = avgdamageheat[avgdamageheat$group == "44/45",], avg_damage~Accession+Treatment, FUN= function(i)mean(i))
head(heatavgdamage_tp4)

table(heatavgdamage_tp4$Accession, heatavgdamage_tp4$Treatment)
heatavgdamage_tp4[heatavgdamage_tp4 == "WW"] <- "heat"
heatavgdamage_tp4[heatavgdamage_tp4 == "WL"] <- "heat.drought"

# make df with only WW data from heat and control temps
heatWW_avgdamage_tp4 <- heatavgdamage_tp4[heatavgdamage_tp4$Treatment == "heat",]

droughtWW_avgdamage_tp4 <- avgdamage_tp4[avgdamage_tp4$Treatment == "control",]

WW_avgdamage_tp4 <- rbind(heatWW_avgdamage_tp4, droughtWW_avgdamage_tp4)

#check to make sure all accessions have data for both treatments on last day
table(WW_avgdamage_tp4$Accession, WW_avgdamage_tp4$Treatment)

damage_dif2 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgdamage_tp4, WW_avgdamage_tp4$Accession), function(i)diff(i$avg_damage)))), c("damage_dif2"))
damage_dif2$Accession <- row.names(damage_dif2)
head(damage_dif2)
row.names(damage_dif2) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgdamage_tp4 <- heatavgdamage_tp4[heatavgdamage_tp4$Treatment == "heat.drought",]

control_and_heatWL_avgdamage_tp4 <- rbind(heatWL_avgdamage_tp4, droughtWW_avgdamage_tp4)

#check to make sure all accessions have data for both treatments on last day
table(control_and_heatWL_avgdamage_tp4$Accession, control_and_heatWL_avgdamage_tp4$Treatment)

damage_dif3 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgdamage_tp4, control_and_heatWL_avgdamage_tp4$Accession), function(i)diff(i$avg_damage)))), c("damage_dif3"))
damage_dif3$Accession <- row.names(damage_dif3)
head(damage_dif3)
row.names(damage_dif3) <- NULL

# create phenotype (Y) file for damage
temporary <- join(damage_dif, damage_dif2, by="Accession")
damage_phenotype_tp4 <- join(temporary, damage_dif3, by="Accession")
damage_phenotype_tp4 <- drop_na(damage_phenotype_tp4)

#organize it
damage_phenotype_tp4 <- dplyr::select(damage_phenotype_tp4, Accession, damage_dif, damage_dif2, damage_dif3)
# rename columns
colnames(damage_phenotype_tp4)[2] <- "Drought"
colnames(damage_phenotype_tp4)[3] <- "Heat"
colnames(damage_phenotype_tp4)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(damage_phenotype_tp4, "./damage_phenotype_tp4_GAPIT.txt")


#### After making phenotype file once, can start here:
#read in phenotype file
damage_phenotype_tp4 <- read.table("./damage_phenotype_tp4_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(damage_phenotype_tp4)

#check phenotype file
hist(damage_phenotype_tp4$Drought)
mean(damage_phenotype_tp4$Drought)
range(damage_phenotype_tp4$Drought)
sd(damage_phenotype_tp4$Drought)
which(is.na(damage_phenotype_tp4$Drought))

hist(damage_phenotype_tp4$Heat)
mean(damage_phenotype_tp4$Heat)
range(damage_phenotype_tp4$Heat)
sd(damage_phenotype_tp4$Heat)
which(is.na(damage_phenotype_tp4$Heat))

hist(damage_phenotype_tp4$Heat.Drought)
mean(damage_phenotype_tp4$Heat.Drought)
range(damage_phenotype_tp4$Heat.Drought)
sd(damage_phenotype_tp4$Heat.Drought)
which(is.na(damage_phenotype_tp4$Heat.Drought))

### Create hapmap subset with only accessions with data for this analysis
#change data to "long" format to be able to sort by Accession
hapmap_og$rn <- seq_len(nrow(hapmap_og))
hapmap_og_melt <- reshape2::melt(hapmap_og, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_og_melt)[13] <- "Accession"

#keep only accessions in this analysis
hapmap_og_melt$Accession <- as.character(hapmap_og_melt$Accession)
in_gwas <- unique(damage_phenotype_tp4$Accession)
hapmap_og_melt$compare.geno <- hapmap_og_melt$Accession %in% in_gwas
hapmap_melt <- filter(hapmap_og_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Accession))
hapmap_melt = dplyr::select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Accession)
hapmap <- dplyr::select(hapmap, -rn)

# save hapmap_og as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp4")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("https://zzlab.net/GAPIT/previous/gapit_functions_20220411.txt")

setwd("/Users/eludwig/Google Drive/My Drive/brachy-gapit/damage_dif_tp4")
#run GAPIT
damage_gapit <- GAPIT(
  Y = damage_phenotype_tp4,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU"),
  Multiple_analysis = TRUE
)

