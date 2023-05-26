rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(ggpubr)


### read in fresh weights
setwd("/Users/eludwig/Desktop/brachy/")
ctl_drought_fw <- read.csv("brachy_drought_fwdw.csv")
colnames(ctl_drought_fw)[1] <- "plantbarcode"
colnames(ctl_drought_fw)[2] <- "Accession"

heat_heatdrought_fw <- read.csv("brachy_heatdrought_fwdw.csv")
colnames(heat_heatdrought_fw)[1] <- "plantbarcode"
colnames(heat_heatdrought_fw)[2] <- "Accession"

### read in area values from images and keep only data from last day time point
setwd("/Users/eludwig/Desktop/brachy/")
DFs <- read.csv("brachy_control_VIS_SV_single_data_compiled_currated.csv")
DFs$Treatment[grep("100", DFs$Treatment.2)] <- "Control"
DFs$Treatment[grep("20", DFs$Treatment.2)] <- "Drought"
DFs_final <- filter(DFs, DAP == 44)
DFs_avg <- group_by(DFs_final, Accession, plantbarcode, Treatment, Replicate) %>% summarise(avg_area = mean(cor_area))

heatDFs <- read.csv("brachy_heat_VIS_SV_single_data_compiled_currated.csv")
heatDFs$Treatment[grep("100", heatDFs$Treatment.2)] <- "Heat"
heatDFs$Treatment[grep("20", heatDFs$Treatment.2)] <- "Heat-drought"
heatDFs_final <- filter(heatDFs, DAP == 45)
heatDFs_avg <- group_by(heatDFs_final, Accession, plantbarcode, Treatment, Replicate) %>% summarise(avg_area = mean(cor_area))


# make combined df with area and fw values and keep only plants with both measurements
ctl_drought_combined <- inner_join(ctl_drought_fw, DFs_avg, by = c("plantbarcode", "Accession"))
length(unique(ctl_drought_combined$Accession))

heat_heatdrought_combined <- inner_join(heat_heatdrought_fw, heatDFs_avg, by = c("plantbarcode", "Accession"))
length(unique(heat_heatdrought_combined$Accession))
heat_heatdrought_combined <- select(heat_heatdrought_combined, -dw)

all_area_fw <- rbind(ctl_drought_combined, heat_heatdrought_combined)
length(unique(all_area_fw$Accession))

# remove two plants where fw was entered incorrectly (in images plants look very large but fw is negative or almost 0)
all_area_fw <- subset(all_area_fw, Accession != "BdTR3B" | fw != -0.0014 | Treatment != "Control" | Replicate != 1) 
all_area_fw <- subset(all_area_fw, Accession != "BdTR5N" | fw != 0.0985 | Treatment != "Control" | Replicate != 2) 

# create scatterplot
sp <- ggscatter(all_area_fw, x = "fw", y = "avg_area",
                title = "PlantCV-calculated plant area vs measured\nfresh weight",
                xlab = "Fresh Weight (g)",
                ylab = expression("Plant Area (cm"^2~")"),
                add = "reg.line",  # Add regression line
                add.params = list(color = "blue", fill = "blue"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = .5, label.y = 45)
