library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(forcats)
library(RColorBrewer)
library(gplots)
library(plotly)
library(FactoMineR)
library(factoextra)
library(missMDA)
library(reshape2)
library(lubridate)
library(corrplot)
library(car)
#library(raster)
library("WorldClimTiles")
library("data.table")
library(lme4)
library(parallel)
library(scales)
library(qqman)
library(multtest)
library(compiler)
library(scatterplot3d)
#library(MASS)

#library(lattice)
#library(lmomco)
#library(MASS, warn.conflicts = FALSE)
#library(nlme, warn.conflicts = FALSE)
#library(mvtnorm, warn.conflicts = FALSE)
#library(grid, warn.conflicts = FALSE)


#setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-drought/")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-heatdrought/")

#####################################################################
# read in the control/drought multi-value data #
#####################################################################

#read in hue data csv files
z500hue <- read.csv("drought_SV_z500_hue.csv")
z1500hue <- read.csv("drought_SV_z1500_hue.csv")
z2500hue <- read.csv("drought_SV_z2500_hue.csv")
z3500hue <- read.csv("drought_SV_z3500_hue.csv")

visz500 <- z500hue
visz1500 <- z1500hue
visz2500 <- z2500hue
visz3500 <- z3500hue

#read in nir data csv files
z500nir <- read.csv("drought_SV_z500_nir.csv")
z1500nir <- read.csv("drought_SV_z1500_nir.csv")
z2500nir <- read.csv("drought_SV_z2500_nir.csv")
z3500nir <- read.csv("drought_SV_z3500_nir.csv")

nirz500 <- z500nir
nirz1500 <- z1500nir
nirz2500 <- z2500nir
nirz3500 <- z3500nir


#filter data
visz500 <- select(visz500, image, plantbarcode, timestamp, value, label)
visz1500 <- select(visz1500, image, plantbarcode, timestamp, value, label)
visz2500 <- select(visz2500, image, plantbarcode, timestamp, value, label)
visz3500 <- select(visz3500, image, plantbarcode, timestamp, value, label)
nirz500 <- select(nirz500, image, plantbarcode, timestamp, trait, other, value, label)
nirz1500 <- select(nirz1500, image, plantbarcode, timestamp, trait, other, value, label)
nirz2500 <- select(nirz2500, image, plantbarcode, timestamp, trait, other, value, label)
nirz3500 <- select(nirz3500, image, plantbarcode, timestamp, trait, other, value, label)

#separate date and time
TS <- c("date","time") #make a vector for the column names of each part of the timestamp
visz500 <- separate(visz500, col = timestamp, into = c("date","time"), sep = " ") #separate each component of timestamp
visz1500 <- separate(visz1500, col = timestamp, into = c("date","time"), sep = " ")
visz2500 <- separate(visz2500, col = timestamp, into = c("date","time"), sep = " ")
visz3500 <- separate(visz3500, col = timestamp, into = c("date","time"), sep = " ")
nirz500 <- separate(nirz500, col = timestamp, into = c("date","time"), sep = " ")
nirz1500 <- separate(nirz1500, col = timestamp, into = c("date","time"), sep = " ")
nirz2500 <- separate(nirz2500, col = timestamp, into = c("date","time"), sep = " ")
nirz3500 <- separate(nirz3500, col = timestamp, into = c("date","time"), sep = " ")

#calculate days after first image
startdate <- as.POSIXct("2015-02-23")
visz500$days = NA
visz1500$days = NA
visz2500$days = NA
visz3500$days = NA
visz500$days = as.integer(difftime(visz500$date, startdate, units = "days"))
visz1500$days = as.integer(difftime(visz1500$date, startdate, units = "days"))
visz2500$days = as.integer(difftime(visz2500$date, startdate, units = "days"))
visz3500$days = as.integer(difftime(visz3500$date, startdate, units = "days"))
nirz500$days = NA
nirz1500$days = NA
nirz2500$days = NA
nirz3500$days = NA
nirz500$days = as.integer(difftime(nirz500$date, startdate, units = "days"))
nirz1500$days = as.integer(difftime(nirz1500$date, startdate, units = "days"))
nirz2500$days = as.integer(difftime(nirz2500$date, startdate, units = "days"))
nirz3500$days = as.integer(difftime(nirz3500$date, startdate, units = "days"))

#copy filtered dfs in case something bad happens
visz500.1 <- visz500
visz1500.1 <- visz1500
visz2500.1 <- visz2500
visz3500.1 <- visz3500
nirz500.1 <- nirz500
nirz1500.1 <- nirz1500
nirz2500.1 <- nirz2500
nirz3500.1 <- nirz3500

#read in barcode csv file and keep relevant info
barcodes <- read.csv("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-drought/brachy-drought-barcodes.csv")
barcodes = select(barcodes, Barcode, Genotype.ID, Treatment.2, Replicate)
colnames(barcodes)[2] <- "Genotype"

#merge barcodes with data file for VIS
DF0 <- filter(visz3500, visz3500$days == "1")
DF0 = merge(DF0, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF0 = drop_na(DF0)
DF1 <- filter(visz3500, visz3500$days == "2")
DF1 = merge(DF1, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF1 = drop_na(DF1)
DF2 <- filter(visz3500, visz3500$days == "3")
DF2 = merge(DF2, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF2 = drop_na(DF2)
DF3 <- filter(visz3500, visz3500$days == "4")
DF3 = merge(DF3, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF3 = drop_na(DF3)
DF4 <- filter(visz3500, visz3500$days == "5")
DF4 = merge(DF4, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF4 = drop_na(DF4)
DF5 <- filter(visz3500, visz3500$days == "6")
DF5 = merge(DF5, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF5 = drop_na(DF5)
DF6 <- filter(visz3500, visz3500$days == "7")
DF6 = merge(DF6, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF6 = drop_na(DF6)
DF7 <- filter(visz3500, visz3500$days == "8")
DF7 = merge(DF7, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF7 = drop_na(DF7)
DF8 <- filter(visz3500, visz3500$days == "9")
DF8 = merge(DF8, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF8 = drop_na(DF8)
DF9 <- filter(visz3500, visz3500$days == "10")
DF9 = merge(DF9, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF9 = drop_na(DF9)
DF10 <- filter(visz2500, visz2500$days == "11")
DF10 = merge(DF10, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF10 = drop_na(DF10)
DF11 <- filter(visz2500, visz2500$days == "12")
DF11 = merge(DF11, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF11 = drop_na(DF11)
DF12 <- filter(visz2500, visz2500$days == "13")
DF12 = merge(DF12, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF12 = drop_na(DF12)
DF13 <- filter(visz2500, visz2500$days == "14")
DF13 = merge(DF13, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF13 = drop_na(DF13)
DF14 <- filter(visz2500, visz2500$days == "15")
DF14 = merge(DF14, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF14 = drop_na(DF14)
DF15 <- filter(visz2500, visz2500$days == "16")
DF15 = merge(DF15, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF15 = drop_na(DF15)
DF16 <- filter(visz1500, visz1500$days == "17")
DF16 = merge(DF16, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF16 = drop_na(DF16)
DF17 <- filter(visz1500, visz1500$days == "18")
DF17 = merge(DF17, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF17 = drop_na(DF17)
DF18 <- filter(visz1500, visz1500$days == "19")
DF18 = merge(DF18, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF18 = drop_na(DF18)
DF19 <- filter(visz1500, visz1500$days == "20")
DF19 = merge(DF19, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF19 = drop_na(DF19)
DF20 <- filter(visz1500, visz1500$days == "21")
DF20 = merge(DF20, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF20 = drop_na(DF20)
DF21 <- filter(visz1500, visz1500$days == "22")
DF21 = merge(DF21, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF21 = drop_na(DF21)
DF22 <- filter(visz1500, visz1500$days == "23")
DF22 = merge(DF22, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF22 = drop_na(DF22)
DF23 <- filter(visz500, visz500$days == "24")
DF23 = merge(DF23, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF23 = drop_na(DF23)
DF24 <- filter(visz500, visz500$days == "25")
DF24 = merge(DF24, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF24 = drop_na(DF24)
DF25 <- filter(visz500, visz500$days == "26")
DF25 = merge(DF25, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF25 = drop_na(DF25)
DF26 <- filter(visz500, visz500$days == "27")
DF26 = merge(DF26, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF26 = drop_na(DF26)
DF27 <- filter(visz500, visz500$days == "28")
DF27 = merge(DF27, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF27 = drop_na(DF27)
DF28 <- filter(visz500, visz500$days == "29")
DF28 = merge(DF28, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF28 = drop_na(DF28)
DF29 <- filter(visz500, visz500$days == "30")
DF29 = merge(DF29, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF29 = drop_na(DF29)
DF30 <- filter(visz500, visz500$days == "31")
DF30 = merge(DF30, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF30 = drop_na(DF30)
DF31 <- filter(visz500, visz500$days == "32")
DF31 = merge(DF31, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF31 = drop_na(DF31)

#merge all dfs and then remove
DF <- bind_rows(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31)
rm(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31)

colnames(DF)[7] <- "Genotype"

#save the csv file of the compiled VIS data
write.csv(DF, "./brachy_VIS_SV_multi_data_compiled_currated.csv", row.names = FALSE)


#merge barcodes with data file for NIR
DF0 <- filter(nirz3500, nirz3500$days == "1")
DF0 = merge(DF0, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF0 = drop_na(DF0)
DF1 <- filter(nirz3500, nirz3500$days == "2")
DF1 = merge(DF1, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF1 = drop_na(DF1)
DF2 <- filter(nirz3500, nirz3500$days == "3")
DF2 = merge(DF2, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF2 = drop_na(DF2)
DF3 <- filter(nirz3500, nirz3500$days == "4")
DF3 = merge(DF3, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF3 = drop_na(DF3)
DF4 <- filter(nirz3500, nirz3500$days == "5")
DF4 = merge(DF4, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF4 = drop_na(DF4)
DF5 <- filter(nirz3500, nirz3500$days == "6")
DF5 = merge(DF5, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF5 = drop_na(DF5)
DF6 <- filter(nirz3500, nirz3500$days == "7")
DF6 = merge(DF6, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF6 = drop_na(DF6)
DF7 <- filter(nirz3500, nirz3500$days == "8")
DF7 = merge(DF7, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF7 = drop_na(DF7)
DF8 <- filter(nirz3500, nirz3500$days == "9")
DF8 = merge(DF8, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF8 = drop_na(DF8)
DF9 <- filter(nirz3500, nirz3500$days == "10")
DF9 = merge(DF9, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF9 = drop_na(DF9)
DF10 <- filter(nirz2500, nirz2500$days == "11")
DF10 = merge(DF10, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF10 = drop_na(DF10)
DF11 <- filter(nirz2500, nirz2500$days == "12")
DF11 = merge(DF11, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF11 = drop_na(DF11)
DF12 <- filter(nirz2500, nirz2500$days == "13")
DF12 = merge(DF12, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF12 = drop_na(DF12)
DF13 <- filter(nirz2500, nirz2500$days == "14")
DF13 = merge(DF13, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF13 = drop_na(DF13)
DF14 <- filter(nirz2500, nirz2500$days == "15")
DF14 = merge(DF14, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF14 = drop_na(DF14)
DF15 <- filter(nirz2500, nirz2500$days == "16")
DF15 = merge(DF15, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF15 = drop_na(DF15)
DF16 <- filter(nirz1500, nirz1500$days == "17")
DF16 = merge(DF16, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF16 = drop_na(DF16)
DF17 <- filter(nirz1500, nirz1500$days == "18")
DF17 = merge(DF17, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF17 = drop_na(DF17)
DF18 <- filter(nirz1500, nirz1500$days == "19")
DF18 = merge(DF18, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF18 = drop_na(DF18)
DF19 <- filter(nirz1500, nirz1500$days == "20")
DF19 = merge(DF19, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF19 = drop_na(DF19)
DF20 <- filter(nirz1500, nirz1500$days == "21")
DF20 = merge(DF20, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF20 = drop_na(DF20)
DF21 <- filter(nirz1500, nirz1500$days == "22")
DF21 = merge(DF21, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF21 = drop_na(DF21)
DF22 <- filter(nirz1500, nirz1500$days == "23")
DF22 = merge(DF22, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF22 = drop_na(DF22)
DF23 <- filter(nirz500, nirz500$days == "24")
DF23 = merge(DF23, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF23 = drop_na(DF23)
DF24 <- filter(nirz500, nirz500$days == "25")
DF24 = merge(DF24, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF24 = drop_na(DF24)
DF25 <- filter(nirz500, nirz500$days == "26")
DF25 = merge(DF25, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF25 = drop_na(DF25)
DF26 <- filter(nirz500, nirz500$days == "27")
DF26 = merge(DF26, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF26 = drop_na(DF26)
DF27 <- filter(nirz500, nirz500$days == "28")
DF27 = merge(DF27, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF27 = drop_na(DF27)
DF28 <- filter(nirz500, nirz500$days == "29")
DF28 = merge(DF28, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF28 = drop_na(DF28)
DF29 <- filter(nirz500, nirz500$days == "30")
DF29 = merge(DF29, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF29 = drop_na(DF29)
DF30 <- filter(nirz500, nirz500$days == "31")
DF30 = merge(DF30, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF30 = drop_na(DF30)
DF31 <- filter(nirz500, nirz500$days == "32")
DF31 = merge(DF31, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
DF31 = drop_na(DF31)

#merge all dfs and then remove
DFNIR <- bind_rows(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31)
rm(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31)

#save the csv file of the compiled data
write.csv(DFNIR, "./brachy_NIR_SV_multi_data_compiled_currated.csv", row.names = FALSE)




#####################################################################
# read in the control/drought single-value data #
#####################################################################
#read in single value trait csv files
z500s <- read.csv("drought_SV_z500-single-value-traits.csv")
z1500s <- read.csv("drought_SV_z1500-single-value-traits.csv")
z2500s <- read.csv("drought_SV_z2500-single-value-traits.csv")
z3500s <- read.csv("drought_SV_z3500-single-value-traits.csv")


#copy single value dfs in case something happens (so don't have to read in again)
copyz500s <- z500s
copyz1500s <- z1500s
copyz2500s <- z2500s
copyz3500s <- z3500s

#subset VIS data
visz500s = z500s[z500s$imgtype == 'VIS',]
visz1500s = z1500s[z1500s$imgtype == 'VIS',]
visz2500s = z2500s[z2500s$imgtype == 'VIS',]
visz3500s = z3500s[z3500s$imgtype == 'VIS',]

#####################################################################
# Control/drought zoom correction #
#####################################################################
zoom.lm <- lm(zoom.camera ~ zoom, data = data.frame(zoom = c(1, 6000), zoom.camera = c(1, 6)))
# Download data for a reference object imaged at different zoom levels
if (!file.exists('zoom_calibration_data.txt')) {
  download.file('http://files.figshare.com/2084101/zoom_calibration_data.txt',
                'zoom_calibration_data.txt')
}
# Read zoom calibration data
z.data <- read.table(file = "zoom_calibration_data.txt", sep = "\t", header = TRUE)
# Calculate px per cm
z.data$px_cm <- z.data$length_px / z.data$length_cm
# Calculate area for each row
z.data$area_cm <- ifelse(z.data$reference == z.data$reference[[1]], (13.2*13.2), (13.2*3.7))
# Calculate px**2 per cm**2
z.data$px2_cm2 <- z.data$area_px / z.data$area_cm
# Convert LemnaTec zoom units to camera zoom units
z.data$zoom.camera <- predict(object = zoom.lm, newdata = z.data)
# Zoom correction for area
area.coef <- coef(nls(log(px2_cm2) ~ log(a * exp(b * zoom.camera)),
                      z.data, start = c(a = 1, b = 0.01)))
area.coef <- data.frame(a = area.coef[1], b = area.coef[2])
area.nls <- nls(px2_cm2 ~ a * exp(b * zoom.camera),
                data = z.data, start = c(a = area.coef$a, b = area.coef$b))
# Zoom correction for height
height.coef <- coef(nls(log(px_cm) ~ log(a * exp(b * zoom.camera)),
                        z.data, start = c(a = 1, b = 0.01)))
height.coef <- data.frame(a = area.coef[1], b = area.coef[2])
height.nls <- nls(px_cm ~ a * exp(b * zoom.camera),
                  data = z.data, start = c(a = area.coef$a, b = area.coef$b))

# change the zoom vector to the zooms in your dataset
exp_zooms = data.frame(zoom=c(500, 1500, 2500, 3500))
exp_zooms$zoom.camera = predict(object = zoom.lm, newdata = exp_zooms)
# this is the one I used:
exp_zooms$px2_cm2 = predict(object = area.nls, newdata = exp_zooms)
exp_zooms$px_cm = predict(object = height.nls, newdata = exp_zooms)

###################zoom correction
# at z500
visz500s <- transform(visz500s, cor_area = area / 1229.355)
visz500s <- transform(visz500s, cor_area_above_reference = area_above_reference / 1229.355)
visz500s <- transform(visz500s, cor_area_below_reference = area_below_reference / 1229.355)
visz500s <- transform(visz500s, cor_width = width / 33.95782)
visz500s <- transform(visz500s, cor_height = height / 33.95782)
visz500s <- transform(visz500s, cor_height_above_reference = height_above_reference / 33.95782)
visz500s <- transform(visz500s, cor_height_below_reference = height_below_reference / 33.95782)
visz500s <- transform(visz500s, cor_convex_hull_area = convex_hull_area / 1229.355)
visz500s <- transform(visz500s, cor_perimeter = perimeter / 33.95782)
visz500s <- transform(visz500s, cor_longest_path = longest_path / 33.95782)
visz500s <- transform(visz500s, cor_ellipse_major_axis = ellipse_major_axis / 33.95782)
visz500s <- transform(visz500s, cor_ellipse_minor_axis = ellipse_minor_axis / 33.95782)

# at z1500
visz1500s <- transform(visz1500s, cor_area = area / 2610.580)
visz1500s <- transform(visz1500s, cor_area_above_reference = area_above_reference / 2610.580)
visz1500s <- transform(visz1500s, cor_area_below_reference = area_below_reference / 2610.580)
visz1500s <- transform(visz1500s, cor_width = width / 49.41139)
visz1500s <- transform(visz1500s, cor_height = height / 49.41139)
visz1500s <- transform(visz1500s, cor_height_above_reference = height_above_reference / 49.41139)
visz1500s <- transform(visz1500s, cor_height_below_reference = height_below_reference / 49.41139)
visz1500s <- transform(visz1500s, cor_convex_hull_area = convex_hull_area / 2610.580)
visz1500s <- transform(visz1500s, cor_perimeter = perimeter / 49.41139)
visz1500s <- transform(visz1500s, cor_longest_path = longest_path / 49.41139)
visz1500s <- transform(visz1500s, cor_ellipse_major_axis = ellipse_major_axis / 49.41139)
visz1500s <- transform(visz1500s, cor_ellipse_minor_axis = ellipse_minor_axis / 49.41139)

###################zoom correction
# at z2500
visz2500s <- transform(visz2500s, cor_area = area / 5543.662)
visz2500s <- transform(visz2500s, cor_area_above_reference = area_above_reference / 5543.662)
visz2500s <- transform(visz2500s, cor_area_below_reference = area_below_reference / 5543.662)
visz2500s <- transform(visz2500s, cor_width = width / 71.89759)
visz2500s <- transform(visz2500s, cor_height = height / 71.89759)
visz2500s <- transform(visz2500s, cor_height_above_reference = height_above_reference / 71.89759)
visz2500s <- transform(visz2500s, cor_height_below_reference = height_below_reference / 71.89759)
visz2500s <- transform(visz2500s, cor_convex_hull_area = convex_hull_area / 5543.662)
visz2500s <- transform(visz2500s, cor_perimeter = perimeter / 71.89759)
visz2500s <- transform(visz2500s, cor_longest_path = longest_path / 71.89759)
visz2500s <- transform(visz2500s, cor_ellipse_major_axis = ellipse_major_axis / 71.89759)
visz2500s <- transform(visz2500s, cor_ellipse_minor_axis = ellipse_minor_axis / 71.89759)

# at z3500
visz3500s <- transform(visz3500s, cor_area = area / 11772.169)
visz3500s <- transform(visz3500s, cor_area_above_reference = area_above_reference / 11772.169)
visz3500s <- transform(visz3500s, cor_area_below_reference = area_below_reference / 11772.169)
visz3500s <- transform(visz3500s, cor_width = width / 104.61683)
visz3500s <- transform(visz3500s, cor_height = height / 104.61683)
visz3500s <- transform(visz3500s, cor_height_above_reference = height_above_reference / 104.61683)
visz3500s <- transform(visz3500s, cor_height_below_reference = height_below_reference / 104.61683)
visz3500s <- transform(visz3500s, cor_convex_hull_area = convex_hull_area / 11772.169)
visz3500s <- transform(visz3500s, cor_perimeter = perimeter / 104.61683)
visz3500s <- transform(visz3500s, cor_longest_path = longest_path / 104.61683)
visz3500s <- transform(visz3500s, cor_ellipse_major_axis = ellipse_major_axis / 104.61683)
visz3500s <- transform(visz3500s, cor_ellipse_minor_axis = ellipse_minor_axis / 104.61683)

########################################
# Control/drought data filtering and organizing #
########################################

#filter data
visz500s <- select(visz500s, image, plantbarcode, timestamp, hue_circular_mean, hue_circular_std, hue_median, percent_unhealthy, cor_area, cor_width, cor_height, cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, cor_area_below_reference, percent_area_below_reference, cor_convex_hull_area, solidity, cor_perimeter, cor_longest_path, convex_hull_vertices, cor_ellipse_major_axis, cor_ellipse_minor_axis, ellipse_angle, ellipse_eccentricity)
visz1500s <- select(visz1500s, image, plantbarcode, timestamp, hue_circular_mean, hue_circular_std, hue_median, percent_unhealthy, cor_area, cor_width, cor_height, cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, cor_area_below_reference, percent_area_below_reference, cor_convex_hull_area, solidity, cor_perimeter, cor_longest_path, convex_hull_vertices, cor_ellipse_major_axis, cor_ellipse_minor_axis, ellipse_angle, ellipse_eccentricity)
visz2500s <- select(visz2500s, image, plantbarcode, timestamp, hue_circular_mean, hue_circular_std, hue_median, percent_unhealthy, cor_area, cor_width, cor_height, cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, cor_area_below_reference, percent_area_below_reference, cor_convex_hull_area, solidity, cor_perimeter, cor_longest_path, convex_hull_vertices, cor_ellipse_major_axis, cor_ellipse_minor_axis, ellipse_angle, ellipse_eccentricity)
visz3500s <- select(visz3500s, image, plantbarcode, timestamp, hue_circular_mean, hue_circular_std, hue_median, percent_unhealthy, cor_area, cor_width, cor_height, cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, cor_area_below_reference, percent_area_below_reference, cor_convex_hull_area, solidity, cor_perimeter, cor_longest_path, convex_hull_vertices, cor_ellipse_major_axis, cor_ellipse_minor_axis, ellipse_angle, ellipse_eccentricity)

#separate date and time
TS <- c("date","time") #make a vector for the column names of each part of the timestamp
visz500s <- separate(visz500s, col = timestamp, into = c("date","time"), sep = " ") #separate each component of timestamp
visz1500s <- separate(visz1500s, col = timestamp, into = c("date","time"), sep = " ")
visz2500s <- separate(visz2500s, col = timestamp, into = c("date","time"), sep = " ")
visz3500s <- separate(visz3500s, col = timestamp, into = c("date","time"), sep = " ")

#calculate days after first image
startdate <- as.POSIXct("2015-02-23")
visz500s$days = NA
visz1500s$days = NA
visz2500s$days = NA
visz3500s$days = NA
visz500s$days = as.integer(difftime(visz500s$date, startdate, units = "days"))
visz1500s$days = as.integer(difftime(visz1500s$date, startdate, units = "days"))
visz2500s$days = as.integer(difftime(visz2500s$date, startdate, units = "days"))
visz3500s$days = as.integer(difftime(visz3500s$date, startdate, units = "days"))

#copy filtered dfs in case something bad happens
visz500.1s <- visz500s
visz1500.1s <- visz1500s
visz2500.1s <- visz2500s
visz3500.1s <- visz3500s

#read in barcode csv file and keep relevant info
barcodes <- read.csv("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-drought/brachy-drought-barcodes.csv")
barcodes = select(barcodes, Barcode, Genotype, Treatment.2, Replicate)
colnames(barcodes)[2] <- "Genotype"

#merge barcodes with data file
DF0 <- filter(visz3500s, visz3500s$days == "1")
DF0 = merge(DF0, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF0 = drop_na(DF0)
DF1 <- filter(visz3500s, visz3500s$days == "2")
DF1 = merge(DF1, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF1 = drop_na(DF1)
DF2 <- filter(visz3500s, visz3500s$days == "3")
DF2 = merge(DF2, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF2 = drop_na(DF2)
DF3 <- filter(visz3500s, visz3500s$days == "4")
DF3 = merge(DF3, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF3 = drop_na(DF3)
DF4 <- filter(visz3500s, visz3500s$days == "5")
DF4 = merge(DF4, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF4 = drop_na(DF4)
DF5 <- filter(visz3500s, visz3500s$days == "6")
DF5 = merge(DF5, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF5 = drop_na(DF5)
DF6 <- filter(visz3500s, visz3500s$days == "7")
DF6 = merge(DF6, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF6 = drop_na(DF6)
DF7 <- filter(visz3500s, visz3500s$days == "8")
DF7 = merge(DF7, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF7 = drop_na(DF7)
DF8 <- filter(visz3500s, visz3500s$days == "9")
DF8 = merge(DF8, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF8 = drop_na(DF8)
DF9 <- filter(visz3500s, visz3500s$days == "10")
DF9 = merge(DF9, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF9 = drop_na(DF9)
DF10 <- filter(visz2500s, visz2500s$days == "11")
DF10 = merge(DF10, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF10 = drop_na(DF10)
DF11 <- filter(visz2500s, visz2500s$days == "12")
DF11 = merge(DF11, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF11 = drop_na(DF11)
DF12 <- filter(visz2500s, visz2500s$days == "13")
DF12 = merge(DF12, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF12 = drop_na(DF12)
DF13 <- filter(visz2500s, visz2500s$days == "14")
DF13 = merge(DF13, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF13 = drop_na(DF13)
DF14 <- filter(visz2500s, visz2500s$days == "15")
DF14 = merge(DF14, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF14 = drop_na(DF14)
DF15 <- filter(visz2500s, visz2500s$days == "16")
DF15 = merge(DF15, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF15 = drop_na(DF15)
DF16 <- filter(visz1500s, visz1500s$days == "17")
DF16 = merge(DF16, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF16 = drop_na(DF16)
DF17 <- filter(visz1500s, visz1500s$days == "18")
DF17 = merge(DF17, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF17 = drop_na(DF17)
DF18 <- filter(visz1500s, visz1500s$days == "19")
DF18 = merge(DF18, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF18 = drop_na(DF18)
DF19 <- filter(visz1500s, visz1500s$days == "20")
DF19 = merge(DF19, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF19 = drop_na(DF19)
DF20 <- filter(visz1500s, visz1500s$days == "21")
DF20 = merge(DF20, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF20 = drop_na(DF20)
DF21 <- filter(visz1500s, visz1500s$days == "22")
DF21 = merge(DF21, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF21 = drop_na(DF21)
DF22 <- filter(visz1500s, visz1500s$days == "23")
DF22 = merge(DF22, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF22 = drop_na(DF22)
DF23 <- filter(visz500s, visz500s$days == "24")
DF23 = merge(DF23, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF23 = drop_na(DF23)
DF24 <- filter(visz500s, visz500s$days == "25")
DF24 = merge(DF24, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF24 = drop_na(DF24)
DF25 <- filter(visz500s, visz500s$days == "26")
DF25 = merge(DF25, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF25 = drop_na(DF25)
DF26 <- filter(visz500s, visz500s$days == "27")
DF26 = merge(DF26, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF26 = drop_na(DF26)
DF27 <- filter(visz500s, visz500s$days == "28")
DF27 = merge(DF27, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF27 = drop_na(DF27)
DF28 <- filter(visz500s, visz500s$days == "29")
DF28 = merge(DF28, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF28 = drop_na(DF28)
DF29 <- filter(visz500s, visz500s$days == "30")
DF29 = merge(DF29, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF29 = drop_na(DF29)
DF30 <- filter(visz500s, visz500s$days == "31")
DF30 = merge(DF30, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF30 = drop_na(DF30)
DF31 <- filter(visz500s, visz500s$days == "32")
DF31 = merge(DF31, barcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF31 = drop_na(DF31)

#merge all dfs and then remove
DFs <- bind_rows(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31)
rm(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31)

colnames(DFs)[28] <- "Genotype"
length(unique(DFs$Genotype))

#read in locations csv file and keep relevant info
locations <- read.csv("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-drought/brachy-locations.csv")
locations = select(locations, Genotype, Collection_location, Latitude, Longitude, Elevation)

#merge locations file with data file
DF0 <- filter(DFs, DFs$days == "1")
DF0 = merge(DF0, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF0 = drop_na(DF0)
DF1 <- filter(DFs, DFs$days == "2")
DF1 = merge(DF1, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF1 = drop_na(DF1)
DF2 <- filter(DFs, DFs$days == "3")
DF2 = merge(DF2, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF2 = drop_na(DF2)
DF3 <- filter(DFs, DFs$days == "4")
DF3 = merge(DF3, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF3 = drop_na(DF3)
DF4 <- filter(DFs, DFs$days == "5")
DF4 = merge(DF4, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF4 = drop_na(DF4)
DF5 <- filter(DFs, DFs$days == "6")
DF5 = merge(DF5, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF5 = drop_na(DF5)
DF6 <- filter(DFs, DFs$days == "7")
DF6 = merge(DF6, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF6 = drop_na(DF6)
DF7 <- filter(DFs, DFs$days == "8")
DF7 = merge(DF7, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF7 = drop_na(DF7)
DF8 <- filter(DFs, DFs$days == "9")
DF8 = merge(DF8, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF8 = drop_na(DF8)
DF9 <- filter(DFs, DFs$days == "10")
DF9 = merge(DF9, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF9 = drop_na(DF9)
DF10 <- filter(DFs, DFs$days == "11")
DF10 = merge(DF10, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF10 = drop_na(DF10)
DF11 <- filter(DFs, DFs$days == "12")
DF11 = merge(DF11, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF11 = drop_na(DF11)
DF12 <- filter(DFs, DFs$days == "13")
DF12 = merge(DF12, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF12 = drop_na(DF12)
DF13 <- filter(DFs, DFs$days == "14")
DF13 = merge(DF13, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF13 = drop_na(DF13)
DF14 <- filter(DFs, DFs$days == "15")
DF14 = merge(DF14, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF14 = drop_na(DF14)
DF15 <- filter(DFs, DFs$days == "16")
DF15 = merge(DF15, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF15 = drop_na(DF15)
DF16 <- filter(DFs, DFs$days == "17")
DF16 = merge(DF16, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF16 = drop_na(DF16)
DF17 <- filter(DFs, DFs$days == "18")
DF17 = merge(DF17, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF17 = drop_na(DF17)
DF18 <- filter(DFs, DFs$days == "19")
DF18 = merge(DF18, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF18 = drop_na(DF18)
DF19 <- filter(DFs, DFs$days == "20")
DF19 = merge(DF19, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF19 = drop_na(DF19)
DF20 <- filter(DFs, DFs$days == "21")
DF20 = merge(DF20, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF20 = drop_na(DF20)
DF21 <- filter(DFs, DFs$days == "22")
DF21 = merge(DF21, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF21 = drop_na(DF21)
DF22 <- filter(DFs, DFs$days == "23")
DF22 = merge(DF22, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF22 = drop_na(DF22)
DF23 <- filter(DFs, DFs$days == "24")
DF23 = merge(DF23, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF23 = drop_na(DF23)
DF24 <- filter(DFs, DFs$days == "25")
DF24 = merge(DF24, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF24 = drop_na(DF24)
DF25 <- filter(DFs, DFs$days == "26")
DF25 = merge(DF25, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF25 = drop_na(DF25)
DF26 <- filter(DFs, DFs$days == "27")
DF26 = merge(DF26, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF26 = drop_na(DF26)
DF27 <- filter(DFs, DFs$days == "28")
DF27 = merge(DF27, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF27 = drop_na(DF27)
DF28 <- filter(DFs, DFs$days == "29")
DF28 = merge(DF28, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF28 = drop_na(DF28)
DF29 <- filter(DFs, DFs$days == "30")
DF29 = merge(DF29, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF29 = drop_na(DF29)
DF30 <- filter(DFs, DFs$days == "31")
DF30 = merge(DF30, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF30 = drop_na(DF30)
DF31 <- filter(DFs, DFs$days == "32")
DF31 = merge(DF31, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF31 = drop_na(DF31)

#merge all dfs and then remove
DFs <- bind_rows(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31)
rm(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31)

length(unique(DFs$Genotype))
#DFs <- select(DFs, image, plantbarcode, Genotype, date, time, days, Treatment.2, Replicate, Collection_location, Latitude, Longitude, Elevation, hue_circular_mean, hue_circular_std, hue_median, percent_unhealthy, cor_area, cor_width, cor_height, cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, cor_area_below_reference, percent_area_below_reference)

DFs <- select(DFs, image, plantbarcode, Genotype, date, time, days, Treatment.2, Replicate, Collection_location, Latitude, Longitude, Elevation, hue_circular_mean, hue_circular_std, hue_median, percent_unhealthy, cor_area, cor_width, cor_height, cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, cor_area_below_reference, percent_area_below_reference, cor_convex_hull_area, solidity, cor_perimeter, cor_longest_path, convex_hull_vertices, cor_ellipse_major_axis, cor_ellipse_minor_axis, ellipse_angle, ellipse_eccentricity)


DFs$imgdays = NA
DFs$imgdays[grep("1", DFs$days)] <- "1"
DFs$imgdays[grep("2", DFs$days)] <- "1"
DFs$imgdays[grep("3", DFs$days)] <- "2"
DFs$imgdays[grep("4", DFs$days)] <- "2"
DFs$imgdays[grep("5", DFs$days)] <- "3"
DFs$imgdays[grep("6", DFs$days)] <- "3"
DFs$imgdays[grep("7", DFs$days)] <- "4"
DFs$imgdays[grep("8", DFs$days)] <- "4"
DFs$imgdays[grep("9", DFs$days)] <- "5"
DFs$imgdays[grep("10", DFs$days)] <- "5"
DFs$imgdays[grep("11", DFs$days)] <- "6"
DFs$imgdays[grep("12", DFs$days)] <- "6"
DFs$imgdays[grep("13", DFs$days)] <- "7"
DFs$imgdays[grep("14", DFs$days)] <- "7"
DFs$imgdays[grep("15", DFs$days)] <- "8"
DFs$imgdays[grep("16", DFs$days)] <- "8"
DFs$imgdays[grep("17", DFs$days)] <- "9"
DFs$imgdays[grep("18", DFs$days)] <- "9"
DFs$imgdays[grep("19", DFs$days)] <- "10"
DFs$imgdays[grep("20", DFs$days)] <- "10"
DFs$imgdays[grep("21", DFs$days)] <- "11"
DFs$imgdays[grep("22", DFs$days)] <- "11"
DFs$imgdays[grep("23", DFs$days)] <- "12"
DFs$imgdays[grep("24", DFs$days)] <- "12"
DFs$imgdays[grep("25", DFs$days)] <- "13"
DFs$imgdays[grep("26", DFs$days)] <- "13"
DFs$imgdays[grep("27", DFs$days)] <- "14"
DFs$imgdays[grep("28", DFs$days)] <- "14"
DFs$imgdays[grep("29", DFs$days)] <- "15"
DFs$imgdays[grep("30", DFs$days)] <- "15"

DFs$Treatment = NA
DFs$Treatment[grep("20", DFs$Treatment.2)] <- "WL"
DFs$Treatment[grep("100", DFs$Treatment.2)] <- "WW"


# multiply the % unhealthy pixels by 100 to convert from fraction to percent
DFs[,16] <- DFs[,16]*100

#save the csv file of the compiled VIS data
write.csv(DFs, "./brachy_VIS_SV_single_data_compiled_currated.csv", row.names = FALSE)

#keep only genotypes that are in both experiments and have complete data
DFs$compare.geno <- DFs$Genotype %in% list3
DFs <- filter(DFs, compare.geno == TRUE)
length(unique(DFs$Genotype))
DFs = select(DFs, -compare.geno)

#write csv with data from only genotypes that are in both experiments
write.csv(DFs, "./brachy_control_VIS_SV_single_data_limited_genotypes.csv", row.names = FALSE)


#####################################################################
# Read in the heat/heat-drought multi-value data #
#####################################################################

#read in hue data csv files
heatz2500hue <- read.csv("heatdrought_SV_z2500_hue.csv")
heatz3500hue <- read.csv("heatdrought_SV_z3500_hue.csv")

visz2500heat <- heatz2500hue
visz3500heat <- heatz3500hue

#filter data
visz2500heat <- select(visz2500heat, image, plantbarcode, timestamp, value, label)
visz3500heat <- select(visz3500heat, image, plantbarcode, timestamp, value, label)

#separate date and time
TS <- c("date","time") #make a vector for the column names of each part of the timestamp
visz2500heat <- separate(visz2500heat, col = timestamp, into = c("date","time"), sep = " ")
visz3500heat <- separate(visz3500heat, col = timestamp, into = c("date","time"), sep = " ")

#change date format to R standard for z2500
visz2500heat$date <- as.Date(visz2500heat$date, format="%m/%d/%y")

# split both z3500 dfs in half to separate out 2nd half that has incorrect date format
newDF <- visz3500heat[1:3963251,]
newDF2 <- visz3500heat[3963252:7926502,]

# Change date format to R standard for z3500
newDF$date <- as.Date(newDF$date)
newDF2$date <- as.Date(newDF2$date, format="%m/%d/%y")

#recombine z3500 dfs with correct date format 
newDF$date <- as.character(newDF$date)
newDF2$date <- as.character(newDF2$date)
visz3500heat <- bind_rows(newDF, newDF2)

visz2500heat$date <- as.character(visz2500heat$date)
visz3500heat$date <- as.character(visz3500heat$date)

#calculate days after first image
startdate1 <- as.POSIXct("2014-08-12")
visz2500heat$days = NA
visz3500heat$days = NA
visz2500heat$days = as.integer(difftime(visz2500heat$date, startdate1, units = "days"))
visz3500heat$days = as.integer(difftime(visz3500heat$date, startdate1, units = "days"))

#copy filtered dfs in case something bad happens
visz2500heat.1 <- visz2500heat
visz3500heat.1 <- visz3500heat

#read in barcode csv file and keep relevant info
heatbarcodes <- read.csv("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-heatdrought/brachy-heatdrought-barcodes.csv")
heatbarcodes = select(heatbarcodes, Barcode, Genotype, Treatment.2, Replicate)
colnames(heatbarcodes)[2] <- "Genotype"

#merge barcodes with data file for VIS
DF0 <- filter(visz3500heat, visz3500heat$days == "1")
DF0 = merge(DF0, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF0 = drop_na(DF0)
DF1 <- filter(visz3500heat, visz3500heat$days == "2")
DF1 = merge(DF1, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF1 = drop_na(DF1)
DF2 <- filter(visz3500heat, visz3500heat$days == "3")
DF2 = merge(DF2, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF2 = drop_na(DF2)
DF3 <- filter(visz3500heat, visz3500heat$days == "4")
DF3 = merge(DF3, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF3 = drop_na(DF3)
DF4 <- filter(visz3500heat, visz3500heat$days == "5")
DF4 = merge(DF4, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF4 = drop_na(DF4)
DF5 <- filter(visz3500heat, visz3500heat$days == "6")
DF5 = merge(DF5, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF5 = drop_na(DF5)
DF6 <- filter(visz3500heat, visz3500heat$days == "7")
DF6 = merge(DF6, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF6 = drop_na(DF6)
DF7 <- filter(visz3500heat, visz3500heat$days == "8")
DF7 = merge(DF7, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF7 = drop_na(DF7)
DF8 <- filter(visz3500heat, visz3500heat$days == "9")
DF8 = merge(DF8, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF8 = drop_na(DF8)
DF9 <- filter(visz3500heat, visz3500heat$days == "10")
DF9 = merge(DF9, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF9 = drop_na(DF9)
DF10 <- filter(visz2500heat, visz2500heat$days == "11")
DF10 = merge(DF10, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF10 = drop_na(DF10)
DF11 <- filter(visz2500heat, visz2500heat$days == "12")
DF11 = merge(DF11, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF11 = drop_na(DF11)
DF12 <- filter(visz2500heat, visz2500heat$days == "13")
DF12 = merge(DF12, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF12 = drop_na(DF12)
DF13 <- filter(visz2500heat, visz2500heat$days == "14")
DF13 = merge(DF13, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF13 = drop_na(DF13)
DF14 <- filter(visz2500heat, visz2500heat$days == "15")
DF14 = merge(DF14, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF14 = drop_na(DF14)
DF15 <- filter(visz2500heat, visz2500heat$days == "16")
DF15 = merge(DF15, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF15 = drop_na(DF15)
DF16 <- filter(visz2500heat, visz2500heat$days == "17")
DF16 = merge(DF16, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF16 = drop_na(DF16)
DF17 <- filter(visz2500heat, visz2500heat$days == "18")
DF17 = merge(DF17, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF17 = drop_na(DF17)
DF18 <- filter(visz2500heat, visz2500heat$days == "19")
DF18 = merge(DF18, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF18 = drop_na(DF18)
DF19 <- filter(visz2500heat, visz2500heat$days == "20")
DF19 = merge(DF19, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF19 = drop_na(DF19)
DF20 <- filter(visz2500heat, visz2500heat$days == "21")
DF20 = merge(DF20, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF20 = drop_na(DF20)
DF21 <- filter(visz2500heat, visz2500heat$days == "22")
DF21 = merge(DF21, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF21 = drop_na(DF21)
DF22 <- filter(visz2500heat, visz2500heat$days == "23")
DF22 = merge(DF22, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF22 = drop_na(DF22)
DF23 <- filter(visz2500heat, visz2500heat$days == "24")
DF23 = merge(DF23, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF23 = drop_na(DF23)
DF24 <- filter(visz2500heat, visz2500heat$days == "25")
DF24 = merge(DF24, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF24 = drop_na(DF24)
DF25 <- filter(visz2500heat, visz2500heat$days == "26")
DF25 = merge(DF25, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF25 = drop_na(DF25)
DF26 <- filter(visz2500heat, visz2500heat$days == "27")
DF26 = merge(DF26, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF26 = drop_na(DF26)
DF27 <- filter(visz2500heat, visz2500heat$days == "28")
DF27 = merge(DF27, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF27 = drop_na(DF27)
DF28 <- filter(visz2500heat, visz2500heat$days == "29")
DF28 = merge(DF28, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF28 = drop_na(DF28)
DF29 <- filter(visz2500heat, visz2500heat$days == "30")
DF29 = merge(DF29, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF29 = drop_na(DF29)
DF30 <- filter(visz2500heat, visz2500heat$days == "31")
DF30 = merge(DF30, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF30 = drop_na(DF30)
DF31 <- filter(visz2500heat, visz2500heat$days == "32")
DF31 = merge(DF31, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF31 = drop_na(DF31)
DF32 <- filter(visz2500heat, visz2500heat$days == "33")
DF32 = merge(DF32, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF32 = drop_na(DF32)
DF33 <- filter(visz2500heat, visz2500heat$days == "34")
DF33 = merge(DF33, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF33 = drop_na(DF33)
DF34 <- filter(visz2500heat, visz2500heat$days == "35")
DF34 = merge(DF34, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF34 = drop_na(DF34)

#merge all dfs and then remove
heatDF <- bind_rows(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31,DF32,DF33,DF34)
rm(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31,DF32,DF33,DF34)

colnames(heatDF)[8] <- "Genotype"
length(unique(heatDF$Genotype))

if (heatDF[1,"Genotype"] == locations[1,"Genotype"]){print("true")}else {print("false")}

{print(locations[1,'Genotype'])}
str(locations)

#read in locations csv file and keep relevant info
locations <- read.csv("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-heatdrought/brachy-locations.csv")
locations = select(locations, Genotype, Collection_location, Latitude, Longitude, Elevation)

#merge locations file with data file
DF0 <- filter(heatDF, heatDF$days == "1")
DF0 = merge(DF0, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF0 = drop_na(DF0)
DF1 <- filter(heatDF, heatDF$days == "2")
DF1 = merge(DF1, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF1 = drop_na(DF1)
DF2 <- filter(heatDF, heatDF$days == "3")
DF2 = merge(DF2, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF2 = drop_na(DF2)
DF3 <- filter(heatDF, heatDF$days == "4")
DF3 = merge(DF3, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF3 = drop_na(DF3)
DF4 <- filter(heatDF, heatDF$days == "5")
DF4 = merge(DF4, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF4 = drop_na(DF4)
DF5 <- filter(heatDF, heatDF$days == "6")
DF5 = merge(DF5, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF5 = drop_na(DF5)
DF6 <- filter(heatDF, heatDF$days == "7")
DF6 = merge(DF6, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF6 = drop_na(DF6)
DF7 <- filter(heatDF, heatDF$days == "8")
DF7 = merge(DF7, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF7 = drop_na(DF7)
DF8 <- filter(heatDF, heatDF$days == "9")
DF8 = merge(DF8, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF8 = drop_na(DF8)
DF9 <- filter(heatDF, heatDF$days == "10")
DF9 = merge(DF9, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF9 = drop_na(DF9)
DF10 <- filter(heatDF, heatDF$days == "11")
DF10 = merge(DF10, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF10 = drop_na(DF10)
DF11 <- filter(heatDF, heatDF$days == "12")
DF11 = merge(DF11, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF11 = drop_na(DF11)
DF12 <- filter(heatDF, heatDF$days == "13")
DF12 = merge(DF12, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF12 = drop_na(DF12)
DF13 <- filter(heatDF, heatDF$days == "14")
DF13 = merge(DF13, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF13 = drop_na(DF13)
DF14 <- filter(heatDF, heatDF$days == "15")
DF14 = merge(DF14, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF14 = drop_na(DF14)
DF15 <- filter(heatDF, heatDF$days == "16")
DF15 = merge(DF15, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF15 = drop_na(DF15)
DF16 <- filter(heatDF, heatDF$days == "17")
DF16 = merge(DF16, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF16 = drop_na(DF16)
DF17 <- filter(heatDF, heatDF$days == "18")
DF17 = merge(DF17, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF17 = drop_na(DF17)
DF18 <- filter(heatDF, heatDF$days == "19")
DF18 = merge(DF18, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF18 = drop_na(DF18)
DF19 <- filter(heatDF, heatDF$days == "20")
DF19 = merge(DF19, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF19 = drop_na(DF19)
DF20 <- filter(heatDF, heatDF$days == "21")
DF20 = merge(DF20, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF20 = drop_na(DF20)
DF21 <- filter(heatDF, heatDF$days == "22")
DF21 = merge(DF21, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF21 = drop_na(DF21)
DF22 <- filter(heatDF, heatDF$days == "23")
DF22 = merge(DF22, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF22 = drop_na(DF22)
DF23 <- filter(heatDF, heatDF$days == "24")
DF23 = merge(DF23, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF23 = drop_na(DF23)
DF24 <- filter(heatDF, heatDF$days == "25")
DF24 = merge(DF24, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF24 = drop_na(DF24)
DF25 <- filter(heatDF, heatDF$days == "26")
DF25 = merge(DF25, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF25 = drop_na(DF25)
DF26 <- filter(heatDF, heatDF$days == "27")
DF26 = merge(DF26, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF26 = drop_na(DF26)
DF27 <- filter(heatDF, heatDF$days == "28")
DF27 = merge(DF27, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF27 = drop_na(DF27)
DF28 <- filter(heatDF, heatDF$days == "29")
DF28 = merge(DF28, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF28 = drop_na(DF28)
DF29 <- filter(heatDF, heatDF$days == "30")
DF29 = merge(DF29, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF29 = drop_na(DF29)
DF30 <- filter(heatDF, heatDF$days == "31")
DF30 = merge(DF30, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF30 = drop_na(DF30)
DF31 <- filter(heatDF, heatDF$days == "32")
DF31 = merge(DF31, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF31 = drop_na(DF31)
DF32 <- filter(heatDF, heatDF$days == "33")
DF32 = merge(DF32, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF32 = drop_na(DF32)
DF33 <- filter(heatDF, heatDF$days == "34")
DF33 = merge(DF33, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF33 = drop_na(DF33)
DF34 <- filter(heatDF, heatDF$days == "35")
DF34 = merge(DF34, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF34 = drop_na(DF34)

#merge all dfs and then remove
heatDF <- bind_rows(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31)
rm(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31)

length(unique(heatDF$Genotype.ID))

heatDF <- select(heatDF, image, plantbarcode, Genotype, date, time, days, Treatment.2, Replicate, value, label)

#save the csv file of the compiled VIS data
write.csv(heatDF, "./brachy_heat_VIS_SV_multi_data_compiled_currated.csv", row.names = FALSE)

#For heatDF keep only genotypes in both experiments
heatDF$compare.geno <- heatDF$Genotype %in% list1
heatDF <- filter(heatDF, compare.geno == TRUE)
length(unique(heatDF$Genotype))
heatDF = select(heatDF, -compare.geno)

#write csv with data from only genotypes that are in both experiments
write.csv(heatDF, "./brachy_heat_VIS_SV_multi_data_limited_genotypes.csv", row.names = FALSE)



#####################################################################
# Read in the heat/heat-drought single-value data #
#####################################################################
#read in single value trait csv files
heatz2500s <- read.csv("heatdrought_SV_z2500-single-value-traits.csv")
heatz3500s <- read.csv("heatdrought_SV_z3500-single-value-traits.csv")

#subset VIS data
heatz2500s = heatz2500s[heatz2500s$imgtype == 'VIS',]
heatz3500s = heatz3500s[heatz3500s$imgtype == 'VIS',]

#####################################################################
# Heat/heat-drought zoom correction #
#####################################################################
zoom.lm <- lm(zoom.camera ~ zoom, data = data.frame(zoom = c(1, 6000), zoom.camera = c(1, 6)))
# Download data for a reference object imaged at different zoom levels
if (!file.exists('zoom_calibration_data.txt')) {
  download.file('http://files.figshare.com/2084101/zoom_calibration_data.txt',
                'zoom_calibration_data.txt')
}
# Read zoom calibration data
z.data <- read.table(file = "zoom_calibration_data.txt", sep = "\t", header = TRUE)
# Calculate px per cm
z.data$px_cm <- z.data$length_px / z.data$length_cm
# Calculate area for each row
z.data$area_cm <- ifelse(z.data$reference == z.data$reference[[1]], (13.2*13.2), (13.2*3.7))
# Calculate px**2 per cm**2
z.data$px2_cm2 <- z.data$area_px / z.data$area_cm
# Convert LemnaTec zoom units to camera zoom units
z.data$zoom.camera <- predict(object = zoom.lm, newdata = z.data)
# Zoom correction for area
area.coef <- coef(nls(log(px2_cm2) ~ log(a * exp(b * zoom.camera)),
                      z.data, start = c(a = 1, b = 0.01)))
area.coef <- data.frame(a = area.coef[1], b = area.coef[2])
area.nls <- nls(px2_cm2 ~ a * exp(b * zoom.camera),
                data = z.data, start = c(a = area.coef$a, b = area.coef$b))
# Zoom correction for height
height.coef <- coef(nls(log(px_cm) ~ log(a * exp(b * zoom.camera)),
                        z.data, start = c(a = 1, b = 0.01)))
height.coef <- data.frame(a = area.coef[1], b = area.coef[2])
height.nls <- nls(px_cm ~ a * exp(b * zoom.camera),
                  data = z.data, start = c(a = area.coef$a, b = area.coef$b))

# change the zoom vector to the zooms in your dataset
exp_zooms = data.frame(zoom=c(500, 1500, 2500, 3500))
exp_zooms$zoom.camera = predict(object = zoom.lm, newdata = exp_zooms)
# this is the one I used:
exp_zooms$px2_cm2 = predict(object = area.nls, newdata = exp_zooms)
exp_zooms$px_cm = predict(object = height.nls, newdata = exp_zooms)

###################zoom correction
# at z2500
heatz2500s <- transform(heatz2500s, cor_area = area / 5543.662)
heatz2500s <- transform(heatz2500s, cor_area_above_reference = area_above_reference / 5543.662)
heatz2500s <- transform(heatz2500s, cor_area_below_reference = area_below_reference / 5543.662)
heatz2500s <- transform(heatz2500s, cor_width = width / 71.89759)
heatz2500s <- transform(heatz2500s, cor_height = height / 71.89759)
heatz2500s <- transform(heatz2500s, cor_height_above_reference = height_above_reference / 71.89759)
heatz2500s <- transform(heatz2500s, cor_height_below_reference = height_below_reference / 71.89759)
heatz2500s <- transform(heatz2500s, cor_convex_hull_area = convex_hull_area / 5543.662)
heatz2500s <- transform(heatz2500s, cor_perimeter = perimeter / 71.89759)
heatz2500s <- transform(heatz2500s, cor_longest_path = longest_path / 71.89759)
heatz2500s <- transform(heatz2500s, cor_ellipse_major_axis = ellipse_major_axis / 71.89759)
heatz2500s <- transform(heatz2500s, cor_ellipse_minor_axis = ellipse_minor_axis / 71.89759)

# at z3500
heatz3500s <- transform(heatz3500s, cor_area = area / 11772.169)
heatz3500s <- transform(heatz3500s, cor_area_above_reference = area_above_reference / 11772.169)
heatz3500s <- transform(heatz3500s, cor_area_below_reference = area_below_reference / 11772.169)
heatz3500s <- transform(heatz3500s, cor_width = width / 104.61683)
heatz3500s <- transform(heatz3500s, cor_height = height / 104.61683)
heatz3500s <- transform(heatz3500s, cor_height_above_reference = height_above_reference / 104.61683)
heatz3500s <- transform(heatz3500s, cor_height_below_reference = height_below_reference / 104.61683)
heatz3500s <- transform(heatz3500s, cor_convex_hull_area = convex_hull_area / 11772.169)
heatz3500s <- transform(heatz3500s, cor_perimeter = perimeter / 104.61683)
heatz3500s <- transform(heatz3500s, cor_longest_path = longest_path / 104.61683)
heatz3500s <- transform(heatz3500s, cor_ellipse_major_axis = ellipse_major_axis / 104.61683)
heatz3500s <- transform(heatz3500s, cor_ellipse_minor_axis = ellipse_minor_axis / 104.61683)


####################################
# Heat/heat-drought data filtering and organizing #
##################################
#filter data
heatz2500s <- select(heatz2500s, image, plantbarcode, timestamp, hue_circular_mean, hue_circular_std, hue_median, percent_unhealthy, cor_area, cor_width, cor_height, cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, cor_area_below_reference, percent_area_below_reference, cor_convex_hull_area, solidity, cor_perimeter, cor_longest_path, convex_hull_vertices, cor_ellipse_major_axis, cor_ellipse_minor_axis, ellipse_angle, ellipse_eccentricity)
heatz3500s <- select(heatz3500s, image, plantbarcode, timestamp, hue_circular_mean, hue_circular_std, hue_median, percent_unhealthy, cor_area, cor_width, cor_height, cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, cor_area_below_reference, percent_area_below_reference, cor_convex_hull_area, solidity, cor_perimeter, cor_longest_path, convex_hull_vertices, cor_ellipse_major_axis, cor_ellipse_minor_axis, ellipse_angle, ellipse_eccentricity)

#separate date and time
TS <- c("date","time") #make a vector for the column names of each part of the timestamp
heatz2500s <- separate(heatz2500s, col = timestamp, into = c("date","time"), sep = " ")
heatz3500s <- separate(heatz3500s, col = timestamp, into = c("date","time"), sep = " ")

#change date format to R standard for z2500
heatz2500s$date <- as.Date(heatz2500s$date, format="%m/%d/%y")

# split both z3500 dfs in half to separate out 2nd half that has incorrect date format
newDFs <- heatz3500s[1:22204,]
newDF2s <- heatz3500s[22205:44408,]

# Change date format to R standard for z3500
newDFs$date <- as.Date(newDFs$date)
newDF2s$date <- as.Date(newDF2s$date, format="%m/%d/%y")

#recombine z3500 dfs with correct date format 
newDFs$date <- as.character(newDFs$date)
newDF2s$date <- as.character(newDF2s$date)
heatz3500s <- bind_rows(newDFs, newDF2s)

heatz2500s$date <- as.character(heatz2500s$date)
heatz3500s$date <- as.character(heatz3500s$date)

#calculate days after first image
startdate1 <- as.POSIXct("2014-08-12")
heatz2500s$days = NA
heatz3500s$days = NA
heatz2500s$days = as.integer(difftime(heatz2500s$date, startdate1, units = "days"))
heatz3500s$days = as.integer(difftime(heatz3500s$date, startdate1, units = "days"))

#copy filtered dfs in case something bad happens
heatz2500s.1 <- heatz2500s
heatz3500s.1 <- heatz3500s

#merge barcodes with data file for VIS
DF0 <- filter(heatz3500s, heatz3500s$days == "1")
DF0 = merge(DF0, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF0 = drop_na(DF0)
DF1 <- filter(heatz3500s, heatz3500s$days == "2")
DF1 = merge(DF1, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF1 = drop_na(DF1)
DF2 <- filter(heatz3500s, heatz3500s$days == "3")
DF2 = merge(DF2, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF2 = drop_na(DF2)
DF3 <- filter(heatz3500s, heatz3500s$days == "4")
DF3 = merge(DF3, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF3 = drop_na(DF3)
DF4 <- filter(heatz3500s, heatz3500s$days == "5")
DF4 = merge(DF4, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF4 = drop_na(DF4)
DF5 <- filter(heatz3500s, heatz3500s$days == "6")
DF5 = merge(DF5, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF5 = drop_na(DF5)
DF6 <- filter(heatz3500s, heatz3500s$days == "7")
DF6 = merge(DF6, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF6 = drop_na(DF6)
DF7 <- filter(heatz3500s, heatz3500s$days == "8")
DF7 = merge(DF7, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF7 = drop_na(DF7)
DF8 <- filter(heatz3500s, heatz3500s$days == "9")
DF8 = merge(DF8, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF8 = drop_na(DF8)
DF9 <- filter(heatz3500s, heatz3500s$days == "10")
DF9 = merge(DF9, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF9 = drop_na(DF9)
DF10 <- filter(heatz2500s, heatz2500s$days == "11")
DF10 = merge(DF10, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF10 = drop_na(DF10)
DF11 <- filter(heatz2500s, heatz2500s$days == "12")
DF11 = merge(DF11, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF11 = drop_na(DF11)
DF12 <- filter(heatz2500s, heatz2500s$days == "13")
DF12 = merge(DF12, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF12 = drop_na(DF12)
DF13 <- filter(heatz2500s, heatz2500s$days == "14")
DF13 = merge(DF13, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF13 = drop_na(DF13)
DF14 <- filter(heatz2500s, heatz2500s$days == "15")
DF14 = merge(DF14, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF14 = drop_na(DF14)
DF15 <- filter(heatz2500s, heatz2500s$days == "16")
DF15 = merge(DF15, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF15 = drop_na(DF15)
DF16 <- filter(heatz2500s, heatz2500s$days == "17")
DF16 = merge(DF16, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF16 = drop_na(DF16)
DF17 <- filter(heatz2500s, heatz2500s$days == "18")
DF17 = merge(DF17, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF17 = drop_na(DF17)
DF18 <- filter(heatz2500s, heatz2500s$days == "19")
DF18 = merge(DF18, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF18 = drop_na(DF18)
DF19 <- filter(heatz2500s, heatz2500s$days == "20")
DF19 = merge(DF19, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF19 = drop_na(DF19)
DF20 <- filter(heatz2500s, heatz2500s$days == "21")
DF20 = merge(DF20, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF20 = drop_na(DF20)
DF21 <- filter(heatz2500s, heatz2500s$days == "22")
DF21 = merge(DF21, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF21 = drop_na(DF21)
DF22 <- filter(heatz2500s, heatz2500s$days == "23")
DF22 = merge(DF22, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF22 = drop_na(DF22)
DF23 <- filter(heatz2500s, heatz2500s$days == "24")
DF23 = merge(DF23, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF23 = drop_na(DF23)
DF24 <- filter(heatz2500s, heatz2500s$days == "25")
DF24 = merge(DF24, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF24 = drop_na(DF24)
DF25 <- filter(heatz2500s, heatz2500s$days == "26")
DF25 = merge(DF25, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF25 = drop_na(DF25)
DF26 <- filter(heatz2500s, heatz2500s$days == "27")
DF26 = merge(DF26, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF26 = drop_na(DF26)
DF27 <- filter(heatz2500s, heatz2500s$days == "28")
DF27 = merge(DF27, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF27 = drop_na(DF27)
DF28 <- filter(heatz2500s, heatz2500s$days == "29")
DF28 = merge(DF28, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF28 = drop_na(DF28)
DF29 <- filter(heatz2500s, heatz2500s$days == "30")
DF29 = merge(DF29, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF29 = drop_na(DF29)
DF30 <- filter(heatz2500s, heatz2500s$days == "31")
DF30 = merge(DF30, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF30 = drop_na(DF30)
DF31 <- filter(heatz2500s, heatz2500s$days == "32")
DF31 = merge(DF31, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF31 = drop_na(DF31)
DF32 <- filter(heatz2500s, heatz2500s$days == "33")
DF32 = merge(DF32, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF32 = drop_na(DF32)
DF33 <- filter(heatz2500s, heatz2500s$days == "34")
DF33 = merge(DF33, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF33 = drop_na(DF33)
DF34 <- filter(heatz2500s, heatz2500s$days == "35")
DF34 = merge(DF34, heatbarcodes, by.x = "plantbarcode", by.y = "Barcode", all = FALSE )
#DF34 = drop_na(DF34)

#merge all dfs and then remove
heatDFs <- bind_rows(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31,DF32,DF33,DF34)
rm(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31,DF32,DF33,DF34)

#colnames(heatDFs)[19] <- "Genotype"
length(unique(heatDFs$Genotype))

copy <- heatDFs

#read in locations csv file and keep relevant info
locations <- read.csv("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-heatdrought/brachy-locations.csv")
locations = select(locations, Genotype, Collection_location, Latitude, Longitude, Elevation)

print(locations$Genotype)

#merge locations file with data file
DF0 <- filter(heatDFs, heatDFs$days == "1")
DF0 = merge(DF0, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF0 = drop_na(DF0)
DF1 <- filter(heatDFs, heatDFs$days == "2")
DF1 = merge(DF1, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF1 = drop_na(DF1)
DF2 <- filter(heatDFs, heatDFs$days == "3")
DF2 = merge(DF2, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF2 = drop_na(DF2)
DF3 <- filter(heatDFs, heatDFs$days == "4")
DF3 = merge(DF3, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF3 = drop_na(DF3)
DF4 <- filter(heatDFs, heatDFs$days == "5")
DF4 = merge(DF4, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF4 = drop_na(DF4)
DF5 <- filter(heatDFs, heatDFs$days == "6")
DF5 = merge(DF5, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF5 = drop_na(DF5)
DF6 <- filter(heatDFs, heatDFs$days == "7")
DF6 = merge(DF6, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF6 = drop_na(DF6)
DF7 <- filter(heatDFs, heatDFs$days == "8")
DF7 = merge(DF7, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF7 = drop_na(DF7)
DF8 <- filter(heatDFs, heatDFs$days == "9")
DF8 = merge(DF8, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF8 = drop_na(DF8)
DF9 <- filter(heatDFs, heatDFs$days == "10")
DF9 = merge(DF9, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF9 = drop_na(DF9)
DF10 <- filter(heatDFs, heatDFs$days == "11")
DF10 = merge(DF10, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF10 = drop_na(DF10)
DF11 <- filter(heatDFs, heatDFs$days == "12")
DF11 = merge(DF11, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF11 = drop_na(DF11)
DF12 <- filter(heatDFs, heatDFs$days == "13")
DF12 = merge(DF12, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF12 = drop_na(DF12)
DF13 <- filter(heatDFs, heatDFs$days == "14")
DF13 = merge(DF13, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF13 = drop_na(DF13)
DF14 <- filter(heatDFs, heatDFs$days == "15")
DF14 = merge(DF14, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF14 = drop_na(DF14)
DF15 <- filter(heatDFs, heatDFs$days == "16")
DF15 = merge(DF15, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF15 = drop_na(DF15)
DF16 <- filter(heatDFs, heatDFs$days == "17")
DF16 = merge(DF16, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF16 = drop_na(DF16)
DF17 <- filter(heatDFs, heatDFs$days == "18")
DF17 = merge(DF17, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF17 = drop_na(DF17)
DF18 <- filter(heatDFs, heatDFs$days == "19")
DF18 = merge(DF18, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF18 = drop_na(DF18)
DF19 <- filter(heatDFs, heatDFs$days == "20")
DF19 = merge(DF19, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF19 = drop_na(DF19)
DF20 <- filter(heatDFs, heatDFs$days == "21")
DF20 = merge(DF20, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF20 = drop_na(DF20)
DF21 <- filter(heatDFs, heatDFs$days == "22")
DF21 = merge(DF21, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF21 = drop_na(DF21)
DF22 <- filter(heatDFs, heatDFs$days == "23")
DF22 = merge(DF22, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF22 = drop_na(DF22)
DF23 <- filter(heatDFs, heatDFs$days == "24")
DF23 = merge(DF23, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF23 = drop_na(DF23)
DF24 <- filter(heatDFs, heatDFs$days == "25")
DF24 = merge(DF24, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF24 = drop_na(DF24)
DF25 <- filter(heatDFs, heatDFs$days == "26")
DF25 = merge(DF25, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF25 = drop_na(DF25)
DF26 <- filter(heatDFs, heatDFs$days == "27")
DF26 = merge(DF26, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF26 = drop_na(DF26)
DF27 <- filter(heatDFs, heatDFs$days == "28")
DF27 = merge(DF27, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF27 = drop_na(DF27)
DF28 <- filter(heatDFs, heatDFs$days == "29")
DF28 = merge(DF28, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF28 = drop_na(DF28)
DF29 <- filter(heatDFs, heatDFs$days == "30")
DF29 = merge(DF29, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF29 = drop_na(DF29)
DF30 <- filter(heatDFs, heatDFs$days == "31")
DF30 = merge(DF30, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF30 = drop_na(DF30)
DF31 <- filter(heatDFs, heatDFs$days == "32")
DF31 = merge(DF31, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF31 = drop_na(DF31)
DF32 <- filter(heatDFs, heatDFs$days == "33")
DF32 = merge(DF32, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF32 = drop_na(DF32)
DF33 <- filter(heatDFs, heatDFs$days == "34")
DF33 = merge(DF33, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF33 = drop_na(DF33)
DF34 <- filter(heatDFs, heatDFs$days == "35")
DF34 = merge(DF34, locations, by.x = "Genotype", by.y = "Genotype", all = FALSE )
#DF34 = drop_na(DF34)

#merge all dfs and then remove
heatDFs <- bind_rows(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31,DF32,DF33,DF34)
rm(DF0,DF1,DF2,DF3,DF4,DF5,DF6,DF7,DF8,DF9,DF10,DF11,DF12,DF13,DF14,DF15,DF16,DF17,DF18,DF19,DF20,DF21,DF22,DF23,DF24,DF25,DF26,DF27,DF28,DF29,DF30,DF31,DF32,DF33,DF34)

length(unique(heatDFs$Genotype))

heatDFs <- select(heatDFs, image, plantbarcode, Genotype, date, time, days, Treatment.2, Replicate, Collection_location, Latitude, Longitude, Elevation, hue_circular_mean, hue_circular_std, hue_median, percent_unhealthy, cor_area, cor_width, cor_height, cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, cor_area_below_reference, percent_area_below_reference, cor_convex_hull_area, solidity, cor_perimeter, cor_longest_path, convex_hull_vertices, cor_ellipse_major_axis, cor_ellipse_minor_axis, ellipse_angle, ellipse_eccentricity)

#add image day (imgdays) column
heatDFs$imgdays = NA
heatDFs$imgdays[grep("1", heatDFs$days)] <- "1"
heatDFs$imgdays[grep("2", heatDFs$days)] <- "1"
heatDFs$imgdays[grep("3", heatDFs$days)] <- "2"
heatDFs$imgdays[grep("4", heatDFs$days)] <- "2"
heatDFs$imgdays[grep("5", heatDFs$days)] <- "3"
heatDFs$imgdays[grep("6", heatDFs$days)] <- "3"
heatDFs$imgdays[grep("7", heatDFs$days)] <- "4"
heatDFs$imgdays[grep("8", heatDFs$days)] <- "4"
heatDFs$imgdays[grep("9", heatDFs$days)] <- "5"
heatDFs$imgdays[grep("10", heatDFs$days)] <- "5"
heatDFs$imgdays[grep("11", heatDFs$days)] <- "6"
heatDFs$imgdays[grep("12", heatDFs$days)] <- "6"
heatDFs$imgdays[grep("13", heatDFs$days)] <- "7"
heatDFs$imgdays[grep("14", heatDFs$days)] <- "7"
heatDFs$imgdays[grep("15", heatDFs$days)] <- "8"
heatDFs$imgdays[grep("16", heatDFs$days)] <- "8"
heatDFs$imgdays[grep("17", heatDFs$days)] <- "9"
heatDFs$imgdays[grep("18", heatDFs$days)] <- "9"
heatDFs$imgdays[grep("19", heatDFs$days)] <- "10"
heatDFs$imgdays[grep("20", heatDFs$days)] <- "10"
heatDFs$imgdays[grep("21", heatDFs$days)] <- "11"
heatDFs$imgdays[grep("22", heatDFs$days)] <- "11"
heatDFs$imgdays[grep("23", heatDFs$days)] <- "12"
heatDFs$imgdays[grep("24", heatDFs$days)] <- "12"
heatDFs$imgdays[grep("25", heatDFs$days)] <- "13"
heatDFs$imgdays[grep("26", heatDFs$days)] <- "13"
heatDFs$imgdays[grep("27", heatDFs$days)] <- "14"
heatDFs$imgdays[grep("28", heatDFs$days)] <- "14"
heatDFs$imgdays[grep("29", heatDFs$days)] <- "15"
heatDFs$imgdays[grep("30", heatDFs$days)] <- "15"
heatDFs$imgdays[grep("31", heatDFs$days)] <- "16"
heatDFs$imgdays[grep("32", heatDFs$days)] <- "16"
heatDFs$imgdays[grep("33", heatDFs$days)] <- "17"
heatDFs$imgdays[grep("34", heatDFs$days)] <- "17"

heatDFs$Treatment = NA
heatDFs$Treatment[grep("20", heatDFs$Treatment.2)] <- "WL"
heatDFs$Treatment[grep("100", heatDFs$Treatment.2)] <- "WW"

# multiply the % unhealthy pixels by 100 to convert from fraction to percent
heatDFs[,16] <- heatDFs[,16]*100

#save the csv file of the compiled VIS data
write.csv(heatDFs, "./brachy_heat_VIS_SV_single_data_compiled_currated.csv", row.names = FALSE)

#keep only genotypes that are in both experiments and have complete data
heatDFs$compare.geno <- heatDFs$Genotype %in% list3
heatDFs <- filter(heatDFs, compare.geno == TRUE)
length(unique(heatDFs$Genotype))
heatDFs = select(heatDFs, -compare.geno)

#write csv with data from only genotypes that are in both experiments
write.csv(heatDFs, "./brachy_heat_VIS_SV_single_data_limited_genotypes.csv", row.names = FALSE)

###############################################################################################################################
### After first time start here ###
###############################################################################################################################



###########################################
# heatmap of percentage of unhealthy plant for control and drought #
###########################################
avg_unhealthy <- group_by(drought_127, Genotype, Treatment, imgdays, .drop=FALSE) %>% summarise(avg_percent_unhealthy = mean(percent_unhealthy))
avg_unhealthy <- ungroup(avg_unhealthy)

ctl_unhealthy <- filter(avg_unhealthy, Treatment == "WW") #make a DF with all images in the control treatment
drt_unhealthy <- filter (avg_unhealthy, Treatment == "WL") #make a DF with all images in the drought treamtment
ctl_unhealthy$imgdays <- as.numeric(ctl_unhealthy$imgdays)
drt_unhealthy$imgdays <- as.numeric(drt_unhealthy$imgdays)

#create a heatmap with a dendrogram
plot_unhealthy <- select(ctl_unhealthy, Genotype, imgdays, avg_percent_unhealthy) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plot_unhealthy$avg_percent_unhealthy <- as.numeric(plot_unhealthy$avg_percent_unhealthy) #make sure the ratio is numeric
spread_unhealthy <- spread(plot_unhealthy, key = imgdays, value = avg_percent_unhealthy) #this flips the dataframe so the days are across the top and the genotypes are down the first row
row.names_unhealthy <- spread_unhealthy[1:127,1] # save the row names
row.names_unhealthy <- t(row.names_unhealthy)
row.names_unhealthy <- as.character(row.names_unhealthy)
spread_unhealthy <- spread_unhealthy[,2:16] # save the DF without the first column
plot.mtx_unhealthy <- as.matrix(spread_unhealthy) #make the DF a matrix
rownames(plot.mtx_unhealthy) <- c(row.names_unhealthy) #add back row names
colnames(plot.mtx_unhealthy) <- as.numeric(c(1:15)) #make the column names the days
# plot the heatmap
#heatmap_height <- heatmap(plot.mtx_height,Colv=NA)
# plot the interactive heatmap
color.palette = colorRampPalette(c("blue","lightyellow","darkorange1"),space="rgb")
col_breaks = c(seq(0,12,length=100), # for blue
               seq(13,25,length=100),  # for yellow
               seq(26,100,length=100)) # for orange
pdf(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/unhealthy_heatmap_control.pdf",width = 10,height = 10)
png(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/unhealthy_heatmap_control.png",width = 900,height = 1000)
heatmapunhealthy <- heatmap.2(plot.mtx_unhealthy, 
                              Colv=FALSE,
                              Rowv = FALSE,
                              dendrogram = "none",
                              xlab = "Imaging Day",
                              ylab = "Genotype",
                              main = "Percent of Unhealthy Tissue in Plants in Control Conditions",
                              col = color.palette,
                              breaks = col_breaks,
                              #margin= c(5,5),
                              srtCol =0,
                              adjCol = c(0.5,.5),
                              #adjRow = c(0.1,0.5),
                              offsetRow = 0,
                              key.ylab = NA,
                              key.xlab = "Percent Stressed",
                              keysize=2,
                              lhei=c(2, 15),
                              lwid = c(2,6),
                              density.info = "none",
                              trace = "none") # plot the heatmap
dev.off()



#create a heatmap with a dendrogram
plot_unhealthy1 <- select(drt_unhealthy, Genotype, imgdays, avg_percent_unhealthy) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plot_unhealthy1$avg_percent_unhealthy <- as.numeric(plot_unhealthy1$avg_percent_unhealthy) #make sure the ratio is numeric
spread_unhealthy1 <- spread(plot_unhealthy1, key = imgdays, value = avg_percent_unhealthy) #this flips the dataframe so the days are across the top and the genotypes are down the first row
row.names_unhealthy1 <- spread_unhealthy1[1:127,1] # save the row names
row.names_unhealthy1 <- t(row.names_unhealthy1)
row.names_unhealthy1 <- as.character(row.names_unhealthy1)
spread_unhealthy1 <- spread_unhealthy1[,2:16] # save the DF without the first column
plot.mtx_unhealthy1 <- as.matrix(spread_unhealthy1) #make the DF a matrix
rownames(plot.mtx_unhealthy1) <- c(row.names_unhealthy1) #add back row names
colnames(plot.mtx_unhealthy1) <- as.numeric(c(1:15)) #make the column names the days
# plot the heatmap
#heatmap_height <- heatmap(plot.mtx_height,Colv=NA)
# plot the interactive heatmap
color.palette = colorRampPalette(c("blue","lightyellow","darkorange1"),space="rgb")
col_breaks = c(seq(0,12,length=100), # for blue
               seq(13,25,length=100),  # for yellow
               seq(26,100,length=100)) # for orange
pdf(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/unhealthy_heatmap_drought.pdf",width = 10,height = 10)
png(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/unhealthy_heatmap_drought.png",width = 900,height = 1000)
heatmapunhealthy1 <- heatmap.2(plot.mtx_unhealthy1, 
                               Colv=FALSE,
                               Rowv = FALSE,
                               dendrogram = "none",
                               xlab = "Imaging Day",
                               ylab = "Genotype",
                               main = "Percent of Unhealthy Tissue in Plants in Drought Conditions",
                               col = color.palette,
                               breaks = col_breaks,
                               #margin= c(5,5),
                               srtCol =0,
                               adjCol = c(0.5,.5),
                               #adjRow = c(0.1,0.5),
                               key.ylab = NA,
                               key.xlab = "Percent",
                               keysize=2,
                               lhei=c(2, 15),
                               lwid = c(2,6),
                               density.info = "none",
                               trace = "none") # plot the heatmap
dev.off()

spread_unhealthy1$compare.geno <- spread_unhealthy1$Genotype %in% list_turkey_110
spread_unhealthy1_turkey <- filter(spread_unhealthy1, compare.geno == TRUE)
length(unique(spread_unhealthy1_turkey$Genotype))
spread_unhealthy1_turkey = select(spread_unhealthy1_turkey, -compare.geno)
spread_unhealthy1 = select(spread_unhealthy1, -compare.geno)

spread_unhealthy1_turkey <- select(spread_unhealthy1_turkey, Genotype, 16)
top_spread_unhealthy1_turkey <- spread_unhealthy1_turkey %>% top_n(27)
bottom_spread_unhealthy1_turkey <- spread_unhealthy1_turkey %>% top_n(-27)

test <- top_spread_unhealthy1_turkey
colnames(test)[2] <- "percent"


# Create a function to print squares of numbers in sequence.
new.function <- function(a) {
  for(i in 1:a) {
    b <- i^2
    print(b)
  }
}

# Call the function new.function supplying 6 as an argument.
new.function(6)




test.function <- function(dataframe, cutoff){
  percents <- dataframe$percent
  x <- sort(percents, decreasing = TRUE)
  y <- x[x >= cutoff] 
  print(y)
}

test.function(test, 20)




###########################################
# biomass heatmaps of control and drought conditions #
###########################################

#calculate average area and create avg df
avgDFs <- group_by(drought_127, Genotype, Treatment, imgdays) %>% summarise(avg_area = mean(cor_area))
avgDFs <- ungroup(avgDFs)
#remove imgday 2 (data super messy)
avgDFs <- avgDFs[!(avgDFs$imgdays=="2"),]

#heatmap of biomass in control
area <- filter(avgDFs, Treatment == "WW") #make a DF with all images in the control treatment
area$imgdays <- as.numeric(area$imgdays)

#create a heatmap with a dendrogram
plotDF <- select(area, Genotype, imgdays, avg_area) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDF$avg_area <- as.numeric(plotDF$avg_area) #make sure the area is numeric
spreadDF <- spread(plotDF, key = imgdays, value = avg_area) #this flips the dataframe so the days are across the top and the genotypes are down the first row
#spreadDFheat <- spreadDFheat[-c(24, 25, 31, 63, 69, 101), ]
row.names <- spreadDF[1:127,1] # save the row names
row.names <- t(row.names)
row.names <- as.character(row.names)
spreadDF <- spreadDF[,2:15] # save the DF without the first column
plot.mtx <- as.matrix(spreadDF) #make the DF a matrix
rownames(plot.mtx) <- c(row.names) #add back row names
colnames(plot.mtx) <- as.numeric(c(1,3:15)) #make the column names the days
# plot the heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,2,length=100), # for orange
               seq(3,4,length=100),  # for yellow
               seq(5,40,length=100)) # for blue
pdf(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/area_heatmap_control.pdf",width = 10,height = 10)
png(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/area_heatmap_control.png",width = 900,height = 1000)
heatmapctl <- heatmap.2(plot.mtx, 
                        Colv=FALSE,
                        Rowv = FALSE,
                        dendrogram = "none",
                        xlab = "Imaging Day",
                        ylab = "Genotype",
                        main = "Heatmap of Average Plant Area in Control",
                        col = color.palette,
                        breaks = col_breaks,
                        #margin= c(5,5),
                        srtCol =0,
                        adjCol = c(0.5,.5),
                        #adjRow = c(0.1,0.5),
                        key.ylab = NA,
                        key.xlab = "Area Value",
                        keysize=2,
                        lhei=c(2, 15),
                        lwid = c(2,6),
                        density.info = "none",
                        trace = "none") # plot the heatmap
dev.off()


#biomass heatmap drought
areadrt <- filter(avgDFs, Treatment == "WL") #make a DF with all images in the control treatment
areadrt$imgdays <- as.numeric(areadrt$imgdays)

#create a heatmap with a dendrogram
plotDF1 <- select(areadrt, Genotype, imgdays, avg_area) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDF1$avg_area <- as.numeric(plotDF1$avg_area) #make sure the area is numeric
spreadDF1 <- spread(plotDF1, key = imgdays, value = avg_area) #this flips the dataframe so the days are across the top and the genotypes are down the first row
#spreadDFheat <- spreadDFheat[-c(24, 25, 31, 63, 69, 101), ]
row.names1 <- spreadDF1[1:127,1] # save the row names
row.names1 <- t(row.names1)
row.names1 <- as.character(row.names1)
spreadDF1 <- spreadDF1[,2:15] # save the DF without the first column
plot.mtx1 <- as.matrix(spreadDF1) #make the DF a matrix
rownames(plot.mtx1) <- c(row.names1) #add back row names
colnames(plot.mtx1) <- as.numeric(c(1,3:15)) #make the column names the days
# plot the heatmap
#heatmap <- heatmap(plot.mtx,Colv=NA)
# plot the interactive heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,2,length=100), # for orange
               seq(3,4,length=100),  # for yellow
               seq(5,40,length=100)) # for blue
pdf(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/area_heatmap_drought.pdf",width = 8,height = 8)
png(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/area_heatmap_drought.png",width = 900,height = 1000)
heatmapdrt <- heatmap.2(plot.mtx1, 
                        Colv=FALSE,
                        Rowv = FALSE,
                        dendrogram = "none",
                        xlab = "Imaging Day",
                        ylab = "Genotype",
                        main = "Heatmap of Average Plant Area in Drought",
                        col = color.palette,
                        breaks = col_breaks,
                        #margin= c(5,5),
                        srtCol =0,
                        adjCol = c(0.5,.5),
                        #adjRow = c(0.1,0.5),
                        key.ylab = NA,
                        key.xlab = "Area Value",
                        keysize=2,
                        lhei=c(2, 15),
                        lwid = c(2,6),
                        density.info = "none",
                        trace = "none") # plot the heatmap
dev.off()

spreadDF1$compare.geno <- spreadDF1$Genotype %in% list_turkey_110
spreadDF1_turkey <- filter(spreadDF1, compare.geno == TRUE)
length(unique(spreadDF1_turkey$Genotype))
spreadDF1_turkey = select(spreadDF1_turkey, -compare.geno)
spreadDF1 = select(spreadDF1, -compare.geno)

spreadDF1_turkey <- select(spreadDF1_turkey, Genotype, 15)
top_spreadDF1_turkey <- spreadDF1_turkey %>% top_n(27)
bottom_spreadDF1_turkey <- spreadDF1_turkey %>% top_n(-27)


###########################################
# height heatmaps of control and drought conditions #
###########################################

#calculate average area and create avg df
avgDFs_height <- group_by(drought_127, Genotype, Treatment, imgdays) %>% summarise(avg_height = mean(cor_height_above_reference))
avgDFs_height <- ungroup(avgDFs_height)

#remove imgday 2 (data super messy)
avgDFs_height <- avgDFs_height[!(avgDFs_height$imgdays=="2"),]

#heatmap of biomass in control
height <- filter(avgDFs_height, Treatment == "WW") #make a DF with all images in the control treatment
height$imgdays <- as.numeric(height$imgdays)

#create a heatmap with a dendrogram
plotDFheight <- select(height, Genotype, imgdays, avg_height) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDFheight$avg_height <- as.numeric(plotDFheight$avg_height) #make sure the area is numeric
spreadDFheight <- spread(plotDFheight, key = imgdays, value = avg_height) #this flips the dataframe so the days are across the top and the genotypes are down the first row
#spreadDFheat <- spreadDFheat[-c(24, 25, 31, 63, 69, 101), ]
row.namesheight <- spreadDFheight[1:127,1] # save the row names
row.namesheight <- t(row.namesheight)
row.namesheight <- as.character(row.namesheight)
spreadDFheight <- spreadDFheight[,2:15] # save the DF without the first column
plot.mtxheight <- as.matrix(spreadDFheight) #make the DF a matrix
rownames(plot.mtxheight) <- c(row.namesheight) #add back row names
colnames(plot.mtxheight) <- as.numeric(c(1,3:15)) #make the column names the days
# plot the heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,4,length=20), # for orange
               seq(5,7,length=20),  # for yellow
               seq(8,20,length=20)) # for blue
pdf(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/height_heatmap_control.pdf",width = 10,height = 10)
png(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/height_heatmap_control.png",width = 900,height = 1000)
heatmapctlheight <- heatmap.2(plot.mtxheight, 
                              Colv=FALSE,
                              Rowv = FALSE,
                              dendrogram = "none",
                              xlab = "Imaging Day",
                              ylab = "Genotype",
                              main = "Average Plant Height in Control",
                              col = color.palette,
                              breaks = col_breaks,
                              #margin= c(5,5),
                              srtCol =0,
                              adjCol = c(0.5,.5),
                              #adjRow = c(0.1,0.5),
                              key.ylab = NA,
                              key.xlab = "Area Value",
                              keysize=2,
                              lhei=c(2, 15),
                              lwid = c(2,6),
                              density.info = "none",
                              trace = "none") # plot the heatmap
dev.off()


#height heatmap drought
heightdrt <- filter(avgDFs_height, Treatment == "WL") #make a DF with all images in the control treatment
heightdrt$imgdays <- as.numeric(heightdrt$imgdays)

#create a heatmap with a dendrogram
plotDFheight1 <- select(heightdrt, Genotype, imgdays, avg_height) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDFheight1$avg_height <- as.numeric(plotDFheight1$avg_height) #make sure the area is numeric
spreadDFheight1 <- spread(plotDFheight1, key = imgdays, value = avg_height) #this flips the dataframe so the days are across the top and the genotypes are down the first row
#spreadDFheat <- spreadDFheat[-c(24, 25, 31, 63, 69, 101), ]
row.namesheight1 <- spreadDFheight1[1:127,1] # save the row names
row.namesheight1 <- t(row.namesheight1)
row.namesheight1 <- as.character(row.namesheight1)
spreadDFheight1 <- spreadDFheight1[,2:15] # save the DF without the first column
plot.mtxheight1 <- as.matrix(spreadDFheight1) #make the DF a matrix
rownames(plot.mtxheight1) <- c(row.namesheight1) #add back row names
colnames(plot.mtxheight1) <- as.numeric(c(1,3:15)) #make the column names the days
# plot the heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,4,length=20), # for orange
               seq(5,7,length=20),  # for yellow
               seq(8,20,length=20)) # for blue
pdf(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/height_heatmap_drought.pdf",width = 10,height = 10)
png(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/height_heatmap_drought.png",width = 900,height = 1000)
heatmapdrtheight <- heatmap.2(plot.mtxheight1, 
                              Colv=FALSE,
                              Rowv = FALSE,
                              dendrogram = "none",
                              xlab = "Imaging Day",
                              ylab = "Genotype",
                              main = "Average Plant Height in Drought",
                              col = color.palette,
                              breaks = col_breaks,
                              #margin= c(5,5),
                              srtCol =0,
                              adjCol = c(0.5,.5),
                              #adjRow = c(0.1,0.5),
                              key.ylab = NA,
                              key.xlab = "Area Value",
                              keysize=2,
                              lhei=c(2, 15),
                              lwid = c(2,6),
                              density.info = "none",
                              trace = "none") # plot the heatmap
dev.off()




spreadDFheight1$compare.geno <- spreadDFheight1$Genotype %in% list_turkey_110
spreadDFheight1_turkey <- filter(spreadDFheight1, compare.geno == TRUE)
length(unique(spreadDFheight1_turkey$Genotype))
spreadDFheight1_turkey = select(spreadDFheight1_turkey, -compare.geno)
spreadDFheight1 = select(spreadDFheight1, -compare.geno)

spreadDFheight1_turkey <- select(spreadDFheight1_turkey, Genotype, 15)
top_spreadDFheight1_turkey <- spreadDFheight1_turkey %>% top_n(27)
bottom_spreadDFheight1_turkey <- spreadDFheight1_turkey %>% top_n(-27)





###########################################
# heatmap of percentage of unhealthy plant for heat and heat-drought #
###########################################
#heatmap of heat unhealthy
heat_127$percent_unhealthy <- as.numeric(heat_127$percent_unhealthy)
avg_unhealthyheat <- group_by(heat_127, Genotype, Treatment, imgdays) %>% summarise(avg_percent_unhealthy = mean(percent_unhealthy, na.rm=TRUE))
avg_unhealthyheat <- ungroup(avg_unhealthyheat)

ctl_unhealthyheat <- filter(avg_unhealthyheat, Treatment == "WW") #make a DF with all images in the control treatment
drt_unhealthyheat <- filter (avg_unhealthyheat, Treatment == "WL") #make a DF with all images in the drought treamtment
ctl_unhealthyheat$imgdays <- as.numeric(ctl_unhealthyheat$imgdays)
drt_unhealthyheat$imgdays <- as.numeric(drt_unhealthyheat$imgdays)

#remove imgday 3 (data super messy)
ctl_unhealthyheat <- ctl_unhealthyheat[!(ctl_unhealthyheat$imgdays=="3"),]

#create a heatmap with a dendrogram
plot_unhealthyheat <- select(ctl_unhealthyheat, Genotype, imgdays, avg_percent_unhealthy) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plot_unhealthyheat$avg_percent_unhealthy <- as.numeric(plot_unhealthyheat$avg_percent_unhealthy) #make sure the ratio is numeric
spread_unhealthyheat <- spread(plot_unhealthyheat, key = imgdays, value = avg_percent_unhealthy) #this flips the dataframe so the days are across the top and the genotypes are down the first row
row.names_unhealthyheat <- spread_unhealthyheat[1:127,1] # save the row names
row.names_unhealthyheat <- t(row.names_unhealthyheat)
row.names_unhealthyheat <- as.character(row.names_unhealthyheat)
spread_unhealthyheat <- spread_unhealthyheat[,2:17] # save the DF without the first column
plot.mtx_unhealthyheat <- as.matrix(spread_unhealthyheat) #make the DF a matrix
rownames(plot.mtx_unhealthyheat) <- c(row.names_unhealthyheat) #add back row names
colnames(plot.mtx_unhealthyheat) <- as.numeric(c(2:17)) #make the column names the days
color.palette = colorRampPalette(c("blue","lightyellow","darkorange1"),space="rgb")
col_breaks = c(seq(0,12,length=100), # for blue
               seq(13,25,length=100),  # for yellow
               seq(26,100,length=100)) # for orange
pdf(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/unhealthy_heatmap_heat.pdf",width = 10,height = 10)
png(file="/Users/ellaludwig/Desktop/test_unhealthy_heatmap_heat.png",width = 900,height = 1000)
heatmapunhealthyheat <- heatmap.2(plot.mtx_unhealthyheat, 
                                  Colv=FALSE,
                                  Rowv = FALSE,
                                  dendrogram = "none",
                                  xlab = "Imaging Day",
                                  ylab = "Genotype",
                                  main = "Percent of Unhealthy Tissue in Plants in Heat",
                                  col = color.palette,
                                  breaks = col_breaks,
                                  #margin= c(5,5),
                                  srtCol =0,
                                  adjCol = c(0.5,.5),
                                  #adjRow = c(0.1,0.5),
                                  offsetRow = 0,
                                  key.ylab = NA,
                                  key.xlab = "Percent Stressed",
                                  keysize=2,
                                  lhei=c(2, 15),
                                  lwid = c(2,6),
                                  density.info = "none",
                                  trace = "none") # plot the heatmap
dev.off()
heatmapunhealthyheat

spread_unhealthyheat$compare.geno <- spread_unhealthyheat$Genotype %in% list_turkey_110
spread_unhealthyheat_turkey <- filter(spread_unhealthyheat, compare.geno == TRUE)
length(unique(spread_unhealthyheat_turkey$Genotype))
spread_unhealthyheat_turkey = select(spread_unhealthyheat_turkey, -compare.geno)
spread_unhealthyheat = select(spread_unhealthyheat, -compare.geno)

spread_unhealthyheat_turkey <- select(spread_unhealthyheat_turkey, Genotype, 17)
top_spread_unhealthyheat_turkey <- spread_unhealthyheat_turkey %>% top_n(27)
bottom_spread_unhealthyheat_turkey <- spread_unhealthyheat_turkey %>% top_n(-27)

topheatmapint_heat <- heatmaply(plot.mtx_unhealthyheat,
                        row_dend_left = FALSE,
                        Colv = FALSE,
                        xlab = "Imaging Day",
                        ylab = "Genotype",
                        main = "Percent of Unhealthy Tissue in Plants in Heat",
                        dendrogram='none',
                        scale='none',
                        col=color.palette,
                        breaks = col_breaks,
                        trace='none',
                        symbreaks=TRUE,
                        cex.main=0.75,
                        keysize=1.4)
heatmapint_heat



#create a heatmap for heat & drought unhealthy
plot_unhealthyheat1 <- select(drt_unhealthyheat, Genotype, imgdays, avg_percent_unhealthy) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plot_unhealthyheat1$avg_percent_unhealthy <- as.numeric(plot_unhealthyheat1$avg_percent_unhealthy) #make sure the ratio is numeric
spread_unhealthyheat1 <- spread(plot_unhealthyheat1, key = imgdays, value = avg_percent_unhealthy) #this flips the dataframe so the days are across the top and the genotypes are down the first row
row.names_unhealthyheat1 <- spread_unhealthyheat1[1:127,1] # save the row names
row.names_unhealthyheat1 <- t(row.names_unhealthyheat1)
row.names_unhealthyheat1 <- as.character(row.names_unhealthyheat1)
spread_unhealthyheat1 <- spread_unhealthyheat1[,2:18] # save the DF without the first column
plot.mtx_unhealthyheat1 <- as.matrix(spread_unhealthyheat1) #make the DF a matrix
rownames(plot.mtx_unhealthyheat1) <- c(row.names_unhealthyheat1) #add back row names
colnames(plot.mtx_unhealthyheat1) <- as.numeric(c(1:17)) #make the column names the days
color.palette = colorRampPalette(c("blue","lightyellow","darkorange1"),space="rgb")
col_breaks = c(seq(0,12,length=100), # for blue
               seq(13,25,length=100),  # for yellow
               seq(26,100,length=100)) # for orange
pdf(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/unhealthy_heatmap_heat&drought.pdf",width = 10,height = 10)
png(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/unhealthy_heatmap_heat&drought.png",width = 900,height = 1000)
heatmapunhealthyheat1 <- heatmap.2(plot.mtx_unhealthyheat1, 
                                   Colv=FALSE,
                                   Rowv = FALSE,
                                   dendrogram = "none",
                                   xlab = "Imaging Day",
                                   ylab = "Genotype",
                                   main = "Percent of Unhealthy Tissue in Plants in Heat & Drought",
                                   col = color.palette,
                                   breaks = col_breaks,
                                   #margin= c(5,5),
                                   srtCol =0,
                                   adjCol = c(0.5,.5),
                                   #adjRow = c(0.1,0.5),
                                   offsetRow = 0,
                                   key.ylab = NA,
                                   key.xlab = "Percent Stressed",
                                   keysize=2,
                                   lhei=c(2, 15),
                                   lwid = c(2,6),
                                   density.info = "none",
                                   trace = "none") # plot the heatmap
dev.off()

spread_unhealthyheat1$compare.geno <- spread_unhealthyheat1$Genotype %in% list_turkey_110
spread_unhealthyheat1_turkey <- filter(spread_unhealthyheat1, compare.geno == TRUE)
length(unique(spread_unhealthyheat1_turkey$Genotype))
spread_unhealthyheat1_turkey = select(spread_unhealthyheat1_turkey, -compare.geno)
spread_unhealthyheat1 = select(spread_unhealthyheat1, -compare.geno)

spread_unhealthyheat1_turkey <- select(spread_unhealthyheat1_turkey, Genotype, 18)
top_spread_unhealthyheat1_turkey <- spread_unhealthyheat1_turkey %>% top_n(27)
bottom_spread_unhealthyheat1_turkey <- spread_unhealthyheat1_turkey %>% top_n(-27)

# drought: top_spread_unhealthy1_turkey, bottom_spread_unhealthy1_turkey
# heat: top_spread_unhealthyheat_turkey, bottom_spread_unhealthyheat_turkey
# heat-drought: top_spread_unhealthyheat1_turkey, bottom_spread_unhealthyheat1_turkey

###see which accessions overlap that perform poorly across stress treatments
# compare heat to heat-drought
list5 <- unique(top_spread_unhealthyheat1_turkey$Genotype)
top_spread_unhealthyheat_turkey$compare.geno <- top_spread_unhealthyheat_turkey$Genotype %in% list5
overlap_stressed_top <- filter(top_spread_unhealthyheat_turkey, compare.geno == TRUE)
length(unique(overlap_stressed_top$Genotype))
overlap_stressed_top = select(overlap_stressed_top, -compare.geno)
list6 <- unique(overlap_stressed_top$Genotype)
#compare heat and drought
list5 <- unique(top_spread_unhealthyheat_turkey$Genotype)
top_spread_unhealthy1_turkey$compare.geno <- top_spread_unhealthy1_turkey$Genotype %in% list5
overlap_stressed_top <- filter(top_spread_unhealthy1_turkey, compare.geno == TRUE)
length(unique(overlap_stressed_top$Genotype))
overlap_stressed_top = select(overlap_stressed_top, -compare.geno)
list9 <- unique(overlap_stressed_top$Genotype)
#compare heat-drought and drought
list5 <- unique(top_spread_unhealthyheat1_turkey$Genotype)
top_spread_unhealthy1_turkey$compare.geno <- top_spread_unhealthy1_turkey$Genotype %in% list5
overlap_stressed_top <- filter(top_spread_unhealthy1_turkey, compare.geno == TRUE)
length(unique(overlap_stressed_top$Genotype))
overlap_stressed_top = select(overlap_stressed_top, -compare.geno)
list10 <- unique(overlap_stressed_top$Genotype)
#compare heat/heat-drought to drought to find overall overlap
top_spread_unhealthy1_turkey$compare.geno <- top_spread_unhealthy1_turkey$Genotype %in% list6
overlap_stressed_top <- filter(top_spread_unhealthy1_turkey, compare.geno == TRUE)



### only turkish accessions!! see which accessions overlap that perform well
# compare heat to heat-drought
list7 <- unique(bottom_spread_unhealthyheat1_turkey$Genotype)
bottom_spread_unhealthyheat_turkey$compare.geno <- bottom_spread_unhealthyheat_turkey$Genotype %in% list7
overlap_stressed_bottom <- filter(bottom_spread_unhealthyheat_turkey, compare.geno == TRUE)
length(unique(overlap_stressed_bottom$Genotype))
overlap_stressed_bottom = select(overlap_stressed_bottom, -compare.geno)
list8 <- unique(overlap_stressed_bottom$Genotype)
#compare heat and drought
list7 <- unique(bottom_spread_unhealthyheat_turkey$Genotype)
bottom_spread_unhealthy1_turkey$compare.geno <- bottom_spread_unhealthy1_turkey$Genotype %in% list7
overlap_stressed_bottom <- filter(bottom_spread_unhealthy1_turkey, compare.geno == TRUE)
length(unique(overlap_stressed_bottom$Genotype))
list11 <- unique(overlap_stressed_bottom$Genotype)
#compare heat-drought and drought
list7 <- unique(bottom_spread_unhealthyheat1_turkey$Genotype)
bottom_spread_unhealthy1_turkey$compare.geno <- bottom_spread_unhealthy1_turkey$Genotype %in% list7
overlap_stressed_bottom <- filter(bottom_spread_unhealthy1_turkey, compare.geno == TRUE)
length(unique(overlap_stressed_bottom$Genotype))
list12 <- unique(overlap_stressed_bottom$Genotype)
#compare heat/heat-drought to drought to find overall overlap
bottom_spread_unhealthy1_turkey$compare.geno <- bottom_spread_unhealthy1_turkey$Genotype %in% list8
overlap_stressed_bottom <- filter(bottom_spread_unhealthy1_turkey, compare.geno == TRUE)



spread_unhealthyheat1 <- select(spread_unhealthyheat1, Genotype, 18)
top_spread_unhealthyheat1 <- spread_unhealthyheat1 %>% top_n(31)
bottom_spread_unhealthyheat1 <- spread_unhealthyheat1 %>% top_n(-31)

spread_unhealthyheat <- select(spread_unhealthyheat, Genotype, 17)
top_spread_unhealthyheat <- spread_unhealthyheat %>% top_n(31)
bottom_spread_unhealthyheat <- spread_unhealthyheat %>% top_n(-31)

spread_unhealthy1 <- select(spread_unhealthy1, Genotype, 16)
top_spread_unhealthy1 <- spread_unhealthy1 %>% top_n(31)
bottom_spread_unhealthy1 <- spread_unhealthy1 %>% top_n(-31)


##### all accessions!! 
###see which accessions overlap that perform poorly across stress treatments
# compare heat to heat-drought
list35 <- unique(top_spread_unhealthyheat1$Genotype)
top_spread_unhealthyheat$compare.geno <- top_spread_unhealthyheat$Genotype %in% list35
all_overlap_stressed_top <- filter(top_spread_unhealthyheat, compare.geno == TRUE)
length(unique(all_overlap_stressed_top$Genotype))
all_overlap_stressed_top = select(all_overlap_stressed_top, -compare.geno)
list36 <- unique(all_overlap_stressed_top$Genotype)
#compare heat and drought
list35 <- unique(top_spread_unhealthyheat$Genotype)
top_spread_unhealthy1$compare.geno <- top_spread_unhealthy1$Genotype %in% list35
all_overlap_stressed_top <- filter(top_spread_unhealthy1, compare.geno == TRUE)
length(unique(all_overlap_stressed_top$Genotype))
all_overlap_stressed_top = select(all_overlap_stressed_top, -compare.geno)
list39 <- unique(all_overlap_stressed_top$Genotype)
#compare heat-drought and drought
list35 <- unique(top_spread_unhealthyheat1$Genotype)
top_spread_unhealthy1$compare.geno <- top_spread_unhealthy1$Genotype %in% list35
all_overlap_stressed_top <- filter(top_spread_unhealthy1, compare.geno == TRUE)
length(unique(all_overlap_stressed_top$Genotype))
all_overlap_stressed_top = select(all_overlap_stressed_top, -compare.geno)
list40 <- unique(all_overlap_stressed_top$Genotype)
#compare heat/heat-drought to drought to find overall overlap
top_spread_unhealthy1$compare.geno <- top_spread_unhealthy1$Genotype %in% list36
all_overlap_stressed_top <- filter(top_spread_unhealthy1, compare.geno == TRUE)



###see which accessions overlap that perform well
# compare heat to heat-drought
list37 <- unique(bottom_spread_unhealthyheat1$Genotype)
bottom_spread_unhealthyheat$compare.geno <- bottom_spread_unhealthyheat$Genotype %in% list37
all_overlap_stressed_bottom <- filter(bottom_spread_unhealthyheat, compare.geno == TRUE)
length(unique(all_overlap_stressed_bottom$Genotype))
all_overlap_stressed_bottom = select(all_overlap_stressed_bottom, -compare.geno)
list38 <- unique(all_overlap_stressed_bottom$Genotype)
#compare heat and drought
list37 <- unique(bottom_spread_unhealthyheat$Genotype)
bottom_spread_unhealthy1$compare.geno <- bottom_spread_unhealthy1$Genotype %in% list37
all_overlap_stressed_bottom <- filter(bottom_spread_unhealthy1, compare.geno == TRUE)
length(unique(all_overlap_stressed_bottom$Genotype))
list41 <- unique(all_overlap_stressed_bottom$Genotype)
#compare heat-drought and drought
list37 <- unique(bottom_spread_unhealthyheat1$Genotype)
bottom_spread_unhealthy1$compare.geno <- bottom_spread_unhealthy1$Genotype %in% list37
all_overlap_stressed_bottom <- filter(bottom_spread_unhealthy1, compare.geno == TRUE)
length(unique(all_overlap_stressed_bottom$Genotype))
list42 <- unique(all_overlap_stressed_bottom$Genotype)
#compare heat/heat-drought to drought to find overall overlap
bottom_spread_unhealthy1$compare.geno <- bottom_spread_unhealthy1$Genotype %in% list38
all_overlap_stressed_bottom <- filter(bottom_spread_unhealthy1, compare.geno == TRUE)



###########################################
# biomass heatmaps of heat and heat-drought conditions #
###########################################
#biomass heatmap for heat
#calculate average area and create avg df
avgheatDFs <- group_by(heat_127, Genotype, Treatment, imgdays) %>% summarise(avg_area = mean(cor_area))
avgheatDFs <- ungroup(avgheatDFs)

heatarea <- filter(avgheatDFs, Treatment == "WW") #make a DF with all images in the control treatment
heatarea$imgdays <- as.numeric(heatarea$imgdays)

#remove imgday 3 (data super messy)
heatarea <- heatarea[!(heatarea$imgdays=="3"),]

#create a heatmap with a dendrogram
plotDFheat <- select(heatarea, Genotype, imgdays, avg_area) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDFheat$avg_area <- as.numeric(plotDFheat$avg_area) #make sure the area is numeric
spreadDFheat <- spread(plotDFheat, key = imgdays, value = avg_area) #this flips the dataframe so the days are across the top and the genotypes are down the first row
#spreadDFheat <- spreadDFheat[-c(24, 25, 31, 63, 69, 101), ]
row.namesheat <- spreadDFheat[1:127,1] # save the row names
row.namesheat <- t(row.namesheat)
row.namesheat <- as.character(row.namesheat)
spreadDFheat <- spreadDFheat[,2:17] # save the DF without the first column
plot.mtxheat <- as.matrix(spreadDFheat) #make the DF a matrix
rownames(plot.mtxheat) <- c(row.namesheat) #add back row names
colnames(plot.mtxheat) <- as.numeric(c(2:17)) #make the column names the days
# plot the heatmap
#heatmap <- heatmap(plot.mtx,Colv=NA)
# plot the interactive heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,2,length=100), # for orange
               seq(3,4,length=100),  # for yellow
               seq(5,40,length=100)) # for blue
pdf(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/area_heatmap_heat.pdf",width = 10,height = 10)
png(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/area_heatmap_heat.png",width = 900,height = 1000)
png(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/heatmap_genotypes.png",width = 900,height = 1200)

heatmapheat <- heatmap.2(plot.mtxheat, 
                         Colv=FALSE,
                         Rowv = FALSE,
                         dendrogram = "none",
                         xlab = "Imaging Day",
                         ylab = "Genotype",
                         main = "Heatmap of Average Plant Area in Heat",
                         col = color.palette,
                         breaks = col_breaks,
                         #margin= c(5,5),
                         srtCol =0,
                         adjCol = c(0.5,.5),
                         #adjRow = c(0.1,0.5),
                         key.ylab = NA,
                         key.xlab = "Area Value",
                         keysize=2,
                         lhei=c(2, 15),
                         lwid = c(2,6),
                         density.info = "none",
                         trace = "none") # plot the heatmap
dev.off()


spreadDFheat$compare.geno <- spreadDFheat$Genotype %in% list_turkey_110
spreadDFheat_turkey <- filter(spreadDFheat, compare.geno == TRUE)
length(unique(spreadDFheat_turkey$Genotype))
spreadDFheat_turkey = select(spreadDFheat_turkey, -compare.geno)
spreadDFheat = select(spreadDFheat, -compare.geno)

spreadDFheat_turkey <- select(spreadDFheat_turkey, Genotype, 17)
top_spreadDFheat_turkey <- spreadDFheat_turkey %>% top_n(27)
bottom_spreadDFheat_turkey <- spreadDFheat_turkey %>% top_n(-27)

#heatmap for biomass in heat and drought
heatdrtarea <- filter (avgheatDFs, Treatment == "WL") #make a DF with all area means in the drought treatment
heatdrtarea$imgdays <- as.numeric(heatdrtarea$imgdays)

#create a heatmap with a dendrogram
plotDFheatdrt <- select(heatdrtarea, Genotype, imgdays, avg_area) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDFheatdrt$avg_area <- as.numeric(plotDFheatdrt$avg_area) #make sure the area is numeric
spreadDFheatdrt <- spread(plotDFheatdrt, key = imgdays, value = avg_area) #this flips the dataframe so the days are across the top and the genotypes are down the first row
#spreadDFheatdrt <- spreadDFheatdrt[-c(24, 25, 31, 63, 69, 101), ]
row.namesheatdrt <- spreadDFheatdrt[1:127,1] # save the row names
row.namesheatdrt <- t(row.namesheatdrt)
row.namesheatdrt <- as.character(row.namesheatdrt)
spreadDFheatdrt <- spreadDFheatdrt[,2:18] # save the DF without the first column
plot.mtxheatdrt <- as.matrix(spreadDFheatdrt) #make the DF a matrix
rownames(plot.mtxheatdrt) <- c(row.namesheatdrt) #add back row names
colnames(plot.mtxheatdrt) <- as.numeric(c(1:17)) #make the column names the days
# plot the heatmap
#heatmap <- heatmap(plot.mtx,Colv=NA)
# plot the interactive heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,2,length=100), # for orange
               seq(3,4,length=100),  # for yellow
               seq(5,40,length=100)) # for blue
pdf(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/area_heatmap_heat&drought.pdf",width = 10,height = 10)
png(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/area_heatmap_heat&drought.png",width = 900,height = 1000)
heatmapheatdrt <- heatmap.2(plot.mtxheatdrt, 
                            Colv=FALSE,
                            Rowv = FALSE,
                            dendrogram = "none",
                            xlab = "Imaging Day",
                            ylab = "Genotype",
                            main = "Heatmap of Average Plant Area in Heat & Drought",
                            col = color.palette,
                            breaks = col_breaks,
                            #margin= c(5,5),
                            srtCol =0,
                            adjCol = c(0.5,.5),
                            #adjRow = c(0.1,0.5),
                            key.ylab = NA,
                            key.xlab = "Area Value",
                            keysize=2,
                            lhei=c(2, 15),
                            lwid = c(2,6),
                            density.info = "none",
                            trace = "none") # plot the heatmap
dev.off()

spreadDFheatdrt$compare.geno <- spreadDFheatdrt$Genotype %in% list_turkey_110
spreadDFheatdrt_turkey <- filter(spreadDFheatdrt, compare.geno == TRUE)
length(unique(spreadDFheatdrt_turkey$Genotype))
spreadDFheatdrt_turkey = select(spreadDFheatdrt_turkey, -compare.geno)
spreadDFheatdrt = select(spreadDFheatdrt, -compare.geno)

spreadDFheatdrt_turkey <- select(spreadDFheatdrt_turkey, Genotype, 18)
top_spreadDFheatdrt_turkey <- spreadDFheatdrt_turkey %>% top_n(27)
bottom_spreadDFheatdrt_turkey <- spreadDFheatdrt_turkey %>% top_n(-27)


# drought: top_spreadDF1_turkey, bottom_spreadDF1_turkey
# heat: top_spreadDFheat_turkey, bottom_spreadDFheat_turkey
# heat-drought: top_spreadDFheatdrt_turkey, bottom_spreadDFheatdrt_turkey

#see which accessions overlap that perform well across stress treatments
# compare heat to heat-drought
list13 <- unique(top_spreadDFheatdrt_turkey$Genotype)
top_spreadDFheat_turkey$compare.geno <- top_spreadDFheat_turkey$Genotype %in% list13
overlap_biomass_top <- filter(top_spreadDFheat_turkey, compare.geno == TRUE)
length(unique(overlap_biomass_top$Genotype))
overlap_biomass_top = select(overlap_biomass_top, -compare.geno)
list14 <- unique(overlap_biomass_top$Genotype)
#compare heat and drought
list13 <- unique(top_spreadDF1_turkey$Genotype)
top_spreadDFheat_turkey$compare.geno <- top_spreadDFheat_turkey$Genotype %in% list13
overlap_biomass_top <- filter(top_spreadDFheat_turkey, compare.geno == TRUE)
length(unique(overlap_biomass_top$Genotype))
overlap_biomass_top = select(overlap_biomass_top, -compare.geno)
list15 <- unique(overlap_biomass_top$Genotype)
#compare heat-drought and drought
list13 <- unique(top_spreadDFheatdrt_turkey$Genotype)
top_spreadDF1_turkey$compare.geno <- top_spreadDF1_turkey$Genotype %in% list13
overlap_biomass_top <- filter(top_spreadDF1_turkey, compare.geno == TRUE)
length(unique(overlap_biomass_top$Genotype))
overlap_biomass_top = select(overlap_biomass_top, -compare.geno)
list16 <- unique(overlap_biomass_top$Genotype)
# compare heat/heat-drought to drought for overall overlap
top_spreadDF1_turkey$compare.geno <- top_spreadDF1_turkey$Genotype %in% list14
overlap_biomass_top <- filter(top_spreadDF1_turkey, compare.geno == TRUE)
length(unique(overlap_biomass_top$Genotype))


#see which accessions overlap that perform poorly across stress treatments
# compare heat to heat-drought
list17 <- unique(bottom_spreadDFheatdrt_turkey$Genotype)
bottom_spreadDFheat_turkey$compare.geno <- bottom_spreadDFheat_turkey$Genotype %in% list17
overlap_biomass_bottom <- filter(bottom_spreadDFheat_turkey, compare.geno == TRUE)
length(unique(overlap_biomass_bottom$Genotype))
overlap_biomass_bottom = select(overlap_biomass_bottom, -compare.geno)
list18 <- unique(overlap_biomass_bottom$Genotype)
#compare heat and drought
list17 <- unique(bottom_spreadDF1_turkey$Genotype)
bottom_spreadDFheat_turkey$compare.geno <- bottom_spreadDFheat_turkey$Genotype %in% list17
overlap_biomass_bottom <- filter(bottom_spreadDFheat_turkey, compare.geno == TRUE)
length(unique(overlap_biomass_bottom$Genotype))
overlap_biomass_bottom = select(overlap_biomass_bottom, -compare.geno)
list19 <- unique(overlap_biomass_bottom$Genotype)
#compare heat-drought and drought
list17 <- unique(bottom_spreadDFheatdrt_turkey$Genotype)
bottom_spreadDF1_turkey$compare.geno <- bottom_spreadDF1_turkey$Genotype %in% list17
overlap_biomass_bottom <- filter(bottom_spreadDF1_turkey, compare.geno == TRUE)
length(unique(overlap_biomass_bottom$Genotype))
overlap_biomass_bottom = select(overlap_biomass_bottom, -compare.geno)
list20 <- unique(overlap_biomass_bottom$Genotype)
# compare heat/heat-drought to drought for overall overlap
bottom_spreadDF1_turkey$compare.geno <- bottom_spreadDF1_turkey$Genotype %in% list18
overlap_biomass_bottom <- filter(bottom_spreadDF1_turkey, compare.geno == TRUE)
length(unique(overlap_biomass_bottom$Genotype))



spreadDF1 <- select(spreadDF1, Genotype, 15)
top_spreadDF1 <- spreadDF1 %>% top_n(31)
bottom_spreadDF1 <- top_spreadDF1 %>% top_n(-31)

spreadDFheat <- select(spreadDFheat, Genotype, 17)
top_spreadDFheat <- spreadDFheat %>% top_n(31)
bottom_spreadDFheat <- spreadDFheat %>% top_n(-31)

spreadDFheatdrt <- select(spreadDFheatdrt, Genotype, 18)
top_spreadDFheatdrt <- spreadDFheatdrt %>% top_n(31)
bottom_spreadDFheatdrt <- spreadDFheatdrt %>% top_n(-31)


# drought: top_spreadDF1_turkey, bottom_spreadDF1_turkey
# heat: top_spreadDFheat_turkey, bottom_spreadDFheat_turkey
# heat-drought: top_spreadDFheatdrt_turkey, bottom_spreadDFheatdrt_turkey

#see which accessions overlap that perform well across stress treatments
# compare heat to heat-drought
list43 <- unique(top_spreadDFheatdrt$Genotype)
top_spreadDFheat$compare.geno <- top_spreadDFheat$Genotype %in% list43
all_overlap_biomass_top <- filter(top_spreadDFheat, compare.geno == TRUE)
length(unique(all_overlap_biomass_top$Genotype))
all_overlap_biomass_top = select(all_overlap_biomass_top, -compare.geno)
list44 <- unique(all_overlap_biomass_top$Genotype)
#compare heat and drought
list43 <- unique(top_spreadDF1$Genotype)
top_spreadDFheat$compare.geno <- top_spreadDFheat$Genotype %in% list43
all_overlap_biomass_top <- filter(top_spreadDFheat, compare.geno == TRUE)
length(unique(all_overlap_biomass_top$Genotype))
all_overlap_biomass_top = select(all_overlap_biomass_top, -compare.geno)
list45 <- unique(all_overlap_biomass_top$Genotype)
#compare heat-drought and drought
list43 <- unique(top_spreadDFheatdrt$Genotype)
top_spreadDF1$compare.geno <- top_spreadDF1$Genotype %in% list43
all_overlap_biomass_top <- filter(top_spreadDF1, compare.geno == TRUE)
length(unique(all_overlap_biomass_top$Genotype))
all_overlap_biomass_top = select(all_overlap_biomass_top, -compare.geno)
list46 <- unique(all_overlap_biomass_top$Genotype)
# compare heat/heat-drought to drought for overall overlap
top_spreadDF1$compare.geno <- top_spreadDF1$Genotype %in% list44
all_overlap_biomass_top <- filter(top_spreadDF1, compare.geno == TRUE)
length(unique(all_overlap_biomass_top$Genotype))


#see which accessions overlap that perform poorly across stress treatments
# compare heat to heat-drought
list47 <- unique(bottom_spreadDFheatdrt$Genotype)
bottom_spreadDFheat$compare.geno <- bottom_spreadDFheat$Genotype %in% list47
all_overlap_biomass_bottom <- filter(bottom_spreadDFheat, compare.geno == TRUE)
length(unique(all_overlap_biomass_bottom$Genotype))
all_overlap_biomass_bottom = select(all_overlap_biomass_bottom, -compare.geno)
list48 <- unique(all_overlap_biomass_bottom$Genotype)
#compare heat and drought
list47 <- unique(bottom_spreadDF1$Genotype)
bottom_spreadDFheat$compare.geno <- bottom_spreadDFheat$Genotype %in% list47
all_overlap_biomass_bottom <- filter(bottom_spreadDFheat, compare.geno == TRUE)
length(unique(all_overlap_biomass_bottom$Genotype))
all_overlap_biomass_bottom = select(all_overlap_biomass_bottom, -compare.geno)
list49 <- unique(all_overlap_biomass_bottom$Genotype)
#compare heat-drought and drought
list47 <- unique(bottom_spreadDFheatdrt$Genotype)
bottom_spreadDF1$compare.geno <- bottom_spreadDF1$Genotype %in% list47
all_overlap_biomass_bottom <- filter(bottom_spreadDF1, compare.geno == TRUE)
length(unique(all_overlap_biomass_bottom$Genotype))
all_overlap_biomass_bottom = select(all_overlap_biomass_bottom, -compare.geno)
list50 <- unique(all_overlap_biomass_bottom$Genotype)
# compare heat/heat-drought to drought for overall overlap
bottom_spreadDF1$compare.geno <- bottom_spreadDF1$Genotype %in% list48
all_overlap_biomass_bottom <- filter(bottom_spreadDF1, compare.geno == TRUE)
length(unique(all_overlap_biomass_bottom$Genotype))




###########################################
# height heatmaps of heat and heat-drought conditions #
###########################################

#calculate average area and create avg df
avgheatDFs_height <- group_by(heat_127, Genotype, Treatment, imgdays) %>% summarise(avg_height = mean(cor_height_above_reference))
avgheatDFs_height <- ungroup(avgheatDFs_height)

#height heatmap for heat
heatheight <- filter(avgheatDFs_height, Treatment == "WW") #make a DF with all images in the control treatment
heatheight$imgdays <- as.numeric(heatheight$imgdays)

#remove imgday 3 (data super messy)
heatheight <- heatheight[!(heatheight$imgdays=="3"),]

#create a heatmap 
plotDFheatheight <- select(heatheight, Genotype, imgdays, avg_height) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDFheatheight$avg_height <- as.numeric(plotDFheatheight$avg_height) #make sure the area is numeric
spreadDFheatheight <- spread(plotDFheatheight, key = imgdays, value = avg_height) #this flips the dataframe so the days are across the top and the genotypes are down the first row
#spreadDFheat <- spreadDFheat[-c(24, 25, 31, 63, 69, 101), ]
row.namesheatheight <- spreadDFheatheight[1:127,1] # save the row names
row.namesheatheight <- t(row.namesheatheight)
row.namesheatheight <- as.character(row.namesheatheight)
spreadDFheatheight <- spreadDFheatheight[,2:17] # save the DF without the first column
plot.mtxheatheight <- as.matrix(spreadDFheatheight) #make the DF a matrix
rownames(plot.mtxheatheight) <- c(row.namesheatheight) #add back row names
colnames(plot.mtxheatheight) <- as.numeric(c(2:17)) #make the column names the days
# plot the heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,4,length=20), # for orange
               seq(5,7,length=20),  # for yellow
               seq(8,20,length=20)) # for blue
pdf(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/height_heatmap_heat.pdf",width = 10,height = 10)
png(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/height_heatmap_heat.png",width = 900,height = 1000)
heatmapheatheight <- heatmap.2(plot.mtxheatheight, 
                               Colv=FALSE,
                               Rowv = FALSE,
                               dendrogram = "none",
                               xlab = "Imaging Day",
                               ylab = "Genotype",
                               main = "Average Plant Height in Heat",
                               col = color.palette,
                               breaks = col_breaks,
                               #margin= c(5,5),
                               srtCol =0,
                               adjCol = c(0.5,.5),
                               #adjRow = c(0.1,0.5),
                               key.ylab = NA,
                               key.xlab = "Area Value",
                               keysize=2,
                               lhei=c(2, 15),
                               lwid = c(2,6),
                               density.info = "none",
                               trace = "none") # plot the heatmap
dev.off()

spreadDFheatheight$compare.geno <- spreadDFheatheight$Genotype %in% list_turkey_110
spreadDFheatheight_turkey <- filter(spreadDFheatheight, compare.geno == TRUE)
length(unique(spreadDFheatheight_turkey$Genotype))
spreadDFheatheight_turkey = select(spreadDFheatheight_turkey, -compare.geno)
spreadDFheatheight = select(spreadDFheatheight, -compare.geno)

spreadDFheatheight_turkey <- select(spreadDFheatheight_turkey, Genotype, 17)
top_spreadDFheatheight_turkey <- spreadDFheatheight_turkey %>% top_n(27)
bottom_spreadDFheatheight_turkey <- spreadDFheatheight_turkey %>% top_n(-27)


#heatmap for height in heat and drought
heatdrtheight <- filter (avgheatDFs_height, Treatment == "WL") #make a DF with all area means in the drought treatment
heatdrtheight$imgdays <- as.numeric(heatdrtheight$imgdays)

#create a heatmap with a dendrogram
plotDFheatdrtheight <- select(heatdrtheight, Genotype, imgdays, avg_height) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDFheatdrtheight$avg_height <- as.numeric(plotDFheatdrtheight$avg_height) #make sure the area is numeric
spreadDFheatdrtheight <- spread(plotDFheatdrtheight, key = imgdays, value = avg_height) #this flips the dataframe so the days are across the top and the genotypes are down the first row
#spreadDFheatdrt <- spreadDFheatdrt[-c(24, 25, 31, 63, 69, 101), ]
row.namesheatdrtheight <- spreadDFheatdrtheight[1:127,1] # save the row names
row.namesheatdrtheight <- t(row.namesheatdrtheight)
row.namesheatdrtheight <- as.character(row.namesheatdrtheight)
spreadDFheatdrtheight <- spreadDFheatdrtheight[,2:18] # save the DF without the first column
plot.mtxheatdrtheight <- as.matrix(spreadDFheatdrtheight) #make the DF a matrix
rownames(plot.mtxheatdrtheight) <- c(row.namesheatdrtheight) #add back row names
colnames(plot.mtxheatdrtheight) <- as.numeric(c(1:17)) #make the column names the days
# plot the heatmap
#heatmap <- heatmap(plot.mtx,Colv=NA)
# plot the interactive heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,4,length=20), # for orange
               seq(5,7,length=20),  # for yellow
               seq(8,20,length=20)) # for blue
pdf(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/height_heatmap_heat&drought.pdf",width = 10,height = 10)
png(file="/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/height_heatmap_heat&drought.png",width = 900,height = 1000)
heatmapheatdrtheight <- heatmap.2(plot.mtxheatdrtheight, 
                                  Colv=FALSE,
                                  Rowv = FALSE,
                                  dendrogram = "none",
                                  xlab = "Imaging Day",
                                  ylab = "Genotype",
                                  main = "Average Plant Height in Heat & Drought",
                                  col = color.palette,
                                  breaks = col_breaks,
                                  #margin= c(5,5),
                                  srtCol =0,
                                  adjCol = c(0.5,.5),
                                  #adjRow = c(0.1,0.5),
                                  key.ylab = NA,
                                  key.xlab = "Area Value",
                                  keysize=2,
                                  lhei=c(2, 15),
                                  lwid = c(2,6),
                                  density.info = "none",
                                  trace = "none") # plot the heatmap
dev.off()



spreadDFheatdrtheight$compare.geno <- spreadDFheatdrtheight$Genotype %in% list_turkey_110
spreadDFheatdrtheight_turkey <- filter(spreadDFheatdrtheight, compare.geno == TRUE)
length(unique(spreadDFheatdrtheight_turkey$Genotype))
spreadDFheatdrtheight_turkey = select(spreadDFheatdrtheight_turkey, -compare.geno)
spreadDFheatdrtheight = select(spreadDFheatdrtheight, -compare.geno)

spreadDFheatdrtheight_turkey <- select(spreadDFheatdrtheight_turkey, Genotype, 18)
top_spreadDFheatdrtheight_turkey <- spreadDFheatdrtheight_turkey %>% top_n(27)
bottom_spreadDFheatdrtheight_turkey <- spreadDFheatdrtheight_turkey %>% top_n(-27)


# drought: top_spreadDFheight1_turkey, bottom_spreadDFheight1_turkey
# heat: top_spreadDFheatheight_turkey, bottom_spreadDFheatheight_turkey
# heat-drought: top_spreadDFheatdrtheight_turkey, bottom_spreadDFheatdrtheight_turkey

#see which accessions overlap that perform well across stress treatments
# compare heat to heat-drought
list21 <- unique(top_spreadDFheatdrtheight_turkey$Genotype)
top_spreadDFheatheight_turkey$compare.geno <- top_spreadDFheatheight_turkey$Genotype %in% list21
overlap_height_top <- filter(top_spreadDFheatheight_turkey, compare.geno == TRUE)
length(unique(overlap_height_top$Genotype))
overlap_height_top = select(overlap_height_top, -compare.geno)
list22 <- unique(overlap_height_top$Genotype)
#compare heat and drought
list21 <- unique(top_spreadDFheight1_turkey$Genotype)
top_spreadDFheatheight_turkey$compare.geno <- top_spreadDFheatheight_turkey$Genotype %in% list21
overlap_height_top <- filter(top_spreadDFheatheight_turkey, compare.geno == TRUE)
length(unique(overlap_height_top$Genotype))
overlap_height_top = select(overlap_height_top, -compare.geno)
list23 <- unique(overlap_height_top$Genotype)
#compare heat-drought and drought
list21 <- unique(top_spreadDFheatdrtheight_turkey$Genotype)
top_spreadDFheight1_turkey$compare.geno <- top_spreadDFheight1_turkey$Genotype %in% list21
overlap_height_top <- filter(top_spreadDFheight1_turkey, compare.geno == TRUE)
length(unique(overlap_height_top$Genotype))
overlap_height_top = select(overlap_height_top, -compare.geno)
list24 <- unique(overlap_height_top$Genotype)
# compare heat/heat-drought to drought for overall overlap
top_spreadDFheight1_turkey$compare.geno <- top_spreadDFheight1_turkey$Genotype %in% list22
overlap_height_top <- filter(top_spreadDFheight1_turkey, compare.geno == TRUE)
length(unique(overlap_height_top$Genotype))


#see which accessions overlap that perform poorly across stress treatments
# compare heat to heat-drought
list21 <- unique(bottom_spreadDFheatdrtheight_turkey$Genotype)
bottom_spreadDFheatheight_turkey$compare.geno <- bottom_spreadDFheatheight_turkey$Genotype %in% list21
overlap_height_bottom <- filter(bottom_spreadDFheatheight_turkey, compare.geno == TRUE)
length(unique(overlap_height_bottom$Genotype))
overlap_height_bottom = select(overlap_height_bottom, -compare.geno)
list25 <- unique(overlap_height_bottom$Genotype)
#compare heat and drought
list21 <- unique(bottom_spreadDFheight1_turkey$Genotype)
bottom_spreadDFheatheight_turkey$compare.geno <- bottom_spreadDFheatheight_turkey$Genotype %in% list21
overlap_height_bottom <- filter(bottom_spreadDFheatheight_turkey, compare.geno == TRUE)
length(unique(overlap_height_bottom$Genotype))
overlap_height_bottom = select(overlap_height_bottom, -compare.geno)
list26 <- unique(overlap_height_bottom$Genotype)
#compare heat-drought and drought
list21 <- unique(bottom_spreadDFheatdrtheight_turkey$Genotype)
bottom_spreadDFheight1_turkey$compare.geno <- bottom_spreadDFheight1_turkey$Genotype %in% list21
overlap_height_bottom <- filter(bottom_spreadDFheight1_turkey, compare.geno == TRUE)
length(unique(overlap_height_bottom$Genotype))
overlap_height_bottom = select(overlap_height_bottom, -compare.geno)
list27 <- unique(overlap_height_bottom$Genotype)
# compare heat/heat-drought to drought for overall overlap
bottom_spreadDFheight1_turkey$compare.geno <- bottom_spreadDFheight1_turkey$Genotype %in% list25
overlap_height_bottom <- filter(bottom_spreadDFheight1_turkey, compare.geno == TRUE)
length(unique(overlap_height_bottom$Genotype))

# Get overall list of accessions that perform well for ANY phenotype
overlap_any_top <- merge(overlap_biomass_top, overlap_height_top, by = "Genotype", all = TRUE)
overlap_any_top <- merge(overlap_any_top, overlap_stressed_bottom, by = "Genotype", all = TRUE)
overlap_any_top <- select(overlap_any_top, Genotype)

# Get overall list of accessions that perform poorly for ANY phenotype
overlap_any_bottom <- merge(overlap_biomass_bottom, overlap_height_bottom, by = "Genotype", all = TRUE)
overlap_any_bottom <- merge(overlap_any_bottom, overlap_stressed_top, by = "Genotype", all = TRUE)
overlap_any_bottom <- select(overlap_any_bottom, Genotype)

# Get overall list of accessions that perform well for ALL phenotypes
overlap_all_top <- merge(overlap_biomass_top, overlap_height_top, by = "Genotype", all = FALSE)
overlap_all_top <- merge(overlap_all_top, overlap_stressed_bottom, by = "Genotype", all = FALSE)
overlap_all_top <- select(overlap_all_top, Genotype)

# Get overall list of accessions that perform poorly for ALL phenotypes
overlap_all_bottom <- merge(overlap_biomass_bottom, overlap_height_bottom, by = "Genotype", all = FALSE)
overlap_all_bottom <- merge(overlap_all_bottom, overlap_stressed_top, by = "Genotype", all = FALSE)
overlap_all_bottom <- select(overlap_all_bottom, Genotype)




spreadDFheight1 <- select(spreadDFheight1, Genotype, 15)
top_spreadDFheight1 <- spreadDFheight1 %>% top_n(31)
bottom_spreadDFheight1 <- top_spreadDFheight1 %>% top_n(-31)

spreadDFheatheight <- select(spreadDFheatheight, Genotype, 17)
top_spreadDFheatheight <- spreadDFheatheight %>% top_n(31)
bottom_spreadDFheatheight <- spreadDFheatheight %>% top_n(-31)

spreadDFheatdrtheight <- select(spreadDFheatdrtheight, Genotype, 18)
top_spreadDFheatdrtheight <- spreadDFheatdrtheight %>% top_n(31)
bottom_spreadDFheatdrtheight <- spreadDFheatdrtheight %>% top_n(-31)


# drought: top_spreadDFheight1_turkey, bottom_spreadDFheight1_turkey
# heat: top_spreadDFheatheight_turkey, bottom_spreadDFheatheight_turkey
# heat-drought: top_spreadDFheatdrtheight_turkey, bottom_spreadDFheatdrtheight_turkey

#see which accessions overlap that perform well across stress treatments
# compare heat to heat-drought
list51 <- unique(top_spreadDFheatdrtheight$Genotype)
top_spreadDFheatheight$compare.geno <- top_spreadDFheatheight$Genotype %in% list51
all_overlap_height_top <- filter(top_spreadDFheatheight, compare.geno == TRUE)
length(unique(all_overlap_height_top$Genotype))
all_overlap_height_top = select(all_overlap_height_top, -compare.geno)
list52 <- unique(all_overlap_height_top$Genotype)
#compare heat and drought
list51 <- unique(top_spreadDFheight1$Genotype)
top_spreadDFheatheight$compare.geno <- top_spreadDFheatheight$Genotype %in% list51
all_overlap_height_top <- filter(top_spreadDFheatheight, compare.geno == TRUE)
length(unique(all_overlap_height_top$Genotype))
all_overlap_height_top = select(all_overlap_height_top, -compare.geno)
list53 <- unique(all_overlap_height_top$Genotype)
#compare heat-drought and drought
list51 <- unique(top_spreadDFheatdrtheight$Genotype)
top_spreadDFheight1$compare.geno <- top_spreadDFheight1$Genotype %in% list51
all_overlap_height_top <- filter(top_spreadDFheight1, compare.geno == TRUE)
length(unique(all_overlap_height_top$Genotype))
all_overlap_height_top = select(all_overlap_height_top, -compare.geno)
list54 <- unique(all_overlap_height_top$Genotype)
# compare heat/heat-drought to drought for overall overlap
top_spreadDFheight1$compare.geno <- top_spreadDFheight1$Genotype %in% list52
all_overlap_height_top <- filter(top_spreadDFheight1, compare.geno == TRUE)
length(unique(all_overlap_height_top$Genotype))


#see which accessions overlap that perform poorly across stress treatments
# compare heat to heat-drought
list51 <- unique(bottom_spreadDFheatdrtheight$Genotype)
bottom_spreadDFheatheight$compare.geno <- bottom_spreadDFheatheight$Genotype %in% list51
all_overlap_height_bottom <- filter(bottom_spreadDFheatheight, compare.geno == TRUE)
length(unique(all_overlap_height_bottom$Genotype))
all_overlap_height_bottom = select(all_overlap_height_bottom, -compare.geno)
list55 <- unique(all_overlap_height_bottom$Genotype)
#compare heat and drought
list51 <- unique(bottom_spreadDFheight1$Genotype)
bottom_spreadDFheatheight$compare.geno <- bottom_spreadDFheatheight$Genotype %in% list51
all_overlap_height_bottom <- filter(bottom_spreadDFheatheight, compare.geno == TRUE)
length(unique(all_overlap_height_bottom$Genotype))
all_overlap_height_bottom = select(all_overlap_height_bottom, -compare.geno)
list56 <- unique(all_overlap_height_bottom$Genotype)
#compare heat-drought and drought
list51 <- unique(bottom_spreadDFheatdrtheight$Genotype)
bottom_spreadDFheight1$compare.geno <- bottom_spreadDFheight1$Genotype %in% list51
all_overlap_height_bottom <- filter(bottom_spreadDFheight1, compare.geno == TRUE)
length(unique(all_overlap_height_bottom$Genotype))
all_overlap_height_bottom = select(all_overlap_height_bottom, -compare.geno)
list57 <- unique(all_overlap_height_bottom$Genotype)
# compare heat/heat-drought to drought for overall overlap
bottom_spreadDFheight1$compare.geno <- bottom_spreadDFheight1$Genotype %in% list55
all_overlap_height_bottom <- filter(bottom_spreadDFheight1, compare.geno == TRUE)
length(unique(all_overlap_height_bottom$Genotype))

# Get overall list of accessions that perform well for ANY phenotype
all_overlap_any_top <- merge(all_overlap_biomass_top, all_overlap_height_top, by = "Genotype", all = TRUE)
all_overlap_any_top <- merge(all_overlap_any_top, all_overlap_stressed_bottom, by = "Genotype", all = TRUE)
all_overlap_any_top <- select(all_overlap_any_top, Genotype)

# Get overall list of accessions that perform poorly for ANY phenotype
all_overlap_any_bottom <- merge(all_overlap_biomass_bottom, all_overlap_height_bottom, by = "Genotype", all = TRUE)
all_overlap_any_bottom <- merge(all_overlap_any_bottom, all_overlap_stressed_top, by = "Genotype", all = TRUE)
all_overlap_any_bottom <- select(all_overlap_any_bottom, Genotype)

# Get overall list of accessions that perform well for ALL phenotypes
all_overlap_all_top <- merge(all_overlap_biomass_top, all_overlap_height_top, by = "Genotype", all = FALSE)
all_overlap_all_top <- merge(all_overlap_all_top, all_overlap_stressed_bottom, by = "Genotype", all = FALSE)
all_overlap_all_top <- select(all_overlap_all_top, Genotype)

# Get overall list of accessions that perform poorly for ALL phenotypes
all_overlap_all_bottom <- merge(all_overlap_biomass_bottom, all_overlap_height_bottom, by = "Genotype", all = FALSE)
all_overlap_all_bottom <- merge(all_overlap_all_bottom, all_overlap_stressed_top, by = "Genotype", all = FALSE)
all_overlap_all_bottom <- select(all_overlap_all_bottom, Genotype)





#compare genotypes in control and heat treatments and keep only the ones in both
newDF <- merge(avgDFs, avgheatDFs, by.x = "Genotype", by.y = "Genotype", all = FALSE)
length(unique(newDF$Genotype))
list1 <- unique(newDF$Genotype)
avgDFs$compare.geno <- avgDFs$Genotype %in% list1
avgheatDFs$compare.geno <- avgheatDFs$Genotype %in% list1
avgDFs <- filter(avgDFs, compare.geno == TRUE)
avgheatDFs <- filter(avgheatDFs, compare.geno == TRUE)
length(unique(avgDFs$Genotype))

#remove -compare.geno column from dataframes
avgDFs = select(avgDFs, -compare.geno)
avgheatDFs = select(avgheatDFs, -compare.geno)




heatDFs1 <- heatDFs

heatDFs1 <- select(heatDFs1, Latitude, Longitude, Elevation, hue_circular_mean, cor_area, cor_height_above_reference, cor_width, cor_height, cor_height_below_reference, cor_area_above_reference, cor_area_below_reference)
heatDFs1$Elevation <- as.numeric(heatDFs1$Elevation)






###############################################
### geographical & climate data ###
###############################################

#download package from github
#devtools::install_github("kapitzas/WorldClimTiles")
#require package (also need raster package)
require("WorldClimTiles")

#use country GADM boundaries to see which tiles we need to download from WorldClim (Spain and Turkey)
#Level 0 indicates country level boundaries, higher number is more specific
esp <- getData("GADM", country = "ESP", level = 0)
tilenames <- tile_name(esp)
tilenames

tur <- getData("GADM", country = "TUR", level = 0)
tilenames2 <- tile_name(tur)
tilenames2

#create list of tiles that you need and want to download and then merge
tilenames_merged <- c("15","16","17")
tilenames_merged <- as.character(tilenames_merged)

#download tiles from WorldClim
tiles <- tile_get(tilenames_merged, "bio")

#merge tiles (LONG STEP, takes about 22 mins per layer per tile so about 13 hrs total)
merged <- tile_merge(tiles)

#copy merged tiles so don't have to repeat merging step in case something goes wrong
merged_copy <- merged

#Plot one raster layer of merged tiles to check it
plot(merged[[1]], main="Annual Mean Temperature")

#save merged raster to disk 
writeRaster(merged, filename = "mergedtiles.tif", format="GTiff", overwrite=TRUE)

#read in merged raster
merged <- stack("mergedtiles.tif")

# Define the extent of the study region we want
e <- extent(c(-7, 45, 35, 47))
# crop worldclim data to this extent
climate <- crop(merged, e)

#plot one raster layer to check cropping
plot(climate[[1]], main="Annual Mean Temperature")

#Convert climate from rasterstack format to dataframe format
climate <- as.data.frame(climate, xy=TRUE)
#climate = drop_na(climate)

#Combine climate data with location data
locations = drop_na((locations))
#climate_sub <- climate[climate$x %between% c(,) & climate$y %between% c(,),]
#Calculate the euclidian distance between the locations and climate data points
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

#bind locations df with back_df
locations_new <- cbind(locations,back_df)
head(locations_new)
filtered <- locations_new[locations_new$euclid < sqrt(2),]

#make a histogram of the euclidian distance between collection locations and climate data points
hist(back_df$euclid)

#keep only genotypes in locations_new that are in both experiments
newDF <- merge(avgDFs, avgheatDFs, by.x = "Genotype", by.y = "Genotype", all = FALSE)
length(unique(newDF$Genotype))
list1 <- unique(newDF$Genotype)
locations_new$compare.geno <- locations_new$Genotype %in% list1
locations_new <- filter(locations_new, compare.geno == TRUE)
length(unique(locations_new$Genotype))
locations_new = select(locations_new, -compare.geno)

#Split collection_location column to have country alone in column
locations_new2 <- locations_new
locations_new$Collection_location <- as.character(locations_new$Collection_location)
locations_new <- separate(data= locations_new, col=Collection_location, into= c("A", "B", "C", "Country"), fill = "left", sep = "\\, ")

# Save csv with location data of final 127 genotypes
locations_new$compare.geno <- locations_new$Genotype %in% list3
locations_127 <- filter(locations_new, compare.geno == TRUE)
length(unique(locations_127$Genotype))
locations_127 = select(locations_127, -compare.geno)
write.csv(locations_127, "./brachy_final_locations.csv", row.names = FALSE)


#Make correlation plot with climate variables
bio <- select(back_df, -x,-y,-diff_x,-diff_y,-euclid)
bio_cor <- cor(bio)
pdf(file= "/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/climatecorrelationplot.pdf", width = 6, height = 7)
png(file= "/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/climatecorrelationplot.png", width = 1200, height = 1400)
corrplot(bio_cor, method = "circle", type = "lower", order = "FPC",
         tl.col="black", tl.srt=45, tl.cex = 3,cl.cex = 3,
         title= "Correlation Plot of Climate Variables")

#climatecorr
dev.off()


# Do PCA with climate data
bio.pca = PCA(locations_new[,11:29], scale.unit=TRUE, ncp=5, graph=T)

pca_df = data.frame(
  "Genotype" = locations_new$Genotype,
  "Country" = locations_new$Country,
  "Latitude" = locations_new$Latitude,
  "Longitude" = locations_new$Longitude,
  "PC1" = bio.pca$ind$coord[, 1],
  "PC2" = bio.pca$ind$coord[, 2],
  "PC3" = bio.pca$ind$coord[, 3],
  "PC4" = bio.pca$ind$coord[, 4])

percentofvariance = data.frame(bio.pca$eig) 
percentofvariance$PC <- row.names(percentofvariance)
row.names(percentofvariance) <- NULL
percentofvariance$PC <- row.names(percentofvariance)
percentofvariance$PC <- as.factor(percentofvariance$PC)
percentofvariance$percentage.of.variance <- as.numeric(percentofvariance$percentage.of.variance)

#plot scree plot
ggplot(percentofvariance, aes(PC, percentage.of.variance, group=1))+
  ggtitle("Scree Plot")+
  geom_line()+
  geom_point()

#calculate loadings for climate PCA
# loadings are significant when > sqrt(1/19) which = 0.229415733870562 (19 is # of climate variables in PCA)
climate.loadings <- sweep(bio.pca$var$coord,2,sqrt(bio.pca$eig[1:ncol(bio.pca$var$coord),1]),FUN="/")



# limit pca_df to only genotypes with all data
pca_df$compare.geno <- pca_df$Genotype %in% list3
pca_df <- filter(pca_df, compare.geno == TRUE)
length(unique(pca_df$Genotype))
pca_df = select(pca_df, -compare.geno)

# make colorblind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# To use for fills, add
scale_fill_manual(values=cbPalette)
# To use for line and point colors, add
scale_colour_manual(values=cbPalette)


climate_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Country))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  theme_minimal()+
  theme(panel.grid = element_blank(), panel.border = element_rect(fill="transparent"))+
  scale_colour_manual(values=cbPalette)+
  stat_ellipse(geom = "polygon", aes(fill = Country), alpha= 0.2, show.legend = FALSE, level = 0.95)+
  xlab("PC 1 (50.0%)")+
  ylab("PC 2 (22.7%)")+
  ggtitle("PCA of Accession Collection Location Climate")+
  geom_point()
climate_pca
ggsave(plot = climate_pca,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/climatepca.png", width = 7.25, height = 6, dpi = 300)
ggsave(plot = climate_pca,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/climatepca.pdf", width = 7.25, height = 6, dpi = 300)


# Graph PCA colored by collection country
fviz_pca_ind(bio.pca,
             geom.ind = "point", # show points only (nbut not "text")
             mean.point = FALSE,
             pointsize = 2.5,
             #palette = "jco",
             ggtheme = theme_gray(),
             col.ind = locations_new$Country, # color by groups
             legend.title = "Country"
)


# Graph PCA colored by collection elevation
locations_new$Elevation <- as.numeric(locations_new$Elevation) #make sure elevation is numeric

fviz_pca_ind(bio.pca,
             geom.ind = "point", # show points only
             mean.point = FALSE,
             pointsize = 2.5,
             ggtheme = theme_gray(),
             col.ind = locations_new$Elevation, # color by elevation
             legend.title = "Elevation"
) +
  scale_color_gradient2(low = "white", mid = "blue", high = "red", midpoint = 1600)



############################################
# make map of population collection locations
#############################################


library(tidyverse)
library(sf)
#library(leaflet)
library(mapview)
library(ggspatial)



# limit locations_new to only genotypes with all data
locations_new$compare.geno <- locations_new$Genotype %in% list3
locations_new <- filter(locations_new, compare.geno == TRUE)
length(unique(locations_new$Genotype))
locations_new = select(locations_new, -compare.geno)

mapview(locations_new, xcol = "Longitude", ycol = "Latitude", crs = 4269, grid = FALSE,legend = TRUE, col)


world_shape <- st_read("ne_10m_admin_0_countries.shp")
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
ggsave(plot = map,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/collectionlocationmap.png", width = 6, height = 4, dpi = 300)
ggsave(plot = map,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/collectionlocationmap.pdf", width = 6, height = 4, dpi = 300)





############################################################
### PCA of Population Structure ###
############################################################
# Read in GBS structure data
gbs_pca <- read.table("brachypca.eigenvec", sep = " ", stringsAsFactors = FALSE)
# Change ARC1 to Arc1 to match locations genotype format
gbs_pca[23, 2] = "Arc1"

# See what PCs 1 and 2 for GBS look like
plot(gbs_pca[,3],gbs_pca[,4])

#keep only genotypes in gbs_pca that are in both experiments
list2 <- unique(locations_new$Genotype)
gbs_pca$compare.geno <- gbs_pca$V2 %in% list2
gbs_pca <- filter(gbs_pca, compare.geno == TRUE)
length(unique(gbs_pca$V2))
gbs_pca = select(gbs_pca, -compare.geno)

#keep only genotypes in locations_new that are in both experiments AND have GBS data
list2 <- unique(gbs_pca$V2)
locations_new$compare.geno <- locations_new$Genotype %in% list2
locations_new <- filter(locations_new, compare.geno == TRUE)
length(unique(locations_new$Genotype))
locations_new = select(locations_new, -compare.geno)

gbs_pca$compare.geno <- gbs_pca$V2 %in% list3

# Check where in PCA the genotypes are that aren't in 127 with all data to make sure not all clustered together
ggplot(gbs_pca, aes(V3,V4))+
  geom_point(aes(color=compare.geno))

# Rename Genotype column in GBS data
colnames(gbs_pca)[2] <- "Genotype"

# Join GBS PCA data with climate PCA data
joined_pca_df <- join(pca_df, gbs_pca[,2:5], by="Genotype")

head(joined_pca_df)

# Plot population structure PCA
pop_pca <- ggplot(joined_pca_df, aes(x = V3, y = V4, color = Country))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  theme_minimal()+
  theme(panel.grid = element_blank(), panel.border = element_rect(fill="transparent"))+
  scale_colour_manual(values=cbPalette)+
  stat_ellipse(geom = "polygon", aes(fill = Country), alpha= 0.2, show.legend = FALSE, level = 0.95)+
  xlab("PC 1")+
  ylab("PC 2")+
  ggtitle("PCA of Population Structure")+
  geom_text(label=joined_pca_df$Genotype, nudge_x = .01,check_overlap = T)+
  geom_point()
pop_pca
ggsave(plot = pop_pca,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/new_poppca.png", width = 7.25, height = 6, dpi = 300)
ggsave(plot = pop_pca,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/new_poppca.pdf", width = 7.25, height = 6, dpi = 300)

ggplotly(pop_pca)

###############################################################
### Height GWAS/Heritability ### control WW and WL
###############################################################


# Keep only genotypes that have WW and WL data on last day (comes to 127)
list3 <- unique(heatarea_dif$Genotype)
DFs$compare.geno <- DFs$Genotype %in% list3
drought_127 <- filter(DFs, compare.geno == TRUE)
length(unique(drought_127$Genotype))
drought_127 = select(drought_127, -compare.geno)
DFs = select(DFs, -compare.geno)

library(lme4)

#Broad sense heritability section
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
dat2 <- drought_127[drought_127$imgdays == 15,]
colnames(dat2)[13] <- "area"
colnames(dat2)[14] <- "width"
colnames(dat2)[16] <- "height_above_reference"
colnames(dat2)[17] <- "height_below_reference"
colnames(dat2)[19] <- "percent_area_above_reference"
colnames(dat2)[21] <- "percent_area_below_reference"
tail(dat2)
shapes <- c("hue_circular_mean", "hue_circular_std", "hue_median", "percent_unhealthy",
            "area", "width", "height_above_reference", "height_below_reference", 
            "percent_area_above_reference", "percent_area_below_reference")
H3 <- c()
for(e in shapes){
  model2 <- lmer(eval(parse(text=e))~(1|Genotype)+(1|Treatment)+(1|Genotype:Treatment),data = dat2)
  re2<-as.numeric(VarCorr(model2))
  res2<-attr(VarCorr(model2), "sc")^2
  interaction.var2 <- re2[1]
  genotype.var2<-re2[2]
  treatment.var2<-re2[3]
  tot.var2<-sum(re2,res2)
  unexp2 <- 1-sum(re2)/sum(re2,res2)
  h3 <- c((genotype.var2/tot.var2),
          (treatment.var2/tot.var2),
          (interaction.var2/tot.var2),
          unexp2)
  H3 <- rbind(H3,h3)
}
H3 <- data.frame(H3,row.names = shapes)
H3$Shape <- rownames(H3)
rownames(H3) <- NULL
colnames(H3) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
H3$Shape <-  ordered(H3$Shape,levels=H3$Shape[order(H3$Unexplained)])
H3_melt <- melt(H3,id=c("Shape"))
H3_melt$variable <- ordered(H3_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(H3_melt)

library(scales)
p2 <- ggplot(data=H3_melt,aes(Shape,value*100))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Variance Explained between Control WW and WL")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y= element_text(size = 12),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill=NA, size=1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
p2
ggsave("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/variance_explained_controlWW&WL.png",width=5.5,height=5.54,plot = p2, dpi = 300)
ggsave("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/variance_explained_controlWW&WL.pdf",width=7.17,height=5.54,plot = p2, dpi = 300)


# GWAS Section
avgheight <- group_by(drought_127, Genotype, Treatment, imgdays) %>% summarise(avg_height = mean(cor_height_above_reference))
avgheight$imgdays <- as.numeric(avgheight$imgdays)

avg_day15 <- aggregate(data = avgheight[avgheight$imgdays == 15,], avg_height~Genotype+Treatment, FUN= function(i)mean(i))
head(avg_day15)

table(avg_day15$Genotype, avg_day15$Treatment)

height_dif <- setNames(data.frame(do.call("rbind",lapply(split(avg_day15, avg_day15$Genotype), function(i)diff(i$avg_height)))), c("height_dif2"))
height_dif$Genotype <- row.names(height_dif)
head(height_dif)
row.names(height_dif) <- NULL

joined_pca_df2 <- join(joined_pca_df, height_dif, by="Genotype", type="inner")
head(joined_pca_df2)

data <- data.table::fread("populations.structure",skip = 1)
my_names <- as.character(read.table("populations.structure",nrows = 1,sep="\t"))

#rename genotype column for gbs data
colnames(data)[1] <- "Genotype"
data[45, 1] = "Arc1"
data[46, 1] = "Arc1"
data = subset(data, select = -c(V2))

# Keep only final 127 genotypes with all data)
data$compare.geno <- data$Genotype %in% list3
snp_data <- filter(data, compare.geno == TRUE)
length(unique(snp_data$Genotype))
snp_data = select(snp_data, -compare.geno)
data = select(data, -compare.geno)

joined_pca_df_final <- join(joined_pca_df2, data, by="Genotype", type="inner")
write.csv(joined_pca_df_final, "./brachy_heat_final_joined_gbs_data.csv", row.names = FALSE)
joined_pca_df_final[1:5,1:20]

library(lme4)
str(joined_pca_df_final)
joined_pca_df_final$`96_83`
joined_pca_df_final$snp <- joined_pca_df_final[,"96_83"]

mod2 <- lmer(data= joined_pca_df_final, height_dif2~0+as.factor(snp)+ (1|V3)+(1|V4)+(1|V5)+(1|PC1))

my_test2 <- Anova(mod2)
my_test2$`Pr(>Chisq)`

library(parallel)
snp_names <- colnames(joined_pca_df_final)[13:ncol(joined_pca_df_final)]

system.time({my_gwas2 <- mclapply(snp_names, mc.cores = 4, function(i){
  tryCatch({
    joined_pca_df_final$snp <- joined_pca_df_final[,i]
    suppressWarnings(mod2 <- lmer(data= joined_pca_df_final, height_dif2~as.factor(snp)+ (1|V3)+(1|V4)+(1|V5)+(1|PC1)))
    suppressWarnings(my_test2 <- Anova(mod2))
    my_test2$`Pr(>Chisq)`
  },warning=function(war){},error=function(err){})
})
})

names(my_gwas2) <- snp_names

#copy my_gwas2 so don't have to re-run gwas
my_gwas2_copy <- my_gwas2

my_gwas2[unlist(lapply(my_gwas2, function(i) !is.null(i)))]

plot(-log10(p.adjust(unlist(my_gwas2),method="bonferroni")),pch=16,col="grey20",main="WW and WL in Control",xlab="SNP",ylab="-log10(P)")
abline(h=-log10(0.00001), col="orange")

# convert to dataframe
my_gwas2_copy.df <- ldply(unlist(my_gwas2_copy), data.frame)

#remove any row names and rename columns
row.names(my_gwas2_copy.df) <- NULL
colnames(my_gwas2_copy.df)[1] <- "ID"
colnames(my_gwas2_copy.df)[2] <- "Value"


### Genomic Inflation Factor Correction ###
chisq <- qchisq(1-my_gwas2_copy.df$Value,1)
lambda <- median(chisq)/qchisq(0.5,1)
my_gwas2_copy.df$P_new <- pchisq(chisq/lambda,df=1, lower.tail = F)


# read in file with chromosome info
chrominfo <- read.csv("populations.plink.csv")

# make sure ID column in both data.frames is character
chrominfo$ID <- as.character(chrominfo$ID)
chrominfo$POS <- as.numeric(chrominfo$POS)

my_gwas2_copy.df$ID <- as.character(my_gwas2_copy.df$ID)
my_gwas2_copy.df$Value <- as.numeric(my_gwas2_copy.df$Value)


# join gwas df with chromosome info
gwas_joined <- join(chrominfo, my_gwas2_copy.df, by="ID", type="right")

#remove centromere snp
gwas_joined <- gwas_joined[!(gwas_joined$CHROM=="Bd1_centromere_containing_Bradi1g41430"),]

# Rename chromosomes to numbers
gwas_joined[gwas_joined=="Bd1"]<-"1"
gwas_joined[gwas_joined=="Bd2"]<-"2"
gwas_joined[gwas_joined=="Bd3"]<-"3"
gwas_joined[gwas_joined=="Bd4"]<-"4"
gwas_joined[gwas_joined=="Bd5"]<-"5"

gwas_joined$CHROM <- as.numeric(gwas_joined$CHROM)
gwas_joined$P_new <- as.numeric(gwas_joined$P_new)

gwas_joined = drop_na(gwas_joined)
unique(gwas_joined$CHROM)

gwas_joined$p.adj <- p.adjust(gwas_joined$Value, method = "bonferroni")

library(qqman)
gwas_joined = drop_na(gwas_joined)
# Make the Manhattan plot on the gwas_joined dataset
manhattan(gwas_joined, chr="CHROM", bp="POS", snp="ID", p="Value", annotatePval = .01 )
# Plot manhattan of genomic inflation corrected p-values
manhattan(gwas_joined, chr="CHROM", bp="POS", snp="ID", p="P_new", annotatePval = .01 )
# Plot manhattan of bonferroni corrected genomic inflation corrected p-values
manhattan(gwas_joined, chr="CHROM", bp="POS", snp="ID", p="p.adj", annotatePval = .01 )
gwas_joined[-log10(gwas_joined$p.adj)>15,]
# plot QQ plots of all 3 versions
qq(gwas_joined$Value)
qq(gwas_joined$P_new)
qq(gwas_joined$p.adj)



don <- gwas_joined %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_joined, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( POScum=POS+tot)

axisdf = don %>% group_by(CHROM) %>% summarize(center=(max(POScum) + min(POScum) ) / 2 )

manhattan <- ggplot(don, aes(x=POScum, y=-log10(p.adj))) +
  # Show all points
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black", "grey"), 5 )) +
  # custom axis:
  scale_x_continuous(label = axisdf$CHROM, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  expand_limits(y = c(0, 65))+
  xlab("Chromosome")+
  #significance line
  geom_hline(yintercept=-log10(5e-8),color='red')+
  geom_hline(yintercept=-log10(5e-6),color='blue')+
  ggtitle("GWAS Comparing Control to WL")+
  # Customize the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
  )
manhattan
ggsave(plot = manhattan,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/new_manhattanplot.png", width = 6, height = 4, dpi = 300)
ggsave(plot = manhattan,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/new_manhattanplot.pdf", width = 6, height = 4, dpi = 300)



###############################################################
### Height GWAS/Heritability ### control WW and heat WW
###############################################################

# make df with only WW data from heat and control temps
heatWW <- heat_127[heat_127$Treatment == "WW",]
heatWW <- heatWW[heatWW$imgdays == 17,]
heatWW[heatWW=="WW"]<-"heat"
heatWW = subset(heatWW, select = -c(image, plantbarcode, date, time, days, Treatment.2, Collection_location, Latitude, Longitude, Elevation))

droughtWW <- drought_127[drought_127$Treatment == "WW",]
droughtWW <- droughtWW[droughtWW$imgdays == 15,]
droughtWW[droughtWW=="WW"]<-"control"
droughtWW = subset(droughtWW, select = -c(image, plantbarcode, date, time, days, Treatment.2))

WW <- rbind(heatWW, droughtWW)


#Broad sense heritability section
lmer(data = heat_127[heat_127$imgdays == 17,], cor_area~1+(1|Genotype)+(1|Treatment)+(1|Genotype:Treatment))
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
dat3 <- WW
colnames(dat3)[7] <- "area"
colnames(dat3)[8] <- "width"
colnames(dat3)[10] <- "height_above_reference"
colnames(dat3)[11] <- "height_below_reference"
colnames(dat3)[13] <- "percent_area_above_reference"
colnames(dat3)[15] <- "percent_area_below_reference"
tail(dat3)
shapes <- c("hue_circular_mean", "hue_circular_std", "hue_median", "percent_unhealthy", "area", "width", "height_above_reference", "height_below_reference", "percent_area_above_reference", "percent_area_below_reference")
H4 <- c()
for(e in shapes){
  #e <- "cor_area"
  model3 <- lmer(eval(parse(text=e))~(1|Genotype)+(1|Treatment)+(1|Genotype:Treatment),data = dat3)
  re3<-as.numeric(VarCorr(model3))
  res3<-attr(VarCorr(model3), "sc")^2
  interaction.var3 <- re3[1]
  genotype.var3<-re3[2]
  treatment.var3<-re3[3]
  tot.var3<-sum(re3,res3)
  unexp3 <- 1-sum(re3)/sum(re3,res3)
  h4 <- c((genotype.var3/tot.var3),
          (treatment.var3/tot.var3),
          (interaction.var3/tot.var3),
          unexp3)
  H4 <- rbind(H4,h4)
}
H4 <- data.frame(H4,row.names = shapes)
H4$Shape <- rownames(H4)
rownames(H4) <- NULL
colnames(H4) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
H4$Shape <-  ordered(H4$Shape,levels=H4$Shape[order(H4$Unexplained)])
H4_melt <- melt(H4,id=c("Shape"))
H4_melt$variable <- ordered(H4_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(H4_melt)
library(scales)
p3 <- ggplot(data=H4_melt,aes(Shape,value*100))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Variance Explained between Control and Heat WW")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y= element_text(size = 12),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill=NA, size=1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
p3
ggsave("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/variance_explained_control&HeatWW.png",width=5.5,height=5.54,plot = p3, dpi = 300)
ggsave("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/variance_explained_control&HeatWW.pdf",width=7.17,height=5.54,plot = p3, dpi = 300)


#### GWAS Section

# make df with only WW data from heat and control temps
heatWW_avg <- heatavg_day17[heatavg_day17$Treatment == "WW",]
heatWW_avg[heatWW_avg=="WW"]<-"heat"

droughtWW_avg <- avg_day15[avg_day15$Treatment == "WW",]
droughtWW_avg[droughtWW_avg=="WW"]<-"control"

WW_avg <- rbind(heatWW_avg, droughtWW_avg)

#check to make sure all genotypes have data for both treatments on last day
table(WW_avg$Genotype, WW_avg$Treatment)

height_dif3 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avg, WW_avg$Genotype), function(i)diff(i$avg_height)))), c("height_dif3"))
height_dif3$Genotype <- row.names(height_dif3)
head(height_dif3)
row.names(height_dif3) <- NULL

joined_pca_df3 <- join(joined_pca_df2, height_dif3, by="Genotype", type="inner")
head(joined_pca_df3)

#read in gbs data
#data <- data.table::fread("populations.structure",skip = 1)
#my_names <- as.character(read.table("populations.structure",nrows = 1,sep="\t"))

#rename genotype column
#colnames(data)[1] <- "Genotype"
#data[45, 1] = "Arc1"
#data[46, 1] = "Arc1"
#data = subset(data, select = -c(V2))

joined_pca_df_final <- join(joined_pca_df3, data, by="Genotype", type="inner")
write.csv(joined_pca_df_final, "./brachy_heat_final_joined_gbs_data.csv", row.names = FALSE)
joined_pca_df_final[1:5,1:20]

library(lme4)
#str(joined_pca_df_final)
#joined_pca_df_final$`96_83`
joined_pca_df_final$snp <- joined_pca_df_final[,"96_83"]

mod3 <- lmer(data= joined_pca_df_final, height_dif3~0+as.factor(snp)+ (1|V3)+(1|V4)+(1|V5)+(1|PC1))

my_test3 <- Anova(mod3)
my_test3$`Pr(>Chisq)`

library(parallel)

snp_names <- colnames(joined_pca_df_final)[14:ncol(joined_pca_df_final)]

system.time({my_gwas3 <- mclapply(snp_names, mc.cores = 4, function(i){
  tryCatch({
    joined_pca_df_final$snp <- joined_pca_df_final[,i]
    suppressWarnings(mod3 <- lmer(data= joined_pca_df_final, height_dif3~as.factor(snp)+ (1|V3)+(1|V4)+(1|V5)+(1|PC1)))
    suppressWarnings(my_test3 <- Anova(mod3))
    my_test3$`Pr(>Chisq)`
  },warning=function(war){},error=function(err){})
})
})

#copy my_gwas3 so don't have to re-run gwas
my_gwas3_copy <- my_gwas3


names(my_gwas3) <- snp_names
my_gwas3[unlist(lapply(my_gwas3, function(i) !is.null(i)))]

#quick manhattan plot to see what results look like
plot(-log10(p.adjust(unlist(my_gwas3),method="bonferroni")),pch=16,col="grey20",main="WW in Control and Heat",xlab="SNP",ylab="-log10(P)")
abline(h=-log10(0.00001), col="orange")


# convert to dataframe
my_gwas3_copy.df <- ldply(unlist(my_gwas3_copy), data.frame)

#remove any row names and rename columns
row.names(my_gwas3_copy.df) <- NULL
colnames(my_gwas3_copy.df)[1] <- "ID"
colnames(my_gwas3_copy.df)[2] <- "Value"


### Genomic Inflation Factor Correction ###
chisq3 <- qchisq(1-my_gwas3_copy.df$Value,1)
lambda3 <- median(chisq3)/qchisq(0.5,1)
my_gwas3_copy.df$P_new <- pchisq(chisq3/lambda3, df=1,lower.tail=FALSE)


# read in file with chromosome info
#chrominfo <- read.csv("populations.plink.csv")

# make sure ID column in both data.frames is character
#chrominfo$ID <- as.character(chrominfo$ID)
#chrominfo$POS <- as.numeric(chrominfo$POS)

my_gwas3_copy.df$ID <- as.character(my_gwas3_copy.df$ID)
my_gwas3_copy.df$Value <- as.numeric(my_gwas3_copy.df$Value)


# join gwas df with chromosome info
gwas3_joined <- join(chrominfo, my_gwas3_copy.df, by="ID", type="left")

#remove centromere snp
gwas3_joined <- gwas3_joined[!(gwas3_joined$CHROM=="Bd1_centromere_containing_Bradi1g41430"),]

# Rename chromosomes to numbers
gwas3_joined[gwas3_joined=="Bd1"]<-"1"
gwas3_joined[gwas3_joined=="Bd2"]<-"2"
gwas3_joined[gwas3_joined=="Bd3"]<-"3"
gwas3_joined[gwas3_joined=="Bd4"]<-"4"
gwas3_joined[gwas3_joined=="Bd5"]<-"5"

gwas3_joined$CHROM <- as.numeric(gwas3_joined$CHROM)

unique(gwas3_joined$CHROM)

gwas3_joined$p.adj <- p.adjust(gwas3_joined$Value, method = "bonferroni")

library(qqman)
gwas3_joined = drop_na(gwas3_joined)
# Make the Manhattan plot on the gwas_joined dataset
manhattan(gwas3_joined, chr="CHROM", bp="POS", snp="ID", p="Value", annotatePval = .01 )
# Plot manhattan of genomic inflation corrected p-values
manhattan(gwas3_joined, chr="CHROM", bp="POS", snp="ID", p="P_new", annotatePval = .01 )
# Plot manhattan of bonferroni corrected genomic inflation corrected p-values
manhattan(gwas3_joined, chr="CHROM", bp="POS", snp="ID", p="p.adj", annotatePval = .01 )
# plot QQ plots of all 3 versions
qq(gwas3_joined$Value)
qq(gwas3_joined$P_new)
qq(gwas3_joined$p.adj)



don3 <- gwas3_joined %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas3_joined, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( POScum=POS+tot)

axisdf3 = don3 %>% group_by(CHROM) %>% summarize(center=(max(POScum) + min(POScum) ) / 2 )

manhattan3 <- ggplot(don3, aes(x=POScum, y=-log10(p.adj))) +
  # Show all points
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black", "grey"), 5 )) +
  # custom axis:
  scale_x_continuous(label = axisdf3$CHROM, breaks= axisdf3$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  expand_limits(y = c(0, 30))+
  xlab("Chromosome")+
  #significance line
  #geom_hline(yintercept=-log10(5e-8),color='red')+
  geom_hline(yintercept=-log10(5e-6),color='blue')+
  ggtitle("GWAS Comparing Control to WL")+
  # Customize the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
  )
manhattan3
#ggsave(plot = manhattan3,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/new_manhattanplot.png", width = 6, height = 4, dpi = 300)
#ggsave(plot = manhattan3,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/new_manhattanplot.pdf", width = 6, height = 4, dpi = 300)

###############################################################
### Height GWAS/Heritability ### control WW and heat WL
###############################################################

# make df with only WW data from control
heatWL <- heat_127[heat_127$Treatment == "WL",]
heatWL <- heatWL[heatWL$imgdays == 17,]
heatWL[heatWL=="WL"]<-"heat WL"
heatWL = subset(heatWL, select = -c(image, plantbarcode, date, time, days, Treatment.2, Collection_location, Latitude, Longitude, Elevation))

# make df with only WL data from heat
droughtWW <- drought_127[drought_127$Treatment == "WW",]
droughtWW <- droughtWW[droughtWW$imgdays == 15,]
droughtWW[droughtWW=="WW"]<-"control"
droughtWW = subset(droughtWW, select = -c(image, plantbarcode, date, time, days, Treatment.2))

control_and_heatWL <- rbind(heatWL, droughtWW)


#Broad sense heritability section
#lmer(data = heat_127[heat_127$imgdays == 17,], cor_height_above_reference~1+(1|Genotype)+(1|Treatment)+(1|Genotype:Treatment))
#*************************************************************************************************
# R^2 and Partial Correlations
#*************************************************************************************************
dat4 <- control_and_heatWL
colnames(dat4)[7] <- "area"
colnames(dat4)[8] <- "width"
colnames(dat4)[10] <- "height_above_reference"
colnames(dat4)[11] <- "height_below_reference"
colnames(dat4)[13] <- "percent_area_above_reference"
colnames(dat4)[15] <- "percent_area_below_reference"
tail(dat3)
shapes <- c("hue_circular_mean", "hue_circular_std", "hue_median", "percent_unhealthy", "area", "width", "height_above_reference", "height_below_reference", "percent_area_above_reference", "percent_area_below_reference")
H5 <- c()
for(e in shapes){
  #e <- "cor_area"
  model4 <- lmer(eval(parse(text=e))~(1|Genotype)+(1|Treatment)+(1|Genotype:Treatment),data = dat4)
  re4<-as.numeric(VarCorr(model4))
  res4<-attr(VarCorr(model4), "sc")^2
  interaction.var4 <- re4[1]
  genotype.var4<-re4[2]
  treatment.var4<-re4[3]
  tot.var4<-sum(re4,res4)
  unexp4 <- 1-sum(re4)/sum(re4,res4)
  h5 <- c((genotype.var4/tot.var4),
          (treatment.var4/tot.var4),
          (interaction.var4/tot.var4),
          unexp4)
  H5 <- rbind(H5,h5)
}
H5 <- data.frame(H5,row.names = shapes)
H5$Shape <- rownames(H5)
rownames(H5) <- NULL
colnames(H5) <- c("Genotype","Treatment","Interaction","Unexplained","Shape")
H5$Shape <-  ordered(H5$Shape,levels=H5$Shape[order(H5$Unexplained)])
H5_melt <- melt(H5,id=c("Shape"))
H5_melt$variable <- ordered(H5_melt$variable,levels=c("Unexplained","Genotype","Treatment","Interaction"))
head(H5_melt)
library(scales)
p4 <- ggplot(data=H5_melt,aes(Shape,value*100))+
  geom_bar(stat = "identity",aes(fill=variable))+
  scale_fill_manual(values = c("gray60",muted("blue",l=35,c=100),"orange","purple"))+
  ylab("Variance Explained (%)")+
  xlab("Element")+
  ggtitle("Variance Explained between Control WW and Heat WL")+
  theme(plot.title = element_text(size = 20))+
  theme_bw()+
  theme(strip.background=element_rect(fill="gray50"),
        strip.text.x=element_text(size=12,color="white"),
        strip.text.y=element_text(size=12,color="white"))+
  theme(axis.text = element_text(size = 12),
        axis.title.y= element_text(size = 12),
        axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"),
        plot.margin=unit(c(0.1,0.25,0.25,0.48), "cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill=NA, size=1,linetype = 1))+
  theme(legend.position = "right")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
p4
ggsave("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/variance_explained_controlWW&HeatWL.png",width=5.5,height=5.54,plot = p4, dpi = 300)
ggsave("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/variance_explained_controlWW&HeatWL.pdf",width=7.17,height=5.54,plot = p4, dpi = 300)


# GWAS Section
# make df with only WW data from heat and control temps
heatWL_avg <- heatavg_day17[heatavg_day17$Treatment == "WL",]
heatWL_avg[heatWL_avg=="WL"]<-"heat WL"

droughtWW_avg <- avg_day15[avg_day15$Treatment == "WW",]
droughtWW_avg[droughtWW_avg=="WW"]<-"control"

control_and_heatWL_avg <- rbind(heatWL_avg, droughtWW_avg)

table(control_and_heatWL_avg$Genotype, control_and_heatWL_avg$Treatment)

height_dif4 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avg, control_and_heatWL_avg$Genotype), function(i)diff(i$avg_height)))), c("height_dif4"))
height_dif4$Genotype <- row.names(height_dif4)
head(height_dif4)
row.names(height_dif4) <- NULL

joined_pca_df4 <- join(joined_pca_df3, height_dif4, by="Genotype", type="inner")
head(joined_pca_df4)

#data <- data.table::fread("populations.structure",skip = 1)
#my_names <- as.character(read.table("populations.structure",nrows = 1,sep="\t"))
#head(data[1:5,1:10])

#rename genotype column
#colnames(data)[1] <- "Genotype"
#data[45, 1] = "Arc1"
#data[46, 1] = "Arc1"
#data = subset(data, select = -c(V2))

joined_pca_df_final <- join(joined_pca_df4, data, by="Genotype", type="inner")
write.csv(joined_pca_df_final, "./brachy_heat_final_joined_gbs_data.csv", row.names = FALSE)
joined_pca_df_final[1:5,1:20]

#climate_asso <- joined_pca_df_final[,c("Genotype", "PC1", "PC2", "PC3", "PC4", "height_dif", "height_dif2", "height_dif3")]
#head(climate_asso)
#climate_asso <- climate_asso[!duplicated(climate_asso),]
#with(climate_asso, cor.test(PC2,height_dif2,method = "spearman"))

library(lme4)
#str(joined_pca_df_final)
#joined_pca_df_final$`96_83`
joined_pca_df_final$snp <- joined_pca_df_final[,"96_83"]

mod4 <- lmer(data= joined_pca_df_final, height_dif4~0+as.factor(snp)+ (1|V3)+(1|V4)+(1|V5)+(1|PC1))

my_test4 <- Anova(mod4)
my_test4$`Pr(>Chisq)`

library(parallel)
head(joined_pca_df_final)
snp_names <- colnames(joined_pca_df_final)[15:ncol(joined_pca_df_final)]

system.time({my_gwas4 <- mclapply(snp_names, mc.cores = 4, function(i){
  tryCatch({
    joined_pca_df_final$snp <- joined_pca_df_final[,i]
    suppressWarnings(mod4 <- lmer(data= joined_pca_df_final, height_dif4~as.factor(snp)+ (1|V3)+(1|V4)+(1|V5)+(1|PC1)))
    suppressWarnings(my_test4 <- Anova(mod4))
    my_test4$`Pr(>Chisq)`
  },warning=function(war){},error=function(err){})
})
})

names(my_gwas4) <- snp_names

#copy my_gwas4 so don't have to re-run gwas
my_gwas4_copy <- my_gwas4


my_gwas4[unlist(lapply(my_gwas4, function(i) !is.null(i)))]

#quick manhattan plot to see what results look like
plot(-log10(p.adjust(unlist(my_gwas4),method="fdr")),pch=16,col="grey20",main="WW in Control and WL in Heat",xlab="SNP",ylab="-log10(P)")
abline(h=-log10(0.00001), col="orange")


# convert to dataframe
my_gwas4_copy.df <- ldply(unlist(my_gwas4_copy), data.frame)

#remove any row names and rename columns
row.names(my_gwas4_copy.df) <- NULL
colnames(my_gwas4_copy.df)[1] <- "ID"
colnames(my_gwas4_copy.df)[2] <- "Value"

### Genomic Inflation Factor Correction ###
chisq4 <- qchisq(1-my_gwas4_copy.df$Value,1)
lambda4 <- median(chisq4)/qchisq(0.5,1)
my_gwas4_copy.df$P_new <- pchisq(chisq4/lambda4, df=1,lower.tail=FALSE)


# read in file with chromosome info
#chrominfo <- read.csv("populations.plink.csv")

# make sure ID column in both data.frames is character
#chrominfo$ID <- as.character(chrominfo$ID)
#chrominfo$POS <- as.numeric(chrominfo$POS)

my_gwas4_copy.df$ID <- as.character(my_gwas4_copy.df$ID)
my_gwas4_copy.df$Value <- as.numeric(my_gwas4_copy.df$Value)


# join gwas df with chromosome info
gwas4_joined <- join(chrominfo, my_gwas4_copy.df, by="ID", type="left")

#remove centromere snp
gwas4_joined <- gwas4_joined[!(gwas4_joined$CHROM=="Bd1_centromere_containing_Bradi1g41430"),]

# Rename chromosomes to numbers
gwas4_joined[gwas4_joined=="Bd1"]<-"1"
gwas4_joined[gwas4_joined=="Bd2"]<-"2"
gwas4_joined[gwas4_joined=="Bd3"]<-"3"
gwas4_joined[gwas4_joined=="Bd4"]<-"4"
gwas4_joined[gwas4_joined=="Bd5"]<-"5"

gwas4_joined$CHROM <- as.numeric(gwas4_joined$CHROM)

unique(gwas4_joined$CHROM)

gwas4_joined$p.adj <- p.adjust(gwas4_joined$Value, method = "bonferroni")

library(qqman)
gwas4_joined = drop_na(gwas4_joined)
# Make the Manhattan plot on the gwas_joined dataset
manhattan(gwas4_joined, chr="CHROM", bp="POS", snp="ID", p="Value", annotatePval = .01 )
# Plot manhattan of genomic inflation corrected p-values
manhattan(gwas4_joined, chr="CHROM", bp="POS", snp="ID", p="P_new", annotatePval = .01 )
# Plot manhattan of bonferroni corrected genomic inflation corrected p-values
manhattan(gwas4_joined, chr="CHROM", bp="POS", snp="ID", p="p.adj", annotatePval = .01 )
# plot QQ plots of all 3 versions
qq(gwas4_joined$Value)
qq(gwas4_joined$P_new)
qq(gwas4_joined$p.adj)


don4 <- gwas4_joined %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas4_joined, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( POScum=POS+tot)

axisdf4 = don4 %>% group_by(CHROM) %>% summarize(center=(max(POScum) + min(POScum) ) / 2 )

manhattan4 <- ggplot(don4, aes(x=POScum, y=-log10(p.adj))) +
  # Show all points
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black", "grey"), 5 )) +
  # custom axis:
  scale_x_continuous(label = axisdf4$CHROM, breaks= axisdf4$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  expand_limits(y = c(0, 30))+
  xlab("Chromosome")+
  #significance line
  geom_hline(yintercept=-log10(5e-8),color='red')+
  geom_hline(yintercept=-log10(5e-6),color='blue')+
  ggtitle("GWAS Comparing Control to Heat WL")+
  # Customize the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
  )
manhattan4
ggsave(plot = manhattan,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.png", width = 6, height = 4, dpi = 300)
ggsave(plot = manhattan,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.pdf", width = 6, height = 4, dpi = 300)


################################################################
# Combine all Height GWAS results together to make one manhattan plot
################################################################

allgwas <- data.frame("ID" = unique(c(gwas_joined$ID, gwas3_joined$ID, gwas4_joined$ID)), stringsAsFactors = FALSE)
head(allgwas)

allgwas <- join(allgwas, gwas_joined, by="ID", type="left")
allgwas = select(allgwas, -P_new, -p.adj)
allgwas <- join(allgwas, gwas3_joined[,c("ID","Value")], by="ID", type="left")
allgwas <- join(allgwas, gwas4_joined[,c("ID","Value")], by="ID", type="left")

colnames(allgwas)[4:6] <- c("Control.Drought", "Control.Heat", "Control.HeatDrought")

allgwas_melt <- reshape2::melt(allgwas, id= c("ID", "CHROM", "POS"))
head(allgwas_melt)

allgwas_melt <- allgwas_melt[!is.na(allgwas_melt$CHROM),]
allgwas_melt$p.adj <- p.adjust(allgwas_melt$value, method = "bonferroni")

#Check Q-Q Plot
qq(allgwas_melt$value)

str(allgwas_melt)
sum(is.na(allgwas_melt$POS))

don_all <- allgwas_melt %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(allgwas_melt, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( POScum=POS+tot)

#don_all_filtered <- don_all[don_all$p.adj<5e-10,]
#don_all_filtered = drop_na(don_all_filtered)

#write.csv(don_all_filtered, "./brachy_combined_filtered_gwas_results.csv", row.names = FALSE)
#write.csv(don_all, "./brachy_combined_gwas_results.csv", row.names = FALSE)

axisdf_all = don_all %>% group_by(CHROM) %>% summarize(center=(max(POScum) + min(POScum) ) / 2 )


manhattan_all <- ggplot(don_all, aes(x=POScum, y=-log10(p.adj))) +
  # Show all points
  geom_point(aes(color=variable), alpha=0.8, size=1.3)+
  geom_point(data = na.omit(don_all[don_all$CHROM== 1&don_all$p.adj>5e-10,]),color="grey60")+
  geom_point(data = na.omit(don_all[don_all$CHROM== 2&don_all$p.adj>5e-10,]),color="grey80")+
  geom_point(data = na.omit(don_all[don_all$CHROM== 5&don_all$p.adj>5e-10,]),color="grey60")+
  geom_point(data = na.omit(don_all[don_all$CHROM== 4&don_all$p.adj>5e-10,]),color="grey80")+
  geom_point(data = na.omit(don_all[don_all$CHROM== 3&don_all$p.adj>5e-10,]),color="grey60")+
  #geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  #scale_color_manual(values = rep(c("black", "grey"), 5 )) +
  # custom axis:
  scale_x_continuous(label = axisdf_all$CHROM, breaks= axisdf_all$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  expand_limits(y = c(0, 65))+
  xlab("Chromosome")+
  #significance line
  geom_hline(yintercept=-log10(5e-8),color='blue')+
  geom_hline(yintercept=-log10(5e-10),color='red')+
  ggtitle("GWAS of Height across Stress Treatments")+
  # Customize the theme:
  theme_bw() +
  theme( 
    legend.position="right",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
  )
manhattan_all

ggsave(plot = manhattan_all,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/new_heightmanhattanplot_adjusted.png", width = 6, height = 4, dpi = 300)
ggsave(plot = manhattan_all,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/new_heightmanhattanplot.pdf", width = 6, height = 4, dpi = 300)

#see how many significant snps from each comparison
don_all
don_all2 <- filter(don_all, (p.adj < 5e-10))
filter(don_all2, (variable == "Control.Drought")) # 6 SNPs (9 after increased maf)
filter(don_all2, (variable == "Control.Heat")) # 35 SNPs (32 after increased maf)
filter(don_all2, (variable == "Control.HeatDrought")) # 89 SNPs (98 after increased maf)









################################################################
### Biomass GWAS ### Control WW and WL
################################################################
# GWAS Section
avgarea <- group_by(drought_127, Genotype, Treatment, imgdays) %>% summarise(avg_area = mean(cor_area))
avgarea$imgdays <- as.numeric(avgarea$imgdays)

avgarea_day15 <- aggregate(data = avgarea[avgarea$imgdays == 15,], avg_area~Genotype+Treatment, FUN= function(i)mean(i))
head(avgarea_day15)

table(avgarea_day15$Genotype, avgarea_day15$Treatment)

area_dif <- setNames(data.frame(do.call("rbind",lapply(split(avgarea_day15, avgarea_day15$Genotype), function(i)diff(i$avg_area)))), c("area_dif"))
area_dif$Genotype <- row.names(area_dif)
head(area_dif)
row.names(area_dif) <- NULL

joined_pca_df_area <- join(joined_pca_df, area_dif, by="Genotype", type="inner")
head(joined_pca_df_area)

data <- data.table::fread("populations.structure",skip = 1)
my_names <- as.character(read.table("populations.structure",nrows = 1,sep="\t"))

#rename genotype column for gbs data
colnames(data)[1] <- "Genotype"
data[45, 1] = "Arc1"
data[46, 1] = "Arc1"
data = subset(data, select = -c(V2))

# Keep only final 127 genotypes with all data)
data$compare.geno <- data$Genotype %in% list3
snp_data <- filter(data, compare.geno == TRUE)
length(unique(snp_data$Genotype))
snp_data = select(snp_data, -compare.geno)
data = select(data, -compare.geno)

joined_pca_df_area_final <- join(joined_pca_df_area, data, by="Genotype", type="inner")
#write.csv(joined_pca_df_area_final, "./brachy_heat_final_joined_gbs_data.csv", row.names = FALSE)
joined_pca_df_area_final[1:5,1:20]

library(lme4)
str(joined_pca_df_area_final)
joined_pca_df_area_final$`96_83`
joined_pca_df_area_final$snp <- joined_pca_df_final[,"96_83"]

mod5 <- lmer(data= joined_pca_df_area_final, area_dif~0+as.factor(snp)+ (1|V3)+(1|V4)+(1|V5)+(1|PC1))

my_test5 <- Anova(mod5)
my_test5$`Pr(>Chisq)`

library(parallel)
snp_names <- colnames(joined_pca_df_area_final)[14:ncol(joined_pca_df_area_final)]

system.time({my_gwasarea <- mclapply(snp_names, mc.cores = 4, function(i){
  tryCatch({
    joined_pca_df_area_final$snp <- joined_pca_df_area_final[,i]
    suppressWarnings(mod5 <- lmer(data= joined_pca_df_area_final, area_dif~as.factor(snp)+ (1|V3)+(1|V4)+(1|V5)+(1|PC1)))
    suppressWarnings(my_test5 <- Anova(mod5))
    my_test5$`Pr(>Chisq)`
  },warning=function(war){},error=function(err){})
})
})

names(my_gwasarea) <- snp_names

my_gwasarea[unlist(lapply(my_gwasarea, function(i) !is.null(i)))]

plot(-log10(p.adjust(unlist(my_gwasarea),method="bonferroni")),pch=16,col="grey20",main="WW and WL in Control",xlab="SNP",ylab="-log10(P)")
abline(h=-log10(0.00001), col="orange")

#copy my_gwas2 so don't have to re-run gwas
my_gwasarea_copy <- my_gwasarea

# convert to dataframe
my_gwasarea.df <- ldply(unlist(my_gwasarea_copy), data.frame)

#remove any row names and rename columns
row.names(my_gwasarea.df) <- NULL
colnames(my_gwasarea.df)[1] <- "ID"
colnames(my_gwasarea.df)[2] <- "Value"


### Genomic Inflation Factor Correction ###
chisq_area <- qchisq(1-my_gwasarea.df$Value,1)
lambda_area <- median(chisq_area)/qchisq(0.5,1)
my_gwasarea.df$P_new <- pchisq(chisq_area/lambda_area,df=1, lower.tail = F)


# read in file with chromosome info
#chrominfo <- read.csv("populations.plink.csv")

# make sure ID column in both data.frames is character
#chrominfo$ID <- as.character(chrominfo$ID)
#chrominfo$POS <- as.numeric(chrominfo$POS)

my_gwasarea.df$ID <- as.character(my_gwasarea.df$ID)
my_gwasarea.df$Value <- as.numeric(my_gwasarea.df$Value)


# join gwas df with chromosome info
gwas_area_joined <- join(chrominfo, my_gwasarea.df, by="ID", type="left")

#remove centromere snp
gwas_area_joined <- gwas_area_joined[!(gwas_area_joined$CHROM=="Bd1_centromere_containing_Bradi1g41430"),]

# Rename chromosomes to numbers
gwas_area_joined[gwas_area_joined=="Bd1"]<-"1"
gwas_area_joined[gwas_area_joined=="Bd2"]<-"2"
gwas_area_joined[gwas_area_joined=="Bd3"]<-"3"
gwas_area_joined[gwas_area_joined=="Bd4"]<-"4"
gwas_area_joined[gwas_area_joined=="Bd5"]<-"5"

gwas_area_joined$CHROM <- as.numeric(gwas_area_joined$CHROM)

unique(gwas_area_joined$CHROM)

gwas_area_joined$p.adj <- p.adjust(gwas_area_joined$P_new, method = "bonferroni")

library(qqman)
gwas_area_joined = drop_na(gwas_area_joined)
# Make the Manhattan plot on the gwas_joined dataset
manhattan(gwas_area_joined, chr="CHROM", bp="POS", snp="ID", p="Value", annotatePval = .01 )
# Plot manhattan of genomic inflation corrected p-values
manhattan(gwas_area_joined, chr="CHROM", bp="POS", snp="ID", p="P_new", annotatePval = .01 )
# Plot manhattan of bonferroni corrected genomic inflation corrected p-values
manhattan(gwas_area_joined, chr="CHROM", bp="POS", snp="ID", p="p.adj", annotatePval = .01 )
gwas_area_joined[-log10(gwas_area_joined$p.adj)>15,]
# plot QQ plots of all 3 versions
qq(gwas_area_joined$Value)
qq(gwas_area_joined$P_new)
qq(gwas_area_joined$p.adj)


donarea <- gwas_area_joined %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_area_joined, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( POScum=POS+tot)

axisdfarea = donarea %>% group_by(CHROM) %>% summarize(center=(max(POScum) + min(POScum) ) / 2 )

manhattanarea <- ggplot(donarea, aes(x=POScum, y=-log10(Value))) +
  # Show all points
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black", "grey"), 5 )) +
  # custom axis:
  scale_x_continuous(label = axisdfarea$CHROM, breaks= axisdfarea$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  expand_limits(y = c(0, 32))+
  xlab("Chromosome")+
  #significance line
  geom_hline(yintercept=-log10(5e-8),color='red')+
  geom_hline(yintercept=-log10(5e-6),color='blue')+
  ggtitle("GWAS Comparing Biomass Control to WL")+
  # Customize the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
  )
manhattanarea
#ggsave(plot = manhattanarea,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.png", width = 6, height = 4, dpi = 300)
#ggsave(plot = manhattanarea,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.pdf", width = 6, height = 4, dpi = 300)


################################################################
### Biomass GWAS ### Control WW and Heat WW
################################################################


#### GWAS Section

avgareaheat <- group_by(heat_127, Genotype, Treatment, imgdays) %>% summarise(avg_area = mean(cor_area))
avgareaheat$imgdays <- as.numeric(avgareaheat$imgdays)

heatavgarea_day17 <- aggregate(data = avgareaheat[avgareaheat$imgdays == 17,], avg_area~Genotype+Treatment, FUN= function(i)mean(i))
head(heatavgarea_day17)

table(heatavgarea_day17$Genotype, heatavgarea_day17$Treatment)

# make df with only WW data from heat and control temps
heatWW_avgarea <- heatavgarea_day17[heatavgarea_day17$Treatment == "WW",]
heatWW_avgarea[heatWW_avgarea=="WW"]<-"heat"

droughtWW_avgarea <- avgarea_day15[avgarea_day15$Treatment == "WW",]
droughtWW_avgarea[droughtWW_avgarea=="WW"]<-"control"

WW_avgarea <- rbind(heatWW_avgarea, droughtWW_avgarea)

#check to make sure all genotypes have data for both treatments on last day
table(WW_avgarea$Genotype, WW_avgarea$Treatment)

area_dif2 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgarea, WW_avgarea$Genotype), function(i)diff(i$avg_area)))), c("area_dif2"))
area_dif2$Genotype <- row.names(area_dif2)
head(area_dif2)
row.names(area_dif2) <- NULL

joined_pca_df_area2 <- join(joined_pca_df_area, area_dif2, by="Genotype", type="inner")
head(joined_pca_df_area2)

#data <- data.table::fread("populations.structure",skip = 1)
#my_names <- as.character(read.table("populations.structure",nrows = 1,sep="\t"))

#rename genotype column for gbs data
#colnames(data)[1] <- "Genotype"
#data[45, 1] = "Arc1"
#data[46, 1] = "Arc1"
#data = subset(data, select = -c(V2))

# Keep only final 127 genotypes with all data)
#data$compare.geno <- data$Genotype %in% list3
#snp_data <- filter(data, compare.geno == TRUE)
#length(unique(snp_data$Genotype))
#snp_data = select(snp_data, -compare.geno)
#data = select(data, -compare.geno)

joined_pca_df_area_final <- join(joined_pca_df_area2, data, by="Genotype", type="inner")
#write.csv(joined_pca_df_area_final, "./brachy_heat_final_joined_gbs_data.csv", row.names = FALSE)
joined_pca_df_area_final[1:5,1:20]

library(lme4)
str(joined_pca_df_area_final)
joined_pca_df_area_final$`96_83`
joined_pca_df_area_final$snp <- joined_pca_df_final[,"96_83"]

mod6 <- lmer(data= joined_pca_df_area_final, area_dif2~0+as.factor(snp) + (1|V3)+(1|V4)+(1|V5)+(1|PC1))

my_test6 <- Anova(mod6)
my_test6$`Pr(>Chisq)`

library(parallel)
#snp_names <- colnames(joined_pca_df_area_final)[15:ncol(joined_pca_df_area_final)]

system.time({my_gwasarea2 <- mclapply(snp_names, mc.cores = 4, function(i){
  tryCatch({
    joined_pca_df_area_final$snp <- joined_pca_df_area_final[,i]
    suppressWarnings(mod6 <- lmer(data= joined_pca_df_area_final, area_dif2~as.factor(snp)+ (1|V3)+(1|V4)+(1|V5)+(1|PC1)))
    suppressWarnings(my_test6 <- Anova(mod6))
    my_test6$`Pr(>Chisq)`
  },warning=function(war){},error=function(err){})
})
})

names(my_gwasarea2) <- snp_names

my_gwasarea2[unlist(lapply(my_gwasarea2, function(i) !is.null(i)))]

plot(-log10(p.adjust(unlist(my_gwasarea2),method="fdr")),pch=16,col="grey20",main="WW and WL in Control",xlab="SNP",ylab="-log10(P)")
abline(h=-log10(0.00001), col="orange")

#copy gwas so don't have to re-run gwas
my_gwasarea2_copy <- my_gwasarea2

# convert to dataframe
my_gwasarea2.df <- ldply(unlist(my_gwasarea2_copy), data.frame)

#remove any row names and rename columns
row.names(my_gwasarea2.df) <- NULL
colnames(my_gwasarea2.df)[1] <- "ID"
colnames(my_gwasarea2.df)[2] <- "Value"


# read in file with chromosome info
#chrominfo <- read.csv("populations.plink.csv")

# make sure ID column in both data.frames is character
chrominfo$ID <- as.character(chrominfo$ID)
chrominfo$POS <- as.numeric(chrominfo$POS)

my_gwasarea2.df$ID <- as.character(my_gwasarea2.df$ID)
my_gwasarea2.df$Value <- as.numeric(my_gwasarea2.df$Value)


# join gwas df with chromosome info
gwas_area_joined2 <- join(chrominfo, my_gwasarea2.df, by="ID", type="left")

#remove centromere snp
gwas_area_joined2 <- gwas_area_joined2[!(gwas_area_joined2$CHROM=="Bd1_centromere_containing_Bradi1g41430"),]

# Rename chromosomes to numbers
gwas_area_joined2[gwas_area_joined2=="Bd1"]<-"1"
gwas_area_joined2[gwas_area_joined2=="Bd2"]<-"2"
gwas_area_joined2[gwas_area_joined2=="Bd3"]<-"3"
gwas_area_joined2[gwas_area_joined2=="Bd4"]<-"4"
gwas_area_joined2[gwas_area_joined2=="Bd5"]<-"5"

gwas_area_joined2$CHROM <- as.numeric(gwas_area_joined2$CHROM)

unique(gwas_area_joined2$CHROM)

library(qqman)
gwas_area_joined2 = drop_na(gwas_area_joined2)
# Make the Manhattan plot on the gwas_joined dataset
manhattan(gwas_area_joined2, chr="CHROM", bp="POS", snp="ID", p="Value", annotatePval = .01 )
#ggtitle("GWAS Comparing Control to WL")
qq(gwas_area_joined2$Value)


donarea2 <- gwas_area_joined2 %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_area_joined2, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( POScum=POS+tot)

axisdfarea2 = donarea2 %>% group_by(CHROM) %>% summarize(center=(max(POScum) + min(POScum) ) / 2 )

manhattanarea2 <- ggplot(donarea2, aes(x=POScum, y=-log10(Value))) +
  # Show all points
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black", "grey"), 5 )) +
  # custom axis:
  scale_x_continuous(label = axisdfarea2$CHROM, breaks= axisdfarea2$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  expand_limits(y = c(0, 32))+
  xlab("Chromosome")+
  #significance line
  geom_hline(yintercept=-log10(5e-8),color='red')+
  geom_hline(yintercept=-log10(5e-6),color='blue')+
  ggtitle("GWAS Comparing Biomass Control to WL")+
  # Customize the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
  )
manhattanarea2
#ggsave(plot = manhattanarea2,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.png", width = 6, height = 4, dpi = 300)
#ggsave(plot = manhattanarea2,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.pdf", width = 6, height = 4, dpi = 300)


################################################################
### Biomass GWAS ### Control WW and Heat WL
################################################################



# make df with only avg WL data from heat
heatWL_avgarea <- heatavgarea_day17[heatavgarea_day17$Treatment == "WL",]
heatWL_avgarea[heatWL_avgarea=="WL"]<-"heat"

#droughtWW_avgarea <- avgarea_day15[avgarea_day15$Treatment == "WW",]
#droughtWW_avgarea[droughtWW_avgarea=="WW"]<-"control"

control_and_heatWL_avgarea <- rbind(heatWL_avgarea, droughtWW_avgarea)

#check to make sure all genotypes have data for both treatments on last day
table(control_and_heatWL_avgarea$Genotype, control_and_heatWL_avgarea$Treatment)

area_dif3 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgarea, control_and_heatWL_avgarea$Genotype), function(i)diff(i$avg_area)))), c("area_dif3"))
area_dif3$Genotype <- row.names(area_dif3)
head(area_dif3)
row.names(area_dif3) <- NULL

joined_pca_df_area3 <- join(joined_pca_df_area2, area_dif3, by="Genotype", type="inner")
head(joined_pca_df_area3)

#data <- data.table::fread("populations.structure",skip = 1)
#my_names <- as.character(read.table("populations.structure",nrows = 1,sep="\t"))

#rename genotype column for gbs data
#colnames(data)[1] <- "Genotype"
#data[45, 1] = "Arc1"
#data[46, 1] = "Arc1"
#data = subset(data, select = -c(V2))

# Keep only final 127 genotypes with all data)
#data$compare.geno <- data$Genotype %in% list3
#snp_data <- filter(data, compare.geno == TRUE)
#length(unique(snp_data$Genotype))
#snp_data = select(snp_data, -compare.geno)
#data = select(data, -compare.geno)

joined_pca_df_area_final <- join(joined_pca_df_area3, data, by="Genotype", type="inner")
#write.csv(joined_pca_df_area_final, "./brachy_heat_final_joined_gbs_data.csv", row.names = FALSE)
joined_pca_df_area_final[1:5,1:20]

library(lme4)
str(joined_pca_df_area_final)
joined_pca_df_area_final$`96_83`
joined_pca_df_area_final$snp <- joined_pca_df_final[,"96_83"]

mod7 <- lmer(data= joined_pca_df_area_final, area_dif3~0+as.factor(snp) + (1|V3)+(1|V4)+(1|V5)+(1|PC1))

my_test7 <- Anova(mod7)
my_test7$`Pr(>Chisq)`

library(parallel)
#snp_names <- colnames(joined_pca_df_area_final)[15:ncol(joined_pca_df_area_final)]

system.time({my_gwasarea3 <- mclapply(snp_names, mc.cores = 4, function(i){
  tryCatch({
    joined_pca_df_area_final$snp <- joined_pca_df_area_final[,i]
    suppressWarnings(mod7 <- lmer(data= joined_pca_df_area_final, area_dif3~as.factor(snp)+ (1|V3)+(1|V4)+(1|V5)+(1|PC1)))
    suppressWarnings(my_test7 <- Anova(mod7))
    my_test7$`Pr(>Chisq)`
  },warning=function(war){},error=function(err){})
})
})

names(my_gwasarea3) <- snp_names

my_gwasarea3[unlist(lapply(my_gwasarea3, function(i) !is.null(i)))]

plot(-log10(p.adjust(unlist(my_gwasarea3),method="fdr")),pch=16,col="grey20",main="WW and WL in Control",xlab="SNP",ylab="-log10(P)")
abline(h=-log10(0.00001), col="orange")

#copy gwas so don't have to re-run gwas
my_gwasarea3_copy <- my_gwasarea3

# convert to dataframe
my_gwasarea3.df <- ldply(unlist(my_gwasarea3_copy), data.frame)

#remove any row names and rename columns
row.names(my_gwasarea3.df) <- NULL
colnames(my_gwasarea3.df)[1] <- "ID"
colnames(my_gwasarea3.df)[2] <- "Value"


# read in file with chromosome info
#chrominfo <- read.csv("populations.plink.csv")

# make sure ID column in both data.frames is character
chrominfo$ID <- as.character(chrominfo$ID)
chrominfo$POS <- as.numeric(chrominfo$POS)

my_gwasarea3.df$ID <- as.character(my_gwasarea3.df$ID)
my_gwasarea3.df$Value <- as.numeric(my_gwasarea3.df$Value)


# join gwas df with chromosome info
gwas_area_joined3 <- join(chrominfo, my_gwasarea3.df, by="ID", type="left")

#remove centromere snp
gwas_area_joined3 <- gwas_area_joined3[!(gwas_area_joined3$CHROM=="Bd1_centromere_containing_Bradi1g41430"),]

# Rename chromosomes to numbers
gwas_area_joined3[gwas_area_joined3=="Bd1"]<-"1"
gwas_area_joined3[gwas_area_joined3=="Bd2"]<-"2"
gwas_area_joined3[gwas_area_joined3=="Bd3"]<-"3"
gwas_area_joined3[gwas_area_joined3=="Bd4"]<-"4"
gwas_area_joined3[gwas_area_joined3=="Bd5"]<-"5"

gwas_area_joined3$CHROM <- as.numeric(gwas_area_joined3$CHROM)

unique(gwas_area_joined3$CHROM)

library(qqman)
gwas_area_joined3 = drop_na(gwas_area_joined3)
# Make the Manhattan plot on the gwas_joined dataset
manhattan(gwas_area_joined3, chr="CHROM", bp="POS", snp="ID", p="Value", annotatePval = .01 )
#ggtitle("GWAS Comparing Control to Heat WL")
qq(gwas_area_joined3$Value)


donarea3 <- gwas_area_joined3 %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_area_joined3, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( POScum=POS+tot)

axisdfarea3 = donarea3 %>% group_by(CHROM) %>% summarize(center=(max(POScum) + min(POScum) ) / 2 )

manhattanarea3 <- ggplot(donarea3, aes(x=POScum, y=-log10(Value))) +
  # Show all points
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black", "grey"), 5 )) +
  # custom axis:
  scale_x_continuous(label = axisdfarea3$CHROM, breaks= axisdfarea3$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  expand_limits(y = c(0, 32))+
  xlab("Chromosome")+
  #significance line
  geom_hline(yintercept=-log10(5e-8),color='red')+
  geom_hline(yintercept=-log10(5e-6),color='blue')+
  ggtitle("GWAS Comparing Biomass Control to Heat WL")+
  # Customize the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
  )
manhattanarea3
#ggsave(plot = manhattanarea2,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.png", width = 6, height = 4, dpi = 300)
#ggsave(plot = manhattanarea2,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.pdf", width = 6, height = 4, dpi = 300)



################################################################
### Combine all Biomass GWAS results together to make one manhattan plot ###
################################################################

allareagwas <- data.frame("ID" = unique(c(gwas_area_joined$ID, gwas_area_joined2$ID, gwas_area_joined3$ID)), stringsAsFactors = FALSE)
head(allareagwas)

allareagwas <- join(allareagwas, gwas_area_joined, by="ID", type="left")
allareagwas <- join(allareagwas, gwas_area_joined2[,c("ID","Value")], by="ID", type="left")
allareagwas <- join(allareagwas, gwas_area_joined3[,c("ID","Value")], by="ID", type="left")

colnames(allareagwas)[4:6] <- c("Control.Drought", "Control.Heat", "Control.HeatDrought")

allareagwas_melt <- reshape2::melt(allareagwas, id= c("ID", "CHROM", "POS"))
head(allareagwas_melt)

allareagwas_melt <- allareagwas_melt[!is.na(allareagwas_melt$CHROM),]
allareagwas_melt$p.adj <- p.adjust(allareagwas_melt$value, method = "bonferroni")

#Check Q-Q Plot
qq(allareagwas_melt$p.adj)

str(allareagwas_melt)
sum(is.na(allareagwas_melt$PO))

all <- as.data.frame(all, xy=TRUE)

donarea_all <- allareagwas_melt %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(allareagwas_melt, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( POScum=POS+tot)

donarea_all_filtered <- donarea_all[donarea_all$p.adj<5e-1,]
donarea_all_filtered = drop_na(donarea_all_filtered)

write.csv(donarea_all_filtered, "./brachy_combined_filtered_biomass_gwas_results.csv", row.names = FALSE)
write.csv(donarea_all, "./brachy_combined_biomass_gwas_results.csv", row.names = FALSE)

axisareadf_all = donarea_all %>% group_by(CHROM) %>% summarize(center=(max(POScum) + min(POScum) ) / 2 )


manhattanarea_all <- ggplot(donarea_all, aes(x=POScum, y=-log10(p.adj))) +
  # Show all points
  geom_point(aes(color=variable), alpha=0.8, size=1.3)+
  geom_point(data = na.omit(donarea_all[donarea_all$CHROM== 1&donarea_all$p.adj>5e-10,]),color="grey60")+
  geom_point(data = na.omit(donarea_all[donarea_all$CHROM== 2&donarea_all$p.adj>5e-10,]),color="grey80")+
  geom_point(data = na.omit(donarea_all[donarea_all$CHROM== 3&donarea_all$p.adj>5e-10,]),color="grey60")+
  geom_point(data = na.omit(donarea_all[donarea_all$CHROM== 4&donarea_all$p.adj>5e-10,]),color="grey80")+
  geom_point(data = na.omit(donarea_all[donarea_all$CHROM== 5&donarea_all$p.adj>5e-10,]),color="grey60")+
  #geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  #scale_color_manual(values = rep(c("black", "grey"), 5 )) +
  # custom axis:
  scale_x_continuous(label = axisdf_all$CHROM, breaks= axisdf_all$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  expand_limits(y = c(0, 60))+
  xlab("Chromosome")+
  #significance line
  geom_hline(yintercept=-log10(5e-8),color='blue')+
  geom_hline(yintercept=-log10(5e-10),color='red')+
  ggtitle("GWAS of Biomass across Stress Treatments")+
  # Customize the theme:
  theme_bw() +
  theme( 
    legend.position="right",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
  )
manhattanarea_all

ggsave(plot = manhattanarea_all,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/areamanhattanplot.png", width = 6, height = 4, dpi = 300)
ggsave(plot = manhattanarea_all,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/areamanhattanplot.pdf", width = 6, height = 4, dpi = 300)

#see how many significant snps from each comparison
donarea_all
donarea_all2 <- filter(donarea_all, (p.adj < 5e-10))
filter(donarea_all2, (variable == "Control.Drought")) #4 SNPs
filter(donarea_all2, (variable == "Control.Heat")) # 0 SNPs
filter(donarea_all2, (variable == "Control.HeatDrought")) #3 SNPs











#####

################################################################
### make dataframes of unhealthy/healthy phenotype comparison between control and stress treatments ###
################################################################
####
#control and drought
####
avgstress <- group_by(drought_127, Genotype, Treatment, imgdays) %>% summarise(avg_unhealthy = mean(percent_unhealthy))
avgstress$imgdays <- as.numeric(avgstress$imgdays)

avgstress_day15 <- aggregate(data = avgstress[avgstress$imgdays == 15,], avg_unhealthy~Genotype+Treatment, FUN= function(i)mean(i))
head(avgstress_day15)

table(avgstress_day15$Genotype, avgstress_day15$Treatment)

#make df with only WL data from drought experiment (control temp)
droughtWL_avgstress <- avgstress_day15[avgstress_day15$Treatment == "WL",]
droughtWL_avgstress[droughtWL_avgstress=="WL"]<-"drought"

stress_dif <- setNames(data.frame(do.call("rbind",lapply(split(avgstress_day15, avgstress_day15$Genotype), function(i)diff(i$avg_unhealthy)))), c("stress_dif"))
stress_dif$Genotype <- row.names(stress_dif)
head(stress_dif)
row.names(stress_dif) <- NULL


####
#control and heat
####
heat_127$percent_unhealthy <- as.numeric(heat_127$percent_unhealthy)
rmna_heat_127 = drop_na(heat_127) 
avgstressheat <- group_by(rmna_heat_127, Genotype, Treatment, imgdays) %>% summarise(avg_unhealthy = mean(percent_unhealthy))
avgstressheat$imgdays <- as.numeric(avgstressheat$imgdays)

heatavgstress_day17 <- aggregate(data = avgstressheat[avgstressheat$imgdays == 17,], avg_unhealthy~Genotype+Treatment, FUN= function(i)mean(i))
head(heatavgstress_day17)

table(heatavgstress_day17$Genotype, heatavgstress_day17$Treatment)

# make df with only WW data from heat and control temps
heatWW_avgstress <- heatavgstress_day17[heatavgstress_day17$Treatment == "WW",]
heatWW_avgstress[heatWW_avgstress=="WW"]<-"heat"

droughtWW_avgstress <- avgstress_day15[avgstress_day15$Treatment == "WW",]
droughtWW_avgstress[droughtWW_avgstress=="WW"]<-"control"

WW_avgstress <- rbind(heatWW_avgstress, droughtWW_avgstress)

#check to make sure all genotypes have data for both treatments on last day
table(WW_avgstress$Genotype, WW_avgstress$Treatment)

stress_dif2 <- setNames(data.frame(do.call("rbind",lapply(split(WW_avgstress, WW_avgstress$Genotype), function(i)diff(i$avg_unhealthy)))), c("stress_dif2"))
stress_dif2$Genotype <- row.names(stress_dif2)
head(stress_dif2)
row.names(stress_dif2) <- NULL

####
# control and heat-drought
####
# make df with only avg WL data from heat
heatWL_avgstress <- heatavgstress_day17[heatavgstress_day17$Treatment == "WL",]
heatWL_avgstress[heatWL_avgstress=="WL"]<-"heat"

control_and_heatWL_avgstress <- rbind(heatWL_avgstress, droughtWW_avgstress)

#check to make sure all genotypes have data for both treatments on last day
table(control_and_heatWL_avgstress$Genotype, control_and_heatWL_avgstress$Treatment)

stress_dif3 <- setNames(data.frame(do.call("rbind",lapply(split(control_and_heatWL_avgstress, control_and_heatWL_avgstress$Genotype), function(i)diff(i$avg_unhealthy)))), c("stress_dif3"))
stress_dif3$Genotype <- row.names(stress_dif3)
head(stress_dif3)
row.names(stress_dif3) <- NULL


################################################################
### Read in hapmap and get in correct format ###
################################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT")
# read in hapmap file
hapmap <- read.csv("./brachy_hapmap.csv")
#rename Genotype names to match other files
colnames(hapmap)[21] <- "Adi-13"
colnames(hapmap)[22] <- "Adi-14"
colnames(hapmap)[23] <- "Adi-15"
colnames(hapmap)[24] <- "Adi-17"
colnames(hapmap)[25] <- "Adi-18"
colnames(hapmap)[26] <- "Adi-1"
colnames(hapmap)[27] <- "Adi-2"
colnames(hapmap)[28] <- "Adi-3"
colnames(hapmap)[29] <- "Adi-4"
colnames(hapmap)[30] <- "Adi-5"
colnames(hapmap)[31] <- "Adi-7"
colnames(hapmap)[32] <- "Adi-8"
colnames(hapmap)[33] <- "Adi-9"
colnames(hapmap)[34] <- "Arc1"
colnames(hapmap)[37] <- "Bd1-1"
colnames(hapmap)[38] <- "Bd21-0"
colnames(hapmap)[39] <- "Bd21-3"
colnames(hapmap)[134] <- "Bis-1"
colnames(hapmap)[136] <- "Gaz-2"
colnames(hapmap)[137] <- "Gaz-3"
colnames(hapmap)[138] <- "Gaz-4"
colnames(hapmap)[139] <- "Gaz-5"
colnames(hapmap)[140] <- "Gaz-7"
colnames(hapmap)[141] <- "Gaz-8"
colnames(hapmap)[143] <- "Kah-2"
colnames(hapmap)[145] <- "Kah-6"
colnames(hapmap)[146] <- "Koz-2"
colnames(hapmap)[154] <- "Tek-4"
colnames(hapmap)[155] <- "Tek-5"

# rename columns to remove any random symbols
colnames(hapmap)[1] <- "rs"
colnames(hapmap)[6] <- "assembly"
#fix snp column format
hapmap <- separate(hapmap, rs, c("rs1","rs2"))
hapmap <- unite(hapmap, rs, c(rs1, rs2))

#change data to "long" format to be able to sort by genotype
hapmap$rn <- seq_len(nrow(hapmap))
hapmap_melt <- reshape2::melt(hapmap, id= c("rs","alleles", "chrom", "pos","strand", "assembly","center", "protLSID","assayLSID", "panelLSID","QCcode", "rn"))
colnames(hapmap_melt)[13] <- "Genotype"

#keep only final 127 genotypes
hapmap_melt$Genotype <- as.character(hapmap_melt$Genotype)
hapmap_melt$compare.geno <- hapmap_melt$Genotype %in% list3
hapmap_melt <- filter(hapmap_melt, compare.geno == TRUE)
length(unique(hapmap_melt$Genotype))
hapmap_melt = select(hapmap_melt, -compare.geno)

#change data back to "wide" format for GAPIT input
hapmap <- reshape2::dcast(hapmap_melt, rn+rs+alleles+chrom+pos+strand+assembly+center+protLSID+assayLSID+panelLSID+QCcode ~ Genotype)
hapmap <- select(hapmap, -rn)

#remove centromere snp
hapmap <- hapmap[!(hapmap$chrom=="Bd1_centromere_containing_Bradi1g41430"),]

# Rename chromosomes to numbers
hapmap[hapmap=="BD1"]<-"1"
hapmap[hapmap=="BD2"]<-"2"
hapmap[hapmap=="BD3"]<-"3"
hapmap[hapmap=="BD4"]<-"4"
hapmap[hapmap=="BD5"]<-"5"

# save hapmap as tsv for GAPIT input
write_tsv(hapmap, "./hapmap_GAPIT.tsv")

### manually open "hapmap_GAPIT.tsv" in Excel and save as "tab-delimited text file"

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT")
#read in hapmap file
hapmap <- read.table("./hapmap_GAPIT.txt", header = FALSE)

# make file with Population Structure PCs
pop_structure <- select(gbs_pca, Genotype, V3, V4, V5, V6)
write.table(pop_structure, "./popstructurePC_GAPIT.txt")

#read in pop structure file
pop_structure <- read.table("./popstructurePC_GAPIT.txt", header = TRUE)
pop_structure$V3 <- as.numeric(pop_structure$V3)
pop_structure$V4 <- as.numeric(pop_structure$V4)
pop_structure$V5 <- as.numeric(pop_structure$V5)
pop_structure$V6 <- as.numeric(pop_structure$V6)

# make file with climate PCs
climate_cov <- select(pca_df, Genotype, PC1, PC2, PC3, PC4)
write.table(climate_cov, "./popstructurePC_GAPIT.txt")

str(climate_cov)

#read in pop structure file
pop_structure <- read.table("./popstructurePC_GAPIT.txt", header = TRUE)
pop_structure$V3 <- as.numeric(pop_structure$V3)
pop_structure$V4 <- as.numeric(pop_structure$V4)
pop_structure$V5 <- as.numeric(pop_structure$V5)
pop_structure$V6 <- as.numeric(pop_structure$V6)

#####


### First tests
################################################################
### GAPIT GWAS (Control and Drought)
################################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/drought")

# create phenotype (Y) file for control and drought
temporary <- join(height_dif, area_dif, by="Genotype")
ctldrt_phenotype <- join(temporary, stress_dif, by="Genotype")

#organize it
ctldrt_phenotype <- select(ctldrt_phenotype, Genotype, height_dif2, area_dif, stress_dif)
# rename columns
colnames(ctldrt_phenotype)[2] <- "Height"
colnames(ctldrt_phenotype)[3] <- "Area"
colnames(ctldrt_phenotype)[4] <- "Stress"
# save file for GAPIT input
write.table(ctldrt_phenotype, "./control_drought_phenotype_GAPIT.txt")
#read in phenotype file
ctldrt_phenotype <- read.table("./control_drought_phenotype_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(ctldrt_phenotype)
#ctldrt_phenotype$Genotype <- as.character(ctldrt_phenotype$Genotype)
#ctldrt_phenotype$Height <- as.numeric(ctldrt_phenotype$Height)
#ctldrt_phenotype$Area <- as.numeric(ctldrt_phenotype$Area)
#ctldrt_phenotype$Stress <- as.numeric(ctldrt_phenotype$Stress)

#check phenotype file
str(ctldrt_phenotype)
hist(ctldrt_phenotype$Height)
mean(ctldrt_phenotype$Height)
range(ctldrt_phenotype$Height)
sd(ctldrt_phenotype$Height)
which(is.na(ctldrt_phenotype$Height))

hist(ctldrt_phenotype$Area)
mean(ctldrt_phenotype$Area)
range(ctldrt_phenotype$Area)
sd(ctldrt_phenotype$Area)
which(is.na(ctldrt_phenotype$Area))

hist(ctldrt_phenotype$Stress)
mean(ctldrt_phenotype$Stress)
range(ctldrt_phenotype$Stress)
sd(ctldrt_phenotype$Stress)
which(is.na(ctldrt_phenotype$Stress))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

str(hapmap)
str(ctldrt_phenotype)
str(pop_structure)
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/drought")

#run GAPIT
ctldrt_gapit <- GAPIT(
  Y = ctldrt_phenotype,
  G = hapmap,
  PCA.total = 0,
  model = c("BLINK", "FarmCPU"),
  #model = c("MLM", "MLMM", "FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
)

################################################################
### GAPIT GWAS (Control and Heat)
################################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/heat")

# create phenotype (Y) file for control and drought
temporary <- join(height_dif3, area_dif2, by="Genotype")
ctlheat_phenotype <- join(temporary, stress_dif2, by="Genotype")

#organize it
ctlheat_phenotype <- select(ctlheat_phenotype, Genotype, height_dif3, area_dif2, stress_dif2)
# rename columns
colnames(ctlheat_phenotype)[2] <- "Height"
colnames(ctlheat_phenotype)[3] <- "Area"
colnames(ctlheat_phenotype)[4] <- "Stress"
# save file for GAPIT input
write.table(ctlheat_phenotype, "./control_heat_phenotype_GAPIT.txt")
#read in phenotype file
ctlheat_phenotype <- read.table("./control_heat_phenotype_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(ctlheat_phenotype)
#ctlheat_phenotype$Genotype <- as.character(ctlheat_phenotype$Genotype)
#ctlheat_phenotype$Height <- as.numeric(ctlheat_phenotype$Height)
#ctlheat_phenotype$Area <- as.numeric(ctlheat_phenotype$Area)
#ctlheat_phenotype$Stress <- as.numeric(ctlheat_phenotype$Stress)

#check phenotype file
hist(ctlheat_phenotype$Height)
mean(ctlheat_phenotype$Height)
range(ctlheat_phenotype$Height)
sd(ctlheat_phenotype$Height)
which(is.na(ctlheat_phenotype$Height))

hist(ctlheat_phenotype$Area)
mean(ctlheat_phenotype$Area)
range(ctlheat_phenotype$Area)
sd(ctlheat_phenotype$Area)
which(is.na(ctlheat_phenotype$Area))

hist(ctlheat_phenotype$Stress)
mean(ctlheat_phenotype$Stress)
range(ctlheat_phenotype$Stress)
sd(ctlheat_phenotype$Stress)
which(is.na(ctlheat_phenotype$Stress))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

#str(hapmap)
str(ctlheat_phenotype)
#str(pop_structure)

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/heat")

#run GAPIT
ctlheat_gapit <- GAPIT(
  Y = ctlheat_phenotype,
  G = hapmap,
  PCA.total = 0,
  model = c("BLINK", "FarmCPU"),
  #model = c("MLM", "MLMM", "FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
)


################################################################
### GAPIT GWAS (Control and Heatdrought)
################################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/heatdrought")

# create phenotype (Y) file for control and drought
temporary <- join(height_dif4, area_dif3, by="Genotype")
ctlheatdrought_phenotype <- join(temporary, stress_dif3, by="Genotype")

#organize it
ctlheatdrought_phenotype <- select(ctlheatdrought_phenotype, Genotype, height_dif4, area_dif3, stress_dif3)
# rename columns
colnames(ctlheatdrought_phenotype)[2] <- "Height"
colnames(ctlheatdrought_phenotype)[3] <- "Area"
colnames(ctlheatdrought_phenotype)[4] <- "Stress"
# save as tsv for GAPIT input
write.table(ctlheatdrought_phenotype, "./control_heatdrought_phenotype_GAPIT.txt")
#read in phenotype file
ctlheatdrought_phenotype <- read.table("./control_heatdrought_phenotype_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(ctlheatdrought_phenotype)
#ctlheatdrought_phenotype$Genotype <- as.character(ctlheatdrought_phenotype$Genotype)
#ctlheatdrought_phenotype$Height <- as.numeric(ctlheatdrought_phenotype$Height)
#ctlheatdrought_phenotype$Area <- as.numeric(ctlheatdrought_phenotype$Area)
#ctlheatdrought_phenotype$Stress <- as.numeric(ctlheatdrought_phenotype$Stress)

#check phenotype file
hist(ctlheatdrought_phenotype$Height)
mean(ctlheatdrought_phenotype$Height)
range(ctlheatdrought_phenotype$Height)
sd(ctlheatdrought_phenotype$Height)
which(is.na(ctlheatdrought_phenotype$Height))

hist(ctlheatdrought_phenotype$Area)
mean(ctlheatdrought_phenotype$Area)
range(ctlheatdrought_phenotype$Area)
sd(ctlheatdrought_phenotype$Area)
which(is.na(ctlheatdrought_phenotype$Area))

hist(ctlheatdrought_phenotype$Stress)
mean(ctlheatdrought_phenotype$Stress)
range(ctlheatdrought_phenotype$Stress)
sd(ctlheatdrought_phenotype$Stress)
which(is.na(ctlheatdrought_phenotype$Stress))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

#str(hapmap)
str(ctlheatdrought_phenotype)
#str(pop_structure)

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/heatdrought")

#run GAPIT
ctlheatdrought_gapit <- GAPIT(
  Y = ctlheatdrought_phenotype,
  G = hapmap,
  PCA.total = 0,
  model = c("BLINK", "FarmCPU"),
  #model = c("MLM", "MLMM", "FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
)

#####



# ***EARLY*** Comparison GWAS with GAPIT (organized by trait)
###############################################################
### GAPIT GWAS for Height (difference)
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/height_dif_early")

# create phenotype (Y) file for control and drought
temporary <- join(height_dif, height_dif3, by="Genotype")
height_phenotype <- join(temporary, height_dif4, by="Genotype")

#organize it
height_phenotype <- select(height_phenotype, Genotype, height_dif2, height_dif3, height_dif4)
# rename columns
colnames(height_phenotype)[2] <- "Drought"
colnames(height_phenotype)[3] <- "Heat"
colnames(height_phenotype)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(height_phenotype, "./height_phenotype_GAPIT.txt")
#read in phenotype file
height_phenotype <- read.table("./height_phenotype_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(height_phenotype)

#check phenotype file
hist(height_phenotype$Drought)
mean(height_phenotype$Drought)
range(height_phenotype$Drought)
sd(height_phenotype$Drought)
which(is.na(height_phenotype$Drought))

hist(height_phenotype$Heat)
mean(height_phenotype$Heat)
range(height_phenotype$Heat)
sd(height_phenotype$Heat)
which(is.na(height_phenotype$Heat))

hist(height_phenotype$Heat.Drought)
mean(height_phenotype$Heat.Drought)
range(height_phenotype$Heat.Drought)
sd(height_phenotype$Heat.Drought)
which(is.na(height_phenotype$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/height_dif_early")
#run GAPIT
height_gapit <- GAPIT(
  Y = height_phenotype,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
  kinship.cluster = "average",
  kinship.group = "Mean",
)

###############################################################
### GAPIT GWAS for Area (difference)
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/area_dif_early")

# create phenotype (Y) file for control and drought
temporary <- join(area_dif, area_dif2, by="Genotype")
area_phenotype <- join(temporary, area_dif3, by="Genotype")

#organize it
area_phenotype <- select(area_phenotype, Genotype, area_dif, area_dif2, area_dif3)
# rename columns
colnames(area_phenotype)[2] <- "Drought"
colnames(area_phenotype)[3] <- "Heat"
colnames(area_phenotype)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(area_phenotype, "./area_phenotype_GAPIT.txt")
#read in phenotype file
area_phenotype <- read.table("./area_phenotype_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(area_phenotype)

#check phenotype file
hist(area_phenotype$Drought)
mean(area_phenotype$Drought)
range(area_phenotype$Drought)
sd(area_phenotype$Drought)
which(is.na(area_phenotype$Drought))

hist(area_phenotype$Heat)
mean(area_phenotype$Heat)
range(area_phenotype$Heat)
sd(area_phenotype$Heat)
which(is.na(area_phenotype$Heat))

hist(area_phenotype$Heat.Drought)
mean(area_phenotype$Heat.Drought)
range(area_phenotype$Heat.Drought)
sd(area_phenotype$Heat.Drought)
which(is.na(area_phenotype$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/area_dif_early")
#run GAPIT
area_gapit <- GAPIT(
  Y = area_phenotype,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE
)


###############################################################
### GAPIT GWAS for Stress (difference)
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/stress_dif_early")

# create phenotype (Y) file for control and drought
temporary <- join(stress_dif, stress_dif2, by="Genotype")
stress_phenotype <- join(temporary, stress_dif3, by="Genotype")

#organize it
stress_phenotype <- select(stress_phenotype, Genotype, stress_dif, stress_dif2, stress_dif3)
# rename columns
colnames(stress_phenotype)[2] <- "Drought"
colnames(stress_phenotype)[3] <- "Heat"
colnames(stress_phenotype)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(stress_phenotype, "./stress_phenotype_GAPIT.txt")
#read in phenotype file
stress_phenotype <- read.table("./stress_phenotype_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(stress_phenotype)

#check phenotype file
hist(stress_phenotype$Drought)
mean(stress_phenotype$Drought)
range(stress_phenotype$Drought)
sd(stress_phenotype$Drought)
which(is.na(stress_phenotype$Drought))

hist(stress_phenotype$Heat)
mean(stress_phenotype$Heat)
range(stress_phenotype$Heat)
sd(stress_phenotype$Heat)
which(is.na(stress_phenotype$Heat))

hist(stress_phenotype$Heat.Drought)
mean(stress_phenotype$Heat.Drought)
range(stress_phenotype$Heat.Drought)
sd(stress_phenotype$Heat.Drought)
which(is.na(stress_phenotype$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/stress_dif_early")
#run GAPIT
stress_gapit <- GAPIT(
  Y = stress_phenotype,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE
)

#####


# ***EARLY*** Individual GWAS with GAPIT (organized by trait)
###############################################################
### GAPIT GWAS for Height
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/height_early")

#make df with only WL data from drought experiment (control temp)
droughtWL_avg <- avg_day15[avg_day15$Treatment == "WL",]
droughtWL_avg[droughtWL_avg=="WL"]<-"drought"

heatWL_avg[heatWL_avg=="heat WL"]<-"heat.drought"

# rename columns
colnames(droughtWW_avg)[3] <- "Control"
colnames(droughtWL_avg)[3] <- "Drought"
colnames(heatWW_avg)[3] <- "Heat"
colnames(heatWL_avg)[3] <- "Heat.Drought"

# create phenotype (Y) file for control and drought
temporary <- join(droughtWW_avg, droughtWL_avg, by="Genotype")
temporary2 <- join(temporary, heatWW_avg, by="Genotype")
height_phenotype2 <- join(temporary2, heatWL_avg, by="Genotype")

#organize it
height_phenotype2 <- select(height_phenotype2, Genotype, Control, Drought, Heat, Heat.Drought)

# save file for GAPIT input
write.table(height_phenotype2, "./height_phenotype2_GAPIT.txt")
#read in phenotype file
height_phenotype2 <- read.table("./height_phenotype2_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(height_phenotype2)

#check phenotype file
hist(height_phenotype2$Control)
mean(height_phenotype2$Control)
range(height_phenotype2$Control)
sd(height_phenotype2$Control)
which(is.na(height_phenotype2$Control))

hist(height_phenotype2$Drought)
mean(height_phenotype2$Drought)
range(height_phenotype2$Drought)
sd(height_phenotype2$Drought)
which(is.na(height_phenotype2$Drought))

hist(height_phenotype2$Heat)
mean(height_phenotype2$Heat)
range(height_phenotype2$Heat)
sd(height_phenotype2$Heat)
which(is.na(height_phenotype2$Heat))

hist(height_phenotype2$Heat.Drought)
mean(height_phenotype2$Heat.Drought)
range(height_phenotype2$Heat.Drought)
sd(height_phenotype2$Heat.Drought)
which(is.na(height_phenotype2$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/height_early")
#run GAPIT
height_gapit2 <- GAPIT(
  Y = height_phenotype2,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
)



###############################################################
### GAPIT GWAS for Area
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/area_early")

#make df with only WL data from drought experiment (control temp)
droughtWL_avgarea <- avgarea_day15[avgarea_day15$Treatment == "WL",]
droughtWL_avgarea[droughtWL_avgarea=="WL"]<-"drought"

heatWL_avgarea[heatWL_avgarea=="heat WL"]<-"heat.drought"

# rename columns
colnames(droughtWW_avgarea)[3] <- "Control"
colnames(droughtWL_avgarea)[3] <- "Drought"
colnames(heatWW_avgarea)[3] <- "Heat"
colnames(heatWL_avgarea)[3] <- "Heat.Drought"

# create phenotype (Y) file for control and drought
temporary <- join(droughtWW_avgarea, droughtWL_avgarea, by="Genotype")
temporary2 <- join(temporary, heatWW_avgarea, by="Genotype")
area_phenotype2 <- join(temporary2, heatWL_avgarea, by="Genotype")

#organize it
### CHECK COLUMN NAMES
area_phenotype2 <- select(area_phenotype2, Genotype, Control, Drought, Heat, Heat.Drought)

# save file for GAPIT input
write.table(area_phenotype2, "./area_phenotype2_GAPIT.txt")
#read in phenotype file
area_phenotype2 <- read.table("./area_phenotype2_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(area_phenotype2)

#check phenotype file
hist(area_phenotype2$Control)
mean(area_phenotype2$Control)
range(area_phenotype2$Control)
sd(area_phenotype2$Control)
which(is.na(area_phenotype2$Control))

hist(area_phenotype2$Drought)
mean(area_phenotype2$Drought)
range(area_phenotype2$Drought)
sd(area_phenotype2$Drought)
which(is.na(area_phenotype2$Drought))

hist(area_phenotype2$Heat)
mean(area_phenotype2$Heat)
range(area_phenotype2$Heat)
sd(area_phenotype2$Heat)
which(is.na(area_phenotype2$Heat))

hist(area_phenotype2$Heat.Drought)
mean(area_phenotype2$Heat.Drought)
range(area_phenotype2$Heat.Drought)
sd(area_phenotype2$Heat.Drought)
which(is.na(area_phenotype2$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/area_early")
#run GAPIT
area_gapit2 <- GAPIT(
  Y = area_phenotype2,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
)


###############################################################
### GAPIT GWAS for Stress
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/stress_early")

#make df with only WL data from drought experiment (control temp)
droughtWL_avgstress <- avgstress_day15[avgstress_day15$Treatment == "WL",]
droughtWL_avgstress[droughtWL_avgstress=="WL"]<-"drought"

heatWL_avgstress[heatWL_avgstress=="heat WL"]<-"heat.drought"

# rename columns
colnames(droughtWW_avgstress)[3] <- "Control"
colnames(droughtWL_avgstress)[3] <- "Drought"
colnames(heatWW_avgstress)[3] <- "Heat"
colnames(heatWL_avgstress)[3] <- "Heat.Drought"

# create phenotype (Y) file for control and drought
temporary <- join(droughtWW_avgstress, droughtWL_avgstress, by="Genotype")
temporary2 <- join(temporary, heatWW_avgstress, by="Genotype")
stress_phenotype2 <- join(temporary2, heatWL_avgstress, by="Genotype")

#organize it
### CHECK COLUMN NAMES 
stress_phenotype2 <- select(stress_phenotype2, Genotype, Control, Drought, Heat, Heat.Drought)

# save file for GAPIT input
write.table(stress_phenotype2, "./stress_phenotype2_GAPIT.txt")
#read in phenotype file
stress_phenotype2 <- read.table("./stress_phenotype2_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(stress_phenotype2)

#check phenotype file
hist(stress_phenotype2$Control)
mean(stress_phenotype2$Control)
range(stress_phenotype2$Control)
sd(stress_phenotype2$Control)
which(is.na(stress_phenotype2$Control))

hist(stress_phenotype2$Drought)
mean(stress_phenotype2$Drought)
range(stress_phenotype2$Drought)
sd(stress_phenotype2$Drought)
which(is.na(stress_phenotype2$Drought))

hist(stress_phenotype2$Heat)
mean(stress_phenotype2$Heat)
range(stress_phenotype2$Heat)
sd(stress_phenotype2$Heat)
which(is.na(stress_phenotype2$Heat))

hist(stress_phenotype2$Heat.Drought)
mean(stress_phenotype2$Heat.Drought)
range(stress_phenotype2$Heat.Drought)
sd(stress_phenotype2$Heat.Drought)
which(is.na(stress_phenotype2$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/stress_early")
#run GAPIT
stress_gapit2 <- GAPIT(
  Y = stress_phenotype2,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU","BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
)



#####


# ***MIDDLE*** Comparison GWAS with GAPIT (organized by trait)
###############################################################
### GAPIT GWAS for Height (difference)
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/height_dif_mid")

# create phenotype (Y) file for control and drought
temporary <- join(height_dif, height_dif3, by="Genotype")
height_phenotype <- join(temporary, height_dif4, by="Genotype")

#organize it
height_phenotype <- select(height_phenotype, Genotype, height_dif2, height_dif3, height_dif4)
# rename columns
colnames(height_phenotype)[2] <- "Drought"
colnames(height_phenotype)[3] <- "Heat"
colnames(height_phenotype)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(height_phenotype, "./height_phenotype_GAPIT.txt")
#read in phenotype file
height_phenotype <- read.table("./height_phenotype_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(height_phenotype)

#check phenotype file
hist(height_phenotype$Drought)
mean(height_phenotype$Drought)
range(height_phenotype$Drought)
sd(height_phenotype$Drought)
which(is.na(height_phenotype$Drought))

hist(height_phenotype$Heat)
mean(height_phenotype$Heat)
range(height_phenotype$Heat)
sd(height_phenotype$Heat)
which(is.na(height_phenotype$Heat))

hist(height_phenotype$Heat.Drought)
mean(height_phenotype$Heat.Drought)
range(height_phenotype$Heat.Drought)
sd(height_phenotype$Heat.Drought)
which(is.na(height_phenotype$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/height_dif_mid")
#run GAPIT
height_gapit <- GAPIT(
  Y = height_phenotype,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
  kinship.cluster = "average",
  kinship.group = "Mean",
)

###############################################################
### GAPIT GWAS for Area (difference)
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/area_dif_mid")

# create phenotype (Y) file for control and drought
temporary <- join(area_dif, area_dif2, by="Genotype")
area_phenotype <- join(temporary, area_dif3, by="Genotype")

#organize it
area_phenotype <- select(area_phenotype, Genotype, area_dif, area_dif2, area_dif3)
# rename columns
colnames(area_phenotype)[2] <- "Drought"
colnames(area_phenotype)[3] <- "Heat"
colnames(area_phenotype)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(area_phenotype, "./area_phenotype_GAPIT.txt")
#read in phenotype file
area_phenotype <- read.table("./area_phenotype_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(area_phenotype)

#check phenotype file
hist(area_phenotype$Drought)
mean(area_phenotype$Drought)
range(area_phenotype$Drought)
sd(area_phenotype$Drought)
which(is.na(area_phenotype$Drought))

hist(area_phenotype$Heat)
mean(area_phenotype$Heat)
range(area_phenotype$Heat)
sd(area_phenotype$Heat)
which(is.na(area_phenotype$Heat))

hist(area_phenotype$Heat.Drought)
mean(area_phenotype$Heat.Drought)
range(area_phenotype$Heat.Drought)
sd(area_phenotype$Heat.Drought)
which(is.na(area_phenotype$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/area_dif_mid")
#run GAPIT
area_gapit <- GAPIT(
  Y = area_phenotype,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE
)


###############################################################
### GAPIT GWAS for Stress (difference)
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/stress_dif_mid")

# create phenotype (Y) file for control and drought
temporary <- join(stress_dif, stress_dif2, by="Genotype")
stress_phenotype <- join(temporary, stress_dif3, by="Genotype")

#organize it
stress_phenotype <- select(stress_phenotype, Genotype, stress_dif, stress_dif2, stress_dif3)
# rename columns
colnames(stress_phenotype)[2] <- "Drought"
colnames(stress_phenotype)[3] <- "Heat"
colnames(stress_phenotype)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(stress_phenotype, "./stress_phenotype_GAPIT.txt")
#read in phenotype file
stress_phenotype <- read.table("./stress_phenotype_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(stress_phenotype)

#check phenotype file
hist(stress_phenotype$Drought)
mean(stress_phenotype$Drought)
range(stress_phenotype$Drought)
sd(stress_phenotype$Drought)
which(is.na(stress_phenotype$Drought))

hist(stress_phenotype$Heat)
mean(stress_phenotype$Heat)
range(stress_phenotype$Heat)
sd(stress_phenotype$Heat)
which(is.na(stress_phenotype$Heat))

hist(stress_phenotype$Heat.Drought)
mean(stress_phenotype$Heat.Drought)
range(stress_phenotype$Heat.Drought)
sd(stress_phenotype$Heat.Drought)
which(is.na(stress_phenotype$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/stress_dif_mid")
#run GAPIT
stress_gapit <- GAPIT(
  Y = stress_phenotype,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE
)

#####


# ***MIDDLE*** Individual GWAS with GAPIT (organized by trait)
###############################################################
### GAPIT GWAS for Height
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/height_mid")

#make df with only WL data from drought experiment (control temp)
droughtWL_avg <- avg_day15[avg_day15$Treatment == "WL",]
droughtWL_avg[droughtWL_avg=="WL"]<-"drought"

heatWL_avg[heatWL_avg=="heat WL"]<-"heat.drought"

# rename columns
colnames(droughtWW_avg)[3] <- "Control"
colnames(droughtWL_avg)[3] <- "Drought"
colnames(heatWW_avg)[3] <- "Heat"
colnames(heatWL_avg)[3] <- "Heat.Drought"

# create phenotype (Y) file for control and drought
temporary <- join(droughtWW_avg, droughtWL_avg, by="Genotype")
temporary2 <- join(temporary, heatWW_avg, by="Genotype")
height_phenotype2 <- join(temporary2, heatWL_avg, by="Genotype")

#organize it
height_phenotype2 <- select(height_phenotype2, Genotype, Control, Drought, Heat, Heat.Drought)

# save file for GAPIT input
write.table(height_phenotype2, "./height_phenotype2_GAPIT.txt")
#read in phenotype file
height_phenotype2 <- read.table("./height_phenotype2_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(height_phenotype2)

#check phenotype file
hist(height_phenotype2$Control)
mean(height_phenotype2$Control)
range(height_phenotype2$Control)
sd(height_phenotype2$Control)
which(is.na(height_phenotype2$Control))

hist(height_phenotype2$Drought)
mean(height_phenotype2$Drought)
range(height_phenotype2$Drought)
sd(height_phenotype2$Drought)
which(is.na(height_phenotype2$Drought))

hist(height_phenotype2$Heat)
mean(height_phenotype2$Heat)
range(height_phenotype2$Heat)
sd(height_phenotype2$Heat)
which(is.na(height_phenotype2$Heat))

hist(height_phenotype2$Heat.Drought)
mean(height_phenotype2$Heat.Drought)
range(height_phenotype2$Heat.Drought)
sd(height_phenotype2$Heat.Drought)
which(is.na(height_phenotype2$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/height_mid")
#run GAPIT
height_gapit2 <- GAPIT(
  Y = height_phenotype2,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
)



###############################################################
### GAPIT GWAS for Area
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/area_mid")

#make df with only WL data from drought experiment (control temp)
droughtWL_avgarea <- avgarea_day15[avgarea_day15$Treatment == "WL",]
droughtWL_avgarea[droughtWL_avgarea=="WL"]<-"drought"

heatWL_avgarea[heatWL_avgarea=="heat WL"]<-"heat.drought"

# rename columns
colnames(droughtWW_avgarea)[3] <- "Control"
colnames(droughtWL_avgarea)[3] <- "Drought"
colnames(heatWW_avgarea)[3] <- "Heat"
colnames(heatWL_avgarea)[3] <- "Heat.Drought"

# create phenotype (Y) file for control and drought
temporary <- join(droughtWW_avgarea, droughtWL_avgarea, by="Genotype")
temporary2 <- join(temporary, heatWW_avgarea, by="Genotype")
area_phenotype2 <- join(temporary2, heatWL_avgarea, by="Genotype")

#organize it
### CHECK COLUMN NAMES
area_phenotype2 <- select(area_phenotype2, Genotype, Control, Drought, Heat, Heat.Drought)

# save file for GAPIT input
write.table(area_phenotype2, "./area_phenotype2_GAPIT.txt")
#read in phenotype file
area_phenotype2 <- read.table("./area_phenotype2_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(area_phenotype2)

#check phenotype file
hist(area_phenotype2$Control)
mean(area_phenotype2$Control)
range(area_phenotype2$Control)
sd(area_phenotype2$Control)
which(is.na(area_phenotype2$Control))

hist(area_phenotype2$Drought)
mean(area_phenotype2$Drought)
range(area_phenotype2$Drought)
sd(area_phenotype2$Drought)
which(is.na(area_phenotype2$Drought))

hist(area_phenotype2$Heat)
mean(area_phenotype2$Heat)
range(area_phenotype2$Heat)
sd(area_phenotype2$Heat)
which(is.na(area_phenotype2$Heat))

hist(area_phenotype2$Heat.Drought)
mean(area_phenotype2$Heat.Drought)
range(area_phenotype2$Heat.Drought)
sd(area_phenotype2$Heat.Drought)
which(is.na(area_phenotype2$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/area_mid")
#run GAPIT
area_gapit2 <- GAPIT(
  Y = area_phenotype2,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
)


###############################################################
### GAPIT GWAS for Stress
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/stress_mid")

#make df with only WL data from drought experiment (control temp)
droughtWL_avgstress <- avgstress_day15[avgstress_day15$Treatment == "WL",]
droughtWL_avgstress[droughtWL_avgstress=="WL"]<-"drought"

heatWL_avgstress[heatWL_avgstress=="heat WL"]<-"heat.drought"

# rename columns
colnames(droughtWW_avgstress)[3] <- "Control"
colnames(droughtWL_avgstress)[3] <- "Drought"
colnames(heatWW_avgstress)[3] <- "Heat"
colnames(heatWL_avgstress)[3] <- "Heat.Drought"

# create phenotype (Y) file for control and drought
temporary <- join(droughtWW_avgstress, droughtWL_avgstress, by="Genotype")
temporary2 <- join(temporary, heatWW_avgstress, by="Genotype")
stress_phenotype2 <- join(temporary2, heatWL_avgstress, by="Genotype")

#organize it
### CHECK COLUMN NAMES 
stress_phenotype2 <- select(stress_phenotype2, Genotype, Control, Drought, Heat, Heat.Drought)

# save file for GAPIT input
write.table(stress_phenotype2, "./stress_phenotype2_GAPIT.txt")
#read in phenotype file
stress_phenotype2 <- read.table("./stress_phenotype2_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(stress_phenotype2)

#check phenotype file
hist(stress_phenotype2$Control)
mean(stress_phenotype2$Control)
range(stress_phenotype2$Control)
sd(stress_phenotype2$Control)
which(is.na(stress_phenotype2$Control))

hist(stress_phenotype2$Drought)
mean(stress_phenotype2$Drought)
range(stress_phenotype2$Drought)
sd(stress_phenotype2$Drought)
which(is.na(stress_phenotype2$Drought))

hist(stress_phenotype2$Heat)
mean(stress_phenotype2$Heat)
range(stress_phenotype2$Heat)
sd(stress_phenotype2$Heat)
which(is.na(stress_phenotype2$Heat))

hist(stress_phenotype2$Heat.Drought)
mean(stress_phenotype2$Heat.Drought)
range(stress_phenotype2$Heat.Drought)
sd(stress_phenotype2$Heat.Drought)
which(is.na(stress_phenotype2$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/stress_mid")
#run GAPIT
stress_gapit2 <- GAPIT(
  Y = stress_phenotype2,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU","BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
)



#####



# ***LAST DAY*** Comparison GWAS with GAPIT (organized by trait)
###############################################################
### GAPIT GWAS for Height (difference)
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/height_dif_final")

# create phenotype (Y) file for control and drought
temporary <- join(height_dif, height_dif3, by="Genotype")
height_phenotype <- join(temporary, height_dif4, by="Genotype")

#organize it
height_phenotype <- select(height_phenotype, Genotype, height_dif2, height_dif3, height_dif4)
# rename columns
colnames(height_phenotype)[2] <- "Drought"
colnames(height_phenotype)[3] <- "Heat"
colnames(height_phenotype)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(height_phenotype, "./height_phenotype_GAPIT.txt")
#read in phenotype file
height_phenotype <- read.table("./height_phenotype_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(height_phenotype)

#check phenotype file
hist(height_phenotype$Drought)
mean(height_phenotype$Drought)
range(height_phenotype$Drought)
sd(height_phenotype$Drought)
which(is.na(height_phenotype$Drought))

hist(height_phenotype$Heat)
mean(height_phenotype$Heat)
range(height_phenotype$Heat)
sd(height_phenotype$Heat)
which(is.na(height_phenotype$Heat))

hist(height_phenotype$Heat.Drought)
mean(height_phenotype$Heat.Drought)
range(height_phenotype$Heat.Drought)
sd(height_phenotype$Heat.Drought)
which(is.na(height_phenotype$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/height_dif_final")
#run GAPIT
height_gapit <- GAPIT(
  Y = height_phenotype,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
  kinship.cluster = "average",
  kinship.group = "Mean",
)

###############################################################
### GAPIT GWAS for Area (difference)
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/area_dif_final")

# create phenotype (Y) file for control and drought
temporary <- join(area_dif, area_dif2, by="Genotype")
area_phenotype <- join(temporary, area_dif3, by="Genotype")

#organize it
area_phenotype <- select(area_phenotype, Genotype, area_dif, area_dif2, area_dif3)
# rename columns
colnames(area_phenotype)[2] <- "Drought"
colnames(area_phenotype)[3] <- "Heat"
colnames(area_phenotype)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(area_phenotype, "./area_phenotype_GAPIT.txt")
#read in phenotype file
area_phenotype <- read.table("./area_phenotype_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(area_phenotype)

#check phenotype file
hist(area_phenotype$Drought)
mean(area_phenotype$Drought)
range(area_phenotype$Drought)
sd(area_phenotype$Drought)
which(is.na(area_phenotype$Drought))

hist(area_phenotype$Heat)
mean(area_phenotype$Heat)
range(area_phenotype$Heat)
sd(area_phenotype$Heat)
which(is.na(area_phenotype$Heat))

hist(area_phenotype$Heat.Drought)
mean(area_phenotype$Heat.Drought)
range(area_phenotype$Heat.Drought)
sd(area_phenotype$Heat.Drought)
which(is.na(area_phenotype$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/area_dif_final")
#run GAPIT
area_gapit <- GAPIT(
  Y = area_phenotype,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE
)


###############################################################
### GAPIT GWAS for Stress (difference)
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/stress_dif_final")

# create phenotype (Y) file for control and drought
temporary <- join(stress_dif, stress_dif2, by="Genotype")
stress_phenotype <- join(temporary, stress_dif3, by="Genotype")

#organize it
stress_phenotype <- select(stress_phenotype, Genotype, stress_dif, stress_dif2, stress_dif3)
# rename columns
colnames(stress_phenotype)[2] <- "Drought"
colnames(stress_phenotype)[3] <- "Heat"
colnames(stress_phenotype)[4] <- "Heat.Drought"
# save file for GAPIT input
write.table(stress_phenotype, "./stress_phenotype_GAPIT.txt")
#read in phenotype file
stress_phenotype <- read.table("./stress_phenotype_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(stress_phenotype)

#check phenotype file
hist(stress_phenotype$Drought)
mean(stress_phenotype$Drought)
range(stress_phenotype$Drought)
sd(stress_phenotype$Drought)
which(is.na(stress_phenotype$Drought))

hist(stress_phenotype$Heat)
mean(stress_phenotype$Heat)
range(stress_phenotype$Heat)
sd(stress_phenotype$Heat)
which(is.na(stress_phenotype$Heat))

hist(stress_phenotype$Heat.Drought)
mean(stress_phenotype$Heat.Drought)
range(stress_phenotype$Heat.Drought)
sd(stress_phenotype$Heat.Drought)
which(is.na(stress_phenotype$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/stress_dif_final")
#run GAPIT
stress_gapit <- GAPIT(
  Y = stress_phenotype,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE
)
  
#####


# ***LAST DAY*** Individual GWAS with GAPIT (organized by trait)
###############################################################
### GAPIT GWAS for Height
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/height_final")

#make df with only WL data from drought experiment (control temp)
droughtWL_avg <- avg_day15[avg_day15$Treatment == "WL",]
droughtWL_avg[droughtWL_avg=="WL"]<-"drought"

heatWL_avg[heatWL_avg=="heat WL"]<-"heat.drought"

# rename columns
colnames(droughtWW_avg)[3] <- "Control"
colnames(droughtWL_avg)[3] <- "Drought"
colnames(heatWW_avg)[3] <- "Heat"
colnames(heatWL_avg)[3] <- "Heat.Drought"

# create phenotype (Y) file for control and drought
temporary <- join(droughtWW_avg, droughtWL_avg, by="Genotype")
temporary2 <- join(temporary, heatWW_avg, by="Genotype")
height_phenotype2 <- join(temporary2, heatWL_avg, by="Genotype")

#organize it
height_phenotype2 <- select(height_phenotype2, Genotype, Control, Drought, Heat, Heat.Drought)

# save file for GAPIT input
write.table(height_phenotype2, "./height_phenotype2_GAPIT.txt")
#read in phenotype file
height_phenotype2 <- read.table("./height_phenotype2_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(height_phenotype2)

#check phenotype file
hist(height_phenotype2$Control)
mean(height_phenotype2$Control)
range(height_phenotype2$Control)
sd(height_phenotype2$Control)
which(is.na(height_phenotype2$Control))

hist(height_phenotype2$Drought)
mean(height_phenotype2$Drought)
range(height_phenotype2$Drought)
sd(height_phenotype2$Drought)
which(is.na(height_phenotype2$Drought))

hist(height_phenotype2$Heat)
mean(height_phenotype2$Heat)
range(height_phenotype2$Heat)
sd(height_phenotype2$Heat)
which(is.na(height_phenotype2$Heat))

hist(height_phenotype2$Heat.Drought)
mean(height_phenotype2$Heat.Drought)
range(height_phenotype2$Heat.Drought)
sd(height_phenotype2$Heat.Drought)
which(is.na(height_phenotype2$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/height_final")
#run GAPIT
height_gapit2 <- GAPIT(
  Y = height_phenotype2,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
)

  

###############################################################
### GAPIT GWAS for Area
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/area_final")

#make df with only WL data from drought experiment (control temp)
droughtWL_avgarea <- avgarea_day15[avgarea_day15$Treatment == "WL",]
droughtWL_avgarea[droughtWL_avgarea=="WL"]<-"drought"

heatWL_avgarea[heatWL_avgarea=="heat WL"]<-"heat.drought"

# rename columns
colnames(droughtWW_avgarea)[3] <- "Control"
colnames(droughtWL_avgarea)[3] <- "Drought"
colnames(heatWW_avgarea)[3] <- "Heat"
colnames(heatWL_avgarea)[3] <- "Heat.Drought"

# create phenotype (Y) file for control and drought
temporary <- join(droughtWW_avgarea, droughtWL_avgarea, by="Genotype")
temporary2 <- join(temporary, heatWW_avgarea, by="Genotype")
area_phenotype2 <- join(temporary2, heatWL_avgarea, by="Genotype")

#organize it
### CHECK COLUMN NAMES
area_phenotype2 <- select(area_phenotype2, Genotype, Control, Drought, Heat, Heat.Drought)

# save file for GAPIT input
write.table(area_phenotype2, "./area_phenotype2_GAPIT.txt")
#read in phenotype file
area_phenotype2 <- read.table("./area_phenotype2_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(area_phenotype2)

#check phenotype file
hist(area_phenotype2$Control)
mean(area_phenotype2$Control)
range(area_phenotype2$Control)
sd(area_phenotype2$Control)
which(is.na(area_phenotype2$Control))

hist(area_phenotype2$Drought)
mean(area_phenotype2$Drought)
range(area_phenotype2$Drought)
sd(area_phenotype2$Drought)
which(is.na(area_phenotype2$Drought))

hist(area_phenotype2$Heat)
mean(area_phenotype2$Heat)
range(area_phenotype2$Heat)
sd(area_phenotype2$Heat)
which(is.na(area_phenotype2$Heat))

hist(area_phenotype2$Heat.Drought)
mean(area_phenotype2$Heat.Drought)
range(area_phenotype2$Heat.Drought)
sd(area_phenotype2$Heat.Drought)
which(is.na(area_phenotype2$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/area_final")
#run GAPIT
area_gapit2 <- GAPIT(
  Y = area_phenotype2,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU", "BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
)


###############################################################
### GAPIT GWAS for Stress
###############################################################
setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/stress_final")

#make df with only WL data from drought experiment (control temp)
droughtWL_avgstress <- avgstress_day15[avgstress_day15$Treatment == "WL",]
droughtWL_avgstress[droughtWL_avgstress=="WL"]<-"drought"

heatWL_avgstress[heatWL_avgstress=="heat WL"]<-"heat.drought"

# rename columns
colnames(droughtWW_avgstress)[3] <- "Control"
colnames(droughtWL_avgstress)[3] <- "Drought"
colnames(heatWW_avgstress)[3] <- "Heat"
colnames(heatWL_avgstress)[3] <- "Heat.Drought"

# create phenotype (Y) file for control and drought
temporary <- join(droughtWW_avgstress, droughtWL_avgstress, by="Genotype")
temporary2 <- join(temporary, heatWW_avgstress, by="Genotype")
stress_phenotype2 <- join(temporary2, heatWL_avgstress, by="Genotype")

#organize it
### CHECK COLUMN NAMES 
stress_phenotype2 <- select(stress_phenotype2, Genotype, Control, Drought, Heat, Heat.Drought)

# save file for GAPIT input
write.table(stress_phenotype2, "./stress_phenotype2_GAPIT.txt")
#read in phenotype file
stress_phenotype2 <- read.table("./stress_phenotype2_GAPIT.txt", header = TRUE)
#make sure columns are in correct format
str(stress_phenotype2)

#check phenotype file
hist(stress_phenotype2$Control)
mean(stress_phenotype2$Control)
range(stress_phenotype2$Control)
sd(stress_phenotype2$Control)
which(is.na(stress_phenotype2$Control))

hist(stress_phenotype2$Drought)
mean(stress_phenotype2$Drought)
range(stress_phenotype2$Drought)
sd(stress_phenotype2$Drought)
which(is.na(stress_phenotype2$Drought))

hist(stress_phenotype2$Heat)
mean(stress_phenotype2$Heat)
range(stress_phenotype2$Heat)
sd(stress_phenotype2$Heat)
which(is.na(stress_phenotype2$Heat))

hist(stress_phenotype2$Heat.Drought)
mean(stress_phenotype2$Heat.Drought)
range(stress_phenotype2$Heat.Drought)
sd(stress_phenotype2$Heat.Drought)
which(is.na(stress_phenotype2$Heat.Drought))

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-GAPIT/stress_final")
#run GAPIT
stress_gapit2 <- GAPIT(
  Y = stress_phenotype2,
  G = hapmap,
  PCA.total = 0,
  model = c("FarmCPU","BLINK"),
  Multiple_analysis = TRUE,
  kinship.algorithm = "Loiselle",
)



#####




# Turkish GWAS
################################################################
### Turkish Climate PCA ###
################################################################
# Do PCA with climate data with just Turkish lines
locations_new$compare.geno <- locations_new$Genotype %in% list_turkey_110
turkey_climate <- filter(locations_new, compare.geno == TRUE)
length(unique(turkey_climate$Genotype))
turkey_climate = select(turkey_climate, -compare.geno)
locations_new = select(locations_new, -compare.geno)

turkey_climate_pca = PCA(turkey_climate[,11:29], scale.unit=TRUE, ncp=5, graph=T)

turkey_climate_pca_df = data.frame(
  "Genotype" = turkey_climate$Genotype,
  "Country" = turkey_climate$Country,
  "Latitude" = turkey_climate$Latitude,
  "Longitude" = turkey_climate$Longitude,
  "PC1" = turkey_climate_pca$ind$coord[, 1],
  "PC2" = turkey_climate_pca$ind$coord[, 2],
  "PC3" = turkey_climate_pca$ind$coord[, 3],
  "PC4" = turkey_climate_pca$ind$coord[, 4])

turkey_percentofvariance = data.frame(turkey_climate_pca$eig) 


turkey_pca <- ggplot(turkey_climate_pca_df, aes(x = PC1, y = PC2, color = Country))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  theme_minimal()+
  theme(panel.grid = element_blank(), panel.border = element_rect(fill="transparent"))+
  scale_colour_manual(values=cbPalette)+
  stat_ellipse(geom = "polygon", aes(fill = Country), alpha= 0.2, show.legend = FALSE, level = 0.95)+
  xlab("PC 1 (57.6%)")+
  ylab("PC 2 (21.9%)")+
  ggtitle("PCA of Accession Collection Location Climate")+
  geom_point()
turkey_pca
############################################################
### PCA of Turkish Group Population Structure ###
############################################################

# Read in GBS structure data
turkey_gbs_pca <- read.table("turkeypca.eigenvec", sep = " ", stringsAsFactors = FALSE)

# Quick look of what PCs 1 and 2 for GBS look like
plot(turkey_gbs_pca[,3],turkey_gbs_pca[,4])


# Rename Genotype column in GBS data
colnames(turkey_gbs_pca)[2] <- "Genotype"

# Join GBS PCA data with climate PCA data
turkey_joined_pca_df <- join(turkey_climate_pca_df, turkey_gbs_pca[,2:5], by="Genotype")

head(turkey_joined_pca_df)


# Plot population structure PCA for Turkish sub-pop
turkeypop_pca <- ggplot(turkey_joined_pca_df, aes(x = V3, y = V4, color = Country))+
  geom_hline(yintercept = 0, lty = 2)+
  geom_vline(xintercept = 0, lty = 2)+
  theme_minimal()+
  theme(panel.grid = element_blank(), panel.border = element_rect(fill="transparent"))+
  scale_colour_manual(values=cbPalette)+
  #stat_ellipse(geom = "polygon", aes(fill = Country), alpha= 0.2, show.legend = FALSE, level = 0.95)+
  xlab("PC 1")+
  ylab("PC 2")+
  ggtitle("PCA of Population Structure (within Turkish Sub-population)")+
  #geom_text(label=turkey_joined_pca_df$Genotype, nudge_x = .01, check_overlap = T)+
  geom_point()
turkeypop_pca
ggsave(plot = turkeypop_pca,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/turkeypoppca.png", width = 7.25, height = 6, dpi = 300)
ggsave(plot = turkeypop_pca,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/turkeypoppca.pdf", width = 7.25, height = 6, dpi = 300)

#ggplotly(turkeypop_pca)


################################################################
### Turkey Height GWAS ### Control WW and WL
################################################################

#Split collection_location column to have country alone in column
DFs$Collection_location <- as.character(DFs$Collection_location)
DFs <- separate(data= DFs, col=Collection_location, into= c("A", "B", "C", "Country"), fill = "left", sep = "\\, ")

# Read in txt file with list of accessions in "Turkish" population
turkey_accessions <- read.delim("pop_fst_turkey.txt", header = FALSE, sep = "\t")

list_turkey <- unique(turkey_accessions$V1)
DFs$compare.geno <- DFs$Genotype %in% list_turkey
turkey <- filter(DFs, compare.geno == TRUE)
length(unique(turkey$Genotype))
turkey = select(turkey, -compare.geno)
DFs = select(DFs, -compare.geno)

list_turkey_110 <- unique(turkey$Genotype)


# GWAS Section
turkey_avgheight <- group_by(turkey, Genotype, Treatment, imgdays) %>% summarise(avg_height = mean(cor_height_above_reference))
turkey_avgheight$imgdays <- as.numeric(turkey_avgheight$imgdays)

turkey_avg_day15 <- aggregate(data = turkey_avgheight[turkey_avgheight$imgdays == 15,], avg_height~Genotype+Treatment, FUN= function(i)mean(i))
head(turkey_avg_day15)

table(turkey_avg_day15$Genotype, turkey_avg_day15$Treatment)

turkey_height_dif <- setNames(data.frame(do.call("rbind",lapply(split(turkey_avg_day15, turkey_avg_day15$Genotype), function(i)diff(i$avg_height)))), c("height_dif"))
turkey_height_dif$Genotype <- row.names(turkey_height_dif)
head(turkey_height_dif)
row.names(height_dif) <- NULL


turkey_pca_df <- join(turkey_joined_pca_df, turkey_height_dif, by="Genotype", type="inner")
head(turkey_pca_df)

turkey_pca_df_final <- join(turkey_pca_df, data, by="Genotype", type="inner")
#write.csv(joined_pca_df_final, "./brachy_heat_final_joined_gbs_data.csv", row.names = FALSE)
turkey_pca_df_final[1:5,1:20]

library(lme4)
str(turkey_pca_df_final)
turkey_pca_df_final$`96_83`
turkey_pca_df_final$snp <- turkey_pca_df_final[,"96_83"]

turkey_mod <- lmer(data= turkey_pca_df_final, height_dif~0+as.factor(snp)+(1|V3)+(1|V4)+(1|V5)+(1|PC1))

turkey_test <- Anova(turkey_mod)
turkey_test$`Pr(>Chisq)`

library(parallel)
snp_names <- colnames(turkey_pca_df_final)[13:ncol(turkey_pca_df_final)]

system.time({turkey_gwas <- mclapply(snp_names, mc.cores = 4, function(i){
  tryCatch({
    turkey_pca_df_final$snp <- turkey_pca_df_final[,i]
    suppressWarnings(turkey_mod <- lmer(data= turkey_pca_df_final, height_dif~0+as.factor(snp)+(1|V3)+(1|V4)+(1|V5)+(1|PC1)))
    suppressWarnings(turkey_test <- Anova(turkey_mod))
    turkey_test$`Pr(>Chisq)`
    },warning=function(war){},error=function(err){})
})
})

names(turkey_gwas) <- snp_names

#copy so don't have to re-run gwas
turkey_gwas_copy <- turkey_gwas

turkey_gwas[unlist(lapply(turkey_gwas, function(i) !is.null(i)))]

plot(-log10(p.adjust(unlist(turkey_gwas),method="bonferroni")),pch=16,col="grey20",main="Turkey WW and WL in Control",xlab="SNP",ylab="-log10(P)")
abline(h=-log10(0.00001), col="orange")



# convert to dataframe
turkey_gwas.df <- ldply(unlist(turkey_gwas_copy), data.frame)

#remove any row names and rename columns
row.names(turkey_gwas.df) <- NULL
colnames(turkey_gwas.df)[1] <- "ID"
colnames(turkey_gwas.df)[2] <- "Value"

### Genomic Inflation Factor Correction ###
chisq_turkey <- qchisq(1-turkey_gwas.df$Value,1)
lambda_turkey <- median(chisq_turkey)/qchisq(0.5,1)
turkey_gwas.df$P_new <- pchisq(chisq_turkey/lambda_turkey, df=1,lower.tail=FALSE)


# read in file with chromosome info
#chrominfo <- read.csv("populations.plink.csv")

# make sure ID column in both data.frames is character
#chrominfo$ID <- as.character(chrominfo$ID)
#chrominfo$POS <- as.numeric(chrominfo$POS)

turkey_gwas.df$ID <- as.character(turkey_gwas.df$ID)
turkey_gwas.df$Value <- as.numeric(turkey_gwas.df$Value)

# join gwas df with chromosome info
turkey_gwas_joined <- join(chrominfo, turkey_gwas.df, by="ID", type="left")

#remove centromere snp
turkey_gwas_joined <- turkey_gwas_joined[!(turkey_gwas_joined$CHROM=="Bd1_centromere_containing_Bradi1g41430"),]

# Rename chromosomes to numbers
turkey_gwas_joined[turkey_gwas_joined=="Bd1"]<-"1"
turkey_gwas_joined[turkey_gwas_joined=="Bd2"]<-"2"
turkey_gwas_joined[turkey_gwas_joined=="Bd3"]<-"3"
turkey_gwas_joined[turkey_gwas_joined=="Bd4"]<-"4"
turkey_gwas_joined[turkey_gwas_joined=="Bd5"]<-"5"

turkey_gwas_joined$CHROM <- as.numeric(turkey_gwas_joined$CHROM)

unique(turkey_gwas_joined$CHROM)

turkey_gwas_joined$p.adj <- p.adjust(turkey_gwas_joined$P_new, method = "bonferroni")

library(qqman)
turkey_gwas_joined = drop_na(turkey_gwas_joined)
# Make the Manhattan plot on the gwas_joined dataset
manhattan(turkey_gwas_joined, chr="CHROM", bp="POS", snp="ID", p="Value", annotatePval = .01 )
# Plot manhattan of genomic inflation corrected p-values
manhattan(turkey_gwas_joined, chr="CHROM", bp="POS", snp="ID", p="P_new", annotatePval = .01 )
# Plot manhattan of bonferroni corrected genomic inflation corrected p-values
manhattan(turkey_gwas_joined, chr="CHROM", bp="POS", snp="ID", p="p.adj", annotatePval = .01 )
# plot QQ plots of all 3 versions
qq(turkey_gwas_joined$Value)
qq(turkey_gwas_joined$P_new)
qq(turkey_gwas_joined$p.adj)


turkey_don <- turkey_gwas_joined %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(turkey_gwas_joined, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( POScum=POS+tot)

turkey_axisdf = turkey_don %>% group_by(CHROM) %>% summarize(center=(max(POScum) + min(POScum) ) / 2 )

turkey_manhattan <- ggplot(turkey_don, aes(x=POScum, y=-log10(P_new))) +
  # Show all points
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black", "grey"), 5 )) +
  # custom axis:
  scale_x_continuous(label = turkey_axisdf$CHROM, breaks= turkey_axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  expand_limits(y = c(0, 15))+
  xlab("Chromosome")+
  #significance line
  geom_hline(yintercept=-log10(5e-8),color='red')+
  geom_hline(yintercept=-log10(5e-6),color='blue')+
  ggtitle("GWAS Comparing Control to WL (only Turkish Accessions)")+
  # Customize the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
  )
turkey_manhattan
ggsave(plot = manhattan,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.png", width = 6, height = 4, dpi = 300)
ggsave(plot = manhattan,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.pdf", width = 6, height = 4, dpi = 300)




################################################################
### Turkey Height GWAS ### Heat WW and Control
################################################################
#Split collection_location column to have country alone in column
heatDFs$Collection_location <- as.character(heatDFs$Collection_location)
heatDFs <- separate(data= heatDFs, col=Collection_location, into= c("A", "B", "C", "Country"), fill = "left", sep = "\\, ")

# Read in txt file with list of accessions in "Turkish" population
#turkey_accessions <- read.delim("pop_fst_turkey.txt", header = FALSE, sep = "\t")

#list_turkey <- unique(turkey_accessions$V1)
heatDFs$compare.geno <- heatDFs$Genotype %in% list_turkey
turkey_heat <- filter(heatDFs, compare.geno == TRUE)
length(unique(turkey_heat$Genotype))
turkey_heat = select(turkey_heat, -compare.geno)
heatDFs = select(heatDFs, -compare.geno)


turkey_heat_avgheight <- group_by(turkey_heat, Genotype, Treatment, imgdays) %>% summarise(avg_height = mean(cor_height_above_reference))
turkey_heat_avgheight$imgdays <- as.numeric(turkey_heat_avgheight$imgdays)

turkey_heat_avg_day17 <- aggregate(data = turkey_heat_avgheight[turkey_heat_avgheight$imgdays == 17,], avg_height~Genotype+Treatment, FUN= function(i)mean(i))
head(turkey_heat_avg_day17)

table(turkey_heat_avg_day17$Genotype, turkey_heat_avg_day17$Treatment)

# make df with only WW data from heat and control temps
turkey_heatWW_avg <- turkey_heat_avg_day17[turkey_heat_avg_day17$Treatment == "WW",]
turkey_heatWW_avg[turkey_heatWW_avg=="WW"]<-"heat"

turkey_control_avg <- turkey_avg_day15[turkey_avg_day15$Treatment == "WW",]
turkey_control_avg[turkey_control_avg=="WW"]<-"control"

turkey_WW <- rbind(turkey_heatWW_avg, turkey_control_avg)

table(turkey_WW$Genotype, turkey_WW$Treatment)

# GWAS Section
turkey_height_dif2 <- setNames(data.frame(do.call("rbind",lapply(split(turkey_WW, turkey_WW$Genotype), function(i)diff(i$avg_height)))), c("height_dif2"))
turkey_height_dif2$Genotype <- row.names(turkey_height_dif2)
head(turkey_height_dif2)
row.names(turkey_height_dif2) <- NULL

turkey_pca_df2 <- join(turkey_pca_df, turkey_height_dif2, by="Genotype", type="inner")
head(turkey_pca_df2)


turkey_pca_df_final <- join(turkey_pca_df2, data, by="Genotype", type="inner")
#write.csv(joined_pca_df_final, "./brachy_heat_final_joined_gbs_data.csv", row.names = FALSE)
turkey_pca_df_final[1:5,1:20]

library(lme4)
str(turkey_pca_df_final)
turkey_pca_df_final$`96_83`
turkey_pca_df_final$snp <- turkey_pca_df_final[,"96_83"]


turkey_mod2 <- lmer(data= turkey_pca_df_final, height_dif2~0+as.factor(snp)+(1|V3)+(1|V4)+(1|V5)+(1|PC1))

turkey_test2 <- Anova(turkey_mod2)
turkey_test2$`Pr(>Chisq)`

library(parallel)
snp_names <- colnames(turkey_pca_df_final)[14:ncol(turkey_pca_df_final)]

system.time({turkey_gwas2 <- mclapply(snp_names, mc.cores = 4, function(i){
  tryCatch({
    turkey_pca_df_final$snp <- turkey_pca_df_final[,i]
    suppressWarnings(turkey_mod2 <- lmer(data= turkey_pca_df_final, height_dif2~as.factor(snp)+(1|V3)+(1|V4)+(1|V5)+(1|PC1)))
    suppressWarnings(turkey_test2 <- Anova(turkey_mod2))
    turkey_test2$`Pr(>Chisq)`
  },warning=function(war){},error=function(err){})
})
})

names(turkey_gwas2) <- snp_names

#copy gwas data so don't have to re-run gwas
turkey_gwas2_copy <- turkey_gwas2


turkey_gwas2[unlist(lapply(turkey_gwas2, function(i) !is.null(i)))]
plot(-log10(p.adjust(unlist(turkey_gwas2),method="bonferroni")),pch=16,col="grey20",main="Turkey Control and Heat WW",xlab="SNP",ylab="-log10(P)")
abline(h=-log10(0.00001), col="orange")

# convert to dataframe
turkey_gwas2.df <- ldply(unlist(turkey_gwas2_copy), data.frame)

#remove any row names and rename columns
row.names(turkey_gwas2.df) <- NULL
colnames(turkey_gwas2.df)[1] <- "ID"
colnames(turkey_gwas2.df)[2] <- "Value"

### Genomic Inflation Factor Correction ###
chisq_turkey2 <- qchisq(1-turkey_gwas2.df$Value,1)
lambda_turkey2 <- median(chisq_turkey2)/qchisq(0.5,1)
turkey_gwas2.df$P_new <- pchisq(chisq_turkey2/lambda_turkey2, df=1,lower.tail=FALSE)

# read in file with chromosome info
#chrominfo <- read.csv("populations.plink.csv")

# make sure ID column in both data.frames is character
#chrominfo$ID <- as.character(chrominfo$ID)
#chrominfo$POS <- as.numeric(chrominfo$POS)

turkey_gwas2.df$ID <- as.character(turkey_gwas2.df$ID)
turkey_gwas2.df$Value <- as.numeric(turkey_gwas2.df$Value)

# join gwas df with chromosome info
turkey_gwas2_joined <- join(chrominfo, turkey_gwas2.df, by="ID", type="left")

#remove centromere snp
turkey_gwas2_joined <- turkey_gwas2_joined[!(turkey_gwas2_joined$CHROM=="Bd1_centromere_containing_Bradi1g41430"),]

# Rename chromosomes to numbers
turkey_gwas2_joined[turkey_gwas2_joined=="Bd1"]<-"1"
turkey_gwas2_joined[turkey_gwas2_joined=="Bd2"]<-"2"
turkey_gwas2_joined[turkey_gwas2_joined=="Bd3"]<-"3"
turkey_gwas2_joined[turkey_gwas2_joined=="Bd4"]<-"4"
turkey_gwas2_joined[turkey_gwas2_joined=="Bd5"]<-"5"

turkey_gwas2_joined$CHROM <- as.numeric(turkey_gwas2_joined$CHROM)

unique(turkey_gwas2_joined$CHROM)

turkey_gwas2_joined$p.adj <- p.adjust(turkey_gwas2_joined$P_new, method = "bonferroni")

library(qqman)
turkey_gwas2_joined = drop_na(turkey_gwas2_joined)
# Make the Manhattan plot on the gwas_joined dataset
manhattan(turkey_gwas2_joined, chr="CHROM", bp="POS", snp="ID", p="Value", annotatePval = .01 )
# Plot manhattan of genomic inflation corrected p-values
manhattan(turkey_gwas2_joined, chr="CHROM", bp="POS", snp="ID", p="P_new", annotatePval = .01 )
# Plot manhattan of bonferroni corrected genomic inflation corrected p-values
manhattan(turkey_gwas2_joined, chr="CHROM", bp="POS", snp="ID", p="p.adj", annotatePval = .01 )
# plot QQ plots of all 3 versions
qq(turkey_gwas2_joined$Value)
qq(turkey_gwas2_joined$P_new)
qq(turkey_gwas2_joined$p.adj)


turkey_don2 <- turkey_gwas2_joined %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(turkey_gwas2_joined, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( POScum=POS+tot)

turkey_axisdf2 = turkey_don2 %>% group_by(CHROM) %>% summarize(center=(max(POScum) + min(POScum) ) / 2 )

turkey_manhattan2 <- ggplot(turkey_don2, aes(x=POScum, y=-log10(Value))) +
  # Show all points
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black", "grey"), 5 )) +
  # custom axis:
  scale_x_continuous(label = turkey_axisdf2$CHROM, breaks= turkey_axisdf2$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  expand_limits(y = c(0, 100))+
  xlab("Chromosome")+
  #significance line
  geom_hline(yintercept=-log10(5e-8),color='red')+
  geom_hline(yintercept=-log10(5e-6),color='blue')+
  ggtitle("GWAS Comparing Control to Heat (only Turkish Accessions)")+
  # Customize the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
  )
turkey_manhattan2
ggsave(plot = manhattan,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.png", width = 6, height = 4, dpi = 300)
ggsave(plot = manhattan,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.pdf", width = 6, height = 4, dpi = 300)


################################################################
### Turkey Height GWAS ### Heat WL and Control
################################################################
# make df with only WL data from heat
turkey_heatWL_avg <- turkey_heat_avg_day17[turkey_heat_avg_day17$Treatment == "WL",]
turkey_heatWL_avg[turkey_heatWL_avg=="WL"]<-"heatdrought"

#combine WL heat with control data and check that there is data for every accession
turkey_WLctl <- rbind(turkey_heatWL_avg, turkey_control_avg)
table(turkey_WW$Genotype, turkey_WW$Treatment)

# GWAS Section
turkey_height_dif3 <- setNames(data.frame(do.call("rbind",lapply(split(turkey_WLctl, turkey_WLctl$Genotype), function(i)diff(i$avg_height)))), c("height_dif3"))
turkey_height_dif3$Genotype <- row.names(turkey_height_dif3)
head(turkey_height_dif3)
row.names(turkey_height_dif3) <- NULL

turkey_pca_df3 <- join(turkey_pca_df2,turkey_height_dif3, by="Genotype", type="inner")
head(turkey_pca_df3)

turkey_pca_df_final <- join(turkey_pca_df3, data, by="Genotype", type="inner")
#write.csv(joined_pca_df_final, "./brachy_heat_final_joined_gbs_data.csv", row.names = FALSE)
turkey_pca_df_final[1:5,1:20]

library(lme4)
str(turkey_pca_df_final)
turkey_pca_df_final$`96_83`
turkey_pca_df_final$snp <- turkey_pca_df_final[,"96_83"]

turkey_mod3 <- lmer(data= turkey_pca_df_final, height_dif3~0+as.factor(snp)+(1|V3)+(1|V4)+(1|V5)+(1|PC1))

turkey_test3 <- Anova(turkey_mod3)
turkey_test3$`Pr(>Chisq)`

library(parallel)
snp_names <- colnames(turkey_pca_df_final)[15:ncol(turkey_pca_df_final)]

system.time({turkey_gwas3 <- mclapply(snp_names, mc.cores = 4, function(i){
  tryCatch({
    turkey_pca_df_final$snp <- turkey_pca_df_final[,i]
    suppressWarnings(turkey_mod3 <- lmer(data= turkey_pca_df_final, height_dif3~as.factor(snp)+(1|V3)+(1|V4)+(1|V5)+(1|PC1)))
    suppressWarnings(turkey_test3 <- Anova(turkey_mod3))
    turkey_test3$`Pr(>Chisq)`
  },warning=function(war){},error=function(err){})
})
})

names(turkey_gwas3) <- snp_names

#copy gwas data so don't have to re-run gwas
turkey_gwas3_copy <- turkey_gwas3

turkey_gwas3[unlist(lapply(turkey_gwas3, function(i) !is.null(i)))]

plot(-log10(p.adjust(unlist(turkey_gwas3),method="bonferroni")),pch=16,col="grey20",main="Turkey Control and Heat WL",xlab="SNP",ylab="-log10(P)")
abline(h=-log10(0.00001), col="orange")

# convert to dataframe
turkey_gwas3.df <- ldply(unlist(turkey_gwas3_copy), data.frame)

#remove any row names and rename columns
row.names(turkey_gwas3.df) <- NULL
colnames(turkey_gwas3.df)[1] <- "ID"
colnames(turkey_gwas3.df)[2] <- "Value"

### Genomic Inflation Factor Correction ###
chisq_turkey3 <- qchisq(1-turkey_gwas3.df$Value,1)
lambda_turkey3 <- median(chisq_turkey3)/qchisq(0.5,1)
turkey_gwas3.df$P_new <- pchisq(chisq_turkey3/lambda_turkey3, df=1,lower.tail=FALSE)

# read in file with chromosome info
#chrominfo <- read.csv("populations.plink.csv")

# make sure ID column in both data.frames is character
#chrominfo$ID <- as.character(chrominfo$ID)
#chrominfo$POS <- as.numeric(chrominfo$POS)

turkey_gwas3.df$ID <- as.character(turkey_gwas3.df$ID)
turkey_gwas3.df$Value <- as.numeric(turkey_gwas3.df$Value)

# join gwas df with chromosome info
turkey_gwas3_joined <- join(chrominfo, turkey_gwas3.df, by="ID", type="left")

#remove centromere snp
turkey_gwas3_joined <- turkey_gwas3_joined[!(turkey_gwas3_joined$CHROM=="Bd1_centromere_containing_Bradi1g41430"),]

# Rename chromosomes to numbers
turkey_gwas3_joined[turkey_gwas3_joined=="Bd1"]<-"1"
turkey_gwas3_joined[turkey_gwas3_joined=="Bd2"]<-"2"
turkey_gwas3_joined[turkey_gwas3_joined=="Bd3"]<-"3"
turkey_gwas3_joined[turkey_gwas3_joined=="Bd4"]<-"4"
turkey_gwas3_joined[turkey_gwas3_joined=="Bd5"]<-"5"

turkey_gwas3_joined$CHROM <- as.numeric(turkey_gwas3_joined$CHROM)

unique(turkey_gwas3_joined$CHROM)

turkey_gwas3_joined$p.adj <- p.adjust(turkey_gwas3_joined$Value, method = "bonferroni")

library(qqman)
turkey_gwas3_joined = drop_na(turkey_gwas3_joined)
# Make the Manhattan plot on the gwas_joined dataset
manhattan(turkey_gwas3_joined, chr="CHROM", bp="POS", snp="ID", p="Value", annotatePval = .01 )
# Plot manhattan of genomic inflation corrected p-values
manhattan(turkey_gwas3_joined, chr="CHROM", bp="POS", snp="ID", p="P_new", annotatePval = .01 )
# Plot manhattan of bonferroni corrected genomic inflation corrected p-values
manhattan(turkey_gwas3_joined, chr="CHROM", bp="POS", snp="ID", p="p.adj", annotatePval = .01 )
# plot QQ plots of all 3 versions
qq(turkey_gwas3_joined$Value)
qq(turkey_gwas3_joined$P_new)
qq(turkey_gwas3_joined$p.adj)


turkey_don3 <- turkey_gwas3_joined %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(turkey_gwas3_joined, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( POScum=POS+tot)

turkey_axisdf3 = turkey_don3 %>% group_by(CHROM) %>% summarize(center=(max(POScum) + min(POScum) ) / 2 )

turkey_manhattan3 <- ggplot(turkey_don3, aes(x=POScum, y=-log10(p.adj))) +
  # Show all points
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black", "grey"), 5 )) +
  # custom axis:
  scale_x_continuous(label = turkey_axisdf3$CHROM, breaks= turkey_axisdf3$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  expand_limits(y = c(0, 25))+
  xlab("Chromosome")+
  #significance line
  geom_hline(yintercept=-log10(5e-8),color='red')+
  geom_hline(yintercept=-log10(5e-6),color='blue')+
  ggtitle("GWAS Comparing Control to WL only Turkish Accessions")+
  # Customize the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
  )
turkey_manhattan3
ggsave(plot = manhattan,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.png", width = 6, height = 4, dpi = 300)
ggsave(plot = manhattan,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/manhattanplot.pdf", width = 6, height = 4, dpi = 300)

################################################################
### Combine all Turkish Height GWAS ###
################################################################
allturkeygwas <- data.frame("ID" = unique(c(turkey_gwas_joined$ID, turkey_gwas2_joined$ID, turkey_gwas3_joined$ID)), stringsAsFactors = FALSE)
head(allturkeygwas)

allturkeygwas <- join(allturkeygwas, turkey_gwas_joined, by="ID", type="left")
allturkeygwas = select(allturkeygwas, -Value, -p.adj)
allturkeygwas <- join(allturkeygwas, turkey_gwas2_joined[,c("ID","P_new")], by="ID", type="left")
allturkeygwas <- join(allturkeygwas, turkey_gwas3_joined[,c("ID","P_new")], by="ID", type="left")

colnames(allturkeygwas)[4:6] <- c("Control.Drought", "Control.Heat", "Control.HeatDrought")

allturkeygwas_melt <- reshape2::melt(allturkeygwas, id= c("ID", "CHROM", "POS"))
head(allturkeygwas_melt)

allturkeygwas_melt <- allturkeygwas_melt[!is.na(allturkeygwas_melt$CHROM),]
allturkeygwas_melt$p.adj <- p.adjust(allturkeygwas_melt$value, method = "bonferroni")

#Check Q-Q Plot
qq(allturkeygwas_melt$value)

str(allturkeygwas_melt)
sum(is.na(allturkeygwas_melt$POS))

turkey_don_all <- allturkeygwas_melt %>% 
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(allturkeygwas_melt, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, POS) %>%
  mutate( POScum=POS+tot)

#turkey_don_all_filtered <- turkey_don_all[don_all$p.adj<5e-10,]
#turkey_don_all_filtered = drop_na(turkey_don_all_filtered)

#write.csv(don_all_filtered, "./brachy_combined_filtered_gwas_results.csv", row.names = FALSE)
#write.csv(don_all, "./brachy_combined_gwas_results.csv", row.names = FALSE)

turkey_axisdf_all = turkey_don_all %>% group_by(CHROM) %>% summarize(center=(max(POScum) + min(POScum) ) / 2 )


manhattan_all_turkey <- ggplot(turkey_don_all, aes(x=POScum, y=-log10(value))) +
  # Show all points
  geom_point(aes(color=variable), alpha=0.8, size=1.3)+
  geom_point(data = na.omit(turkey_don_all[turkey_don_all$CHROM== 1&turkey_don_all$value>5e-10,]),color="grey60")+
  geom_point(data = na.omit(turkey_don_all[turkey_don_all$CHROM== 2&turkey_don_all$value>5e-10,]),color="grey80")+
  geom_point(data = na.omit(turkey_don_all[turkey_don_all$CHROM== 3&turkey_don_all$value>5e-10,]),color="grey60")+
  geom_point(data = na.omit(turkey_don_all[turkey_don_all$CHROM== 4&turkey_don_all$value>5e-10,]),color="grey80")+
  geom_point(data = na.omit(turkey_don_all[turkey_don_all$CHROM== 5&turkey_don_all$value>5e-10,]),color="grey60")+
  #geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  #scale_color_manual(values = rep(c("black", "grey"), 5 )) +
  # custom axis:
  scale_x_continuous(label = turkey_axisdf_all$CHROM, breaks= turkey_axisdf_all$center ) +
  scale_y_continuous(expand = c(0, 0) ) +  # remove space between plot area and x axis
  expand_limits(y = c(0, 15))+
  xlab("Chromosome")+
  #significance line
  geom_hline(yintercept=-log10(5e-8),color='blue')+
  geom_hline(yintercept=-log10(5e-10),color='red')+
  ggtitle("GWAS of Height across Stress Treatments (only Turkish Accessions)")+
  # Customize the theme:
  theme_bw() +
  theme( 
    legend.position="right",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x.bottom = element_line(),
    axis.line.y.left = element_line(),
  )
manhattan_all_turkey

ggsave(plot = manhattan_all_turkey,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/heightmanhattanplot_turkey_cor.png", width = 6, height = 4, dpi = 300)
ggsave(plot = manhattan_all_turkey,"/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-figures/heightmanhattanplot_turkey_cor.pdf", width = 6, height = 4, dpi = 300)

#see how many significant snps from each comparison
turkey_don_all
turkey_don_all2 <- filter(turkey_don_all, (value < 5e-10))
filter(turkey_don_all2, (variable == "Control.Drought")) # 17 SNPs (2 after correction)
filter(turkey_don_all2, (variable == "Control.Heat")) # 110 SNPs (27 after correction)
filter(turkey_don_all2, (variable == "Control.HeatDrought")) # 335 SNPs (117 after correction)





