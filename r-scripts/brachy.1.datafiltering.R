rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.


library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
#library(readr)


#####################################################################
# read in the control/drought single-value data #
#####################################################################

#setwd("/Users/ellaludwig/Desktop/Danforth-Center/brachy/brachy-drought/")
setwd("/Users/eludwig/Desktop/brachy/")

#read in single value trait csv files
z500s <- read.csv("drought_SV_z500-single-value-traits.csv")
z1500s <- read.csv("drought_SV_z1500-single-value-traits.csv")
z2500s <- read.csv("drought_SV_z2500-single-value-traits.csv")
z3500s <- read.csv("drought_SV_z3500-single-value-traits.csv")

z500s = subset(z500s, select = -c(percent_unhealthy))
z1500s = subset(z1500s, select = -c(percent_unhealthy))
z2500s = subset(z2500s, select = -c(percent_unhealthy))
z3500s = subset(z3500s, select = -c(percent_unhealthy))

#remove duplicates of all rows due to re-running plantcv and rows being added to the end of original file instead of file being overwritten
z500s <- z500s %>% distinct(plantbarcode, timestamp, frame, imgtype, .keep_all = TRUE)
z1500s <- z1500s %>% distinct(plantbarcode, timestamp, frame, imgtype, .keep_all = TRUE)
z2500s <- z2500s %>% distinct(plantbarcode, timestamp, frame, imgtype, .keep_all = TRUE)
z3500s <- z3500s %>% distinct(plantbarcode, timestamp, frame, imgtype, .keep_all = TRUE)

#subset VIS data
visz500s = z500s[z500s$imgtype == 'VIS',]
visz1500s = z1500s[z1500s$imgtype == 'VIS',]
visz2500s = z2500s[z2500s$imgtype == 'VIS',]
visz3500s = z3500s[z3500s$imgtype == 'VIS',]

setwd("/Users/eludwig/Desktop/brachy/")
#read in single value trait csv files
z500s_damage <- read.csv("drought_SV_z500_damage-single-value-traits.csv")
z1500s_damage <- read.csv("drought_SV_z1500_damage-single-value-traits.csv")
z2500s_damage <- read.csv("drought_SV_z2500_damage-single-value-traits.csv")
z3500s_damage <- read.csv("drought_SV_z3500_damage-single-value-traits.csv")

#subset VIS data
visz500s_damage = z500s_damage[z500s_damage$imgtype == 'VIS',]
visz1500s_damage = z1500s_damage[z1500s_damage$imgtype == 'VIS',]
visz2500s_damage = z2500s_damage[z2500s_damage$imgtype == 'VIS',]
visz3500s_damage = z3500s_damage[z3500s_damage$imgtype == 'VIS',]

visz500s <- inner_join(visz500s, visz500s_damage, by = c("camera", "imgtype", "zoom", "exposure", "gain", "frame", "lifter", 
                                                            "timestamp", "id", "plantbarcode", "treatment", "cartag", "measurementlabel", "other", "image"))
visz1500s <- inner_join(visz1500s, visz1500s_damage, by = c("camera", "imgtype", "zoom", "exposure", "gain", "frame", "lifter", 
                   "timestamp", "id", "plantbarcode", "treatment", "cartag", "measurementlabel", "other", "image"))
visz2500s <- inner_join(visz2500s, visz2500s_damage, by = c("camera", "imgtype", "zoom", "exposure", "gain", "frame", "lifter", 
                                                            "timestamp", "id", "plantbarcode", "treatment", "cartag", "measurementlabel", "other", "image"))
visz3500s <- inner_join(visz3500s, visz3500s_damage, by = c("camera", "imgtype", "zoom", "exposure", "gain", "frame", "lifter", 
                                                            "timestamp", "id", "plantbarcode", "treatment", "cartag", "measurementlabel", "other", "image"))

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
visz500s <- transform(visz500s, cor_unhealthy_area = unhealthy_area / 1229.355)
visz500s <- transform(visz500s, cor_total_plant_area = total_plant_area / 1229.355)

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
visz1500s <- transform(visz1500s, cor_unhealthy_area = unhealthy_area / 2610.580)
visz1500s <- transform(visz1500s, cor_total_plant_area = total_plant_area / 2610.580)

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
visz2500s <- transform(visz2500s, cor_unhealthy_area = unhealthy_area / 5543.662)
visz2500s <- transform(visz2500s, cor_total_plant_area = total_plant_area / 5543.662)

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
visz3500s <- transform(visz3500s, cor_unhealthy_area = unhealthy_area / 11772.169)
visz3500s <- transform(visz3500s, cor_total_plant_area = total_plant_area / 11772.169)

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

DFs <- rbind(visz3500s, visz2500s, visz1500s, visz500s)

#filter data
DFs <- select(DFs, image, plantbarcode, timestamp, hue_circular_mean, hue_circular_std, hue_median, cor_area, cor_width, cor_height, 
              cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, 
              cor_area_below_reference, percent_area_below_reference, cor_convex_hull_area, solidity, cor_perimeter, cor_longest_path, 
              convex_hull_vertices, cor_ellipse_major_axis, cor_ellipse_minor_axis, ellipse_angle, ellipse_eccentricity, cor_unhealthy_area, 
              cor_total_plant_area)
DFs <- filter(DFs, plantbarcode != "Cc00AA000000") #remove any images of the color card

#separate date and time
TS <- c("date","time") #make a vector for the column names of each part of the timestamp
DFs <- separate(DFs, col = timestamp, into = c("date","time"), sep = " ") #separate each component of timestamp

#read in barcode csv file with lane information to calculate DAP (since split into two groups for planting and imaging)
barcodes_even <- read.csv("/Users/eludwig/Desktop/brachy/B002dr_032515_barcodes_statefile_order_even.csv", header= FALSE)
colnames(barcodes_even)[1] <- "Barcode"

even <- (unique(barcodes_even$Barcode))
DFs$compare.barcode <- DFs$plantbarcode %in% even
even <- filter(DFs, compare.barcode == TRUE)
even$lane_group <- "even"
even <- select(even, -compare.barcode)

barcodes_odd <- read.csv("/Users/eludwig/Desktop/brachy/B002dr_032515_barcodes_statefile_order_odd.csv", header= FALSE)
colnames(barcodes_odd)[1] <- "Barcode"

odd <- (unique(barcodes_odd$Barcode))
DFs$compare.barcode2 <- DFs$plantbarcode %in% odd
odd <- filter(DFs, compare.barcode2 == TRUE)
odd$lane_group <- "odd"
odd <- select(odd, -compare.barcode, -compare.barcode2)

#check that there isn't overlap between odd and even groups (overlap should be 0)
test <- as.data.frame(unique(even$plantbarcode))
colnames(test)[1] <- "Barcode"
test2 <- as.data.frame(unique(odd$plantbarcode))
colnames(test2)[1] <- "Barcode"
overlap <- merge(test, test2, by.x = "Barcode", by.y = "Barcode", all = FALSE)

#calculate days after planting
plantdate_even <- as.Date("2015-02-08")
plantdate_odd <- as.Date("2015-02-09")

even$DAP = NA
odd$DAP = NA

even$DAP = as.integer(difftime(even$date, plantdate_even, units = "days"))
odd$DAP = as.integer(difftime(odd$date, plantdate_odd, units = "days"))

#recombine data that is now labeled as odd or even lanes with calculated DAP
DFs <- rbind(odd, even)
DFs <- DFs[order(DFs$DAP),]


#read in barcode csv file and keep relevant info
barcodes <- read.csv("/Users/eludwig/Desktop/brachy/brachy-drought-barcodes.csv")
barcodes = select(barcodes, Barcode, Genotype.ID, Treatment.2, Replicate)
colnames(barcodes)[1] <- "plantbarcode"
colnames(barcodes)[2] <- "Accession"

DFs <- left_join(DFs, barcodes, by = "plantbarcode")
length(unique(DFs$Accession))
DFs <- DFs[order(DFs$DAP, DFs$time, DFs$Accession),]

#read in locations csv file and keep relevant info
locations <- read.csv("brachy-locations.csv")
locations = select(locations, Genotype, Collection_location, Latitude, Longitude, Elevation)
colnames(locations)[1] <- "Accession"

DFs <- left_join(DFs, locations, by = "Accession")
length(unique(DFs$Accession))

DFs$Treatment = NA
DFs$Treatment[grep("20", DFs$Treatment.2)] <- "WL"
DFs$Treatment[grep("100", DFs$Treatment.2)] <- "WW"

DFs$percent_damage = NA
DFs$percent_damage <- ((DFs$cor_unhealthy_area/DFs$cor_total_plant_area)*100)

DFs <- select(DFs, image, plantbarcode, Accession, date, time, DAP, Treatment, Treatment.2, Replicate, Collection_location, Latitude, 
                  Longitude, Elevation, hue_circular_mean, hue_circular_std, hue_median, cor_area, cor_width, cor_height, 
                  cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, 
                  cor_area_below_reference, percent_area_below_reference, cor_convex_hull_area, solidity, cor_perimeter, 
                  cor_longest_path, convex_hull_vertices, cor_ellipse_major_axis, cor_ellipse_minor_axis, ellipse_angle, 
                  ellipse_eccentricity, cor_unhealthy_area, cor_total_plant_area, percent_damage)

#remove a few images with issues (plant is missing)
DFs[DFs=='NA'] <- NA
DFs <- DFs[! is.na(DFs$hue_circular_mean),]

# a few imaging jobs were interrupted during experiment so they extended to next day, combine them with previous day for ease of analysis
DFs$DAP[grep("33", DFs$DAP)] <- "32"
DFs$DAP[grep("37", DFs$DAP)] <- "36"



#save the csv file of the compiled VIS data
write.csv(DFs, "./brachy_control_VIS_SV_single_data_compiled_currated.csv", row.names = FALSE)



#####################################################################
# Read in the heat/heat-drought single-value data #
#####################################################################

setwd("/Users/eludwig/Desktop/brachy/")
#read in single value trait csv files
heatz2500s <- read.csv("heatdrought_SV_z2500-single-value-traits.csv")
heatz3500s <- read.csv("heatdrought_SV_z3500-single-value-traits.csv")

heatz2500s = subset(heatz2500s, select = -c(percent_unhealthy))
heatz3500s = subset(heatz3500s, select = -c(percent_unhealthy))

#remove duplicates of all rows due to re-running plantcv and rows being added to the end of original file instead of file being overwritten
heatz2500s <- heatz2500s %>% distinct(plantbarcode, timestamp, frame, imgtype, .keep_all = TRUE)
heatz3500s <- heatz3500s[44405:88808,] #drops too many doing it with distinct
rownames(heatz3500s) <- NULL

#subset VIS data
heatvisz2500s = heatz2500s[heatz2500s$imgtype == 'VIS',]
heatvisz3500s = heatz3500s[heatz3500s$imgtype == 'VIS',]


#read in damage trait csv files
heatz2500s_damage <- read.csv("heatdrought_SV_z2500_damage-single-value-traits.csv")
heatz3500s_damage <- read.csv("heatdrought_SV_z3500_damage-single-value-traits.csv")

#subset VIS data
heatvisz2500s_damage = heatz2500s_damage[heatz2500s_damage$imgtype == 'VIS',]
heatvisz3500s_damage = heatz3500s_damage[heatz3500s_damage$imgtype == 'VIS',]


heatvisz2500s <- inner_join(heatvisz2500s, heatvisz2500s_damage, by = c("camera", "imgtype", "zoom", "exposure", "gain", "frame", "lifter", 
                                                            "timestamp", "id", "plantbarcode", "treatment", "cartag", "measurementlabel", "other", "image"))
heatvisz3500s <- inner_join(heatvisz3500s, heatvisz3500s_damage, by = c("camera", "imgtype", "zoom", "exposure", "gain", "frame", "lifter", 
                                                            "timestamp", "id", "plantbarcode", "treatment", "cartag", "measurementlabel", "other", "image"))


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
heatvisz2500s <- transform(heatvisz2500s, cor_area = area / 5543.662)
heatvisz2500s <- transform(heatvisz2500s, cor_area_above_reference = area_above_reference / 5543.662)
heatvisz2500s <- transform(heatvisz2500s, cor_area_below_reference = area_below_reference / 5543.662)
heatvisz2500s <- transform(heatvisz2500s, cor_unhealthy_area = unhealthy_area / 5543.662)
heatvisz2500s <- transform(heatvisz2500s, cor_total_plant_area = total_plant_area / 5543.662)

heatvisz2500s <- transform(heatvisz2500s, cor_width = width / 71.89759)
heatvisz2500s <- transform(heatvisz2500s, cor_height = height / 71.89759)
heatvisz2500s <- transform(heatvisz2500s, cor_height_above_reference = height_above_reference / 71.89759)
heatvisz2500s <- transform(heatvisz2500s, cor_height_below_reference = height_below_reference / 71.89759)
heatvisz2500s <- transform(heatvisz2500s, cor_convex_hull_area = convex_hull_area / 5543.662)
heatvisz2500s <- transform(heatvisz2500s, cor_perimeter = perimeter / 71.89759)
heatvisz2500s <- transform(heatvisz2500s, cor_longest_path = longest_path / 71.89759)
heatvisz2500s <- transform(heatvisz2500s, cor_ellipse_major_axis = ellipse_major_axis / 71.89759)
heatvisz2500s <- transform(heatvisz2500s, cor_ellipse_minor_axis = ellipse_minor_axis / 71.89759)

# at z3500
heatvisz3500s <- transform(heatvisz3500s, cor_area = area / 11772.169)
heatvisz3500s <- transform(heatvisz3500s, cor_area_above_reference = area_above_reference / 11772.169)
heatvisz3500s <- transform(heatvisz3500s, cor_area_below_reference = area_below_reference / 11772.169)
heatvisz3500s <- transform(heatvisz3500s, cor_unhealthy_area = unhealthy_area / 11772.169)
heatvisz3500s <- transform(heatvisz3500s, cor_total_plant_area = total_plant_area / 11772.169)

heatvisz3500s <- transform(heatvisz3500s, cor_width = width / 104.61683)
heatvisz3500s <- transform(heatvisz3500s, cor_height = height / 104.61683)
heatvisz3500s <- transform(heatvisz3500s, cor_height_above_reference = height_above_reference / 104.61683)
heatvisz3500s <- transform(heatvisz3500s, cor_height_below_reference = height_below_reference / 104.61683)
heatvisz3500s <- transform(heatvisz3500s, cor_convex_hull_area = convex_hull_area / 11772.169)
heatvisz3500s <- transform(heatvisz3500s, cor_perimeter = perimeter / 104.61683)
heatvisz3500s <- transform(heatvisz3500s, cor_longest_path = longest_path / 104.61683)
heatvisz3500s <- transform(heatvisz3500s, cor_ellipse_major_axis = ellipse_major_axis / 104.61683)
heatvisz3500s <- transform(heatvisz3500s, cor_ellipse_minor_axis = ellipse_minor_axis / 104.61683)


#####################################################################
# Heat/heat-drought data filtering and organizing #
#####################################################################

heatDFs <- rbind(heatvisz3500s, heatvisz2500s)

#filter data
heatDFs <- select(heatDFs, image, plantbarcode, timestamp, hue_circular_mean, hue_circular_std, hue_median, cor_area, cor_width, cor_height, 
              cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, 
              cor_area_below_reference, percent_area_below_reference, cor_convex_hull_area, solidity, cor_perimeter, cor_longest_path, 
              convex_hull_vertices, cor_ellipse_major_axis, cor_ellipse_minor_axis, ellipse_angle, ellipse_eccentricity, cor_unhealthy_area, 
              cor_total_plant_area)
heatDFs <- filter(heatDFs, plantbarcode != "Cc00AA000000") #remove any images of the color card

#separate date and time
TS <- c("date","time") #make a vector for the column names of each part of the timestamp
heatDFs <- separate(heatDFs, col = timestamp, into = c("date","time"), sep = " ") #separate each component of timestamp

#change date format to R standard 
heatDFs$date <- as.Date(heatDFs$date, format="%m/%d/%y")

#read in barcode csv file with lane information to calculate DAP (since split into two groups for planting and imaging)
heatbarcodes_even <- read.csv("/Users/eludwig/Desktop/brachy/B001dr_081214_barcodes_statefile_order_even.csv", header= FALSE)
colnames(heatbarcodes_even)[1] <- "Barcode"

heat_even <- (unique(heatbarcodes_even$Barcode))
heatDFs$compare.barcode <- heatDFs$plantbarcode %in% heat_even
heat_even <- filter(heatDFs, compare.barcode == TRUE)
heat_even$lane_group <- "even"
heat_even <- select(heat_even, -compare.barcode)

heatbarcodes_odd <- read.csv("/Users/eludwig/Desktop/brachy/B001dr_081214_barcodes_statefile_order_odd.csv", header= FALSE)
colnames(heatbarcodes_odd)[1] <- "Barcode"

heat_odd <- (unique(heatbarcodes_odd$Barcode))
heatDFs$compare.barcode2 <- heatDFs$plantbarcode %in% heat_odd
heat_odd <- filter(heatDFs, compare.barcode2 == TRUE)
heat_odd$lane_group <- "odd"
heat_odd <- select(heat_odd, -compare.barcode, -compare.barcode2)

#check that there isn't overlap between odd and even groups (overlap should be 0)
test <- as.data.frame(unique(heat_even$plantbarcode))
colnames(test)[1] <- "Barcode"
test2 <- as.data.frame(unique(heat_odd$plantbarcode))
colnames(test2)[1] <- "Barcode"
overlap <- merge(test, test2, by.x = "Barcode", by.y = "Barcode", all = FALSE)

#calculate days after planting
heatplantdate_even <- as.Date("2014-07-28")
heatplantdate_odd <- as.Date("2014-07-29")

heat_even$DAP = NA
heat_odd$DAP = NA

heat_even$DAP = as.integer(difftime(heat_even$date, heatplantdate_even, units = "days"))
heat_odd$DAP = as.integer(difftime(heat_odd$date, heatplantdate_odd, units = "days"))

#recombine data that is now labeled as odd or even lanes with calculated DAP
heatDFs <- rbind(heat_odd, heat_even)
heatDFs <- heatDFs[order(heatDFs$DAP),]


#read in barcode csv file and keep relevant info
heatbarcodes <- read.csv("/Users/eludwig/Desktop/brachy/brachy-heatdrought-barcodes.csv")
heatbarcodes = select(heatbarcodes, Barcode, Genotype.ID, Treatment.2, Replicate)
colnames(heatbarcodes)[1] <- "plantbarcode"
colnames(heatbarcodes)[2] <- "Accession"

heatDFs <- left_join(heatDFs, heatbarcodes, by = "plantbarcode")
length(unique(heatDFs$Accession))
heatDFs <- heatDFs[order(heatDFs$DAP, heatDFs$time, heatDFs$Accession),]

#read in locations csv file and keep relevant info
locations <- read.csv("brachy-locations.csv")
locations = select(locations, Genotype, Collection_location, Latitude, Longitude, Elevation)
colnames(locations)[1] <- "Accession"

heatDFs <- left_join(heatDFs, locations, by = "Accession")
length(unique(heatDFs$Accession))

heatDFs$Treatment = NA
heatDFs$Treatment[grep("20", heatDFs$Treatment.2)] <- "WL"
heatDFs$Treatment[grep("100", heatDFs$Treatment.2)] <- "WW"

heatDFs$percent_damage = NA
heatDFs$percent_damage <- ((heatDFs$cor_unhealthy_area/heatDFs$cor_total_plant_area)*100)

heatDFs <- select(heatDFs, image, plantbarcode, Accession, date, time, DAP, Treatment, Treatment.2, Replicate, Collection_location, Latitude, 
                  Longitude, Elevation, hue_circular_mean, hue_circular_std, hue_median, cor_area, cor_width, cor_height, 
                  cor_height_above_reference, cor_height_below_reference, cor_area_above_reference, percent_area_above_reference, 
                  cor_area_below_reference, percent_area_below_reference, cor_convex_hull_area, solidity, cor_perimeter, 
                  cor_longest_path, convex_hull_vertices, cor_ellipse_major_axis, cor_ellipse_minor_axis, ellipse_angle, 
                  ellipse_eccentricity, cor_unhealthy_area, cor_total_plant_area, percent_damage)

#remove a few images with issues (plant is missing)
heatDFs[heatDFs=='NA'] <- NA
heatDFs <- heatDFs[! is.na(heatDFs$hue_circular_mean),]

# a few imaging jobs were interrupted during experiment so they extended to next day, combine them with previous day for ease of analysis
heatDFs$DAP[grep("28", heatDFs$DAP)] <- "27"
heatDFs$DAP[grep("32", heatDFs$DAP)] <- "31"

#save the csv file of the compiled VIS data
write.csv(heatDFs, "./brachy_heat_VIS_SV_single_data_compiled_currated.csv", row.names = FALSE)

#####################################################################
# Filtering dataframes to include only accessions in both experiments #
#####################################################################

#compare accessions in control and heat treatments and keep only the ones in both
geno_control <- as.data.frame(unique(DFs$Accession))
colnames(geno_control)[1] <- "Accession"
geno_heat <- as.data.frame(unique(heatDFs$Accession))
colnames(geno_heat)[1] <- "Accession"
bothDF <- merge(geno_control, geno_heat, by.x = "Accession", by.y = "Accession", all = FALSE)
length(unique(bothDF$Accession))
in_both <- unique(bothDF$Accession)

DFs$compare.geno <- DFs$Accession %in% in_both
heatDFs$compare.geno <- heatDFs$Accession %in% in_both
DFs_132 <- filter(DFs, compare.geno == TRUE)
heatDFs_132 <- filter(heatDFs, compare.geno == TRUE)
DFs_132 <- select(DFs_132, -compare.geno)
heatDFs_132 <- select(heatDFs_132, -compare.geno)
DFs <- select(DFs, -compare.geno)
heatDFs <- select(heatDFs, -compare.geno)

#write csv with data from only genotypes that are in both experiments
write.csv(DFs_132, "./brachy_control_single_data_132_accessions.csv", row.names = FALSE)

#write csv with data from only genotypes that are in both experiments
write.csv(heatDFs_132, "./brachy_heat_single_data_132_accessions.csv", row.names = FALSE)
