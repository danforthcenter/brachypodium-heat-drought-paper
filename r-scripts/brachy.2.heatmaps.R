rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(gplots)


setwd("/Users/eludwig/Desktop/brachy/")

DFs_132 <- read.csv("brachy_control_single_data_132_accessions.csv")

#remove days 46 and 47 (data incomplete because half the plants were harvested)
DFs_132 <- DFs_132[!(DFs_132$DAP=="46"),]
DFs_132 <- DFs_132[!(DFs_132$DAP=="47"),]

heatDFs_132 <- read.csv("brachy_heat_single_data_132_accessions.csv")

###########################################
# heatmap of percentage of unhealthy plant for control and drought #
###########################################
avg_damage <- group_by(DFs_132, Accession, Treatment, DAP) %>% summarise(avg_percent_damage = mean(percent_damage, na.rm = TRUE))
avg_damage <- ungroup(avg_damage)

ctl_damage <- filter(avg_damage, Treatment == "WW") #make a DF with all images in the control treatment
drt_damage <- filter (avg_damage, Treatment == "WL") #make a DF with all images in the drought treamtment
ctl_damage$DAP <- as.numeric(ctl_damage$DAP)
drt_damage$DAP <- as.numeric(drt_damage$DAP)

#create a heatmap with a dendrogram
plot_damage <- select(ctl_damage, Accession, DAP, avg_percent_damage) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plot_damage$avg_percent_damage <- as.numeric(plot_damage$avg_percent_damage) #make sure the ratio is numeric
spread_damage <- spread(plot_damage, key = DAP, value = avg_percent_damage) #this flips the dataframe so the days are across the top and the Accessions are down the first row
row.names_damage <- spread_damage[1:132,1] # save the row names
row.names_damage <- t(row.names_damage)
row.names_damage <- as.character(row.names_damage)
spread_damage <- spread_damage[,2:16] # save the DF without the first column
plot.mtx_damage <- as.matrix(spread_damage) #make the DF a matrix
rownames(plot.mtx_damage) <- c(row.names_damage) #add back row names
# plot the heatmap
#heatmap_height <- heatmap(plot.mtx_height,Colv=NA)
# plot the interactive heatmap
color.palette = colorRampPalette(c("blue","lightyellow","darkorange1"),space="rgb")
col_breaks = c(seq(0,15,length=100), # for blue
               seq(16,25,length=100),  # for yellow
               seq(26,100,length=100)) # for orange
pdf(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/damage_heatmap_control.pdf",width = 12,height = 18)
png(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/damage_heatmap_control.png",width = 900,height = 1200)
heatmapunhealthy <- heatmap.2(plot.mtx_damage, 
                              Colv=FALSE,
                              Rowv = FALSE,
                              dendrogram = "none",
                              xlab = "Days After Planting (DAP)",
                              ylab = "Accession",
                              main = "Percent of Damaged Tissue in Plants in Control Conditions",
                              col = color.palette,
                              breaks = col_breaks,
                              #na.color = "black",
                              #margin= c(5,5),
                              srtCol =0,
                              adjCol = c(0.5,.5),
                              #adjRow = c(0.1,0.5),
                              offsetRow = 0,
                              key.ylab = NA,
                              key.xlab = "Percent Damage",
                              keysize=2,
                              lhei=c(2, 15),
                              lwid = c(2,6),
                              density.info = "none",
                              trace = "none") # plot the heatmap
dev.off()



#create a heatmap with a dendrogram
plot_damage1 <- select(drt_damage, Accession, DAP, avg_percent_damage) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plot_damage1$avg_percent_damage <- as.numeric(plot_damage1$avg_percent_damage) #make sure the ratio is numeric
spread_damage1 <- spread(plot_damage1, key = DAP, value = avg_percent_damage) #this flips the dataframe so the days are across the top and the Accessions are down the first row
row.names_damage1 <- spread_damage1[1:132,1] # save the row names
row.names_damage1 <- t(row.names_damage1)
row.names_damage1 <- as.character(row.names_damage1)
spread_damage1 <- spread_damage1[,2:16] # save the DF without the first column
plot.mtx_damage1 <- as.matrix(spread_damage1) #make the DF a matrix
rownames(plot.mtx_damage1) <- c(row.names_damage1) #add back row names
# plot the heatmap
#heatmap_height <- heatmap(plot.mtx_height,Colv=NA)
# plot the interactive heatmap
color.palette = colorRampPalette(c("blue","lightyellow","darkorange1"),space="rgb")
col_breaks = c(seq(0,15,length=100), # for blue
               seq(16,25,length=100),  # for yellow
               seq(26,100,length=100)) # for orange
pdf(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/damage_heatmap_drought.pdf",width = 12,height = 18)
png(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/damage_heatmap_drought.png",width = 900,height = 1200)
heatmapunhealthy1 <- heatmap.2(plot.mtx_damage1, 
                               Colv=FALSE,
                               Rowv = FALSE,
                               dendrogram = "none",
                               xlab = "Days After Planting (DAP)",
                               ylab = "Accession",
                               main = "Percent of Unhealthy Tissue in Plants in Drought Conditions",
                               col = color.palette,
                               breaks = col_breaks,
                               #na.color = "black",
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


###########################################
# heatmap of percentage of unhealthy plant for heat and heat-drought #
###########################################
#heatmap of heat unhealthy
heatDFs_132$percent_damage <- as.numeric(heatDFs_132$percent_damage)
avg_damageheat <- group_by(heatDFs_132, Accession, Treatment, DAP) %>% summarise(avg_percent_damage = mean(percent_damage, na.rm=TRUE))
avg_damageheat <- ungroup(avg_damageheat)

ctl_damageheat <- filter(avg_damageheat, Treatment == "WW") #make a DF with all images in the control treatment
drt_damageheat <- filter (avg_damageheat, Treatment == "WL") #make a DF with all images in the drought treamtment

missing1 <- filter(ctl_damageheat, Accession == "Bar2")
missing2 <- filter(ctl_damageheat, Accession == "BdTR1H")
missing3 <- filter(ctl_damageheat, Accession == "BdTR3S")
missing <- rbind(missing1, missing2)
missing <- rbind(missing, missing3)
missing$Treatment[grep("WW", missing$Treatment)] <- "WL"
missing$avg_percent_damage = NA
drt_damageheat <- merge(drt_damageheat, missing, by = c("Accession", "Treatment", "DAP", "avg_percent_damage"), all= TRUE)
length(unique(drt_damageheat$Accession))
ctl_damageheat$DAP <- as.numeric(ctl_damageheat$DAP)
drt_damageheat$DAP <- as.numeric(drt_damageheat$DAP)
#remove imgday 3 (data super messy)
ctl_damageheat <- ctl_damageheat[!(ctl_damageheat$DAP=="3"),]

#create a heatmap with a dendrogram
plot_damageheat <- select(ctl_damageheat, Accession, DAP, avg_percent_damage) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plot_damageheat$avg_percent_damage <- as.numeric(plot_damageheat$avg_percent_damage) #make sure the ratio is numeric
spread_damageheat <- spread(plot_damageheat, key = DAP, value = avg_percent_damage) #this flips the dataframe so the days are across the top and the Accessions are down the first row
row.names_damageheat <- spread_damageheat[1:132,1] # save the row names
row.names_damageheat <- t(row.names_damageheat)
row.names_damageheat <- as.character(row.names_damageheat)
spread_damageheat <- spread_damageheat[,2:19] # save the DF without the first column
plot.mtx_damageheat <- as.matrix(spread_damageheat) #make the DF a matrix
rownames(plot.mtx_damageheat) <- c(row.names_damageheat) #add back row names
color.palette = colorRampPalette(c("blue","lightyellow","darkorange1"),space="rgb")
col_breaks = c(seq(0,15,length=100), # for blue
               seq(16,25,length=100),  # for yellow
               seq(26,100,length=100)) # for orange
pdf(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/damage_heatmap_heat.pdf",width = 12,height = 18)
png(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/damage_heatmap_heat.png",width = 900,height = 1200)
heatmapunhealthyheat <- heatmap.2(plot.mtx_damageheat, 
                                  Colv=FALSE,
                                  Rowv = FALSE,
                                  dendrogram = "none",
                                  xlab = "Days After Planting (DAP)",
                                  ylab = "Accession",
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
#heatmapunhealthyheat



#create a heatmap for heat & drought unhealthy
plot_damageheat1 <- select(drt_damageheat, Accession, DAP, avg_percent_damage) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
length(unique(avg_damageheat$Accession))
plot_damageheat1$avg_percent_damage <- as.numeric(plot_damageheat1$avg_percent_damage) #make sure the ratio is numeric
spread_damageheat1 <- spread(plot_damageheat1, key = DAP, value = avg_percent_damage) #this flips the dataframe so the days are across the top and the Accessions are down the first row
row.names_damageheat1 <- spread_damageheat1[1:132,1] # save the row names
row.names_damageheat1 <- t(row.names_damageheat1)
row.names_damageheat1 <- as.character(row.names_damageheat1)
spread_damageheat1 <- spread_damageheat1[,2:19] # save the DF without the first column
plot.mtx_damageheat1 <- as.matrix(spread_damageheat1) #make the DF a matrix
rownames(plot.mtx_damageheat1) <- c(row.names_damageheat1) #add back row names
color.palette = colorRampPalette(c("blue","lightyellow","darkorange1"),space="rgb")
col_breaks = c(seq(0,15,length=100), # for blue
               seq(16,25,length=100),  # for yellow
               seq(26,100,length=100)) # for orange
pdf(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/damage_heatmap_heat&drought.pdf",width = 12,height = 18)
png(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/damage_heatmap_heat&drought.png",width = 900,height = 1200)
heatmapunhealthyheat1 <- heatmap.2(plot.mtx_damageheat1, 
                                   Colv=FALSE,
                                   Rowv = FALSE,
                                   dendrogram = "none",
                                   xlab = "Days After Planting (DAP)",
                                   ylab = "Accession",
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





###########################################
# biomass heatmaps of control and drought conditions #
###########################################



#calculate average area and create avg df
avgDFs <- group_by(DFs_132, Accession, Treatment, DAP) %>% summarise(avg_area = mean(cor_area))
avgDFs <- ungroup(avgDFs)


#heatmap of biomass in control
area <- filter(avgDFs, Treatment == "WW") #make a DF with all images in the control treatment
area$DAP <- as.numeric(area$DAP)

#create a heatmap with a dendrogram
plotDF <- select(area, Accession, DAP, avg_area) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDF$avg_area <- as.numeric(plotDF$avg_area) #make sure the area is numeric
spreadDF <- spread(plotDF, key = DAP, value = avg_area) #this flips the dataframe so the days are across the top and the Accessions are down the first row
row.names <- spreadDF[1:132,1] # save the row names
row.names <- t(row.names)
row.names <- as.character(row.names)
spreadDF <- spreadDF[,2:16] # save the DF without the first column
plot.mtx <- as.matrix(spreadDF) #make the DF a matrix
rownames(plot.mtx) <- c(row.names) #add back row names
# plot the heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,5,length=40), # for orange
               seq(6,9,length=40),  # for yellow
               seq(10,40,length=40)) # for blue
pdf(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/area_heatmap_control.pdf",width = 12,height = 18)
png(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/area_heatmap_control.png",width = 900,height = 1200)
heatmapctl <- heatmap.2(plot.mtx, 
                        Colv=FALSE,
                        Rowv = FALSE,
                        dendrogram = "none",
                        xlab = "Days After Planting (DAP)",
                        ylab = "Accession",
                        main = "Heatmap of Average Plant Area in Control",
                        col = color.palette,
                        breaks = col_breaks,
                        #na.color = "black",
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
areadrt$DAP <- as.numeric(areadrt$DAP)

#create a heatmap with a dendrogram
plotDF1 <- select(areadrt, Accession, DAP, avg_area) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDF1$avg_area <- as.numeric(plotDF1$avg_area) #make sure the area is numeric
spreadDF1 <- spread(plotDF1, key = DAP, value = avg_area) #this flips the dataframe so the days are across the top and the Accessions are down the first row
row.names1 <- spreadDF1[1:132,1] # save the row names
row.names1 <- t(row.names1)
row.names1 <- as.character(row.names1)
spreadDF1 <- spreadDF1[,2:16] # save the DF without the first column
plot.mtx1 <- as.matrix(spreadDF1) #make the DF a matrix
rownames(plot.mtx1) <- c(row.names1) #add back row names
# plot the heatmap
#heatmap <- heatmap(plot.mtx,Colv=NA)
# plot the interactive heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,5,length=40), # for orange
               seq(6,9,length=40),  # for yellow
               seq(10,40,length=40)) # for blue
pdf(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/area_heatmap_drought.pdf",width = 12,height = 18)
png(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/area_heatmap_drought.png",width = 900,height = 1200)
heatmapdrt <- heatmap.2(plot.mtx1, 
                        Colv=FALSE,
                        Rowv = FALSE,
                        dendrogram = "none",
                        xlab = "Days After Planting (DAP)",
                        ylab = "Accession",
                        main = "Heatmap of Average Plant Area in Drought",
                        col = color.palette,
                        breaks = col_breaks,
                        #na.color = "black",
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



###########################################
# biomass heatmaps of heat and heat-drought conditions #
###########################################
#biomass heatmap for heat
#calculate average area and create avg df
avgheatDFs <- group_by(heatDFs_132, Accession, Treatment, DAP) %>% summarise(avg_area = mean(cor_area))
avgheatDFs <- ungroup(avgheatDFs)

heatarea <- filter(avgheatDFs, Treatment == "WW") #make a DF with all images in the control treatment
heatarea$DAP <- as.numeric(heatarea$DAP)


#create a heatmap with a dendrogram
plotDFheat <- select(heatarea, Accession, DAP, avg_area) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDFheat$avg_area <- as.numeric(plotDFheat$avg_area) #make sure the area is numeric
spreadDFheat <- spread(plotDFheat, key = DAP, value = avg_area) #this flips the dataframe so the days are across the top and the Accessions are down the first row
row.namesheat <- spreadDFheat[1:132,1] # save the row names
row.namesheat <- t(row.namesheat)
row.namesheat <- as.character(row.namesheat)
spreadDFheat <- spreadDFheat[,2:19] # save the DF without the first column
plot.mtxheat <- as.matrix(spreadDFheat) #make the DF a matrix
rownames(plot.mtxheat) <- c(row.namesheat) #add back row names
# plot the heatmap
#heatmap <- heatmap(plot.mtx,Colv=NA)
# plot the interactive heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,5,length=40), # for orange
               seq(6,9,length=40),  # for yellow
               seq(10,40,length=40)) # for blue
pdf(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/area_heatmap_heat.pdf",width = 12,height = 18)
png(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/area_heatmap_heat.png",width = 900,height = 1200)
heatmapheat <- heatmap.2(plot.mtxheat, 
                         Colv=FALSE,
                         Rowv = FALSE,
                         dendrogram = "none",
                         xlab = "Days After Planting (DAP)",
                         ylab = "Accession",
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




#heatmap for biomass in heat and drought
heatdrtarea <- filter (avgheatDFs, Treatment == "WL") #make a DF with all area means in the drought treatment
heatdrtarea$DAP <- as.numeric(heatdrtarea$DAP)

missing1 <- filter(heatarea, Accession == "Bar2")
missing2 <- filter(heatarea, Accession == "BdTR1H")
missing3 <- filter(heatarea, Accession == "BdTR3S")
missing <- rbind(missing1, missing2)
missing <- rbind(missing, missing3)
missing$Treatment[grep("WW", missing$Treatment)] <- "WL"
missing$avg_area = NA
heatdrtarea <- merge(heatdrtarea, missing, by = c("Accession", "Treatment", "DAP", "avg_area"), all= TRUE)


#create a heatmap with a dendrogram
plotDFheatdrt <- select(heatdrtarea, Accession, DAP, avg_area) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDFheatdrt$avg_area <- as.numeric(plotDFheatdrt$avg_area) #make sure the area is numeric
spreadDFheatdrt <- spread(plotDFheatdrt, key = DAP, value = avg_area) #this flips the dataframe so the days are across the top and the Accessions are down the first row
row.namesheatdrt <- spreadDFheatdrt[1:132,1] # save the row names
row.namesheatdrt <- t(row.namesheatdrt)
row.namesheatdrt <- as.character(row.namesheatdrt)
spreadDFheatdrt <- spreadDFheatdrt[,2:19] # save the DF without the first column
plot.mtxheatdrt <- as.matrix(spreadDFheatdrt) #make the DF a matrix
rownames(plot.mtxheatdrt) <- c(row.namesheatdrt) #add back row names
# plot the heatmap
#heatmap <- heatmap(plot.mtx,Colv=NA)
# plot the interactive heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,5,length=40), # for orange
               seq(6,9,length=40),  # for yellow
               seq(10,40,length=40)) # for blue
pdf(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/area_heatmap_heat&drought.pdf",width = 12,height = 18)
png(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/area_heatmap_heat&drought.png",width = 900,height = 1200)
heatmapheatdrt <- heatmap.2(plot.mtxheatdrt, 
                            Colv=FALSE,
                            Rowv = FALSE,
                            dendrogram = "none",
                            xlab = "Days After Planting (DAP)",
                            ylab = "Accession",
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



###########################################
# height heatmaps of control and drought conditions #
###########################################

#calculate average area and create avg df
avgDFs_height <- group_by(DFs_132, Accession, Treatment, DAP) %>% summarise(avg_height = mean(cor_height_above_reference))
avgDFs_height <- ungroup(avgDFs_height)

#heatmap of biomass in control
height <- filter(avgDFs_height, Treatment == "WW") #make a DF with all images in the control treatment
height$DAP <- as.numeric(height$DAP)

#create a heatmap with a dendrogram
plotDFheight <- select(height, Accession, DAP, avg_height) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDFheight$avg_height <- as.numeric(plotDFheight$avg_height) #make sure the area is numeric
spreadDFheight <- spread(plotDFheight, key = DAP, value = avg_height) #this flips the dataframe so the days are across the top and the Accessions are down the first row
row.namesheight <- spreadDFheight[1:132,1] # save the row names
row.namesheight <- t(row.namesheight)
row.namesheight <- as.character(row.namesheight)
spreadDFheight <- spreadDFheight[,2:16] # save the DF without the first column
plot.mtxheight <- as.matrix(spreadDFheight) #make the DF a matrix
rownames(plot.mtxheight) <- c(row.namesheight) #add back row names
# plot the heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,6,length=20), # for orange
               seq(7,9,length=20),  # for yellow
               seq(10,20,length=20)) # for blue
pdf(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/height_heatmap_control.pdf",width = 12,height = 18)
png(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/height_heatmap_control.png",width = 900,height = 1200)
heatmapctlheight <- heatmap.2(plot.mtxheight, 
                              Colv=FALSE,
                              Rowv = FALSE,
                              dendrogram = "none",
                              xlab = "Days After Planting (DAP)",
                              ylab = "Accession",
                              main = "Average Plant Height in Control",
                              col = color.palette,
                              breaks = col_breaks,
                              #na.color = "black",
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
heightdrt$DAP <- as.numeric(heightdrt$DAP)

#create a heatmap with a dendrogram
plotDFheight1 <- select(heightdrt, Accession, DAP, avg_height) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDFheight1$avg_height <- as.numeric(plotDFheight1$avg_height) #make sure the area is numeric
spreadDFheight1 <- spread(plotDFheight1, key = DAP, value = avg_height) #this flips the dataframe so the days are across the top and the Accessions are down the first row
row.namesheight1 <- spreadDFheight1[1:132,1] # save the row names
row.namesheight1 <- t(row.namesheight1)
row.namesheight1 <- as.character(row.namesheight1)
spreadDFheight1 <- spreadDFheight1[,2:16] # save the DF without the first column
plot.mtxheight1 <- as.matrix(spreadDFheight1) #make the DF a matrix
rownames(plot.mtxheight1) <- c(row.namesheight1) #add back row names
# plot the heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,6,length=20), # for orange
               seq(7,9,length=20),  # for yellow
               seq(10,20,length=20)) # for blue
pdf(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/height_heatmap_drought.pdf",width = 12,height = 18)
png(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/height_heatmap_drought.png",width = 900,height = 1200)
heatmapdrtheight <- heatmap.2(plot.mtxheight1, 
                              Colv=FALSE,
                              Rowv = FALSE,
                              dendrogram = "none",
                              xlab = "Days After Planting (DAP)",
                              ylab = "Accession",
                              #main = "Average Plant Height in Drought",
                              col = color.palette,
                              breaks = col_breaks,
                              #na.color = "black",
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



###########################################
# height heatmaps of heat and heat-drought conditions #
###########################################

#calculate average area and create avg df
avgheatDFs_height <- group_by(heatDFs_132, Accession, Treatment, DAP) %>% summarise(avg_height = mean(cor_height_above_reference))
avgheatDFs_height <- ungroup(avgheatDFs_height)

#height heatmap for heat
heatheight <- filter(avgheatDFs_height, Treatment == "WW") #make a DF with all images in the control treatment
heatheight$DAP <- as.numeric(heatheight$DAP)

#create a heatmap 
plotDFheatheight <- select(heatheight, Accession, DAP, avg_height) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDFheatheight$avg_height <- as.numeric(plotDFheatheight$avg_height) #make sure the area is numeric
spreadDFheatheight <- spread(plotDFheatheight, key = DAP, value = avg_height) #this flips the dataframe so the days are across the top and the Accessions are down the first row
row.namesheatheight <- spreadDFheatheight[1:132,1] # save the row names
row.namesheatheight <- t(row.namesheatheight)
row.namesheatheight <- as.character(row.namesheatheight)
spreadDFheatheight <- spreadDFheatheight[,2:19] # save the DF without the first column
plot.mtxheatheight <- as.matrix(spreadDFheatheight) #make the DF a matrix
rownames(plot.mtxheatheight) <- c(row.namesheatheight) #add back row names
# plot the heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,6,length=20), # for orange
               seq(7,9,length=20),  # for yellow
               seq(10,20,length=20)) # for blue
pdf(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/height_heatmap_heat.pdf",width = 12,height = 18)
png(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/height_heatmap_heat.png",width = 900,height = 1200)
heatmapheatheight <- heatmap.2(plot.mtxheatheight, 
                               Colv=FALSE,
                               Rowv = FALSE,
                               dendrogram = "none",
                               xlab = "Days After Planting (DAP)",
                               ylab = "Accession",
                               #main = "Average Plant Height in Heat",
                               col = color.palette,
                               breaks = col_breaks,
                               ##na.color = "black",
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




#heatmap for height in heat and drought
heatdrtheight <- filter (avgheatDFs_height, Treatment == "WL") #make a DF with all area means in the drought treatment
heatdrtheight$DAP <- as.numeric(heatdrtheight$DAP)

missing1 <- filter(heatheight, Accession == "Bar2")
missing2 <- filter(heatheight, Accession == "BdTR1H")
missing3 <- filter(heatheight, Accession == "BdTR3S")
missing <- rbind(missing1, missing2, missing3)
missing$Treatment[grep("WW", missing$Treatment)] <- "WL"
missing$avg_height = NA
heatdrtheight <- merge(heatdrtheight, missing, by = c("Accession", "Treatment", "DAP", "avg_height"), all= TRUE)

#create a heatmap with a dendrogram
plotDFheatdrtheight <- select(heatdrtheight, Accession, DAP, avg_height) #select only the data you want in the heatmap (x axis, y axis, and the ratio)
plotDFheatdrtheight$avg_height <- as.numeric(plotDFheatdrtheight$avg_height) #make sure the area is numeric
spreadDFheatdrtheight <- spread(plotDFheatdrtheight, key = DAP, value = avg_height) #this flips the dataframe so the days are across the top and the Accessions are down the first row
row.namesheatdrtheight <- spreadDFheatdrtheight[1:132,1] # save the row names
row.namesheatdrtheight <- t(row.namesheatdrtheight)
row.namesheatdrtheight <- as.character(row.namesheatdrtheight)
spreadDFheatdrtheight <- spreadDFheatdrtheight[,2:19] # save the DF without the first column
plot.mtxheatdrtheight <- as.matrix(spreadDFheatdrtheight) #make the DF a matrix
rownames(plot.mtxheatdrtheight) <- c(row.namesheatdrtheight) #add back row names
# plot the heatmap
#heatmap <- heatmap(plot.mtx,Colv=NA)
# plot the interactive heatmap
color.palette = colorRampPalette(c("darkorange1","lightyellow","blue"),space="rgb")
col_breaks = c(seq(0,5,length=20), # for orange
               seq(6,8,length=20),  # for yellow
               seq(9,20,length=20)) # for blue
pdf(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/height_heatmap_heat&drought.pdf",width = 12,height = 18)
png(file="/Users/eludwig/Google Drive/My Drive/brachy_project/brachy-figures/height_heatmap_heat&drought.png",width = 900,height = 1200)
heatmapheatdrtheight <- heatmap.2(plot.mtxheatdrtheight, 
                                  Colv=FALSE,
                                  Rowv = FALSE,
                                  dendrogram = "none",
                                  xlab = "Days After Planting (DAP)",
                                  ylab = "Accession",
                                  #main = "Average Plant Height in Heat & Drought",
                                  col = color.palette,
                                  breaks = col_breaks,
                                  ##na.color = "black",
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


