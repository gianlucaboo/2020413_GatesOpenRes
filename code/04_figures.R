####SET UP THE WORKING ENVIRONMENT####
#load libraries
library(raster); library(rgdal); library(tcltk); 
library(tmap);library(spatstat); library(sf); library(ggplot2); 
library(grid); library(rgeos); library(osmdata);
library(dplyr); library(purrr); library(gridsample)

####LOAD AND PROCESS DATA####
#select working environment
setwd(tk_choose.dir(default = getwd(), caption = "Select working directory!"))
getwd()


####FIGURE 2###
#The study area 
#load boundaries
country <- st_combine(st_read("sample_design/data/in/cod_bcr_admin_0.shp"))
provinces <- st_read("sample_design/data/in/cod_bcr_admin_1.shp")
provinces <- subset(provinces, admin1Type == "Urban")
provinces_inset <- st_as_sfc(st_bbox(provinces))
rivers <-  st_crop(st_read(paste0("/Users/", unname(Sys.info()[8]), "/Dropbox/congo-democratic-republic-latest-free/gis_osm_water_a_free_1.shp")), provinces)
places <- st_crop(st_read(paste0("/Users/", unname(Sys.info()[8]), "/Dropbox/congo-democratic-republic-latest-free/gis_osm_places_free_1.shp")), provinces)
places_city <- places[places$fclass == "city", drop=F]
places_town <- places[places$fclass == "town", drop=F]
places_capital <- places[places$fclass == "national_capital", drop=F]
#download.file("https://storage.googleapis.com/earthenginepartners-hansen/GFC-2017-v1.5/Hansen_GFC-2017-v1.5_treecover2000_00N_010E.tif", paste0("/Users/", unname(Sys.info()[8]), "/Downloads/hansen.tif"), quiet = FALSE, method="auto", cacheOK=T)
hansen <- raster(paste0("/Users/", unname(Sys.info()[8]), "/Downloads/hansen.tif"))              
hansen <- mask(crop(hansen, as(provinces, "Spatial")), as(provinces, "Spatial"))

tm_shape(provinces, is.master=T, unit="km") +
  tm_borders(col="black", lwd=0.5) +
  tm_shape(brick(hansen))+
  tm_raster(col="hansen", palette="YlGn", legend.show = FALSE) +
  tm_shape(rivers)+
  tm_fill("blue") +
tm_shape(places_capital)+
  tm_symbols(size=0.2, col="black", shape=21)+
  tm_text("name", ymod=-0.4)+
  tm_shape(places_city)+
  tm_symbols(size=0.2, col="black", shape=1)+
  tm_text("name", ymod=-0.4)+
tm_shape(places_town)+
  tm_symbols(size=0.2, col="black", shape=2)+
  tm_text("name", ymod=-0.4) +
tm_shape(provinces, is.master=T, unit="km") +
  tm_borders(col="black", lwd=0.5) +
tm_scale_bar(breaks=50, color.dark = "black", color.light = "black", position=c("left", "bottom"))+
  tm_layout(frame = F)

print(
  tm_shape(country) +
    tm_borders(col = "black", lw=3) +
    tm_shape(provinces_inset) +
    tm_borders(col = "black", lw=3) +
    tm_layout(legend.show = FALSE, bg.color = NA, frame = FALSE),
  vp=viewport(x= 0.15, y= 0.8, width= 0.3, height= 0.3))

###FIGURE 3###
#The gridded sampling frame (with zoom-in in rural and urban location) 

#Settled pixels
settled_area <- raster("sample_design/data/in/cod_2prov_100m_ornl_2016_settlmentbinary.tif")
settled_inset <- rep(list(raster("sample_design/data/in/cod_2prov_100m_ornl_2016_settlmentbinary.tif")), 4)
inset <- st_read(dsn = paste0("/Users/", unname(Sys.info()[8]), "/Dropbox/inserts.shp"))

for (i in 1:4){
  settled_inset[[i]] <- mask(settled_area, as(inset[i, drop=F], "Spatial"))
}

for (i in 1:4){
  settled_inset[[i]]  <-crop(settled_inset[[i]], as(inset[i, drop=F], "Spatial"))
  settled_inset[[i]][ settled_inset[[i]] == 0] <- NA
}

#Settlement layer
settlement_area <- st_read(dsn = paste0("/Users/", unname(Sys.info()[8]), "/Dropbox/settlement_province.shp"))

settlement_inset <- list()
for (i in 1:4){
  settlement_inset[[i]] <- settlement_area[st_intersects(settlement_area, inset[i, drop=F], sparse=F), drop=F]
}

tmap_options(max.raster = c(plot = 14879410, view = 14879410))
tm_shape(brick(settled_area)) +
  tm_raster(col="cod_2prov_100m_ornl_2016_settlmentbinary", pal="gray50", legend.show = FALSE) +
  tm_shape(provinces, is.master=T, unit="km") +
  tm_borders(col="black", lwd=0.5) +
  tm_shape(inset)+
  tm_borders(col="black", lwd=3) +
  tm_scale_bar(breaks=50, color.dark = "black", color.light = "black", position=c("left", "bottom"))+
  tm_layout(frame = F)
view_x <- rev(c(0.15, 0.38, 0.84, 0.61))
view_y <- rev(c(0.8,0.8, 0.15,0.15))

for (i in 1:4) {
  print(
    tm_shape(inset[i, drop=F]) +
      tm_borders(col = "black", lw=03) +
      tm_shape(brick(settled_inset[[i]])) +
      tm_raster(col="cod_2prov_100m_ornl_2016_settlmentbinary", pal="gray50", colorNA="white", legend.show = FALSE)+
      tm_shape(settlement_inset[[i]]) +
      tm_fill("black") +
      tm_shape(inset[i, drop=F]) +
      tm_borders(col = "black", lw=03) +
      tm_layout(legend.show = FALSE, bg.color = NA, frame = FALSE),
    vp=viewport(x= view_x[i], y= view_y[i], width= 0.3, height= 0.3))
}

###FIGURE 4###
###Variance explained for the k-means and map of strata###

#load pca table
covariate_table_short_kmeans <- read.csv("sample_design/data/out/kmeans/covariate_table_short_kmeans.csv")
covariate_table_short_kmeans <- covariate_table_short_kmeans[2:ncol(covariate_table_short_kmeans)]

#within sum of sqaure plot for kmeans 1:10
nc <- 10
seed <- 1
wss <- (nrow(covariate_table_short_kmeans)-1)*sum(apply(covariate_table_short_kmeans,2,var))

for (i in 1:nc){
  set.seed(seed)
  centers_selected <- sample(1:nrow(covariate_table_short_kmeans), i)
  wss[i] <- sum(kmeans(covariate_table_short_kmeans, centers=covariate_table_short_kmeans[centers_selected,], nstart = 20, iter.max = 200, algorithm="Lloyd")$withinss)}
wss <- wss*100/max(wss)

plot(1:nc, wss, type="b", pch=16, main="Within-cluster sum of squares reduction [%]", xlab="Number of clusters", ylab="", axes=T, family="serif")
axis(1, at=1:nc, label=1:nc, tick=F, family="serif")
axis(2, at=wss[c(1, 3, 5, 8)], labels=paste(round(wss[c(1, 3, 5, 8)],1), "%"), tick=F, las=2, family="serif")
points(1, wss[1], pch=1, bkg="white", cex=1.5)
points(3, wss[3], pch=1, bkg="white", cex=1.5)
points(5, wss[5], pch=1, bkg="white", cex=1.5)
points(8, wss[8], pch=1, bkg="white", cex=1.5)

#load kmeans 1:16
kmeans_list <- list.files("sample_design/data/out/kmeans/", pattern=".tif$")
kmeans <- stack(paste0("sample_design/data/out/kmeans/" , kmeans_list))
kmeans <- stack(kmeans$kmeans_3, kmeans$kmeans_5, kmeans$kmeans_8)
kmeans <- mask(kmeans, as(provinces, "Spatial"))

#color palette — http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=8
pal <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf")
#STRATA 3
strata3_inset <- list()
for (i in 1:4){
  strata3_inset[[i]] <- mask(kmeans$kmeans_3, as(inset[i, drop=F], "Spatial"))
}

for (i in 1:4){
  strata3_inset[[i]]  <-crop(strata3_inset[[i]], as(inset[i, drop=F], "Spatial"))
  strata3_inset[[i]][ strata3_inset[[i]] == 0] <- NA
}

tmap_options(max.raster = c(plot = 14879410, view = 14879410))
tm_shape(brick(kmeans$kmeans_3)) +
  tm_raster(col="kmeans_3", palette=pal[c(1:3)], n=3, style="cat")+
  tm_shape(provinces, is.master=T, unit="km") +
  tm_borders(col="black", lwd=0.5) +
  tm_shape(inset)+
  tm_borders(col="black", lwd=3) +
  tm_scale_bar(breaks=50, color.dark = "black", color.light = "black", position=c("left", "bottom"))+
  tm_layout(frame = F)
view_x <- rev(c(0.15, 0.38, 0.84, 0.61))
view_y <- rev(c(0.8,0.8, 0.15,0.15))

for (i in 1:4) {
  pal_insert <- sort(unique(strata3_inset[[i]][]))
  print(
    tm_shape(inset[i, drop=F]) +
      tm_borders(col = "black", lw=03) +
      tm_shape(brick(strata3_inset[[i]])) +
      tm_raster(col="kmeans_3", palette=pal[c(1:3)][pal_insert], n=3, style="cat")+
      tm_shape(inset[i, drop=F]) +
      tm_borders(col = "black", lw=03) +
      tm_layout(legend.show = FALSE, bg.color = NA, frame = FALSE),
    vp=viewport(x= view_x[i], y= view_y[i], width= 0.3, height= 0.3))
}

#STRATA 5#
strata5_inset <- list()
for (i in 1:4){
  strata5_inset[[i]] <- mask(kmeans$kmeans_5, as(inset[i, drop=F], "Spatial"))
}

for (i in 1:4){
  strata5_inset[[i]]  <-crop(strata5_inset[[i]], as(inset[i, drop=F], "Spatial"))
  strata5_inset[[i]][ strata5_inset[[i]] == 0] <- NA
}

tmap_options(max.raster = c(plot = 14879410, view = 14879410))
tm_shape(brick(kmeans$kmeans_5)) +
  tm_raster(col="kmeans_5", palette=pal[c(1, 4, 3, 2, 7)], n=5, style="cat") +
  tm_shape(provinces, is.master=T, unit="km") +
  tm_borders(col="black", lwd=0.5) +
  tm_shape(inset)+
  tm_borders(col="black", lwd=3) +
  tm_scale_bar(breaks=50, color.dark = "black", color.light = "black", position=c("left", "bottom"))+
  tm_layout(frame = F)
view_x <- rev(c(0.15, 0.38, 0.84, 0.61))
view_y <- rev(c(0.8,0.8, 0.15,0.15))

for (i in 1:4) {
  pal_insert <- sort(unique(strata5_inset[[i]][]))
  print(
    tm_shape(inset[i, drop=F]) +
      tm_borders(col = "black", lw=03) +
      tm_shape(brick(strata5_inset[[i]])) +
      tm_raster(col="kmeans_5", palette=pal[c(1, 4, 3, 2, 7)][pal_insert], n=5, style="cat") +
      tm_shape(inset[i, drop=F]) +
      tm_borders(col = "black", lw=03) +
      tm_layout(legend.show = FALSE, bg.color = NA, frame = FALSE),
    vp=viewport(x= view_x[i], y= view_y[i], width= 0.3, height= 0.3))
}

#STRATA 8#
strata8_inset <- list()
for (i in 1:4){
  strata8_inset[[i]] <- mask(kmeans$kmeans_8, as(inset[i, drop=F], "Spatial"))
}

for (i in 1:4){
  strata8_inset[[i]]  <-crop(strata8_inset[[i]], as(inset[i, drop=F], "Spatial"))
  strata8_inset[[i]][ strata5_inset[[i]] == 0] <- NA
}

tmap_options(max.raster = c(plot = 14879410, view = 14879410))
tm_shape(brick(kmeans$kmeans_8)) +
  tm_raster(col="kmeans_8", palette=pal[c(1, 8, 4, 2, 6, 3, 5, 7)], n=8, style="cat")+
  tm_shape(provinces, is.master=T, unit="km") +
  tm_borders(col="black", lwd=0.5) +
  tm_shape(inset)+
  tm_borders(col="black", lwd=3) +
  tm_scale_bar(breaks=50, color.dark = "black", color.light = "black", position=c("left", "bottom"))+
  tm_layout(frame = F)
view_x <- rev(c(0.15, 0.38, 0.84, 0.61))
view_y <- rev(c(0.8,0.8, 0.15,0.15))

for (i in 1:4) {
  pal_insert <- sort(unique(strata8_inset[[i]][]))
  print(
    tm_shape(inset[i, drop=F]) +
      tm_borders(col = "black", lw=03) +
      tm_shape(brick(strata8_inset[[i]])) +
      tm_raster(col="kmeans_8", palette=pal[c(1, 8, 4, 2, 6, 3, 5, 7)][pal_insert], n=8, style="cat")+
      tm_shape(inset[i, drop=F]) +
      tm_borders(col = "black", lw=03) +
      tm_layout(legend.show = FALSE, bg.color = NA, frame = FALSE),
    vp=viewport(x= view_x[i], y= view_y[i], width= 0.3, height= 0.3))
}

#tmap_options(max.raster = c(plot = 14879410, view = 14879410))
strata3 <- list()
for (i in 1:4) {
  pal_insert <- sort(unique(strata3_inset[[i]][]))
  strata3[[i]] <- 
    tm_layout(frame=F, fontfamily="serif", legend.title.fontface=2)+
    tm_shape(brick(strata3_inset[[i]])) +
    tm_raster(col="kmeans_3", palette=pal[c(1:3)][pal_insert], n=3, style="cat")+
    tm_shape(inset[i, drop=F]) +
    tm_borders(col = "black", lw=03)
}

barplot(kmeans$kmeans_3, col=pal[c(1:3)], border=NA, axes=F, xlab="", ylab="", family="serif", horiz=T)
text(unname(freq(kmeans$kmeans_3)[c(1:3),2]), unname(freq(kmeans$kmeans_3)[c(1:3),1]), round(unname(freq(kmeans$kmeans_3)[c(1:3),2])*100/sum(unname(freq(kmeans$kmeans_3)[c(1:3),2])), 2), family="serif") 

strata5 <- list()
for (i in 1:4) {
  pal_insert <- sort(unique(strata5_inset[[i]][]))
  strata5[[i]] <- 
    tm_layout(frame=F, fontfamily="serif", legend.title.fontface=2)+
    tm_shape(brick(strata5_inset[[i]])) +
    tm_raster(col="kmeans_5", palette=pal[c(1, 4, 3, 2, 7)][pal_insert], n=5, style="cat")+
    tm_shape(inset[i, drop=F]) +
    tm_borders(col = "black", lw=03)
}

barplot(kmeans$kmeans_5, col=pal[c(1, 4, 3, 2, 7)], border=NA, axes=F, xlab="", ylab="", family="serif", horiz=T)
text(unname(freq(kmeans$kmeans_5)[c(1:5),2]), unname(freq(kmeans$kmeans_5)[c(1:5),1]), round(unname(freq(kmeans$kmeans_5)[c(1:5),2])*100/sum(unname(freq(kmeans$kmeans_5)[c(1:5),2])), 2), family="serif") 

strata8 <- list()
for (i in 1:4) {
  pal_insert <- sort(unique(strata8_inset[[i]][]))
  strata8[[i]] <- 
    tm_layout(frame=F, fontfamily="serif", legend.title.fontface=2)+
    tm_shape(brick(strata8_inset[[i]])) +
    tm_raster(col="kmeans_8", palette=pal[c(1, 8, 4, 2, 6, 3, 5, 7)][pal_insert], n=8, style="cat")+
    tm_shape(inset[i, drop=F]) +
    tm_borders(col = "black", lw=03)
}

barplot(kmeans$kmeans_8, col=pal[c(1, 8, 4, 2, 6, 3, 5, 7)], border=NA, axes=F, xlab="", ylab="", family="serif", horiz=T)
text(unname(freq(kmeans$kmeans_8)[c(1:8),2]), unname(freq(kmeans$kmeans_8)[c(1:8),1]), round(unname(freq(kmeans$kmeans_8)[c(1:8),2])*100/sum(unname(freq(kmeans$kmeans_8)[c(1:8),2])), 2), family="serif") 

tmap_arrange(list(strata3[[4]],strata3[[3]], strata3[[1]], strata3[[2]],
                  strata5[[4]],strata5[[3]], strata5[[1]], strata5[[2]],
                  strata8[[4]],strata8[[3]], strata8[[1]], strata8[[2]]), ncol=4, nrow=3)

###FIGURE 5###
#ECDF population per strata
population_strata_table_short <- read.csv("sample_design/data/out/sampling/population_strata_table_short.csv")
population_strata_table_short <- population_strata_table_short[3:4]
population_strata_table_short$sampling_covariate_3_kmeans <- as.factor(population_strata_table_short$sampling_covariate_3_kmeans)



boxplot(sampling_population~sampling_covariate_3_kmeans,data=population_strata_table_short, main="xxx", 
        xlab="xxx", ylab="xxx")

ggplot(population_strata_table_short, aes(x=sampling_covariate_3_kmeans, y=sampling_population)) + 
  geom_violin()

par(mfrow=c(1,1))
for (i in 1:3){
  plot(ecdf(population_strata_table_short[which(population_strata_table_short$sampling_covariate_3_kmeans==i),]$sampling_population), 
       col=pal[i], xlab="Population per cell", main="ECDF per stratum", add=ifelse(i==1, "F", "T"), family="serif", axes=T)
  legend(900, 0.2, legend=c("Stratum 1", "Stratum 2", "Stratum 3"), col=pal[1:3], lty=1, cex=0.9)
}
axis(1, at=c(1, max(population_strata_table_short[which(population_strata_table_short$sampling_covariate_3_kmeans==1),]$sampling_population)), tick=F, las=2, family="serif")
axis(2, at=c(0,1), tick=F, family="serif")

###FIGURE 6###
sampleSize <- c(1:1000)
#WECDF for sample per strata for different sample sizes
par(mfrow=c(1,3))
for (i in 1:3)  {
  strata_frame <- population_strata_table_short[which(population_strata_table_short$sampling_covariate_3_kmeans==i),]
  strata_frame$sampling_probability <- strata_frame$sampling_population/sum(strata_frame$sampling_population)
  plot(ecdf(strata_frame$sampling_population), col="black", lty=2, main=paste("Stratum", i), xlab="Population per cell")
  for (n in sampleSize){
    set.seed(1)
    #extract samples of different sizes from the selected stratum
    samp_frame <-strata_frame[sample(nrow(strata_frame), size=n, prob=strata_frame$sampling_probability, replace=F),]
    a <- ewcdf(samp_frame$sampling_population, 1/samp_frame$sampling_probability, normalise=T)
    m <- matrix(c(sort(rep(c(0, samp_frame$sampling_population), 2))[-1], 
                  sort(rep(c(0,a(samp_frame$sampling_population)), 2))[-length(samp_frame$sampling_population)*2]), 
                ncol=2, 
                byrow = F)
    lines(m, col=alpha(pal[i], 0.1))
  }
}

###FIGURE 7###
#Plot of sample size vs kolomogorov smirnov.
load("sample_design/data/out/sampling/ks_strata.RData")
for (i in 1:length(ks_list)){
  ks_list[[i]]$iteartion_mean <- rowMeans(ks_list[[i]][,c(3:ncol(ks_list[[i]]))])
  ks_list[[i]]$size_strata <- ks_list[[i]]$sample_size/nrow(population_strata_table_short[which(population_strata_table_short$sampling_covariate_3_kmeans==i),])
}

for (i in 1:length(ks_list)){
  ks_list[[i]]$improvement <- (ks_list[[i]]$iteartion_mean)*100/max(ks_list[[i]]$iteartion_mean)
}
par(mfrow=c(1,1))

plot(ks_list[[1]]$sample_size, ks_list[[3]]$iteartion_mean, type="l", col=pal[3], main="Kolmogorov–Smirnov distance reduction [%]", xlab="Sample Size", ylab="", axes=T, family="serif")
lines(ks_list[[2]]$sample_size, ks_list[[2]]$iteartion_mean, type="l", col=pal[2])
lines(ks_list[[3]]$sample_size, ks_list[[1]]$iteartion_mean, type="l", col=pal[1])

abline(h=0.1, lty=2)
abline(h=0.2, lty=2)

axis(1, at=c(0,
             ks_list[[3]][which(abs(ks_list[[3]]$iteartion_mean-0.15)==min(abs(ks_list[[3]]$iteartion_mean-0.15))),]$sample_size,
             ks_list[[1]][which(abs(ks_list[[1]]$iteartion_mean-0.15)==min(abs(ks_list[[1]]$iteartion_mean-0.15))),]$sample_size,
             ks_list[[2]][which(abs(ks_list[[2]]$iteartion_mean-0.15)==min(abs(ks_list[[2]]$iteartion_mean-0.15))),]$sample_size,
             500,1000), 
     label=c(0,
             ks_list[[3]][which(abs(ks_list[[3]]$iteartion_mean-0.15)==min(abs(ks_list[[3]]$iteartion_mean-0.15))),]$sample_size,
             ks_list[[1]][which(abs(ks_list[[1]]$iteartion_mean-0.15)==min(abs(ks_list[[1]]$iteartion_mean-0.15))),]$sample_size,
             ks_list[[2]][which(abs(ks_list[[2]]$iteartion_mean-0.15)==min(abs(ks_list[[2]]$iteartion_mean-0.15))),]$sample_size,
             500,1000), tick=F, family="serif")
axis(2, at=c(round(min(ks_list[[3]]$iteartion_mean), 2), 0.10, 0.15, 0.20, round(max(ks_list[[1]]$iteartion_mean), 2)), labels=c(round(min(ks_list[[3]]$iteartion_mean), 2), 0.1, 0.15, 0.2, round(max(ks_list[[1]]$iteartion_mean),2)), tick=F, las=2, family="serif")


points(ks_list[[1]][which(abs(ks_list[[1]]$iteartion_mean-0.15)==min(abs(ks_list[[1]]$iteartion_mean-0.15))),]$sample_size, 0.15, pch=1, cex=1.5)
points(ks_list[[2]][which(abs(ks_list[[2]]$iteartion_mean-0.15)==min(abs(ks_list[[2]]$iteartion_mean-0.15))),]$sample_size, 0.15, pch=1, cex=1.5)
points(ks_list[[3]][which(abs(ks_list[[3]]$iteartion_mean-0.15)==min(abs(ks_list[[3]]$iteartion_mean-0.15))),]$sample_size, 0.15, pch=1, cex=1.5)

###GRID SMAPLING
#load kmeans 1:16
kmeans_list <- list.files("sample_design/data/out/kmeans/", pattern=".tif$")
kmeans <- stack(paste0("sample_design/data/out/kmeans/" , kmeans_list))
#kmeans <- mask(kmeans, as(provinces, "Spatial"))
kmeans_3 <- layerize(kmeans$kmeans_3)
NAvalue(kmeans_3) <- 0
urban <- reclassify(kmeans_3, c(0,4,1))
NAvalue(urban) <- 0
pop <- raster("sample_design/data/in/cod_2prov_100m_worldpop_2016_ppp_v2.tif")

kmeans_3_pop <- mask(pop, kmeans_3)
sampleSize <- c(139, 83, 171)

sampleSize_tot <- sum(139, 83, 171)

sampleSize_perc <- 116.713/sampleSize_tot*sampleSize
  
sample_strata <- list()

for (i in 1:3){
sample_strata[i] <- gs_sample(population_raster=kmeans_3_pop[[i]], 
                              strata_raster=kmeans_3[[i]], 
                              urban_raster=urban[[i]], 
                              cfg_hh_per_stratum=sampleSize[i],
                              cfg_hh_per_urban=1,
                              cfg_hh_per_rural=1,
                              cfg_pop_per_psu=0,
                              cfg_sample_rururb=F,
                              cfg_sample_spatial=F,
                              cfg_desired_cell_size=1,
                              cfg_max_psu_size=0.1,
                              cfg_min_pop_per_cell=0.0001,
                              cfg_psu_growth=F,
                              cfg_random_number=1,
                              output_path=paste0("sample_design/data/out/psu_",i,"_strata.shp"),
                              sample_name=paste0("psu_",i,"_strata")
                              )

}

sample_strata[1]



for (i in 1:3){
  sample_strata[[i]]@data$probStrata <- nrow(sample_strata[[1]]@data)/sample_strata[[i]]@data$str_cells
  sample_strata[[i]]@data$probPop <- sample_strata[[i]]@data$psu_pop/sample_strata[[i]]@data$str_pop
  sample_strata[[i]]@data$weight <- (1/sample_strata[[i]]@data$probStrata)*(1/sample_strata[[i]]@data$probPop)
  sample_strata[[i]]@data$strata <- as.numeric(i)
}

sample_strata_all <- raster::union( sample_strata[[1]],  sample_strata[[2]])
sample_strata_all <- raster::union(sample_strata_all, sample_strata[[3]])

pal <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf")

tmap_options(max.raster = c(plot = 14879410, view = 14879410))
tm_shape(brick(settled_area)) +
  tm_raster(col="cod_2prov_100m_ornl_2016_settlmentbinary", pal="gray50", legend.show = FALSE) +
  tm_shape(provinces, is.master=T, unit="km") +
  tm_borders(col="black", lwd=0.5) +
  tm_shape(sample_strata_all)+
  tm_bubbles("weight", "strata", border.col = "white", border.lwd=0.25, palette=pal[c(1,3,2)]) +
    tm_scale_bar(breaks=50, color.dark = "black", color.light = "black", position=c("left", "bottom"))+
  tm_layout(frame = F)



  i, "_strata"))
  plot(psu_polygons[[i]])
  plot(strata_R[[i]], add = T)
