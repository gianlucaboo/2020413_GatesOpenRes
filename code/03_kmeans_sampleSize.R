####SET UP THE WORKING ENVIRONMENT####
gc()
#load libraries
library(rdrop2); library(raster); library(tcltk); library(parallel); library(foreach); library(doParallel); 
library(ggridges); library(Hmisc); library(spatstat); library(matrixStats)


####LOAD AND PROCESS DATA####
#select working environment
setwd(tclvalue(tkchooseDirectory()))
getwd()

dir.create("sample_design/data/out/sampling")

#load the pca data
covariate_kmeans_list <- list.files("sample_design/data/out/kmeans", pattern=".tif$")

#read in the raster data
covariate_kmeans <- stack(paste0("sample_design/data/out/kmeans/", covariate_kmeans_list))
rm(covariate_kmeans_list)

#Extract WorldPop population estimates within settelement area as vector
covariate_3_kmeans <- covariate_kmeans$kmeans_3
#plot(covariate_3_kmeans)
covariate_3_kmeans_reclass <- matrix(c(c(1,2,3), c(1,3,2)), ncol = 2)
plot(covariate_3_kmeans)

names(covariate_3_kmeans_reclass) <- c("is", "becomes")
covariate_3_kmeans <- reclassify(covariate_3_kmeans, covariate_3_kmeans_reclass, right=NA)
plot(covariate_3_kmeans)

rm(covariate_kmeans, covariate_3_kmeans_reclass)
settled <- reclassify(covariate_3_kmeans, c(0, 10, 1))
population <- raster("sample_design/data/in/cod_2prov_100m_worldpop_2016_ppp_v2.tif")

population <- raster::mask(population, settled)
population_strata <- stack(population, covariate_3_kmeans)
names(population_strata) <- c("sampling_population", "sampling_covariate_3_kmeans")
rm(settled)
writeRaster(population_strata, filename=paste0("sample_design/data/out/sampling/", names(population_strata)), bylayer=TRUE,format="GTiff", overwrite=T)
rm(covariate_3_kmeans, population)

#create table
population_strata_table_long <- data.frame(getValues(population_strata))
cell_id <- 1:nrow(population_strata_table_long)
population_strata_table_long <- cbind(cell_id, population_strata_table_long)
population_strata_table_short <- population_strata_table_long[complete.cases(population_strata_table_long),]

write.csv(population_strata_table_short, "sample_design/data/out/sampling/population_strata_table_short.csv")

####VISUALIZE POPULATION DENSITIES WITHIN SETTLEMENT AREA####
#histogram of population densities within settelement area and classes breaks
boxplot(sampling_population~as.factor(sampling_covariate_3_kmeans),data=population_strata_table_short, xlab="Strata", ylab="Population per pixel")

#density plot
#ggplot(population_strata_table_short, aes(x = sampling_population, y = as.factor(sampling_covariate_3_kmeans))) + geom_density_ridges()

####ASSESS SAMPLE SIZE####
#startify the population estimates within settelement area
population_strata_table_short$sampling_covariate_3_kmeans <- as.character(population_strata_table_short$sampling_covariate_3_kmeans)

population_strata_list <- c()
for (i in sort(unique(population_strata_table_short$sampling_covariate_3_kmeans))) {
  population_strata_list[[i]] <- population_strata_table_short[grepl(paste(i), population_strata_table_short$sampling_covariate_3_kmeans),]
}

#compute propability proportional to population size
for (i in sort(unique(population_strata_table_short$sampling_covariate_3_kmeans))) {
  population_strata_list[[i]]$sampling_probability <- population_strata_list[[i]]$sampling_population/sum(population_strata_list[[i]]$sampling_population)
}

rm(population_strata_table_long, population_strata_table_short, population_strata, cell_id, i)

#smirnv kormogorv test with weighted ecdf
ks_weighted <- function(strata, populationStratified, maxSize, maxIteration){
  #create dataframe to store values
  ks_stat_frame <- data.frame(matrix(NA, nrow =as.numeric(maxSize), ncol = as.numeric(maxIteration)+2))
  names(ks_stat_frame) <- c("strata", "sample_size", paste0("iteration_", 1:maxIteration))
  
  #extract data frame for selected stratum
  true_frame <- populationStratified[[as.numeric(strata)]][,c(2,4)]
  true_frame <- true_frame[order(true_frame$sampling_population),]
  #compute the ecdf for the population
  f_true_pop <- ecdf(true_frame$sampling_population)
  
for (i in 1:maxIteration){
  set.seed(i)
  for (n in 1:maxSize){
  #extract samples of different sizes from the selected stratum
    samp_frame <-true_frame[sample(nrow(true_frame), size=n, prob=true_frame$sampling_probability, replace=F),]
    samp_frame <- samp_frame[order(samp_frame$sampling_population),]
    f_samp_pop <- ewcdf(samp_frame$sampling_population, 1/samp_frame$sampling_probability, normalise=T)
    ks_stat <- max(abs(f_true_pop(samp_frame$sampling_population) - f_samp_pop(samp_frame$sampling_population)))
    rm(samp_frame, f_samp_pop)
    ks_stat_frame[n, 1] <- as.numeric(strata)
    ks_stat_frame[n, 2] <- as.numeric(n)
    ks_stat_frame[n,2+i] <- ks_stat
  }
}
rm(f_true_pop)
ks_stat_frame
}

#parallel computation
if(!file.exists("sample_design/data/out/sampling/ks_strata.RData")[1]){
cores <- detectCores()
cl <- makeCluster(3)
registerDoParallel(cl)

ks_list <- list()
ks_list <- foreach (s = c(1:3), .packages = c("spatstat")) %dopar% {
  ks_weighted(strata=s, populationStratified=population_strata_list, maxSize=1000, maxIteration=1000)
}

save(ks_list, file="sample_design/data/out/sampling/ks_strata.RData")} else{
load("sample_design/data/out/sampling/ks_strata.RData")
rm(ks_weighted)
}


for (i in 1:length(ks_list)){
  ks_list[[i]]$iteration_mean <- rowMeans(ks_list[[i]][,c(3:ncol(ks_list[[i]]))])
  ks_list[[i]]$iteration_sd <- rowSds(as.matrix(ks_list[[i]][,c(3:ncol(ks_list[[i]]))]))
  ks_list[[i]]$performance <- ks_list[[i]]$iteration_mean / ks_list[[i]]$iteration_sd
  ks_list[[i]]$std_sd <- ks_list[[i]]$iteration_sd / ks_list[[i]]$iteration_mean
  
  }
  