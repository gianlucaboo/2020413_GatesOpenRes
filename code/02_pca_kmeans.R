####SET UP THE WORKING ENVIRONMENT####
gc()
#load libraries
library(rdrop2); library(raster); library(factoextra); library(tcltk); library(cluster); library(clValid); library(readr)

#select working environment
setwd(tclvalue(tkchooseDirectory()))
getwd()

dir.create("sample_design/data/out/kmeans")

####LOAD AND PROCESS DATA####
#load the pca data
covariate_pca_list <- list.files("sample_design/data/out/pca", pattern=".tif$")

#read in the raster data
covariate_pca <- stack(paste0("sample_design/data/out/pca/" , covariate_pca_list))
rm(covariate_pca_list)

#Extract standardized covariate values
covariate_table_long_PCA_scores <- getValues(covariate_pca)
rm(covariate_pca)
cell_id <- 1:nrow(covariate_table_long_PCA_scores)
covariate_table_long_PCA_scores <- cbind(cell_id, covariate_table_long_PCA_scores)

#Discard NA values
covariate_table_short_PCA_scores <- data.frame(covariate_table_long_PCA_scores[complete.cases(covariate_table_long_PCA_scores),])
covariate_table_short_kmeans <- covariate_table_short_PCA_scores[,-1]
rm(covariate_table_long_PCA_scores)

####K MEANS####
covariate_table_short_kmeans_cluster <- list()
for (i in 1:16){
  set.seed(1)
  centers_selected <- sample(1:nrow(covariate_table_short_kmeans), i)
  covariate_table_short_kmeans_cluster[[i]] <- kmeans(covariate_table_short_kmeans, centers=covariate_table_short_kmeans[centers_selected,], nstart = 20, iter.max = 200, algorithm="Lloyd")
}

#extract cluster classes for each cluster number
covariate_table_short_kmeans_cluster_classes <- list()
for (i in 1:length(covariate_table_short_kmeans_cluster)){
  covariate_table_short_kmeans_cluster_classes[[i]] <- covariate_table_short_kmeans_cluster[[i]]$cluster
}
covariate_table_short_kmeans_cluster_classes <- as.data.frame(t(matrix(unlist(covariate_table_short_kmeans_cluster_classes), ncol=length(covariate_table_short_kmeans_cluster_classes[[1]]), byrow =TRUE)))
covariate_table_short_kmeans_cluster_classes <- cbind(covariate_table_short_PCA_scores[,1], covariate_table_short_kmeans_cluster_classes)
names(covariate_table_short_kmeans_cluster_classes) <- c("cell_id", paste0("kmeans_", 1:(ncol(covariate_table_short_kmeans_cluster_classes)-1)))

#map cluster cluster classes for each cluster number
covariate_table_long_kmeans_cluster_classes <- merge(data.frame(cell_id), covariate_table_short_kmeans_cluster_classes, by="cell_id", all.x=T, sort=T)
covariate_table_long_kmeans_cluster_classes <- covariate_table_long_kmeans_cluster_classes[,-1]
rm(covariate_table_short_kmeans_cluster_classes, covariate_table_short_PCA_scores)

ornlSettlement <- raster("sample_design/data/in/cod_2prov_100m_acled_2016_distanceconflict_pca_num.tif")
covariate_kmeans <- stack(replicate(ncol(covariate_table_long_kmeans_cluster_classes), ornlSettlement))
rm(ornlSettlement)

for (i in 1:ncol(covariate_table_long_kmeans_cluster_classes)) {
  values(covariate_kmeans[[i]]) <- as.vector(unlist(covariate_table_long_kmeans_cluster_classes[,i]))   
}
names(covariate_kmeans) <- names(covariate_table_long_kmeans_cluster_classes)

###WRITE RASTER###
writeRaster(covariate_kmeans, filename=paste0("sample_design/data/out/kmeans/", names(covariate_kmeans)), bylayer=TRUE,format="GTiff", overwrite=TRUE)

#cluster assessment
wssplot <- function(data, nc=16, seed=1){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 1:nc){
    set.seed(seed)
    centers_selected <- sample(1:nrow(covariate_table_short_kmeans), i)
    wss[i] <- sum(kmeans(covariate_table_short_kmeans, centers=covariate_table_short_kmeans[centers_selected,], nstart = 20, iter.max = 200, algorithm="Lloyd")$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of clusters",
       ylab="Within groups sum of squares")}

dev.off()
setEPS()
postscript("sample_design/graphic/covariate_kmeans.eps", family="Helvetica")
wssplot(covariate_table_short_kmeans)
plot()
dev.off()

write.csv(covariate_table_short_kmeans, "sample_design/data/out/kmeans/covariate_table_short_kmeans.csv")
rm(covariate_table_short_kmeans, covariate_table_long_kmeans_cluster_classes, covariate_table_short_kmeans_cluster, cell_id, wssplot, centers_selected)
