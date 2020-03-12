####SET UP THE WORKING ENVIRONMENT####
gc()

#load libraries
library(rdrop2); library(raster); library(factoextra); library(PCAmixdata); library(tcltk); library(readr); library(rgdal)

#select working environment

setwd(tclvalue(tkchooseDirectory()))
getwd()

#create folder structure
dir.create("sample_design")
dir.create("sample_design/data")
dir.create("sample_design/data/in")
dir.create("sample_design/data/out")
dir.create("sample_design/data/out/masked")
dir.create("sample_design/data/out/pca")
dir.create("sample_design/graphic")

####LOAD AND PROCESS DATA####
##download the data from dropbox if not there
if (length(list.files("sample_design/data/in/", pattern=".tif")) == 0) {
data_list <- drop_dir("/sampling_paper_data", recursive=T, include_media_info=F, include_deleted=F)
data_list <- data_list[which(data_list$.tag == "file"),]

for (i in 1:nrow(data_list)){
drop_download(data_list$path_lower[i], local_path="sample_design/data/in", overwrite=F) 
}
}
#load the raster data, harmonize them, and mask them using the admin boundaries
boundaries <- readOGR("sample_design/data/in", "cod_bcr_admin_1")

raster_list <- list.files("sample_design/data/in/", pattern=".tif", recursive = TRUE)
raster <- stack(paste0("sample_design/data/in/", raster_list))
raster <- mask(raster, boundaries)
raster <- mask(raster, raster$cod_2prov_100m_ornl_2016_settlmentbinary)

#extract covariate data
covariate_list <- grepl("_pca_", raster_list)
covariate <- raster[[which(covariate_list)]]

#reclassify covariates
covariate$cod_2prov_100m_osm_2016_distanceroads_pca_num <- reclassify(covariate$cod_2prov_100m_osm_2016_distanceroads_pca_num, rcl=cbind(250, 260, NA))
covariate$cod_2prov_100m_viirs_2016_annualintens_pca_num <- reclassify(covariate$cod_2prov_100m_viirs_2016_annualintens_pca_num, rcl=cbind(60, 600, 60))

#extract numeric and categoric covariates for PCA
covariate_numeric <- dropLayer(covariate, c("cod_2prov_100m_esacci_2015_reclassified_pca_cat", "cod_2prov_100m_ghs_2015_smod_pca_cat"))
covariate_categoric <- subset(covariate, c("cod_2prov_100m_esacci_2015_reclassified_pca_cat", "cod_2prov_100m_ghs_2015_smod_pca_cat"))

covariate_categoric_dummy <- list()
for (i in 1:length(covariate_categoric@layers)){
  covariate_categoric_dummy[[i]] <- layerize(covariate_categoric[[i]], falseNA=F)
}

#combine numeric and categoric covariates
covariate_categoric <- stack(covariate_categoric_dummy)
covariate <- stack(covariate_numeric, covariate_categoric)
rm(covariate_categoric_dummy, covariate_numeric, covariate_categoric, covariate_list, raster_list)

#mask covariate by settlement layer
covariate_mask <- mask(covariate, raster$cod_2prov_100m_ornl_2016_settlmentbinary)

#drop covariate wih very few observation within the settlement layer
covariate_mask <- dropLayer(covariate_mask, c("X140", "X160", "X200", "X210"))
writeRaster(covariate_mask, filename=paste0("sample_design/data/out/masked/masked_", names(covariate_mask)), bylayer=TRUE,format="GTiff", overwrite=T)
rm(covariate)

####EXTRACT DATA AND PC ANALYSIS####
#Extract covariate values
covariate_table_long <- data.frame(getValues(covariate_mask))
rm(covariate_mask)
cell_id <- 1:nrow(covariate_table_long)
covariate_table_long <- cbind(cell_id, covariate_table_long)

#select observations and variables to keep
covariate_table_short <-  covariate_table_long[complete.cases(covariate_table_long),]
write_csv(data.frame(covariate_table_short), "sample_design/data/out/masked/masked_covariate_table_short.csv")
rm(covariate_table_long)
covariate_table_short_PCA <- covariate_table_short[,-1]

#perform the PCA
covariate_table_short_PCA <- prcomp(covariate_table_short_PCA, scale=TRUE, center=TRUE)

#cumulative explained variance
figure_PCA <- cbind(c(1:16), c(cumsum(covariate_table_short_PCA$sdev^2/sum(covariate_table_short_PCA$sdev^2)*100)))
dev.off()
setEPS()
postscript("sample_design/graphic/covaraite_PCA.eps", family="Helvetica")
plot(figure_PCA[,1], figure_PCA[,2], type = "b", ylab="Cumulative explained variance (%)", xlab="Principal Components")
points(figure_PCA[9,1], figure_PCA[9,2], col="red", cex = 1)
text(x=figure_PCA[9,1], y=figure_PCA[9,2], labels = paste(round(figure_PCA[9,2], digits=2)), pos = 2)
dev.off()
rm(figure_PCA)

#compute standard deviation of each principal component
std_dev <- covariate_table_short_PCA$sdev
pr_var <- std_dev^2

#proportion of variance explained
prop_varex <- pr_var/sum(pr_var)
sum(prop_varex[c(1:9)])

#extarct first seven independent variables
covariate_table_short_PCA_scores <- as.data.frame(covariate_table_short_PCA$x)[1:9]

dev.off()
setEPS()
postscript("sample_design/graphic/covaraite_2PCA.eps", family="Helvetica")
fviz_pca_var(covariate_table_short_PCA, col.ind = "#00AFBB",
             repel = FALSE)
dev.off()

#map cluster based on standardized covariate
covariate_table_short_PCA_scores <- cbind(covariate_table_short[,1], covariate_table_short_PCA_scores)
colnames(covariate_table_short_PCA_scores) <- c("cell_id", "PCA_1", "PCA_2", "PCA_3", "PCA_4", "PCA_5", "PCA_6", "PCA_7", "PCA_8", "PCA_9")
covariate_table_long_PCA_scores <- merge(data.frame(cell_id), covariate_table_short_PCA_scores, by="cell_id", all.x=T, sort=T, incomparables = NULL)
write.csv(covariate_table_short_PCA_scores, "sample_design/data/out/pca/covariate_table_short_pca_score.csv")
rm(covariate_table_short, covariate_table_short_PCA, covariate_table_short_PCA_scores, pr_var, prop_varex, std_dev, i)

covariate_table_long_PCA_scores <- covariate_table_long_PCA_scores[,-1]


covariate_pca <- stack(replicate(ncol(covariate_table_long_PCA_scores), raster$cod_2prov_100m_ornl_2016_settlmentbinary))
for (i in 1:ncol(covariate_table_long_PCA_scores)) {
  values(covariate_pca[[i]]) <- as.vector(unlist(covariate_table_long_PCA_scores[,i]))   
}
names(covariate_pca) <- tolower(paste0(names(covariate_table_long_PCA_scores)))

###WRITE RASTER###
writeRaster(covariate_pca, filename=paste0("sample_design/data/out/pca/", names(covariate_pca)), bylayer=TRUE,format="GTiff", overwrite=T)
rm(boundaries, cell_id, i, data_list, raster)
