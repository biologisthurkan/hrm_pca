############################################################################
############################################################################
###                                                                      ###
###                            HRM PCA SCRIPT:                           ###
###                       DEVELOPED BY HAKAN DUMAN                       ###
###                          HKND1977@GMAIL.COM                          ###
###                        EDITED BY KAAN HÜRKAN                         ###
###                       KAAN.HURKAN@IGDIR.EDU.TR                       ###
###                           IGDIR UNIVERSITY                           ###
###               DEPARTMENT OF AGRICULTURAL BIOTECHNOLOGY               ###
###                              5 JULY 2020                             ###
###                                                                      ###
############################################################################
############################################################################

###########################################################################
###                          Dependent Packages                         ###
###########################################################################
library(tidyverse) # data editing
library(psych) # pca
library(plot3D) # 3d graphics
library(factoextra) #  clustering graphics
library(mclust) # clustering analysis
#- use install.packages('packet name') to install packages  -


###########################################################################
###                                Data Import                          ###
###########################################################################
#------------ read.csv with file directory ------------
data <- read.csv(file = "data/primer1.csv",header = TRUE)
#---------- identify data as string ---------
str(data)
dependentVariable <- data[,1]

###########################################################################
###                               HRM Range                             ###
###########################################################################
maxLim <- 80.85 #upper limit
minLim <- 76.86 #lower limit
#-------------------------- Filtering the temperature -------------------------
normdata  <- as_tibble(data) %>% 
  filter(Temperature<maxLim,Temperature>minLim)

###########################################################################
###                   Turning the Variables to Percentage                   ###
###########################################################################
normdata <- apply(normdata[,-1],2,function(x) x/max(x))


###########################################################################
###                            Median Curve                            ###
###########################################################################
medcurve <- apply(normdata, 1, median)

###########################################################################
###                 Normalisation of the Melting Curve                  ###
###########################################################################
for(i in 1:ncol(normdata)){
  normdata[,i]=normdata[,i]-medcurve}
head(normdata)

###########################################################################
###                 Principal Component Analysis (PCA)                  ###
###########################################################################
fitPCA <-prcomp(normdata, scale = FALSE) 
plot(fitPCA, type="line",  main="PCs of HRM") # Pareto plot

dataPCA <- as.data.frame(fitPCA$rotation[,1:3])

###########################################################################
###                              Clustering                               ###
###########################################################################
fitCluster <- Mclust(dataPCA)
#----------------------------  Best model name ----------------------------
fitCluster$modelName
#--------------------------  Best cluster count  --------------------------
fitCluster$G
#-------------------- Cluster probability(first 6 row) --------------------
head(fitCluster$z)
#-------------------------------- Clusters ---------------------------------
fitCluster$classification

###########################################################################
###                      Clustering Visualisation                       ###
###########################################################################

fviz_mclust(fitCluster, "BIC", palette = "jco") 
#--------------------------- Clustering graphic ---------------------------
fviz_mclust(fitCluster, 
            "classification", 
            geom = "point", 
            ellipse.type = "t",
            pointsize = 1.5, 
            palette = "jco") 
#-------------------------- Uncertainty graphic----------------------------
fviz_mclust(fitCluster, 
            "uncertainty", 
            palette = "jco")

###########################################################################
###                          ggplot2 Graphics                           ###
###########################################################################

dataPCA <- as_tibble(dataPCA) 
dataPCA <- dataPCA %>% 
  mutate(Genotypes = factor(fitCluster$classification))

dataPCA %>% ggplot(aes(x=PC1,y=PC2,color=Genotypes,shape=Genotypes))+
  geom_point(size = 5) + theme_classic(base_size = 20)

dataPCA %>% ggplot(aes(x=PC2,y=PC3,color=Genotypes,shape=Genotypes))+
  geom_point(size = 5) + theme_classic(base_size = 20)

dataPCA %>% ggplot(aes(x=PC1,y=PC3,color=Genotypes,shape=Genotypes))+
  geom_point(size = 5) + theme_classic(base_size = 20)

#-------------------------------- 3D Graphic ------------------------------
points3D(x=dataPCA$PC1,y=dataPCA$PC2,z=dataPCA$PC3, 
         colvar=fitCluster$classification, 
         xlab = 'PC1', ylab = 'PC2', zlab = 'PC3', pch = 16,
         colkey = F, col=c('black', 'red', 'green')
)

identify(dataPCA[,1:3], labels = rownames(dataPCA))

###########################################################################
###                       Exporting the Graphics                        ###
###########################################################################
#------ pdf,jpeg,png formats are available -------
pdf("PC1 vs PC2.pdf")
plot(dataPCA[,1:2], bg=fitCluster$classification,
     pch=21,xlab='PC1', ylab='PC2')
dev.off()

pdf("PC2 vs PC3.pdf")
plot(dataPCA[,2:3], bg=fitCluster$classification,
     pch=21,xlab='PC2', ylab='PC3')
dev.off()

pdf("PC1 vs PC3.pdf")
plot(dataPCA[,c(1,3)], bg=fitCluster$classification,
     pch=21,xlab='PC1', ylab='PC3')
dev.off()

pdf("3D.pdf")
points3D(x=dataPCA$PC1,y=dataPCA$PC2,z=dataPCA$PC3, 
         colvar=fitCluster$classification, 
         xlab = 'PC1', ylab = 'PC2', zlab = 'PC3', pch = 16,
         colkey = F, col=c('black', 'red', 'green')
)
dev.off()

###########################################################################
###                   Exporting the Clustering Data                     ###
###########################################################################
write.csv(fitCluster$classification, 
          file = "variant call.csv", sep=" ")

