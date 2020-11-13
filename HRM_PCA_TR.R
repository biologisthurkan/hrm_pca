###########################################################################
###########################################################################
###                                                                     ###
###                           HRM PCA SCRIPT:                           ###
###                      GELISTIREN BY HAKAN DUMAN                      ###
###                          HKND1977@GMAIL.COM                         ###
###                        DUZENLEYEN KAAN HÜRKAN                       ###
###                      KAAN.HURKAN@IGDIR.EDU.TR                       ###
###                          IGDIR ÜNIVERSITESI                         ###
###                    TARIMSAL BIYOTEKNOLOJI BOLUMU                    ###
###                            5 TEMMUZ 2020                            ###
###                                                                     ###
###########################################################################
###########################################################################

###########################################################################
###                          Gerekli Paketler                           ###
###########################################################################
library(tidyverse) # veri duzenlemesi icin
library(psych) # pca icin
library(plot3D) # 3d grafikler
library(factoextra) #  clustering grafikleri icin
library(mclust) # kumeleme analizi
#- yuklu olmayan paketler icin install.packages('paket ismi') kullanilir  -


###########################################################################
###                         Verinin yuklenmesi                          ###
###########################################################################
#------------ read.csv komutuna dosya yoluyla berbaber yazilir ------------
veri <- read.csv(file = "data/primer12.csv",header = TRUE)
#---------- verisetinin yapisal durumu icin str komutu kullanilir ---------
str(veri)
bagimliDegisken <- veri[,1]

###########################################################################
###                            Erime bolgesi                            ###
###########################################################################
maxLim <- 80.85 #ust limit
minLim <- 76.86 #alt limit
#---------------------------- Filtreleme islemi ---------------------------
normveri  <- as_tibble(veri) %>% 
  filter(Temperature<maxLim,Temperature>minLim)

###########################################################################
###                 Her bir degiskeni yuzdelige cevirme                 ###
###########################################################################
#---- apply verilen fonksiyonu butun satirlara veya sutunlara uygular ----
#------------ eger 2. parametre 1 ise satirlar,  2 ise sutunlar -----------
normveri <- apply(normveri[,-1],2,function(x) x/max(x))


###########################################################################
###                            Medyan egrisi                            ###
###########################################################################
medcurve <- apply(normveri, 1, median)

###########################################################################
### Medyan egrileri yardimiyla yuzde egrisinin normallesitirlmesi ###
###########################################################################
for(i in 1:ncol(normveri)){
  normveri[,i]=normveri[,i]-medcurve}
head(normveri)

###########################################################################
###                     Temel Bilesen Analizi (PCA)                     ###
###########################################################################
fitPCA <-prcomp(normveri, scale = FALSE) 
plot(fitPCA, type="line",  main="PCs of HRM") # Pareto plot
#-------- anlamli olan PC 1:3 olarak secilip yeni veri olusturuldu --------
veriPCA <- as.data.frame(fitPCA$rotation[,1:3])

###########################################################################
###                              Kumeleme                               ###
###########################################################################
fitCluster <- Mclust(veriPCA)
#---------------------------- En iyi model adi ----------------------------
fitCluster$modelName
#-------------------------- En uygun kume sayisi --------------------------
fitCluster$G
#--------------------- Kume olasiliklari (ilk 6 satir) ---------------------
head(fitCluster$z)
#-------------------------------- Kumeler ---------------------------------
fitCluster$classification

###########################################################################
###                         Kumeleme Gorselleri                         ###
###########################################################################

#--- BIC'ye gore en uygun kume sayilarinin belirlenmesini gosteren grafik ---
fviz_mclust(fitCluster, "BIC", palette = "jco") 
#---------------------------- Kumeleme Grafigi ----------------------------
fviz_mclust(fitCluster, 
            "classification", 
            geom = "point", 
            ellipse.type = "t",
            pointsize = 1.5, 
            palette = "jco") 
#-------------------------- Belirsizlik Grafigi ---------------------------
fviz_mclust(fitCluster, 
            "uncertainty", 
            palette = "jco")

###########################################################################
###                         ggplot2 grafikleri                          ###
###########################################################################

veriPCA <- as_tibble(veriPCA) 
veriPCA <- veriPCA %>% 
  mutate(Genotypes = factor(fitCluster$classification))

veriPCA %>% ggplot(aes(x=PC1,y=PC2,color=Genotypes,shape=Genotypes))+
  geom_point(size = 5) + theme_classic(base_size = 20)

veriPCA %>% ggplot(aes(x=PC2,y=PC3,color=Genotypes,shape=Genotypes))+
  geom_point(size = 5) + theme_classic(base_size = 20)

veriPCA %>% ggplot(aes(x=PC1,y=PC3,color=Genotypes,shape=Genotypes))+
  geom_point(size = 5) + theme_classic(base_size = 20)

#-------------------------------- 3D grafik -------------------------------
points3D(x=veriPCA$PC1,y=veriPCA$PC2,z=veriPCA$PC3, 
         colvar=fitCluster$classification, 
         xlab = 'PC1', ylab = 'PC2', zlab = 'PC3', pch = 16,
         colkey = F, col=c('black', 'red', 'green')
)

identify(veriPCA[,1:3], labels = rownames(veriPCA))

###########################################################################
###                   Grafiklerin Disari Aktarilmasi                    ###
###########################################################################
#------ pdf yerine jpeg,png gibi farkli dosya turleri de secilebilir -------
pdf("PC1 vs PC2.pdf")
plot(veriPCA[,1:2], bg=fitCluster$classification,
     pch=21,xlab='PC1', ylab='PC2')
dev.off()

pdf("PC2 vs PC3.pdf")
plot(veriPCA[,2:3], bg=fitCluster$classification,
     pch=21,xlab='PC2', ylab='PC3')
dev.off()

pdf("PC1 vs PC3.pdf")
plot(veriPCA[,c(1,3)], bg=fitCluster$classification,
     pch=21,xlab='PC1', ylab='PC3')
dev.off()

pdf("3D.pdf")
points3D(x=veriPCA$PC1,y=veriPCA$PC2,z=veriPCA$PC3, 
         colvar=fitCluster$classification, 
         xlab = 'PC1', ylab = 'PC2', zlab = 'PC3', pch = 16,
         colkey = F, col=c('black', 'red', 'green')
)
dev.off()

###########################################################################
###                 Kumeleme Verisinin Disari Alinmasi                  ###
###########################################################################
write.csv(fitCluster$classification, 
          file = "variant call.csv", sep=" ")

