###########################################################################
###########################################################################
###                                                                     ###
###                           HRM PCA SCRIPT:                           ###
###                      GELI�TIREN BY HAKAN DUMAN                      ###
###                          HKND1977@GMAIL.COM                         ###
###                        D�ZENLEYEN KAAN H�RKAN                       ###
###                      KAAN.HURKAN@IGDIR.EDU.TR                       ###
###                          I�D�R �NIVERSITESI                         ###
###                    TAR�MSAL BIYOTEKNOLOJI B�L�M�                    ###
###                            5 TEMMUZ 2020                            ###
###                                                                     ###
###########################################################################
###########################################################################

###########################################################################
###                          Gerekli Paketler                           ###
###########################################################################
library(tidyverse) # veri düzenlemesi için
library(psych) # pca için
library(plot3D) # 3d grafikler
library(factoextra) #  clustering grafikleri için
library(mclust) # kümeleme analizi
#- yüklü olmayan paketler için install.packages('paket ismi') kullanılır  -


###########################################################################
###                         Verinin yüklenmesi                          ###
###########################################################################
#------------ read.csv komutuna dosya yoluyla berbaber yazılır ------------
veri <- read.csv(file = "data/primer12.csv",header = TRUE)
#---------- verisetinin yapısal durumu için str komutu kullanılır ---------
str(veri)
bagimliDegisken <- veri[,1]

###########################################################################
###                            Erime bölgesi                            ###
###########################################################################
maxLim <- 80.85 #üst limit
minLim <- 76.86 #alt limit
#---------------------------- Filtreleme işlemi ---------------------------
normveri  <- as_tibble(veri) %>% 
  filter(Temperature<maxLim,Temperature>minLim)

###########################################################################
###                 Her bir Değişkeni Yüzdeliğe çevirme                 ###
###########################################################################
#---- apply verilen fonksiyonu bütün satırtlara veya sütunlara uygular ----
#------------ eğer 2. parametre 1 ise satırlar,  2 ise sütunlar -----------
normveri <- apply(normveri[,-1],2,function(x) x/max(x))


###########################################################################
###                            Medyan eğrisi                            ###
###########################################################################
medcurve <- apply(normveri, 1, median)

###########################################################################
###Medyan eğrileri yardımıyla yüzde erime eğirisinin normalleştirilmesi ###
###########################################################################
for(i in 1:ncol(normveri)){
  normveri[,i]=normveri[,i]-medcurve}
head(normveri)

###########################################################################
###                 Principal component analysis (PCA)                  ###
###########################################################################
fitPCA <-prcomp(normveri, scale = FALSE) 
plot(fitPCA, type="line",  main="PCs of HRM") # Pareto plot
#-------- anlamlı olan PC 1:3 olarak seçilip yeni veri oluşturuldu --------
veriPCA <- as.data.frame(fitPCA$rotation[,1:3])

###########################################################################
###                              Kümeleme                               ###
###########################################################################
fitCluster <- Mclust(veriPCA)
#---------------------------- En iyi model adı ----------------------------
fitCluster$modelName
#-------------------------- En uygun küme sayısı --------------------------
fitCluster$G
#--------------------- Küme olasılıkları(ilk 6 satır) ---------------------
head(fitCluster$z)
#-------------------------------- Kümeler ---------------------------------
fitCluster$classification

###########################################################################
###                         Kümeleme Görselleri                         ###
###########################################################################

#--- BIC'ye göre en uygun küme sayısının belirlenmesini gösteren grafik ---
fviz_mclust(fitCluster, "BIC", palette = "jco") 
#---------------------------- Kümeleme Grafiği ----------------------------
fviz_mclust(fitCluster, 
            "classification", 
            geom = "point", 
            ellipse.type = "t",
            pointsize = 1.5, 
            palette = "jco") 
#-------------------------- uncertainty grafiği ---------------------------
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
###                   Grafiklerin dışarı aktarılması                    ###
###########################################################################
#------ pdf yerine jpeg,png gibi farklı dosya türleride seçilebilir -------
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
###                 Kümeleme Verisinin dışarı alınması                  ###
###########################################################################
write.csv(fitCluster$classification, 
          file = "variant call.csv", sep=" ")

