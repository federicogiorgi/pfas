## This is the code for the correlation plots (Supplementary Figure S1,Supplementary Figure S2,Figure 2A).
## In order to verify if exposure to PFAS molecules stimulate correlated responses in distinct species,
## the signature matrix (Supplementary File S2) was used for all plots and was calculated the Pearson correlation.

library(openxlsx)
library(corrplot)
library(Cairo)

signatureMatrix=read.xlsx('SupplementaryFileS3.xlsx',rowNames=T)


# Supplementary Figure S1 : correlation plot of all contrast of signature matrix (110 contrasts - 7 species)

col=c(rep('darkblue',6),rep('blue',3),rep('dodgerblue2',4),rep('orangered',2),rep('deeppink',16),
      rep('#AF0040',37),rep('green4',14),rep('orange2',2),rep('darkorchid1',12),
      rep('turquoise3',8),rep('lightsteelblue4',6))


CairoPNG(width=5000,height=5000,filename='Supplementary Figure S1.png',bg='white',units='px',res=95,dpi='auto')
corrplot(cor(signatureMatrix,use='pairwise.complete.obs'),method='ellipse',type='lower',
         addCoef.col='black',tl.cex=2,tl.col=col,tl.srt=45,cl.cex=3.2,font=2,
         col=colorRampPalette(c('midnightblue',"blue3",'mediumblue','royalblue3','deepskyblue',
                                'white','orangered','firebrick1','red2',"red3",'#660000'))(1000))
legend(70,112,legend=c('M.musculus (Attema et al. 2022)','M.musculus (Heintz et al. 2022)',
                       'M.musculus (Pfohl et al. 2021)','H.sapiens (Imir et al. 2021)',
                       'H.sapiens (Rowan-Carroll et al. 2021)','H.sapiens (Reardon et al. 2021)',
                       'C.elegans (Feng et al. 2022)','D.rerio (Dasgupta et al. 2020)',
                       'G.morhua (Khan et al. 2021)','M.salmoides (Collí-Dulá et al. 2016)',
                       'P.promelas (Rodríguez-Jorquera et al. 2019)'),
       pch=15,col=c('darkblue','blue','dodgerblue2','orangered','deeppink','#AF0040','green4',
                    'orange2','darkorchid1','turquoise3','lightsteelblue4'),
                    cex=5,text.font=3,ncol=1)  
dev.off()


# Supplementary Figure S2 : correlation plot between H. sapiens and M. musculus

col=c(rep('darkblue',6),rep('blue',3),rep('dodgerblue2',4),rep('orangered',2),rep('deeppink',16),rep('#AF0040',37))

CairoPNG(width=7000,height=7000,filename='Supplementary Figure S2.png',bg='white',units='px',res=100,dpi='auto')
corrplot(cor(signatureMatrix[,1:68],use='pairwise.complete.obs'),method='ellipse',type='lower',
         addCoef.col='black',tl.cex=3.4,tl.col=col,tl.srt=45,cl.cex=5,number.cex=1.2,number.font=2,font=2,
         col=colorRampPalette(c('midnightblue',"blue3",'mediumblue','royalblue3','deepskyblue','white',
                                 'orangered','firebrick1','red2',"red3",'#660000'))(1000))
legend(45,73,legend=c('M.musculus (Attema et al. 2022)','M.musculus (Heintz et al. 2022)',
                       'M.musculus (Pfohl et al. 2021)','H.sapiens (Imir et al. 2021)',
                       'H.sapiens (Rowan-Carroll et al. 2021)','H.sapiens (Reardon et al. 2021)'),
       pch=15,col=c('darkblue','blue','dodgerblue2','orangered','deeppink','#AF0040'),cex=7,text.font=3,ncol=1)
dev.off()


# Figure 2A : correlation plot between H. sapiens and M. musculus exposed to PFOA

cor.PFOA=cor(signatureMatrix[,c(1:2,4:5,20:23,40:42,60,61)],use='pairwise.complete.obs',method='pearson')
index=which(cor.PFOA==0.36,arr.ind=T)
xleft=index[2]-0.5
ybottom=index[8]-0.5
xright=index[2]+0.5
ytop=index[8]+0.5

col=c(rep('darkblue',4),rep('deeppink',4),rep('#AF0040',5))

CairoPNG(width=4000,height=4000,filename='Figure 2A.png',bg='white',units='px',res=130,dpi='auto')
corrplot(cor.PFOA,method='ellipse',type='lower',
         addCoef.col='black',tl.cex=3.4,tl.col=col,tl.srt=45,cl.cex=3.2,number.cex=2.5,number.font=2,
         col=colorRampPalette(c('midnightblue',"blue3",'mediumblue','royalblue3','deepskyblue','white',
                                 'orangered','firebrick1','red2',"red3",'#660000'))(1000),font=2)
rect(xleft=1.5,ybottom=5.5,xright=2.5,ytop=6.5,border='black',lwd=10)                                
legend(7,16,legend=c('M.musculus (Attema et al. 2022)','H.sapiens (Rowan-Carroll et al. 2021)',
       'H.sapiens (Reardon et al. 2021)'),pch=15,
       col=c('darkblue','deeppink','#AF0040'),cex=4,text.font=3,ncol=1)
dev.off()


sessionInfo()

# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)

# Matrix products: default

# locale:
#   [1] LC_COLLATE=Italian_Italy.utf8  LC_CTYPE=Italian_Italy.utf8    LC_MONETARY=Italian_Italy.utf8
# [4] LC_NUMERIC=C                   LC_TIME=Italian_Italy.utf8    

# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] Cairo_1.6-0      corrplot_0.92    openxlsx_4.2.5.2

# loaded via a namespace (and not attached):
#   [1] compiler_4.2.2  tools_4.2.2     rstudioapi_0.14 Rcpp_1.0.10     stringi_1.7.12  zip_2.2.2
