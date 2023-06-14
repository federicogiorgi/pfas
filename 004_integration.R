# This is the code for the integration of signature matrix and Figure 4, followed by the code for Figure 5 and 
# Supplementary figure S5.

library(openxlsx)
library(corto)
library(ggplot2)
library(ggrepel)
library(dplyr)

signatureMatrix=read.xlsx('SupplementaryFileS3.xlsx',rowNames=T)


## Step 1: integrate the signature scores for each species through stouffer method in order to
## overcome the uneven representation of genes

# M. musculus 
integration.mouse=as.data.frame(apply(signatureMatrix[,1:13],1,function(x){
  mat=x[!is.na(x)]
  stouffer(mat)
}))

# H.sapiens 
integration.human=as.data.frame(apply(signatureMatrix[,14:68],1,function(x){
  mat=x[!is.na(x)]
  stouffer(mat)
}))

# C. elegans 
integration.Celegans=as.data.frame(apply(signatureMatrix[,69:82],1,function(x){
  mat=x[!is.na(x)]
  stouffer(mat)
}))

# D. rerio 
integration.zebrafish=as.data.frame(apply(signatureMatrix[,83:84],1,function(x){
  mat=x[!is.na(x)]
  stouffer(mat)
}))

# G. morhua
integration.Gmorhua=as.data.frame(apply(signatureMatrix[,85:96],1,function(x){
  mat=x[!is.na(x)]
  stouffer(mat)
}))

# M. salmoides
integration.Msalmoides=as.data.frame(apply(signatureMatrix[,97:104],1,function(x){
  mat=x[!is.na(x)]
  stouffer(mat)
}))

# P. promelas
integration.Ppromelas=as.data.frame(apply(signatureMatrix[,105:110],1,function(x){
  mat=x[!is.na(x)]
  stouffer(mat)
}))


## Step 2: create a matrix containing stouffer coefficients of 7 species

mat=cbind(integration.mouse,integration.human,integration.Celegans,integration.zebrafish,
                  integration.Gmorhua,integration.Msalmoides,integration.Ppromelas)
colnames(mat)=c('integration.mouse','integration.human','integration.Celegans','integration.zebrafish',
                  'integration.Gmorhua','integration.Msalmoides','integration.Ppromelas')


## Step 3: retention of genes present in almost 6 species

mat[mat=='NaN']=NA
mat=mat[rowSums(is.na(mat))<2,]


## Step 4: integrate the stouffer coefficients of 7 species through stouffer method

stouffer_coefficient=as.data.frame(apply(mat,1,function(x){
  mat=x[!is.na(x)]
  stouffer(mat)
}))


## Step 5: calculate sandard deviation on stouffer coefficients of 7 species 

sd=as.data.frame(apply(mat,1,function(x){
  mat=x[!is.na(x)]
  sd(mat)
}))


# Step 6: create the matrix with integrated scores and select the genes basing on stouffer coefficient and sd

matIntegration=cbind(stouffer_coefficient,sd)
colnames(matIntegration)=c('stouffer_coefficient','sd')

matIntegration$genes=rownames(matIntegration)

matIntegration$expression='NS'
matIntegration$expression[matIntegration$sd<10 & matIntegration$stouffer_coefficient<(-5)]=
  'Downregulation & Low Sd'
matIntegration$expression[matIntegration$sd<10 & matIntegration$stouffer_coefficient>5]=
  'Upregulation & Low Sd'
matIntegration$expression[matIntegration$sd>10 & matIntegration$stouffer_coefficient<(-10)]=
  'Downregulation & High Sd'
matIntegration$expression[matIntegration$sd>10 & matIntegration$stouffer_coefficient>10]=
  'Upregulation & High Sd'

matIntegration$label=NA
matIntegration$label[matIntegration$expression!='NS']=matIntegration$genes[matIntegration$expression!='NS']


# Step 7: Figure 4 - plot with the integrated scores across datasets

png('Figure 4.png',w=5000,h=3000,res=420)
ggplot(data=matIntegration,aes(stouffer_coefficient,sd,color=expression,label=label)) + 
  geom_point(shape=19,size=4) +
  scale_colour_manual(values=c('Downregulation & Low Sd'='royalblue',
                               'Upregulation & Low Sd'='firebrick1',
                               'Downregulation & High Sd'='darkturquoise',
                               'Upregulation & High Sd'='orange'),
                      breaks=c('Downregulation & Low Sd','Upregulation & Low Sd',
                               'Downregulation & High Sd','Upregulation & High Sd')) +
  geom_text_repel(box.padding=0.5,max.overlaps=Inf,colour='black') + 
  labs(title='Integrated signatures',color='Expression') +
  xlab(label='Integrated signature score across datasets') + 
  ylab(label='Standard deviation') + 
  theme_bw()
dev.off()


# Step 8: Figure 5 - line plot with expression levels of genes showing high stouffer coefficient and low sd

up=filter(matIntegration,expression=='Upregulation & Low Sd')
dn=filter(matIntegration,expression=='Downregulation & Low Sd')

signUp=filter(signatureMatrix,rownames(signatureMatrix)%in%rownames(up))
signDn=filter(signatureMatrix,rownames(signatureMatrix)%in%rownames(dn))
genes=rbind(signUp,signDn)

genes[is.na(genes)]=0


col1=c('#7D0025','#6D0026','#8B0020','#990316','#9E1021','#9F0032','#8E063B','#AC1425','#B71632',
       '#C91837','#A70B00','#B41500','#C01F00','#CC2900','#F5191C','#E42237','#FA4E5A','#CF2E2D',
       '#E23D00','#ED6132','#FF6E5A','#EC5A00','#E88600','#EDA92B','#FFDC00')
col2=c('#0D1E34','#10253F','#142D4A','#1C3C61','#00366C','#023FA5','#244B7A','#00467D','#1C518F',
       '#2C5B92','#005691','#344FA5','#356CAC','#4C5FA9','#376CC3','#0067A7','#0F71AF','#0078BD',
       '#5F6DAE','#707CB4','#3488CB','#0084A7','#00A0BB')
 
png('Figure 5.png',w=5000,h=3000,res=420)
matplot(t(genes),type='l',col=c(col1,col2),lty=1,lwd=2,xlab='Contrasts',ylab='Signature',xaxt='n')
abline(h=c(1.30103,-1.30103),col='gray22',lty=5)
axis(1,at=c(7,42,76,84,91,101,107),labels=c('M.musculus','H.sapiens','C.elegans','D.rerio',
                                            'G.morhua','M. salmoides','P. promelas'),
     las=2,font=3,cex.axis=1,gap.axis=0.1,mgp=c(3,0.2,0),tick=F)
abline(v=c(1,14,69,83,85,97,105),col='black',lty=3)
legend('topright',legend=c(rownames(genes)),lty=1,lwd=2,col=c(col1,col2),cex=0.5,text.font=3,ncol=6)
dev.off() 


# Step 9: Supplementary Figure S5 - line plot with expression levels of outliers genes 
# showing high stouffer coefficient and high sd

up2=filter(matIntegration,expression=='Upregulation & High Sd')
dn2=filter(matIntegration,expression=='Downregulation & High Sd')

signUp2=filter(signatureMatrix,rownames(signatureMatrix)%in%rownames(up2))
signDn2=filter(signatureMatrix,rownames(signatureMatrix)%in%rownames(dn2))
genes2=rbind(signUp2,signDn2)

genes2[is.na(genes2)]=0


col3=c('#CC5522','#D75C20','#E26320','#EA6D20','#EF7924','#F2862C','#F49335','#F5A040','#F6AC50')
col4=c('#006064','#00838F','#0097A7','#00ACC1','#00BCD4','#26C6DA','#4DD0E1','#80DEEA')

png('Supplementary Figure S5.png',w=5000,h=3000,res=420)
matplot(t(genes2),type='l',col=c(col3,col4),lty=1,lwd=2,ylab='Signature',xlab='Contrasts',xaxt='n')
abline(h=c(1.30103,-1.30103),col='gray22',lty=5)
axis(1,at=c(7,42,76,84,91,101,107),labels=c('M.musculus','H.sapiens','C.elegans','D.rerio',
                                            'G.morhua','M. salmoides','P. promelas'),
     las=2,font=3,cex.axis=1,gap.axis=0.1,mgp=c(3,0.2,0),tick=F)
abline(v=c(1,14,69,83,85,97,105),col='black',lty=3)
legend('topright',legend=c(rownames(genes2)),lty=1,lwd=2,col=c(col3,col4),cex=0.7,text.font=3,ncol=4)
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
#   [1] dplyr_1.1.0      ggrepel_0.9.3    ggplot2_3.4.1    corto_1.2.0      openxlsx_4.2.5.2

# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.10        rstudioapi_0.14    magrittr_2.0.3     munsell_0.5.0      tidyselect_1.2.0   colorspace_2.1-0  
# [7] R6_2.5.1           rlang_1.0.6        pbapply_1.7-0      fansi_1.0.4        caTools_1.18.2     tools_4.2.2       
# [13] grid_4.2.2         parallel_4.2.2     gtable_0.3.1       plotrix_3.8-2      KernSmooth_2.23-20 utf8_1.2.2        
# [19] cli_3.6.0          withr_2.5.0        gtools_3.9.4       tibble_3.1.8       lifecycle_1.0.3    zip_2.2.2         
# [25] farver_2.1.1       vctrs_0.5.2        bitops_1.0-7       glue_1.6.2         labeling_0.4.2     stringi_1.7.12    
# [31] compiler_4.2.2     pillar_1.8.1       scales_1.2.1       generics_0.1.3     gplots_3.1.3       pkgconfig_2.0.3       