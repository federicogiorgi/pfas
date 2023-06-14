## This is the code for Figure 2B and Supplementary Figure S3

library(openxlsx)
library(dplyr)
library(corto)
library(RColorBrewer)
library(NLP)
library(wordcloud)
library(tm)

signatureMatrix=read.xlsx('Supplementary File S2.xlsx',rowNames=T)


## Figure 2B : Scatter plot between two contrasts of mouse and human datasets, 
## Attema et al. 2022 & Rowan-Carroll et al. 2021, showing the highest Pearson correlation (0.36) in Figure 2A

# Step 1: extract the contrasts from signature matrix (Supplementary File S2) and 
# create the vectors for the scatterplot

mat=subset(signatureMatrix,select=c('Mm_AT_WT_PFOA0.3','Hs_RO_D14_PFOA20'))
mat=na.omit(mat)

mouse=setNames(mat$Mm_AT_WT_PFOA0.3,paste0(rownames(mat)))
human=setNames(mat$Hs_RO_D14_PFOA20,paste0(rownames(mat)))


# Step 2: Select genes driving the correlation between human and mouse basing on p-value threshold ≤ 0.001

thr=-log10(0.001)
mouse.up=mouse[mouse>=thr]
mouse.dn=mouse[mouse<=(-thr)]
human.up=human[human>=thr]
human.dn=human[human<=(-thr)]

up=intersect(names(human.up),names(mouse.up))
dn=intersect(names(human.dn),names(mouse.dn))
genes=c(up,dn)


# Step 3: make the scatter plot

scatter2=function (x, y, method = "pearson", threshold = 0.01, showLine = TRUE, 
                   grid = TRUE, bgcol = FALSE, pch = 16, subtitle = NULL, extendXlim = FALSE, 
                   ci = FALSE,cex=1, ...) 
{
  common <- intersect(names(x), names(y))
  x <- x[common]
  y <- y[common]
  if (!extendXlim) {
    plot(x, y, pch = pch, ...)
  }
  else {
    plot(x, y, pch = pch, xlim = 1.2 * c(min(x), max(x)), 
         ...)
  }
  cc <- cor.test(x, y, method = method)
  ccp <- signif(cc$p.value, 3)
  cccor <- signif(cc$estimate, 3)
  if (is.null(subtitle)) {
    if (ccp < 0.01) {
      vv <- format(ccp, scientific = TRUE)
      v1 <- gsub("e.+", "", vv)
      v2 <- gsub(".+e", "", vv)
      v2 <- gsub("-0+", "-", v2)
      v2 <- gsub("\\+0", "+", v2)
      v2 <- gsub("\\++", "", v2)
      bq <- as.expression(bquote("CC=" ~ .(cccor) ~ " (p=" ~ 
                                   .(v1) ~ x ~ 10^.(v2) ~ ")"))
    }
    else {
      bq <- mtext(paste0("CC=", cccor, " (p=", ccp, ")"), 
                  cex = 0.7)
    }
    mtext(bq, cex = 0.7)
  }
  else {
    mtext(subtitle, cex = 0.7)
  }
  if (bgcol) {
    if (cccor >= 0) {
      bgcol <- "#FF000033"
    }
    else {
      bgcol <- "#0000FF33"
    }
    if (ccp > threshold) {
      bgcol <- "#FFFFFF00"
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], 
         col = bgcol)
  }
  if (grid) {
    grid(col = "gray10")
  }
  if (showLine) {
    lm1 <- lm(y ~ x)
    abline(lm1$coef,col='midnightblue',lwd=2)
  }
  if (ci) {
    lm2 <- lm(y ~ x)
    newx = data.frame(x = seq(min(x), max(x), length.out = length(x)))
    confInterval = predict(lm2, newdata = data.frame(x = newx), 
                           interval = "confidence", level = 0.95)
    matlines(newx$x, confInterval[, 2:3], col = "dodgerblue1", 
             lty = 2,lwd=1)
  }
}

textplot2=function (x, y, words, cex = 1, new = TRUE, show.lines = TRUE, pointcolor='red',line.width=1,
                    line.col='gray69',rstep=0.5,padding=' ',...) 
{
  if (new) 
    plot(x, y, type = "n", ...)
  words<-paste0(padding,words,padding)
  points(x, y, pch = 16, col = pointcolor, cex = 1)
  lay <- wordlayout(x, y, words, cex, ...)
  if (show.lines) {
    for (i in 1:length(x)) {
      xl <- lay[i, 1]
      yl <- lay[i, 2]
      w <- lay[i, 3]
      h <- lay[i, 4]
      if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] > 
          yl + h) {
        points(x[i], y[i], pch = 16, col = pointcolor, cex = 1)
        nx <- xl + 0.5 * w
        ny <- yl + 0.5 * h
        lines(c(x[i], nx), c(y[i], ny), col = line.col,lwd=line.width)
      }
    }
  }
  text(lay[, 1] + 0.5 * lay[, 3], lay[, 2] + 0.5 * lay[, 4], 
       words, cex = cex, ...)
}

png('Figure 2B.png',w=1750,h=1750,res=300)
scatter2(human,mouse,method='pearson',xlab=substitute(paste(italic('H. sapiens '),'PFOA 20µM (Day 14)')),
        ylab=substitute(paste(italic('M. musculus '),'PFOA 0.3 mg/kg body weight/day')),
        main=substitute(paste('Scatterplot ',italic('H. sapiens '),'vs ',italic('M. musculus'))),
        col='gray85',xlim=c(-20,50),ci=T,cex=1.5)
textplot2(human[genes],mouse[genes],words=genes,new=F,cex=0.85,col='black',show.lines=T,
          pointcolor='maroon1',font=2)
dev.off()



## Supplementary Figure S3 : scatterplot between zebrafish and human contrasts showing the highest 
# Pearson correlation (0.2) among these two species

# Step 1: extract the contrasts from signature matrix (Supplementary File S2) and 
# create the vectors for the scatterplot

mat2=subset(signatureMatrix,select=c('Hs_RO_D10_PFBS20','Dr_DA_24hpf_PFOSA12.5'))
mat2=na.omit(mat2)

human2=setNames(mat2$Hs_RO_D10_PFBS20,rownames(mat2))
Drerio=setNames(mat2$Dr_DA_24hpf_PFOSA12.5,rownames(mat2))


# Step 2: Select genes driving the correlation by p-value threshold ≤ 0.05

thr=-log10(0.03)
human2.up=human2[human2>=thr] 
human2.dn=human2[human2<=(-thr)]
Drerio.up=Drerio[Drerio>=thr]
Drerio.dn=Drerio[Drerio<=(-thr)]

up2=intersect(names(human2.up),names(Drerio.up))
dn2=intersect(names(human2.dn),names(Drerio.dn))
genes2=c(up2,dn2)


# Step 3: make the scatter plot

png('Supplementary Figure S3.png',w=1750,h=1750,res=300)
scatter2(human2,Drerio,method='pearson',xlab=substitute(paste(italic('H. sapiens '),'PFBS 20µM (Day 10)')),
        ylab=substitute(paste(italic('D. rerio '),'PFOSA 12.5µM (24 hpf)')),
        main=substitute(paste('Scatterplot ',italic('H. sapiens '),'vs ',italic('D. rerio'))),
        col='gray85',xlim=c(-20,5),ci=T,cex=1.5)
textplot2(human2[genes2],Drerio[genes2],words=genes2,new=F,cex=0.85,col='black',show.lines=T,
          pointcolor='maroon1',font=2)
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
#   [1] tm_0.7-11          wordcloud_2.6      NLP_0.2-1          RColorBrewer_1.1-3 corto_1.2.0        dplyr_1.1.0       
# [7] openxlsx_4.2.5.2  

# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.10        rstudioapi_0.14    xml2_1.3.3         magrittr_2.0.3     tidyselect_1.2.0   R6_2.5.1          
# [7] rlang_1.0.6        pbapply_1.7-0      fansi_1.0.4        caTools_1.18.2     tools_4.2.2        parallel_4.2.2    
# [13] plotrix_3.8-2      KernSmooth_2.23-20 utf8_1.2.2         cli_3.6.0          gtools_3.9.4       tibble_3.1.8      
# [19] lifecycle_1.0.3    zip_2.2.2          vctrs_0.5.2        bitops_1.0-7       slam_0.1-50        glue_1.6.2        
# [25] stringi_1.7.12     compiler_4.2.2     pillar_1.8.1       generics_0.1.3     gplots_3.1.3       pkgconfig_2.0.3

