# This is the code for fgsea analysis carried out on signature matrix (Supplementary File S2), in order to identify 
# the pathways mostly up- and down-regulated across species, and the code for Figure 6.

library(openxlsx)
library(msigdbr)
library(fgsea)
library(corto)
library(dplyr)
library(tidyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)

signatureMatrix=read.xlsx('SupplementaryFileS3.xlsx',rowNames=T)


# Step 1: Prepare msigdbr list

fname="mlist.rda"
if(!file.exists(fname)){
  as.data.frame(msigdbr_species())
  as.data.frame(msigdbr_collections())
  # C2 KEGG
  mdf=msigdbr(species="Homo sapiens",category="C2",subcategory="CP:KEGG")
  mlist=mdf %>% split(x=.$gene_symbol,f=.$gs_name)
  # C2 Wiki Pathways
  mdf=msigdbr(species="Homo sapiens",category="C2",subcategory="CP:WIKIPATHWAYS")
  mlist=c(mlist,mdf %>% split(x=.$gene_symbol,f=.$gs_name))
  # C5 Gene Ontology (excluding HPO)
  mdf=msigdbr(species="Homo sapiens",category="C5",subcategory=c("GO:BP"))
  mlist=c(mlist,mdf %>% split(x=.$gene_symbol,f=.$gs_name))
  mdf=msigdbr(species="Homo sapiens",category="C5",subcategory=c("GO:CC"))
  mlist=c(mlist,mdf %>% split(x=.$gene_symbol,f=.$gs_name))
  mdf=msigdbr(species="Homo sapiens",category="C5",subcategory=c("GO:MF"))
  mlist=c(mlist,mdf %>% split(x=.$gene_symbol,f=.$gs_name))
  # H Hallmark pathways
  mdf=msigdbr(species="Homo sapiens",category="H")
  mlist=c(mlist,mdf %>% split(x=.$gene_symbol,f=.$gs_name))
  save(mlist,file=fname)
}else{load(fname)}


# Step 2: perform fgsea for each contrast of signature matrix

fgseas=list()
for(i in colnames(signatureMatrix)){
    message('Doing ',i)
    sign=signatureMatrix[,i]   
    names(sign)=rownames(signatureMatrix)
    sign=sign[!is.na(sign)]
    set.seed(1)
    sign=sign+rnorm(length(sign),sd=0.0001)
    fgseaResults=fgsea(pathways=mlist,stats=sign,minSize=5,maxSize=Inf,nproc=8)
    fgseaResults=fgseaResults[order(fgseaResults$pval),]
    fgseas[[i]]=fgseaResults
}


# Step 3: create a matrix with NES of pathways for each contrast of signature matrix and fill the matrix
# with Normalized Enrichment Score

pathways=c()
for(i in 1:length(fgseas)){
  pathways=union(pathways,fgseas[[i]]$pathway)
}

gseamat=matrix(NA,nrow=length(pathways),ncol=length(fgseas))
rownames(gseamat)=pathways
colnames(gseamat)=names(fgseas)

for(i in 1:length(fgseas)){
  nes=setNames(fgseas[[i]]$NES,fgseas[[i]]$pathway)
  gseamat[,i]=nes[rownames(gseamat)]
}


# Step 4: calculate the p-value and the p-adjusted of the pathways through the integrated scores

pvalmat=as.data.frame(apply(gseamat,1,function(x){
  mat=x[!is.na(x)]
  z2p(stouffer(mat))
}))
colnames(pvalmat)='pval'

pvalmat$padj=p.adjust(pvalmat$pval,method='BH')


# Step 5: select most up- and down-regulated pathways across species

mat=replace(gseamat,is.na(gseamat),0)

up=apply(mat,1,function(x){
  sum(x>0)
})
pathup=names(sort(up,decreasing=TRUE)[1:20]) # Mostly upregulated pathways

dn=apply(mat,1,function(x){
  sum(x<0)
})
pathdn=names(sort(dn,decreasing=TRUE)[1:20]) # Mostly downregulated pathways

topgsea=as.data.frame(mat[c(pathup,pathdn),])


# Step 6: select the p-value of top up- and down-regulated pathways and order pathways by p-value

pvalmatUp=filter(pvalmat,rownames(pvalmat)%in%pathup)
pvalmatUp=arrange(pvalmatUp,padj)

pvalmatDn=filter(pvalmat,rownames(pvalmat)%in%pathdn)
pvalmatDn=arrange(pvalmatDn,desc(padj))

pvalmat=rbind(pvalmatUp,pvalmatDn)


# Step 7: reorder the top differentially expressed pathways (topgsea data frame) by p-adjusted 

order=rownames(pvalmat)
topgsea$path=rownames(topgsea)
mat2=topgsea %>% mutate(path=factor(path,levels=order)) %>% arrange(path)
topgsea=mat2[,-111]


# Step 8: modify names of pathways

name=rownames(topgsea)
name1=gsub('_',' ',name)
name2=str_to_title(name1)
name3=gsub('Wp ','WP:',name2)
name4=gsub('Gobp ','GOBP:',name3)
name5=gsub('Gocc ','GOCC:',name4)
name6=gsub('Gomf ','GOMF:',name5)
name7=gsub('GOBP:Dna','GOBP:DNA',name6)

rownames(topgsea)=name7
rownames(pvalmat)=name7



### Figure 6: heatmap of top differetially expressed pathways across species
# Step 9: create the annotation of the heatmap

Sp=colnames(topgsea)
colSp=c(rep('darkblue',13),rep('deeppink',55),rep('green3',14),rep('orangered',2),
        rep('darkmagenta',12),rep('turquoise3',8),rep('goldenrod1',6))
names(colSp)=Sp


# Step 10: superscript the p-adjusted

p=pvalmat$padj
p=signif(p,3)
padj=gsub('e.+','x10',p) 
padj=as.character(padj)

pvalmat=separate(pvalmat,padj,c('pad','exp'),sep='e',remove=F)
exp=pvalmat$exp

Padj=mapply(function(x,y){
  as.expression(bquote(.(x)^.(y)))
},padj,exp)


# Step 11: create the heatmap

gseaHeat=Heatmap(as.matrix(topgsea),col=colorRamp2(c(-2,0,2),c('navy','white','red3')),
                 cluster_rows=F,cluster_columns=F,show_column_names=F,
                 bottom_annotation=columnAnnotation(Species=Sp,col=list(Species=colSp),
                                                    show_annotation_name=T,show_legend=F,
                                                    annotation_name_gp=gpar(fontsize=18,fontface='bold'),
                                                    gap=unit(2, "points"),height=unit(1,'cm')),
                 left_annotation=rowAnnotation(padjusted=anno_text(Padj,which='row',
                                                                   gp=gpar(fontsize=20,col='black')),
                                               show_annotation_name=F,show_legend=F),
                 show_heatmap_legend=F,width=unit(20, "cm"),height=nrow(as.matrix(topgsea))*unit(1,'cm'),
                 row_names_gp=gpar(fontsize=20))


# Step 12: create heatmap legends

lgdHeat=Legend(col_fun=colorRamp2(c(-2,0,2),c('navy','white','red3')),at=c(-2,-1,0,1,2),
               title='Normilized Enrichment Score',title_gp=gpar(fontsize=20,fontface='bold'),
               legend_width=unit(15,'cm'),direction='horizontal',title_position='topcenter',
               labels_gp=gpar(col='black',fontsize=15,title_gap=unit(3,'mm')),gap=unit(3,"mm"))

col=c('darkblue','deeppink','green3','orangered','darkmagenta','turquoise3','goldenrod1')

lgdAnno=Legend(labels=c("M. musculus","H. sapiens","C. elegans",'D. rerio','G. morhua',
                        'M. salmoides','P. promelas'),title='Species',
               title_gp=gpar(fontsize=20,fontface='bold'),
               labels_gp=gpar(col='black',font=3,fontsize=20),legend_gp=gpar(fill=col),
               grid_height=unit(4,'mm'),grid_width=unit(4, "mm"),title_gap=unit(3,'mm'),ncol=7)

legends=packLegend(lgdAnno,lgdHeat,gap=unit(1,'cm'),direction='vertical')


# Step 13: plot the heatmap

png('Figure 6.png',w=10500,h=5500,res=300)
draw(gseaHeat,heatmap_legend_list=legends,heatmap_legend_side='top')
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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] ComplexHeatmap_2.14.0 circlize_0.4.15       stringr_1.5.0         tidyr_1.3.0           dplyr_1.1.0          
# [6] corto_1.2.0           fgsea_1.24.0          msigdbr_7.5.1         openxlsx_4.2.5.2     

# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.10         lattice_0.20-45     snow_0.4-4          png_0.1-8           gtools_3.9.4       
# [6] digest_0.6.31       foreach_1.5.2       utf8_1.2.2          R6_2.5.1            stats4_4.2.2       
# [11] ggplot2_3.4.1       pillar_1.8.1        gplots_3.1.3        GlobalOptions_0.1.2 rlang_1.0.6        
# [16] rstudioapi_0.14     data.table_1.14.8   S4Vectors_0.36.1    GetoptLong_1.0.5    Matrix_1.5-3       
# [21] BiocParallel_1.32.5 munsell_0.5.0       compiler_4.2.2      pkgconfig_2.0.3     BiocGenerics_0.44.0
# [26] shape_1.4.6         tidyselect_1.2.0    tibble_3.1.8        IRanges_2.32.0      codetools_0.2-19   
# [31] matrixStats_0.63.0  fansi_1.0.4         withr_2.5.0         crayon_1.5.2        bitops_1.0-7       
# [36] gtable_0.3.1        lifecycle_1.0.3     magrittr_2.0.3      scales_1.2.1        KernSmooth_2.23-20 
# [41] zip_2.2.2           cli_3.6.0           stringi_1.7.12      pbapply_1.7-0       doParallel_1.0.17  
# [46] generics_0.1.3      vctrs_0.5.2         cowplot_1.1.1       fastmatch_1.1-3     RColorBrewer_1.1-3 
# [51] rjson_0.2.21        Cairo_1.6-0         iterators_1.0.14    tools_4.2.2         glue_1.6.2         
# [56] purrr_1.0.1         plotrix_3.8-2       parallel_4.2.2      clue_0.3-64         babelgene_22.9     
# [61] colorspace_2.1-0    cluster_2.1.4       caTools_1.18.2