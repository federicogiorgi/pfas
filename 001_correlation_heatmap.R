## This is the code for Figure 1:
## Correlation heatmap of signature matrix (Supplementary File S2) containing humanized genes
## In order to verify if exposure to PFAS molecules stimulate correlated responses in distinct species,
## Pearson correlation was calculated on signature matrix

library(openxlsx)
library(circlize)
library(ComplexHeatmap)


signatureMatrix=read.xlsx('SupplementaryFileS3.xlsx',rowNames=T)

cormat=cor(signatureMatrix,use='pairwise.complete.obs',method='pearson')


# Step 1: create the annotations of heatmap

Sp=colnames(cormat)
col1=c(rep('darkblue',13),rep('deeppink',55),rep('green3',14),rep('orangered',2),
             rep('darkmagenta',12),rep('turquoise3',8),rep('goldenrod1',6))
names(col1)=Sp

tissue=c(rep('Liver',13),rep('Cancer xenograft',2),rep('Liver',53),rep('Whole body',14),
         rep('Embryo tissues',2),rep('Reproductive system',16),rep('Liver',7),rep('Blood',3))
col2=c(rep('dodgerblue',13),rep('gray17',2),rep('dodgerblue',53),rep('antiquewhite1',14),
       rep('olivedrab1',2),rep('plum2',16),rep('dodgerblue',7),rep('indianred2',3))
names(col2)=tissue


# Step 2: create the heatmap

CorHeatmap=Heatmap(cormat,col=colorRamp2(c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),
                                         c('midnightblue',"blue3",'mediumblue','royalblue3','deepskyblue',
                                           'white','orangered','firebrick1','red2',"red3",'#660000')),
                   cluster_rows=F,cluster_columns=F,show_row_names=F,show_column_names=F,
                   left_annotation=rowAnnotation(Species=Sp,col=list(Species=col1),
                                                 show_annotation_name=F,show_legend=F),
                   top_annotation=columnAnnotation(Tissues=tissue,Species=Sp,
                                                   col=list(Tissues=col2,Species=col1),
                                                   show_annotation_name=T,
                                                   annotation_name_gp=gpar(fontsize=20,
                                                                           fontface='bold'),
                                                   show_legend=F,gap=unit(2,'mm')),
                   show_heatmap_legend=F)


# Step 3: create the legends for the heatmap

lgdHeatmap=Legend(col_fun=colorRamp2(c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),
                                     c('midnightblue',"blue3",'mediumblue','royalblue3','deepskyblue',
                                       'white','orangered','firebrick1','red2',"red3",'#660000')),
                  at=c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),title='Pearson correlation',
                  title_gp=gpar(fontsize=20,fontface='bold'),legend_width=unit(20,'cm'),
                  direction='horizontal',title_position='topcenter',
                  labels_gp=gpar(col='black',fontsize=20),title_gap=unit(3,'mm'))


col3=c('darkblue','deeppink','green3','orangered','darkmagenta','turquoise3','goldenrod1')
lgdSpecies=Legend(labels=c("M. musculus","H. sapiens","C. elegans",'D. rerio','G. morhua',
                           'M. salmoides','P. promelas'),
                  title='Species',title_gp=gpar(fontsize=20,fontface='bold'),
                  labels_gp=gpar(col='black',font=3,fontsize=22),
                  legend_gp=gpar(fill=col3),grid_height=unit(1,'cm'),grid_width=unit(5, "mm"),
                  title_gap=unit(3,'mm'))


col4=c('dodgerblue','gray17','antiquewhite1','olivedrab1','plum2','indianred2')
lgdTissue=Legend(labels=c('Liver','Cancer xenograft','Whole body','Embryo tissues',
                          'Reproductive system','Blood'),
                 title='Tissues',title_gp=gpar(fontsize=20,fontface='bold'),
                 labels_gp=gpar(col='black',fontsize=22),
                 legend_gp=gpar(fill=col4),grid_height=unit(1,'cm'),grid_width=unit(5,"mm"),
                 title_gap=unit(3,'mm'))

legends=packLegend(lgdSpecies,lgdTissue,gap=unit(1.2,'cm'))


# Step 4: plot the heatmap

png('Figure 1.png',w=4400,h=4000,res=300)
draw(CorHeatmap,heatmap_legend_list=lgdHeatmap,heatmap_legend_side='bottom',
     annotation_legend_list=legends)
dev.off()


sessionInfo()

# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)

# Matrix products: default

# locale:
#  [1] LC_COLLATE=Italian_Italy.utf8  LC_CTYPE=Italian_Italy.utf8    LC_MONETARY=Italian_Italy.utf8
# [4] LC_NUMERIC=C                   LC_TIME=Italian_Italy.utf8    

# attached base packages:
#  [1] grid      stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#  [1] ComplexHeatmap_2.14.0 circlize_0.4.15       openxlsx_4.2.5.2     

# loaded via a namespace (and not attached):
#  [1] zip_2.2.2           Rcpp_1.0.10         RColorBrewer_1.1-3  pillar_1.8.1        compiler_4.2.2     
# [6] iterators_1.0.14    tools_4.2.2         digest_0.6.31       clue_0.3-64         lifecycle_1.0.3    
# [11] tibble_3.1.8        gtable_0.3.1        pkgconfig_2.0.3     png_0.1-8           rlang_1.0.6        
# [16] foreach_1.5.2       cli_3.6.0           rstudioapi_0.14     parallel_4.2.2      cluster_2.1.4      
# [21] dplyr_1.1.0         S4Vectors_0.36.1    IRanges_2.32.0      generics_0.1.3      GlobalOptions_0.1.2
# [26] vctrs_0.5.2         stats4_4.2.2        tidyselect_1.2.0    glue_1.6.2          R6_2.5.1           
# [31] GetoptLong_1.0.5    fansi_1.0.4         ggplot2_3.4.1       magrittr_2.0.3      BiocGenerics_0.44.0
# [36] scales_1.2.1        codetools_0.2-19    matrixStats_0.63.0  shape_1.4.6         colorspace_2.1-0   
# [41] utf8_1.2.2          stringi_1.7.12      munsell_0.5.0       doParallel_1.0.17   rjson_0.2.21       
# [46] crayon_1.5.2        Cairo_1.6-0 



