# This is the code for Figure 7. As network for the metabolites prediction on human datasets, we used the network of
# Cavicchioli et al. 2022 (PMID: 35409231). Supplementary File S3 contains the normalized enrichment scores of 147 
# predicted metabolites for 55 PFAS exposure contrasts of human datasets.

library(openxlsx)
library(corto)
library(dplyr)
library(ggplot2)

metab.matrix=read.xlsx('Supplementary File S3.xlsx',rowNames=T)


# Step 1: integration to identify top up- and down-regulated metabolites

integration=as.data.frame(apply(metab.matrix,1,stouffer))
colnames(integration)='stouffer_coefficient'

integration$sd=apply(metab.matrix,1,sd)

top_up=arrange(integration,desc(stouffer_coefficient))[1:10,]


top_down=arrange(integration,stouffer_coefficient)[1:10,]
top_down=arrange(top_down,desc(stouffer_coefficient))

topmetab=rbind(top_up,top_down)
topmetab$metabolite=as.factor(rownames(topmetab))
topmetab$Effect=c(rep('Upregulated',10),rep('Downregulated',10))
topmetab$pvalue=z2p(topmetab$stouffer_coefficient)


# Step 2: Figure 7 - barplot of top metabolites

order=topmetab$metabolite

png('Figure 7.png',w=2300,h=1700,res=420)
barplot=ggplot(topmetab)
barplot + geom_bar(aes(x=factor(metabolite,levels=order),y=stouffer_coefficient,fill=Effect),
                   stat='identity',position='dodge') + 
  geom_errorbar(aes(x=factor(metabolite,levels=order),ymin=stouffer_coefficient-sd,
                    ymax=stouffer_coefficient+sd),width=0.4,colour=c(rep('red3',10),rep('navy',10)),
                alpha=1,linewidth=0.5) +
  geom_hline(yintercept=p2z(1e-10),linetype='dashed',color='gray') + # pvalue threshold 
  geom_hline(yintercept=-p2z(1e-10),linetype='dashed',color='gray') +
  scale_fill_manual(values=c('navy','red3')) + 
  labs(title='Top metabolites of human datasets') + 
  xlab(label='Metabolites') + 
  ylab(label='Integrated Normalized Enrichment Score') +
  theme_bw() +
  theme(axis.text.x.bottom=element_text(angle=45,vjust=1,hjust=1,colour='black'))
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
#   [1] ggplot2_3.4.1    dplyr_1.1.0      corto_1.2.0      openxlsx_4.2.5.2

# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.10        rstudioapi_0.14    magrittr_2.0.3     munsell_0.5.0      tidyselect_1.2.0   colorspace_2.1-0  
# [7] R6_2.5.1           rlang_1.0.6        pbapply_1.7-0      fansi_1.0.4        caTools_1.18.2     tools_4.2.2       
# [13] grid_4.2.2         parallel_4.2.2     gtable_0.3.1       plotrix_3.8-2      KernSmooth_2.23-20 utf8_1.2.2        
# [19] cli_3.6.0          withr_2.5.0        gtools_3.9.4       tibble_3.1.8       lifecycle_1.0.3    zip_2.2.2         
# [25] farver_2.1.1       vctrs_0.5.2        bitops_1.0-7       glue_1.6.2         labeling_0.4.2     stringi_1.7.12    
# [31] compiler_4.2.2     pillar_1.8.1       scales_1.2.1       generics_0.1.3     gplots_3.1.3       pkgconfig_2.0.3



