#########EdgeR#########

library(edgeR)
library(xlsx)
#Import count matrix
counts<-read.delim('./fish_strains_project2_raw_counts.txt')
rownames(counts)<-counts$id
id_ref <- counts[,1:2]
counts <- counts[,4:(ncol(counts)-1)]
# counts0 <- counts

# NA v SAT
study <- '7dpci'
counts <- counts[,c(22:24,30:32)] # 7dpci

#Define group design
sampname <- colnames(counts)
meta<-data.frame('Group'=strsplit2(sampname,'_')[,1],
                 'timepoint'=strsplit2(sampname,'_')[,2],
                 'rep'=strsplit2(sampname,'_')[,3])
# meta0<-meta
design<-model.matrix(~0+Group,data = meta)

#Filter out low expression genes
dgelist<-DGEList(counts = counts,group = meta$Group)
keep<-filterByExpr(dgelist,group = meta$Group)
dge<-dgelist[keep,]

# Normalization and fit model
dgelist_norm <- calcNormFactors(dge, method = 'TMM')
dge <- estimateGLMCommonDisp(dgelist_norm, design)
dge <- estimateGLMTrendedDisp(dge,design)
dge <- estimateGLMTagwiseDisp(dge,design)
fit <- glmFit(dge, design, robust = TRUE)

# Make contrast
comp <- "NAvSAT"

cont<-makeContrasts(GroupNA-GroupSAT,levels = design)

lrt <- topTags(glmLRT(fit,contrast = cont), n = nrow(dgelist$counts))
write.table(lrt, paste0(comp,'_',study,'_glmLRT.txt'), sep = '\t', col.names = NA, quote = FALSE)

library(ggplot2)
# MAplot
lrt<-read.delim(paste0(comp,'_',study,'_glmLRT.txt'),row.names = 1)
plot(lrt$logFC,lrt$logCPM,pch=16)
title(paste0(comp,'_',study))
plot.new()
hist(lrt$PValue,main="")
title(paste0(comp,'_',study))

lrt$id <- rownames(lrt)
lrt <- merge(lrt,id_ref,by = 'id')
lrt <- lrt[order(lrt$logFC,decreasing = T),]
data<-lrt

###Volcano###
# Select DEG
down<- intersect(which(data$FDR<=0.1) , which(data$logFC<=(-log2(1.5))))
up<- intersect(which(data$FDR<=0.1) , which(data$logFC>=log2(1.5)))

downlrt <- lrt[down,]
downlrt <- downlrt[order(downlrt$logFC),]
uplrt <- lrt[up,]

write.csv(downlrt,paste0('downlrt_',comp,'_',study,'.ONLY.csv'))
write.csv(uplrt,paste0('uplrt_',comp,'_',study,'.ONLY.csv'))

### EnrichR
library(enrichR)
setEnrichrSite("FishEnrichr")
dbs_list <- listEnrichrDbs()
dbs <- c("GO_Biological_Process_2018", "WikiPathways_2018","KEGG_2019")

comp='NAvSAT'
study='7dpci'

downlrt <- read.csv(paste0('downlrt_',comp,'_',study,'.ONLY.csv'))
uplrt <- read.csv(paste0('uplrt_',comp,'_',study,'.ONLY.csv'))

deglrt <- list("down"=downlrt$symbol,"up"=uplrt$symbol)

library(Cairo)
enriched <- list()
for(ls in c("down","up")){
  enriched[[ls]] <- enrichr(deglrt[[ls]], dbs)

  for(db in dbs){
    enr <- enriched[[ls]]
    p <- plotEnrich(enr[[db]], showTerms = 20, numChar = 40,
                    y = "Ratio", orderBy = "P.Value",
                    title=paste0(comp,'_',ls,'_',study),
                    xlab=db
                    )
    pdf(paste0(comp,'_',ls,'_',study,'_',db,'_enrichR.pdf'),height = 5,width = 6)
    print(p)
    dev.off()

    CairoPNG(paste0(comp,'_',ls,'_',study,'_',db,'_enrichR.png'),height = 5,width = 6,units="in",res = 300)
    print(p)
    dev.off()
  }
}

saveRDS(enr,paste0('enr_',comp,'_',ls,'_',study,'_',db,'_enrichR.rds'))
kegg_7dpci_up <- enriched$up$KEGG_2019
write.table(kegg_7dpci_up,'kegg_7dpci_up.txt')

