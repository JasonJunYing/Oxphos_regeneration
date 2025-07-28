#########EdgeR#########

library(edgeR)
library(xlsx)
library(pheatmap)
#Import count matrix
counts<-read.csv('./Files/scores_gimme_1000X.csv',row.names = 1)
counts$id <- rownames(counts)
counts <- as.data.frame(lapply(counts[,1:(ncol(counts)-1)], round, 3),row.names = counts$id)
counts0 <- counts

### Heatmap ###
data <- counts0

# 7dpci only
data <- data[,grepl('7Dpci',colnames(data))]

TCA_reactions <- read.delim('Reactions for subsystem tricarboxylic_acid_cycle_and_glyoxylate_dicarboxylate_metabolism.tsv')
TCA_reactions <- TCA_reactions[!duplicated(TCA_reactions$Genes),]
data <- data[TCA_reactions$Reaction.ID,]
data <- data[!is.na(data[,1]),]

termname <- "TCA"
ls <- ""

sampname <- colnames(data)
meta<-data.frame('Strain'=strsplit2(sampname,'_')[,1],
                 'timepoint'=strsplit2(sampname,'_')[,2],
                 'rep'=strsplit2(sampname,'_')[,3])

annotation_row <- data.frame(Group = meta$Group)
annotation_row <- meta[,1:2]
row.names(annotation_row) <- colnames(data)

# Get reaction info
id_anno$anno <- paste0(id_anno$id,'_',id_anno$name)
data$id <- rownames(data)
data <- merge(data,id_anno,by = 'id')

rownames(data)<-data$anno
data <- data[,-1]
data <- data[,-ncol(data)]
data <- data[,-ncol(data)]
data <- data[,-ncol(data)]

data <- data[rowSums(data)>0,]

annotation_color <- list(timepoint = c(
  "7Dpci"="orange"
  )
  )
p<-pheatmap(data,
            scale = "row",
            color = colorRampPalette(c("blue", "white","red"))(100),
            annotation = annotation_row,
            annotation_colors = annotation_color,
            cluster_rows = T,
            cluster_cols = F,
            show_colnames = T, # Parameters could be changed if necessary
)

pdf(paste0('./figures/',comp,'_',ls,'_',study,'_',termname,'_heatmap.pdf'),height = 6,width = 12)
print(p)
dev.off()

library(Cairo)
CairoPNG(paste0('./figures/',comp,'_',ls,'_',study,'_',termname,'_heatmap.png'),res = 300,unit="in",height = 6,width = 12)
print(p)
dev.off()
