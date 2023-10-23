
library(ggplot2)
library(DESeq2)




##Counts / Meta
meta <-read.table("Z:/Projects/Project Management/iPSC-CM U01 RNA sequencing/CEU/Counts/deseq2_design_fc1_to_7_to_charles.txt",row.names=1,header=T, sep="\t")
counts_table <- read.table("Z:/Projects/Project Management/iPSC-CM U01 RNA sequencing/CEU/Counts/ceu_u01_fc1_to_7_rnaseq_counts.txt", sep="\t", row.names = 1, header=T, stringsAsFactors = F)
colnames(counts_table) <- rownames(meta)
##filter data
meta[which(rownames(meta)%in%c("c1407_Unstim","c1287_Stim")),]
idx=which(rownames(meta)%in%c("c1407_Unstim","c1287_Stim"))
meta <- meta[-idx,]
counts_table=counts_table[,-idx]
##expression analysis
counts_table <- counts_table[apply(counts_table,1,median)>=10,]
dat <- varianceStabilizingTransformation(as.matrix(counts_table))


#### PCA #####
pca <- prcomp(t(dat), scale=F)
percentVar <- pca$sdev^2/sum(pca$sdev^2)
pc1<-1; pc2<-2
d <- data.frame(PCA=pca$x[,pc1], PCB=pca$x[,pc2], meta)
colnames(d)
ggplot(data=d, aes(x=PCA, y=PCB, colour=Condition)) +
  theme_bw()+
  geom_hline(yintercept = 0, colour = "gray65") + 
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(size=7) + 
  #scale_color_manual(values=cols) +
  xlab(paste0("PC1 (", round((percentVar[pc1]*100),2), "%)")) + 
  ylab(paste0("PC2 (", round((percentVar[pc2]*100),2), "%)")) + 
  theme(legend.title=element_blank(), 
        legend.text = element_text(size=15),
        axis.title=element_text(size=15))
View(pca$x)



#c1287_Stim, c1407_Unstim
