## Estimating PC Association using ANOVA ## 


#### Read-In Expression Data ####

#### DB Data ####
# counts_table <- read.table("Z:/Projects/Project Management/DB Matrix/Data/AmpliSeq_Counts_Oct19/data_for_charles/raw_counts.tsv", sep="\t", header=T, row.names = 1)
# meta_data <- read.table("Z:/Projects/Project Management/DB Matrix/Data/AmpliSeq_Counts_Oct19/data_for_charles/deseq2_design.tsv", sep="\t", header=T, row.names = 1)
# Input = list();
# Input$count_data <- counts_table
# Input$meta_data <- meta_data
# colnames(Input$count_data) <- rownames(Input$meta_data)
# Input_filtered <- Input
# samples <- rownames(meta_data)
# Input_filtered$meta_data <- Input$meta_data[samples,]
# #for (i in colnames(Input_filtered$meta_data)) {
#   #print(class(Input_filtered$meta_data[[i]]))
# #  Input_filtered$meta_data[[i]] <-factor(Input_filtered$meta_data[[i]])
# #}
# Input_filtered$count_data <- Input$count_data[,samples]
# Input_filtered$count_data <- Input_filtered$count_data[apply(Input_filtered$count_data,1,median)>=10, ]
# sampleinfo <- Input_filtered$meta_data
# sampleinfo <- DataFrame(sampleinfo)
# countdata <- Input_filtered$count_data
# dds_from_counts <- DESeqDataSetFromMatrix(countData = countdata, colData = sampleinfo, design = ~1)
# x <- DESeq2::varianceStabilizingTransformation(dds_from_counts)



#### DRC Data ####
counts_table <- read.table("Z:/Projects/Project Management/Cytotoxicity/Data/TKI_UWGS/run1_2_3_4/raw_counts_nwgc_run1_2_3_4_combined_genenames.counts", sep="\t", header=T, row.names=1)
meta_data <- read.table("Z:/Projects/Project Management/Cytotoxicity/Data/TKI_UWGS/run1_2_3_4/deseq2_design.txt", sep="\t", header=T, row.names=1)

#counts_table <- read.table("Z:/Projects/Project Management/Cytotoxicity/Data/TKI_UWGS/run1_to_5_combined/raw_counts_nwgc_run1_2_3_4_5_combined_genenames.counts", sep="\t", header=T, row.names=1)
#meta_data <- read.table("Z:/Projects/Project Management/Cytotoxicity/Data/TKI_UWGS/run1_to_5_combined/deseq2_design_1_5.txt", sep="\t", header=T, row.names=1)

#counts_table <- read.table("Z:/Projects/Project Management/iPSC-CM U01 RNA sequencing/Raw data/Counts_byPA/FC1_8/raw_counts.counts", sep="\t", header=T, row.names=1)
#meta_data <- read.table("Z:/Projects/Project Management/iPSC-CM U01 RNA sequencing/Raw data/Counts_byPA/FC1_8/deseq2_design_for_Anova.txt", sep="\t", header=T, row.names=1)

Input = list();
Input$count_data <- counts_table
Input$meta_data <- meta_data
colnames(Input$count_data) <- rownames(Input$meta_data)
Input_filtered <- Input
#samples <- rownames(meta_data)
#samples <- rownames(meta_data[which((meta_data$Project == "TKI") & (meta_data$SampleID != 'c1060_DMSO2') & (meta_data$Cohort != 'iCell')),])
samples <- rownames(meta_data[which((meta_data$Project == "DRC" & meta_data$Condition=="Unstim")),])
Input_filtered$meta_data <- Input$meta_data[samples,]
Input_filtered$count_data <- Input$count_data[,samples]
Input_filtered$count_data <- Input_filtered$count_data[apply(Input_filtered$count_data,1,median)>=0,]
sampleinfo <- Input_filtered$meta_data
sampleinfo <- DataFrame(sampleinfo)
countdata <- Input_filtered$count_data
dds_from_counts <- DESeqDataSetFromMatrix(countData = countdata, colData = sampleinfo, design = ~1)
x <- varianceStabilizingTransformation(dds_from_counts)

#### Retain Top N Most Variable Genes ####
rv <- rowVars(assay(x))
select <- order(rv, decreasing = TRUE)[seq_len(min(length(rv), length(rv)))]
pca <- prcomp(t(assay(x)[select, ]), scale = F)

#### Run ANOVA ####
meta <- Input_filtered$meta_data
pc_anova <- list()
prohoc_all <- list()
prohoc_rownames <- list()
max.pca <- min(min(dim(pca$x)), 25)
for (pc in 1:max.pca) {
  for (var in colnames(meta)) {
    #skip if theres only one group or if the # of groups equals the # of samples
    if (!is.factor(meta$Project) | 
        length(unique(meta[[var]])) == 1 | 
        length(unique(meta[[var]])) == nrow(meta)) {
      next
    }
    else {
      #run ANOVA (PC~MetaVariable) and save p value
      s <- summary(aov(as.numeric(pca$x[,pc])~meta[[var]]))
      thsd <- TukeyHSD(aov(as.numeric(pca$x[,pc])~meta[[var]]))
      prohoc <- data.frame(thsd$`meta[[var]]`)
      prohoc_all[[var]][pc] <- prohoc["p.adj"]
      pc_anova[[var]][pc] <- s[[1]][["Pr(>F)"]][1]
      prohoc_rownames[[var]] <- rownames(prohoc["p.adj"])
    }
  }
}

#### Create Table of p Values #####
pc_anova <- as.data.frame(pc_anova)
rownames(pc_anova) <- paste0("PC", rownames(pc_anova))
View(pc_anova)

# ##### Linear Model for Continous Variables ######
# pc_lm <- list()
# for (pc in 1:25) {
#   print(paste0("Currently on...PC", pc))
#   for (var in colnames(meta_data)[c(6,11,12,13,14,15,16,21,22)]) {
#     pc_lm[[var]][pc] <- summary(lm(as.numeric(pca$x[,pc])~as.numeric(meta_data[[var]])))[[4]][2,4] }
# }
# pc_lm <- as.data.frame(pc_lm)
# rownames(pc_lm) <- paste0("PC", rownames(pc_lm))
# View(pc_lm)


#### Verify with PCA ####
intgroup <- c("Gender")
percentVar <- pca$sdev^2/sum(pca$sdev^2)
intgroup.df <- GenomicRanges::as.data.frame(colData(x)[, intgroup, drop=FALSE])
group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
#group <- reorder(group, new.order=c(8,6,4,2,7,5,3,1))
pc1 <- 3
pc2 <- 4
d <- data.frame(PCA = pca$x[, pc1], PCB = pca$x[, pc2], 
                group = group, intgroup.df, 
                names = colnames(x))
ggplot(data = d, aes_string(x = "PCA", y = "PCB", colour = group, label = "names")) +
  geom_hline(yintercept = 0, colour = "gray65") + 
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(size=5) + 
  xlab(paste0("PC",pc1,": ", round(percentVar[pc1] * 100), "% variance")) + 
  ylab(paste0("PC",pc2,": ", round(percentVar[pc2] * 100), "% variance"))
  #+scale_color_brewer(type="seq", direction=1, palette=1)

