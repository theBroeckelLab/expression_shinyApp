install.packages("openxlsx")
library(openxlsx)

## RUN 22
run22=read.table("Z:/Projects/Project Management/Cytotoxicity/Data/TKI_UWGS/run1_to_22_combined/run22_counts.counts", header=T, sep="\t", row.names = 1)
run1.to.21=read.table("Z:/Projects/Project Management/Cytotoxicity/Data/TKI_UWGS/run1_to_21_combined/counts_combined_run1_to_21.counts", header=F, sep="\t", row.names = 1)
colnames(run22)=sub(".", "c", colnames(run22))

#identical genes?
identical(rownames(run22), rownames(run1.to.21))
#bind data
run1.to.22=cbind(run1.to.21, run22)
write.table(run1.to.22, "Z:/Projects/Project Management/Cytotoxicity/Data/TKI_UWGS/run1_to_22_combined/run22_counts.counts", sep="\t", row.names=T, col.names=F)

#u01 phenos
phenos=read.xlsx("Z:/Projects/Project Management/Analysis/12.16.2019/Data Analysis/Data and Gene Sets/PhenotypeData/Copy of hiPSC_ID_master_list_041315_fromCharles_PA_110316_JB_12302019.xlsx")

#pull data
cat(colnames(run22), sep="\n")
cat(sapply(strsplit(colnames(run22), "_"), `[`, 1), sep="\n")
cat(sapply(strsplit(colnames(run22), "_"), `[`, 2), sep="\n")
cat(sapply(strsplit(colnames(run22), "_"), `[`, 3), sep="\n")

cat(phenos$sex[match(sapply(strsplit(colnames(run22), "_"), `[`, 1), paste0("c", phenos$id.cdi))], sep="\n")
cat(phenos$lvmht27[match(sapply(strsplit(colnames(run22), "_"), `[`, 1), paste0("c", phenos$id.cdi))], sep="\n")
cat(phenos$lvmht27.nw[match(sapply(strsplit(colnames(run22), "_"), `[`, 1), paste0("c", phenos$id.cdi))], sep="\n")
cat(phenos$race[match(sapply(strsplit(colnames(run22), "_"), `[`, 1), paste0("c", phenos$id.cdi))], sep="\n")
