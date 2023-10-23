cmrg=read.xlsx("W:/Clinical Working/SVI Project/Oxford Nanopore Technologies/Nanopore data from Maple/MCW-SVI-0034/mandelker_genes.xlsx")
human=useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_ids=cmrg$gene
View(listAttributes(human))
out=getBM(attributes=c('hgnc_symbol', 'chromosome_name','start_position','end_position', 'description'),
          filters = 'hgnc_symbol',
          values = gene_ids,
          mart = human)

write.xlsx(cmrg_towrite, rowNames=T, file="Z:/Projects/Project Management/Cytotoxicity/Data/TKI_UWGS/run1_to_15_combined_and_DOX/DOX/DE_results/Male_vs_Female_DMSO_resultsTable.xlsx")