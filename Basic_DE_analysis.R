#### Attempt One to make a basic DE analysis script
#### Praful Aggarwal 4/20/2015 - 8:47AM
#### 0.1

## prepare the output folder and a single object that contains the countdata and the metadata
custom_shape_palette <<- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)





#####################
## Shiny Front End ##
#####################

##**## STEP1: COMBINE META AND COUNTS
consolidate_genes_samples <- function(counts_table = NULL, meta_data = NULL){
  print("STEP1: COMBINE META AND COUNTS")
  
  Input = list();
  
  if (is.null(counts_table) || is.null(meta_data)){
    stop("Please provide both the counts file and the meta-data file.")
  }
  else{
    Input$count_data <- counts_table
    Input$meta_data <- meta_data
    colnames(Input$count_data) <- rownames(Input$meta_data)
  }
  
  input_object <<- Input
  
  if (dim(input_object$count_data)[2] != dim(input_object$meta_data)[1]) {
    stop("Check Input Files! Number of Samples in Design Matrix doesn't Match Number of Samples in Counts Table") }
  
  return(input_object)
}


##**## STEP 2: QC NORMALIZATION
read_transformation <- function(x, type="vst"){
  start <- Sys.time()
  print("START STEP 2: QC NORMALIZATION")
  
  type_vals <- c("rlog","vst")
  
  if (type == "rlog"){
    rtransformed <- DESeq2::rlog(x) } 
  else if (type == "vst"){
    rtransformed <- DESeq2::varianceStabilizingTransformation(x) }
    #x <- x+1
    #rtransformed <- DESeq2::varianceStabilizingTransformation(x)
    #rtransformed[which(is.infinite(as.matrix(rtransformed)))] <- 0}
  
  print("END STEP 2: QC NORMALIZATION")
  
  end <- Sys.time()
  print(paste0("**Runtime of Normalization: ", round(difftime(end,start,units="mins"),4), " Minutes"))
  return(rtransformed)
}


##**## STEP 3: RUN DESEQ NORMALIZATION
run_DESeq2 <- function(Input){
  
  print("START STEP 3: RUN DUMMY DESEQ NORMALIZATION")
  

  sampleinfo <- Input$meta_data
  sampleinfo <- DataFrame(sampleinfo)
  countdata <- Input$count_data
  
  sampleinfo[sapply(sampleinfo, is.character)]=lapply(sampleinfo[sapply(sampleinfo, is.character)], as.factor)
  dds <- DESeqDataSetFromMatrix(countData = countdata, colData = sampleinfo, design = ~1)
  
  print("END STEP 3: RUN DUMMY DESEQ NORMALIZATION")
  return(dds)
}


##**## STEP 4: FILTER BY MIN MEDIAN
filter_gene_sample_list <- function(Input, minMedian = 10,  samples = rownames(raw_metadata)){
  print("START STEP 4: FILTER BY MIN MEDIAN")
  Input_filtered <- Input

  Input_filtered$meta_data <- Input$meta_data[samples,]
  
  for (i in 1:ncol(Input_filtered$meta_data)) {
    if (is.factor(Input_filtered$meta_data[,i])) {
      Input_filtered$meta_data[,i] <- factor(Input_filtered$meta_data[,i]) }
  }

  
  ##**##** TKI/DRC-Specific Code ##**##**
  ## Reorder ET1 Conditions or KI Doses 
  # if ("Project" %in% names(Input_filtered$meta_data)) {
  #   if (all(Input_filtered$meta_data$Project == "DRC")) {
  #     Input_filtered$meta_data$Condition <- factor(Input_filtered$meta_data$Condition, 
  #          levels=c("Unstim","640fM","3.2pM","16pM","80pM","400pM","2nM","10nM")) 
  #     Input_filtered$meta_data$Condition <- factor(Input_filtered$meta_data$Condition)
  #   }
  #   else if (all(Input_filtered$meta_data$Project == "TKI")) {
  #     Input_filtered$meta_data$Dosage <- factor(Input_filtered$meta_data$Dosage, 
  #                                   levels=c("DMSO","1x","3x","3x-r","10x","10x-r"))
  #     Input_filtered$meta_data$Dosage <- factor(Input_filtered$meta_data$Dosage)
  #   }
  # }
  ##**##**#########################**##**
  
  Input_filtered$count_data <- Input$count_data[,samples]
  Input_filtered$count_data <- Input_filtered$count_data[apply(Input_filtered$count_data,1,median)>=minMedian, ]
  GeneNames_Global <<- rownames(Input_filtered$count_data)
  #print(length(GeneNames_Global))
  
  print("END STEP 4: FILTER BY MIN MEDIAN")
  return(Input_filtered)
  
}

##**## STEP 5: FILTER BY USER INPUT (AND GENERATE FILTER UI)
samples_to_analyze <- function(input, x, val){

  print("START STEP 5: FILTER BY USER INPUT")
  
  sample_list <- rownames(input)[which(!input[[x]] %in% val)]
  
  print("END STEP 5: FILTER BY USER INPUT")
  return(sample_list)
}

##**## STEP 5.2: GENERATE UI FOR DYNAMIC VAR FILTERING
make_ui <- function(x, var) {
  levs <- unique(x)
  selectInput(var, var, choices=levs, selected="", multiple = TRUE)
}


#############################################
## Side Panel Options/Downloads  ############
#############################################
### Batch effect removal for clustering 
remove_batcheffect <- function(x, metadata, factors){
  
  if(length(factors)==1){
    factorU1 <- factors
    factorU2 <- ""
  } else if(length(factors)==2){
    factorU1 <- factors[1]
    factorU2 <- factors[2]
  }
  if (factorU2 == ""){
    batchremoved <- removeBatchEffect(assay(x), batch=metadata[,factorU1])
  } else {
    batchremoved <- removeBatchEffect(assay(x),batch=metadata[,factorU1], batch2=metadata[,factorU2])
  }
  
  return(batchremoved)
  
}

### PCA rotation for ALL GENES
pca_rotation <- function(x, PC = "PC1"){
  pca_for_rotation <- prcomp(t(assay(x)))
  #pca_for_rotation <- prcomp(t(x))
  pc_rotation <- sort(abs(pca_for_rotation$rotation[,PC]), decreasing = TRUE)
  write.table(pca_for_rotation$x, "./PCA_Coordinates.txt", col.names = T, row.names = T, sep="\t")
  return(pc_rotation)
}

pca_rotation_batch_removed <- function(x, PC = "PC1"){
  pca_for_rotation <- prcomp(t(x))
  pc_rotation <- sort(abs(pca_for_rotation$rotation[,PC]), decreasing = TRUE)
  write.table(pca_for_rotation$x, "./PCA_Coordinates.txt", col.names = T, row.names = T, sep="\t")
  return(pc_rotation)
}


### PCA rotation for CUSTOM GENES 
custom_pca_rotation <- function(x, PC = "PC1", geneset){
  if (length(grep(",", geneset)>0)) {
    geneset <- unlist(strsplit(geneset, ",")) }
  
  if (length(grep(" ", geneset)>0)) {
    geneset <- unlist(strsplit(geneset, " ")) }
  
  select <- which(rownames(x) %in% geneset)
  pca_for_rotation <- prcomp(t(assay(x[select, ])))
  #pca_for_rotation <- prcomp(t(assay(x)))
  write.table(pca_for_rotation$x, "./PCA_Coordinates.txt", col.names = T, row.names = T, sep="\t")
  pc_rotation <- sort(abs(pca_for_rotation$rotation[,PC]), decreasing = TRUE)

  return(pc_rotation)
}

custom_pca_rotation_batch_removed <- function(x, PC = "PC1", geneset){
  if (length(grep(",", geneset)>0)) {
    geneset <- unlist(strsplit(geneset, ",")) }
  
  if (length(grep(" ", geneset)>0)) {
    geneset <- unlist(strsplit(geneset, " ")) }
  
  select <- which(rownames(x) %in% geneset)
  pca_for_rotation <- prcomp(t(x[select, ]))
  #pca_for_rotation <- prcomp(t(x))
  pc_rotation <- sort(abs(pca_for_rotation$rotation[,PC]), decreasing = TRUE)
  return(pc_rotation)
}



### PCA loading for ALL GENES
pca_loadings <- function(x, ntop){
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca_for_rotation <- prcomp(t(assay(x)[select, ]))
  #pca <- prcomp(t(assay(x)[select, ]), scale = F)
  return(pca_for_rotation$x)
}
pca_loadings_batch_removed <- function(x, ntop){
  rv <- rowVars(x)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca_for_rotation <- prcomp(t(x[select, ]))
  return(pca_for_rotation$x)
}

### PCA loadigs for CUSTOM GENES 
custom_pca_loadings <- function(x, geneset){
  if (length(grep(",", geneset)>0)) {
    geneset <- unlist(strsplit(geneset, ",")) }
  
  if (length(grep(" ", geneset)>0)) {
    geneset <- unlist(strsplit(geneset, " ")) }
  
  select <- which(rownames(x) %in% geneset)
  pca_for_rotation <- prcomp(t(assay(x[select, ])))
  return(pca_for_rotation$x)
}

custom_pca_loadings_batch_removed <- function(x, geneset){
  if (length(grep(",", geneset)>0)) {
    geneset <- unlist(strsplit(geneset, ",")) }
  
  if (length(grep(" ", geneset)>0)) {
    geneset <- unlist(strsplit(geneset, " ")) }
  
  select <- which(rownames(x) %in% geneset)
  pca_for_rotation <- prcomp(t(x[select, ]))
  return(pca_for_rotation$x)
}


########################
## QC Plots ############
#########################

################ NORMALIZATION BOXPLOTS ####################
## For rlog/vst
plot_read_transformation_boxplot <- function(x,y,type){
  par(mfrow=c(3,1))
  boxplot(counts(x), main = "Raw Read Counts", xaxt='n')
  boxplot(log2(counts(x)+1), main = "log2 Normalized Counts", xaxt='n')
  if (type == "vst") {
    boxplot(assay(y), main = "vst Normalized Read Counts", xaxt='n') }
  else if (type == "rlog") {
    boxplot(assay(y), main = "rlog Normalized Read Counts", xaxt='n') }
}

########### PLURITEST #############
run_pluritest <- function(x, intgroup, annotation, returnData = FALSE){
  
  mydata_5 <- as.data.frame(assay(x))
  mydata_5$Symbol <- rownames(mydata_5)
  mydata_6 <- annotation
  
  fData3 <- merge(mydata_5, mydata_6, by="Symbol")
  lengthfdata3 <- length(fData3)-1
  fData4 <- as.matrix(fData3[,2:lengthfdata3])
  rownames(fData4) <- fData3$Illumina_ID
  working.lumi <- ExpressionSet(assayData = fData4)
  sel<-match(rownames(W15),rownames(fData4))
  H15.new<-predictH(exprs(working.lumi[sel,][!is.na(sel),]),W15[!is.na(sel),])
  coef.pt<-c(-1.267095e+02,4.567437e-03, 4.377068e-03,1.043193e-03)
  
  s.new <-drop(coef.pt[1] + coef.pt[2:4]%*%H15.new[c(1,14,13),])
  working.lumi <- lumiN(working.lumi, method =
                          "rankinvariant")
  H12.new<-predictH(exprs(working.lumi[sel,][!is.na(sel),]),
                    W12[!is.na(sel),])
  rss.new<-apply((exprs(working.lumi[sel,][!is.na(sel),])-W12[!is.na(sel),]%*%H12.new)^2,2,sum)
  RMSE.new<-sqrt(rss.new/sum(!is.na(sel)))
  novel.new<-apply((exprs(working.lumi[sel,][!is.na(sel),])-
                      W12[!is.na(sel),]%*%H12.new)^8,2,sum)
  novel.new<-(novel.new/sum(!is.na(sel)))^(1/8)
  mydata <- as.data.frame(cbind(s.new, novel.new))
  
  intgroup.df <- GenomicRanges::as.data.frame(colData(x)[, intgroup, drop=FALSE])
  #print(intgroup.df)
  group <- factor(apply(intgroup.df, 1, paste, collapse = ":"))
  d <- data.frame(Pscore = mydata$s.new, Nscore = mydata$novel.new, group = group, intgroup.df, names = colnames(x))
  
  if (returnData) {
    return(mydata)
  }
  #print(mydata)
  
  #print(intgroup[1])
  
  DoxCategory_CI = c('Resistant' = "midnightblue", 'Sensitive' = "red2")
  
  #print(DoxCategory_CI)
  # 
  if (length(intgroup)==1){
      ggplot(data = d, aes_string(x = "Nscore", y = "Pscore", color = intgroup[1], label = "names"))+
      #geom_text(colour = "black", alpha=0.8, size=7, hjust=0, vjust=0, position = "jitter") +
      geom_point(size = 12) +
      #scale_color_manual(values = c("High Fractality (HF)" = "#F8766D", "No Topography (NT)"= "#00BA38", "Round Topography (RT)" = "#619CFF",
      #                              "hiPSC" = "#D39200", "Immortalized Human Podocytes" = "#93AA00", "Induced Human Podocytes" = "#00C19F",
      #                              "Nephron Progenitor" = "#00B9E3", "Sorted Human Adult Podocytes" = "#DB72FB",
      #                              "Sorted Human Non-Podocytes" = "#FF61C3"))+
      theme(plot.title=element_text(size=20, face="bold"),
            #plot.title = element_blank(),
            #axis.title = element_blank(),
            axis.title=element_text(size=30, face="bold"),
            axis.text= element_text(size=25),
            legend.title=element_text(size=25, face="bold"),
            legend.text = element_text(size=30))+xlab("Novelty")+
      ylab("Pluripotency Score") +
      labs(title = "PluriTest Output")

  }
  
  # if (length(intgroup)==1){
  #   ggplot(data = d, aes_string(x = "Nscore", y = "Pscore", color = intgroup[1], label = "names"))+
  #     #geom_text(colour = "black", alpha=0.8, size=7, hjust=0, vjust=0, position = "jitter") + 
  #     geom_point(size = 12) + 
  #     scale_color_manual(values = c("High Fractality (HF)" = "#F8766D", "No Topography (NT)"= "#00BA38", "Round Topography (RT)" = "#619CFF",
  #                                   "hiPSC" = "#D39200", "Immortalized Human Podocytes" = "#93AA00", "Induced Human Podocytes" = "#00C19F",
  #                                   "Nephron Progenitor" = "#00B9E3", "Sorted Human Adult Podocytes" = "#DB72FB", 
  #                                   "Sorted Human Non-Podocytes" = "#FF61C3"))+
  #     theme(plot.title=element_text(size=20, face="bold"),
  #           axis.title=element_text(size=30, face="bold"),
  #           axis.text= element_text(size=25), legend.position = "none")+xlab("Novelty")+ 
  #     ylab("Pluripotency Score") + 
  #     labs(title = "PluriTest Output")    
  #   
  #} 
  # 
  # if (length(intgroup)==1){
  #   ggplot(data = d, aes_string(x = "Nscore", y = "Pscore", color = intgroup[1], label = "names"))+
  #     geom_text(colour = "black", alpha=0.8, size=5, hjust=0, vjust=0, position = "jitter") + 
  #     geom_point(size = 10) +
  #     theme(axis.text=element_text(size=10))+
  #     xlab("Novelty score")+ 
  #     ylab("Pluripotency score") + 
  #     labs(title = "PluriTest Output")    
  # } 
  
  else if (length(intgroup) == 2){
    #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
    required_shapes <- custom_shape_palette[1:length(unique(d[,intgroup[2]]))]
    ggplot(data = d, aes_string(x = "Nscore", y = "Pscore", color = intgroup[1], shape = intgroup[2], label = "names"))+
      geom_text(colour = "black", alpha = 0.8, size = 7, hjust=0, vjust = 0, position = "jitter") +
      scale_shape_manual(values=required_shapes) +
      geom_point(size = 10) +
      xlab("Novelty score")+ 
      ylab("Pluripotency score") + 
      theme(axis.text=element_text(size=10))+
      labs(title = "PluriTest Output")
  } 
  
  else if (length(intgroup)==3){
    #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
    required_shapes <- custom_shape_palette[1:length(unique(d[,intgroup[2]]))]
    ggplot(data=d,aes_string(x="Nscore",y="Pscore",color=intgroup[1],shape=intgroup[2],alpha=intgroup[3],label="names"))+
      geom_text(colour = "black", alpha = 0.8, size = 5, hjust=0, vjust = 0, position = "jitter") + 
      scale_shape_manual(values=required_shapes) +
      geom_point(size = 10) +
      xlab("Novelty score")+ 
      ylab("Pluripotency score") + 
      theme(axis.text=element_text(size=10))+
      labs(title = "PluriTest Output")
  }
  
}


########## 2D PCA PLOT ####################
plot_PCA_prcomp <- function(x, intgroup, PCa = 1, PCb = 2, ntop, returnData = FALSE, main="PCA plot"){
     
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
	pca <- prcomp(t(assay(x)[select, ]), scale = F)
	percentVar <- pca$sdev^2/sum(pca$sdev^2)

   intgroup.df <- GenomicRanges::as.data.frame(colData(x)[, intgroup, drop=FALSE])
   group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
   d <- data.frame(PCA = pca$x[, PCa], PCB = pca$x[, PCb], 
                   group = group, intgroup.df, 
                   names = colnames(x))
                   #names = colData(x)[["CellLine"]])


    if (length(intgroup)==1) {
     ggplot(data = d, aes_string(x = "PCA", y = "PCB", 
                color = intgroup[1], label = "names"))+
            # scale_color_manual(values = c("High Fractality (HF)" = "#00BA38", "No Topography (NT)"= "#00A9FF", "Round Topography (RT)" = "#DB72FB",
            #                                "hiPSC" = "#D39200", "Immortalized Human Podocytes" = "#93AA00", "Induced Human Podocytes" = "#00C19F",
            #                                "Nephron Progenitor" = "#619CFF", "Sorted Human Adult Podocytes" = "#FF61C3", 
            #                              "Sorted Human Non-Podocytes" = "#FF61C3"))+
      theme(plot.title=element_text(size=20, face="bold"),
            axis.title=element_text(size=30, face="bold"),
            axis.text = element_text(size=25),
            legend.title = element_text(size=25, face = "bold"),
            legend.text = element_text(size=30))+#legend.text = element_text(size=30), legend.position = c(0.75,0.85))
      geom_hline(yintercept = 0, colour = "gray65") + 
      geom_vline(xintercept = 0, colour = "gray65") +
      xlim(min(d$PCA)-2, max(d$PCA)+2)+
      #geom_text(colour="black", alpha=0.8, size=5, hjust=(-0.25), vjust=0, position="jitter") + 
      geom_point(size=12) + 
      xlab(paste0("PC",PCa,": ", round(percentVar[PCa] * 100), "% variance")) + 
      ylab(paste0("PC",PCb,": ", round(percentVar[PCb] * 100), "% variance")) + 
      labs(title=main)
     
     
   } else if (length(intgroup) == 2) {
     #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
     required_shapes <- custom_shape_palette[1:dim(unique(d[intgroup[2]]))[1]]
     ggplot(data = d, aes_string(x = "PCA", y = "PCB", 
                color = intgroup[1],shape = intgroup[2], label = "names")) +
       theme(plot.title=element_text(size=20, face="bold"),
             axis.title=element_text(size=30, face="bold"),
             axis.text = element_text(size = 25),
             legend.title=element_text(size=25),
             legend.text = element_text(size=30))+
      geom_hline(yintercept = 0, colour = "gray65") + 
      geom_vline(xintercept = 0, colour = "gray65") +
      #geom_text(colour = "black", alpha = 0.8, size = 5, hjust=0, vjust=0, position = "jitter") + 
      geom_point(size=12) + 
      scale_shape_manual(values=required_shapes) +
      xlab(paste0("PC",PCa,": ", round(percentVar[PCa] * 100), "% variance")) +
      ylab(paste0("PC",PCb,": ", round(percentVar[PCb] * 100), "% variance")) +
      labs(title=main) #+ scale_color_manual(values=c("#F8766D","#7CAE00","#00BFC4"))
      #+scale_color_brewer()
     
     
   } else if (length(intgroup) == 3) {
     #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
     required_shapes <- custom_shape_palette[1:dim(unique(d[intgroup[2]]))[1]]
     ggplot(data = d, aes_string(x = "PCA", y = "PCB", 
                color = intgroup[1], shape = intgroup[2], 
                alpha = intgroup[3], label = "names")) +
       theme(plot.title=element_text(size=20, face="bold"),
             axis.title=element_text(size=20, face="bold"),
             legend.title=element_text(size=20),
             legend.text = element_text(size=20))+
       geom_hline(yintercept = 0, colour = "gray65") + 
      geom_vline(xintercept = 0, colour = "gray65") +
      geom_text(colour="black", alpha=0.8, size=5, hjust=0, 
                vjust=0, position = "jitter") + 
      geom_point(size=10) + 
      scale_shape_manual(values=required_shapes) +
      xlab(paste0("PC",PCa,": ", round(percentVar[PCa] * 100), "% variance")) + 
      ylab(paste0("PC",PCb,": ", round(percentVar[PCb] * 100), "% variance")) + 
      scale_alpha_discrete(range = c(0.35,0.9)) + 
      labs(title=main)
     
   }
}



############# 2D PCA Batch Removed ####################
plot_PCA_prcomp_batch_removed <- function(x,y, intgroup, PCa=1, PCb=2, ntop="all", returnData=FALSE, main="PCA plot"){
  
  rv <- rowVars(x)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(x[select, ]))
  
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(y)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- GenomicRanges::as.data.frame(colData(y)[, intgroup, drop=FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PCA = pca$x[, PCa], PCB = pca$x[, PCb], 
                  group = group, intgroup.df, 
                  names = colnames(x))
                  #names = colData(y)[["CellLine"]])
  if (returnData) {
    attr(d, "percentVaxr") <- percentVar[PCa:PCb]
    return(d)
  }
  
  if (length(intgroup)==1) {
    ggplot(data=d, aes_string(x="PCA", y="PCB", color=intgroup[1], label="names")) +
      theme(plot.title=element_text(size=20, face="bold"),
            axis.title=element_text(size=20, face="bold"),
            legend.title=element_text(size=20),
            legend.text = element_text(size=20))+
      geom_hline(yintercept = 0, colour = "gray65") + 
      geom_vline(xintercept = 0, colour = "gray65") + 
      xlim(min(d$PCA)-2, max(d$PCA)+2)+
      #geom_text(colour="black", alpha=0.8, size=5, hjust=(-0.25), vjust=0, position="jitter") + 
      geom_point(size=10) + 
      xlab(paste0("PC",PCa,": ", round(percentVar[PCa] * 100), "% variance")) +
      ylab(paste0("PC",PCb,": ", round(percentVar[PCb] * 100),"% variance")) + 
      labs(title=main)
    
    
  } else if (length(intgroup) == 2) {
    #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
    required_shapes <- custom_shape_palette[1:dim(unique(d[intgroup[2]]))[1]]
    ggplot(data=d, aes_string(x="PCA", y="PCB", color=intgroup[1],shape=intgroup[2], label = "names")) +
      theme(plot.title=element_text(size=20, face="bold"),
            axis.title=element_text(size=30, face="bold"),
            axis.text=element_text(size=25),
            legend.title=element_text(size=25),
            legend.text = element_text(size=30))+
      geom_hline(yintercept = 0, colour = "gray65") + 
      geom_vline(xintercept = 0, colour = "gray65") +
      #geom_text(colour = "black", alpha = 0.8, size = 5, hjust=0, vjust=0, position = "jitter") + 
      geom_point(size=12) + 
      scale_shape_manual(values=required_shapes) +
      xlab(paste0("PC",PCa,": ", round(percentVar[PCa] * 100), "% variance")) +
      ylab(paste0("PC",PCb,": ", round(percentVar[PCb] * 100), "% variance")) +
      labs(title=main) #+ scale_color_manual(values=c("#F8766D","#7CAE00","#00BFC4"))
    
    
  } else if (length(intgroup) == 3) {
    #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
    required_shapes <- custom_shape_palette[1:dim(unique(d[intgroup[2]]))[1]]
    ggplot(data=d,aes_string(x="PCA", y="PCB", color=intgroup[1], shape=intgroup[2], alpha=intgroup[3], label="names"))+
      theme(plot.title=element_text(size=20, face="bold"),
            axis.title=element_text(size=20, face="bold"),
            legend.title=element_text(size=20),
            legend.text = element_text(size=20))+
      geom_hline(yintercept = 0, colour = "gray65") + 
      geom_vline(xintercept = 0, colour = "gray65") +
      geom_text(colour = "black", alpha = 0.8, size = 5, hjust=0, vjust=0, position = "jitter") + 
      geom_point(size=10) + 
      scale_shape_manual(values=required_shapes) +
      xlab(paste0("PC",PCa,": ", round(percentVar[PCa] * 100),"% variance")) + 
      ylab(paste0("PC",PCb,": ", round(percentVar[PCb] * 100),"% variance")) + 
      scale_alpha_discrete(range = c(0.35,0.9)) + 
      labs(title=main)
    
  }
}


############# PC ANOVA TABLE #####################
pc_anova_table <- function(x, ntop) {
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]), scale = F)
  
  ## Factor catagorical data in metadata ##
  for (c in 1:ncol(colData(x))) {if (is.character(colData(x)[,c])) {colData(x)[,c]=factor(colData(x)[,c])}}
  
  #### Run ANOVA ####
  pc_anova<-list(); prohoc_all<-list(); prohoc_rownames<-list()
  max.pca <- min(min(dim(pca$x)), 10)
  for (pc in 1:max.pca) {
    for (var in colnames(colData(x))) {
      if (!is.factor(colData(x)[[var]]) | 
          length(unique(colData(x)[[var]])) == 1 | 
          length(unique(colData(x)[[var]])) == nrow(colData(x)))  { next }
      else {
        s <- summary(aov(as.numeric(pca$x[,pc])~colData(x)[[var]]))
        thsd <- TukeyHSD(aov(as.numeric(pca$x[,pc])~colData(x)[[var]]))
        prohoc <- data.frame(thsd$`colData(x)[[var]]`)
        prohoc_all[[var]][pc] <- prohoc["p.adj"]
        pc_anova[[var]][pc] <- s[[1]][["Pr(>F)"]][1]
        prohoc_rownames[[var]] <- rownames(prohoc["p.adj"])
      }
    }
  }
  
  #### Create Table of p Values #####
  pc_anova <- as.data.frame(pc_anova)
  rownames(pc_anova) <- paste0("PC", rownames(pc_anova))
  return(pc_anova)
}


############# PC ANOVA TABLE Batch Removed #####################
pc_anova_table_batch_removed <- function(x, y, ntop) {
  rv <- rowVars(x)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(x[select, ]), scale = F)
  
  #### Run ANOVA ####
  pc_anova<-list(); prohoc_all<-list(); prohoc_rownames<-list()
  max.pca <- min(min(dim(pca$x)), 10)
  for (pc in 1:max.pca) {
    for (var in colnames(colData(y))) {
      if (
          !is.factor(y[[var]]) | 
          length(unique(y[[var]])) == 1 | 
          length(unique(y[[var]])) == nrow(y))  { next }
      else {
        s <- summary(aov(as.numeric(pca$x[,pc])~y[[var]]))
        thsd <- TukeyHSD(aov(as.numeric(pca$x[,pc])~y[[var]]))
        #prohoc <- data.frame(thsd$`y[[var]]`)
        #prohoc_all[[var]][pc] <- prohoc["p.adj"]
        pc_anova[[var]][pc] <- s[[1]][["Pr(>F)"]][1]
        #prohoc_rownames[[var]] <- rownames(prohoc["p.adj"])
      }
    }
  }
  
  #### Create Table of p Values #####
  pc_anova <- as.data.frame(pc_anova)
  rownames(pc_anova) <- paste0("PC", rownames(pc_anova))
  return(pc_anova)
}



############# 3D PCA PLOT #####################
plot3DPCA <- function(x, intgroup, PCa=1, PCb=2, PCc=3, ntop="all", returnData=FALSE, main = "3D PCA plot"){
  
  ## preform pca
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(x)[select, ]), scale = F)
    
  if (!all(intgroup %in% names(colData(x)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  ## create pca table
  pca_tbl <- rbind(pca$rotation, 
                   "Individual.%" = (pca$sdev^2 / sum(pca$sdev^2)*100),
                   "Cumulative.%" = (cumsum(pca$sdev^2/sum(pca$sdev^2)*100)))
  var <- (pca$sdev^2 / sum(pca$sdev^2)*100)
  
  ## create df
  intgroup.df <- GenomicRanges::as.data.frame(colData(x)[, intgroup,drop=F])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  
  d <- data.frame(PCA = pca$x[, PCa], 
                  PCB = pca$x[, PCb], 
                  PCC = pca$x[,PCc], 
                  group = group, intgroup.df, names = colnames(x))
  
  ## store user input factors
  ig1 <- intgroup[1]
  ig2 <- intgroup[2]
  ig3 <- intgroup[3]
  
  ## if 1 factor
  if (length(intgroup) == 1) {
    plot_ly(d, x = ~PCA, y = ~PCB, z = ~PCC, 
            color=x@colData[ig1][[1]], type="scatter3d", mode="markers") %>% 
      layout(scene = list(
        xaxis = list(title = paste('PC',PCa,': ',round(var[PCa],2), '%')), 
        yaxis = list(title = paste('PC',PCb,': ',round(var[PCb],2), '%')),
        zaxis = list(title = paste('PC',PCc,': ',round(var[PCc],2), '%'))))
  }
  
  ## if 2 factors (cant do >2)
  else if (length(intgroup) > 1) {
    plot_ly(d, x = ~PCA, y = ~PCB, z = ~PCC, 
            color=x@colData[ig1][[1]], symbol=x@colData[ig2][[1]], 
            colors = brewer.pal(6,"Set3"),symbols = c(16,15,18,1,0,5), 
            mode = "markers", type="scatter3d") %>% 
    layout(scene = list(
      xaxis = list(title = paste('PC',PCa,': ',round(var[PCa],2), '%')), 
      yaxis = list(title = paste('PC',PCb,': ',round(var[PCb],2), '%')),
      zaxis = list(title = paste('PC',PCc,': ',round(var[PCc],2), '%'))))
  }
}



############# 3D PCA PLOT Batch Removed #####################
plot3DPCA_batch_removed <- function(x,y, intgroup, PCa = 1, PCb = 2, PCc = 3, ntop = "all", returnData = FALSE, main = "3D PCA plot"){
    
  ## preform pca
  rv <- rowVars(x)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(x[select, ]), scale = F)

  ## create pca table
  pca_tbl <- rbind(pca$rotation, 
                   "Individual.%" = (pca$sdev^2 / sum(pca$sdev^2)*100),
                   "Cumulative.%" = (cumsum(pca$sdev^2/sum(pca$sdev^2)*100)))
  var <- (pca$sdev^2 / sum(pca$sdev^2)*100)
  
  ## create df
  intgroup.df <- GenomicRanges::as.data.frame(colData(y)[, intgroup,drop=F])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  
  d <- data.frame(PCA = pca$x[, PCa], 
                  PCB = pca$x[, PCb], 
                  PCC = pca$x[,PCc], 
                  group = group, intgroup.df, names = colnames(x))
  
  ## store user input factors
  ig1 <- intgroup[1]
  ig2 <- intgroup[2]
  ig3 <- intgroup[3]
  
  ## if 1 factor
  if (length(intgroup) == 1) {
    plot_ly(d, x = ~PCA, y = ~PCB, z = ~PCC, 
            #color=x@colData[ig1][[1]], type="scatter3d") %>% 
            color=y@colData[ig1][[1]], type="scatter3d") %>% 
      layout(scene = list(
        xaxis = list(title = paste('PC',PCa,': ',round(var[PCa],2), '%')), 
        yaxis = list(title = paste('PC',PCb,': ',round(var[PCb],2), '%')),
        zaxis = list(title = paste('PC',PCc,': ',round(var[PCc],2), '%'))))
  }
  
  ## if 2 factors (cant do >2)
  else if (length(intgroup) > 1) {
    plot_ly(d, x = ~PCA, y = ~PCB, z = ~PCC, 
            #color=x@colData[ig1][[1]], symbol=x@colData[ig2][[1]],
            color=y@colData[ig1][[1]], symbol=y@colData[ig2][[1]],
            colors = brewer.pal(6,"Set3"),symbols = c(16,15,18,1,0,5), 
            mode = "markers", type="scatter3d") %>% 
      layout(scene = list(
        xaxis = list(title = paste('PC',PCa,': ',round(var[PCa],2), '%')), 
        yaxis = list(title = paste('PC',PCb,': ',round(var[PCb],2), '%')),
        zaxis = list(title = paste('PC',PCc,': ',round(var[PCc],2), '%'))))
  }
}



########## DENDROGRAM ############
plot_hclust <- function(x, dist_method = "euclidean", hclust_method = "ward.D2", main="Hierarchical Clustering Plot", intgroup, ntop){
  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  dissimilarity <- 1 - abs(cor(assay(x[select,]), method="spearman"))
  correlation <- cor(assay(x[select,]), method="spearman")
  distance = as.dist(dissimilarity)
  h <- hclust(dist(t(assay(x[select,])),method=dist_method), hclust_method)

  ig1 <- intgroup[1]
  ig2 <- intgroup[2]

  if (length(intgroup) == 1) {
    dendr <- dendro_data(h, type="rectangle") # convert for ggplot
    idxs <- match(as.character(dendr[["labels"]]$label), rownames(as.data.frame(colData(x))))
    dendr[["labels"]]$Data.Set <- as.data.frame(colData(x))[[ig1]][idxs]
    print(dendr)
    ggplot() + 
      geom_segment(data=dendr$segments, aes(x=x, y=y, xend=xend, yend=yend)) + 
      geom_text(data=dendr$labels, aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
      geom_point(data=dendr$labels, aes(x, y, color=Data.Set), size=4) +
      ylim(-max(dendr$segments$y)*(1/5), max(dendr$segments$y))+
      #coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + ylab("")+
      #ggtitle("All Genes (n=18,851)")+
      theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.line.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            plot.title=element_text(size=20, face="bold", hjust=0.5),
            legend.title=element_text(size=20),
            legend.text = element_text(size=20),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
  }
  
  else if (length(intgroup) == 2) {
    #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
    required_shapes <- custom_shape_palette[1:length(unique(colData(x)[,ig2]))]
    dendr <- dendro_data(h, type="rectangle") # convert for ggplot
    idxs <- match(as.character(dendr[["labels"]]$label), rownames(as.data.frame(colData(x))))
    dendr[["labels"]]$Data.Set <- as.data.frame(colData(x))[[ig1]][idxs]
    dendr[["labels"]]$Data.Set2 <- as.data.frame(colData(x))[[ig2]][idxs]
    ggplot() +
      geom_segment(data=dendr$segments, aes(x=x, y=y, xend=xend, yend=yend)) +
      geom_text(data=dendr$labels, aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
      geom_point(data=dendr$labels, aes(x, y, color=Data.Set, shape=Data.Set2), size=4) +
      ylim(-max(dendr$segments$y)*(1/5), max(dendr$segments$y))+
      scale_shape_manual(values=required_shapes) +
      #coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + ylab("")+
      #ggtitle("All Genes (n=18,851)")+
      theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.line.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            plot.title=element_text(size=20, face="bold", hjust=0.5),
            legend.title=element_text(size=20),
            legend.text = element_text(size=20),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
  }
  
}


########## Batch Removed DENDROGRAM ############
plot_hclust_batch_removed <- function(x, y, dist_method="euclidean", hclust_method="ward.D2", main="Hierarchical Clustering Plot Batch Removed",intgroup, ntop){
  
  rv <- rowVars(x)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  dissimilarity <- 1 - abs(cor(x[select,], method="spearman"))
  correlation <- cor(x[select,], method="spearman")
  distance = as.dist(dissimilarity)
  h <- hclust(dist(t(x[select,]),method=dist_method), hclust_method)
  
  ig1 <- intgroup[1]
  ig2 <- intgroup[2]
  
  if (length(intgroup) == 1) {
    dendr <- dendro_data(h, type="rectangle") # convert for ggplot
    idxs <- match(as.character(dendr[["labels"]]$label), rownames(as.data.frame(colData(y))))
    dendr[["labels"]]$Data.Set <- as.data.frame(colData(y))[[ig1]][idxs]
    #print(dendr[["labels"]])
    ggplot() + 
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
      geom_text(data=label(dendr), aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
      geom_point(data=label(dendr), aes(x, y, color=Data.Set), size=4) +
      ylim(-max(segment(dendr)$y)*(1/5), max(segment(dendr)$y))+
      #coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + ylab("")+
      #ggtitle("All Genes (n=18,851)")+
      theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.line.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            plot.title=element_text(size=20, face="bold", hjust=0.5),
            legend.title=element_text(size=20),
            legend.text = element_text(size=20),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
  }
  
  else if (length(intgroup) == 2) {
    #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
    required_shapes <- custom_shape_palette[1:length(unique(colData(y)[,ig2]))]
    dendr <- dendro_data(h, type="rectangle") # convert for ggplot
    idxs <- match(as.character(dendr[["labels"]]$label), rownames(as.data.frame(colData(y))))
    dendr[["labels"]]$Data.Set <- as.data.frame(colData(y))[[ig1]][idxs]
    dendr[["labels"]]$Data.Set2 <- as.data.frame(colData(y))[[ig2]][idxs]
    ggplot() +
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
      geom_text(data=label(dendr), aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
      geom_point(data=label(dendr), aes(x, y, color=Data.Set, shape=Data.Set2), size=4) +
      ylim(-max(segment(dendr)$y)*(1/5), max(segment(dendr)$y))+
      scale_shape_manual(values=required_shapes) +
      #coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + ylab("")+
      #ggtitle("All Genes (n=18,851)")+
      theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.line.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            plot.title=element_text(size=20, face="bold", hjust=0.5),
            legend.title=element_text(size=20),
            legend.text = element_text(size=20),            
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
  }
}


########## DENDROGRAM - PVCLUST ############
plot_pvclust <- function(x, dist_method="euclidean", hclust_method="ward.D2", intgroup, geneset=NA, ntop, nboots){
  
  if(geneset=="") {
    rv <- rowVars(assay(x))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  }
  
  else if (geneset!="") {
    if (length(grep(",", geneset)>0)) { geneset <- unlist(strsplit(geneset, ",")) }
    if (length(grep(" ", geneset)>0)) { geneset <- unlist(strsplit(geneset, " ")) }
    select <- which(rownames(assay(x)) %in% geneset)
  }
  
  dissimilarity <- 1 - abs(cor(assay(x[select,]), method="spearman"))
  correlation <- cor(assay(x[select,]), method="spearman")
  distance = as.dist(dissimilarity)
  
  h <- hclust(dist(t(assay(x[select,])),method=dist_method), hclust_method)
  #set.seed(86)
  pv <- pvclust(data=assay(x)[select,], method.dist=dist_method, method.hclust=hclust_method, nboot=nboots)
  
  
  ig1 <- intgroup[1]
  ig2 <- intgroup[2]
  
  if (length(intgroup) == 1) {
    TestpvGLOBAL <<- pv
    dendr <- dendro_data(h, type="rectangle") # convert for ggplot
    segs <- dendr$segments[which(dendr$segments$x!=dendr$segments$xend),]
    segs <- segs[which(!duplicated(segs$y)),]
    segs <- segs[order(segs$y),]
    segs$edge.no<-1:nrow(segs); segs$au<-round(pv$edges$au*100,0); segs$bp<-round(pv$edges$bp*100,0)
    idxs <- match(as.character(dendr[["labels"]]$label), rownames(as.data.frame(colData(x))))
    dendr[["labels"]]$Data.Set <- as.data.frame(colData(x))[[ig1]][idxs]
    ggplot() + 
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), alpha=.3) + 
      geom_text(data=label(dendr), aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
      geom_text(data=segs, aes(x, y, label=au), color="red", vjust=-.2, hjust=1.1, size=3) +
      geom_text(data=segs, aes(x, y, label=bp), color="darkgreen", vjust=-.2, hjust=-.1, size=3) +
      geom_point(data=label(dendr), aes(x, y, color=Data.Set), size=4) +
      ylim(-max(segment(dendr)$y)*(1/5), max(segment(dendr)$y))+
      theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.line.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            plot.title=element_text(size=20, face="bold", hjust=0.5),
            legend.title=element_text(size=20),
            legend.text = element_text(size=20),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
  }
  
  else if (length(intgroup) == 2) {
    required_shapes <- custom_shape_palette[1:length(unique(colData(x)[,ig2]))]
    dendr <- dendro_data(h, type="rectangle") # convert for ggplot
    segs <- dendr$segments[which(dendr$segments$x!=dendr$segments$xend),]
    segs <- segs[which(!duplicated(segs$y)),]
    segs <- segs[order(segs$y),]
    segs$edge.no<-1:nrow(segs); segs$au<-round(pv$edges$au*100,0); segs$bp<-round(pv$edges$bp*100,0)
    idxs <- match(as.character(dendr[["labels"]]$label), rownames(as.data.frame(colData(x))))
    dendr[["labels"]]$Data.Set <- as.data.frame(colData(x))[[ig1]][idxs]
    dendr[["labels"]]$Data.Set2 <- as.data.frame(colData(x))[[ig2]][idxs]
    ggplot() +
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), alpha=.3) +
      geom_text(data=label(dendr), aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
      geom_text(data=segs, aes(x, y, label=au), color="red", vjust=-.2, hjust=-.1, size=3) +
      geom_text(data=segs, aes(x, y, label=bp), color="darkgreen", vjust=-.2, hjust=1.1, size=3) +
      geom_point(data=label(dendr), aes(x, y, color=Data.Set, shape=Data.Set2), size=4) +
      ylim(-max(segment(dendr)$y)*(1/5), max(segment(dendr)$y))+
      scale_shape_manual(values=required_shapes) +
      theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.line.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            plot.title=element_text(size=20, face="bold", hjust=0.5),
            legend.title=element_text(size=20),
            legend.text = element_text(size=20),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
  }
}

########## DENDROGRAM - PVCLUST Batch Removed ############
plot_pvclust_batch_removed <- function(x, y, dist_method="euclidean", hclust_method="ward.D2", intgroup, geneset=NA, ntop, nboots){
  
  if(geneset=="") {
    rv <- rowVars(x)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  }
  
  else if (geneset!="") {
    if (length(grep(",", geneset)>0)) { geneset <- unlist(strsplit(geneset, ",")) }
    if (length(grep(" ", geneset)>0)) { geneset <- unlist(strsplit(geneset, " ")) }
    select <- which(rownames(x) %in% geneset)
  }
  
  dissimilarity <- 1 - abs(cor(x[select,], method="spearman"))
  correlation <- cor(x[select,], method="spearman")
  distance = as.dist(dissimilarity)
  
  h <- hclust(dist(t(x[select,]),method=dist_method), hclust_method)
  pv <- pvclust(data=x[select,], method.dist=dist_method, method.hclust=hclust_method, nboot=nboots)
  
  ig1 <- intgroup[1]
  ig2 <- intgroup[2]
  
  if (length(intgroup) == 1) {
    dendr <- dendro_data(h, type="rectangle") # convert for ggplot
    segs <- dendr$segments[which(dendr$segments$x!=dendr$segments$xend),]
    segs <- segs[which(!duplicated(segs$y)),]
    segs <- segs[order(segs$y),]
    segs$edge.no<-1:nrow(segs); segs$au<-round(pv$edges$au*100,0); segs$bp<-round(pv$edges$bp*100,0)
    idxs <- match(as.character(dendr[["labels"]]$label), rownames(as.data.frame(colData(y))))
    dendr[["labels"]]$Data.Set <- as.data.frame(colData(y))[[ig1]][idxs]
    ggplot() + 
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), alpha=.3) + 
      geom_text(data=label(dendr), aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
      geom_text(data=segs, aes(x, y, label=au), color="red", vjust=-.2, hjust=-.1, size=3) +
      geom_text(data=segs, aes(x, y, label=bp), color="darkgreen", vjust=-.2, hjust=1.1, size=3) +
      geom_point(data=label(dendr), aes(x, y, color=Data.Set), size=4) +
      ylim(-max(segment(dendr)$y)*(1/5), max(segment(dendr)$y))+
      theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.line.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            plot.title=element_text(size=20, face="bold", hjust=0.5),
            legend.title=element_text(size=20),
            legend.text = element_text(size=20),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
  }
  
  else if (length(intgroup) == 2) {
    required_shapes <- custom_shape_palette[1:length(unique(colData(y)[,ig2]))]
    dendr <- dendro_data(h, type="rectangle") # convert for ggplot
    segs <- dendr$segments[which(dendr$segments$x!=dendr$segments$xend),]
    segs <- segs[which(!duplicated(segs$y)),]
    segs <- segs[order(segs$y),]
    segs$edge.no<-1:nrow(segs); segs$au<-round(pv$edges$au*100,0); segs$bp<-round(pv$edges$bp*100,0)
    idxs <- match(as.character(dendr[["labels"]]$label), rownames(as.data.frame(colData(y))))
    dendr[["labels"]]$Data.Set <- as.data.frame(colData(y))[[ig1]][idxs]
    dendr[["labels"]]$Data.Set2 <- as.data.frame(colData(y))[[ig2]][idxs]
    ggplot() +
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), alpha=.3) +
      geom_text(data=label(dendr), aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
      geom_text(data=segs, aes(x, y, label=au), color="red", vjust=-.2, hjust=-.1, size=3) +
      geom_text(data=segs, aes(x, y, label=bp), color="darkgreen", vjust=-.2, hjust=1.1, size=3) +
      geom_point(data=label(dendr), aes(x, y, color=Data.Set, shape=Data.Set2), size=4) +
      ylim(-max(segment(dendr)$y)*(1/5), max(segment(dendr)$y))+
      scale_shape_manual(values=required_shapes) +
      theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.line.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            plot.title=element_text(size=20, face="bold", hjust=0.5),
            legend.title=element_text(size=20),
            legend.text = element_text(size=20),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
  }
  
}






########## Heatmaps ############
ann_colors = list(
#   Phenotype = c(Normal = "darkgreen", Hypertrophy = "darkred", Unknown = "gray"),
#   CellLine = c(Nor_A = "darkolivegreen3", Nor_B = "darkseagreen3", Nor_C="green3",
#                Hyp_A = "firebrick1", Hyp_B = "indianred", Hyp_C = "orangered", FHM = "gray"),
  Gender = c(Male = "blue", Female = "pink", Unknown = "ghostwhite"),
  Phenotype = c(Normal = "green", Hypertrophy = "brown", Unknown = "ghostwhite"),
  CellLine = c(Nor_A = "limegreen", Nor_B = "darkseagreen3", Nor_C="seagreen",
               Hyp_D = "burlywood3", Hyp_E = "bisque", Hyp_F = "darksalmon", FHM = "ghostwhite"),
  Condition = c(ET1 = "dodgerblue4", Unstim = "darkslategray"))


########## Heatmap Top 100 Most Var ############
plot_gene_cluster100 <- function(x, n = 100, main="Most variable genes", intgroup){
  topVarGenes <- head(order(-rowVars(assay(x))), n)
  mat <- assay(x)[ topVarGenes, ]
  mat <- mat - rowMeans(mat)
  ig1 = intgroup[1]
  ig2 = intgroup[2]
  
  n <- length(unique(colData(x)[,ig1]))
  cols <- list(hcl(h=seq(15, 375, length=n+1), l=65, c=100)[1:n])
  names(cols)[1] <- ig1
  names(cols[[1]]) <- levels(as.factor(colData(x)[,ig1]))

  if (length(intgroup) == 1) {
    df <- as.data.frame(colData(x)[,ig1])
    colnames(df) <- ig1
    rownames(df) <- rownames(colData(x))
    pheatmap(mat, annotation_col=df, annotation_colors=cols)
    #pheatmap(mat, annotation_col=df)
    }
  else {
    df <- as.data.frame(colData(x)[,c(ig1, ig2)])
    pheatmap(mat, annotation_col=df, annotation_colors=cols)
    #pheatmap(mat, annotation_col=df)
    }
}

########## Heatmap Top 100 Batch Removed Most Var ############
plot_gene_cluster100_batch_removed <- function(x, y, n=100, intgroup, main="Most variable genes after Batch effect removal"){
  
  topVarGenes <- head(order(-rowVars(x)), n)
  mat <- x[ topVarGenes, ]
  mat <- mat - rowMeans(mat)
  ig1 = intgroup[1]
  ig2 = intgroup[2]
  
  n <- length(unique(colData(y)[,ig1]))
  cols <- list(hcl(h=seq(15, 375, length=n+1), l=65, c=100)[1:n])
  names(cols)[1] <- ig1
  names(cols[[1]]) <- levels(as.factor(colData(y)[,ig1]))

  if (length(intgroup) == 1) {
    df <- as.data.frame(colData(y)[,ig1])
    colnames(df) <- ig1
    rownames(df) <- rownames(colData(y))
    pheatmap(mat, annotation_col = df, annotation_colors = cols)
    #pheatmap(mat, annotation_col=df)   
    }
  else {
    df <- as.data.frame(colData(y)[,c(ig1, ig2)])
    pheatmap(mat, annotation_col = df, annotation_colors = cols)
    #pheatmap(mat, annotation_col=df)   
    }
}
  



########## Heatmap Top 500 Most Var ############
plot_gene_cluster500 <- function(x, n = 500, main="Most variable genes", intgroup){
  topVarGenes <- head(order(-rowVars(assay(x))), n)
  mat <- assay(x)[ topVarGenes, ]
  mat <- mat - rowMeans(mat)
  ig1 = intgroup[1]
  ig2 = intgroup[2]
  
  n <- length(unique(colData(x)[,ig1]))
  cols <- list(hcl(h=seq(15, 375, length=n+1), l=65, c=100)[1:n])
  names(cols)[1] <- ig1
  names(cols[[1]]) <- levels(as.factor(colData(x)[,ig1]))

  if (length(intgroup) == 1) {
    df <- as.data.frame(colData(x)[,ig1])
    colnames(df) <- ig1
    rownames(df) <- rownames(colData(x))
    pheatmap(mat, annotation_col = df, fontsize_row = 7, annotation_colors = cols)
    #pheatmap(mat, annotation_col=df, fontsize_row=7)   
    }
  else {
    df <- as.data.frame(colData(x)[,c(ig1, ig2)])
    pheatmap(mat, annotation_col = df, fontsize_row = 7, annotation_colors=cols)
    #pheatmap(mat, annotation_col = df, fontsize_row = 7)
  }
}

########## Heatmap Top 500 Batch Removed Most Var ############
plot_gene_cluster500_batch_removed <- function(x, y, n=500, intgroup, main="Most variable genes after Batch effect removal"){
  topVarGenes <- head(order(-rowVars(x)), n)
  mat <- x[ topVarGenes, ]
  mat <- mat - rowMeans(mat)
  ig1 = intgroup[1]
  ig2 = intgroup[2]
  
  n <- length(unique(colData(y)[,ig1]))
  cols <- list(hcl(h=seq(15, 375, length=n+1), l=65, c=100)[1:n])
  names(cols)[1] <- ig1
  names(cols[[1]]) <- levels(as.factor(colData(y)[,ig1]))
  
  if (length(intgroup) == 1) {
    df <- as.data.frame(colData(y)[,ig1])
    colnames(df) <- ig1
    rownames(df) <- rownames(colData(y))
    pheatmap(mat, annotation_col = df, fontsize_row = 7, annotation_colors = cols)
    #pheatmap(mat, annotation_col = df, fontsize_row = 7)
  }
  else {
    df <- as.data.frame(colData(y)[,c(ig1, ig2)])
    pheatmap(mat, annotation_col = df, fontsize_row = 7, annotation_colors = cols)
    #pheatmap(mat, annotation_col = df, fontsize_row = 7)
  }
}



########## Heatmap HouseKeeping Genes ############
plot_housekeeping <- function(x){
  housekeepinggenes <- c("GAPDH","B2M","C1orf43","CHMP2A","EMC7","GPI","PSMB2","PSMB4","RAB7A","REEP5","SNRPD3","VCP","VPS29")
  hm.rld <- round(t(assay(x)[ which(rownames(assay(x)) %in% housekeepinggenes), ]),1)
  #hm.rld <- round(t(counts(x, normalized=TRUE)[ housekeepinggenes, ]),1)
  heatmap.2(t(assay(x)[ which(rownames(assay(x)) %in% housekeepinggenes), ]), 
            trace="none", cellnote = hm.rld, notecex = 0.8, 
            notecol = "black", margin = c(5,15))
}



##### QC GENES #######
plot_qcgenes <- function(x, intgroup){
  qc_genes <- c("TNC","VCAN","PCOLCE","FOS","SEMA5A","SNAI2","SERPINE2","COL6A2","THBS1","MYC",
                "PODXL","KLF4","SOX2","LOX","DCN","FN1","FBN1","GREM1","COL3A1")
  mat <- assay(x)[ which(rownames(assay(x)) %in% qc_genes), ]
  mat <- mat - rowMeans(mat)
  ig1 = intgroup[1]
  ig2 = intgroup[2]
  if (length(intgroup) == 1) {
    df <- as.data.frame(colData(x)[,ig1])
    colnames(df) <- ig1
    rownames(df) <- rownames(colData(x))
    pheatmap(mat, annotation_col = df)
  }
  else {
    df <- as.data.frame(colData(x)[,c(ig1, ig2)])
    pheatmap(mat, annotation_col = df)
  }
}

plot_qcgenes_batchremoved <- function(x,y,intgroup){
  qc_genes <- c("TNC","VCAN","PCOLCE","FOS","SEMA5A","SNAI2","SERPINE2","COL6A2","THBS1","MYC",
                "PODXL","KLF4","SOX2","LOX","DCN","FN1","FBN1","GREM1","COL3A1")
  mat <- x[ which(rownames(x) %in% qc_genes), ]
  mat <- mat - rowMeans(mat)
  ig1 = intgroup[1]
  ig2 = intgroup[2]
  if (length(intgroup) == 1) {
    df <- as.data.frame(colData(y)[,ig1])
    colnames(df) <- ig1
    rownames(df) <- rownames(colData(y))
    pheatmap(mat, annotation_col = df)
  }
  else {
    df <- as.data.frame(colData(y)[,c(ig1, ig2)])
    pheatmap(mat, annotation_col = df)
  }
}





##################################
#### HEATMAP Gene Sets ###########
###################################
get_genesets <- function(setname) {
  
  if (setname=="Cardiac Hypertrophy Marker Genes") {
    return(c("ACTA1","ACTA2","ACTB","ACTC1","ACTG1","ATP2A2","CORIN","FOS","GATA4","MYH3","MYH4","MYH6","MYH7","MYH7B","MYH9","NFATC3","NFATC4","NPPA","NPPB","NIN","GSK3B")) }
  
  if (setname=="HyperGEN Linkage Genes") {
    return(c('ABHD1','ACCS','AIM1L','ANGEL1','ARFRP1','ARHGAP25','ARHGEF3','ARSK','ATF7IP2','ATP8A1','C14orf105','C15orf62','C2orf73','C5orf52','C6orf89','C4orf41','CA3','CCR3','CNTN5','COL28A1','CRYBG3','CSPG5','CX3CR1','CYLC2','DCLK1','DDRGK1','DIXDC1','DPRX','DSCAML1','DUOXA1','EBI3','EDN1','EWSR1','FAM136A','FAM151A','FAM183A','FAM219A','FGFR1','FRMD3','GCA','GGNBP2','GLTSCR1','GPD1','GUCA2A','HAUS6','HAX1','HOXA1','HOXD9','HSPA12B','HYAL1','IGFBP2','ITGA10','ITIH4','KCNK2','KLF14','LPA','LPAR2','LRP5L','MED13L','MSRB2','MYRIP','NAT6','NDUFB5','NEURL3','NLRP3','NPRL2','NT5C1B-RDH14;RDH14','NTAN1','OBP2B','OR10A3','OR1I1','OR52A1','OR52D1','OR9A4','PCED1B','PDE12','PDE1B','PDS5A','PEX11G','PHOSPHO1','PIPOX','PLA2G7','PLEKHB2','PLEKHH1','PLXDC2','PPM1M','PPP5C','PRR29','PSG7','PTBP2','PTGER3','PYGO1','RAMP1','REPIN1','RPL14','RSRC1','SALL2','SBNO2','SERPINB11','SETD3','SF3B1','SFXN5','SHISA3','SLC27A6','SLC35C1','SLC39A7','SLC4A8','SP140','STAB1','TBCK','TCFL5','TCIRG1','TCL1A','TMC1','TNFRSF4','TNNI2','TRIM60','TUBGCP3','UBXN2A','UCP2','WSB2','ZBTB44','ZNF200','ZNF644','ZNF727','ZNF778','ZNF83','ZSWIM1')) }
  
  if (setname=="Colagen Genes") {
    return(c('COL10A1','COL11A1','COL11A2','COL12A1','COL13A1','COL14A1','COL15A1','COL16A1','COL17A1',
             'COL18A1','COL19A1','COL1A1','COL1A2','COL20A1','COL21A1','COL22A1','COL23A1','COL24A1','COL25A1','COL27A1','COL28A1','COL2A1','COL3A1','COL4A1','COL4A2','COL4A3','COL4A3BP','COL4A4','COL4A5',
'COL4A6','COL5A1','COL5A2','COL5A3','COL6A1','COL6A2','COL6A3','COL6A4P2','COL6A5','COL6A6','COL7A1'
,'COL8A1','COL8A2','COL9A1','COL9A2','COL9A3')) }
  
  if (setname=="Laminin Genes") {
    return(c('LAMA1','LAMA2','LAMA3','LAMA4','LAMA5','LAMB1','LAMB2','LAMB2P1','LAMB3','LAMB4','LAMC1','LAMC2','LAMC3','LAMP1','LAMP2','LAMP3','LAMP5'))}
  
  if (setname=="Integrin Genes") {
    return(c('ITGA1','ITGA10','ITGA11','ITGA2','ITGA2B','ITGA3','ITGA4','ITGA5','ITGA6','ITGA7','ITGA8','ITGA9','ITGAD','ITGAE','ITGAL','ITGAM','ITGAV','ITGAX','ITGB1','ITGB1BP1','ITGB1BP2','ITGB2','ITGB3','ITGB3BP','ITGB4','ITGB5','ITGB6','ITGB7','ITGB8','ITGBL1'))  }
  
  if (setname=="MAPK Genes") {
    return(c('MAP2K4P1','MAP2K1','MAP2K2','MAP2K3','MAP2K4','MAP2K5','MAP2K6','MAP2K7','MAPK1'))  }
  
  if (setname=="PIK3 Genes") {
    return(c('PIK3AP1','PIK3C2A','PIK3C2B','PIK3C2G','PIK3C3','PIK3CA','PIK3CB','PIK3CD','PIK3CG','PIK3IP1','PIK3R1','PIK3R2','PIK3R3','PIK3R4','PIK3R5','PIK3R6')) }
    
  if (setname=="Rac/Ras Genes") {
    return(c('PTK2','PXN','RAC1','RAC2','RAC3','RAF1','RASA1','RASA2','RASA3','RASA4','TLN1','TLN2','VCL','FBN1','FBN2','FBN3','FN1','GRB2')) }
  
  if (setname=="ECM Genes") {
    return(c('COL10A1','COL11A1','COL11A2','COL12A1','COL13A1','COL14A1','COL15A1','COL16A1','COL17A1','COL18A1','COL19A1','COL1A1','COL1A2','COL20A1','COL21A1','COL22A1','COL23A1','COL24A1','COL25A1','COL27A1','COL28A1','COL2A1','COL3A1','COL4A1','COL4A2','COL4A3','COL4A3BP','COL4A4','COL4A5','COL4A6','COL5A1','COL5A2','COL5A3','COL6A1','COL6A2','COL6A3','COL6A4P2','COL6A5','COL6A6','COL7A1','COL8A1','COL8A2','COL9A1','COL9A2','COL9A3','FBN1','FBN2','FBN3','FN1','GRB2','ITGA1','ITGA10','ITGA11','ITGA2','ITGA2B','ITGA3','ITGA4','ITGA5','ITGA6','ITGA7','ITGA8','ITGA9','ITGAD','ITGAE','ITGAL','ITGAM','ITGAV',
  'ITGAX','ITGB1','ITGB1BP1','ITGB1BP2','ITGB2','ITGB3','ITGB3BP','ITGB4','ITGB5','ITGB6','ITGB7','ITGB8','ITGBL1','LAMA1','LAMA2','LAMA3','LAMA4','LAMA5','LAMB1','LAMB2','LAMB2P1','LAMB3','LAMB4','LAMC1','LAMC2','LAMC3','LAMP1','LAMP2','LAMP3','LAMP5','MAP2K4P1','MAP2K1','MAP2K2','MAP2K3','MAP2K4','MAP2K5','MAP2K6','MAP2K7','MAPK1','PIK3AP1','PIK3C2A','PIK3C2B','PIK3C2G','PIK3C3','PIK3CA','PIK3CB'
,'PIK3CD','PIK3CG','PIK3IP1','PIK3R1','PIK3R2','PIK3R3','PIK3R4','PIK3R5','PIK3R6','PTK2','PXN','RAC1','RAC2','RAC3','RAF1','RASA1','RASA2','RASA3','RASA4','TLN1','TLN2','VCL','COL6A4P1')) }
  
  if (setname=="QC Marker Genes") {
    return(c("TNC","VCAN","PCOLCE","FOS","SEMA5A","SNAI2","SERPINE2","COL6A2","THBS1","MYC","PODXL","KLF4","SOX2","LOX","DCN","FN1","FBN1","GREM1","COL3A1"))  }
  
  if (setname == "C.Gu VR Genes") {
    return(c('A4GALT','ABCA2','AC004057.1','AC007292.3','AC009531.2','AC074091.13','AC093818.1','AC093850.2','ACSS1','ACVRL1','ADAMTS16','ADM2','AKT1','AMD1','AMPD2','ANAPC10','ANKMY2','ANO10','ANP32E','AP000640.10','AP001062.7','AP2S1','APOBEC2','ARAP1','ARHGAP32','ATL3','ATP2A2','ATP2C1','ATP6V0E1','AURKAIP1','B3GAT2','BANF1','BMP5','BRI3','C3','C6orf120','C7orf60','CACNA1G-AS1','CACNA2D3','CAPZB','CASKIN2','CCDC124','CCDC81','CCL4','CDC37','CDH1','CDYL2','CEP85L','CFLAR','CHCHD10','CHCHD6','CHL1','CNR1','COL4A2','COLEC12','COMT','CPLX2','CRIM1','CRLF1','CTC-471F3.4','CTC-523E23.1','CUEDC1','CXCR2','CXorf38','CYP1B1','CYP27C1','DACT2','DCC','DDR1','DDX50P1','DNAJB14','DNM3OS','DPP4','DYRK4','EBF1','EDEM2','EEF1A1P9','EFEMP2','EIF3J','EIF3KP1','EPHA4','EPHA6','EPM2AIP1','EPN1','EPS8L1','ERC1','EREG','EXOC5P1','FAM115C','FAM26F','FAP','FAT2','FENDRR','FEZ2','FOXL1','FSTL4','FTH1','GABRA4','GABRA5','GADD45GIP1','GALNT2','GAS5','GATA3','GCAT','GGA1','GJA1P1','GPAA1','GPC1','GUSB','HCN1','HLA-DRB5','HNRNPA3','HOXA5','HOXC4','HOXD3','HPN','HPR','HRC','IL15','IL1RL1','ISLR2','ITIH3','JPH4','KCNB2','KCNIP1','KCNJ16','KCNK12','KHSRPP1','KIAA2013','KLHL23','KLHL30','KRAS','LASP1','LGALSL','LIFR','LINC00662','LINC00689','LINGO2','LINS','LMOD3','LMTK3','LPCAT2','LRFN2','LRRC73','M6PR','MAGEF1','MED30','MESTP1','MFGE8','MGAT5B','MKLN1','MOV10L1','MRPL47','MRPL52','MRPL54','MT-ND4','MTHFR','MTRF1','MTRNR2L11','MYEF2','MYPOP','NAA25','NBPF20','NCAN','NCR3LG1','NECAP1','NGFR','NOP56','NPFFR2' ,'NPM1','NQO2','NT5C1A','OGN','OPN4','OPRK1','OR8S1','OST4','OSTC','P4HA3','PAPD4','PCSK1','PDCD5','PDCL3','PEAR1','PKD1','PKD1L2','PLEKHH1','POLR2I','POLR2L','PPAPDC1A','PPDPF','PRDX2P4','PRDX6','PREPL','PRRG3','PSKH1','PTPRN2','PTPRZ1','RAB3IP','RANBP3L','RBM28','RCN3','REXO1','RFX4','RHOJ','RIMKLA','RN7SKP130','RN7SL18P','RN7SL197P','RN7SL2','RN7SL424P','RN7SL483P','RN7SL635P','RN7SL711P','RP1-122P22.2','RP11-1212A22.1','RP11-121L10.3','RP11-174O3.1','RP11-254B13.3','RP11-301G21.1','RP11-380D23.1','RP11-385J23.1','RP11-39C10.1','RP11-422J8.1','RP11-428G5.2','RP11-554D14.2','RP11-64D22.2','RP11-69L16.3','RP11-834C11.4','RP5-857K21.7','RPF1','RPL11','RPL12','RPL13P12','RPL18A','RPL23A','RPL24P4','RPL31P2','RPL9','RPS15A','RPS24','RPS3A','RPS3AP26','RPS3AP6','RSL24D1','RUFY3','S100A13','S100A6','SART1','SC5D','SCARNA6','SCARNA7','SCD5','SEMA5B','SETX','SGCA','SGSM1','SLC16A12','SLC22A18','SLC25A5P2','SLC7A2','SMYD1','SNAPC3','SNHG1','SNRPG','SOX2','SPEG','SPON2','SPRYD7','SRCRB4D','SRGN','SSTR1','ST8SIA1','SYT17','TAF1D','TAGLN2P1','TIMM13','TKT','TMEM130','TMEM151B','TMEM200A','TMEM259','TMEM64','TOMM7','TPT1','TRMT6','TRPC6','TXN','UPK1A','UQCRB','UTP11L','WDR74','WI2-1896O14.1','YIPF5','ZBTB10','ZIC2','ZNF124','ZNF320','ZNF444','ZNF611','ZNF701','ZNRF2P2')) }
  
  if (setname == "C.Gu Pure DE Genes") {
    return(c('ABCA6','AC073283.4','ACPP','ADAM21','ADCY8','AJAP1','AMHR2','AMPH','ANGPT1','ANKRD2','AP001258.4','AREG','ARTN','BMP10','BTBD8','C10orf82','C16orf89','C19orf57','C1orf127','C5','CACNG4','CALN1','CD200','CDC45','CENPH','CHAF1B','CLDN4','CLSPN','CMPK2','CNRIP1','CTC-558O2.1','DRD1','DSCC1','DUSP4','EXO1','FGF12-AS3','FOSL1','GDF6','GINS3','GPR61','HELLS','HFE','HIST1H3J','HK2','IFIH1','IL17B','IMPA2','KB-173C10.1','KCNN3','KRT18P38','LHX1','LIG1','LINC00595','LINC00648','LINC00707','LINC00954','LIPG','MAP2K6','MCM10','METTL7A','MFSD4','MRAP2','NECAB1','NFE2','NPY4R','OPRL1','ORC1','P2RY6','PAK6','PDE1A','PLCB2','PLCD4','PLSCR4','POF1B','PRIM1','PTGS1','RAD51','RAD51B','RIBC1','RMI2','RP11-13K12.1','RP11-161D15.3','RP11-366M4.3','RP11-400K9.4','RP11-554I8.2','RP11-647K16.1','RPS6KA5','SEPP1','SERPINB2','SFRP5','SLC16A9','SLC9A9','SYT10','TMEM179','TMEM180','TMEM254-AS1','UBASH3B','UHRF1','ZIC1','ZNF815P','AC138623.1','ACAN','ACOXL','ACTN3','AGAP4','AK5','AP000350.5','AP006222.2','ATP1A4','BEGAIN','BST2','C6orf223','CECR1','COL6A4P2','CPNE7','CSPG4P13','CTA-390C10.10','CTB-47B11.3','CTD-2291D10.4','CTD-2651B20.3','CYP4V2','DOCK10','EBF2','EEF1GP1','EIF1AY','ELAVL3','F10','GIPC3','GLDC','GSTT2B','GYG2P1','HCG4P3','HEYL','HPR','HSD3B2','KALP','KLHL1','KRT8P8','LGALS4','LINC00960','LL0XNC01-116E7.4','MAB21L1','MAB21L2','MTATP6P1','MTATP8P2','OR7D2','OSTN','PCDHB8','PCSK1','PLEKHG7','PPIAP29','PRKY','PRSS45','PRSS46','PRSS50','PSORS1C1','RBM44','RELL1','RFTN2','RN7SL471P','RP11-165H4.2','RP11-21I10.2','RP11-249L21.4','RP11-256I23.3','RP11-326I11.4','RP11-356B19.11','RP11-368L12.1','RP11-379P15.1','RP11-428C19.4','RP11-428C19.5','RP11-455F5.3','RP11-567P17.1','RP11-706O15.5','RP11-734I18.1','RP11-893F2.14','RP11-94I2.1','RP3-323K23.3','RPL5P25','RPS4Y1','SCG3','SDPR','SETD9','SHISA9','SLC25A47','SPTLC3','SULT1A1','SUSD4','TDGF1','TDGF1P3','TTTY10','TTTY15','TXLNG2P','ULK4P3','USP9Y','UTY','ZFY','ZNF285','ZNF585B','ZNF680','ZNF730',"hsa-mir-6723")) }
    
  if (setname=="Greber's QC Marker Genes") {
    return(c("SOX2","NANOG","POU5F1","IRX3","MYL7","MYL2","MYL4","MYH6","MYH7","ACTC1"))}
      
  if (setname=="Greber's Atrial Marker Genes") {
    return(c("NPPB","HAMP","MYBPHL","NPPA","PLA2G2A","COMP","TCEAL2","SLPI","DHRS9","HP")) }

  if (setname=="Greber's Ventricular Marker Genes") {
    return(c("DLK1","IRX4","MYL2","XDH","TMEM190","HYLA2","CPNE4","CYP1A1","IRX5","C3orf23"))}
  
  if (setname=="Greber's Atrial / Ventricular Marker Genes") {
    return(c("NPPB","HAMP","MYBPHL","NPPA","PLA2G2A","COMP","TCEAL2","SLPI","DHRS9","HP","DLK1","IRX4","MYL2","XDH","TMEM190","HYLA2","CPNE4","CYP1A1","IRX5","C3orf23")) }
  
  if (setname=="Greber's Early iPSC-CM Marker Genes") {
    return(c("NKX2-5","IRX4","TBX2","COL2A1","ISL1","HAND1","ID2","LEF1","IRS1","MDK")) }
    
  if (setname=="Greber's Late iPSC-CM Marker Genes") {
    return(c("MYH7","MYL2","TNNI3K","HSPB7","PLN","CSRP3","ACTN2","RBM20","TRIM63","CORIN")) }
      
  if (setname=="Greber's Up-Regulated Cardio Marker Genes") {
    return(c("MYL2","TNNI3","ACTN2","MYH7","MYL3","TNNC1","S100A1","TNNT2","MYH11","SORBS1","ATP2A2","PLN","CASQ2","RYR2","RYR3","TRDN","ITPR1","ITPR3","ASPH","HRC","KCNA4","KCNA5","KCNAB1","KCNAB2","KCND2","KCND3","KCNE4","KCNG1","KCNH2","KCNH7","KCNIP2","KCNJ2","KCNJ3","KCNJ5","KCNJ8","KCNK1","KCNQ1","KCNV1","SCN1A","SCN1B","SCN2B","SCN3A","SCN4B","SCB5A","HCN1","HCN4","CACNA1C","CACNA1D","CACNA1H","CACNA1G","CACNA2D1","CACNB2","SLC8A1","TRPC3","TRPC4","TRPC6","CFTR")) }
  
  if (setname=="Wu's Cardio Matration Marker Genes") {
    return(c("TNNI3","ACTN2","MYH7","RYR2","ITPR3","KCNH2","MYH6","CAV3","TLR3","HCN4","MYPN")) }
  
  if (setname=="DB C57 Marker Genes") {
    return(c("CNTF","CNTFR","OXT","OXTR","LEP","LEPR","STAT3","LNPEP","SLC2A4","GLRX3","LDLR","MTTP","SORBS1","NAMPT","GPIHBP1","PTPLB","PTPN11","LIAS","NDUFS6","ENDOG","ATPAF2","AIFM2","PRKAA2","APC","DES","ILK","TTN","MYL2","DAPK3","DMD","GSN","COL1A2","LMNA","LAMA2","CANX","MYOM1","CACNB2","SMAD3","TGFB1","TSP1","MMP1","MMP2","MMP3","MMP8","MMP9","MMP12","MMP13","MMP14","MMP28","TIMP1","TIMP2","TIMP3","TIMP4","ANGPT2")) }
  
      
}





#######################################
#### Custom Gene Set  Tab #############
#######################################

############### 2D PCA Custom Gene Set ###########################
plot_PCA_prcomp_custom <- function(x, intgroup, PCa=1, PCb=2, returnData=F, geneset, selectset, predefined){
  
  if (predefined==F) {
    if (length(grep(",", geneset)>0)) {
      geneset <- unlist(strsplit(geneset, ",")) }
    if (length(grep(" ", geneset)>0)) {
      geneset <- unlist(strsplit(geneset, " ")) }
  }
  
  else if (predefined==T) {
    geneset <- get_genesets(selectset)
  }
  
  select <- which(rownames(x) %in% geneset)
  pca <- prcomp(t(assay(x)[select, ]), scale = F)
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  intgroup.df <- GenomicRanges::as.data.frame(colData(x)[, intgroup, drop=FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  #group <- reorder(group, new.order=c(8,6,4,2,7,5,3,1))
  d <- data.frame(PCA = pca$x[, PCa], PCB = pca$x[, PCb], 
                  group = group, intgroup.df, 
                  names = colnames(x))
  #main=paste0("PC",PCa," vs. PC", PCb," using ",length(geneset)," Genes from Custom Gene Set")
  #main=paste0("PC",PCa," vs. PC", PCb," using ",nrow(assay(x)[select,])," Genes from Custom Gene Set")
  main = "Custom gene PCA plot"
  
  if (length(intgroup)==1) {
    ggplot(data = d, aes_string(x = "PCA", y = "PCB", 
                                color = intgroup[1], label = "names"))+
      geom_hline(yintercept = 0, colour = "gray65") + 
      geom_vline(xintercept = 0, colour = "gray65") +
      #geom_text(colour = "black", alpha = 0.8, size = 5,hjust=0, vjust=0, position = position_nudge(x=0.03)) + 
      geom_point(size=12) + 
      theme(plot.title = element_text(size = 20, face = "bold"),
            axis.title = element_text(size = 30, face = "bold"),
            axis.text = element_text(size = 25),
            legend.title = element_text(size = 25),
            legend.text = element_text(size = 20))+
      xlab(paste0("PC",PCa,": ", round(percentVar[PCa] * 100), "% variance")) + 
      ylab(paste0("PC",PCb,": ", 
                  round(percentVar[PCb] * 100), "% variance")) + 
      labs(title=main)
      #+scale_color_brewer(type="seq", direction=1, palette=1)
    
    
  } else if (length(intgroup) == 2) {
    #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
    required_shapes <- custom_shape_palette[1:dim(unique(d[intgroup[2]]))[1]]
    ggplot(data = d, aes_string(x = "PCA", y = "PCB", 
                          color = intgroup[1],shape = intgroup[2], label = "names")) +
      geom_hline(yintercept = 0, colour = "gray65") + 
      geom_vline(xintercept = 0, colour = "gray65") +
      geom_text(colour = "black", alpha = 0.8, size = 5, 
                hjust=0, vjust=0, position = "jitter") + 
      geom_point(size=10) + 
      scale_shape_manual(values=required_shapes) +
      xlab(paste0("PC",PCa,": ", round(percentVar[PCa] * 100), "% variance")) +
      ylab(paste0("PC",PCb,": ", round(percentVar[PCb] * 100), "% variance")) +
      labs(title=main)
      #+scale_color_brewer(type="seq", direction=1, palette=1)

    
  } else if (length(intgroup) == 3) {
    #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
    required_shapes <- custom_shape_palette[1:dim(unique(d[intgroup[2]]))[1]]
    ggplot(data = d, aes_string(x = "PCA", y = "PCB", 
                          color = intgroup[1], shape = intgroup[2], 
                          alpha = intgroup[3], label = "names")) +
      geom_hline(yintercept = 0, colour = "gray65") + 
      geom_vline(xintercept = 0, colour = "gray65") +
      geom_text(colour="black", alpha=0.8, size=5, hjust=0, 
                vjust=0, position = "jitter") + 
      geom_point(size=10) + 
      scale_shape_manual(values=required_shapes) +
      xlab(paste0("PC",PCa,": ", round(percentVar[PCa] * 100), "% variance")) + 
      ylab(paste0("PC",PCb,": ", round(percentVar[PCb] * 100), "% variance")) + 
      scale_alpha_discrete(range = c(0.35,0.9)) + 
      labs(title=main)
    
  }
}

############### 2D PCA Custom Genes Set Batch Removed ###########################
plot_PCA_prcomp_batch_removed_custom <- function(x,y,intgroup, PCa=1, PCb=2, returnData=F, geneset, selectset, predefined){
  
  if (predefined==F) {
    if (length(grep(",", geneset)>0)) {
      geneset <- unlist(strsplit(geneset, ",")) }
    if (length(grep(" ", geneset)>0)) {
      geneset <- unlist(strsplit(geneset, " ")) }
  }
  
  else if (predefined==T) {
    geneset <- get_genesets(selectset)
  }
  
  select <- which(rownames(x) %in% geneset)
  pca <- prcomp(t(x[select, ]))
  
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(y)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- GenomicRanges::as.data.frame(colData(y)[, intgroup, drop=FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
  d <- data.frame(PCA = pca$x[, PCa], PCB = pca$x[, PCb], 
                  group = group, intgroup.df, 
                  names = colnames(x))
  
  main = paste0("Batch Effect Removed: PC", PCa, " vs. PC", 
                PCb," using ", length(geneset), " Genes from Custom Gene Set")
  
  if (returnData) {
    attr(d, "percentVaxr") <- percentVar[PCa:PCb]
    return(d)
  }
  
  if (length(intgroup)==1) {
    ggplot(data=d, aes_string(x="PCA", y="PCB", color=intgroup[1], label="names")) +
      geom_hline(yintercept = 0, colour = "gray65") + 
      geom_vline(xintercept = 0, colour = "gray65") + 
      geom_text(colour="black", alpha=0.8, size=5, hjust=0, vjust=0, position="jitter") + 
      geom_point(size=10) + 
      xlab(paste0("PC",PCa,": ", round(percentVar[PCa] * 100), "% variance")) +
      ylab(paste0("PC",PCb,": ", round(percentVar[PCb] * 100),"% variance")) + 
      labs(title=main)
    
    
  } else if (length(intgroup) == 2) {
    #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
    required_shapes <- custom_shape_palette[1:dim(unique(d[intgroup[2]]))[1]]
    ggplot(data=d, aes_string(x="PCA", y="PCB", color=intgroup[1],shape=intgroup[2], label = "names")) +
      geom_hline(yintercept = 0, colour = "gray65") + 
      geom_vline(xintercept = 0, colour = "gray65") +
      geom_text(colour = "black", alpha = 0.8, size = 5, hjust=0, vjust=0, position = "jitter") + 
      geom_point(size=10) + 
      scale_shape_manual(values=required_shapes) +
      xlab(paste0("PC",PCa,": ", round(percentVar[PCa] * 100), "% variance")) +
      ylab(paste0("PC",PCb,": ", round(percentVar[PCb] * 100), "% variance")) +
      labs(title=main)
    
    
  } else if (length(intgroup) == 3) {
    #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
    required_shapes <- custom_shape_palette[1:dim(unique(d[intgroup[2]]))[1]]
    ggplot(data = d, aes_string(x = "PCA", y = "PCB", color = intgroup[1], shape = intgroup[2], alpha = intgroup[3], label = "names")) +
      geom_hline(yintercept = 0, colour = "gray65") + 
      geom_vline(xintercept = 0, colour = "gray65") +
      geom_text(colour = "black", alpha = 0.8, size = 5, hjust=0, vjust=0, position = "jitter") + 
      geom_point(size=10) + 
      scale_shape_manual(values=required_shapes) +
      xlab(paste0("PC",PCa,": ", round(percentVar[PCa] * 100),"% variance")) + 
      ylab(paste0("PC",PCb,": ", round(percentVar[PCb] * 100),"% variance")) + 
      scale_alpha_discrete(range = c(0.35,0.9)) + 
      labs(title=main)
    
  }
}


############### Heatmap Custom Gene Set ###########################
plot_ui_geneset <- function(x, intgroup, geneset, selectset, predefined){
  if (predefined==F) {
    if (length(grep(",", geneset)>0)) {
    geneset <- unlist(strsplit(geneset, ",")) }
    if (length(grep(" ", geneset)>0)) {
      geneset <- unlist(strsplit(geneset, " ")) }
  }
  
  else if (predefined==T) { geneset <- get_genesets(selectset) }
  
  mat <- assay(x)[ which(rownames(assay(x)) %in% geneset), ]
  mat <- mat - rowMeans(mat)
  ig1 = intgroup[1]
  ig2 = intgroup[2]
  
  n <- length(unique(colData(x)[,ig1]))
  cols <- list(hcl(h=seq(15, 375, length=n+1), l=65, c=100)[1:n])
  names(cols)[1] <- ig1
  names(cols[[1]]) <- levels(as.factor(colData(x)[,ig1]))
  
  #main=paste0(nrow(mat), " of ", length(geneset), " Custom Genes found")
  main = "Custom gene heatmap"
  
  if (length(intgroup) == 1) {
    df <- as.data.frame(colData(x)[,ig1])
    ann2_colors = list(
      #Phenotype = c('EC-ExoCtrl' = "gray60", 'EC-ExoLVH' = "gray20")
      #Phenotype = c('EC-ExoCtrl' = "midnightblue", 'EC-ExoLVH' = "red2")
      #Group = c('High Fractality (HF)'="#F8766D", 'No Topography (NT)' = "#00BA38", 'Round Topography (RT)' = '#619CFF')
      Group = c("Frozen" = "#BB9D00", "High Fractality (HF)" = "#F8766D", "No Topography (NT)"= "#00BA38", "Round Topography (RT)" = "#619CFF",
                "Induced Human Podocytes" = "#00C19F","Nephron Progenitor" = "#00B9E3", "Sorted Human Adult Podocytes" = "#DB72FB") 
                
    )
    colnames(df) <- ig1
    rownames(df) <- rownames(colData(x))
    pheatmap(mat, annotation_col = df, annotation_colors = cols, main=main, clustering_method = "ward.D2")
    #pheatmap(mat, annotation_col = df, annotation_colors = cols, fontsize = 15, main=main, color = colorRampPalette(c("blue", "black", "yellow"))(25))
    #pheatmap(mat, annotation_col = df, annotation_colors = ann2_colors, fontsize = 15, main=main,
    #          color = colorRampPalette(c("yellow", "black", "blue"))(25))
    #pheatmap(mat, legend = TRUE, annotation_legend = FALSE, show_rownames = T, 
    #        show_colnames = F, annotation_col = df, annotation_colors = ann2_colors, fontsize_col = 12, fontsize_row = 6, angle_col = 45)
    
  }
  else {
    df <- as.data.frame(colData(x)[,c(ig1, ig2)])
    #pheatmap(mat, annotation_col = df, annotation_colors = cols, main=main)
    pheatmap(mat, annotation_col = df, annotation_colors = cols, fontsize = 15, main=main, color = colorRampPalette(c("blue", "black", "yellow"))(25))
    #pheatmap(mat, annotation_col = df, main=main)
  }
}

############### Heatmap Custom Gene Set Batch Removed ###########################
plot_ui_geneset_batchremoved <- function(x, y, intgroup, geneset, selectset, predefined){

  if (predefined==F) {
    if (length(grep(",", geneset)>0)) {
      geneset <- unlist(strsplit(geneset, ",")) }
    if (length(grep(" ", geneset)>0)) {
      geneset <- unlist(strsplit(geneset, " ")) }
  }
  
  else if (predefined==T) {
    geneset <- get_genesets(selectset)
  }
    
  mat <- x[ which(rownames(x) %in% geneset), ]
  mat <- mat - rowMeans(mat)
  ig1 = intgroup[1]
  ig2 = intgroup[2]
  
  n <- length(unique(colData(y)[,ig1]))
  cols <- list(hcl(h=seq(15, 375, length=n+1), l=65, c=100)[1:n])
  names(cols)[1] <- ig1
  names(cols[[1]]) <- levels(as.factor(colData(y)[,ig1]))
  
  main=paste0(nrow(mat), " of ", length(geneset), " Custom Genes found")
  
  if (length(intgroup) == 1) {
    df <- as.data.frame(colData(y)[,ig1])
    colnames(df) <- ig1
    rownames(df) <- rownames(colData(y))
    pheatmap(mat, annotation_col = df, annotation_colors=cols, main=main)
    #pheatmap(mat, annotation_col = df, main=main)
  }
  else {
    df <- as.data.frame(colData(y)[,c(ig1, ig2)])
    pheatmap(mat, annotation_col = df, annotation_colors=cols, main=main)
    #pheatmap(mat, annotation_col = df, main=main)
  }
}


############### Dendrogram Custom Gene Set ###########################
plot_hclust_custom <- function(x, dist_method, hclust_method, intgroup, geneset,  selectset, predefined){
  
  if (predefined==F) {
    if (length(grep(",", geneset)>0)) {
      geneset <- unlist(strsplit(geneset, ",")) }
    if (length(grep(" ", geneset)>0)) {
      geneset <- unlist(strsplit(geneset, " ")) }
  }
  
  else if (predefined==T) {
    geneset <- get_genesets(selectset)
  }
  
  select <- which(rownames(x) %in% geneset)
  dissimilarity <- 1 - abs(cor(assay(x[select,]), method="spearman"))
  correlation <- cor(assay(x[select,]), method="spearman")
  distance = as.dist(dissimilarity)
  h <- hclust(dist(t(assay(x[select,])),method=dist_method), hclust_method)
  
  ig1 <- intgroup[1]
  ig2 <- intgroup[2]
  
  if (length(intgroup) == 1) {
    dendr <- dendro_data(h, type="rectangle") # convert for ggplot
    idxs <- match(as.character(dendr[["labels"]]$label), rownames(as.data.frame(colData(x))))
    dendr[["labels"]]$Data.Set <- as.data.frame(colData(x))[[ig1]][idxs]
    ggplot() + 
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
      geom_text(data=label(dendr), aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
      geom_point(data=label(dendr), aes(x, y, color=Data.Set), size=4) +
      ylim(-max(segment(dendr)$y)*(1/5), max(segment(dendr)$y))+
      #coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + ylab("")+
      #ggtitle("All Genes (n=18,851)")+
      theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.line.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
  }
  
  else if (length(intgroup) == 2) {
    #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
    required_shapes <- custom_shape_palette[1:length(unique(colData(x)[,ig2]))]
    dendr <- dendro_data(h, type="rectangle") # convert for ggplot
    idxs <- match(as.character(dendr[["labels"]]$label), rownames(as.data.frame(colData(x))))
    dendr[["labels"]]$Data.Set <- as.data.frame(colData(x))[[ig1]][idxs]
    dendr[["labels"]]$Data.Set2 <- as.data.frame(colData(x))[[ig2]][idxs]
    ggplot() +
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
      geom_text(data=label(dendr), aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
      geom_point(data=label(dendr), aes(x, y, color=Data.Set, shape=Data.Set2), size=4) +
      ylim(-max(segment(dendr)$y)*(1/5), max(segment(dendr)$y))+
      scale_shape_manual(values=required_shapes) +
      #coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + ylab("")+
      #ggtitle("All Genes (n=18,851)")+
      theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.line.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
  }
}

############### Dendrogram Custom Gene Set Bath Removed ###########################
plot_hclust_custom_batch_removed <- function(x, y, dist_method, hclust_method,intgroup, geneset, selectset, predefined){
  
  if (predefined==F) {
    if (length(grep(",", geneset)>0)) {
      geneset <- unlist(strsplit(geneset, ",")) }
    if (length(grep(" ", geneset)>0)) {
      geneset <- unlist(strsplit(geneset, " ")) }
  }
  
  else if (predefined==T) {
    geneset <- get_genesets(selectset)
  }
  
  select <- which(rownames(x) %in% geneset)
  
  dissimilarity <- 1 - abs(cor(x[select,], method="spearman"))
  correlation <- cor(x[select,], method="spearman")
  distance = as.dist(dissimilarity)
  h <- hclust(dist(t(x[select,]),method=dist_method), hclust_method)
  
  ig1 <- intgroup[1]
  ig2 <- intgroup[2]
  
  if (length(intgroup) == 1) {
    dendr <- dendro_data(h, type="rectangle") # convert for ggplot
    idxs <- match(as.character(dendr[["labels"]]$label), rownames(as.data.frame(colData(y))))
    dendr[["labels"]]$Data.Set <- as.data.frame(colData(y))[[ig1]][idxs]
    ggplot() + 
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
      geom_text(data=label(dendr), aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
      geom_point(data=label(dendr), aes(x, y, color=Data.Set), size=4) +
      ylim(-max(segment(dendr)$y)*(1/5), max(segment(dendr)$y))+
      #coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + ylab("")+
      #ggtitle("All Genes (n=18,851)")+
      theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.line.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
  }
  
  else if (length(intgroup) == 2) {
    #custom_shape_palette <- c(16,17,15,18,1,2,0,5,7,9,10,11,12,13,14)
    required_shapes <- custom_shape_palette[1:length(unique(colData(y)[,ig2]))]
    dendr <- dendro_data(h, type="rectangle") # convert for ggplot
    idxs <- match(as.character(dendr[["labels"]]$label), rownames(as.data.frame(colData(y))))
    dendr[["labels"]]$Data.Set <- as.data.frame(colData(y))[[ig1]][idxs]
    dendr[["labels"]]$Data.Set2 <- as.data.frame(colData(y))[[ig2]][idxs]
    ggplot() +
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
      geom_text(data=label(dendr), aes(x, y, label=label), angle=90, hjust=1.1, size=3) +
      geom_point(data=label(dendr), aes(x, y, color=Data.Set, shape=Data.Set2), size=4) +
      ylim(-max(segment(dendr)$y)*(1/5), max(segment(dendr)$y))+
      scale_shape_manual(values=required_shapes) +
      #coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + ylab("")+
      #ggtitle("All Genes (n=18,851)")+
      theme(axis.line.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.x=element_blank(), axis.title.x=element_blank(),
            axis.line.y=element_blank(), axis.ticks.y=element_blank(),
            axis.text.y=element_blank(), axis.title.y=element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank())
  }
}








##############################
## DE Analysis Tab ###########
##############################

########### Normalizartipn Boxplots ########## 
plot_DESEQ_transformation_boxplot <- function(x){
  if (is.null(x)) { return() }
  par(mfrow=c(3,1))
  boxplot(counts(x), main = "Raw Read Counts", xaxt='n')
  boxplot(log2(counts(x) + 1), main = "log2 Normalized Read Counts", xaxt='n')
  boxplot(counts(x,normalized =TRUE),main = "DESeq Normalized Read Counts", xaxt='n')
  
}


########### Make Results Table #############
final_contrast <- function(x, factorU, levelsU){
  if (is.null(factorU) || is.null(levelsU)) { return() }
  else if ((length(levelsU)<2)) { return() }
  else{ res <- results(x, contrast=c(factorU, levelsU[1], levelsU[2])); return(res) }
}	



########### Generate Design Formula #############
des_form_text <- function(factorU, controlU) {
  if (length(factorU) == 0) { return() }
  else if (length(controlU) == 0) { out<-paste("~", factorU) }
  else if (length(controlU) == 1) { out<-paste("~", controlU[1], "+", factorU) }
  else if (length(controlU) == 2) { out<-paste("~", controlU[1], "+", controlU[2], "+", factorU) }
  else if (length(controlU) == 3) { out<-paste("~", controlU[1], "+", controlU[2], "+", controlU[3], "+", factorU) }
  return(out)
}

########### Text for Design/DE Levels #############
lev_text <- function(factorU, levelU) {
  if (length(factorU) == 0) { return() }
  else {
    if (length(levelU) == 0) { return() }
    else if (length(levelU) == 1) { return(levelU[1]) }
    if (length(levelU) == 2) { return(paste(levelU[1], "vs.", levelU[2]))  }
  }
}

  
########### Preform DE Step #############
de_1 <- function(x, factorU, controlU) {
  foi <- factorU
  if (length(factorU) == 0) { return() }
  
  else {
    if (length(controlU) == 0){ des <- paste("~", foi)  }
    else if (length(controlU) == 1) { des <- paste("~",controlU[1],"+",foi) }
    else if (length(controlU) == 2) {  des <- paste("~",controlU[1],"+",controlU[2],"+",foi) }
    else if (length(controlU) == 3) { des <- paste("~",controlU[1],"+",controlU[2],"+", controlU[3],"+",foi) }
    else if (length(controlU) == 4) { des<-paste("~",controlU[1],"+",controlU[2],"+", controlU[3],"+",controlU[4],"+",foi) }
  }

  design(x) <- formula(des)
  x <- DESeq(x)
  return(x)
}


########### Store Results Table #############
de_res <- function(x, factorU, levelU, padj_in) {
  if (length(levelU) < 2) { return() } 
  else { l1<-levelU[1]; l2<-levelU[2] }
  
  rt <- results(x, contrast=c(factorU, levelU[1], levelU[2]), pAdjustMethod=padj_in)
  return(rt)
}


########### Generate/Store Results Summary #############
return_res_summary <- function(x) {
  if (is.null(x)) { return() }
  else {
  l2fc_vec <- x$log2FoldChange
  tot <- dim(x)[1]
  res_summary <- data.frame(total=tot,
                  outlier=paste0(sum(is.na(x$padj)), " (", round((sum(is.na(x$padj))/tot)*100,2), "%)"),
                  l2fcgt1=paste0(length(which(l2fc_vec > 0)), " (", round((length(which(l2fc_vec > 0))/tot)*100,2), "%)"),
                  l2fclt1=paste0(length(which(l2fc_vec < 0)), " (", round((length(which(l2fc_vec < 0))/tot)*100,2),"%)"),
                  l2fcgt2=paste0(length(which(l2fc_vec > 2)), " (",  round((length(which(l2fc_vec > 2))/tot)*100,2), "%)"),
                  l2fclt2=paste0(length(which(l2fc_vec < -2)), " (", round((length(which(l2fc_vec < -2))/tot)*100,2),"%)"),
                  padjlt05=length(which(x$padj < .05)))
  colnames(res_summary) <- c("Total Genes", "Outliers", "Log2FC > 0 (up)", "Log2FC < 0 (down)", "Log2FC > 2 (up)", "Log2FC < -2 (down)", "Total padj < 0.05")
  return(res_summary)
  }
}


########### Top DE Genes Tables (Up/Down) #############
return_H_topDE <- function(x, n) {
  if (is.null(x)) { return()}
  else{
    if (n==0) { return() }
    else {
      h <- as.data.frame(x[head(order(x$log2FoldChange, decreasing=T), n),])[,c(1,2,3,6)]
      return(h)
    }
  }
}

return_L_topDE <- function(x, n) {
  if (is.null(x)) {  return() }
  else {
    if (n==0) { return() }
    else {
      l <- as.data.frame(x[head(order(x$log2FoldChange,  decreasing=F), n),])[,c(1,2,3,6)]
      return(l)
    }
  }
}



########### DESEq Gene Expression Plot #############
plot_gene <- function(x, factorU, annot_var, gene) {
  if (length(factorU) == 0) { return() }
  
  else {
    if (is.null(annot_var)) {
      plotCounts(x, gene=gene, intgroup=factorU, pch=19,  xlab="Factor of Interest")
    }
    else {
      c1 <- as.numeric(x@colData[[annot_var]])
      c1 <- factor(c1, labels=1:length(unique(c1)))
      plotCounts(x, gene=gene, intgroup=factorU, pch=19,  col=c1, xlab="Factor of Interest")
      legend("topright", legend=unique(x@colData[[annot_var]]), fill=unique(c1),cex=1, bty="n")
    }
  }
}

