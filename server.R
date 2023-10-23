library(shiny)
library(DESeq2)
library(limma)
library(ggplot2)
library(gplots)
#library(vsn)
library(genefilter)
library(Hmisc)
library(corrplot)
library(pvclust)
library(pheatmap)
library(calibrate)
library(lumi)
library(plotly)
library(RColorBrewer)
library(DT)
library(lme4)
library(purrr)
#library(openxlsx)
library(ggdendro)



## Max File Upload Size = 100MB
options(shiny.maxRequestSize=250*2048^2)
if (length(dev.list()>0)){graphics.off()}
load("data/NMETH-BC11016B-pluritestsub.unk", envir= .GlobalEnv)
source("Basic_DE_analysis.R")
illumina_ids <<- read.delim("data/Illumina_PluritestID.txt", header = TRUE)
#phenotypes <<- read.xlsx("Z:/Projects/Stem Cell Projects/U01/PhenotypeData/Copy of hiPSC_ID_master_list_041315_fromCharles_PA_110316.xlsx", sheet=1)



# Define server logic 
shinyServer(
  function(input, output, session) {

    
    #### Read-In Design and Counts Files
    storeDesign <- reactive({ 
      if(is.null(input$designmatrix)) return(NULL)
      design_df <- read.csv(input$designmatrix$datapath, 
                              header=T, sep="\t", row.names = 1)
      assign("raw_metadata", design_df, envir=.GlobalEnv)
      return(design_df) })
    
    storeRawcounts <- reactive({   
      if(is.null(input$rawcounts)) return(NULL)
      rawcounts_df <- read.table(input$rawcounts$datapath, 
                                 header=T, sep="\t", row.names=1)
      assign("raw_countdata", rawcounts_df, envir=.GlobalEnv)
      return(na.omit(rawcounts_df)) })
    
    
    #### Process Samples
    rtransformation <- reactive({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      read_transformation(dds(), type = input$Transform)
    })
    
    dds <- reactive({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      run_DESeq2(filtered_gene_samples())
    })

    filtered_gene_samples <- reactive({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      filter_gene_sample_list(read_data(), minMedian = input$MinMedianCount, samples = samples_TO_analyze())
    })
    
    read_data <- reactive({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      consolidate_genes_samples(counts_table = storeRawcounts(), meta_data = storeDesign())
    })
    
    samples_TO_analyze <- reactive({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      each_var <- map(vars(), ~ samples_to_analyze(read_data()[["meta_data"]], .x, input[[.x]]))
      purrr::reduce(each_var,intersect)
    }) 
    
    
    ## Dynamic UI 
    vars <- reactive(names(read_data()[["meta_data"]]))
    output$filter <- renderUI(
      map(vars(), ~ make_ui(read_data()[["meta_data"]][[.x]], .x))
    )

    #### Pluritest
    pluritest_data <- reactive({
      run_pluritest(rtransformation(), illumina_ids, returnData = TRUE)
    })
    
    #### BE Removal
    batch_removal <- reactive({
         remove_batcheffect(rtransformation(),raw_metadata[samples_TO_analyze(),], input$BEFactor)
     })
    
    #### Run DESEQ DE Analysis
    generate_design <- reactive({
      de_1(dds(), factorU=input$Factor, 
           controlU=input$ControlInteract) 
    })
    
    #### Store DE Results Table
    generate_results <- reactive({
      de_res(generate_design(), factorU=input$Factor, 
             levelU=input$FactorLevels, padj_in=input$padj)
    })
    
    res <- reactive({
      input$DEButton
      if (is.null(input$Factor)) {
        return() }
      else {
        de_1(dds(), factorU=input$Factor, controlU=input$Controls) }
    })
    
    levelsU <- reactive({
      levels(raw_metadata[,input$Factor])
    })
    
    degenes <- reactive({
      toupper(unlist(strsplit(input$GeneName, split=",")))
    })
    
    pcs <- reactive({
      unlist(strsplit(input$PCs, split=","))
    })
    
    custom_pcs <- reactive({
      unlist(strsplit(input$CustomPCs, split=","))
    })
    
    intgrouplevels <- reactive({
      input$HeatmapAnnotation
    })
    
    intgrouplevels_custom_geneset <- reactive({
      input$HeatmapAnnotation_custom_geneset
    })
    
    intgrouplevels_pvclust <- reactive({
      input$HeatmapAnnotation_pvclust
    })
    
    pcrotation <- reactive({
      if (input$BatchEffect == FALSE) {
        pca_rotation(rtransformation(), PC = unlist(input$PCrotation)[1])
      }
      else {
        pca_rotation_batch_removed(batch_removal(),PC = unlist(input$PCrotation)[1])
      }
    })
    
    CustomPcrotation <- reactive({
      if (input$BatchEffect == FALSE) {
        custom_pca_rotation(rtransformation(), PC = unlist(input$CustomPCrotation)[1], 
                            geneset=input$GenesUI)
      }
      else {
        custom_pca_rotation_batch_removed(batch_removal(),PC = unlist(input$CustomPCrotation)[1],
                                          geneset=input$GenesUI)
      }
    })
    
    pcloadings <- reactive({
      if (input$BatchEffect == FALSE) {
        pca_loadings(rtransformation(), ntop=input$PCAGenes)
      }
      else {
        pca_loadings_batch_removed(batch_removal(), ntop=input$PCAGenes)
      }
    })
    
    CustomPcloadings <- reactive({
      if (input$BatchEffect == FALSE) {
        custom_pca_loadings(rtransformation(),
                            geneset=input$GenesUI)
      }
      else {
        custom_pca_loadings_batch_removed(batch_removal(),
                                          geneset=input$GenesUI)
      }
    })
    
    
    
    #### UI Server Outputs ####
    output$FactorLevelsUI <- renderUI({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
       selectInput("FactorLevels", "Select the levels to use for Contrast", 
                   choices = levelsU(),
                   multiple = TRUE,
                   selected = c(levelsU()[1], levelsU()[2]))
      checkboxGroupInput("FactorLevels", 
                         "Select the levels to use for Contrast", 
                         choices = levelsU(),
                         selected = c(levelsU()[1], levelsU()[2]))
    })
    
    output$PCAGenesUI <- renderUI({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      sliderInput("PCAGenes", label = "Number of highly variable genes for PCA and Dendrogram", min = 100, 
                  max = length(rownames(counts(dds()))), value = 500, step = 100)
    })

    output$BEUI <- renderUI({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      selectInput("BEFactor","", 
                  choices = colnames(raw_metadata), multiple = T,
                  selected = colnames(raw_metadata)[1])
      #checkboxGroupInput("BEFactor", "", choices = colnames(raw_metadata))
    })
    
    output$HeatmapAnnotationUI <- renderUI({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      selectInput("HeatmapAnnotation", "Choose Factors for PCA, Dendrogram, and Heatmap Annotation", 
      choices = colnames(input_object$meta_data),
      selected = colnames(input_object$meta_data)[1],
      multiple=TRUE)
    })
    
    output$HeatmapAnnotationUI_custom_geneset <- renderUI({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      selectInput("HeatmapAnnotation_custom_geneset", "Choose Factors for PCA, Dendrogram, and Heatmap Annotation", 
                  choices = colnames(input_object$meta_data),
                  selected = colnames(input_object$meta_data)[1],
                  multiple=TRUE)
    })
    
    output$HeatmapAnnotationUI_pvclust <- renderUI({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      selectInput("HeatmapAnnotation_pvclust", "Choose Factors for Dendrogram Annotation", 
                  choices = colnames(input_object$meta_data),
                  selected = colnames(input_object$meta_data)[1],
                  multiple=TRUE)
    })
    
    output$PCAGenesUI_pvclust <- renderUI({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      sliderInput("PCAGenes_pvclust", label = "Number of highly variable genes for pvclust", min = 500, 
                  max = length(rownames(counts(dds()))), value = 500, step = 500)
    })
    
    
    
    
    output$FactorUI2 <- renderUI({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      checkboxGroupInput("Factor",
                   label = "Choose the factor of interest for DE analysis",
                   choices = colnames(input_object$meta_data))
    })
    
    output$AnnotVarUI <- renderUI({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      checkboxGroupInput("AnnotVar",
                          "Select One Factor for Expression Plot Annotation",
                           choices = colnames(input_object$meta_data),
                           selected = NULL)
    })

    #### OUTPUT QC PLOTS ####
    output$HKGplot <- renderPlot({
      if (input$PlotsButton == F) {
        return() }
      else {
      isolate({plot_housekeeping(rtransformation())}) }
    })

    output$RTBoxPlot <- renderPlot({
      if (input$PlotsButton == F) {
        return() }
      else {
      isolate({plot_read_transformation_boxplot(dds(),
                              rtransformation(),
                              type=input$Transform)}) }
    })


    output$QCGplot <- renderPlot({
      if (input$PlotsButton == F) {
        return() }
      else {
        if (input$BatchEffect == FALSE) {
        isolate({plot_qcgenes(rtransformation(), 
                              intgroup=input$HeatmapAnnotation)})
        }
        else{
        isolate({plot_qcgenes_batchremoved(batch_removal(), 
                                           rtransformation(), 
                                           intgroup=input$HeatmapAnnotation)}) }
      }
    })


    output$PCAplot <- renderPlot({
      if (input$PlotsButton == F) {
        return() }
      else {
        # if (input$BatchEffect == FALSE) {
        # isolate({ plot_PCA_prcomp(rtransformation(), 
        #                 intgroup=intgrouplevels(), 
        #                 PCa=as.numeric(pcs()[1]), PCb=as.numeric(pcs()[2]), 
        #                 ntop=input$PCAGenes, 
        #                 main = paste(input$Transform," based PC", 
        #                              pcs()[1], " vs. PC", pcs()[2], " using ",
        #                              input$PCAGenes, " genes", sep="")) }) }
        if (input$BatchEffect == FALSE) {
          isolate({ plot_PCA_prcomp(rtransformation(), 
                                    intgroup=intgrouplevels(), 
                                    PCa=as.numeric(pcs()[1]), PCb=as.numeric(pcs()[2]), 
                                    ntop=input$PCAGenes) }) }
  
        else{
          isolate({ plot_PCA_prcomp_batch_removed(batch_removal(),rtransformation(),
                          #intgroup=input$HeatmapAnnotation, 
                          intgroup=intgrouplevels(),
                          PCa=as.numeric(pcs()[1]), PCb=as.numeric(pcs()[2]), 
                          ntop=input$PCAGenes, 
                          # main = paste("Batch effect removed ",
                          #              input$Transform," based PC", 
                          #              pcs()[1], " vs. PC", pcs()[2]," using ",
                          #              input$PCAGenes, " genes", sep="")
                          ) }) }
      }
    })
    
    ## PC ANOVA ##
    output$pc_anova <- DT::renderDataTable({
      if (input$PlotsButton == F) { return() }
      else {
        if (input$BatchEffect == FALSE) {
          isolate({ pc_anova_table(rtransformation(), ntop=input$PCAGenes) }) } 
        else if (input$BatchEffect == TRUE) {
          isolate({ pc_anova_table_batch_removed(batch_removal(), rtransformation(), ntop=input$PCAGenes) }) } 
        }
      })


    
    
    output$PCA3Dplot <- renderPlotly({
      if (input$PlotsButton == F) {
        return() }
      else {
        if (input$BatchEffect == FALSE) {
        isolate({ plot3DPCA(rtransformation(), intgroup=input$HeatmapAnnotation, 
                            PCa=as.numeric(pcs()[1]), PCb=as.numeric(pcs()[2]),
                            PCc=as.numeric(pcs()[3]), ntop=input$PCAGenes, 
                            main = paste(input$Transform," based PC", pcs()[1], 
                                " vs. PC", pcs()[2], " vs. PC", pcs()[3],
                                " using ", input$PCAGenes, " genes", sep="")) }) }
        
        else {
        isolate({ plot3DPCA_batch_removed(batch_removal(), rtransformation(), 
                                          intgroup=input$HeatmapAnnotation, 
                            PCa=as.numeric(pcs()[1]), PCb=as.numeric(pcs()[2]),
                            PCc=as.numeric(pcs()[3]), ntop=input$PCAGenes, 
                            main = paste(input$Transform," based PC", pcs()[1], 
                                " vs. PC", pcs()[2]," vs. PC", pcs()[3],
                                " using ", input$PCAGenes, " genes", sep="")) }) }
        }
      
    })


    output$PTplot <- renderPlot({
      if (input$PlotsButton == F) {
        return() }
      else {
      isolate({run_pluritest(rtransformation(), 
                             intgroup=input$HeatmapAnnotation, 
                             illumina_ids)}) }
    })


    output$genecluster100 <- renderPlot({
      if (input$PlotsButton == F) {
        return() }
      else {
        if (input$BatchEffect == FALSE) {
          isolate({plot_gene_cluster100(rtransformation(), 
                                        intgroup=input$HeatmapAnnotation)})
        }
        else {
          isolate({plot_gene_cluster100_batch_removed(batch_removal(), 
                                        rtransformation(), 
                                        intgroup=input$HeatmapAnnotation)}) }
      } 
    })

    output$genecluster500 <- renderPlot({
      if (input$PlotsButton == F) {
        return() }
      else {
        if (input$BatchEffect == FALSE) {
          isolate({plot_gene_cluster500(rtransformation(), 
                                        intgroup=input$HeatmapAnnotation)})
        }
        else {
          isolate({plot_gene_cluster500_batch_removed(batch_removal(), 
                                        rtransformation(), 
                                        intgroup=input$HeatmapAnnotation)}) }
      }
    })
    

    
    output$hclust <- renderPlot({
      if (input$PlotsButton == F) {
        return() }
      else {
         if (input$BatchEffect == FALSE) {
          isolate({plot_hclust(rtransformation(), hclust_method = 'ward.D2',
                        intgroup=input$HeatmapAnnotation, 
                        ntop=input$PCAGenes)})
         }
         else {
           isolate({plot_hclust_batch_removed(batch_removal(), rtransformation(),
                        hclust_method = 'ward.D2', 
                        intgroup=input$HeatmapAnnotation,
                        ntop=input$PCAGenes)}) }
       }
    })  
    
    
    output$pvclust <- renderPlot({
      if (input$pvClustButton == F) {
        return() }
      else {
        if (input$BatchEffect == FALSE) {
          isolate({plot_pvclust(rtransformation(), hclust_method = 'ward.D2',
                               intgroup=input$HeatmapAnnotation_pvclust, geneset=input$GenesUI_pvclust,
                               ntop=input$PCAGenes,  nboots=as.numeric(input$nboots))})
        }
        else {
          isolate({plot_pvclust_batch_removed(batch_removal(), rtransformation(),
                                             hclust_method = 'ward.D2', 
                                             intgroup=input$HeatmapAnnotation_pvclust, geneset=input$GenesUI_pvclust, 
                                             ntop=input$PCAGenes, nboots=as.numeric(input$nboots))}) }
      }
    })  
    
    
    
#### CUstom Gene Set Analysis ####
    
    
    ## Heatmap with CUstom Gene Set
    output$ui_geneset <- renderPlot({
      if (input$CustomGeneButton == F) {
        return() }
    
      else {
        if (input$BatchEffect == FALSE) {
          isolate({plot_ui_geneset(rtransformation(), 
            intgroup=intgrouplevels_custom_geneset(), geneset=input$GenesUI,
            selectset=input$SelectSetUI, predefined=input$PreDefined)})
        }
        else {
          isolate({plot_ui_geneset_batchremoved(batch_removal(), 
           rtransformation(), intgroup=intgrouplevels_custom_geneset(), 
           geneset=input$GenesUI, selectset=input$SelectSetUI, predefined=input$PreDefined)})}
      }
    })
    

    #### Dynamic custom geneset heatmap size #### color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(10)
    plotCount <- reactive({
      if (input$PreDefined==F) {
        if (length(grep(",", input$GenesUI)>0)) {
          geneset <- unique(unlist(strsplit(input$GenesUI, ","))) }
        if (length(grep(" ", input$GenesUI)>0)) {
          geneset <- unique(unlist(strsplit(input$GenesUI, " ")))  }
        px=(400+(sum(geneset%in%GeneNames_Global)^1.5))
        #observeEvent(input$GenesUI, { print(paste0("TEST OBSERVE EVENT: ", px)) })
        return(px)
      }
      else if (input$PreDefined==T) {
        geneset <- get_genesets(input$SelectSetUI)
        px=(400+(sum(geneset%in%GeneNames_Global)^1.5))
        return(px)
        }
    })

    output$plot.ui <- renderUI({
      if (input$CustomGeneButton == F) { return() }
      else isolate({ plotOutput("ui_geneset", width = "100%", height=paste0(plotCount(),"px")) })
      #else isolate({ plotOutput("ui_geneset", width = "100%", height = "1000px") })
    })

    

    
    ## Sample PCA using Custom Gene Set
    output$PCAplot_Geneset <- renderPlot({
      if (input$CustomGeneButton == F) {
        return() }
      else {
        if (input$BatchEffect == FALSE) {
          isolate({ plot_PCA_prcomp_custom(rtransformation(), 
                        intgroup=intgrouplevels_custom_geneset(), 
                        PCa=as.numeric(custom_pcs()[1]),PCb=as.numeric(custom_pcs()[2]),
                        geneset=input$GenesUI, selectset=input$SelectSetUI, 
                        predefined=input$PreDefined) }) }
        
        
        else{
          isolate({ plot_PCA_prcomp_batch_removed_custom(batch_removal(), rtransformation(),
                      intgroup=intgrouplevels_custom_geneset(),
                      PCa=as.numeric(custom_pcs()[1]), PCb=as.numeric(custom_pcs()[2]),
                      geneset=input$GenesUI, selectset=input$SelectSetUI, 
                      predefined=input$PreDefined) }) }
      }
    })
    
    ## Hclust Custom Gene Set
    output$hclust_Geneset <- renderPlot({
      if (input$CustomGeneButton == F) {
        return() }
      else {
        if (input$BatchEffect == FALSE) {
          isolate({plot_hclust_custom(rtransformation(),
                               dist_method="euclidean", 
                               hclust_method="ward.D2",
                               intgroup=intgrouplevels_custom_geneset(), 
                               geneset=input$GenesUI, selectset=input$SelectSetUI, 
                               predefined=input$PreDefined)})
        }
        else {
          isolate({plot_hclust_custom_batch_removed(batch_removal(), rtransformation(),
                                             dist_method="euclidean", 
                                             hclust_method="ward.D2",
                                             intgroup=intgrouplevels_custom_geneset(),
                                             geneset=input$GenesUI, selectset=input$SelectSetUI, 
                                             predefined=input$PreDefined)}) }
      }
    })  
    
    
    
    
    ## raw/normalized QC boxplots
    output$GeneCount_QC <- DT::renderDataTable({
      if (is.null(input$count_table_QC)) {
        return()
      }
      if (input$count_table_QC == "Normalized (rlog/vst)") {
        isolate({assay(rtransformation(), normalized=TRUE)})
      }
      else if (input$count_table_QC == "Raw") {
        #isolate({counts(dds(), normalized=FALSE)})
        isolate({raw_counts})
      }
    }, rownames = T)

    
    
    
    
    ##### Output Design Matrix #####
    output$DM <- DT::renderDataTable({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      else ({input_object$meta_data[samples_TO_analyze(),]})
    }, options=list(rownames = T, pageLength=500))    
    
    
    ##### Output Raw Counts #####
    #output$RC <- shiny::renderTable({
    output$RC <- DT::renderDataTable({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      #else ({raw_countdata[samples_TO_analyze()]})
      else ({ input_object$count_data[samples_TO_analyze()]})
    }, options=list(rownames = T, pageLength=500))  
    
    
    # ##### Output Phenotype Data #####
    # output$Pheno <- DT::renderDataTable({
    #   if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
    #   else ({phenotypes[which(phenotypes$id.cdi %in% 
    #                       substr(input_object$meta_data[samples_TO_analyze(),]$CellLine, 2, 5)),]})
    # }, options=list(rownames = T, pageLength=500))  
  

    
    #### Output DE ANalysis Page ####
    
    ### RAW/NORMALIZED COUNTS for DE
    output$GeneCount <- DT::renderDataTable({
      if (is.null(input$count_table)) {
        return()
      }
      if (input$count_table == "Normalized (DESeq)") {
        isolate({counts(dds(), normalized=TRUE)})
      }
      else if (input$count_table == "Raw") {
        #isolate({counts(dds(), normalized=FALSE)})
        isolate({raw_counts})
      }
    }, rownames = T)
    
    
    ## Print Formula and Levels Text
    output$DesFormula <- renderText({
      input$DEButton
      #isolate({des_form_text(factorU=input$Factor, controlU=input$ControlInteract)})
      des_form_text(factorU=input$Factor, controlU=input$ControlInteract)
    })
    
    output$LevelsText <- renderText({
      #input$DEButton
      #isolate({des_form_text(factorU=input$Factor, controlU=input$ControlInteract)})
      lev_text(factorU=input$Factor, levelU=input$FactorLevels)
    })

    controlsU <- reactive({
      cnts <- c()
      for (i in colnames(raw_metadata)) {
        cnts <- c(cnts, i)
      }
      for (i in 1:length(colnames(raw_metadata))) {
        for (j in 1:(length(colnames(raw_metadata))-1)) {
        cnts <- c(cnts, paste0(colnames(raw_metadata)[i], ":", 
                               colnames(raw_metadata)[j]))
        }
      }
      cnts
    })
    
    output$ControlsUI <- renderUI({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      selectInput("ControlInteract","Choose Control Variables",
      ##checkboxGroupInput("ControlInteract","Choose Control Variables",
                         choices = controlsU(), multiple=T)
                         ##, selected = c(controlsU()[1], controlsU()[2]))
    })
    
    levelsU <- reactive({
      levels(as.factor(raw_metadata[,input$Factor]))
    })
    
    output$FactorLevelsUI <- renderUI({
      if(is.null(input$designmatrix) | is.null(input$rawcounts)) return(NULL)
      #print(levelsU)
      selectInput("FactorLevels", "Choose Two Levels for DE analysis",
      ##checkboxGroupInput("FactorLevels", "Choose Two Levels for DE analysis",
                         choices = levelsU(),
                         multiple=T)
                         ##, selected = c(levelsU()[1], levelsU()[2]))
    })
    

    output$DESEQBoxPlot <- renderPlot({
      input$DEButton
      isolate({plot_DESEQ_transformation_boxplot(generate_design())})
    })

    ## DE Results Table
    output$ResTable <- DT::renderDataTable({
      input$DEButton
      isolate({as.data.frame(generate_results())})
    }, rownames = T)

    
    ## DE Results Summary
    output$res_summary <- shiny::renderTable({
      input$DEButton

      isolate({return_res_summary(generate_results())})
    })
    
    ### Top DE Genes
    output$top_H_DEgenes <- shiny::renderTable({
      input$DEButton
      isolate({return_H_topDE(generate_results(), input$DEGenes)})
    }, rownames = T, digits = -2)
    
    output$top_L_DEgenes <- shiny::renderTable({
      input$DEButton
      isolate({return_L_topDE(generate_results(), input$DEGenes)})
    }, rownames = T, digits = -2)

    ### Expression Plot
    output$plot_gene <- renderPlot({
      input$DEButton
      isolate({plot_gene(generate_design(), input$Factor, input$AnnotVar, input$GeneName)})
    })
    
    
    ### DOWNLOADS ###
    
    output$downloadPCrotation <- downloadHandler(
      filename = function() {paste(input$PCrotation, "_rotations.csv",sep="")},
      content = function(file) {
        write.table(as.data.frame(pcrotation()), file, sep="\t", quote = F)
      }
      )
    
    output$CustomDownloadPCrotation <- downloadHandler(
      filename = function() {paste(input$CustomPCrotation, "_rotations.csv",sep="")},
      content = function(file) {
        write.table(as.data.frame(CustomPcrotation()), file, sep="\t", quote = F)
      }
    )
    
    output$downloadPCloadings <- downloadHandler(
      filename = function() {paste(input$PCloadings, "_loadings.csv",sep="")},
      content = function(file) {
        write.table(as.data.frame(pcloadings()), file, sep="\t", quote = F)
      }
    )
    
    output$CustomDownloadPCloadings <- downloadHandler(
      filename = function() {paste(input$CustomPCloadings, "_loadings.csv",sep="")},
      content = function(file) {
        write.table(as.data.frame(CustomPcloadings()), file, sep="\t", quote = F)
      }
    )
    
    
    output$downloadDEResults <- downloadHandler(
      filename = function() {paste(input$Factor,"_",levels(raw_metadata[,input$Factor])[1],"_vs_",
                                  levels(raw_metadata[,input$Factor])[2],"_deseq2_results.csv", sep=" ")},
      content = function(file) {
        write.csv(as.data.frame(res()), file)
      }
      )
    
    
    output$downloadResultsTable <- downloadHandler(
      filename = function() {
        paste("results_table", Sys.Date(), ".tsv", sep="")
          },
      content = function(file) {
        #write.csv(generate_results(), file, quote=F, sep="\t")
        write.table(generate_results(), file, quote=F, sep="\t", col.names=T)
      }
    )

    output$downloadDESeqNormReads <- downloadHandler(
      filename = "deseq2_normalized_read_counts.csv",
      content = function(file) {
        write.table(as.data.frame(counts(generate_design(), normalized=TRUE)), file, sep="\t", row.names = T, quote = F)
      }
    )
    
    output$downloadDESeqReads <- downloadHandler(
      filename = "deseq_read_counts.csv",
      content = function(file) {
        write.table(as.data.frame(counts(generate_design(), normalized=FALSE)), file, sep="\t", row.names = T, quote = F)
      }
    )
    
    output$downloadDM <- downloadHandler(
      filename = function() {paste("DM",".csv", sep=" ")},
      content = function(file) {
        write.table(as.data.frame(raw_metadata[samples_TO_analyze(),]), file, quote=F, sep="\t", col.names = T, row.names = T)
      }
    )
    
    output$downloadRawReads_QC <- downloadHandler(
      filename = function() {paste("Library_size_normalized_read_counts",".csv", sep=" ")},
      content = function(file) {
        write.csv(as.data.frame(raw_counts), file)
      }
    )

    output$downloadNormReads_transformed <- downloadHandler(
      filename = function() {paste(input$Transform,"_normalized_read_counts",".csv", sep=" ")},
      content = function(file) {
        write.table(assay(rtransformation()), file, sep="\t",
                    row.names = T, quote = F, col.names = T)
      }
    )

    output$downloadPluriTestOutput <- downloadHandler(
      filename = function(){paste("PluriTest_pluripotency_novelty_score_output",".csv",sep=" ")},
      content = function(file){
        write.csv(as.data.frame(pluritest_data()), file)
      })
    
  }
)