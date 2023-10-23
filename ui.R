
library(ggplot2)
library(gplots)
library(pvclust)
library(pheatmap)
library(plotly)
# library(openxlsx)
library(ggdendro)



#### Test conditional UI
#### Two stages : QC and/or DE results

shinyUI(fluidPage(
  titlePanel("Differential Gene Expression Analysis"),
  
  ################ SIDEBAR Panel ###################
  sidebarLayout(
    sidebarPanel(
      
      ###### Upload and Filtering ######
      #tags$hr(),
      h3(strong("File Upload and Filtering")),
      tags$hr(),
      
      fileInput("designmatrix", "Upload Design Matrix File", multiple = FALSE),
      fileInput("rawcounts", "Upload Raw Counts File", multiple = FALSE),  
      #br(),
      #br(),
      
      
        #h3(strong("Set Filtering Parameters")),
        #tags$hr(),
        #br(),

        selectInput("Transform", 
            label = "Select Method for Read Transformation",
            choices = c("rlog","vst"), selected = "vst"),
        br(),

        checkboxInput("BatchEffect", strong("Remove Batch Effect?"), FALSE),
        conditionalPanel(
           condition = "input.BatchEffect == true",
           htmlOutput("BEUI")),
        br(),
      
        
        checkboxInput("Filters", strong("Filter Samples?"), FALSE),
        conditionalPanel(
          condition = "input.Filters == true",
            uiOutput("filter")
        ),
        br(),

        sliderInput("MinMedianCount", label = "Set the Minimum Median Read Depth", 
                    min = 0, max = 100, value = 10, step = 5),
        #tags$hr(),
        br(),
        br(),


        ###### QC Conditional Panel ######
        conditionalPanel(condition="input.tabselected==1",
          br(),
          #tags$hr(),
          h3(strong("Set QC Plotting Parameters")),
          tags$hr(),
          #br(),
                
          htmlOutput("PCAGenesUI"),
          br(),
          
          htmlOutput("HeatmapAnnotationUI"),
          br(),
      
          textInput("PCs","Principal Components to Plot", "1,2,3"),
          br(),
          
          textInput("PCrotation","Principal Component Rotation for Download","PC1"),
          br(),
  
          downloadButton("downloadNormReads_transformed", "Download Normalized Counts"),
          br(),

          downloadButton("downloadPCrotation", "Download PC Rotation"),
          br(),
          
          downloadButton("downloadPCloadings", "Download PC Loadings"),
          br(),
  
          #downloadButton("downloadPluriTestOutput", "Download PluriTest Output"),
          #br(),
          br(),
          
          actionButton("PlotsButton", strong("Update Plots")),
          br()
        ),


        ###### DE Conditional Panel ######
        conditionalPanel(condition="input.tabselected==2",
            br(),
            h3(strong("Differential Expression Analysis Parameters")),
            ## Factor of Interest
            htmlOutput("FactorUI2"),

            # Factor Levels
              htmlOutput("FactorLevelsUI"),

            # Control Variables 
              htmlOutput("ControlsUI"),

            # PADJ MEthod
              selectInput("padj",
                            label = "Choose a Method for Multiple Test Correction",
                            choices = c("holm", "hochberg", "hommel", 
                                        "bonferroni", "BY", "fdr", "none"), 
                          selected = "fdr"
              ),
              br(),
              
              htmlOutput("FactorLevelsUIDE"),
              br(),        
            
              sliderInput("DEGenes", 
                          label = "Number of Genes in Top DE Genes Tables", 
                          min = 0, max = 200, value = 5, step = 5 ),
              br(),

              textInput("GeneName", "Gene for Expression Plot", "NPPB"),
              br(),
            
            htmlOutput("AnnotVarUI"),
            br(),
            downloadButton('downloadDESeqReads',"Download DESeq Counts"), 
            br(),
            downloadButton('downloadDESeqNormReads',"Download DESeq Normalized Counts"), 
            br(),
            downloadButton('downloadResultsTable', "Download Results Table"),
            br(),
            br(),
            
            actionButton("DEButton", strong("Perform/Update DE Analysis"))
        ),

            ###### Custom Geneset Conditional Panel ######
            conditionalPanel(condition="input.tabselected==3",
                 br(),
                 h3(strong("Custom Gene Set Parameters")),
                 
                 #h5("Enter Custom Gene Set"),
                 textInput("GenesUI","Enter Custom Gene Set", value=NA),
                 br(),
                 #checkboxInput("UserEntered", strong("Enter Your Own Gene Set?"), FALSE),
                 #conditionalPanel(
                 #    condition = "input.UserEntered == true",
                 #  textInput("GenesUI", "", "")),
                 checkboxInput("PreDefined", strong("Chose from a List of PreDefined Gene Sets?"), FALSE),
                 conditionalPanel(
                   condition = "input.PreDefined == true",
                   selectInput("SelectSetUI",
                               label = "Chose a PreDefined Gene Set",
                               choices = c("Cardiac Hypertrophy Marker Genes", 
                                           "HyperGEN Linkage Genes", "Colagen Genes", 
                                           "Laminin Genes", "Integrin Genes", 
                                           "MAPK Genes", "PIK3 Genes", 
                                           "Rac/Ras Genes", "ECM Genes", 
                                           "QC Marker Genes", "C.Gu VR Genes", 
                                           "C.Gu Pure DE Genes", "Greber's QC Marker Genes", 
                                           "Greber's Atrial Marker Genes", 
                                           "Greber's Ventricular Marker Genes", 
                                           "Greber's Atrial / Ventricular Marker Genes", 
                                           "Greber's Early iPSC-CM Marker Genes", 
                                           "Greber's Late iPSC-CM Marker Genes", 
                                           "Greber's Up-Regulated Cardio Marker Genes", 
                                           "Wu's Cardio Matration Marker Genes", 
                                           "DB C57 Marker Genes"), 
                               selected = "Cardiac Hypertrophy Marker Genes")),
                  br(),
                 
                 #textInput("GenesUI","Enter Custom Gene Set", ""),
                 #br(),
                 htmlOutput("HeatmapAnnotationUI_custom_geneset"),
                 br(),
                 textInput("CustomPCs","Principal Components to Plot", "1,2"),
                 br(),
                 textInput("CustomPCrotation","Principal Component Rotation for Download","PC1"),
                 br(),
                 downloadButton("CustomDownloadPCrotation", "Download PC Rotation"),
                 br(),
                 downloadButton("CustomDownloadPCloadings", "Download PC Loadings"),
                 br(),
                 br(),
          
                 actionButton("CustomGeneButton", strong("Update Plots"))
            ),
      
      ###### pvclust Panel ######
      conditionalPanel(condition="input.tabselected==4",
                       br(),
                       h3(strong("Set pvclust Parameters")),
                       htmlOutput("PCAGenesUI_pvclust"),
                       br(),
                       checkboxInput("CustomUI_pvclust", strong("Run pvclust with Custom Gene Set?"), FALSE),
                       conditionalPanel(
                         condition = "input.CustomUI_pvclust == true",
                         textInput("GenesUI_pvclust","Enter Custom Gene Set", value=NA)),
                       br(),
                       htmlOutput("HeatmapAnnotationUI_pvclust"),
                       br(),
                       textInput("nboots","Number of Bootstrap Replications", "100"),
                       br(),
                       actionButton("pvClustButton", strong("Update Plots"))
      )
      
    ),
    


################ MAIN Panel ###################
    
mainPanel(
  
  tabsetPanel(

    ###### QC Main ######
    tabPanel(strong("QC Plots"),
            value=1,
        
            h5(strong("Normalized Expression of House Keeping Genes")),
            plotOutput("HKGplot", width = "100%", height = "500px"),
            br(),
        
            h5(strong("PluriTest Plot")),
            plotOutput("PTplot", width = "100%", height = "600px"),
            br(),
            
            h5(strong("Read Normalization Distribution")),
            plotOutput("RTBoxPlot", width = "100%", height = "600px"),
            br(),
            
            h5(strong("Normalized Expression of QC Genes")),
            plotOutput("QCGplot", width = "100%", height = "1000px"),
            br(),
            
            h5(strong("Dendrogram")),      
            plotOutput("hclust", width = "100%", height = "600px"),
            br(),
            
            h5(strong("PCA plot")),      
            plotOutput("PCAplot", width = "100%", height = "700px"),
            DT::dataTableOutput("pc_anova"),
            br(),
            
            h5(strong("PCA 3Dplot (Click to Display)")),      
            #plotOutput("PCA3Dplot", width = "100%", height = "1000px"),
            plotlyOutput("PCA3Dplot"),
            br(),
            
            h5(strong("Heatmap of Top 100 Most Variable Genes")),
            plotOutput("genecluster100", width = "100%", height = "1000px"),
            br(),
            
            h5(strong("Heatmap of 500 Most Variable Genes")),
            plotOutput("genecluster500", width = "100%", height = "3000px"),
            br(),
            br()),
    
    ###### DE Main ######
    tabPanel(strong("DE Analysis"), 
               value=2,  
               h5(strong("Design Formula")),
               verbatimTextOutput("DesFormula"),
               
               h5(strong("Contrast Levels")),
               verbatimTextOutput("LevelsText"),
               br(),
             
               ## Read DESEQ Normalization
               h5(strong("DESeq2 Normalization Distribution")),
               plotOutput("DESEQBoxPlot", width = "100%", height = "600px"),
               br(),
             
               h5(strong("Results Table")),
               DT::dataTableOutput("ResTable"),
               br(),
               
               h5(strong("Results Summary")),
               shiny::tableOutput("res_summary"),
               br(),
               
               h5(strong("Top DE Genes")),
               h6(strong("(Over-Expressed)")),
               shiny::tableOutput("top_H_DEgenes"),
               h6(strong("(Under-Expressed)")),
               shiny::tableOutput("top_L_DEgenes"),
               br(),
               
               h5(strong("Gene Expression Plot")),
               br(),
               plotOutput("plot_gene", width = "100%", height = "400px"),
               
               verbatimTextOutput("GeneDesc"),
               br()),
    

    ###### Custom Gene Set ANalysis #####
    tabPanel(strong("Custom Gene Set"),
             value=3,
             h5(strong("Dendrogram using Custom Gene Set")),      
             plotOutput("hclust_Geneset", width = "100%"),
             br(),
             h5(strong("PCA Plot using Custom Gene Set")),      
             plotOutput("PCAplot_Geneset", width = "70%", height = "600px"),
             br(),
             h5(strong("Heatmap using Custom Gene Set")),
             #plotOutput("ui_geneset", width = "200%", height = "00px"),
             uiOutput("plot.ui"),
             br(),
             br()),
    
    ## pvclust ##
    tabPanel(strong("Run pvclust"),
             value=4,
             h5(strong("pvclust")),      
             plotOutput("pvclust", width = "100%", height = "700px"),
             br(),
             br()),
    
    #### Design Matrix Tab ####
    tabPanel("View Design Matrix",
             h5(strong("Design Matrix")),
             #shiny::tableOutput("DM"),
             DT::dataTableOutput("DM"),
             br()),
    
    
    #### Raw Counts Tab ####
    tabPanel("View Raw Counts",
             h5(strong("Raw Counts")),
             #shiny::tableOutput("RC"),
             DT::dataTableOutput("RC"),
             br()),
    
    # #### Phenotype Tab ####
    # tabPanel("View Phenotypes",
    #          h5(strong("Phenotypes")),
    #          DT::dataTableOutput("Pheno"),
    #          br()),
    
    
    id = "tabselected"
        
    )
      
      
      )
    
  )
))
