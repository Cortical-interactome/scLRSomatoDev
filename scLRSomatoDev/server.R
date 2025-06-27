#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

source('utils.R', local = TRUE)




# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  # Hide the "Loading..." tab once the app finishes loading in a consistent way
  observe({
    if (!is.null(session$clientData$url_search)) {
        hideTab("navbar", "Loading")
        showTab("navbar", "Overview", select = TRUE) 
    }
  })


  output$app_content <- renderUI({

    #####################
    # Overview
    #####################
    
    #####################
    #Datasets
    #####################
    
    ############ Table

      # --------------------------------------------------
      # Description: Render a DataTable output for the Cortex metadata
      # Input:
      #   - colData(GECortex_PostMsub): The metadata associated with the GECortex_PostMsub object
      # Output: A DataTable output with filtering, pagination, and download options
      # --------------------------------------------------
  
    output$Cortex_metadata <- renderDT(
        as.data.frame(colData(GECortex_PostMsub)),
        filter = "top",extensions = c('Buttons'),
        options = list(dom = 'Blfrtip',
                       pageLength=10, lengthMenu=c(10,20,50,100,500),
                       buttons = list(
                         list(
                         extend = 'collection',
                         buttons = list(
                           list(extend = 'csv', filename = "Metadata"),
                           list(extend = 'excel', filename = "Metadata")),
                         text = 'Download' )
                       )
                  )
    )
    
    ################# Data Visualization

      # --------------------------------------------------
      # Description: Update picker inputs for various data filters and selections
      # Input:
      #   - session: The Shiny session object
      #   - inputId: The ID of the picker input to update
      #   - choices: A vector of choices to display in the picker input
      #   - selected: The initially selected value(s) in the picker input
      # Output: Updates the picker inputs in the Shiny app UI
      # --------------------------------------------------
    
    #picker Input Dimension reduction
    observe({
        updatePickerInput(
            session = session,
            inputId = "Dim_red",
            choices = c("UMAP_2D", "UMAP_3D", "tSNE_2D", "tSNE_3D", "PCA_2D", "PCA_3D", "ICA_2D", "ICA_3D"),
            selected = "UMAP_2D"
        )
    })
    
    #picker Input grouping level
    observe({
        updatePickerInput(
            session = session,
            inputId = "DataVis_level",
            choices = c("class", "family","subclass", "supertype", "cell-type","original cell-type", "age","study","platform","RNA-seq method"),
            selected = "cell-type"
        )
    })   
    
    #picker Input class
    observe({
        updatePickerInput(
            session = session,
            inputId = "DataVis_class",
            choices = sort(unique(colData(GECortex_PostMsub)$class_label)),
            selected = sort(unique(colData(GECortex_PostMsub)$class_label))
        )
    }) 
    
    #picker Input family
    observe({
      updatePickerInput(
        session = session,
        inputId = "DataVis_family",
        choices = sort(unique(colData(GECortex_PostMsub)$family_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$family_label))
      )
    })
    
    #picker Input subclass
    observe({
        updatePickerInput(
            session = session,
            inputId = "DataVis_subclass",
            choices = sort(unique(colData(GECortex_PostMsub)$subclass_label)),
            selected = sort(unique(colData(GECortex_PostMsub)$subclass_label))
        )
    }) 
    
    #picker Input supertype
    observe({
      updatePickerInput(
        session = session,
        inputId = "DataVis_supertype",
        choices = sort(unique(colData(GECortex_PostMsub)$supertype_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$supertype_label))
      )
    }) 
    
    #picker Input celltype
    observe({
        updatePickerInput(
            session = session,
            inputId = "DataVis_celltype",
            choices = sort(unique(colData(GECortex_PostMsub)$celltype_label)),
            selected = sort(unique(colData(GECortex_PostMsub)$celltype_label))
        )
    }) 
    
    #picker Input celltype_original
    observe({
        updatePickerInput(
            session = session,
            inputId = "DataVis_celltype_original",
            choices = sort(unique(colData(GECortex_PostMsub)$celltype_original_label)),
            selected = sort(unique(colData(GECortex_PostMsub)$celltype_original_label))
        )
    }) 
    
    #picker Input age
    observe({
      updatePickerInput(
        session = session,
        inputId = "DataVis_age",
        choices = sort(unique(colData(GECortex_PostMsub)$age.at.collection.grp)),
        selected = sort(unique(colData(GECortex_PostMsub)$age.at.collection.grp))
      )
    }) 
    
    
    #picker Input region
    observe({
      updatePickerInput(
        session = session,
        inputId = "DataVis_region",
        choices = sort(unique(colData(GECortex_PostMsub)$region_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$region_label))
      )
    }) 

    
    #picker Input study
    observe({
        updatePickerInput(
            session = session,
            inputId = "DataVis_study",
            choices = sort(unique(colData(GECortex_PostMsub)$study_label)),
            selected = sort(unique(colData(GECortex_PostMsub)$study_label))
        )
    }) 
    
    #picker Input platform
    observe({
      updatePickerInput(
        session = session,
        inputId = "DataVis_platform",
        choices = sort(unique(colData(GECortex_PostMsub)$platform_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$platform_label))
      )
    }) 
    
    #picker Input RNAseq.method
    observe({
      updatePickerInput(
        session = session,
        inputId = "DataVis_RNAseq.method",
        choices = sort(unique(colData(GECortex_PostMsub)$RNAseq.method_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$RNAseq.method_label))
      )
    }) 
    
    
    #picker Input Genes
    observe({
        updatePickerInput(
            session = session,
            inputId = "DataVis_gene",
            choices = sort(unique(rownames(GECortex_PostMsub))),
            selected = c("Neurod2", "Gad1")
        )
    }) 
    
    
    
    
    ########Dim red
    
      # --------------------------------------------------
      # Description: Render a Plotly output for dimensionality reduction of the Cortex data
      # Input:
      #   - GECortex_PostMsub: The SingleCellExperiment object containing the Cortex data
      #   - Various input$ parameters: User-selected filtering and visualization options
      # Output: A Plotly plot showing the dimensionality reduction of the filtered Cortex data
      # --------------------------------------------------

    
    output$Cortex_Dim.red <- renderPlotly({

            # Call the plot_Dim.red function with user-selected filtering and visualization options
            plot_Dim.red(GECortex_PostMsub, grouping.level= input$"DataVis_level", dim.red.type=input$"Dim_red", ident.class=input$"DataVis_class", 
                         ident.family=input$"DataVis_family", ident.subclass=input$"DataVis_subclass",ident.supertype=input$"DataVis_supertype", ident.celltype=input$"DataVis_celltype", 
                         ident.celltype_original=input$"DataVis_celltype_original",ident.age=input$"DataVis_age",ident.study=input$"DataVis_study",
                         ident.region=input$"DataVis_region", ident.RNAseq.method=input$"DataVis_RNAseq.method" ,ident.platform=input$"DataVis_platform",
                         height=input$"DataVis_height")
        
                                        })
    
      # --------------------------------------------------
      # Description: Render a Plotly output for dimensionality reduction plot with gene expression overlay in the Cortex data
      # Input:
      #   - input$"Dimred_button": The action button value (not used in this context, but required to trigger the plot rendering)
      #   - Various input$ parameters: User-selected filtering and visualization options
      # Output: A Plotly plot showing a dimensionality reduction plot with gene expression overlay for the filtered Cortex data
      # --------------------------------------------------

    output$Cortex_Dim.red_ExprGene <- renderPlotly({

        # Read the input button value (not used in this context, but required to trigger the plot rendering)
        input$"Dimred_button"

        # Render the dimensionality reduction plot with gene expression overlay
        plot_Dim.red_ExprGene(GECortex_PostMsub, "Data/logCPM_Cortex_PostMsub.tome", grouping.level= input$"DataVis_level", dim.red.type=input$"Dim_red", 
                              ident.class=input$"DataVis_class", ident.family=input$"DataVis_family", 
                              ident.subclass=input$"DataVis_subclass",ident.supertype=input$"DataVis_supertype", ident.celltype=input$"DataVis_celltype", 
                              ident.celltype_original=input$"DataVis_celltype_original",ident.age=input$"DataVis_age",ident.study=input$"DataVis_study",
                              ident.region=input$"DataVis_region", ident.RNAseq.method=input$"DataVis_RNAseq.method" , ident.platform=input$"DataVis_platform" ,
                              height=input$"DataVis_height",GeneList=isolate(input$"DataVis_gene"), color_gdt1=input$"DataVis_pal1", color_gdt2=input$"DataVis_pal2")
        
                                                 })

    
    ################ Gene Expression
    
    
    ############# Heatmap by grouping level

    # --------------------------------------------------
    # Description: Update picker inputs for various data filters and selections
    # Input:
    #   - session: The Shiny session object
    #   - inputId: The ID of the picker input to update
    #   - choices: A vector of choices to display in the picker input
    #   - selected: The initially selected value(s) in the picker input
    # Output: Updates the picker inputs in the Shiny app UI
    # --------------------------------------------------

    #picker Input Plot type
    observe({
        updatePickerInput(
            session = session,
            inputId = "ExprData_plot",
            choices = c("group_heatmap_plot", "group_dot_plot"),
            selected = "group_dot_plot"
        )
    })
    
    #picker Input grouping level
    observe({
        updatePickerInput(
            session = session,
            inputId = "ExprData_level",
            choices = c("class","family", "subclass", "celltype", "celltype_original"),
            selected = "celltype"
        )
    })   
    
    #picker Input class
    observe({
        updatePickerInput(
            session = session,
            inputId = "ExprData_class",
            choices = sort(unique(colData(GECortex_PostMsub)$class_label)),
            selected = sort(unique(colData(GECortex_PostMsub)$class_label))
        )
    }) 
    
    #picker Input family
    observe({
      updatePickerInput(
        session = session,
        inputId = "ExprData_family",
        choices = sort(unique(colData(GECortex_PostMsub)$family_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$family_label))
      )
    }) 
    
    
    #picker Input subclass
    observe({
        updatePickerInput(
            session = session,
            inputId = "ExprData_subclass",
            choices = sort(unique(colData(GECortex_PostMsub)$subclass_label)),
            selected = sort(unique(colData(GECortex_PostMsub)$subclass_label))
        )
    }) 
    

    #picker Input celltype
    observe({
        updatePickerInput(
            session = session,
            inputId = "ExprData_celltype",
            choices = sort(unique(colData(GECortex_PostMsub)$celltype_label)),
            selected = sort(unique(colData(GECortex_PostMsub)$celltype_label))
        )
    }) 
    
    #picker Input celltype_original
    observe({
        updatePickerInput(
            session = session,
            inputId = "ExprData_celltype_original",
            choices = sort(unique(colData(GECortex_PostMsub)$celltype_original_label)),
            selected = sort(unique(colData(GECortex_PostMsub)$celltype_original_label))
        )
    }) 
    
    #picker Input age
    observe({
      updatePickerInput(
        session = session,
        inputId = "ExprData_age",
        choices = sort(unique(colData(GECortex_PostMsub)$age.at.collection.grp)),
        selected = sort(unique(colData(GECortex_PostMsub)$age.at.collection.grp))
      )
    }) 
    
    #picker Input study
    observe({
      updatePickerInput(
        session = session,
        inputId = "ExprData_study",
        choices = sort(unique(colData(GECortex_PostMsub)$study_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$study_label))
      )
    }) 
    
    #picker Input region
    observe({
      updatePickerInput(
        session = session,
        inputId = "ExprData_region",
        choices = sort(unique(colData(GECortex_PostMsub)$region_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$region_label))
      )
    }) 
    
    #picker Input RNA.seq_method
    observe({
      updatePickerInput(
        session = session,
        inputId = "ExprData_RNAseq.method",
        choices = sort(unique(colData(GECortex_PostMsub)$RNAseq.method_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$RNAseq.method_label))
      )
    }) 
    
    #picker Input platform
    observe({
      updatePickerInput(
        session = session,
        inputId = "ExprData_platform",
        choices = sort(unique(colData(GECortex_PostMsub)$platform_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$platform_label))
      )
    }) 
    
    #picker Input Selected Genes
    observe({
        updatePickerInput(
            session = session,
            inputId = "ExprData_gene",
            choices = sort(unique(rownames(GECortex_PostMsub))),
            selected = c("Neurod2", "Gad1", "Gad2", "Dlx1", "Dlx2", "Trp73")
        )
    }) 
    
    
    
    
    ########Cortex

      # --------------------------------------------------
      # Description: Render a Plotly output for gene expression data of the Cortex data
      # Input:
      #   - GECortex_PostMsub: The SingleCellExperiment object containing the Cortex data
      #   - Various input$ parameters: User-selected filtering and visualization options
      # Output: A Plotly plot showing either a grouped heatmap or a grouped dot plot of the filtered gene expression data
      # --------------------------------------------------

    output$Cortex_ExprData <- renderPlotly({

       # Read the input button value (not used in this context, but required to trigger the plot rendering)
        input$"ExprData_button"

        # Check the selected plot type and render the corresponding plot
        if(!!input$"ExprData_plot"=="group_heatmap_plot"){

        # Render a grouped heatmap plot  
        plotly_heatmap_markers(GECortex_PostMsub,"Data/logCPM_Cortex_PostMsub.tome", grouping= input$"ExprData_level",genes=isolate(input$"ExprData_gene"),
                               colorset=input$"ExprData_pal", ident.class=input$"ExprData_class", ident.family=isolate(input$"ExprData_family"),
                               ident.subclass=isolate(input$"ExprData_subclass"), ident.celltype=isolate(input$"ExprData_celltype"), 
                               ident.celltype_original=isolate(input$"ExprData_celltype_original") ,ident.age=isolate(input$"ExprData_age"),
                               ident.study=isolate(input$"ExprData_study"),ident.region=isolate(input$"ExprData_region"), 
                               ident.RNAseq.method=isolate(input$"ExprData_RNAseq.method"),
                               ident.platform=isolate(input$"ExprData_platform"),
                               height=input$"ExprData_height")
        }
        
        else if(!!input$"ExprData_plot"=="group_dot_plot"){

        # Render a grouped dot plot
        plotly_dot_markers(GECortex_PostMsub,"Data/logCPM_Cortex_PostMsub.tome", grouping= input$"ExprData_level",genes=isolate(input$"ExprData_gene"),
                               colorset=input$"ExprData_pal", ident.class=input$"ExprData_class", ident.family=input$"ExprData_family", 
                               ident.subclass=input$"ExprData_subclass",
                               ident.celltype=input$"ExprData_celltype", 
                               ident.celltype_original=input$"ExprData_celltype_original" ,ident.age=input$"ExprData_age",
                               ident.study=input$"ExprData_study",ident.region=input$"ExprData_region", ident.RNAseq.method=input$"ExprData_RNAseq.method",
                               ident.platform=input$"ExprData_platform",
                               height=input$"ExprData_height")
        }
        
        
    })

    # --------------------------------------------------
    # Description: Render a text output showing the number of selected genes and the number of detected genes in the dataset
    # Input:
    #   - input$"ExprData_button": The action button value (not used in this context, but required to trigger the text rendering)
    #   - input$"ExprData_gene": The vector of selected genes
    # Output: A text output displaying the number of selected genes and the number of detected genes in the dataset
    # --------------------------------------------------
    
    output$ExprData_selectLR <- renderText({
        input$"ExprData_button"
        paste("You have selected", length(isolate(input$"ExprData_gene")), "genes" ,"//", 
              length(rownames(GECortex_PostMsub)[rownames(GECortex_PostMsub)%in%isolate(input$"ExprData_gene")]), "detected in the dataset")
    })
    
    ############# Heatmap by cell-type*Age

      # --------------------------------------------------
      # Description: Update picker inputs for various data filters and selections
      # Input:
      #   - session: The Shiny session object
      #   - inputId: The ID of the picker input to update
      #   - choices: A vector of choices to display in the picker input
      #   - selected: The initially selected value(s) in the picker input
      # Output: Updates the picker inputs in the Shiny app UI
      # --------------------------------------------------
    
    #picker Input grouping level
    observe({
      updatePickerInput(
        session = session,
        inputId = "HmctAge_level",
        choices = c("class","family", "subclass", "celltype", "celltype_original"),
        selected = "celltype"
      )
    })   
    
    
    #picker Input class
    
    observe({
      updatePickerInput(
        session = session,
        inputId = "HmctAge_class",
        choices = sort(unique(colData(GECortex_PostMsub)$class_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$class_label))
      )
    }) 
    
    
    #picker Input class
    observe({
      updatePickerInput(
        session = session,
        inputId = "HmctAge_family",
        choices = sort(unique(colData(GECortex_PostMsub)$family_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$family_label))
      )
    }) 
    
    #picker Input subclass
    observe({
      updatePickerInput(
        session = session,
        inputId = "HmctAge_subclass",
        choices = sort(unique(colData(GECortex_PostMsub)$subclass_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$subclass_label))
      )
    }) 
    
    
    #picker Input celltype
    observe({
      updatePickerInput(
        session = session,
        inputId = "HmctAge_celltype",
        choices = sort(unique(colData(GECortex_PostMsub)$celltype_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$celltype_label))
      )
    }) 
    
    #picker Input celltype_original
    observe({
      updatePickerInput(
        session = session,
        inputId = "HmctAge_celltype_original",
        choices = sort(unique(colData(GECortex_PostMsub)$celltype_original_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$celltype_original_label))
      )
    }) 
    
    #picker Input age
    observe({
      updatePickerInput(
        session = session,
        inputId = "HmctAge_age",
        choices = sort(unique(colData(GECortex_PostMsub)$age.at.collection.grp)),
        selected = sort(unique(colData(GECortex_PostMsub)$age.at.collection.grp))
      )
    }) 
    
    #picker Input study
    observe({
      updatePickerInput(
        session = session,
        inputId = "HmctAge_study",
        choices = sort(unique(colData(GECortex_PostMsub)$study_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$study_label))
      )
    }) 
    
    #picker Input region
    observe({
      updatePickerInput(
        session = session,
        inputId = "HmctAge_region",
        choices = sort(unique(colData(GECortex_PostMsub)$region_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$region_label))
      )
    }) 
    
    #picker Input RNA.seq_method
    observe({
      updatePickerInput(
        session = session,
        inputId = "HmctAge_RNAseq.method",
        choices = sort(unique(colData(GECortex_PostMsub)$RNAseq.method_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$RNAseq.method_label))
      )
    }) 
    
    #picker Input platform
    observe({
      updatePickerInput(
        session = session,
        inputId = "HmctAge_platform",
        choices = sort(unique(colData(GECortex_PostMsub)$platform_label)),
        selected = sort(unique(colData(GECortex_PostMsub)$platform_label))
      )
    }) 
    
    #picker Input Selected Genes
    observe({
      updatePickerInput(
        session = session,
        inputId = "HmctAge_gene",
        choices = sort(unique(rownames(GECortex_PostMsub))),
        selected = c("Neurod2", "Gad1", "Gad2", "Dlx1", "Dlx2", "Trp73")
      )
    }) 
    
    
    
    
    ########Cortex
    
      # --------------------------------------------------
      # Description: Render a Plotly output for gene expression heatmap by age in the Cortex data
      # Input:
      #   - input$"HmctAge_button": The action button value (not used in this context, but required to trigger the plot rendering)
      #   - Various input$ parameters: User-selected filtering and visualization options
      # Output: A Plotly plot showing a gene expression heatmap by age for the filtered Cortex data
      # --------------------------------------------------
    
    output$Cortex_HmctAge <- renderPlotly({

      # Read the input button value (not used in this context, but required to trigger the plot rendering)
      input$"HmctAge_button"

        # Render the gene expression heatmap by age plot
        plotlyHm_Exprgene_byAge(GECortex_PostMsub,"Data/logCPM_Cortex_PostMsub.tome",gene=isolate(input$"HmctAge_gene"),grouping=input$"HmctAge_level",
                               colorset=input$"HmctAge_pal", ident.class=input$"HmctAge_class",ident.family=isolate(input$"HmctAge_family"),
                               ident.subclass=isolate(input$"HmctAge_subclass"), 
                               ident.celltype=isolate(input$"HmctAge_celltype"),ident.age=isolate(input$"HmctAge_age"),
                               ident.study=isolate(input$"HmctAge_study"),ident.region=isolate(input$"HmctAge_region"), 
                               ident.RNAseq.method=isolate(input$"HmctAge_RNAseq.method"),
                               ident.platform=isolate(input$"HmctAge_platform"),
                               height=input$"HmctAge_height")
      
    })
    
    # --------------------------------------------------
    # Description: Render a text output showing the number of selected genes and the number of detected genes in the dataset
    # Input:
    #   - input$"HmctAge_button": The action button value (not used in this context, but required to trigger the text rendering)
    #   - input$"HmctAge_gene": The vector of selected genes
    # Output: A text output displaying the number of selected genes and the number of detected genes in the dataset
    # --------------------------------------------------

    output$HmctAge_selectGene <- renderText({
      input$"HmctAge_button"
      paste("You have selected", length(isolate(input$"HmctAge_gene")), "genes" ,"//", 
            length(rownames(GECortex_PostMsub)[rownames(GECortex_PostMsub)%in%isolate(input$"HmctAge_gene")]), "detected in the dataset")
    })
    

    ############# Pseudo-maturation

      # --------------------------------------------------
      # Description: Update picker inputs for various data filters and selections
      # Input:
      #   - session: The Shiny session object
      #   - inputId: The ID of the picker input to update
      #   - choices: A vector of choices to display in the picker input
      #   - selected: The initially selected value(s) in the picker input
      # Output: Updates the picker inputs in the Shiny app UI
      # --------------------------------------------------
    
    #picker Input celltype
    observe({
      updatePickerInput(
        session = session,
        inputId = "pM_celltype",
        choices = sort(names(Res_pt.list)),
        selected = names(Res_pt.list)
      )
    }) 
    
    #picker Input Selected Genes
    observe({
      updatePickerInput(
        session = session,
        inputId = "pM_gene",
        choices = Gene.pM,
        selected = c("Tubb3", "Dcx", "Dlg4", "Gphn", "Rbfox3", "Syn1")
      )
    }) 
    
    
    ########Plot 
    
      # --------------------------------------------------
      # Description: Render a plotly output for gene expression heatmap and spline plot by cell type in the Res_pt.list data
      # Input:
      #   - input$"pM_button": The action button value (not used in this context, but required to trigger the plot rendering)
      #   - Various input$ parameters: User-selected filtering and visualization options
      # Output: Two plotly plots showing a gene expression heatmap and a spline plot by cell type for the filtered Res_pt.list data
      # --------------------------------------------------

    # Dynamically adjust the output height and width based on the number of selected genes  
    observe({
      output$Hm_pM <- renderPlot({

        # Read the input button value (not used in this context, but required to trigger the plot rendering)
        input$"pM_button"
        
        # Render the gene expression heatmap by cell type plot
        plotGene_Hm_pt(Res_pt.list,celltype=input$"pM_celltype", gene = isolate(input$"pM_gene"), free.scale = TRUE,  ncol = NA, 
                       col_gdt=input$"pM_pal")
      }, 

      # Set the output height and width based on the number of selected genes
      height=input$"pM_height"+10*length(isolate(input$"pM_gene")), width=input$"pM_width"+10*length(isolate(input$"pM_gene")))
      
    })
    
    # Render the gene expression spline plot by cell type
    output$Spline_pM <- renderPlotly({
      input$"pM_button"
        
        plotlyGene_smoothSpline_pt(Res_pt.list,celltype=input$"pM_celltype", gene = isolate(input$"pM_gene"), free.scale = TRUE,  ncol = NA, 
                                   height=input$"pM_height", width=input$"pM_width")
      
    })

    # --------------------------------------------------
    # Description: Render text outputs showing the number of selected genes and the number of detected genes in the dataset for the heatmap and spline plots
    # Input:
    #   - input$"pM_button": The action button value (not used in this context, but required to trigger the text rendering)
    #   - input$"pM_gene": The vector of selected genes
    # Output: Two text outputs displaying the number of selected genes and the number of detected genes in the dataset for the heatmap and spline plots
    # --------------------------------------------------
    
    # Render text output for heatmap plot
    output$pM_selectGene_Hm <- renderText({
      input$"pM_button"
      paste("You have selected", length(isolate(input$"pM_gene")), "genes" ,"//", 
            length(Gene.pM[Gene.pM%in%isolate(input$"pM_gene")]), "detected in the dataset")
    })

    # Render text output for spline plot
    output$pM_selectGene_Spline <- renderText({
      input$"pM_button"
      paste("You have selected", length(isolate(input$"pM_gene")), "genes" ,"//", 
            length(Gene.pM[Gene.pM%in%isolate(input$"pM_gene")]), "detected in the dataset")
    })
    
    ############# Pseudo-layer
    
      # --------------------------------------------------
      # Description: Update picker inputs for various data filters and selections
      # Input:
      #   - session: The Shiny session object
      #   - inputId: The ID of the picker input to update
      #   - choices: A vector of choices to display in the picker input
      #   - selected: The initially selected value(s) in the picker input
      # Output: Updates the picker inputs in the Shiny app UI
      # --------------------------------------------------

    #picker Input celltype
    observe({
      updatePickerInput(
        session = session,
        inputId = "pL_celltype",
        choices = c("IT","ET","Sst","Pvalb","Vip","Lamp5"),
        selected = "IT"
      )
    }) 
    
    #picker Input age
    observe({
      updatePickerInput(
        session = session,
        inputId = "pL_age",
        choices = names(Res_pL.list$IT),
        selected = names(Res_pL.list$IT)
      )
    }) 
    
    #picker Input Selected Genes
    observe({
      updatePickerInput(
        session = session,
        inputId = "pL_gene",
        choices = Gene.pL,
        selected = c( "Cux2", "Rorb", "Hs3st2","Cryab")
      )
    }) 
    
    
    ########Plot 
    
      # --------------------------------------------------
      # Description: Render a plotly output for gene expression heatmap by cell type and age in the Res_pL.list data
      # Input:
      #   - input$"pL_button": The action button value (not used in this context, but required to trigger the plot rendering)
      #   - Various input$ parameters: User-selected filtering and visualization options
      # Output: A plotly plot showing a gene expression heatmap by cell type and age for the filtered Res_pL.list data
      # --------------------------------------------------

    # Dynamically adjust the output height and width based on the number of selected genes
    observe({
      output$Hm_pL <- renderPlot({   
        input$"pL_button"

         # Render the gene expression heatmap by cell type and age plot based on the selected cell type
        if(!!input$"pL_celltype"=="IT"){
          plotGene_Hm_pL(Res_pL.list$IT,age=input$"pL_age", gene = isolate(input$"pL_gene"), free.scale = TRUE,  ncol = NA, 
                         col_gdt=input$"pL_pal")
        }
        else if(!!input$"pL_celltype"=="ET"){
          plotGene_Hm_pL(Res_pL.list$ET,age=input$"pL_age", gene =  isolate(input$"pL_gene"), free.scale = TRUE,  ncol = NA, 
                         col_gdt=input$"pL_pal")
        }
        else if(!!input$"pL_celltype"=="Sst"){
         plotGene_Hm_pL(Res_pL.list$Sst,age=input$"pL_age", gene = isolate(input$"pL_gene"), free.scale = TRUE,  ncol = NA, 
                        col_gdt=input$"pL_pal")
        }
        else if(!!input$"pL_celltype"=="Pvalb"){
          plotGene_Hm_pL(Res_pL.list$Pvalb,age=input$"pL_age", gene = isolate(input$"pL_gene"), free.scale = TRUE,  ncol = NA, 
                         col_gdt=input$"pL_pal")
        }
        else if(!!input$"pL_celltype"=="Vip"){
          plotGene_Hm_pL(Res_pL.list$Vip,age=input$"pL_age", gene = isolate(input$"pL_gene"), free.scale = TRUE,  ncol = NA, 
                         col_gdt=input$"pL_pal")
        }
        else if(!!input$"pL_celltype"=="Lamp5"){
          plotGene_Hm_pL(Res_pL.list$Lamp5,age=input$"pL_age", gene = isolate(input$"pL_gene"), free.scale = TRUE,  ncol = NA, 
                         col_gdt=input$"pL_pal")
        }
        
        # Set the output height and width based on the number of selected genes
      }, height=input$"pL_height"+10*length(isolate(input$"pL_gene")), width=input$"pL_width"+10*length(isolate(input$"pL_gene")))
      
    })
    
    output$Spline_pL <- renderPlotly({
      input$"pL_button"

        # Render the gene expression spline plot by cell type and age based on the selected cell type
        if(!!input$"pL_celltype"=="IT"){
        plotlyGene_smoothSpline_pL(Res_pL.list$IT,age=input$"pL_age", gene =isolate(input$"pL_gene"), free.scale = TRUE,  ncol = NA, 
                                   line.size = 1,  showlegend=F, col_gdt=input$"pL_spline_pal", height=input$"pL_height",width=input$"pL_width")
        }
        else if(!!input$"pL_celltype"=="ET"){
          plotlyGene_smoothSpline_pL(Res_pL.list$ET,age=input$"pL_age", gene = isolate(input$"pL_gene"), free.scale = TRUE,  ncol = NA, 
                                    line.size = 1,  showlegend=F, col_gdt=input$"pL_spline_pal", height=input$"pL_height",width=input$"pL_width")
        }
        else if(!!input$"pL_celltype"=="Sst"){
          plotlyGene_smoothSpline_pL(Res_pL.list$Sst,age=input$"pL_age", gene = isolate(input$"pL_gene"), free.scale = TRUE,  ncol = NA, 
                                    line.size = 1,  showlegend=F, col_gdt=input$"pL_spline_pal", height=input$"pL_height",width=input$"pL_width")
        }
        else if(!!input$"pL_celltype"=="Pvalb"){
          plotlyGene_smoothSpline_pL(Res_pL.list$Pvalb,age=input$"pL_age", gene = isolate(input$"pL_gene"), free.scale = TRUE,  ncol = NA, 
                                     line.size = 1,  showlegend=F, col_gdt=input$"pL_spline_pal", height=input$"pL_height",width=input$"pL_width")
        }
      else if(!!input$"pL_celltype"=="Vip"){
        plotlyGene_smoothSpline_pL(Res_pL.list$Vip,age=input$"pL_age", gene = isolate(input$"pL_gene"), free.scale = TRUE,  ncol = NA, 
                                   line.size = 1,  showlegend=F, col_gdt=input$"pL_spline_pal", height=input$"pL_height",width=input$"pL_width")
      }
      else if(!!input$"pL_celltype"=="Lamp5"){
        plotlyGene_smoothSpline_pL(Res_pL.list$Lamp5,age=input$"pL_age", gene = isolate(input$"pL_gene"), free.scale = TRUE,  ncol = NA, 
                                   line.size = 1,  showlegend=F, col_gdt=input$"pL_spline_pal", height=input$"pL_height",width=input$"pL_width")
      }

})
    
    # --------------------------------------------------
    # Description: Render text outputs showing the number of selected genes and the number of detected genes in the dataset for the heatmap and spline plots based on the selected cell type
    # Input:
    #   - input$"pL_button": The action button value (not used in this context, but required to trigger the text rendering)
    #   - Various input$ parameters: User-selected filtering and visualization options
    # Output: Two text outputs displaying the number of selected genes and the number of detected genes in the dataset for the heatmap and spline plots
    # --------------------------------------------------
 
    
    output$pL_selectGene_Hm <- renderText({
      input$"pL_button"
      
      # Render the text output for the heatmap plot based on the selected cell type
      if(!!input$"pL_celltype"=="IT"){
      paste("You have selected", length(isolate(input$"pL_gene")), "genes" ,"//", 
            length(Gene_IT.pL[Gene_IT.pL%in%isolate(input$"pL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"pL_celltype"=="ET"){
        paste("You have selected", length(isolate(input$"pL_gene")), "genes" ,"//", 
              length(Gene_ET.pL[Gene_ET.pL%in%isolate(input$"pL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"pL_celltype"=="Sst"){
        paste("You have selected", length(isolate(input$"pL_gene")), "genes" ,"//", 
              length(Gene_Sst.pL[Gene_Sst.pL%in%isolate(input$"pL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"pL_celltype"=="Pvalb"){
        paste("You have selected", length(isolate(input$"pL_gene")), "genes" ,"//", 
              length(Gene_Pvalb.pL[Gene_Pvalb.pL%in%isolate(input$"pL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"pL_celltype"=="Vip"){
        paste("You have selected", length(isolate(input$"pL_gene")), "genes" ,"//", 
              length(Gene_Vip.pL[Gene_Vip.pL%in%isolate(input$"pL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"pL_celltype"=="Lamp5"){
        paste("You have selected", length(isolate(input$"pL_gene")), "genes" ,"//", 
              length(Gene_Lamp5.pL[Gene_Lamp5.pL%in%isolate(input$"pL_gene")]), "detected in the dataset")
      }
      
    })
    
    
    
    output$pL_selectGene_Spline <- renderText({
      input$"pL_button"
      
      # Render the text output for the spline plot based on the selected cell type
      if(!!input$"pL_celltype"=="IT"){
        paste("You have selected", length(isolate(input$"pL_gene")), "genes" ,"//", 
              length(Gene_IT.pL[Gene_IT.pL%in%isolate(input$"pL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"pL_celltype"=="ET"){
        paste("You have selected", length(isolate(input$"pL_gene")), "genes" ,"//", 
              length(Gene_ET.pL[Gene_ET.pL%in%isolate(input$"pL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"pL_celltype"=="Sst"){
        paste("You have selected", length(isolate(input$"pL_gene")), "genes" ,"//", 
              length(Gene_Sst.pL[Gene_Sst.pL%in%isolate(input$"pL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"pL_celltype"=="Pvalb"){
        paste("You have selected", length(isolate(input$"pL_gene")), "genes" ,"//", 
              length(Gene_Pvalb.pL[Gene_Pvalb.pL%in%isolate(input$"pL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"pL_celltype"=="Vip"){
        paste("You have selected", length(isolate(input$"pL_gene")), "genes" ,"//", 
              length(Gene_Vip.pL[Gene_Vip.pL%in%isolate(input$"pL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"pL_celltype"=="Lamp5"){
        paste("You have selected", length(isolate(input$"pL_gene")), "genes" ,"//", 
              length(Gene_Lamp5.pL[Gene_Lamp5.pL%in%isolate(input$"pL_gene")]), "detected in the dataset")
      }
      
    })
 

    ############# Transcriptional landscape
    
      # --------------------------------------------------
      # Description: Update picker inputs for various data filters and selections
      # Input:
      #   - session: The Shiny session object
      #   - inputId: The ID of the picker input to update
      #   - choices: A vector of choices to display in the picker input
      #   - selected: The initially selected value(s) in the picker input
      # Output: Updates the picker inputs in the Shiny app UI
      # --------------------------------------------------

    #picker Input celltype
    observe({
      updatePickerInput(
        session = session,
        inputId = "TL_celltype",
        choices = c("IT","ET","Sst","Pvalb","Vip","Lamp5","Other GlutNs","Other GABANs"),
        selected = "IT"
      )
    }) 
    
    #picker Input Selected Genes
    observe({
      updatePickerInput(
        session = session,
        inputId = "TL_gene",
        choices = Gene_TL,
        selected = c( "Cux2", "Rorb", "Hs3st2","Cryab")
      )
    }) 
    
    ########Plot 

      # --------------------------------------------------
      # Description: Render a plot output for the trajectory landscape based on the selected cell type and genes
      # Input:
      #   - input$"TL_button": The action button value (not used in this context, but required to trigger the plot rendering)
      #   - Various input$ parameters: User-selected filtering and visualization options
      # Output: A plot output showing the trajectory landscape based on the selected cell type and genes
      # --------------------------------------------------
    
    observe({
      output$Cortex_TL <- renderPlot({

      # Read the input button value (not used in this context, but required to trigger the plot rendering)
      input$"TL_button"
        
        # Render the trajectory landscape plot based on the selected cell type
        if(!!input$"TL_celltype"=="IT"){
          Plot.Tlandscape(celltype = "IT", col_gdt=input$"TL_pal",gene=isolate(input$"TL_gene"),sig.gene=thr_name.family.list[[c("IT")]])
        }
        else if(!!input$"TL_celltype"=="ET"){
         Plot.Tlandscape(celltype = "ET", col_gdt=input$"TL_pal",gene=isolate(input$"TL_gene"),sig.gene=thr_name.family.list[[c("ET")]])
          
        }
        else if(!!input$"TL_celltype"=="Sst"){
          Plot.Tlandscape(celltype = "Sst", col_gdt=input$"TL_pal",gene=isolate(input$"TL_gene"),sig.gene=thr_name.family.list[[c("Sst")]])
          
        }
        else if(!!input$"TL_celltype"=="Pvalb"){
          Plot.Tlandscape(celltype = "Pvalb", col_gdt=input$"TL_pal",gene=isolate(input$"TL_gene"),sig.gene=thr_name.family.list[[c("Pvalb")]])
          
        }
        else if(!!input$"TL_celltype"=="Vip"){
          Plot.Tlandscape(celltype = "Vip", col_gdt=input$"TL_pal",gene=isolate(input$"TL_gene"),sig.gene=thr_name.family.list[[c("Vip")]])
          
        }
        else if(!!input$"TL_celltype"=="Lamp5"){
          Plot.Tlandscape(celltype = "Lamp5", col_gdt=input$"TL_pal",gene=isolate(input$"TL_gene"),sig.gene=thr_name.family.list[[c("Lamp5")]])
          
        }
        
        else if(!!input$"TL_celltype"=="Other GlutNs"){
          Plot.Tlandscape(celltype = "Other GlutNs", col_gdt=input$"TL_pal",gene=isolate(input$"TL_gene"),sig.gene=thr_name.family.list[[c("Other GlutNs")]])
          
        }
        
        else if(!!input$"TL_celltype"=="Other GABANs"){
          Plot.Tlandscape(celltype = "Other GABANs", col_gdt=input$"TL_pal",gene=isolate(input$"TL_gene"),sig.gene=thr_name.family.list[[c("Other GABANs")]])
          
        }
      
       # Set the output height and width based on the number of selected genes 
    },height=input$"TL_height"+10*length(isolate(input$"TL_gene")), width=input$"TL_width"+10*length(isolate(input$"TL_gene")))

      
    })
    
      # --------------------------------------------------
      # Description: Render a text output showing the number of selected genes and the number of detected genes in the dataset based on the selected cell type
      # Input:
      #   - input$"TL_button": The action button value (not used in this context, but required to trigger the text rendering)
      #   - Various input$ parameters: User-selected filtering and visualization options
      # Output: A text output showing the number of selected genes and the number of detected genes in the dataset based on the selected cell type
      # --------------------------------------------------

    
    #Number of genes selected
    output$TL_selectGene <- renderText({
      input$"TL_button"
      
      # Render the text output based on the selected cell type
      if(!!input$"TL_celltype"=="IT"){
        paste("You have selected", length(isolate(input$"TL_gene")), "genes" ,"//", 
              length(Gene_IT.TL[Gene_IT.TL%in%isolate(input$"TL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"TL_celltype"=="ET"){
        paste("You have selected", length(isolate(input$"TL_gene")), "genes" ,"//", 
              length(Gene_ET.TL[Gene_ET.TL%in%isolate(input$"TL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"TL_celltype"=="Sst"){
        paste("You have selected", length(isolate(input$"TL_gene")), "genes" ,"//", 
              length(Gene_Sst.TL[Gene_Sst.TL%in%isolate(input$"TL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"TL_celltype"=="Pvalb"){
        paste("You have selected", length(isolate(input$"TL_gene")), "genes" ,"//", 
              length(Gene_Pvalb.TL[Gene_Pvalb.TL%in%isolate(input$"TL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"TL_celltype"=="Vip"){
        paste("You have selected", length(isolate(input$"TL_gene")), "genes" ,"//", 
              length(Gene_Vip.TL[Gene_Vip.TL%in%isolate(input$"TL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"TL_celltype"=="Lamp5"){
        paste("You have selected", length(isolate(input$"TL_gene")), "genes" ,"//", 
              length(Gene_Lamp5.TL[Gene_Lamp5.TL%in%isolate(input$"TL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"TL_celltype"=="Other GlutNs"){
        paste("You have selected", length(isolate(input$"TL_gene")), "genes" ,"//", 
              length(Gene_OtherGlutNs.TL[Gene_OtherGlutNs.TL%in%isolate(input$"TL_gene")]), "detected in the dataset")
      }
      
      else if(!!input$"TL_celltype"=="Other GABANs"){
        paste("You have selected", length(isolate(input$"TL_gene")), "genes" ,"//", 
              length(Gene_OtherGABANs.TL[Gene_OtherGABANs.TL%in%isolate(input$"TL_gene")]), "detected in the dataset")
      }
      
      
    })
    
    ######################
    #LR Significant data   
    ######################

     ########### Table ########

    # render a data table in the output table named output$table_LRDB
    output$table_LRDB <- renderDT(
        LRintercellNetworkDB, rownames=F,
        filter = "top", extensions = c('Buttons'),
        options = list(pageLength = 10,
                       lengthMenu=c(10,20,50,100,500),
                       buttons = list(
                         list(
                           extend = 'collection',
                           buttons = list(
                             list(extend = 'csv', filename = "LR Metadata"),
                             list(extend = 'excel', filename = "LR Metadata")),
                           text = 'Download' )
                       )
          )
        
    )
  
    
    ########### Ligand & Receptor family ##########
    
    # Render a Plotly output showing the ligand or receptor families based on the user's selection
    output$LRfamily <- renderPlotly({
    
    if(!!input$"LRfamily_var"=="ligand_family"){
        plot_LigLRDB_family(LRintercellNetworkDB)
    }
    else if(!!input$"LRfamily_var"=="receptor_family"){
        plot_RecLRDB_family(LRintercellNetworkDB)
    }
        
    })

    
    ############## Table

    # --------------------------------------------------
    # Description: Render a DataTable output showing the significant ligand-receptor interactions for the "E18.5-P0" time point
    # Output: A DataTable output with filtering, pagination, and download options
    # --------------------------------------------------

    output$tableP0_sig <- renderDT({
      
        # Extract the desired columns from the Cortex_edges_sig.list$`E18.5-P0` data frame
        Cortex_edges_sig.list$`E18.5-P0`[,c("ligand","ligand_1","ligand_2","ligand_3","receptor","receptor_1","receptor_2","receptor_3","LR_pair","cluster_L","cluster_R",
                                            "interaction","L_score_S_lr","R_score_S_lr","S_inter","S_inter_diff","pathway","rpathway_score","S_intra","SOURCE.class",
                                            "TARGET.class","Ligand.category","Receptor.category","Ligand.family","Receptor.family","directed","stimulation","inhibition",
                                            "conn.score","dev.score","cell.type.score","detection.rate_L1","detection.rate_L2","detection.rate_L3","detection.rate_R1",
                                            "detection.rate_R2","detection.rate_R3","avg.S_inter","value.S_inter.ctscore","Thr.S_inter","Thr.S_inter_diff","age.at.collection",
                                            "numdup.LR_pair.interaction","numdup.LR_pair.cluster_L","numdup.LR_pair.cluster_R")]},
        
        # Enable top filtering and add the Buttons extension
        filter = "top",extensions = c('Buttons'),rownames=F,

        # Configure DataTable options
        options = list(dom = 'Blfrtip',
                       pageLength=10, lengthMenu=c(10,20,50,100,500),
                       buttons = list(
                         list(
                           extend = 'collection',
                           buttons = list(
                             list(extend = 'csv', filename = "E18.5-P0_LR_table"),
                             list(extend = 'excel', filename = "E18.5-P0_LR_table")),
                           text = 'Download' )
                       )
        )
    )
    

    # --------------------------------------------------
    # Description: Render a DataTable output showing the significant ligand-receptor interactions for the "P1-P2" time point
    # Output: A DataTable output with filtering, pagination, and download options
    # --------------------------------------------------

    output$tableP2_sig <- renderDT({

        # Extract the desired columns from the Cortex_edges_sig.list$`P1-P2` data frame
        Cortex_edges_sig.list$`P1-P2`[,c("ligand","ligand_1","ligand_2","ligand_3","receptor","receptor_1","receptor_2","receptor_3","LR_pair","cluster_L","cluster_R",
                                         "interaction","L_score_S_lr","R_score_S_lr","S_inter","S_inter_diff","pathway","rpathway_score","S_intra","SOURCE.class",
                                         "TARGET.class","Ligand.category","Receptor.category","Ligand.family","Receptor.family","directed","stimulation","inhibition",
                                         "conn.score","dev.score","cell.type.score","detection.rate_L1","detection.rate_L2","detection.rate_L3","detection.rate_R1",
                                         "detection.rate_R2","detection.rate_R3","avg.S_inter","value.S_inter.ctscore","Thr.S_inter","Thr.S_inter_diff","age.at.collection",
                                         "numdup.LR_pair.interaction","numdup.LR_pair.cluster_L","numdup.LR_pair.cluster_R")]},
        
        # Enable top filtering and add the Buttons extension
        filter = "top",extensions = c('Buttons'),rownames=F,

        # Configure DataTable options
        options = list(dom = 'Blfrtip',
                       pageLength=10, lengthMenu=c(10,20,50,100,500),
                       buttons = list(
                         list(
                           extend = 'collection',
                           buttons = list(
                             list(extend = 'csv', filename = "P1-P2_LR_table"),
                             list(extend = 'excel', filename = "P1-P2_LR_table")),
                           text = 'Download' )
                       )
        )
    )
    
    # --------------------------------------------------
    # Description: Render a DataTable output showing the significant ligand-receptor interactions for the "P4-P5" time point
    # Output: A DataTable output with filtering, pagination, and download options
    # --------------------------------------------------
    
    output$tableP5_sig <- renderDT({

        # Extract the desired columns from the Cortex_edges_sig.list$`P4-P5` data frame      
        Cortex_edges_sig.list$`P4-P5`[,c("ligand","ligand_1","ligand_2","ligand_3","receptor","receptor_1","receptor_2","receptor_3","LR_pair","cluster_L","cluster_R",
                                         "interaction","L_score_S_lr","R_score_S_lr","S_inter","S_inter_diff","pathway","rpathway_score","S_intra","SOURCE.class",
                                         "TARGET.class","Ligand.category","Receptor.category","Ligand.family","Receptor.family","directed","stimulation","inhibition",
                                         "conn.score","dev.score","cell.type.score","detection.rate_L1","detection.rate_L2","detection.rate_L3","detection.rate_R1",
                                         "detection.rate_R2","detection.rate_R3","avg.S_inter","value.S_inter.ctscore","Thr.S_inter","Thr.S_inter_diff","age.at.collection",
                                         "numdup.LR_pair.interaction","numdup.LR_pair.cluster_L","numdup.LR_pair.cluster_R")]},
        
        # Enable top filtering and add the Buttons extension
        filter = "top",extensions = c('Buttons'),rownames=F,

        # Configure DataTable options
        options = list(dom = 'Blfrtip',
                       pageLength=10, lengthMenu=c(10,20,50,100,500),
                       buttons = list(
                         list(
                           extend = 'collection',
                           buttons = list(
                             list(extend = 'csv', filename = "P4-P5_LR_table"),
                             list(extend = 'excel', filename = "P4-P5_LR_table")),
                           text = 'Download' )
                       )
        )
    )
    
     # --------------------------------------------------
     # Description: Render a DataTable output showing the significant ligand-receptor interactions for the "P8" time point
     # Output: A DataTable output with filtering, pagination, and download options
     # --------------------------------------------------

    output$tableP8_sig <- renderDT({

      # Extract the desired columns from the Cortex_edges_sig.list$`P8` data frame
      Cortex_edges_sig.list$P8[,c("ligand","ligand_1","ligand_2","ligand_3","receptor","receptor_1","receptor_2","receptor_3","receptor_4","LR_pair","cluster_L",
                                  "cluster_R","interaction","L_score_S_lr","R_score_S_lr","S_inter","S_inter_diff","pathway","rpathway_score","S_intra","SOURCE.class",
                                  "TARGET.class","Ligand.category","Receptor.category","Ligand.family","Receptor.family","directed","stimulation","inhibition",
                                  "conn.score","dev.score","cell.type.score","detection.rate_L1","detection.rate_L2","detection.rate_L3","detection.rate_R1",
                                  "detection.rate_R2","detection.rate_R3","detection.rate_R4","avg.S_inter","value.S_inter.ctscore","Thr.S_inter","Thr.S_inter_diff",
                                  "age.at.collection","numdup.LR_pair.interaction","numdup.LR_pair.cluster_L","numdup.LR_pair.cluster_R")]},

      # Enable top filtering and add the Buttons extension                            
      filter = "top",extensions = c('Buttons'),rownames=F,

      # Configure DataTable options
      options = list(dom = 'Blfrtip',
                     pageLength=10, lengthMenu=c(10,20,50,100,500),
                     buttons = list(
                       list(
                         extend = 'collection',
                         buttons = list(
                           list(extend = 'csv', filename = "P8_LR_table"),
                           list(extend = 'excel', filename = "P8_LR_table")),
                         text = 'Download' )
                     )
      )
    )
    
     # --------------------------------------------------
     # Description: Render a DataTable output showing the significant ligand-receptor interactions for the "P16" time point
     # Output: A DataTable output with filtering, pagination, and download options
     # --------------------------------------------------

    output$tableP16_sig <- renderDT({

      # Extract the desired columns from the Cortex_edges_sig.list$`P16` data frame
      Cortex_edges_sig.list$P16[,c("ligand","ligand_1","ligand_2","ligand_3","receptor","receptor_1","receptor_2","receptor_3","LR_pair","cluster_L","cluster_R",
                                   "interaction","L_score_S_lr","R_score_S_lr","S_inter","S_inter_diff","pathway","rpathway_score","S_intra","SOURCE.class",
                                   "TARGET.class","Ligand.category","Receptor.category","Ligand.family","Receptor.family","directed","stimulation","inhibition",
                                   "conn.score","dev.score","cell.type.score","detection.rate_L1","detection.rate_L2","detection.rate_L3","detection.rate_R1",
                                   "detection.rate_R2","detection.rate_R3","avg.S_inter","value.S_inter.ctscore","Thr.S_inter","Thr.S_inter_diff","age.at.collection",
                                   "numdup.LR_pair.interaction","numdup.LR_pair.cluster_L","numdup.LR_pair.cluster_R")]},

      # Enable top filtering and add the Buttons extension                             
      filter = "top",extensions = c('Buttons'),rownames=F,

      # Configure DataTable options
      options = list(dom = 'Blfrtip',
                     pageLength=10, lengthMenu=c(10,20,50,100,500),
                     buttons = list(
                       list(
                         extend = 'collection',
                         buttons = list(
                           list(extend = 'csv', filename = "P16_LR_table"),
                           list(extend = 'excel', filename = "P16_LR_table")),
                         text = 'Download' )
                     )
      )
    )
    
     # --------------------------------------------------
     # Description: Render a DataTable output showing the significant ligand-receptor interactions for the "P30" time point
     # Output: A DataTable output with filtering, pagination, and download options
     # --------------------------------------------------

    output$tableP30_sig <- renderDT({

        # Extract the desired columns from the Cortex_edges_sig.list$`P30` data frame     
        Cortex_edges_sig.list$P30[,c("ligand","ligand_1","ligand_2","ligand_3","receptor","receptor_1","receptor_2","receptor_3","receptor_4","LR_pair","cluster_L",
                                     "cluster_R","interaction","L_score_S_lr","R_score_S_lr","S_inter","S_inter_diff","pathway","rpathway_score","S_intra","SOURCE.class",
                                     "TARGET.class","Ligand.category","Receptor.category","Ligand.family","Receptor.family","directed","stimulation","inhibition",
                                     "conn.score","dev.score","cell.type.score","detection.rate_L1","detection.rate_L2","detection.rate_L3","detection.rate_R1",
                                     "detection.rate_R2","detection.rate_R3","detection.rate_R4","avg.S_inter","value.S_inter.ctscore","Thr.S_inter","Thr.S_inter_diff",
                                     "age.at.collection","numdup.LR_pair.interaction","numdup.LR_pair.cluster_L","numdup.LR_pair.cluster_R")]},
        
        # Enable top filtering and add the Buttons extension
        filter = "top",extensions = c('Buttons'),rownames=F,

        # Configure DataTable options
        options = list(dom = 'Blfrtip',
                       pageLength=10, lengthMenu=c(10,20,50,100,500),
                       buttons = list(
                         list(
                           extend = 'collection',
                           buttons = list(
                             list(extend = 'csv', filename = "P30_LR_table"),
                             list(extend = 'excel', filename = "P30_LR_table")),
                           text = 'Download' )
                       )
        )
    )
      # --------------------------------------------------
      # Description: Render a DataTable output showing the significant ligand-receptor interactions for the "Adult" time point
      # Output: A DataTable output with filtering, pagination, and download options
      # --------------------------------------------------
    
    output$tableAdult_sig <- renderDT({

        # Extract the desired columns from the Cortex_edges_sig.list$`Adult` data frame
        Cortex_edges_sig.list$Adult[,c("ligand","ligand_1","ligand_2","ligand_3","receptor","receptor_1","receptor_2","receptor_3","receptor_4","LR_pair","cluster_L",
                                       "cluster_R","interaction","L_score_S_lr","R_score_S_lr","S_inter","S_inter_diff","pathway","rpathway_score","S_intra","SOURCE.class",
                                       "TARGET.class","Ligand.category","Receptor.category","Ligand.family","Receptor.family","directed","stimulation","inhibition",
                                       "conn.score","dev.score","cell.type.score","detection.rate_L1","detection.rate_L2","detection.rate_L3","detection.rate_R1",
                                       "detection.rate_R2","detection.rate_R3","detection.rate_R4","avg.S_inter","value.S_inter.ctscore","Thr.S_inter","Thr.S_inter_diff",
                                       "age.at.collection","numdup.LR_pair.interaction","numdup.LR_pair.cluster_L","numdup.LR_pair.cluster_R")]},

        # Enable top filtering and add the Buttons extension
        filter = "top",extensions = c('Buttons'), rownames=F,

        # Configure DataTable options
        options = list(dom = 'Blfrtip',
                       pageLength=10, lengthMenu=c(10,20,50,100,500),
                       buttons = list(
                         list(
                           extend = 'collection',
                           buttons = list(
                             list(extend = 'csv', filename = "Adult_LR_table"),
                             list(extend = 'excel', filename = "Adult_LR_table")),
                           text = 'Download' )
                       )
        )
    )
    
    
    ################ Heatmap: Number of interactions

    # --------------------------------------------------
    # Description: Create legends and update picker inputs for various filters in the heatmap of significant ligand-receptor interactions
    # --------------------------------------------------
    
    # Create legends for class and normalized soma depth
    lgdCT = Legend(labels = c("GABAergic", "Glutamatergic"),legend_gp = gpar(fill=c("#D80D00","#1144E3")) ,title = "Class",
                   title_gp = gpar(fontsize = 16, fontface = "bold"), labels_gp = gpar(fontsize = 14.5))
    lgdSoma = Legend(title="Normalized \nsoma depth", at = rev(c(0, 0.2, 0.4, 0.6,0.8,1)), 
                     col_fun = colorRamp2(c(0,1),c("grey95","black")), title_gp = gpar(fontsize = 16, fontface = "bold"), labels_gp = gpar(fontsize = 14.5))
    
    # Update picker inputs for various filters
    observe({
        updatePickerInput(
            session = session,
            inputId = "htNinter_Sending.cluster",
            choices = sort(unique(Cortex_edges_sig.list$Cortex$cluster_L)),
            selected = sort(unique(Cortex_edges_sig.list$Cortex$cluster_L))
        )
    })
    
    observe({
        updatePickerInput(
            session = session,
            inputId = "htNinter_Target.cluster",
            choices = sort(unique(Cortex_edges_sig.list$Cortex$cluster_R)),
            selected = sort(unique(Cortex_edges_sig.list$Cortex$cluster_R))
        )
    })
    
    observe({
        updatePickerInput(
            session = session,
            inputId = "htNinter_SOURCE.class",
            choices = sort(unique(Cortex_edges_sig.list$Cortex$SOURCE.class)),
            selected = sort(unique(Cortex_edges_sig.list$Cortex$SOURCE.class))
        )
    })
    
    observe({
        updatePickerInput(
            session = session,
            inputId = "htNinter_TARGET.class",
            choices = sort(unique(Cortex_edges_sig.list$Cortex$TARGET.class)),
            selected = sort(unique(Cortex_edges_sig.list$Cortex$TARGET.class))
        )
    })

    observe({
        updatePickerInput(
            session = session,
            inputId = "htNinter_LR",
            choices = sort(unique(c(LRintercellNetworkDB$ligand_receptor, Cortex_edges_sig.list$Cortex$LR_pair))),
            selected = sort(unique(c(LRintercellNetworkDB$ligand_receptor, Cortex_edges_sig.list$Cortex$LR_pair)))
        )
    })
    

    observe({
        updatePickerInput(
            session = session,
            inputId = "Cortex_htNinter_Ligand.category",
            choices = sort(unique(unlist(strsplit(LRintercellNetworkDB$ligand_category, fixed=T, split=";")))),
            selected = sort(unique(unlist(strsplit(LRintercellNetworkDB$ligand_category, fixed=T, split=";"))))
        )
    })
    
    observe({
        updatePickerInput(
            session = session,
            inputId = "Cortex_htNinter_Receptor.category",
            choices = sort(unique(unlist(strsplit(LRintercellNetworkDB$receptor_category, fixed=T, split=";")))),
            selected = sort(unique(unlist(strsplit(LRintercellNetworkDB$receptor_category, fixed=T, split=";"))))
        )
    })
    
    observe({
        updatePickerInput(
            session = session,
            inputId = "Cortex_htNinter_Ligand.family",
            choices = sort(unique(unlist(strsplit(LRintercellNetworkDB$ligand_family, fixed=T, split=";")))),
            selected = sort(unique(unlist(strsplit(LRintercellNetworkDB$ligand_family, fixed=T, split=";"))))
        )
    })
    
    observe({
        updatePickerInput(
            session = session,
            inputId = "Cortex_htNinter_Receptor.family",
            choices = sort(unique(unlist(strsplit(LRintercellNetworkDB$receptor_family, fixed=T, split=";")))),
            selected = sort(unique(unlist(strsplit(LRintercellNetworkDB$receptor_family, fixed=T, split=";"))))
        )
    })
    
    observe({
    updatePickerInput(
     session = session,
    inputId = "Cortex_htNinter_pathway",
    choices = sort(unique(Cortex_edges_sig.list$Cortex$pathway)),
    selected = sort(unique(Cortex_edges_sig.list$Cortex$pathway))
    )
    })
    
    observe({
      updatePickerInput(
        session = session,
        inputId = "htNinter_numdup.LR_pair.interaction",
        choices = str_sort(unique(Cortex_edges_sig.list$Cortex$numdup.LR_pair.interaction), numeric=T),
        selected = str_sort(unique(Cortex_edges_sig.list$Cortex$numdup.LR_pair.interaction), numeric=T)
      )
    })
    

    
    
    observe({
      updatePickerInput(
        session = session,
        inputId = "htNinter_numdup.LR_pair.cluster_L",
        choices = str_sort(unique(Cortex_edges_sig.list$Cortex$numdup.LR_pair.cluster_L), numeric=T),
        selected = str_sort(unique(Cortex_edges_sig.list$Cortex$numdup.LR_pair.cluster_L), numeric=T)
      )
    })
    
    observe({
      updatePickerInput(
        session = session,
        inputId = "htNinter_numdup.LR_pair.cluster_R",
        choices = str_sort(unique(Cortex_edges_sig.list$Cortex$numdup.LR_pair.cluster_R), numeric=T),
        selected = str_sort(unique(Cortex_edges_sig.list$Cortex$numdup.LR_pair.cluster_R), numeric=T)
      )
    })
 
    ## P0_wt
    
    # --------------------------------------------------
    # Description: Set key columns for the Cortex_edges_sig.list$`E18.5-P0` data table for efficient merging and subsetting later
    # --------------------------------------------------
   
    setkey(Cortex_edges_sig.list$`E18.5-P0`, cluster_L, cluster_R, SOURCE.class, TARGET.class, LR_pair, numdup.LR_pair.interaction, numdup.LR_pair.cluster_L,
           numdup.LR_pair.cluster_R, pathway, S_inter, S_intra, S_inter_diff,
           Ligand.category1, Ligand.category2, Ligand.category3, Ligand.category4, Ligand.category5, Ligand.category6, Ligand.category7,
           Receptor.category1, Receptor.category2, Receptor.category3, Receptor.category4, Receptor.category5,
           Ligand.family1, Ligand.family2, Ligand.family3, Ligand.family4, Ligand.family5,
           Receptor.family1, Receptor.family2, Receptor.family3, Receptor.family4, Receptor.family5)
   


    # --------------------------------------------------
    #  Description: eventReactive function filters the data table "Cortex_edges_sig.list$'E18.5-P0'" based on user-defined input criteria 
    #  related to ligand and receptor categories and families.
    #  Output : the filtered data as a data frame.
    # --------------------------------------------------

    P0_htNinter_filter <- eventReactive(input$"htNinter_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$`E18.5-P0`[Ligand.category1 %in% req(input$"Cortex_htNinter_Ligand.category")  |
                        Ligand.category2  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category3  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category4  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category5  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category6  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category7  %in% req(input$"Cortex_htNinter_Ligand.category")
      ]
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_htNinter_Receptor.category") |
                        Receptor.category2  %in% req(input$"Cortex_htNinter_Receptor.category") |
                        Receptor.category3  %in% req(input$"Cortex_htNinter_Receptor.category") |
                        Receptor.category4  %in% req(input$"Cortex_htNinter_Receptor.category") |
                        Receptor.category5  %in% req(input$"Cortex_htNinter_Receptor.category") 
      ]
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family3  %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_htNinter_Ligand.family")
      ]
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_htNinter_Receptor.family") |
                      Receptor.family2 %in% req(input$"Cortex_htNinter_Receptor.family") |
                      Receptor.family3 %in% req(input$"Cortex_htNinter_Receptor.family") |
                      Receptor.family4  %in% req(input$"Cortex_htNinter_Receptor.family") |
                      Receptor.family5  %in% req(input$"Cortex_htNinter_Receptor.family")  
      ]
      mydata <- mydata[S_intra >= req(input$"htNinter_S_intra") |
                         is.na(S_intra)==T]
      mydata <- mydata[cluster_L %in% req(input$"htNinter_Sending.cluster") & 
                      cluster_R %in% req(input$"htNinter_Target.cluster") & 
                      SOURCE.class %in% req(input$"htNinter_SOURCE.class") & 
                      TARGET.class %in% req(input$"htNinter_TARGET.class") & 
                      LR_pair %in% req(input$"htNinter_LR") & 
                      numdup.LR_pair.interaction %in% req(input$"htNinter_numdup.LR_pair.interaction") & 
                      numdup.LR_pair.cluster_L %in% req(input$"htNinter_numdup.LR_pair.cluster_L") & 
                      numdup.LR_pair.cluster_R  %in% req(input$"htNinter_numdup.LR_pair.cluster_R") & 
                      pathway  %in% req(input$"Cortex_htNinter_pathway") & 
                      S_inter >= req(input$"htNinter_S_inter") &
                      S_inter_diff >= req(input$"htNinter_S_inter_diff")]

   

      # Convert the data.table back to a data.frame 
      as.data.frame(mydata)
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render an interactive heatmap plot using plotly for the filtered ligand-receptor interaction data
    # Input:
    #   - P0_htNinter_filter(): A reactive function that returns the filtered data frame based on user input criteria
    #   - Cortex_edges_sig.list$`E18.5-P0`: A data frame containing ligand-receptor interaction data for the E18.5-P0 time point
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------

        output$P0_htNinter <- renderPlotly({

          # Check if the filtered data frame is empty
          if (nrow(req(P0_htNinter_filter())) == 0) return(NULL)

          # Generate an interactive heatmap plot using the scSeqComm_heatmaply_cardinality function
          scSeqComm_heatmaply_cardinality(as.data.frame(Cortex_edges_sig.list$`E18.5-P0`),req(P0_htNinter_filter()), y = "cluster_L", x = "cluster_R", hover_anno="Number of interactions",
                                                      ylab = "Ligand-expressing cell types", xlab = "Receptor-expressing cell types",  
                                                      cl.order=order_celltype[order_celltype$celltype_label%in%unique(Cortex_edges_sig.list$`E18.5-P0`$cluster_L),],
                                                      key.title="Number of \ninteractions", hide_colorbar=F, margins=c(0,0,0,0), sizefont=10, 
                                                      main="", scale="none", 
                                                      Sending.cl=req(input$"htNinter_Sending.cluster"),
                                                      Target.cl=req(input$"htNinter_Target.cluster"), SOURCE.CT=req(input$"htNinter_SOURCE.class"),TARGET.CT=req(input$"htNinter_TARGET.class"),
                                                      col.heatmap=req(input$"ht_ninter_pal"), height=req(input$"htNinter_height"))
          
        

    })

    # --------------------------------------------------
    # Description: Create an eventReactive function that filters the LRintercellNetworkDB data frame based on user input criteria
    # Output: A filtered data frame containing rows that match the user input criteria
    # --------------------------------------------------

    LRDB_mouse_Hm <- eventReactive(input$"htNinter_button",{

      # Filter the LRintercellNetworkDB data frame based on user input criteria
      myLRdata <-  LRintercellNetworkDB[which( LRintercellNetworkDB$ligand_receptor %in% req(input$"htNinter_LR")   &
                                                   rownames(LRintercellNetworkDB) %in% rownames(LRintercellNetworkDB[unlist(sapply(input$"Cortex_htNinter_Ligand.category", grep, LRintercellNetworkDB$ligand_category, fixed=T, USE.NAMES = F)),]) &
                                                   rownames(LRintercellNetworkDB) %in% rownames(LRintercellNetworkDB[unlist(sapply(input$"Cortex_htNinter_Receptor.category", grep, LRintercellNetworkDB$receptor_category, fixed=T, USE.NAMES = F)),]) &
                                                   rownames(LRintercellNetworkDB) %in% rownames(LRintercellNetworkDB[unlist(sapply(input$"Cortex_htNinter_Ligand.family", grep, LRintercellNetworkDB$ligand_family, fixed=T, USE.NAMES = F)),]) &
                                                   rownames(LRintercellNetworkDB) %in% rownames(LRintercellNetworkDB[unlist(sapply(input$"Cortex_htNinter_Receptor.family", grep, LRintercellNetworkDB$receptor_family, fixed=T, USE.NAMES = F)),])),]
      
      # Return the filtered data frame
      myLRdata

    },ignoreNULL = FALSE)
        
    # --------------------------------------------------
    # Description: Create a text output that displays the number of selected ligand-receptor pairs and the number of pairs expressed in the E18.5-P0 dataset with the current filters
    # --------------------------------------------------

    output$P0_htNinter_selectLR <- renderText({
        

        paste("Your selection contains", nrow(LRDB_mouse_Hm()), "ligand-receptor pairs", "|", 
              length(unique(P0_htNinter_filter()$LR_pair[which(P0_htNinter_filter()$cluster_L %in% req(input$"htNinter_Sending.cluster")  & 
                                                                       P0_htNinter_filter()$cluster_R %in% req(input$"htNinter_Target.cluster")  &  
                                                                       P0_htNinter_filter()$SOURCE.class %in% req(input$"htNinter_SOURCE.class")  &  
                                                                       P0_htNinter_filter()$TARGET.class %in% req(input$"htNinter_TARGET.class") )])), 
              "expressed in the E18.5-P0 dataset with the current filters")
        
        
        
        
    })
    
    # Create a plot output for the combined legend of the heatmap of ligand-receptor interactions
    output$P0_lgd_ninter_ht <- renderPlot({
        draw(packLegend(lgdCT, lgdSoma))
    })
    
    ## P2
    
    # --------------------------------------------------
    # Description: Create an eventReactive function that filters the Cortex_edges_sig.list$P1-P2 data frame based on user input criteria
    # Output: A filtered data frame containing rows that match the user input criteria
    # --------------------------------------------------

    setkey(Cortex_edges_sig.list$`P1-P2`, cluster_L, cluster_R, SOURCE.class, TARGET.class, LR_pair, numdup.LR_pair.interaction, numdup.LR_pair.cluster_L,
           numdup.LR_pair.cluster_R, pathway, S_inter, S_intra, S_inter_diff,
           Ligand.category1, Ligand.category2, Ligand.category3, Ligand.category4, Ligand.category5, Ligand.category6, Ligand.category7,
           Receptor.category1, Receptor.category2, Receptor.category3, Receptor.category4, Receptor.category5,
           Ligand.family1, Ligand.family2, Ligand.family3, Ligand.family4, Ligand.family5,
           Receptor.family1, Receptor.family2, Receptor.family3, Receptor.family4, Receptor.family5)
    
    # --------------------------------------------------
    # Create an eventReactive function that filters the Cortex_edges_sig.list$P1-P2 data frame based on user input criteria
    # Output: A filtered data frame containing rows that match the user input criteria
    # --------------------------------------------------

    P2_htNinter_filter <- eventReactive(input$"htNinter_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$`P1-P2`[Ligand.category1 %in% req(input$"Cortex_htNinter_Ligand.category")  |
                        Ligand.category2  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category3  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category4  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category5  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category6  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category7  %in% req(input$"Cortex_htNinter_Ligand.category")
      ]
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_htNinter_Receptor.category") 
      ]
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family3  %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_htNinter_Ligand.family")
      ]
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family4  %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family5  %in% req(input$"Cortex_htNinter_Receptor.family")  
      ]
      mydata <- mydata[S_intra >= req(input$"htNinter_S_intra") |
                         is.na(S_intra)==T]
      mydata <- mydata[cluster_L %in% req(input$"htNinter_Sending.cluster") & 
                         cluster_R %in% req(input$"htNinter_Target.cluster") & 
                         SOURCE.class %in% req(input$"htNinter_SOURCE.class") & 
                         TARGET.class %in% req(input$"htNinter_TARGET.class") & 
                         LR_pair %in% req(input$"htNinter_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"htNinter_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"htNinter_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"htNinter_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_htNinter_pathway") & 
                         S_inter >= req(input$"htNinter_S_inter") &
                         S_inter_diff >= req(input$"htNinter_S_inter_diff")]
      
      
      
      # Convert the data.table back to a data.frame 
      as.data.frame(mydata)
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render an interactive heatmap plot using plotly for the filtered ligand-receptor interaction data
    # Input:
    #   - P2_htNinter_filter(): A reactive function that returns the filtered data frame based on user input criteria
    #   - Cortex_edges_sig.list$`P1-P2`: A data frame containing ligand-receptor interaction data for the P1-P2 time point
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------

    output$P2_htNinter <- renderPlotly({
      
      # Check if the filtered data frame is empty
      if (nrow(req(P2_htNinter_filter())) == 0) return(NULL)
      
      # Generate an interactive heatmap plot using the scSeqComm_heatmaply_cardinality function
      scSeqComm_heatmaply_cardinality(as.data.frame(Cortex_edges_sig.list$`P1-P2`),req(P2_htNinter_filter()), y = "cluster_L", x = "cluster_R", hover_anno="Number of interactions",
                                      ylab = "Ligand-expressing cell types", xlab = "Receptor-expressing cell types",  
                                      cl.order=order_celltype[order_celltype$celltype_label%in%unique(Cortex_edges_sig.list$`P1-P2`$cluster_L),],
                                      key.title="Number of \ninteractions", hide_colorbar=F, margins=c(0,0,0,0), sizefont=10, 
                                      main="", scale="none", 
                                      Sending.cl=req(input$"htNinter_Sending.cluster"),
                                      Target.cl=req(input$"htNinter_Target.cluster"), SOURCE.CT=req(input$"htNinter_SOURCE.class"),TARGET.CT=req(input$"htNinter_TARGET.class"),
                                      col.heatmap=req(input$"ht_ninter_pal"), height=req(input$"htNinter_height"))
      
    })
    

    # --------------------------------------------------
    # Description: Create a text output that displays the number of selected ligand-receptor pairs and the number of pairs expressed in the P1-P2 dataset with the current filters
    # --------------------------------------------------

    output$P2_htNinter_selectLR <- renderText({
      
      
      paste("Your selection contains", nrow(LRDB_mouse_Hm()), "ligand-receptor pairs", "|", 
            length(unique(P2_htNinter_filter()$LR_pair[which(P2_htNinter_filter()$cluster_L %in% req(input$"htNinter_Sending.cluster")  & 
                                                               P2_htNinter_filter()$cluster_R %in% req(input$"htNinter_Target.cluster")  &  
                                                               P2_htNinter_filter()$SOURCE.class %in% req(input$"htNinter_SOURCE.class")  &  
                                                               P2_htNinter_filter()$TARGET.class %in% req(input$"htNinter_TARGET.class") )])), 
            "expressed in the P1-P2 dataset with the current filters")
      
      
      
      
    })

    # Create a plot output for the combined legend of the heatmap of ligand-receptor interactions
    output$P2_lgd_ninter_ht <- renderPlot({
      draw(packLegend(lgdCT, lgdSoma))
    })

    #P5

    # --------------------------------------------------
    # Description: Create an eventReactive function that filters the Cortex_edges_sig.list$P4-P5 data frame based on user input criteria
    # Output: A filtered data frame containing rows that match the user input criteria
    # --------------------------------------------------

    setkey(Cortex_edges_sig.list$`P4-P5`, cluster_L, cluster_R, SOURCE.class, TARGET.class, LR_pair, numdup.LR_pair.interaction, numdup.LR_pair.cluster_L,
           numdup.LR_pair.cluster_R, pathway, S_inter, S_intra, S_inter_diff,
           Ligand.category1, Ligand.category2, Ligand.category3, Ligand.category4, Ligand.category5, Ligand.category6, Ligand.category7,
           Receptor.category1, Receptor.category2, Receptor.category3, Receptor.category4, Receptor.category5,
           Ligand.family1, Ligand.family2, Ligand.family3, Ligand.family4, Ligand.family5,
           Receptor.family1, Receptor.family2, Receptor.family3, Receptor.family4, Receptor.family5)
    
    # --------------------------------------------------
    # Create an eventReactive function that filters the Cortex_edges_sig.list$P4-P5 data frame based on user input criteria
    # Output: A filtered data frame containing rows that match the user input criteria
    # --------------------------------------------------
    
    P5_htNinter_filter <- eventReactive(input$"htNinter_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$`P4-P5`[Ligand.category1 %in% req(input$"Cortex_htNinter_Ligand.category")  |
                        Ligand.category2  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category3  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category4  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category5  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category6  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category7  %in% req(input$"Cortex_htNinter_Ligand.category")
      ]
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_htNinter_Receptor.category") 
      ]
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family3  %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_htNinter_Ligand.family")
      ]
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family4  %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family5  %in% req(input$"Cortex_htNinter_Receptor.family")  
      ]
      mydata <- mydata[S_intra >= req(input$"htNinter_S_intra") |
                         is.na(S_intra)==T]
      mydata <- mydata[cluster_L %in% req(input$"htNinter_Sending.cluster") & 
                         cluster_R %in% req(input$"htNinter_Target.cluster") & 
                         SOURCE.class %in% req(input$"htNinter_SOURCE.class") & 
                         TARGET.class %in% req(input$"htNinter_TARGET.class") & 
                         LR_pair %in% req(input$"htNinter_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"htNinter_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"htNinter_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"htNinter_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_htNinter_pathway") & 
                         S_inter >= req(input$"htNinter_S_inter") &
                         S_inter_diff >= req(input$"htNinter_S_inter_diff")]
      
      
      
      # Convert the data.table back to a data.frame 
      as.data.frame(mydata)
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render an interactive heatmap plot using plotly for the filtered ligand-receptor interaction data
    # Input:
    #   - P5_htNinter_filter(): A reactive function that returns the filtered data frame based on user input criteria
    #   - Cortex_edges_sig.list$`P4-P5`: A data frame containing ligand-receptor interaction data for the "P4-P5" time point
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------
    
    output$P5_htNinter <- renderPlotly({
      
      # Check if the filtered data frame is empty
      if (nrow(req(P5_htNinter_filter())) == 0) return(NULL)
      
      # Generate an interactive heatmap plot using the scSeqComm_heatmaply_cardinality function
      scSeqComm_heatmaply_cardinality(as.data.frame(Cortex_edges_sig.list$`P4-P5`),req(P5_htNinter_filter()), y = "cluster_L", x = "cluster_R", hover_anno="Number of interactions",
                                      ylab = "Ligand-expressing cell types", xlab = "Receptor-expressing cell types",  
                                      cl.order=order_celltype[order_celltype$celltype_label%in%unique(Cortex_edges_sig.list$`P4-P5`$cluster_L),],
                                      key.title="Number of \ninteractions", hide_colorbar=F, margins=c(0,0,0,0), sizefont=10, 
                                      main="", scale="none", 
                                      Sending.cl=req(input$"htNinter_Sending.cluster"),
                                      Target.cl=req(input$"htNinter_Target.cluster"), SOURCE.CT=req(input$"htNinter_SOURCE.class"),TARGET.CT=req(input$"htNinter_TARGET.class"),
                                      col.heatmap=req(input$"ht_ninter_pal"), height=req(input$"htNinter_height"))
      
    })
    
    # --------------------------------------------------
    # Description: Create a text output that displays the number of selected ligand-receptor pairs and the number of pairs expressed in the P4-P5 dataset with the current filters
    # --------------------------------------------------

    output$P5_htNinter_selectLR <- renderText({
      
      
      paste("Your selection contains", nrow(LRDB_mouse_Hm()), "ligand-receptor pairs", "|", 
            length(unique(P5_htNinter_filter()$LR_pair[which(P5_htNinter_filter()$cluster_L %in% req(input$"htNinter_Sending.cluster")  & 
                                                               P5_htNinter_filter()$cluster_R %in% req(input$"htNinter_Target.cluster")  &  
                                                               P5_htNinter_filter()$SOURCE.class %in% req(input$"htNinter_SOURCE.class")  &  
                                                               P5_htNinter_filter()$TARGET.class %in% req(input$"htNinter_TARGET.class") )])), 
            "expressed in the P4-P5 dataset with the current filters")
      
      
      
      
    })
    
    # Create a plot output for the combined legend of the heatmap of ligand-receptor interactions
    output$P5_lgd_ninter_ht <- renderPlot({
      draw(packLegend(lgdCT, lgdSoma))
    })
    
    
    #P8
    
     # --------------------------------------------------
     # Description: Create an eventReactive function that filters the Cortex_edges_sig.list$P8 data frame based on user input criteria
     # Output: A filtered data frame containing rows that match the user input criteria
     # --------------------------------------------------

    setkey(Cortex_edges_sig.list$P8, cluster_L, cluster_R, SOURCE.class, TARGET.class, LR_pair, numdup.LR_pair.interaction, numdup.LR_pair.cluster_L,
           numdup.LR_pair.cluster_R, pathway, S_inter, S_intra, S_inter_diff,
           Ligand.category1, Ligand.category2, Ligand.category3, Ligand.category4, Ligand.category5, Ligand.category6, 
           Receptor.category1, Receptor.category2, Receptor.category3, Receptor.category4, Receptor.category5,
           Ligand.family1, Ligand.family2, Ligand.family3, Ligand.family4, Ligand.family5,
           Receptor.family1, Receptor.family2, Receptor.family3, Receptor.family4, Receptor.family5)
    
    # --------------------------------------------------
    # Create an eventReactive function that filters the Cortex_edges_sig.list$P8 data frame based on user input criteria
    # Output: A filtered data frame containing rows that match the user input criteria
    # --------------------------------------------------
    
    P8_htNinter_filter <- eventReactive(input$"htNinter_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$P8[Ligand.category1 %in% req(input$"Cortex_htNinter_Ligand.category")  |
                        Ligand.category2  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category3  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category4  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category5  %in% req(input$"Cortex_htNinter_Ligand.category") |
                        Ligand.category6  %in% req(input$"Cortex_htNinter_Ligand.category") 
      ]
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_htNinter_Receptor.category") 
      ]
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family3  %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_htNinter_Ligand.family")
      ]
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family4  %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family5  %in% req(input$"Cortex_htNinter_Receptor.family")  
      ]
      mydata <- mydata[S_intra >= req(input$"htNinter_S_intra") |
                         is.na(S_intra)==T]
      mydata <- mydata[cluster_L %in% req(input$"htNinter_Sending.cluster") & 
                         cluster_R %in% req(input$"htNinter_Target.cluster") & 
                         SOURCE.class %in% req(input$"htNinter_SOURCE.class") & 
                         TARGET.class %in% req(input$"htNinter_TARGET.class") & 
                         LR_pair %in% req(input$"htNinter_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"htNinter_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"htNinter_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"htNinter_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_htNinter_pathway") & 
                         S_inter >= req(input$"htNinter_S_inter") &
                         S_inter_diff >= req(input$"htNinter_S_inter_diff")]
      
      
      
      # Convert the data.table back to a data.frame 
      as.data.frame(mydata)
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render an interactive heatmap plot using plotly for the filtered ligand-receptor interaction data
    # Input:
    #   - P8_htNinter_filter(): A reactive function that returns the filtered data frame based on user input criteria
    #   - Cortex_edges_sig.list$`P8`: A data frame containing ligand-receptor interaction data for the "P8" time point
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------

    output$P8_htNinter <- renderPlotly({
      
      # Check if the filtered data frame is empty
      if (nrow(req(P8_htNinter_filter())) == 0) return(NULL)
      
      # Generate an interactive heatmap plot using the scSeqComm_heatmaply_cardinality function
      scSeqComm_heatmaply_cardinality(as.data.frame(Cortex_edges_sig.list$P8),req(P8_htNinter_filter()), y = "cluster_L", x = "cluster_R", hover_anno="Number of interactions",
                                      ylab = "Ligand-expressing cell types", xlab = "Receptor-expressing cell types",  
                                      cl.order=order_celltype[order_celltype$celltype_label%in%unique(Cortex_edges_sig.list$P8$cluster_L),],
                                      key.title="Number of \ninteractions", hide_colorbar=F, margins=c(0,0,0,0), sizefont=10, 
                                      main="", scale="none", 
                                      Sending.cl=req(input$"htNinter_Sending.cluster"),
                                      Target.cl=req(input$"htNinter_Target.cluster"), SOURCE.CT=req(input$"htNinter_SOURCE.class"),TARGET.CT=req(input$"htNinter_TARGET.class"),
                                      col.heatmap=req(input$"ht_ninter_pal"), height=req(input$"htNinter_height"))
      
    })
    
     # --------------------------------------------------
     # Description: Create a text output that displays the number of selected ligand-receptor pairs and the number of pairs expressed in the P8 dataset with the current filters
     # --------------------------------------------------
    
    output$P8_htNinter_selectLR <- renderText({
      
      
      paste("Your selection contains", nrow(LRDB_mouse_Hm()), "ligand-receptor pairs", "|", 
            length(unique(P8_htNinter_filter()$LR_pair[which(P8_htNinter_filter()$cluster_L %in% req(input$"htNinter_Sending.cluster")  & 
                                                               P8_htNinter_filter()$cluster_R %in% req(input$"htNinter_Target.cluster")  &  
                                                               P8_htNinter_filter()$SOURCE.class %in% req(input$"htNinter_SOURCE.class")  &  
                                                               P8_htNinter_filter()$TARGET.class %in% req(input$"htNinter_TARGET.class") )])), 
            "expressed in the P8 dataset with the current filters")
      
      
      
      
    })
    
    # Create a plot output for the combined legend of the heatmap of ligand-receptor interactions
    output$P8_lgd_ninter_ht <- renderPlot({
      draw(packLegend(lgdCT, lgdSoma))
    })
    
    
    #P16

      # --------------------------------------------------
      # Description: Create an eventReactive function that filters the Cortex_edges_sig.list$P16 data frame based on user input criteria
      # Output: A filtered data frame containing rows that match the user input criteria
      # --------------------------------------------------
    
    setkey(Cortex_edges_sig.list$P16, cluster_L, cluster_R, SOURCE.class, TARGET.class, LR_pair, numdup.LR_pair.interaction, numdup.LR_pair.cluster_L,
           numdup.LR_pair.cluster_R, pathway, S_inter, S_intra, S_inter_diff,
           Ligand.category1, Ligand.category2, Ligand.category3, Ligand.category4, Ligand.category5, Ligand.category6, 
           Receptor.category1, Receptor.category2, Receptor.category3, Receptor.category4, Receptor.category5,
           Ligand.family1, Ligand.family2, Ligand.family3, Ligand.family4, Ligand.family5,
           Receptor.family1, Receptor.family2, Receptor.family3, Receptor.family4, Receptor.family5)
    
      # --------------------------------------------------
      # Create an eventReactive function that filters the Cortex_edges_sig.list$P16 data frame based on user input criteria
      # Output: A filtered data frame containing rows that match the user input criteria
      # --------------------------------------------------
    
    P16_htNinter_filter <- eventReactive(input$"htNinter_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$P16[Ligand.category1 %in% req(input$"Cortex_htNinter_Ligand.category")  |
                         Ligand.category2  %in% req(input$"Cortex_htNinter_Ligand.category") |
                         Ligand.category3  %in% req(input$"Cortex_htNinter_Ligand.category") |
                         Ligand.category4  %in% req(input$"Cortex_htNinter_Ligand.category") |
                         Ligand.category5  %in% req(input$"Cortex_htNinter_Ligand.category") |
                         Ligand.category6  %in% req(input$"Cortex_htNinter_Ligand.category")
      ]
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_htNinter_Receptor.category") 
      ]
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family3  %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_htNinter_Ligand.family")
      ]
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family4  %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family5  %in% req(input$"Cortex_htNinter_Receptor.family")  
      ]
      mydata <- mydata[S_intra >= req(input$"htNinter_S_intra") |
                         is.na(S_intra)==T]
      mydata <- mydata[cluster_L %in% req(input$"htNinter_Sending.cluster") & 
                         cluster_R %in% req(input$"htNinter_Target.cluster") & 
                         SOURCE.class %in% req(input$"htNinter_SOURCE.class") & 
                         TARGET.class %in% req(input$"htNinter_TARGET.class") & 
                         LR_pair %in% req(input$"htNinter_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"htNinter_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"htNinter_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"htNinter_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_htNinter_pathway") & 
                         S_inter >= req(input$"htNinter_S_inter") &
                         S_inter_diff >= req(input$"htNinter_S_inter_diff")]
      
      # Convert the data.table back to a data.frame 
      as.data.frame(mydata)
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render an interactive heatmap plot using plotly for the filtered ligand-receptor interaction data
    # Input:
    #   - P16_htNinter_filter(): A reactive function that returns the filtered data frame based on user input criteria
    #   - Cortex_edges_sig.list$`P8`: A data frame containing ligand-receptor interaction data for the "P16" time point
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------
    
    output$P16_htNinter <- renderPlotly({
      
      # Check if the filtered data frame is empty
      if (nrow(req(P16_htNinter_filter())) == 0) return(NULL)
      
      # Generate an interactive heatmap plot using the scSeqComm_heatmaply_cardinality function
      scSeqComm_heatmaply_cardinality(as.data.frame(Cortex_edges_sig.list$P16),req(P16_htNinter_filter()), y = "cluster_L", x = "cluster_R", hover_anno="Number of interactions",
                                      ylab = "Ligand-expressing cell types", xlab = "Receptor-expressing cell types",  
                                      cl.order=order_celltype[order_celltype$celltype_label%in%unique(Cortex_edges_sig.list$P16$cluster_L),],
                                      key.title="Number of \ninteractions", hide_colorbar=F, margins=c(0,0,0,0), sizefont=10, 
                                      main="", scale="none", 
                                      Sending.cl=req(input$"htNinter_Sending.cluster"),
                                      Target.cl=req(input$"htNinter_Target.cluster"), SOURCE.CT=req(input$"htNinter_SOURCE.class"),TARGET.CT=req(input$"htNinter_TARGET.class"),
                                      col.heatmap=req(input$"ht_ninter_pal"), height=req(input$"htNinter_height"))
      
    })
    
      # --------------------------------------------------
      # Description: Create a text output that displays the number of selected ligand-receptor pairs and the number of pairs expressed in the P16 dataset with the current filters
      # --------------------------------------------------
    
    output$P16_htNinter_selectLR <- renderText({
      
      
      paste("Your selection contains", nrow(LRDB_mouse_Hm()), "ligand-receptor pairs", "|", 
            length(unique(P16_htNinter_filter()$LR_pair[which(P16_htNinter_filter()$cluster_L %in% req(input$"htNinter_Sending.cluster")  & 
                                                               P16_htNinter_filter()$cluster_R %in% req(input$"htNinter_Target.cluster")  &  
                                                               P16_htNinter_filter()$SOURCE.class %in% req(input$"htNinter_SOURCE.class")  &  
                                                               P16_htNinter_filter()$TARGET.class %in% req(input$"htNinter_TARGET.class") )])), 
            "expressed in the P16 dataset with the current filters")
      
      
      
      
    })
    
    # Create a plot output for the combined legend of the heatmap of ligand-receptor interactions
    output$P16_lgd_ninter_ht <- renderPlot({
      draw(packLegend(lgdCT, lgdSoma))
    })

    
    
    #P30
    
      # --------------------------------------------------
      # Description: Create an eventReactive function that filters the Cortex_edges_sig.list$P30 data frame based on user input criteria
      # Output: A filtered data frame containing rows that match the user input criteria
      # --------------------------------------------------

    setkey(Cortex_edges_sig.list$P30, cluster_L, cluster_R, SOURCE.class, TARGET.class, LR_pair, numdup.LR_pair.interaction, numdup.LR_pair.cluster_L,
           numdup.LR_pair.cluster_R, pathway, S_inter, S_intra, S_inter_diff,
           Ligand.category1, Ligand.category2, Ligand.category3, Ligand.category4, Ligand.category5, Ligand.category6, Ligand.category7,
           Receptor.category1, Receptor.category2, Receptor.category3, Receptor.category4, Receptor.category5,
           Ligand.family1, Ligand.family2, Ligand.family3, Ligand.family4, Ligand.family5,
           Receptor.family1, Receptor.family2, Receptor.family3, Receptor.family4, Receptor.family5)
    
      # --------------------------------------------------
      # Create an eventReactive function that filters the Cortex_edges_sig.list$P30 data frame based on user input criteria
      # Output: A filtered data frame containing rows that match the user input criteria
      # --------------------------------------------------
    
    P30_htNinter_filter <- eventReactive(input$"htNinter_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$P30[Ligand.category1 %in% req(input$"Cortex_htNinter_Ligand.category")  |
                         Ligand.category2  %in% req(input$"Cortex_htNinter_Ligand.category") |
                         Ligand.category3  %in% req(input$"Cortex_htNinter_Ligand.category") |
                         Ligand.category4  %in% req(input$"Cortex_htNinter_Ligand.category") |
                         Ligand.category5  %in% req(input$"Cortex_htNinter_Ligand.category") |
                         Ligand.category6  %in% req(input$"Cortex_htNinter_Ligand.category") |
                         Ligand.category7  %in% req(input$"Cortex_htNinter_Ligand.category")
      ]
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_htNinter_Receptor.category") 
      ]
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family3  %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_htNinter_Ligand.family")
      ]
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family4  %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family5  %in% req(input$"Cortex_htNinter_Receptor.family")  
      ]
      mydata <- mydata[S_intra >= req(input$"htNinter_S_intra") |
                         is.na(S_intra)==T]
      mydata <- mydata[cluster_L %in% req(input$"htNinter_Sending.cluster") & 
                         cluster_R %in% req(input$"htNinter_Target.cluster") & 
                         SOURCE.class %in% req(input$"htNinter_SOURCE.class") & 
                         TARGET.class %in% req(input$"htNinter_TARGET.class") & 
                         LR_pair %in% req(input$"htNinter_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"htNinter_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"htNinter_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"htNinter_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_htNinter_pathway") & 
                         S_inter >= req(input$"htNinter_S_inter") &
                         S_inter_diff >= req(input$"htNinter_S_inter_diff")]
      
      # Convert the data.table back to a data.frame 
      as.data.frame(mydata)
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render an interactive heatmap plot using plotly for the filtered ligand-receptor interaction data
    # Input:
    #   - P30_htNinter_filter(): A reactive function that returns the filtered data frame based on user input criteria
    #   - Cortex_edges_sig.list$`P8`: A data frame containing ligand-receptor interaction data for the "P30" time point
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------
    
    output$P30_htNinter <- renderPlotly({
      
      # Check if the filtered data frame is empty
      if (nrow(req(P30_htNinter_filter())) == 0) return(NULL)
      
      # Generate an interactive heatmap plot using the scSeqComm_heatmaply_cardinality function
      scSeqComm_heatmaply_cardinality(as.data.frame(Cortex_edges_sig.list$P30),req(P30_htNinter_filter()), y = "cluster_L", x = "cluster_R", hover_anno="Number of interactions",
                                      ylab = "Ligand-expressing cell types", xlab = "Receptor-expressing cell types",  
                                      cl.order=order_celltype[order_celltype$celltype_label%in%unique(Cortex_edges_sig.list$P30$cluster_L),],
                                      key.title="Number of \ninteractions", hide_colorbar=F, margins=c(0,0,0,0), sizefont=10, 
                                      main="", scale="none", 
                                      Sending.cl=req(input$"htNinter_Sending.cluster"),
                                      Target.cl=req(input$"htNinter_Target.cluster"), SOURCE.CT=req(input$"htNinter_SOURCE.class"),TARGET.CT=req(input$"htNinter_TARGET.class"),
                                      col.heatmap=req(input$"ht_ninter_pal"), height=req(input$"htNinter_height"))
      
    })
    
      # --------------------------------------------------
      # Description: Create a text output that displays the number of selected ligand-receptor pairs and the number of pairs expressed in the P30 dataset with the current filters
      # --------------------------------------------------
    
    output$P30_htNinter_selectLR <- renderText({
      
      
      paste("Your selection contains", nrow(LRDB_mouse_Hm()), "ligand-receptor pairs", "|", 
            length(unique(P30_htNinter_filter()$LR_pair[which(P30_htNinter_filter()$cluster_L %in% req(input$"htNinter_Sending.cluster")  & 
                                                               P30_htNinter_filter()$cluster_R %in% req(input$"htNinter_Target.cluster")  &  
                                                               P30_htNinter_filter()$SOURCE.class %in% req(input$"htNinter_SOURCE.class")  &  
                                                               P30_htNinter_filter()$TARGET.class %in% req(input$"htNinter_TARGET.class") )])), 
            "expressed in the P30 dataset with the current filters")
      
      
      
      
    })
    
    # Create a plot output for the combined legend of the heatmap of ligand-receptor interactions
    output$P30_lgd_ninter_ht <- renderPlot({
      draw(packLegend(lgdCT, lgdSoma))
    })
    


    
    #Adult
    
      # --------------------------------------------------
      # Description: Create an eventReactive function that filters the Cortex_edges_sig.list$Adult data frame based on user input criteria
      # Output: A filtered data frame containing rows that match the user input criteria
      # --------------------------------------------------

    setkey(Cortex_edges_sig.list$Adult, cluster_L, cluster_R, SOURCE.class, TARGET.class, LR_pair, numdup.LR_pair.interaction, numdup.LR_pair.cluster_L,
           numdup.LR_pair.cluster_R, pathway, S_inter, S_intra, S_inter_diff,
           Ligand.category1, Ligand.category2, Ligand.category3, Ligand.category4, Ligand.category5, Ligand.category6, Ligand.category7,
           Receptor.category1, Receptor.category2, Receptor.category3, Receptor.category4, Receptor.category5,
           Ligand.family1, Ligand.family2, Ligand.family3, Ligand.family4, Ligand.family5,
           Receptor.family1, Receptor.family2, Receptor.family3, Receptor.family4, Receptor.family5)
    
      # --------------------------------------------------
      # Create an eventReactive function that filters the Cortex_edges_sig.list$Adult data frame based on user input criteria
      # Output: A filtered data frame containing rows that match the user input criteria
      # --------------------------------------------------

    Adult_htNinter_filter <- eventReactive(input$"htNinter_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$Adult[Ligand.category1 %in% req(input$"Cortex_htNinter_Ligand.category")  |
                           Ligand.category2  %in% req(input$"Cortex_htNinter_Ligand.category") |
                           Ligand.category3  %in% req(input$"Cortex_htNinter_Ligand.category") |
                           Ligand.category4  %in% req(input$"Cortex_htNinter_Ligand.category") |
                           Ligand.category5  %in% req(input$"Cortex_htNinter_Ligand.category") |
                           Ligand.category6  %in% req(input$"Cortex_htNinter_Ligand.category") |
                           Ligand.category7  %in% req(input$"Cortex_htNinter_Ligand.category")
      ]
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_htNinter_Receptor.category") 
      ]
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family3  %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_htNinter_Ligand.family")
      ]
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family4  %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family5  %in% req(input$"Cortex_htNinter_Receptor.family")  
      ]
      mydata <- mydata[S_intra >= req(input$"htNinter_S_intra") |
                         is.na(S_intra)==T]
      mydata <- mydata[cluster_L %in% req(input$"htNinter_Sending.cluster") & 
                         cluster_R %in% req(input$"htNinter_Target.cluster") & 
                         SOURCE.class %in% req(input$"htNinter_SOURCE.class") & 
                         TARGET.class %in% req(input$"htNinter_TARGET.class") & 
                         LR_pair %in% req(input$"htNinter_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"htNinter_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"htNinter_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"htNinter_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_htNinter_pathway") & 
                         S_inter >= req(input$"htNinter_S_inter") &
                         S_inter_diff >= req(input$"htNinter_S_inter_diff")]
      
      # Convert the data.table back to a data.frame 
      as.data.frame(mydata)
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render an interactive heatmap plot using plotly for the filtered ligand-receptor interaction data
    # Input:
    #   - Adult_htNinter_filter(): A reactive function that returns the filtered data frame based on user input criteria
    #   - Cortex_edges_sig.list$`P8`: A data frame containing ligand-receptor interaction data for the "Adult" time point
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------
    
    output$Adult_htNinter <- renderPlotly({
      
      # Check if the filtered data frame is empty
      if (nrow(req(Adult_htNinter_filter())) == 0) return(NULL)
      
      # Generate an interactive heatmap plot using the scSeqComm_heatmaply_cardinality function
      scSeqComm_heatmaply_cardinality(as.data.frame(Cortex_edges_sig.list$Adult),req(Adult_htNinter_filter()), y = "cluster_L", x = "cluster_R", hover_anno="Number of interactions",
                                      ylab = "Ligand-expressing cell types", xlab = "Receptor-expressing cell types",  
                                      cl.order=order_celltype[order_celltype$celltype_label%in%unique(Cortex_edges_sig.list$Adult$cluster_L),],
                                      key.title="Number of \ninteractions", hide_colorbar=F, margins=c(0,0,0,0), sizefont=10, 
                                      main="", scale="none", 
                                      Sending.cl=req(input$"htNinter_Sending.cluster"),
                                      Target.cl=req(input$"htNinter_Target.cluster"), SOURCE.CT=req(input$"htNinter_SOURCE.class"),TARGET.CT=req(input$"htNinter_TARGET.class"),
                                      col.heatmap=req(input$"ht_ninter_pal"), height=req(input$"htNinter_height"))
      
    })
    
      # --------------------------------------------------
      # Description: Create a text output that displays the number of selected ligand-receptor pairs and the number of pairs expressed in the Adult dataset with the current filters
      # --------------------------------------------------
    
    output$Adult_htNinter_selectLR <- renderText({
      
      
      paste("Your selection contains", nrow(LRDB_mouse_Hm()), "ligand-receptor pairs", "|", 
            length(unique(Adult_htNinter_filter()$LR_pair[which(Adult_htNinter_filter()$cluster_L %in% req(input$"htNinter_Sending.cluster")  & 
                                                               Adult_htNinter_filter()$cluster_R %in% req(input$"htNinter_Target.cluster")  &  
                                                               Adult_htNinter_filter()$SOURCE.class %in% req(input$"htNinter_SOURCE.class")  &  
                                                               Adult_htNinter_filter()$TARGET.class %in% req(input$"htNinter_TARGET.class") )])), 
            "expressed in the Adult dataset with the current filters")
      
      
      
      
    })
    
    # Create a plot output for the combined legend of the heatmap of ligand-receptor interactions
    output$Adult_lgd_ninter_ht <- renderPlot({
      draw(packLegend(lgdCT, lgdSoma))
    })
    
    
    
    
    #Cortex

      # --------------------------------------------------
      # Description: Create an eventReactive function that filters the Cortex_edges_sig.list$Cortex data frame based on user input criteria
      # Output: A filtered data frame containing rows that match the user input criteria
      # --------------------------------------------------

    setkey(Cortex_edges_sig.list$Cortex, cluster_L, cluster_R, SOURCE.class, TARGET.class, LR_pair, numdup.LR_pair.interaction, numdup.LR_pair.cluster_L,
           numdup.LR_pair.cluster_R, pathway, S_inter, S_intra, S_inter_diff,
           Ligand.category1, Ligand.category2, Ligand.category3, Ligand.category4, Ligand.category5, Ligand.category6, Ligand.category7,
           Receptor.category1, Receptor.category2, Receptor.category3, Receptor.category4, Receptor.category5,
           Ligand.family1, Ligand.family2, Ligand.family3, Ligand.family4, Ligand.family5,
           Receptor.family1, Receptor.family2, Receptor.family3, Receptor.family4, Receptor.family5)
    
    # --------------------------------------------------
      # Create an eventReactive function that filters the Cortex_edges_sig.list$Cortex data frame based on user input criteria
      # Output: A filtered data frame containing rows that match the user input criteria
      # --------------------------------------------------
    
    Cortex_htNinter_filter <- eventReactive(input$"htNinter_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$Cortex[Ligand.category1 %in% req(input$"Cortex_htNinter_Ligand.category")  |
                           Ligand.category2  %in% req(input$"Cortex_htNinter_Ligand.category") |
                           Ligand.category3  %in% req(input$"Cortex_htNinter_Ligand.category") |
                           Ligand.category4  %in% req(input$"Cortex_htNinter_Ligand.category") |
                           Ligand.category5  %in% req(input$"Cortex_htNinter_Ligand.category") |
                           Ligand.category6  %in% req(input$"Cortex_htNinter_Ligand.category") |
                           Ligand.category7  %in% req(input$"Cortex_htNinter_Ligand.category")
      ]
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_htNinter_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_htNinter_Receptor.category") 
      ]
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family3  %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_htNinter_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_htNinter_Ligand.family")
      ]
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family4  %in% req(input$"Cortex_htNinter_Receptor.family") |
                         Receptor.family5  %in% req(input$"Cortex_htNinter_Receptor.family")  
      ]
      mydata <- mydata[S_intra >= req(input$"htNinter_S_intra") |
                         is.na(S_intra)==T]
      mydata <- mydata[cluster_L %in% req(input$"htNinter_Sending.cluster") & 
                         cluster_R %in% req(input$"htNinter_Target.cluster") & 
                         SOURCE.class %in% req(input$"htNinter_SOURCE.class") & 
                         TARGET.class %in% req(input$"htNinter_TARGET.class") & 
                         LR_pair %in% req(input$"htNinter_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"htNinter_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"htNinter_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"htNinter_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_htNinter_pathway") & 
                         S_inter >= req(input$"htNinter_S_inter") &
                         S_inter_diff >= req(input$"htNinter_S_inter_diff")]
      
      # Convert the data.table back to a data.frame 
      as.data.frame(mydata)
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render an interactive heatmap plot using plotly for the filtered ligand-receptor interaction data
    # Input:
    #   - Cortex_htNinter_filter(): A reactive function that returns the filtered data frame based on user input criteria
    #   - Cortex_edges_sig.list$Cortex: A data frame containing ligand-receptor interaction data for Cortex
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------

    output$Cortex_htNinter <- renderPlotly({
      
      # Check if the filtered data frame is empty
      if (nrow(req(Cortex_htNinter_filter())) == 0) return(NULL)
      
      # Generate an interactive heatmap plot using the scSeqComm_contigheatmaply_cardinality function
      scSeqComm_contigheatmaply_cardinality(Cortex_edges_sig.list$Cortex,req(Cortex_htNinter_filter()), y = "cluster_L", x = "cluster_R", hover_anno="Number of interactions",
                                      ylab = "Ligand-expressing cell types", xlab = "Receptor-expressing cell types",  
                                      cl.order=order_celltype[order_celltype$celltype_label%in%unique(Cortex_edges_sig.list$Cortex$cluster_L),],
                                      key.title="Number of \ninteractions", hide_colorbar=F, margins=c(0,0,0,0), sizefont=10, 
                                      main="", scale="none", 
                                      Sending.cl=req(input$"htNinter_Sending.cluster"),
                                      Target.cl=req(input$"htNinter_Target.cluster"), SOURCE.CT=req(input$"htNinter_SOURCE.class"),TARGET.CT=req(input$"htNinter_TARGET.class"),
                                      col.heatmap=req(input$"ht_ninter_pal"), height=req(input$"htNinter_height"))
      
    })
    
      # --------------------------------------------------
      # Description: Create a text output that displays the number of selected ligand-receptor pairs and the number of pairs expressed in the Cortex with the current filters
      # --------------------------------------------------

    output$Cortex_htNinter_selectLR <- renderText({
      
      paste("Your selection contains", nrow(LRDB_mouse_Hm()), "ligand-receptor pairs", "|", 
            length(unique(Cortex_htNinter_filter()$LR_pair[which(Cortex_htNinter_filter()$cluster_L %in% req(input$"htNinter_Sending.cluster")  & 
                                                                  Cortex_htNinter_filter()$cluster_R %in% req(input$"htNinter_Target.cluster")  &  
                                                                  Cortex_htNinter_filter()$SOURCE.class %in% req(input$"htNinter_SOURCE.class")  &  
                                                                  Cortex_htNinter_filter()$TARGET.class %in% req(input$"htNinter_TARGET.class") )])), 
            "expressed across all datasets with the current filters")
      
    })
    
    # Create a plot output for the combined legend of the heatmap of ligand-receptor interactions
    output$Cortex_lgd_ninter_ht <- renderPlot({
      draw(packLegend(lgdCT, lgdSoma))
    })

    
    ############## LR interactions

    # --------------------------------------------------
    # Description: update picker inputs for various filters in the ligand-receptor interactions
    # --------------------------------------------------
    
    observe({
        updatePickerInput(
           session = session,
            inputId = "Cortex_LRpairs_SOURCE.class",
            choices = sort(unique(Cortex_edges_sig.list$Cortex$SOURCE.class)),
            selected = sort(unique(Cortex_edges_sig.list$Cortex$SOURCE.class))
        )
    })
    
    observe({
      updatePickerInput(
        session = session,
        inputId = "Cortex_LRpairs_TARGET.class",
        choices = sort(unique(Cortex_edges_sig.list$Cortex$TARGET.class)),
        selected = sort(unique(Cortex_edges_sig.list$Cortex$TARGET.class))
      )
    })
    
    observe({
      updatePickerInput(
        session = session,
        inputId = "Cortex_LRpairs_cluster_L",
        choices = sort(unique(Cortex_edges_sig.list$Cortex$cluster_L)),
        selected = sort(unique(Cortex_edges_sig.list$Cortex$cluster_L))
      )
    })
    
    observe({
      updatePickerInput(
        session = session,
        inputId = "Cortex_LRpairs_cluster_R",
        choices = sort(unique(Cortex_edges_sig.list$Cortex$cluster_R)),
        selected = sort(unique(Cortex_edges_sig.list$Cortex$cluster_R))
      )
    })
    
    observe({
        updatePickerInput(
            session = session,
            inputId = "Cortex_LRpairs_LR",
            choices = sort(unique(c(LRintercellNetworkDB$ligand_receptor, Cortex_edges_sig.list$Cortex$LR_pair))),
            selected = c( "Cbln4 - Dcc", "Cbln4 - Glud1", "Cbln4 - Neo1" )
        )
    })
    
    observe({
      updatePickerInput(
        session = session,
        inputId = "Cortex_LRpairs_pathway",
        choices = sort(unique(Cortex_edges_sig.list$Cortex$pathway)),
        selected = sort(unique(Cortex_edges_sig.list$Cortex$pathway))
      )
    })
    
    observe({
        updatePickerInput(
            session = session,
            inputId = "Cortex_LRpairs_Ligand.category",
            choices = sort(unique(unlist(strsplit(LRintercellNetworkDB$ligand_category, fixed=T, split=";")))),
            selected = sort(unique(unlist(strsplit(LRintercellNetworkDB$ligand_category, fixed=T, split=";"))))
        )
    })
    
    observe({
       updatePickerInput(
            session = session,
            inputId = "Cortex_LRpairs_Receptor.category",
            choices = sort(unique(unlist(strsplit(LRintercellNetworkDB$receptor_category, fixed=T, split=";")))),
            selected = sort(unique(unlist(strsplit(LRintercellNetworkDB$receptor_category, fixed=T, split=";"))))
        )
    })
    
    observe({
        updatePickerInput(
            session = session,
            inputId = "Cortex_LRpairs_Ligand.family",
            choices = sort(unique(unlist(strsplit(LRintercellNetworkDB$ligand_family, fixed=T, split=";")))),
            selected = sort(unique(unlist(strsplit(LRintercellNetworkDB$ligand_family, fixed=T, split=";"))))
        )
    })
    
    observe({
        updatePickerInput(
            session = session,
            inputId = "Cortex_LRpairs_Receptor.family",
            choices = sort(unique(unlist(strsplit(LRintercellNetworkDB$receptor_family, fixed=T, split=";")))),
            selected = sort(unique(unlist(strsplit(LRintercellNetworkDB$receptor_family, fixed=T, split=";"))))
        )
    })
    
    observe({
      updatePickerInput(
        session = session,
        inputId = "Cortex_LRpairs_numdup.LR_pair.interaction",
        choices = str_sort(unique(Cortex_edges_sig.list$Cortex$numdup.LR_pair.interaction), numeric=T),
        selected = str_sort(unique(Cortex_edges_sig.list$Cortex$numdup.LR_pair.interaction), numeric=T)
      )
    })
    
    observe({
      updatePickerInput(
        session = session,
        inputId = "Cortex_LRpairs_numdup.LR_pair.cluster_L",
        choices = str_sort(unique(Cortex_edges_sig.list$Cortex$numdup.LR_pair.cluster_L), numeric=T),
        selected = str_sort(unique(Cortex_edges_sig.list$Cortex$numdup.LR_pair.cluster_L), numeric=T)
      )
    })
    
    observe({
      updatePickerInput(
        session = session,
        inputId = "Cortex_LRpairs_numdup.LR_pair.cluster_R",
        choices = str_sort(unique(Cortex_edges_sig.list$Cortex$numdup.LR_pair.cluster_R), numeric=T),
        selected = str_sort(unique(Cortex_edges_sig.list$Cortex$numdup.LR_pair.cluster_R), numeric=T)
      )
    })

    observe({
      updatePickerInput(
        session = session,
        inputId = "Cortex_LRpairs_x.axis",
        choices = c("cluster_L","cluster_R","interaction","LR_pair", "pathway"),
        selected = "interaction"
      )
    })

      observe({
        updatePickerInput(
          session = session,
          inputId = "Cortex_LRpairs_y.axis",
          choices = c("cluster_L","cluster_R","interaction","LR_pair","pathway"),
          selected = "LR_pair"
      )
      })
    
      observe({
        updatePickerInput(
          session = session,
          inputId = "Cortex_LRpairs_facetgrid.x",
          choices = c("cluster_L","cluster_R","interaction","LR_pair", "pathway","NULL"),
          selected = "cluster_L"
      )
      })

      observe({
        updatePickerInput(
          session = session,
          inputId = "Cortex_LRpairs_facetgrid.y",
          choices = c("cluster_L","cluster_R","interaction","LR_pair", "pathway","NULL"),
          selected = "NULL"
      )
      })
      
      observe({
        updatePickerInput(
          session = session,
          inputId = "Cortex_LRpairs_pt.size",
          choices = c("S_inter","S_intra","S_inter_diff"),
          selected = "S_inter"
        )
      })
      
      observe({
        updatePickerInput(
          session = session,
          inputId = "Cortex_LRpairs_fill",
          choices = c("S_inter","S_intra","S_inter_diff"),
          selected = "S_intra"
        )
      })
    
    #E18.5-P0

      # --------------------------------------------------
      # Description: An eventReactive function  to filter a dataset based on various user-selected inputs and return the filtered dataset. 
      # --------------------------------------------------

    P0_LRpairs_filter <- eventReactive(input$"Cortex_LRpairs_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$`E18.5-P0`[Ligand.category1 %in% req(input$"Cortex_LRpairs_Ligand.category")  |
                        Ligand.category2  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category3  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category4  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category5  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category6  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category7  %in% req(input$"Cortex_LRpairs_Ligand.category")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.category' columns match the user-selected receptor category
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_LRpairs_Receptor.category")
      ]

      # Filter the 'mydata' by rows where any of the 'Ligand.family' columns match the user-selected ligand family
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family3 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_LRpairs_Ligand.family")
      ]

      # Filter the 'mydata' by rows where any of the 'Receptor.family' columns match the user-selected receptor family
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family4 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family5 %in% req(input$"Cortex_LRpairs_Receptor.family")
      ]

      # Filter the 'mydata' by rows where the 'S_intra' column is greater than or equal to the user-selected 'S_intra' value or is NA
      mydata <- mydata[S_intra >= req(input$"Cortex_LRpairs_S_intra") |
                         is.na(S_intra)==T]

      # Perform a final filtering step on 'mydata' based on multiple user-selected inputs related to various columns in the data table                   
      mydata <- mydata[cluster_L %in% req(input$"Cortex_LRpairs_cluster_L") & 
                         cluster_R %in% req(input$"Cortex_LRpairs_cluster_R") & 
                         SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class") & 
                         TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") & 
                         LR_pair %in% req(input$"Cortex_LRpairs_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"Cortex_LRpairs_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_LRpairs_pathway") & 
                         S_inter >= req(input$"Cortex_LRpairs_S_inter") &
                         S_inter_diff >= req(input$"Cortex_LRpairs_S_inter_diff")]
      }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render a heatmap plot using the 'scSeqComm_heatmaply' function for the filtered ligand-receptor interaction data
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------
    
    output$P0_LRpairs <- renderPlotly({

      input$"Cortex_LRsetting_button"

      # Create a heatmap plot using the filtered data from 'P0_LRpairs_filter()' and various user-selected inputs
      scSeqComm_heatmaply(req(P0_LRpairs_filter()), x=isolate(input$"Cortex_LRpairs_x.axis"), y=isolate(input$"Cortex_LRpairs_y.axis"), 
                          fill = isolate(input$"Cortex_LRpairs_fill"), size = isolate(input$"Cortex_LRpairs_pt.size"),  
                          height = isolate(input$"LRpairs_height"),width= isolate(input$"LRpairs_width"),
                          gap_width= isolate(input$"LRpairs_step_width"), gap_height= isolate(input$"LRpairs_step_height"),
                          xlab="", ylab = "",  gdt_pal=isolate(input$"Cortex_LRpairs_pal"),
                          limit_fill = c(0, 1), limit_size = c(0, 1), breaks_size = c(0.2, 
                                                                                      0.5, 0.8), 
                          facet_grid_x = isolate(input$"Cortex_LRpairs_facetgrid.x"), facet_grid_y = isolate(input$"Cortex_LRpairs_facetgrid.y"))
    })

      # --------------------------------------------------
      # Description: This eventReactive function listens to the "Cortex_LRpairs_button" input, which triggers the filtering process
      # Output: The filtered ligand-receptor interaction data
      # --------------------------------------------------
    
    LRDB_mouse_LRpairs <- eventReactive(input$"Cortex_LRpairs_button",{

       # Filter the data from the 'LRintercellNetworkDB' data frame based on various user-selected inputs
      myLRdata <-  LRintercellNetworkDB[which( LRintercellNetworkDB$ligand_receptor %in% req(input$"Cortex_LRpairs_LR")   &

      # Match the 'ligand_receptor' column with the user-selected 'Ligand.category', 'Receptor.category', 'Ligand.family', and 'Receptor.family' inputs
                                                 rownames(LRintercellNetworkDB) %in% rownames(LRintercellNetworkDB[unlist(sapply(input$"Cortex_LRpairs_Ligand.category", grep, LRintercellNetworkDB$ligand_category, fixed=T, USE.NAMES = F)),]) &
                                                 rownames(LRintercellNetworkDB) %in% rownames(LRintercellNetworkDB[unlist(sapply(input$"Cortex_LRpairs_Receptor.category", grep, LRintercellNetworkDB$receptor_category, fixed=T, USE.NAMES = F)),]) &
                                                 rownames(LRintercellNetworkDB) %in% rownames(LRintercellNetworkDB[unlist(sapply(input$"Cortex_LRpairs_Ligand.family", grep, LRintercellNetworkDB$ligand_family, fixed=T, USE.NAMES = F)),]) &
                                                 rownames(LRintercellNetworkDB) %in% rownames(LRintercellNetworkDB[unlist(sapply(input$"Cortex_LRpairs_Receptor.family", grep, LRintercellNetworkDB$receptor_family, fixed=T, USE.NAMES = F)),])),]
      
      # Return the filtered data
      myLRdata
      
    },ignoreNULL = FALSE)
    
      # --------------------------------------------------
      # Description: create a text output that displays information about the user's selection
      # --------------------------------------------------

    output$P0_selectLR <- renderText({
      
      paste("Your selection contains", nrow(LRDB_mouse_LRpairs()), "ligand-receptor pairs", "|", 
            length(unique(P0_LRpairs_filter()$LR_pair[which(P0_LRpairs_filter()$cluster_L %in% req(input$"Cortex_LRpairs_cluster_L")  & 
                                                               P0_LRpairs_filter()$cluster_R %in% req(input$"Cortex_LRpairs_cluster_R")  &  
                                                               P0_LRpairs_filter()$SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class")  &  
                                                               P0_LRpairs_filter()$TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") )])), 
            "expressed in the E18.5-P0 dataset with the current filters")
      
    }) 
    
    #P1-P2

      # --------------------------------------------------
      # Description: An eventReactive function  to filter a dataset based on various user-selected inputs and return the filtered dataset. 
      # --------------------------------------------------

    P2_LRpairs_filter <- eventReactive(input$"Cortex_LRpairs_button",{
      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$`P1-P2`[Ligand.category1 %in% req(input$"Cortex_LRpairs_Ligand.category")  |
                        Ligand.category2  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category3  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category4  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category5  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category6  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category7  %in% req(input$"Cortex_LRpairs_Ligand.category")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.category' columns match the user-selected receptor category
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_LRpairs_Receptor.category")
      ]

      # Further filter the 'mydata' by rows where any of the 'Ligand.family' columns match the user-selected ligand family
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family3 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_LRpairs_Ligand.family")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.family' columns match the user-selected receptor family
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family4 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family5 %in% req(input$"Cortex_LRpairs_Receptor.family")
      ]

      # Filter the 'mydata' by rows where the 'S_intra' column is greater than or equal to the user-selected 'S_intra' value or is NA
      mydata <- mydata[S_intra >= req(input$"Cortex_LRpairs_S_intra") |
                         is.na(S_intra)==T]

      # Perform a final filtering step on 'mydata' based on multiple user-selected inputs related to various columns in the data table                         
      mydata <- mydata[cluster_L %in% req(input$"Cortex_LRpairs_cluster_L") & 
                         cluster_R %in% req(input$"Cortex_LRpairs_cluster_R") & 
                         SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class") & 
                         TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") & 
                         LR_pair %in% req(input$"Cortex_LRpairs_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"Cortex_LRpairs_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_LRpairs_pathway") & 
                         S_inter >= req(input$"Cortex_LRpairs_S_inter") &
                         S_inter_diff >= req(input$"Cortex_LRpairs_S_inter_diff")]
    }, ignoreNULL = FALSE)
    
    
    # --------------------------------------------------
    # Description: Render a heatmap plot using the 'scSeqComm_heatmaply' function for the filtered ligand-receptor interaction data
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------
    
    output$P2_LRpairs <- renderPlotly({

      input$"Cortex_LRsetting_button"
      
      # Render a heatmap plot using the 'scSeqComm_heatmaply' function for the filtered ligand-receptor interaction data
      scSeqComm_heatmaply(req(P2_LRpairs_filter()), x=isolate(input$"Cortex_LRpairs_x.axis"), y=isolate(input$"Cortex_LRpairs_y.axis"), 
                          fill = isolate(input$"Cortex_LRpairs_fill"), size = isolate(input$"Cortex_LRpairs_pt.size"),  
                          height = isolate(input$"LRpairs_height"),width= isolate(input$"LRpairs_width"),
                          gap_width= isolate(input$"LRpairs_step_width"), gap_height= isolate(input$"LRpairs_step_height"),
                          xlab="", ylab = "",  gdt_pal=isolate(input$"Cortex_LRpairs_pal"),
                          limit_fill = c(0, 1), limit_size = c(0, 1), breaks_size = c(0.2, 
                                                                                      0.5, 0.8), 
                          facet_grid_x = isolate(input$"Cortex_LRpairs_facetgrid.x"), facet_grid_y = isolate(input$"Cortex_LRpairs_facetgrid.y"))
      
    })

      # --------------------------------------------------
      # Description: create a text output that displays information about the user's selection
      # --------------------------------------------------
    
    output$P2_selectLR <- renderText({
      
      
      paste("Your selection contains", nrow(LRDB_mouse_LRpairs()), "ligand-receptor pairs", "|", 
            length(unique(P2_LRpairs_filter()$LR_pair[which(P2_LRpairs_filter()$cluster_L %in% req(input$"Cortex_LRpairs_cluster_L")  & 
                                                              P2_LRpairs_filter()$cluster_R %in% req(input$"Cortex_LRpairs_cluster_R")  &  
                                                              P2_LRpairs_filter()$SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class")  &  
                                                              P2_LRpairs_filter()$TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") )])), 
            "expressed in the P1-P2 dataset with the current filters")

    })
    
    #P4-P5

      # --------------------------------------------------
      # Description: An eventReactive function  to filter a dataset based on various user-selected inputs and return the filtered dataset. 
      # --------------------------------------------------

    P5_LRpairs_filter <- eventReactive(input$"Cortex_LRpairs_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$`P4-P5`[Ligand.category1 %in% req(input$"Cortex_LRpairs_Ligand.category")  |
                        Ligand.category2  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category3  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category4  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category5  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category6  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category7  %in% req(input$"Cortex_LRpairs_Ligand.category")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.category' columns match the user-selected receptor category
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_LRpairs_Receptor.category")
      ]

      # Further filter the 'mydata' by rows where any of the 'Ligand.family' columns match the user-selected ligand family
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family3 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_LRpairs_Ligand.family")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.family' columns match the user-selected receptor family
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family4 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family5 %in% req(input$"Cortex_LRpairs_Receptor.family")
      ]

      # Further filter the 'mydata' by rows where the 'S_intra' column is greater than or equal to the user-selected value
      mydata <- mydata[S_intra >= req(input$"Cortex_LRpairs_S_intra") |
                         is.na(S_intra)==T]

      # Perform a final filtering step on 'mydata' based on multiple user-selected inputs related to various columns in the data table                        
      mydata <- mydata[cluster_L %in% req(input$"Cortex_LRpairs_cluster_L") & 
                         cluster_R %in% req(input$"Cortex_LRpairs_cluster_R") & 
                         SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class") & 
                         TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") & 
                         LR_pair %in% req(input$"Cortex_LRpairs_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"Cortex_LRpairs_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_LRpairs_pathway") & 
                         S_inter >= req(input$"Cortex_LRpairs_S_inter") &
                         S_inter_diff >= req(input$"Cortex_LRpairs_S_inter_diff")]
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render a heatmap plot using the 'scSeqComm_heatmaply' function for the filtered ligand-receptor interaction data
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------
    
    output$P5_LRpairs <- renderPlotly({

      input$"Cortex_LRsetting_button"
      
      # Generate an interactive heatmap plot using the 'scSeqComm_heatmaply' function
      scSeqComm_heatmaply(req(P5_LRpairs_filter()), x=isolate(input$"Cortex_LRpairs_x.axis"), y=isolate(input$"Cortex_LRpairs_y.axis"), 
                          fill = isolate(input$"Cortex_LRpairs_fill"), size = isolate(input$"Cortex_LRpairs_pt.size"),  
                          height = isolate(input$"LRpairs_height"),width= isolate(input$"LRpairs_width"),
                          gap_width= isolate(input$"LRpairs_step_width"), gap_height= isolate(input$"LRpairs_step_height"),
                          xlab="", ylab = "",  gdt_pal=isolate(input$"Cortex_LRpairs_pal"),
                          limit_fill = c(0, 1), limit_size = c(0, 1), breaks_size = c(0.2, 
                                                                                      0.5, 0.8), 
                          facet_grid_x = isolate(input$"Cortex_LRpairs_facetgrid.x"), facet_grid_y = isolate(input$"Cortex_LRpairs_facetgrid.y"))
      
    })

      # --------------------------------------------------
      # Description: create a text output that displays information about the user's selection
      # --------------------------------------------------
    
    output$P5_selectLR <- renderText({
      
      paste("Your selection contains", nrow(LRDB_mouse_LRpairs()), "ligand-receptor pairs", "|", 
            length(unique(P5_LRpairs_filter()$LR_pair[which(P5_LRpairs_filter()$cluster_L %in% req(input$"Cortex_LRpairs_cluster_L")  & 
                                                              P5_LRpairs_filter()$cluster_R %in% req(input$"Cortex_LRpairs_cluster_R")  &  
                                                              P5_LRpairs_filter()$SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class")  &  
                                                              P5_LRpairs_filter()$TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") )])), 
            "expressed in the P4-P5 dataset with the current filters")

    })
    
    
    #P8

      # --------------------------------------------------
      # Description: An eventReactive function  to filter a dataset based on various user-selected inputs and return the filtered dataset. 
      # --------------------------------------------------

    P8_LRpairs_filter <- eventReactive(input$"Cortex_LRpairs_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$P8[Ligand.category1 %in% req(input$"Cortex_LRpairs_Ligand.category")  |
                        Ligand.category2  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category3  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category4  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category5  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                        Ligand.category6  %in% req(input$"Cortex_LRpairs_Ligand.category")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.category' columns match the user-selected receptor category
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_LRpairs_Receptor.category")
      ]

      # Further filter the 'mydata' by rows where any of the 'Ligand.family' columns match the user-selected ligand family
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family3 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_LRpairs_Ligand.family")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.family' columns match the user-selected receptor family
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family4 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family5 %in% req(input$"Cortex_LRpairs_Receptor.family")
      ]

      # Further filter the 'mydata' by rows where the 'S_intra' column is greater than or equal to the user-selected 'S_intra' value or is NA
      mydata <- mydata[S_intra >= req(input$"Cortex_LRpairs_S_intra") |
                         is.na(S_intra)==T]

      # Perform a final filtering step on 'mydata' based on multiple user-selected inputs related to various columns in the data table                     
      mydata <- mydata[cluster_L %in% req(input$"Cortex_LRpairs_cluster_L") & 
                         cluster_R %in% req(input$"Cortex_LRpairs_cluster_R") & 
                         SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class") & 
                         TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") & 
                         LR_pair %in% req(input$"Cortex_LRpairs_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"Cortex_LRpairs_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_LRpairs_pathway") & 
                         S_inter >= req(input$"Cortex_LRpairs_S_inter") &
                         S_inter_diff >= req(input$"Cortex_LRpairs_S_inter_diff")]
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render a heatmap plot using the 'scSeqComm_heatmaply' function for the filtered ligand-receptor interaction data
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------
    
    output$P8_LRpairs <- renderPlotly({

      input$"Cortex_LRsetting_button"
      
      # Generate an interactive heatmap plot using the 'scSeqComm_heatmaply' function
      scSeqComm_heatmaply(req(P8_LRpairs_filter()), x=isolate(input$"Cortex_LRpairs_x.axis"), y=isolate(input$"Cortex_LRpairs_y.axis"), 
                          fill = isolate(input$"Cortex_LRpairs_fill"), size = isolate(input$"Cortex_LRpairs_pt.size"),  
                          height = isolate(input$"LRpairs_height"),width= isolate(input$"LRpairs_width"),
                          gap_width= isolate(input$"LRpairs_step_width"), gap_height= isolate(input$"LRpairs_step_height"),
                          xlab="", ylab = "",  gdt_pal=isolate(input$"Cortex_LRpairs_pal"),
                          limit_fill = c(0, 1), limit_size = c(0, 1), breaks_size = c(0.2, 
                                                                                      0.5, 0.8), 
                          facet_grid_x = isolate(input$"Cortex_LRpairs_facetgrid.x"), facet_grid_y = isolate(input$"Cortex_LRpairs_facetgrid.y"))
      
    })

      # --------------------------------------------------
      # Description: create a text output that displays information about the user's selection
      # --------------------------------------------------
    
    output$P8_selectLR <- renderText({
      
      paste("Your selection contains", nrow(LRDB_mouse_LRpairs()), "ligand-receptor pairs", "|", 
            length(unique(P8_LRpairs_filter()$LR_pair[which(P8_LRpairs_filter()$cluster_L %in% req(input$"Cortex_LRpairs_cluster_L")  & 
                                                              P8_LRpairs_filter()$cluster_R %in% req(input$"Cortex_LRpairs_cluster_R")  &  
                                                              P8_LRpairs_filter()$SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class")  &  
                                                              P8_LRpairs_filter()$TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") )])), 
            "expressed in the P8 dataset with the current filters")

    })
  
    
    #P16

      # --------------------------------------------------
      # Description: An eventReactive function  to filter a dataset based on various user-selected inputs and return the filtered dataset. 
      # --------------------------------------------------

    P16_LRpairs_filter <- eventReactive(input$"Cortex_LRpairs_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$P16[Ligand.category1 %in% req(input$"Cortex_LRpairs_Ligand.category")  |
                         Ligand.category2  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                         Ligand.category3  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                         Ligand.category4  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                         Ligand.category5  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                         Ligand.category6  %in% req(input$"Cortex_LRpairs_Ligand.category")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.category' columns match the user-selected receptor category
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_LRpairs_Receptor.category")
      ]

      # Further filter the 'mydata' by rows where any of the 'Ligand.family' columns match the user-selected ligand family
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family3 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_LRpairs_Ligand.family")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.family' columns match the user-selected receptor family
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family4 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family5 %in% req(input$"Cortex_LRpairs_Receptor.family")
      ]

      # Filter the 'mydata' by rows where the 'S_intra' column is greater than or equal to the user-selected 'S_intra' value or is NA
      mydata <- mydata[S_intra >= req(input$"Cortex_LRpairs_S_intra") |
                         is.na(S_intra)==T]

      # Perform a final filtering step on 'mydata' based on multiple user-selected inputs related to various columns in the data table                   
      mydata <- mydata[cluster_L %in% req(input$"Cortex_LRpairs_cluster_L") & 
                         cluster_R %in% req(input$"Cortex_LRpairs_cluster_R") & 
                         SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class") & 
                         TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") & 
                         LR_pair %in% req(input$"Cortex_LRpairs_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"Cortex_LRpairs_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_LRpairs_pathway") & 
                         S_inter >= req(input$"Cortex_LRpairs_S_inter") &
                         S_inter_diff >= req(input$"Cortex_LRpairs_S_inter_diff")]
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render a heatmap plot using the 'scSeqComm_heatmaply' function for the filtered ligand-receptor interaction data
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------
    
    output$P16_LRpairs <- renderPlotly({

      input$"Cortex_LRsetting_button"

      # Create a plotly heatmap plot using the 'scSeqComm_heatmaply' function for the filtered ligand-receptor interaction data
      scSeqComm_heatmaply(req(P16_LRpairs_filter()), x=isolate(input$"Cortex_LRpairs_x.axis"), y=isolate(input$"Cortex_LRpairs_y.axis"), 
                          fill = isolate(input$"Cortex_LRpairs_fill"), size = isolate(input$"Cortex_LRpairs_pt.size"),  
                          height = isolate(input$"LRpairs_height"),width= isolate(input$"LRpairs_width"),
                          gap_width= isolate(input$"LRpairs_step_width"), gap_height= isolate(input$"LRpairs_step_height"),
                          xlab="", ylab = "",  gdt_pal=isolate(input$"Cortex_LRpairs_pal"),
                          limit_fill = c(0, 1), limit_size = c(0, 1), breaks_size = c(0.2, 
                                                                                      0.5, 0.8), 
                          facet_grid_x = isolate(input$"Cortex_LRpairs_facetgrid.x"), facet_grid_y = isolate(input$"Cortex_LRpairs_facetgrid.y"))
      
    })
    
      # --------------------------------------------------
      # Description: create a text output that displays information about the user's selection
      # --------------------------------------------------

    output$P16_selectLR <- renderText({
      
      
      paste("Your selection contains", nrow(LRDB_mouse_LRpairs()), "ligand-receptor pairs", "|", 
            length(unique(P16_LRpairs_filter()$LR_pair[which(P16_LRpairs_filter()$cluster_L %in% req(input$"Cortex_LRpairs_cluster_L")  & 
                                                              P16_LRpairs_filter()$cluster_R %in% req(input$"Cortex_LRpairs_cluster_R")  &  
                                                              P16_LRpairs_filter()$SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class")  &  
                                                              P16_LRpairs_filter()$TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") )])), 
            "expressed in the P16 dataset with the current filters")

    })

    #P30

      # --------------------------------------------------
      # Description: An eventReactive function  to filter a dataset based on various user-selected inputs and return the filtered dataset. 
      # --------------------------------------------------

    P30_LRpairs_filter <- eventReactive(input$"Cortex_LRpairs_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$P30[Ligand.category1 %in% req(input$"Cortex_LRpairs_Ligand.category")  |
                         Ligand.category2  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                         Ligand.category3  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                         Ligand.category4  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                         Ligand.category5  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                         Ligand.category6  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                         Ligand.category7  %in% req(input$"Cortex_LRpairs_Ligand.category")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.category' columns match the user-selected receptor category
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_LRpairs_Receptor.category")
      ]

      # Filter the 'mydata' by rows where any of the 'Ligand.family' columns match the user-selected ligand family
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family3 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_LRpairs_Ligand.family")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.family' columns match the user-selected receptor family
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family4 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family5 %in% req(input$"Cortex_LRpairs_Receptor.family")
      ]

      # Filter the 'mydata' by rows where the 'S_intra' column is greater than or equal to the user-selected 'S_intra' value or is NA
      mydata <- mydata[S_intra >= req(input$"Cortex_LRpairs_S_intra") |
                         is.na(S_intra)==T]

      # Perform a final filtering step on 'mydata' based on multiple user-selected inputs related to various columns in the data table                     
      mydata <- mydata[cluster_L %in% req(input$"Cortex_LRpairs_cluster_L") & 
                         cluster_R %in% req(input$"Cortex_LRpairs_cluster_R") & 
                         SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class") & 
                         TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") & 
                         LR_pair %in% req(input$"Cortex_LRpairs_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"Cortex_LRpairs_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_LRpairs_pathway") & 
                         S_inter >= req(input$"Cortex_LRpairs_S_inter") &
                         S_inter_diff >= req(input$"Cortex_LRpairs_S_inter_diff")]
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render a heatmap plot using the 'scSeqComm_heatmaply' function for the filtered ligand-receptor interaction data
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------
    
    output$P30_LRpairs <- renderPlotly({
      input$"Cortex_LRsetting_button"
      
      # Create a heatmap plot using the filtered data from 'P30_LRpairs_filter()' and various user-selected inputs
      scSeqComm_heatmaply(req(P30_LRpairs_filter()), x=isolate(input$"Cortex_LRpairs_x.axis"), y=isolate(input$"Cortex_LRpairs_y.axis"), 
                          fill = isolate(input$"Cortex_LRpairs_fill"), size = isolate(input$"Cortex_LRpairs_pt.size"),  
                          height = isolate(input$"LRpairs_height"),width= isolate(input$"LRpairs_width"),
                          gap_width= isolate(input$"LRpairs_step_width"), gap_height= isolate(input$"LRpairs_step_height"),
                          xlab="", ylab = "",  gdt_pal=isolate(input$"Cortex_LRpairs_pal"),
                          limit_fill = c(0, 1), limit_size = c(0, 1), breaks_size = c(0.2, 
                                                                                      0.5, 0.8), 
                          facet_grid_x = isolate(input$"Cortex_LRpairs_facetgrid.x"), facet_grid_y = isolate(input$"Cortex_LRpairs_facetgrid.y"))
      
    })
    
      # --------------------------------------------------
      # Description: create a text output that displays information about the user's selection
      # --------------------------------------------------

    output$P30_selectLR <- renderText({
      
      
      paste("Your selection contains", nrow(LRDB_mouse_LRpairs()), "ligand-receptor pairs", "|", 
            length(unique(P30_LRpairs_filter()$LR_pair[which(P30_LRpairs_filter()$cluster_L %in% req(input$"Cortex_LRpairs_cluster_L")  & 
                                                              P30_LRpairs_filter()$cluster_R %in% req(input$"Cortex_LRpairs_cluster_R")  &  
                                                              P30_LRpairs_filter()$SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class")  &  
                                                              P30_LRpairs_filter()$TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") )])), 
            "expressed in the P30 dataset with the current filters")

    })
    
    #Adult

      # --------------------------------------------------
      # Description: An eventReactive function  to filter a dataset based on various user-selected inputs and return the filtered dataset. 
      # --------------------------------------------------

    Adult_LRpairs_filter <- eventReactive(input$"Cortex_LRpairs_button",{

      # Filter the data using data.table syntax
      mydata <- Cortex_edges_sig.list$Adult[Ligand.category1 %in% req(input$"Cortex_LRpairs_Ligand.category")  |
                           Ligand.category2  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                           Ligand.category3  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                           Ligand.category4  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                           Ligand.category5  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                           Ligand.category6  %in% req(input$"Cortex_LRpairs_Ligand.category") |
                           Ligand.category7  %in% req(input$"Cortex_LRpairs_Ligand.category")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.category' columns match the user-selected receptor category
      mydata <- mydata[Receptor.category1 %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category2  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category3  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category4  %in% req(input$"Cortex_LRpairs_Receptor.category") |
                         Receptor.category5  %in% req(input$"Cortex_LRpairs_Receptor.category")
      ]

      # Further filter the 'mydata' by rows where any of the 'Ligand.family' columns match the user-selected ligand family
      mydata <- mydata[Ligand.family1   %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family2 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family3 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family4 %in% req(input$"Cortex_LRpairs_Ligand.family") |
                         Ligand.family5 %in% req(input$"Cortex_LRpairs_Ligand.family")
      ]

      # Further filter the 'mydata' by rows where any of the 'Receptor.family' columns match the user-selected receptor family
      mydata <- mydata[Receptor.family1 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family2 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family3 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family4 %in% req(input$"Cortex_LRpairs_Receptor.family") |
                         Receptor.family5 %in% req(input$"Cortex_LRpairs_Receptor.family")
      ]

      # Further filter the 'mydata' by rows where any of the 'S_intra' columns match the user-selected S_intra
      mydata <- mydata[S_intra >= req(input$"Cortex_LRpairs_S_intra") |
                         is.na(S_intra)==T]
      
      # Perform a final filtering step on 'mydata' based on multiple user-selected inputs related to various columns in the data table 
      mydata <- mydata[cluster_L %in% req(input$"Cortex_LRpairs_cluster_L") & 
                         cluster_R %in% req(input$"Cortex_LRpairs_cluster_R") & 
                         SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class") & 
                         TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") & 
                         LR_pair %in% req(input$"Cortex_LRpairs_LR") & 
                         numdup.LR_pair.interaction %in% req(input$"Cortex_LRpairs_numdup.LR_pair.interaction") & 
                         numdup.LR_pair.cluster_L %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_L") & 
                         numdup.LR_pair.cluster_R  %in% req(input$"Cortex_LRpairs_numdup.LR_pair.cluster_R") & 
                         pathway  %in% req(input$"Cortex_LRpairs_pathway") & 
                         S_inter >= req(input$"Cortex_LRpairs_S_inter") &
                         S_inter_diff >= req(input$"Cortex_LRpairs_S_inter_diff")]
    }, ignoreNULL = FALSE)
    
    # --------------------------------------------------
    # Description: Render a heatmap plot using the 'scSeqComm_heatmaply' function for the filtered ligand-receptor interaction data
    # Output: A plotly heatmap plot displayed in the Shiny app UI
    # --------------------------------------------------
    
    output$Adult_LRpairs <- renderPlotly({

      input$"Cortex_LRsetting_button"

      # Render a heatmap plot using the 'scSeqComm_heatmaply' function for the filtered ligand-receptor interaction data
      scSeqComm_heatmaply(req(Adult_LRpairs_filter()), x=isolate(input$"Cortex_LRpairs_x.axis"), y=isolate(input$"Cortex_LRpairs_y.axis"), 
                          fill = isolate(input$"Cortex_LRpairs_fill"), size = isolate(input$"Cortex_LRpairs_pt.size"),  
                          height = isolate(input$"LRpairs_height"),width= isolate(input$"LRpairs_width"),
                          gap_width= isolate(input$"LRpairs_step_width"), gap_height= isolate(input$"LRpairs_step_height"),
                          xlab="", ylab = "",  gdt_pal=isolate(input$"Cortex_LRpairs_pal"),
                          limit_fill = c(0, 1), limit_size = c(0, 1), breaks_size = c(0.2, 
                                                                                      0.5, 0.8), 
                          facet_grid_x = isolate(input$"Cortex_LRpairs_facetgrid.x"), facet_grid_y = isolate(input$"Cortex_LRpairs_facetgrid.y"))
      
    })
    
      # --------------------------------------------------
      # Description: create a text output that displays information about the user's selection
      # --------------------------------------------------

    output$Adult_selectLR <- renderText({
      
      paste("Your selection contains", nrow(LRDB_mouse_LRpairs()), "ligand-receptor pairs", "|", 
            length(unique(Adult_LRpairs_filter()$LR_pair[which(Adult_LRpairs_filter()$cluster_L %in% req(input$"Cortex_LRpairs_cluster_L")  & 
                                                              Adult_LRpairs_filter()$cluster_R %in% req(input$"Cortex_LRpairs_cluster_R")  &  
                                                              Adult_LRpairs_filter()$SOURCE.class %in% req(input$"Cortex_LRpairs_SOURCE.class")  &  
                                                              Adult_LRpairs_filter()$TARGET.class %in% req(input$"Cortex_LRpairs_TARGET.class") )])), 
            "expressed in the Adult dataset with the current filters")

    })
    
  })
  
})
        
   
