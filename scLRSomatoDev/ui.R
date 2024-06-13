#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


source("utils.R", local = TRUE)



# Define UI for application
shinyUI(

  navbarPage(

    # Define the title of the application and the theme
    title = "scLRSomatoDev", id = "navbar",
    theme = bs_theme(bootswatch = "flatly", version = 5),

    # Define the loading page
    tabPanel(title = "Loading", value = "Loading",
      shinycssloaders::withSpinner(
        uiOutput("app_content"),
        color = "#02E7B9",
        type = 8,
        size = 3,
        caption = div(
          strong("Loading data"), br(), "Please wait about 2 min...",
          style = "font-size:30px;"
        )
      )
    ),
    ###############
    # Overview Tab
    ###############

    tabPanel(
      title = "Overview",

      ########################
      ## Overview Page to present the scLRSomatoDev Shiny App
      ## this page give a short description of the app functionality
      ########################

      ## Header to introduce the  author and the logo
      layout_columns(

        ## set column widths
        col_widths = c(1, 4, 6, 1),

        ## Logo + Project description

        div(),

        div(
          style = "text-align: left; 
          align-content: center;
          align-items: center;",
          
          br(), br(), br(),
          
          ## Welcome message
          h1("Welcome to scLRSomatoDev", style = "text-align: center;"),

          ## Project description
          p("This Shiny App has been made to visualize the ligand-receptor interactions between GABAergic and Glutamatergic cells during Somatosensory Cortex development."),
          p(" You can navigate across the shiny App through the navigation bar Menu located at the top of the window."),
          p("Created and maintained by Rémi MATHIEU (INMED, INSERM, Aix Marseille Univ, France)."),
          p("UI/UX Reworked by Lucas SILVAGNOLI (INMED, INSERM, Aix Marseille Univ, France)."),
          p("last update: 2024-05"),
        ),

        ## Shiny App tab showcase images
        div(
          fluidPage(
            div(img(src = "logo_scLRSomatoDev.png", width = 950, height = 425),
                style = "text-align: center; 
                align-content: center;align-items: center;"))
        ),

        div(),
      ),

      br(),

      ## Value boxes to present some statistics of the shiny app
      layout_columns(

        col_widths = c(2, 4, 4, 2),

        div(),

        ## number of ligand-receptor pairs
        div(
          h4("8789 Ligand-Receptor Pairs",
             style = "text-align: center; 
             align-content: center;align-items: center;"),
        ),

        ## number of time points
        div(
          h4("17 Time Points", 
             style = "text-align: center; 
             align-content: center;align-items: center;"),
        ),

        div(),
      ),

      br(),

      ########################
      ## Gene expression menu
      ##provide descriptive information on the gene expression menu tab
      ## I have to abuse the columns layout to divide content and align properly
      ########################

      div(style = "align-content: center;align-items: center; 
        background-color: #2f3e50; color: #ffffff;",

        ## section title
        div(style = "text-align: center;",

          br(), br(),

          h3("Gene Expression features"),
          p("In each section, you will find an help tab panel 
          that will help you to interact with the different options."),
        ),

        br(), br(),

        ## metadata table Visualization container (image on the right + text)
        layout_columns(

          ## set column widths to divide content
          col_widths = c(2, 4, 4, 2),

          ## left margin
          div(),

          ## Text container
          div(
            br(), br(), br(),
            strong("Metadata table"),
            p("Access and review the metadata associated with the 17 datasets used in this study."),
            style = "text-align: left; padding-left: 20%;",
          ),

          ## Image container
          div(
            fluidPage(
              div(img(src = "metadata_table.png", width = 344, height = 250),
                  style = "text-align: center;")
            )
          ),

          ## right margin
          div(),
        ),

        ## Clustering results visualization container (image on the left + text)
        layout_columns(

          ## set column widths to divide content
          col_widths = c(2, 4, 4, 2),

          ## left margin
          div(),

          ## Image container
          div(
            fluidPage(
              div(img(src = "Clustering_results_L.png",
                      width = 800, height = 250),
                  style = "text-align: center; padding-left: 10%; position: relative; right: 55%;")
            )
          ),

          ## Text container
          div(
            br(), br(), br(),
            strong("Clustering Results and feature plots"),
            p("Visualize the distinct cell-type populations in the dataset using UMAP 
            (or alternative dimensionality reduction methods) and explore the expression profiles of one or two genes."),
            style = "text-align: left; padding-right: 20%;"
          ),

          ## right margin
          div(),
        ),

        ## Gene expression Visualization container (image on the right + text)
        layout_columns(

          ## set column widths to divide content
          col_widths = c(2, 4, 4, 2),

          ## left margin
          div(),

          ## Text container
          div(
            br(), br(), br(),
            strong("Absolute expression"),
            p("Visualize the expression profiles of multiple genes per cell-type across cortical development using a heatmap or Dot plot."),
            style = "text-align: left; padding-left: 20%;"
          ),

          ## Image container
          div(
            fluidPage(
              div(img(src = "Gene_expression_L.png",
                      width = 800, height = 250),
                style = "text-align: center; padding-left: 10%;"
              )
            )
          ),

          ## right margin
          div(),
        ),

        ## Pseudo-maturation visualization container (image on the left + text)
        layout_columns(

          ## set column widths to divide content
          col_widths = c(2, 4, 4, 2),

          ## left margin
          div(),

          ## Image container
          div(
            fluidPage(
              div(img(src = "Pseudo_maturation_D.png",
                      width = 344, height = 250),
                  style = "text-align: center; padding-left: 10%;")
            )
          ),

          ## Text container
          div(
            br(), br(), br(),
            strong("Pseudo-maturation"),
            p("Visualize the dynamics of gene expression along the pseudo-maturation axis."),
            style = "text-align: left; padding-right: 10%;"
          ),

          ## right margin
          div(),
        ),

        ## Pseudo-layer Visualization container (image on the right + text)
        layout_columns(

          ## set column widths to divide content
          col_widths = c(2, 4, 4, 2),

          ## left margin
          div(),

          ## Text container
          div(
            br(), br(), br(),
            strong("Pseudo-layer"),
            p("Visualize the dynamics of gene expression along the pseudo-layer axis."),
            style = "text-align: left; padding-left: 20%;"
          ),

          ## Image container
          div(
            fluidPage(
              div(img(src = "Pseudo_layer_D.png",
                      width = 344, height = 250),
                  style = "text-align: center; padding-left: 10%;")
            )
          ),

          ## right margin
          div(),
        ),

        ## Transcriptional landscape visualization (image on the left + text)
        layout_columns(

          ## set column widths to divide content
          col_widths = c(2, 4, 4, 2),

          ## left margin
          div(),

          ## Image container
          div(
            fluidPage(
              div(img(src = "Transcriptional_landscape_D.png",
                      width = 344, height = 250),
                  style = "text-align: center; padding-left: 10%;")
            )
          ),

          ## Text container
          div(
            br(), br(), br(),
            strong("Transcriptional landscape"),
            p("Visualize the dynamics of gene expression along both the pseudo-maturation and pseudo-layer axes, represented as a 2D map."),
            style = "text-align: left; padding-right: 10%;"
          ),

          ## right margin
          div(),
        ),

        br(), br(),

      ),


      ########################
      ## Ligand-Receptor menu
      ## provide descriptive information on the gene expression menu tab
      ## I have to abuse the columns layout to divide content and align properly
      ########################

      br(), br(),

      h3("Ligand-Receptor features", style = "text-align: center;"),
      p("In each section, you will find an help tab panel 
      that will help you to interact with the different options.", 
        style = "text-align: center;"),

      br(), br(),

      div(

        br(),

        ## LRintercellNetworkDB Visualization container 
        ##(image on the right + text)
        layout_columns(

          ## set column widths to divide content
          col_widths = c(2, 4, 4, 2),

          ## left margin
          div(),

          ## Text container
          div(
            br(), br(), br(), 
            strong("LRintercellNetworkDB"),
            p("Explore the curated database of ligand-receptor pairs."),
            style = "text-align: left;align-content: 
            center;align-items: center; padding-left: 20%;"
          ),

          ## Image container
          div(
            fluidPage(
              div(img(src = "LRintercellNetworkDB.png",
                      width = 344, height = 250),
                  style = "text-align: center;align-content: center;
                align-items: center; padding-right: 10%;")
            )
          ),

          ## right margin
          div(),
        ),

        br(),

        ## Table visualization container (image on the left + text)
        layout_columns(

          ## set column widths to divide content
          col_widths = c(2, 4, 4, 2),

          ## left margin
          div(),

          ## Image container
          div(
            fluidPage(
              div(img(src = "LR_table.png",
                      width = 344, height = 250),
                style = "text-align: center;align-content: 
                center;align-items: center; padding-left: 10%;"
              )
            )
          ),

          ## Text container
          div(
            br(), br(), br(),
            strong("LR Table"),
            p("Access the tables resulting from intercellular and intracellular signaling analyses performed with scSeqComm (Baruzzo et al., 2022) 
            for each developmental age, .ie., E18.5-P0, P1-P2, P4-P5, P8, P16, P30, and Adult."),
            style = "text-align: left;align-content: center;
            align-items: center; padding-right: 20%;"
          ),

          ## right margin
          div(),
        ),

        br(),

        ## Number of interaction Visualization container 
        ## (image on the right + text)
        layout_columns(

          ## set column widths to divide content
          col_widths = c(2, 4, 4, 2),

          ## left margin
          div(),

          ## Text container
          div(
            br(), br(), br(),
            strong("Number of interactions"),
            p("Visualize the number of predicted interactions (ligand-receptor pairs) between each cell-type pair for each developmental age, .i.e. , E18.5-P0, P1-P2, P4-P5, P8, P16, P30, and Adult."),
            style = "text-align: left;align-content: 
            center;align-items: center; padding-left: 20%;"
          ),

          ## Image container
          div(
            fluidPage(
              div(img(src = "Number_of_interaction.png",
                      width = 344, height = 250),
                style = "text-align: center;align-content: center;
                align-items: center; padding-right: 10%;"
              )
            )
          ),

          ## right margin
          div(),
        ),

        br(),

        ## Intercellular/Intracellular signaling visualization container 
        ## (image on the left + text)
        layout_columns(

          ## set column widths to divide content
          col_widths = c(2, 4, 4, 2),

          ## left margin
          div(),

          ## Image container
          div(
            fluidPage(
              div(img(src = "Intercellular_Intracellular_signaling_D.png",
                      width = 344, height = 250),
                  style = "text-align: center;align-content: 
                  center;align-items: center; padding-left: 10%;")
            )
          ),

          ## Text container
          div(
            br(), br(), br(),
            strong("Intercellular/Intracellular signaling"),
            p("Visualize which ligand-receptor pairs are likely to be present between cell-type pairs and the pathways in which they are involved."),
            style = "text-align: left;align-content: 
            center;align-items: center; padding-right: 20%;"
          ),

          ## right margin
          div(),
        ),

        br(), br(), br(),

      ),

    ),

    ###############
    #Datasets Tab
    ###############

    navbarMenu(title = "Gene Expression",
      ######### Metadata table

      # a Page to display the metadata table
      tabPanel(title = "Metadata table",
        fluidPage(
          tabsetPanel(
            tabPanel("Metadata table",
                     dataTableOutput(outputId = "Cortex_metadata")),
            tabPanel("Help",
              htmltools::includeMarkdown("Data/Help_datasets_metadata.Rmd")
            )
          )
        )
      ),
      ############ Data Visualization

      # a Page to display the data visualization (UMAP and tSNE)
      tabPanel(title = "Clustering results and Feature plots",
        fluidPage(
          titlePanel("Visualize clustering results and feature plots"),
          br(),
          sidebarLayout(

            # Sidebar panel for User selection inputs
            sidebarPanel(

              ## sidebar parameters
              width=3,

              ## Title
              h4("Data visualization Explorer:"),

              ## Dimension reduction Method selection
              pickerInput(inputId = "Dim_red", label= "Dimension reduction:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F, width="100%"),

              ## Grouping level selection
              pickerInput(inputId = "DataVis_level", label= "Grouping level:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F, width="100%"),

              ## Gene selection
              pickerInput(inputId = "DataVis_gene", label= "Select 1 or 2 genes:", choices= NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE, "max-options" = 2,"max-options-text" = "No more selectable genes"), multiple=T, width="100%"),

              ## Age selection
              pickerInput(inputId = "DataVis_age", label= "Age:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

              hr(),

              accordion(

                accordion_panel("Grouping Level filters",

                  ## class selection
                  pickerInput(inputId = "DataVis_class", label= "Class:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                  ## family selection
                  pickerInput(inputId = "DataVis_family", label= "Family:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                  ## subclass selection
                  pickerInput(inputId = "DataVis_subclass", label= "Subclass:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                  ## Supertype selection
                  pickerInput(inputId = "DataVis_supertype", label= "Supertype:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                  ## cell-type selection
                  pickerInput(inputId = "DataVis_celltype", label= "Cell-type:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                  ## Original cell-type selection
                  pickerInput(inputId = "DataVis_celltype_original", label= "Original cell-type:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                ),

                accordion_panel("Other filters",

                  ## Region selection
                  pickerInput(inputId = "DataVis_region", label= "Region:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple=T, width="100%"),

                  ## Study selection
                  pickerInput(inputId = "DataVis_study", label= "Study:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                  ## Platform selection
                  pickerInput(inputId = "DataVis_platform", label= "Platform:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                  ## RNA-seq method selection
                  pickerInput(inputId = "DataVis_RNAseq.method", label= "RNA-seq method:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                ),

                accordion_panel("Plot options",

                  ## palette selection for gene1
                  palettePicker(inputId = "DataVis_pal1",
                    label = "Select a color palette for gene1:",
                    choices = list(
                      Reds = brewer_pal(palette = "Reds")(8),
                      Greens = brewer_pal(palette = "Greens")(8),
                      Blues = brewer_pal(palette = "Blues")(8),
                      Purples = brewer_pal(palette = "Purples")(8),
                      Oranges = brewer_pal(palette = "Oranges")(8),
                      OrRd = brewer_pal(palette = "OrRd")(8),
                      YlOrRd = brewer_pal(palette = "YlOrRd")(8),
                      YlOrBr = brewer_pal(palette = "YlOrBr")(8),
                      BuYlRd = brewer_pal(palette = "RdYlBu", direction = -1)(8)
                    ),
                    textColor = rep("black", 9),
                    selected = "Reds"
                  ),

                  ## palette selection for gene2
                  palettePicker(inputId = "DataVis_pal2",
                    label = "Select a color palette for gene2:",
                    choices = list(
                      Reds = brewer_pal(palette = "Reds")(8),
                      Greens = brewer_pal(palette = "Greens")(8),
                      Blues = brewer_pal(palette = "Blues")(8),
                      Purples = brewer_pal(palette = "Purples")(8),
                      Oranges = brewer_pal(palette = "Oranges")(8),
                      OrRd = brewer_pal(palette = "OrRd")(8),
                      YlOrRd = brewer_pal(palette = "YlOrRd")(8),
                      YlOrBr = brewer_pal(palette = "YlOrBr")(8),
                      BuYlRd = brewer_pal(palette = "RdYlBu", direction = -1)(8)
                    ),
                    textColor = rep("black", 9),
                    selected = "Greens"
                  ),

                  ## Dimension reduction plot settings
                  sliderInput(inputId = "DataVis_height", 
                              label = "Set the Dimension reduction plot height:",
                              min = 400, max = 2000, value = 750, step = 50),
                )
              ),

              hr(),

              ## Apply selection button
            fluidPage(
              layout_columns(
                col.widths = c(2, 8, 2),
                div(),
                div(
                style = "align-content: center; align-self: center;",
                actionButton(inputId = "Dimred_button", label = "Apply Selection"),
                ),
                div(),
              )
              
            )
          ),

            ## Tabs Panel for dimension reduction
            mainPanel(
              tabsetPanel(

                ## Tab for grouping level
                tabPanel("Grouping level",

                  ## Render plot
                  shinycssloaders::withSpinner(
                  plotlyOutput(outputId = 'Cortex_Dim.red', height="850px",width="100%"),
                    color = "#02E7B9",
                    type = 3,
                    color.background = "white",
                    size = 1.5,
                    caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)", 
                    style = "font-size:20px;")
                  )
                ),

                ## Tab for gene expression
                tabPanel("Feature plots", br(),

                  ## Render plot
                  shinycssloaders::withSpinner(
                  plotlyOutput(outputId = 'Cortex_Dim.red_ExprGene', height="850px", width="100%"),
                  color = "#02E7B9",
                  type = 3,
                  color.background = "white",
                  size = 1.5,
                  caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                style = "font-size:20px;")
                  )
                ),

                ## Help tab
                tabPanel("Help", htmltools::includeMarkdown("Data/Help_datasets_Data visualization.Rmd"))
              ), width = 9
            )
          )
        )
      ),

      ####################### Gene Expression

      # A Page to display the gene expression
      tabPanel(title = "Absolute Expression",
        fluidPage(
          titlePanel("Visualize Gene Expression"),
          br(),
          tabsetPanel(

            ## Tab for grouped gene expression heatmap and dot plot
            tabPanel("Gene Expression by grouping level",
              sidebarLayout(

                ## Sidebar panel for user input selection
                sidebarPanel(

                  ## Sidebar Parameters
                  width= 3,

                  ## Title
                  h4("Gene Expression Explorer:"),

                  ## Plot type selection
                  pickerInput(inputId = "ExprData_plot", label= "Plot type:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F, width="100%"),

                  ## Grouping level selection
                  pickerInput(inputId = "ExprData_level", label= "Grouping level:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F, width="100%"),

                  ## Gene selection
                  pickerInput(inputId = "ExprData_gene", label= "Select  genes:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE, "max-options"=500, "max-options-text" = "No more selectable genes"), multiple=T, width="100%"),

                  ## Age selection
                  pickerInput(inputId = "ExprData_age", label= "Age:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                  hr(),

                  accordion(

                    accordion_panel( "Grouping Level filters",

                      ## Class selection
                      pickerInput(inputId = "ExprData_class", label= "Class:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## Family selection
                      pickerInput(inputId = "ExprData_family", label= "Family:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## Subclass selection
                      pickerInput(inputId = "ExprData_subclass", label= "Subclass:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## Cell-type selection
                      pickerInput(inputId = "ExprData_celltype", label= "Cell-type:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## Original cell-type selection
                      pickerInput(inputId = "ExprData_celltype_original", label= "Original cell-type:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),
                    ),

                    accordion_panel( "Other filters",

                      ## Study selection
                      pickerInput(inputId = "ExprData_study", label= "Study:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## Region selection
                      pickerInput(inputId = "ExprData_region", label= "Region:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## RNA-seq method selection
                      pickerInput(inputId = "ExprData_RNAseq.method", label= "RNA-seq method:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## Platform selection
                      pickerInput(inputId = "ExprData_platform", label= "Platform:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),
                                        
                    ),

                    accordion_panel("Plot options",

                      ## Palette selection
                      palettePicker(inputId = "ExprData_pal", 
                                    label = "Select a color palette:", 
                                    choices = list(
                                      PuOr = brewer_pal(palette = "PuOr", direction=-1)(11),
                                      GBBr = brewer_pal(palette = "BrBG", direction=-1)(11),
                                      BuYlRd = brewer_pal(palette = "RdYlBu", direction = -1)(11),
                                      BuPhRd = c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                                      GnPhRd = c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                                      BuGyOrRd = c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500"),
                                      BuGyRd= c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                                      Rainbow = rainbow(5,rev=T),
                                      Heat = heat.colors(10, rev=T),
                                      viridis = viridis_pal(option = "viridis")(10),
                                      magma = viridis_pal(option = "magma")(10),
                                      inferno = viridis_pal(option = "inferno")(10),
                                      plasma = viridis_pal(option = "plasma")(10),
                                      cividis = viridis_pal(option = "cividis")(10)
                                    ),
                                    textColor = c(rep("white", 8), "black" ,rep("white", 5) ),selected = "BuGyOrRd"),

                      ## Plot height selection   
                      sliderInput(inputId = "ExprData_height", label = "Set the  plot height:", min = 400, max = 2000, value = 600, step = 50),
                    )
                  ),

                  hr(),

                  ## Apply selection button 
                  
                  fluidPage(
                    layout_columns(
                      col.widths = c(2, 8, 2),
                      div(),
                      div(
                        style = "align-content: center; align-self: center;",
                        actionButton(inputId = "ExprData_button", label = "Apply Selection"),
                      ),
                      div(),
                    )
                  )
                ),

                ## Plot area
                mainPanel((div(style='width:100%;height:800px; overflow-y: scroll; position: relative',

                 ## Render the plot (with Loading Animation)
                  shinycssloaders::withSpinner(
                    plotlyOutput(outputId = "Cortex_ExprData"),
                    color = "#02E7B9",
                    type = 3,
                    color.background = "white",
                    size = 1.5,
                    caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                  style = "font-size:20px;")
                  ), ) ),
                textOutput("ExprData_selectLR")
                )
              )
            ),

            ## tab for gene expression by age
            tabPanel("Gene Expression by grouping level and by age",
              sidebarLayout(

                ## sidebar panel for user input selection
                sidebarPanel(

                  ##Sidebar parameters
                  width = 3,

                  ## Title
                  h4("Gene Expression Explorer:"),

                  ## Grouping level selection
                  pickerInput(inputId = "HmctAge_level", label= "Grouping level:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F, width="100%"),

                  ## Gene selection
                  pickerInput(inputId = "HmctAge_gene", label= "Select  genes:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F, width="100%"),

                  hr(),

                  accordion(

                    accordion_panel("Grouping Level filters",

                      ## Class selection
                      pickerInput(inputId = "HmctAge_class", label= "Class:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## Family selection
                      pickerInput(inputId = "HmctAge_family", label= "Family:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## Subclass selection
                      pickerInput(inputId = "HmctAge_subclass", label= "Subclass:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## Cell-type selection
                      pickerInput(inputId = "HmctAge_celltype", label= "Cell-type:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## Original cell-type selection
                      pickerInput(inputId = "HmctAge_age", label= "Age:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                    ),

                    accordion_panel("Other filters",

                      ## Study selection
                      pickerInput(inputId = "HmctAge_study", label= "Study:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## Region selection
                      pickerInput(inputId = "HmctAge_region", label= "Region:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## RNA-seq method selection
                      pickerInput(inputId = "HmctAge_RNAseq.method", label= "RNA-seq method:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                      ## Platform selection
                      pickerInput(inputId = "HmctAge_platform", label= "Platform:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

                    ),

                    accordion_panel("Plot options",

                      ## Color palette selection
                      palettePicker(
                       inputId = "HmctAge_pal", label = "Select a color palette:", 
                        choices = list(
                          PuOr = brewer_pal(palette = "PuOr", direction=-1)(11),
                          GBBr = brewer_pal(palette = "BrBG", direction=-1)(11),
                          BuYlRd = brewer_pal(palette = "RdYlBu", direction = -1)(11),
                          BuPhRd = c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                          GnPhRd = c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                          BuGyOrRd = c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500"),
                          BuGyRd= c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                          Rainbow = rainbow(5,rev=T),
                          Heat = heat.colors(10, rev=T),
                          viridis = viridis_pal(option = "viridis")(10),
                          magma = viridis_pal(option = "magma")(10),
                          inferno = viridis_pal(option = "inferno")(10),
                          plasma = viridis_pal(option = "plasma")(10),
                          cividis = viridis_pal(option = "cividis")(10)
                        ),
                        textColor = c(rep("white", 8), "black" ,rep("white", 5) ),
                        selected = "BuGyOrRd"
                      ),

                      ##  Plot height selection
                      sliderInput(inputId = "HmctAge_height", label = "Set the  plot height:", min = 400, max = 2000, value = 600, step = 50), 

                    )
                  ),

                  hr(),

                  ## Apply selection button
                  fluidPage(
                    layout_columns(
                      col.widths = c(2, 8, 2),
                      div(),
                      div(
                        style = "align-content: center; align-self: center;",
                        actionButton(inputId = "HmctAge_button", label = "Apply Selections"),
                      ),
                      div(),
                    )
                  )

                ),

                ## Plot Area
                mainPanel((div(style='width:100%;height:850px; overflow-y: scroll; position: relative',

                ## Render the Heatmap plot by Age (with Loading Animation)
                  shinycssloaders::withSpinner(
                    plotlyOutput(outputId = "Cortex_HmctAge"),
                    color = "#02E7B9",
                    type = 3,
                    color.background = "white",
                    size = 1.5,
                    caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                  style = "font-size:20px;")
                  )
                )
                ),
                textOutput("HmctAge_selectGene")
                )
              )
            ),

            ## Help tab
            tabPanel("Help", htmltools::includeMarkdown("Data/Help_datasets_Gene expression.Rmd"))

          )
        )
      ),


      ##################### Pseudo-maturation

      ## A panel for Pseudo-maturation plots
      tabPanel(title = "Pseudo-maturation", 
      fluidPage(
        titlePanel("Visualize Gene Expression along the pseudo-maturation axis"),
        br(),
        sidebarLayout(

          ## sidebar panel for user input selection
          sidebarPanel(

            ## sidebar Parameters
            width = 3,

            ## Title
            h4("Gene Expression Explorer:"),

            ## Cell-type selection
            pickerInput(inputId = "pM_celltype", label= "Cell-type:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

            ## Gene selection
            pickerInput(inputId = "pM_gene", label= "Select genes:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE, "max-options"=50, "max-options-text" = "No more selectable genes"), multiple=T, width="100%"),

            accordion(

              accordion_panel("Plot options",

                ## Color palette selection
                palettePicker(inputId = "pM_pal", label = "Select a color palette:", 
                  choices = list(
                    PuOr = brewer_pal(palette = "PuOr", direction=-1)(11),
                    GBBr = brewer_pal(palette = "BrBG", direction=-1)(11),
                    BuYlRd = brewer_pal(palette = "RdYlBu", direction = -1)(11),
                    BuRd = brewer_pal(palette = "RdBu", direction = -1)(11),
                    BuPhRd = c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                    GnPhRd = c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                    BuGyOrRd = c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500"),
                    BuGyRd= c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                    BuSkyWhPkRd=c("#0000CD", "#87CEEB","#FFFFFF", "#FFC0CB", "#CD0000"),
                    Vik=scico(11, palette="vik")
                  ),

                  textColor = c(rep("white", 10)),

                  selected = "BuSkyWhPkRd"
                ),

                ##  Plot height and width selection
                sliderInput(inputId = "pM_height", label = "Set the plot height:", min = 400, max = 2000, value = 1000, step = 50),

                sliderInput(inputId = "pM_width", label = "Set the plot width:", min = 400, max = 2000, value = 1000, step = 50),

              )
            ),

            hr(),

            ## Apply selection button
            fluidPage(
                    layout_columns(
                      col.widths = c(2, 8, 2),
                      div(),
                      div(
                        style = "align-content: center; align-self: center;",
                        actionButton(inputId = "pM_button", label = "Apply Selection"),
                      ),
                      div(),
                    )
                  )

          ),

          mainPanel(
            tabsetPanel(

              ## Heatmap Plot Tab
              tabPanel("Heatmap",(div

              ## Set the plot display area
              (style='width:100%;height:850px; overflow-y: scroll; position: relative', 

                ## Render the Heatmap plot (with Loading Animation)
                shinycssloaders::withSpinner(
                  plotOutput(outputId = "Hm_pM"),
                  color = "#02E7B9", 
                  type = 3, 
                  color.background = "white", 
                  size = 1.5,
                  caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                style = "font-size:20px;")
                )
              )
              ), 

              ## Render the text output
              textOutput("pM_selectGene_Hm")
              ),

              ## Spline Plot Tab
              tabPanel("Curve",(div(

                ## Set the plot display area
                style='width:100%;height:850px; overflow-y: scroll; position: relative', 

                ## Render the Spline plot (with Loading Animation)
                shinycssloaders::withSpinner(
                  plotlyOutput(outputId = "Spline_pM"),
                  color = "#02E7B9",
                  type = 3,
                  color.background = "white",
                  size = 1.5,
                  caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                style = "font-size:20px;")
                )
              )
            ),

              ## Render the text output
              textOutput("pM_selectGene_Spline")),

              ## Help tab
              tabPanel("Help", htmltools::includeMarkdown("Data/Help_datasets_Pseudo-maturation.Rmd"))
            )
          )
        )
      )
    ),

      ##################### Pseudo-layer

      ## A Page for Pseudo-layer plots
      tabPanel(title="Pseudo-layer", 
        fluidPage(
          titlePanel("Visualize Gene Expression along Pseudo-layer axis"),
          br(),
          sidebarLayout(
            
            ## sidebar panel for user input selection
            sidebarPanel(

            ## Sidebar parameters
            width=3,
                        
            ## Cell-type selection
            pickerInput(inputId = "pL_celltype", label= "Family:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F, width="100%"),

            ## Age selection
            pickerInput(inputId = "pL_age", label= "Age:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE, "max-options"=50, "max-options-text" = "No more selectable genes"), multiple=T, width="100%"),

            ## Genes selection
            pickerInput(inputId = "pL_gene", label= "Select genes:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T, width="100%"),

            accordion(

              accordion_panel("Plot options",

                ## Color palette selection for heatmap plot
                palettePicker( inputId = "pL_pal", label = "Select a color palette for the heatmap:", 
                  choices = list(
                    PuOr = brewer_pal(palette = "PuOr", direction=-1)(11),
                    GBBr = brewer_pal(palette = "BrBG", direction=-1)(11),
                    BuYlRd = brewer_pal(palette = "RdYlBu", direction = -1)(11),
                    BuRd = brewer_pal(palette = "RdBu", direction = -1)(11),
                    BuPhRd = c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                    GnPhRd = c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                    BuGyOrRd = c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500"),
                    BuGyRd= c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                    BuSkyWhPkRd=c("#0000CD", "#87CEEB","#FFFFFF", "#FFC0CB", "#CD0000"),
                    Vik=scico(11, palette="vik")
                    ),

                  textColor = c(rep("white", 10)),
                            
                  selected = "BuSkyWhPkRd"
                ),

                ## Color palette selection for spline plot
                palettePicker( inputId = "pL_spline_pal", label = "Select a color palette for Spline:", 
                  choices = list(
                    Blues = blues9,
                    Reds = brewer_pal(palette="Reds")(9)
                    ),

                  textColor = c(rep("black", 2)),
                            
                  selected = "Blues"
                ),

                ## Plot height and width selection
                sliderInput(inputId = "pL_height", label = "Set the  plot height:", min = 400, max = 2000, value = 750, step = 50),

                sliderInput(inputId = "pL_width", label = "Set the plot width:", min = 400, max = 2000, value = 1000, step = 50),
              )
            ),    

            hr(),

            ## Apply selection button
            fluidPage(
                    layout_columns(
                      col.widths = c(2, 8, 2),
                      div(),
                      div(
                        style = "align-content: center; align-self: center;",
                        actionButton(inputId = "pL_button", label = "Apply Selection"),
                      ),
                      div(),
                    )
                  )
            ),
                    
                    mainPanel(

                        tabsetPanel(
                        
                            ## Heatmap tab
                            tabPanel("Heatmap",(div
                        
                                ## set the plot display area
                                (style='width:100%;height:850px; overflow-y: scroll; position: relative', 
                            
                                 ## Render the heatmap plot (with Loading Animation)
                                shinycssloaders::withSpinner(
                                    plotOutput(outputId = "Hm_pL"),
                                    color = "#02E7B9", 
                                    type = 3, 
                                    color.background = "white", 
                                    size = 1.5,
                                caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                              style = "font-size:20px;")
                                    )
                                )
                            ), 
                            
                            ## Render the text output
                            textOutput("pL_selectGene_Hm")),

                            ## Spline Plot Tab
                            tabPanel("Curve",(div
                        
                            ## Set the plot display area
                            (style='width:100%;height:850px; overflow-y: scroll; position: relative', 
                            
                            ## Render the Spline plot (with Loading Animation)
                            shinycssloaders::withSpinner(
                                plotlyOutput(outputId = "Spline_pL"),
                                color = "#02E7B9", 
                                type = 3, 
                                color.background = "white", 
                                size = 1.5,
                                caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                              style = "font-size:20px;")
                                    )
                                )
                            ), 
                            
                            ## Render the text output
                            textOutput("pL_selectGene_Spline")),

                            ## Help tab
                            tabPanel("Help", htmltools::includeMarkdown("Data/Help_datasets_Pseudo-layer.Rmd"))

                                )
                            )
                        )     
                    )
      ),

      ##################### Transcriptional landscape

      ## A page to display the transcriptional landscape
      tabPanel(title="Transcriptional landscape", 
        fluidPage(
          titlePanel("Visualize Gene transcriptional landscape"),
          br(),
                tabsetPanel(tabPanel("Transcriptional landscape",
                    sidebarLayout(
                        
                        ## sidebar panel for user input selection
                        sidebarPanel(

                            ## Sidebar parameters
                            width=3,

                            ## Title
                            h4("Gene Expression Explorer:"),
                            
                            ## Family selection
                            pickerInput(inputId = "TL_celltype", label= "Family:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F, width="100%"),

                            ## Genes selection
                            pickerInput(inputId = "TL_gene", label= "Select  genes:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE, "max-options"=50, "max-options-text" = "No more selectable genes"), multiple=T, width="100%"),

                            accordion(
                                accordion_panel("Plot options",

                                ## Color palette selection
                                    palettePicker(inputId = "TL_pal", label = "Select a color palette:", 
                                        choices = list(
                                            PuOr = brewer_pal(palette = "PuOr", direction=-1)(11),
                                            GBBr = brewer_pal(palette = "BrBG", direction=-1)(11),
                                            BuYlRd = brewer_pal(palette = "RdYlBu", direction = -1)(11),
                                            BuRd = brewer_pal(palette = "RdBu", direction = -1)(11),
                                            BuPhRd = c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                                            GnPhRd = c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                                            BuGyOrRd = c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500"),
                                            BuGyRd= c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                                            BuSkyWhPkRd=c("#0000CD", "#87CEEB","#FFFFFF", "#FFC0CB", "#CD0000"),
                                            Vik=scico(11, palette="vik")
                                            ),
                                
                                    textColor = c(rep("white", 10)),
                                
                                    selected = "BuRd"
                                        ),
                            
                                ## Plot height and width selection
                            
                                sliderInput(inputId = "TL_height", label = "Set the plot height:", min = 400, max = 2000, value = 750, step = 50),

                                sliderInput(inputId = "TL_width", label = "Set the plot width:", min = 400, max = 2000, value = 1000, step = 50),

                                )
                            ),

                            hr(),

                            ## Apply selection button
                            fluidPage(
                              layout_columns(
                                col.widths = c(2, 8, 2),
                                div(),
                                div(
                                  style = "align-content: center; align-self: center;",
                                  actionButton(inputId = "TL_button", label = "Apply selection"),
                                  ),
                                div(),
                              )
                            )
                        ),
                    
                    ## Plot area
                    mainPanel((div
                    
                        ## Plot area
                        (style='width:100%;height:850px; overflow-y: scroll; position: relative', 
                    
                        ## Render the Heatmap landscape plots (with Loading Animation)
                        shinycssloaders::withSpinner(
                        plotOutput(outputId = 'Cortex_TL'),
                        color = "#02E7B9", 
                        type = 3, 
                        color.background = "white", 
                        size = 1.5,
                        caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)", 
                                      style = "font-size:20px;")
                            )
                        )
                    ),
                        
                        ## Render the text output
                        textOutput("TL_selectGene"))
                                )
                             ),
                    ## Help tab
                    tabPanel("Help", htmltools::includeMarkdown("Data/Help_datasets_transcriptional landscape.Rmd"))
                        )      
                    )
      )             
    ),
    

############################
# LR Significant Data
###########################

    ## A page to display the LR data
    navbarMenu(title = "Ligand-Receptor",
               
        ############ LRintercellNetworkDB ######

        ## A page to display the LRintercellNetworkDB
        tabPanel("LRintercellNetworkDB",
            tabsetPanel(

                ## tab to show how the LRintercellNetworkDB is constructed
                tabPanel("Database construction",fluidPage(div(img(src="Fig LRintercellNetworkDB.png", width=1782, height=1347, align="center"), style="text-align: center;"))),
                
                ## tab to show a table of the LRintercellNetworkDB data
                tabPanel("Table", fluidPage(
                    
                    ## Render the table
                    dataTableOutput(outputId = 'table_LRDB'))),

                ## tab to show an image of the ligand-receptor location
                tabPanel("Ligand & Receptor location", fluidPage(div(img(src="LRintercellNetworkDB_ProtLoc.png", width=1782, height=1347, align="center"), style="text-align: center;"))),

                ## tab to show an image of the ligand-receptor category
                tabPanel("Ligand & Receptor category", fluidPage(div(img(src="Ligand_Recptor_category_LRintercellNetworkDB.png", width=1782, height=1347, align="center"), style="text-align: center;"))),

                ## tab to show a plot of the ligand-receptor families
                tabPanel("Ligand & Receptor family", varSelectInput(inputId = "LRfamily_var", label= "Select Ligand or Receptor family:", data=LRintercellNetworkDB[,7:8]),
                     (div
                     
                        ## Plot area
                        (style='width:100%;height:950px; overflow-y: scroll; position: relative',
                        
                        ## Render the plot
                        shinycssloaders::withSpinner(
                        plotlyOutput(outputId = "LRfamily"),
                        color = "#02E7B9", 
                        type = 3, 
                        color.background = "white", 
                        size = 1.5,
                        caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                      style = "font-size:20px;")
                        )
                        
                                )
                            )
                        )
                    )
                ),
           
        ######### Table ############
        
        ## A page to display the data tables of the LRintercellNetworkDB
        tabPanel("Table", fluidPage(
          titlePanel("Table ligand-receptor atlas"),
          br(),  
            ## Age tabs
            tabsetPanel(
                
                ## E18.5-P0 tab
                tabPanel("E18.5-P0",
                
                    ## Render the E18.5-P0 table
                    dataTableOutput(outputId = 'tableP0_sig')),
                
                ## P1-P2 tab
                tabPanel("P1-P2", 
                
                    ## Render the P1-P2 table
                    dataTableOutput(outputId = 'tableP2_sig')),

                ## P4-P5 tab
                tabPanel("P4-P5", 
                
                    ## Render the P4-P5 table
                    dataTableOutput(outputId = 'tableP5_sig')),

                ## P8 tab
                tabPanel("P8", 
                
                    ## Render the P8 table
                    dataTableOutput(outputId = 'tableP8_sig')),

                ## P16 tab
                tabPanel("P16", 
                
                    ## Render the P16 table
                    dataTableOutput(outputId = 'tableP16_sig')),

                ## P30 tab
                tabPanel("P30", 
                
                    ## Render the P30 table
                    dataTableOutput(outputId = 'tableP30_sig')),
                
                ## Adult tab
                tabPanel("Adult",
                
                    ## Render the Adult table
                    dataTableOutput(outputId = 'tableAdult_sig')),

                ## Help tab
                tabPanel("Help", h1("Columns description"), htmltools::includeMarkdown("Data/Help_LRatlas_table.Rmd"))
                                           
                        )
                    )
                ),
        
        ########### Number of interactions ###############
          
        ## A page to display the number of interactions (heatmap)
        tabPanel("Number of interactions", fluidPage(
          titlePanel("Visualize the number of interactions"),
          br(),
            tabsetPanel(
                
                ## Heatmap tab
                tabPanel("Heatmap: Number of interactions",
                    sidebarLayout(
                        sidebarPanel(

                            ## sidebar parameters
                            width=3,

                            ## Ligand-Receptor Pairs selection
                            pickerInput(inputId = "htNinter_LR", label= "Ligand-Receptor:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple=T ),

                            ## ligand category selection
                                    pickerInput(inputId = "Cortex_htNinter_Ligand.category", label= "Ligand category:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T ),

                            ## Receptor category selection
                            pickerInput(inputId = "Cortex_htNinter_Receptor.category", label= "Receptor category:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T ),

                            ## Sending cell-type selection
                                    pickerInput(inputId = "htNinter_Sending.cluster", label= "Sending cell-type:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T ),

                            ## Target cell-type selection
                            pickerInput(inputId = "htNinter_Target.cluster", label= "Target cell-type:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T ),

                            accordion(
                                
                                ## LR filters tab (user inputs selection)
                                accordion_panel("LR filters",

                                    ## ligand family selection
                                    pickerInput(inputId = "Cortex_htNinter_Ligand.family", label= "Ligand family:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T ),

                                    ## Receptor family selection
                                    pickerInput(inputId = "Cortex_htNinter_Receptor.family", label= "Receptor family:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T ),

                                    ## Pathway selection
                                    pickerInput(inputId = "Cortex_htNinter_pathway", label= "pathway:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),

                                    ## Number of LRpair duplicates overall selection       
                                    pickerInput(inputId = "htNinter_numdup.LR_pair.interaction", label= "Number of LRpair duplicates in overall:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T ),

                                    ## Number of LRpair duplicates by fixing the Source
                                    pickerInput(inputId = "htNinter_numdup.LR_pair.cluster_L", label= "Number of LRpair duplicates by fixing the Source:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T ),

                                    ## Number of LRpair duplicates by fixing the Target
                                    pickerInput(inputId = "htNinter_numdup.LR_pair.cluster_R", label= "Number of LRpair duplicates by fixing the Target:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T ),

                                    ## S_inter score threshold 
                                    sliderInput(inputId = "htNinter_S_inter", label = "Filter LRpair with S_inter greater than or equal to the selected value:", min = 0, max = 1, value = 0, step = 0.05),

                                    ## S_intra score threshold
                                    sliderInput(inputId = "htNinter_S_intra", label = "Filter LRpair with S_intra greater than or equal to the selected value:", min = 0, max = 1, value = 0, step = 0.05),

                                    ## S_inter_diff score threshold
                                    sliderInput(inputId = "htNinter_S_inter_diff", label = "Filter LRpair with S_inter_diff greater than or equal to the selected value:", min = 0, max = 1, value = 0, step = 0.05),

                                            ),

                                ## Cell-type filters tab (user inputs selection)
                                accordion_panel("Cell-type filters",

                                    ## SOURCE class selection
                                    pickerInput(inputId = "htNinter_SOURCE.class", label= "SOURCE class:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple=T ),

                                    ## TARGET class selection
                                    pickerInput(inputId = "htNinter_TARGET.class", label= "TARGET class:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T )
                                        ),
                                
                                ## Plot settings tab
                                accordion_panel("Plot settings",
                                    
                                    ## Palette selection
                                    palettePicker(inputId = "ht_ninter_pal", label = "Select a color palette:", 
                                        choices = list(
                                            Sevilla=c("#FFFFFF", "#D6D6D7", "#C7B0B0","#C98F8F", "#CA6E6E", "#C64D4D", "#9F2D2D", "#761A1A", "#4C1111", "#270A0A", "#000000" ),
                                            Oslo=scico(11,palette = "oslo",direction=-1),
                                            PuOr = brewer_pal(palette = "PuOr", direction=-1)(11),
                                            GBBr = brewer_pal(palette = "BrBG", direction=-1)(11),
                                            BuYlRd = brewer_pal(palette = "RdYlBu", direction = -1)(11),
                                            BuPhRd = c("#00007F", "#0012FF", "#00A3FF","#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                                            GnPhRd = c("#2A363BFF", "#019875FF", "#99B898FF", "#FECEA8FF", "#FF847CFF", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                                            BuGyOrRd = c("#00008b", "#1e90ff", "#cccccc", "#ffa500","#ff4500"),
                                            BuGyRd= c("#00007F", "#0012FF", "#00A3FF", "#cccccc", "#E84A5FFF", "#C0392BFF", "#96281BFF"),
                                            viridis = viridis_pal(option = "viridis")(10),
                                            magma = viridis_pal(option = "magma")(10),
                                            inferno = viridis_pal(option = "inferno")(10),
                                            plasma = viridis_pal(option = "plasma")(10),
                                            cividis = viridis_pal(option = "cividis")(10)
                                                    ),

                                        textColor = c(rep("black", 2), rep("white", 14)),

                                        selected = "Sevilla"
                                                ),

                                    ## Plot height selection
                                    sliderInput(inputId = "htNinter_height", label = "Set the Heatmap height:", min = 400, max = 2000, value = 800, step = 50)
                                        )
                                    ),

                                    hr(),

                                    ## Apply filters button
                                    fluidPage(
                                      layout_columns(
                                        col.widths = c(2, 8, 2),
                                        div(),
                                        div(
                                          style = "align-content: center; align-self: center;",
                                          actionButton(inputId = "htNinter_button",label = "Apply Selection")
                                        ),
                                        div(),
                                      )
                                    )
                                  ),
                
                ## Plot Area
                mainPanel(

                    ## Age Tab set
                    tabsetPanel(
                        
                        ## E18.5-P0 tab
                        tabPanel("E18.5-P0",

                            ## Plot Area
                            fluidRow(

                                ## Heatmap Rendering (with Loading Animation) 
                                column(8, offset = 1, 
                                    shinycssloaders::withSpinner(
                                        plotlyOutput(outputId = 'P0_htNinter', height="850px",width="100%"),
                                        color = "#02E7B9", 
                                        type = 3, 
                                        color.background = "white", 
                                        size = 1.5,
                                        caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                                      style = "font-size:20px;")
                                    )
                                ),

                                ## Text Rendering and Plot Legend
                                column(2, textOutput("P0_htNinter_selectLR"),plotOutput(outputId = 'P0_lgd_ninter_ht', height="600px", width="100%"))
                            )
                        ),

                        ## P1-P2 tab
                        tabPanel("P1-P2",

                            ## Plot Area
                            fluidRow( 
                                
                                ## Heatmap Rendering (with Loading Animation)
                                column(8, offset = 1, 
                                    shinycssloaders::withSpinner(
                                        plotlyOutput(outputId = 'P2_htNinter',height="850px", width="100%"),
                                        color = "#02E7B9", 
                                        type = 3, 
                                        color.background = "white", 
                                        size = 1.5,
                                        caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                                      style = "font-size:20px;")
                                    )
                                ),
                                
                                ## Text Rendering and Plot Legend
                                column(2, textOutput("P2_htNinter_selectLR"),plotOutput(outputId = 'P2_lgd_ninter_ht', height="600px", width="100%"))
                            )
                        ),
                        
                        ## P4-P5 tab
                        tabPanel("P4-P5",

                            ## Plot Area
                            fluidRow( 
                                
                                ## Heatmap Rendering (with Loading Animation)
                                column(8, offset = 1,
                                    shinycssloaders::withSpinner(
                                        plotlyOutput(outputId = 'P5_htNinter',height="850px", width="100%"),
                                        color = "#02E7B9", 
                                        type = 3, 
                                        color.background = "white", 
                                        size = 1.5,
                                        caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)", 
                                                      style = "font-size:20px;")
                                    )
                                ),

                                ## Text Rendering and Plot Legend
                                column(2, textOutput("P5_htNinter_selectLR"), plotOutput(outputId = 'P5_lgd_ninter_ht', height="600px", width="100%"))
                            )
                        ),
                        
                        ## P8 tab
                        tabPanel("P8",

                            ## Plot Area
                            fluidRow(
                                
                                ## Heatmap Rendering (with Loading Animation)
                                column(8, offset = 1,
                                    shinycssloaders::withSpinner(
                                        plotlyOutput(outputId = 'P8_htNinter',height="850px", width="100%"),
                                        color = "#02E7B9", 
                                        type = 3, 
                                        color.background = "white", 
                                        size = 1.5,
                                        caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                                      style = "font-size:20px;")
                                    )    
                                ),

                                ## Text Rendering and Plot Legend
                                column(2, textOutput("P8_htNinter_selectLR"), plotOutput(outputId = 'P8_lgd_ninter_ht', height="600px", width="100%"))
                            )
                        ),
                        
                        ## P16 tab
                        tabPanel("P16",

                            ## Plot Area
                            fluidRow( 
                                
                                ## Heatmap Rendering (with Loading Animation)
                                column(8, offset = 1,
                                    shinycssloaders::withSpinner(
                                        plotlyOutput(outputId = 'P16_htNinter',height="850px", width="100%"),
                                        color = "#02E7B9", 
                                        type = 3, 
                                        color.background = "white", 
                                        size = 1.5,
                                        caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                                      style = "font-size:20px;")
                                    )
                                ),
                                
                                ## Text Rendering and Plot Legend
                                column(2, textOutput("P16_htNinter_selectLR"), plotOutput(outputId = 'P16_lgd_ninter_ht', height="600px", width="100%"))
                            )
                        ),
                        
                        ## P30 tab
                        tabPanel("P30",

                            ## Plot Area
                            fluidRow( 
                                
                                ## Heatmap Rendering (with Loading Animation)
                                column(8, offset = 1,
                                    shinycssloaders::withSpinner(
                                        plotlyOutput(outputId = 'P30_htNinter',height="850px", width="100%"),
                                        color = "#02E7B9", 
                                        type = 3, 
                                        color.background = "white", 
                                        size = 1.5,
                                        caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)", 
                                                      style = "font-size:20px;")
                                    )
                                
                                ),

                                ## Text Rendering and Plot Legend
                                column(2, textOutput("P30_htNinter_selectLR"), plotOutput(outputId = 'P30_lgd_ninter_ht', height="600px", width="100%"))
                            )
                        ),
                        
                        ## Adult tab
                        tabPanel("Adult",

                            ## Plot Area
                            fluidRow( 
                                
                                ## Heatmap Rendering (with Loading Animation)
                                column(8, offset = 1,
                                    shinycssloaders::withSpinner(
                                        plotlyOutput(outputId = 'Adult_htNinter',height="850px", width="100%"),
                                        color = "#02E7B9", 
                                        type = 3, 
                                        color.background = "white", 
                                        size = 1.5,
                                        caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                                      style = "font-size:20px;")
                                    ) 
                                
                                ),
                                
                                ## Text Rendering and Plot Legend
                                column(2, textOutput("Adult_htNinter_selectLR"), plotOutput(outputId = 'Adult_lgd_ninter_ht', height="600px", width="100%"))
                            )
                        ),
                        
                        ## Contiguous Age tab
                        tabPanel("Contiguous from E18.5-P0 to P30",

                            ## Plot Area
                            fluidRow( 
                                
                                ## Heatmap Rendering (with Loading Animation)
                                column(8, offset = 1,
                                    shinycssloaders::withSpinner(
                                        plotlyOutput(outputId = 'Cortex_htNinter',height="850px", width="100%"),
                                        color = "#02E7B9", 
                                        type = 3, 
                                        color.background = "white", 
                                        size = 1.5,
                                        caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)", 
                                                      style = "font-size:20px;")
                                    ) 
                                
                                ),
                                
                                ## Text Rendering and Plot Legend
                                column(2, textOutput("Cortex_htNinter_selectLR"), plotOutput(outputId = 'Cortex_lgd_ninter_ht', height="600px", width="100%"))
                            )
                         )
                    ),
                            
                            width=9
                )
                                               
            )
                                           
        ),
                                 
                ## Help tab           
                tabPanel("Help", htmltools::includeMarkdown("Data/Help_LRatlas_nbintr.Rmd"))
                )
            )
        ),

    ############# LR interactions ##################
        
        ## A page to display the Intercellular/Intracellular signaling plots
        tabPanel("Intercellular/Intracellular signaling",
            
            ## A page to display Ligand-Receptor pairs between cell-type pairs
            fluidPage(
              
              titlePanel("Visualize Ligand-Receptor pairs between cell-type pairs"),
            br(),
            
            sidebarLayout(
                
                ## Sidebar panel for user input
                sidebarPanel(


                    ## Sidebar parameters
                    width=3,

                    ##Ligand-Receptor pair selection                        
                    pickerInput(inputId = "Cortex_LRpairs_LR", label= "Ligand-Receptor:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),
                        
                    ## Source class selection
                    pickerInput(inputId = "Cortex_LRpairs_SOURCE.class", label= "SOURCE class:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),
                        
                    ## Target class selection
                    pickerInput(inputId = "Cortex_LRpairs_TARGET.class", label= "TARGET class:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),

                    ## Sending cell-type selection
                    pickerInput(inputId = "Cortex_LRpairs_cluster_L", label= "Sending cell-type:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE),multiple=T),

                    ## Target cell-type selection
                    pickerInput(inputId = "Cortex_LRpairs_cluster_R", label= "Target cell-type:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),
                    
                    accordion(
                    
                        ## Cell-type & LR filters
                        accordion_panel("LR filters",
                                                        
                            ## Pathway selection
                            pickerInput(inputId = "Cortex_LRpairs_pathway", label= "pathway:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),

                            ## Ligand category selection
                            pickerInput(inputId = "Cortex_LRpairs_Ligand.category", label= "Ligand category:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),

                            ## Receptor category selection
                            pickerInput(inputId = "Cortex_LRpairs_Receptor.category", label= "Receptor category:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),

                            ## Ligand family selection
                            pickerInput(inputId = "Cortex_LRpairs_Ligand.family", label= "Ligand family:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),

                            ## Receptor family selection
                            pickerInput(inputId = "Cortex_LRpairs_Receptor.family", label= "Receptor family:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),

                            ## Number of LRpair duplicates overall selection
                            pickerInput(inputId = "Cortex_LRpairs_numdup.LR_pair.interaction", label= "Number of LRpair duplicates in overall:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),

                            ## Number of LRpair duplicates by fixing the Source
                            pickerInput(inputId = "Cortex_LRpairs_numdup.LR_pair.cluster_L", label= "Number of LRpair duplicates by fixing the Source:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),
                        
                            ## Number of LRpair duplicates by fixing the Target
                            pickerInput(inputId = "Cortex_LRpairs_numdup.LR_pair.cluster_R", label= "Number of LRpair duplicates by fixing the Target:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=T),

                            ## S_inter score threshold
                            sliderInput(inputId = "Cortex_LRpairs_S_inter", label = "Filter LRpair with S_inter greater than or equal to the selected value:", min = 0, max = 1, value = 0, step = 0.05),
                        
                            ## S_intra score threshold
                            sliderInput(inputId = "Cortex_LRpairs_S_intra", label = "Filter LRpair with S_intra greater than or equal to the selected value:", min = 0, max = 1, value = 0, step = 0.05),
                        
                            ## S_inter_diff score threshold
                            sliderInput(inputId = "Cortex_LRpairs_S_inter_diff", label = "Filter LRpair with S_inter_diff greater than or equal to the selected value:", min = 0, max = 1, value = 0, step = 0.05),
                        
                        ## Apply filters button
                        fluidPage(
                          layout_columns(
                            col.widths = c(2, 8, 2),
                            div(),
                            div(
                              style = "align-content: center; align-self: center;",
                              actionButton( inputId = "Cortex_LRpairs_button", label = "Apply Selection")
                            ),
                            div(),
                          )
                        )
                      ),
                    
                    ## Plot settings tab
                    accordion_panel("Plot settings",
                        
                        ## X-axis variable selection
                        pickerInput(inputId = "Cortex_LRpairs_x.axis", label= "Select the variable on the x-axis:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F),

                        ## Y-axis variable selection   
                        pickerInput(inputId = "Cortex_LRpairs_y.axis", label= "Select the variable on the y-axis:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F),

                        ## X-axis Facet selection  
                        pickerInput(inputId = "Cortex_LRpairs_facetgrid.x", label= "Select the variable for the facet on the x-axis:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F),

                        ## Y-axis Facet selection  
                        pickerInput(inputId = "Cortex_LRpairs_facetgrid.y", label= "Select the variable for the facet on the y-axis:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F),

                        ## Point size selection
                        pickerInput(inputId = "Cortex_LRpairs_pt.size", label= "Select the variable for the point size:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F),
                        
                        ## Color gradient selection
                        pickerInput(inputId = "Cortex_LRpairs_fill", label= "Select the variable for the color gradient:", choices=NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE), multiple=F),

                        ## Color palette selection
                        palettePicker( inputId = "Cortex_LRpairs_pal", label = "Select a color palette:", 
                            choices = list(
                                lajolla = scico(11, palette="lajolla"),
                                bilbao = scico(11, palette="bilbao"),
                                acton= scico(11, palette="acton", direction=-1),
                                batlowK = scico(11, palette="batlowK", direction=-1),
                                davos= scico(11, palette="davos", direction=-1),
                                hawaii = scico(11, palette="hawaii", direction=-1),
                                imola= scico(11,palette="imola", direction=-1),
                                lapaz = scico(11,palette="lapaz", direction=-1),
                                nuuk = scico(11,palette="nuuk", direction=-1),
                                tokyo = scico(11,palette="tokyo", direction=-1),
                                turku = scico(11,palette="turku", direction=-1)
                                ),
                                    
                                textColor = c( rep("black", 11)), selected = "lajolla"
                                ),
                        
                        ## Plot height and width selection
                        sliderInput(inputId = "LRpairs_width", label = "Set the plot width:", min = 100, max = 2000, value = 500, step = 50),

                        sliderInput(inputId = "LRpairs_height", label = "Set the plot height:", min =100, max = 2000, value = 500, step = 50),

                        sliderInput(inputId = "LRpairs_step_height", label = "Set the gap between 2 y-values:", min =0, max = 1000, value = 30, step = 10),

                        sliderInput(inputId = "LRpairs_step_width", label = "Set the gap between 2 x-values:", min =0, max = 600, value = 30, step = 10),
                                                                    
                        )
                    ),
                    
                    hr(),

                    ## Apply settings button
                    fluidPage(
                      layout_columns(
                        col.widths = c(2, 8, 2),
                        div(),
                        div(
                          style = "align-content: center; align-self: center;",
                          actionButton( inputId = "Cortex_LRsetting_button", label = "Apply Selection")
                        ),
                        div(),
                      )
                    )  
                  ),
            
            ## Plot panel
            mainPanel(
                tabsetPanel(

                    # E18.5-P0 tab
                    tabPanel("E18.5-P0", (
                        
                        ## Plot area
                        div(style='overflow-y: scroll; position: relative', 
                        
                        ## Plot Rendering (with Loading Animation)
                        shinycssloaders::withSpinner(
                            plotlyOutput(outputId = 'P0_LRpairs', height="850px",width="100%"),
                            color = "#02E7B9", 
                            type = 3, 
                            color.background = "white", 
                            size = 1.5,
                            caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                          style = "font-size:20px;")
                            )
                        )
                    ), 
                        
                        ## Render the text output
                        textOutput("P0_selectLR")),

                    ## P1-P2 tab
                    tabPanel("P1-P2", (
                        
                        ## Plot area
                        div(style='overflow-y: scroll; position: relative', 
                        
                        ## Plot Rendering (with Loading Animation)
                        shinycssloaders::withSpinner(
                            plotlyOutput(outputId = 'P2_LRpairs', height="850px",width="100%"),
                            color = "#02E7B9", 
                            type = 3, 
                            color.background = "white", 
                            size = 1.5,
                            caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                          style = "font-size:20px;")
                            )
                        )
                    ), 
                        
                        ## Render the text output
                        textOutput("P2_selectLR")),

                    ## P4-P5 tab
                    tabPanel("P4-P5", (
                        
                        ## Plot area
                        div(style='overflow-y: scroll; position: relative', 
                        
                        ## Plot Rendering (with Loading Animation)
                        shinycssloaders::withSpinner(
                            plotlyOutput(outputId = 'P5_LRpairs', height="850px",width="100%"),
                            color = "#02E7B9", 
                            type = 3, 
                            color.background = "white", 
                            size = 1.5,
                            caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                          style = "font-size:20px;")
                            )
                        )
                    ), 
                        
                        ## Render the text output
                        textOutput("P5_selectLR")),

                    ## P8 tab
                    tabPanel("P8", (
                        
                        ## Plot area
                        div(style='overflow-y: scroll; position: relative', 
                        
                        ## Plot Rendering (with Loading Animation)
                        shinycssloaders::withSpinner(
                            plotlyOutput(outputId = 'P8_LRpairs', height="850px",width="100%"),
                            color = "#02E7B9", 
                            type = 3, 
                            color.background = "white", 
                            size = 1.5,
                            caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)", 
                                          style = "font-size:20px;")
                            )
                        )
                    ), 
                        
                        ## Render the text output
                        textOutput("P8_selectLR")),

                    ## P16 tab
                    tabPanel("P16", (
                        
                        ## Plot area
                        div(style='overflow-y: scroll; position: relative', 
                        
                        ## Plot Rendering (with Loading Animation)
                        shinycssloaders::withSpinner(
                            plotlyOutput(outputId = 'P16_LRpairs', height="850px",width="100%"),
                            color = "#02E7B9", 
                            type = 3, 
                            color.background = "white", 
                            size = 1.5,
                            caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                          style = "font-size:20px;")
                            )
                        )
                    ), 
                        
                        ## Render the text output
                        textOutput("P16_selectLR")),

                    ## P30 tab
                    tabPanel("P30", (
                        
                        ## Plot area
                        div(style='overflow-y: scroll; position: relative', 
                        
                        ## Plot Rendering (with Loading Animation)
                        shinycssloaders::withSpinner(
                            plotlyOutput(outputId = 'P30_LRpairs', height="850px",width="100%"),
                            color = "#02E7B9", 
                            type = 3, 
                            color.background = "white", 
                            size = 1.5,
                            caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                          style = "font-size:20px;")
                            )
                        )
                    ), 
                        
                        ## Render the text output
                        textOutput("P30_selectLR")),

                    ## Adult tab
                    tabPanel("Adult", (
                        
                        ## Plot area
                        div(style='overflow-y: scroll; position: relative', 
                        
                        ## Plot Rendering
                        shinycssloaders::withSpinner(
                            plotlyOutput(outputId = 'Adult_LRpairs', height="850px",width="100%"),
                            color = "#02E7B9", 
                            type = 3, 
                            color.background = "white", 
                            size = 1.5,
                            caption = div(strong("Rendering Your Plot"), br(), "(Please wait, processing time depends on your machine's specifications)",
                                          style = "font-size:20px;")
                            )
                        )
                    ), 
                        
                        ## Render the text output
                        textOutput("Adult_selectLR")),

                    ## Help tab
                    tabPanel("Help", htmltools::includeMarkdown("Data/Help_LRatlas_CCC.Rmd"))
                            
                            ),
                            
                            width=9
                        )
                    )
                )
            )
        )
    )  
)


