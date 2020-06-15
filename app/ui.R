########
## ui ##
########

library(DT) # For data table
library(plotly)
library(shinycssloaders) # For graph loader
library(msigdbr)
library(shiny)
library(shinydashboard)
library(shinyjs) # Loader at the beginning for annotationhub loading

UIset <- list(headTitle = "Gene sets community network",
              skin_color = "purple"
)


header <- dashboardHeader(
  title = UIset$headTitle
  
) # End of dashboardHeader()


############################################################################################################################################

sidebar <- dashboardSidebar(
  sidebarMenu(
    
    # ------------- File upload box (Start) ------------------
    fileInput(
      inputId = 'upload', 
      label = 'Upload a gene expression file',
      accept = c(
        'text/csv',
        'text/comma-separated-values',
        '.csv'
        )
      ),
    # -------------File upload box (End) ------------------
    
    
    
    
    
    # ------------- Biological filter (Start) -------------
    menuItem(
      "Biological filter", tabName="Biological", icon = icon("kiwi-bird"),
      
      # Select species
      selectInput(
        inputId = "species", 
        label = "Select species", 
        msigdbr_show_species(),
        selected = "Homo sapiens"
        ),
      
      
      # Select gene ID type
      # selectInput(
      #   "idtype", 
      #   "Select gene ID type", 
      #   c("--Non-specified--", "Entrez ID", "Ensembl gene ID", "Gene symbol"), 
      #   selected = "--Non-specified--"
      # ),
      
      
      # Select database
      selectInput(
        "dbtype", 
        "Select database", 
        c("HALLMARK", "KEGG", "Reactome", "GO"), 
        selected = "GO"
      )
      
    ),
    # ------------- Biological filter (End) -------------
    
    
    
    
    
    
    # ------------- Network filter (Start) -------------
    menuItem(
      "Network analysis filter", tabName="Network", icon = icon("project-diagram"),
      
      
      
      selectInput(
        "Connection methods", 
        "Gene sets connect by", 
        c("Jaccard index"),
        selected = "Jaccard index"
      )
      
      
      
    ),
    # ------------- Network filter (End) -------------
    
    
    
    
    # ------------- Submit button (Start) ------------
    submitButton("Start analysis")
    # ------------- Submit button (End) -------------
    
    
    
  ) # End of sidebarMenu()
  
) # End of dashboardSidebar()

############################################################################################################################################

body <- dashboardBody(
  
  fluidRow(

    tabBox(width = 12,

           # -------------  <Tab> Welcome page (Start) -------------
           tabPanel(
             status = "primary",
             title = "Hello World!",
             div(
               h5("Welcome to gene sets network analysis.
                  This is a shiny webapp for gene sets enrichment network analysis.
                  ")
             )
           ),
             
           # -------------  <Tab> Welcome page (End) -------------
           
           
           # -------------  <Tab> Show the input file (Start) -------------
           tabPanel(
             status = "primary",
             title = "Input gene list",
             
             # Loading screen
             useShinyjs(),
             div(
               id = "loading_page",
               h3("Loading AnnotationHub..."),
               h5("Please upload your file after this is finished")
             ),
             hidden(
               div(
                 id = "main_content",
                 h5("AnnotationHub is ready, Please upload your expression file")
               )
             ),
             
             # Show table
             div(
               h3("1. Raw data"),
               h5("This is table shows all the input data.")
             ),
             DT::dataTableOutput("input_table") %>% 
               withSpinner()
             ),
           # -------------  <Tab> Show the input file (End) -------------
           
           
           
           
           
           
           
           # -------------  <Tab> Community network (Start) -------------
           tabPanel(
             status = "success",
             title = "Pathway enrichment testing result",
             # Loading screen
             useShinyjs(),
             div(
               id = "PWF_loading",
               h3("Calculation in progress..."),
               h5("This takes around 50s")
             ),
             hidden(
               div(
                 id = "main_content",
                 h5("Calculation completed")
               )
             ),
             DT::dataTableOutput("de"),
             plotOutput("PWF"),
             downloadButton(outputId = "pwf_download", label = "Download PWF")
             
           ),
           # -------------  <Tab> Community network (End) -------------
           
           
           
           
           # -------------  <Tab> Jaccard index (Start) -------------
           tabPanel(
             status = "success",
             title = "Gene sets jaccard index result",
             # Loading screen
             useShinyjs(),
             div(
               id = "genesets_ji_loading",
               h3("Calculation in progress...")
             ),
             hidden(
               div(
                 id = "main_content",
                 h5("Calculation completed")
               )
             ),
             verbatimTextOutput("ji_dist"), 
             plotOutput("ji_hist"),
             downloadButton(outputId = "ji_hist_download", label = "Download histogram"),
             DT::dataTableOutput("genesets_ji") %>% 
               withSpinner()
           )
           # -------------  <Tab> Community network (End) -------------
      
    )
    
    
    
    
    
    
  ) # End of fluidRow()
) # End of dashboardBody()

ui <- dashboardPage(
  skin = UIset$skin_color, # Theme color as purple
  header,
  sidebar,
  body
)

