library(tidyverse)
library(AnnotationHub) # Mapping Entrez and ensembl
library(msigdbr) # Featch msigDB information
library(GO.db)
library(networkD3) # Network d3 graph
library(car) # Rescale node size
library(parallel) # Parallel running in R
library(neo4r) # Connect neo4j and R
library(RNeo4j)
library(magrittr) # For %>%
library(htmltools) # Network search box
library(goseq) # Enrichment testing
library(shiny)
library(shinydashboard)
library(shinyjs) # Loading screen

goSummaries <- reactive({
  url("https://uofabioinformaticshub.github.io/summaries2GO/data/goSummaries.RDS") %>%
    readRDS() %>%
    mutate(
      gs_name = Term(id),
      gs_name = str_to_upper(gs_name), 
      gs_name = str_replace_all(gs_name, "[ -]", "_"),
      gs_name = paste0("GO_", gs_name)
    )
}
)


############
## server ##
############


server <- shinyServer(function(input, output, session) {
  
  ensDb <- reactive({
    ah <- AnnotationHub(localHub = TRUE)
    ah %>%
      subset(rdataclass == "EnsDb") %>%
      subset(grepl(input$species, species))
    ensDb <- ah[["AH75011"]] #This is ensembl 98
  })

  
  
  
  # Load annotation hub
  ens2entrez <- reactive({
    # Get the mapping between Entrez & ensembl
    ens2entrez <-
      genes(ensDb()) %>% # Get all genes encoded on ensDb
      mcols() %>% # mcols: a DataFrame object containing the metadata columns. Columns cannot be named "seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end", "width", or "element".
      .[c("gene_id", "entrezid", "gene_name")] %>% # Form a tibble with 3 coulumn, gene_id, entreid, gene_name
      as.data.frame() %>%
      as_tibble()
  })


  
  
  # Load data
  topTable <- reactive({
    inFile <- input$upload
    if (is.null(inFile))
      return(NULL)
    load_data()
    topTable <- read_csv(inFile$datapath, col_names = TRUE) %>% 
      inner_join(ens2entrez(), by = c("ID" = "gene_id")) %>%
      dplyr::select(ID, gene_name, everything()) 
    return(topTable)
  })
  
  
  # ------------- File upload (Start) ------------------
  output$input_table <- DT::renderDataTable({
    topTable() %>% 
      DT::datatable() %>% 
      DT::formatRound(
        columns = c("AveExpr", "logFC", "t"),
        digits = 2
      ) %>%
      DT::formatSignif(
        columns = c("P.Value", "adj.P.Val"),
        digits = 3
      )
  })
  # ------------- File upload (End) ------------------
  
  
  # all_go is all go gene sets
  all_go <- reactive({
    msigdbr(species = input$species, category = "C5") %>% 
      left_join(goSummaries()) %>% 
      dplyr::filter(
        shortest_path >3 | terminal_node, # Shortest path must be > 3 OR it's a termina node
        entrez_gene %in% unlist(topTable()$entrezid) # Restrict to genes detected in our dataset
        #entrez_gene %in% unlist(dplyr::filter(topTable,adj.P.Val < 0.05)$entrezid) # Or restrict to the DE genes only
      ) %>% 
      left_join(
        unnest(ens2entrez(), entrezid), 
        by = c("entrez_gene" = "entrezid")
      )
  })

  

  
  ## Pathway enrichment testing ##
  ######################################################################################################################

  
  
  #Group Go set by entrez gene id
    
   
   
    
    # Group Go set by entrez gene id
    goByGene <- reactive({
      all_go() %>%
        split(f = .$gene_id) %>%
        lapply(extract2, "gs_name")
    })

    
    # Group gene by gene set
    goByID <- reactive({
      all_go() %>%
        split(f = .$gs_name) %>%
        lapply(extract2, "gene_id") %>%
        lapply(unique)
    })

    
    
    # Estimate PWF (probability weight function)
    
    # Get the gene length
    gene_length <- reactive({
      ensDb() %>% 
        transcriptLengths() %>% 
        as_tibble() %>% 
        group_by(gene_id) %>% 
        summarise(len = max(tx_len))
    })

    
    # Get the DE gene vector
    de.vector <- reactive({
      topTable() %>% 
        with(
          structure(adj.P.Val < 0.05, names = ID)
        )
    })

    
    # Get the gene length factor
    length.vector <- reactive({
      gene_length() %>% 
        dplyr::filter(gene_id %in% names(de.vector()))%>% 
        with(
          structure(len, names = gene_id)
        )
    })

    
    pwf <- reactive({
      pwf <- nullp(de.vector(), bias.data = length.vector())
      pwf
    })
    
    de.go.set <- reactive({
      de.go.set <- goseq(pwf(), gene2cat = goByGene()) %>%
        as_tibble() %>%
        dplyr::select(category, over_represented_pvalue, numDEInCat, numInCat) %>%
        mutate(fdr = p.adjust(over_represented_pvalue, "fdr")) %>%
        dplyr::filter(fdr < 0.05)
      de.go.set
    })
    
    
    
    
    
    #--------Table DE gene sets list (Start)------------
    output$de <- DT::renderDataTable({
      datatable(data = de.go.set(), 
                # Table title
                caption = htmltools::tags$caption (style = 'text-align: center;', 'Table 2: ', htmltools::em('Differential expression gene sets list')),
                # Filter for each column
                filter = 'top',
                # Download button
                extensions = 'Buttons', 
                options = list(autoWidth = TRUE,
                               searchHighlight = TRUE,
                               dom = 'lBfrtip',
                               scrollX = TRUE,
                               #fixedColumns = list(leftColumns =2, rightColumns = 1),
                               buttons = c('copy', 'csv', 'pdf'),
                               lengthMenu = c(10, 20, 50, -1)
                )
      )
    })
    #--------Table DE gene sets list (End)------------
    
    
    
  
    de.go <- reactive({
      de.go <- de.go.set() %>%
        left_join(all_go(), by = c("category" = "gs_name")) %>%
        subset(., select = -c(gene_id, human_gene_symbol, ontology)) %>%
        unique() %>% 
        arrange(over_represented_pvalue)
    })
    
    
    
    
    
    
    # ------------- Plot PWF (Start) -----------------
    output$PWF <- renderPlot({
      # Loader
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message = 'Calculation in progress',
                   detail = 'This may take a while...')
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.5)
      }
      plotPWF(pwf())
      })
    
    # ------------- Plot PWF (End) ------------------
    
    
    
    
    
    
    # ------------- download PWF button (Start) -----------------
    
    output$pwf_download <- downloadHandler(
      filename = function(){
        # Specify file name
        paste("pwf", "png", sep = ".")
        },
      content = function(file){
        png(file)
        # Create plot
        plotPWF(pwf())
        # Close the device
        dev.off()
      }
    )
    
    # ------------- download PWF (End) -----------------
    
    
    
    
    
    ## Jaccard index calculation ##
    ######################################################################################################################
    
    # Filter nodes
    goList <- reactive({
      goList <- de.go() %>% # Separate genes into groups according to their gene name. The genes are named byreasons to filter by their functionality. This is the reason why separate them by gene name.
        split(f = .$category)
    })
    
    
    # Get a logical vector for each gene set
    genes2de_go <- reactive({
      genes2de_go <- goList() %>%
        lapply(function(x){
          topTable()$gene_name %in% x$gene_symbol
          }) %>%
        as_tibble() %>%
        mutate(gene_id = topTable$ID) %>%
        dplyr::select(gene_id, everything())
    })
    
    
    ji <- reactive({
      # Define the pairs for calculation of jaccard indices
      go_Pairs <- combn(names(goList()), 2) %>% 
        t()
      
      goByID <- goByID() # Save the goByID() content in variable goByID, otherwise goByID[[s1]] will not working
      
      # Parallel processing calculate jaccard index
      mc <- detectCores() - 1 # Keep one core free
      ji <- go_Pairs %>%
        nrow() %>% #Count how many pairs we have
        seq_len() %>% # Generate a sequence from 1:nPairs
        mclapply( # Apply this process to every value in thea %>% gr sequence
          FUN = function(x){
            s1 <- go_Pairs[x, 1] # Set variable s1 as go_Pairs column 1
            s2 <- go_Pairs[x, 2] # Set variable s1 as go_Pairs column 2
            n <- length( intersect(goByID[[s1]], goByID[[s2]]) )
            w <- n / length( union(goByID[[s1]], goByID[[s2]]) )
            tibble(
              Set1 = s1, # Copy content in s1 to the new table column 1
              Set2 = s2, # Copy content in s1 to the new table column 2
              n = n,
              weight = w
            )
          }, 
          mc.cores = mc) %>% # Split the calculation as number of mc cores (which mc=15 in here)
        bind_rows() %>%
        dplyr::filter(!is.na(weight), !is.nan(weight)) %>% 
        .[which(.$weight > 0),] %>% # Select value is greater then 0
        arrange(desc(weight))
    })
    
    # Adding jaccard index explaination
    
    
    
    
    
    
    #--------Table gene sets Jaccard index (Start)------------
    output$genesets_ji <- DT::renderDataTable({
      ji() %>% 
        dplyr::rename("common genes" = n, "Jaccard index" = weight) %>% 
      DT::datatable() %>%
      DT::formatRound(
        columns = "Jaccard index",
        digits = 3
      ) 
        
    })
    #--------Table gene sets Jaccard index (End)------------
    
    
    
    
    
    #--------Gene sets Jaccard index distribution (Start)------------
    output$ji_dist <- renderPrint({
      ji()$weight %>%
        quantile(probs=seq(0, 1, by=0.05))
    })
    #--------Gene sets Jaccard index distribution (End)------------
    
    
    
    
    
    
    #--------Gene sets Jaccard index histogram (Start)------------
    output$ji_hist <- renderPlot({
      ji()$weight %>%
        hist(breaks = 50, main = paste("Histogram of DE gene sets ji>0 & shortest path>3"), xlim = c(0,1))
    })
    
    #--------Gene sets Jaccard index histogram (End)------------
    
    
    
    
    
    
    
    #------------- download Jaccard index histogram button (Start) -----------------
    
    output$ji_hist_download <- downloadHandler(
      filename = function(){
        # Specify file name
        paste("Jaccard index histogram", "png", sep = ".")
      },
      content = function(file){
        png(file)
        # Create plot
        ji()$weight %>%
          hist(breaks = 50, main = paste("Histogram of DE gene sets ji>0 & shortest path>3"), xlim = c(0,1))
        # Close the device
        dev.off()
      }
    )
    
    # ------------- download Jaccard index histogram button (End) -----------------
    
    
    
  } # End of function(input, output) {}
  ) # End of shinyServer()




