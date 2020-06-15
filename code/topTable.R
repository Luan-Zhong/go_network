# Last update: 2020.06.10
# Title: GO DE community in neo4j
# Environment: VM linux 64g
# Major change: 



######################
## Require packages ##
######################

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
goSummaries <- url("https://uofabioinformaticshub.github.io/summaries2GO/data/goSummaries.RDS") %>%
  readRDS() %>%
  mutate(
    gs_name = Term(id),
    gs_name = str_to_upper(gs_name), 
    gs_name = str_replace_all(gs_name, "[ -]", "_"),
    gs_name = paste0("GO_", gs_name)
  )





ah <- AnnotationHub()
ah %>%
  subset(rdataclass == "EnsDb") %>%
  subset(grepl("Mus musculus", species))
ensDb <- ah[["AH75036"]] #This is ensembl 98


# Load in our example data
topTable <- read_csv("data/topTable.csv") %>% 
  dplyr::rename("ID" = "gene_id") 


# Get the mapping between Entrez & ensembl
ens2entrez <- topTable %>% 
  subset(., select = c(ID, entrezid, gene_name))
ens2entrez$entrezid <- as.integer((ens2entrez$entrezid))
ens2entrez <- ens2entrez %>% 
  dplyr::rename("gene_id" = "ID")


#genes(ensDb) %>% # Get all genes encoded on ensDb
#   mcols() %>% # mcols: a DataFrame object containing the metadata columns. Columns cannot be named "seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end", "width", or "element".
#   .[c("gene_id", "entrezid", "gene_name")] %>% # Form a tibble with 3 coulumn, gene_id, entreid, gene_name
#   as.data.frame() %>%
#   as_tibble()


# all_go is all go gene sets
all_go <- msigdbr(species = "Mus musculus", category = "C5") %>% 
  left_join(goSummaries) %>% 
  dplyr::filter(
    shortest_path >3 | terminal_node, # Shortest path must be > 3 OR it's a terminal node
    entrez_gene %in% unlist(topTable$entrezid) # Restrict to genes detected in our dataset
    #entrez_gene %in% unlist(dplyr::filter(topTable,adj.P.Val < 0.05)$entrezid) # Or restrict to the DE genes only
  ) %>% 
  left_join(ens2entrez,
    by = c("entrez_gene" = "entrezid")
    )




## Pathway enrichment testing ##
######################################################################################################################


# Group Go set by entrez gene id
goByGene <- all_go %>%
  split(f = .$gene_id) %>%
  lapply(extract2, "gs_name")

# Group gene by gene set
goByID <- all_go %>%
  split(f = .$gs_name) %>%
  lapply(extract2, "gene_id") %>%
  lapply(unique)


# Estimate PWF (probability weight function)

# Get the gene length
gene_length <- ensDb %>% 
  transcriptLengths() %>% 
  as_tibble() %>% 
  group_by(gene_id) %>% 
  summarise(len = max(tx_len))

# Get the DE gene vector
de.vector <- topTable %>% 
  with(
    structure(PValue < 0.05, names = ID)
  ) 

# Get the gene length factor
length.vector <- gene_length %>% 
  dplyr::filter(gene_id %in% names(de.vector))%>% 
  with(
    structure(len, names = gene_id)
  )

pwf <- nullp(de.vector, bias.data = length.vector)

de.go.set <- goseq(pwf, gene2cat = goByGene) %>% 
  as_tibble() %>% 
  dplyr::select(category, over_represented_pvalue, numDEInCat, numInCat) %>% 
  mutate(fdr = p.adjust(over_represented_pvalue, "fdr")) %>% 
  dplyr::filter(fdr < 0.05)


de.go <- de.go.set %>%
  left_join(all_go, by = c("category" = "gs_name")) %>%
  subset(., select = -c(gene_id, human_gene_symbol, ontology)) %>%
  unique() %>% 
  arrange(over_represented_pvalue)
######################################################################################################################


## Jaccard index calculation ##
######################################################################################################################

# Filter nodes
goList <- de.go %>% # Separate genes into groups according to their gene name. The genes are named byreasons to filter by their functionality. This is the reason why separate them by gene name.
  split(f = .$category)



# Get a logical vector for each gene set
genes2de_go <- goList %>%
  lapply(function(x){
    topTable$gene_name %in% x$gene_symbol 
  }) %>%
  as_tibble() %>%
  mutate(gene_id = topTable$ID) %>%
  dplyr::select(gene_id, everything())


# Define the pairs for calculation of jaccard indices
go_Pairs <- combn(names(goList), 2) %>% 
  t()

# Parallel processing calculate jaccard index
mc <- detectCores() - 1 # Keep one core free
ji <- go_Pairs %>%
  nrow() %>% #Count how many pairs we have
  seq_len() %>% # Generate a sequence from 1:nPairs
  mclapply( # Apply this process to every value 
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

ji$weight %>% quantile(probs=seq(0, 1, by=0.05))
ji$weight %>% .[.>0] %>% quantile(probs=seq(0, 1, by=0.05))

# Histogram
ji$weight %>%
  hist(breaks = 50, main = paste("Distribution of DE gene sets Jaccard index"), xlim = c(0,0.9), ylim = c(1, 1500), xlab='Jaccard index', ylab='Frequency')

# Distribution of jaccard index result
ji$weight %>% 
  quantile(probs=seq(0, 1, by=0.05))


# Count # gene in a a gene set
num_genes <- genes2de_go %>%
  ncol() %>%
  seq_len() %>%
  lapply (function(x){
    genes2de_go[,x] %>%
      .[which(.==TRUE),] %>%
      nrow()}) %>% 
  tibble(names(genes2de_go), .) %>%
  dplyr::rename(., Gene_set ="names(genes2de_go)", numGenes = ".") %>% 
  .[-1,] %>%
  mutate(numGenes = vapply(numGenes, function(x){x[[1]]}, integer(1)))

######################################################################################################################


## Data cleaning & generate csv for neo4j ##
######################################################################################################################

# go_de.csv : all go genesets which difficial expressed
aaa <- all_go %>% 
  subset(select = c(gs_name, id, gs_cat, gs_subcat, species_name, shortest_path, longest_path)) %>% 
  unique()

de.go.set %>% 
  left_join(aaa, by = c("category" = "gs_name")) %>% 
  arrange(over_represented_pvalue) %>% 
  left_join(num_genes, by = c("category"="Gene_set")) %>% 
  write_csv(path = "/home/a1696632/neo4j-community-3.5.14/import/go_de.csv")


# ji.csv 
write_csv(x = ji, path = "/home/a1696632/neo4j-community-3.5.14/import/ji.csv")

######################################################################################################################



## neo4j ##
######################################################################################################################

con <- neo4r::neo4j_api$new(
  url = "http://10.150.8.113:7474", # http://10.150.8.113 is the VM IP address NB: is http:// rather than https://
  user = "neo4j",
  password = "" # No need password
)

con$ping()

#Clear any existing data

"
MATCH (n) 
DETACH DELETE n
" %>%
  call_neo4j(con)

con$get_version()
con$get_labels()
con$get_relationships()

# Import genes & properties (nodes)
# topTable.csv : All genes information 
# Making sure there is no "." in the *.csv rowname



# Import go de gene sets
# go_de.csv :  go de gene sets name
# Absolute path of go_de.csv is /home/a1696632/neo4j-community-3.5.14/import/go_de.csv

"
USING PERIODIC COMMIT 300
LOAD CSV WITH HEADERS FROM 'file:///go_de.csv' AS degenesets
MERGE (gode:GODE{
		gs_name: degenesets.category,
		over_represented_pvalue: degenesets.over_represented_pvalue,
		numDEInCat: degenesets.numDEInCat,
		numInCat: degenesets.numInCat,
		fdr: degenesets.fdr,
		go_id: degenesets.id,
		gs_cat: degenesets.gs_cat,
		gs_subcat: degenesets.gs_subcat,
		species_name: degenesets.species_name,
		numGenes: degenesets.numGenes
        })
WITH gode
MATCH (gode)
RETURN (gode)
" %>%
  call_neo4j(con)



# Import jaccard
# Import the content in ji.csv and store in variable jaccard

"
  USING PERIODIC COMMIT 1000
  LOAD CSV WITH HEADERS FROM 'file:///ji.csv' AS jaccard_index
  WITH jaccard_index
  MATCH (set1:GODE{gs_name:jaccard_index.Set1}), (set2:GODE{gs_name:jaccard_index.Set2})
  MERGE (set1) -[ji:jaccard {jaccard_weight: jaccard_index.weight}]-> (set2)
  SET ji.jaccard_weight=toFloat(ji.jaccard_weight)
" %>%
  call_neo4j(con)

con$get_relationships()


# Community detection

# Louvain

# Return result of Louvain, # of community, and # nodes in each community
louvain <- "
CALL algo.louvain.stream('GODE','jaccard', {iterations:20, write:true, writeProperty:'louvain_01', weightProperty: 'jaccard_weight', direction:'BOTH'})
YIELD nodeId,community
RETURN community, count(*) AS communitySize ORDER BY community DESC
" %>%
  call_neo4j(con)


graph = startGraph("http://10.150.8.113:7474/db/data")

query <- "CALL algo.louvain('GODE','jaccard', {iterations:20, write:true, writeProperty:'louvain_01', weightProperty: 'jaccard_weight', direction:'BOTH'})" %>% 
  cypher(graph, .)


"
CALL algo.triangleCount.stream('GODE', 'jaccard', {concurrency:4})
YIELD nodeId, triangles, coefficient

RETURN algo.asNode(nodeId).gs_name AS name, triangles, coefficient
ORDER BY coefficient DESC
" %>% 
  call_neo4j(con)


#####################
## Auto anootation ##



G <- "MATCH p=()-[r:jaccard]->() RETURN p" %>%
  call_neo4j(con, type = "graph")

# Create dataframe for nodes
G$nodes <- G$nodes %>%
  unnest_nodes(what = "properties")


# Create dataframe for relationships
G$relationships <- G$relationships %>% 
  unnest_relationships() %>%
  dplyr::select(startNode, endNode, type, everything()) %>%
  dplyr::rename(., neo4j_relid = "id")

# Create the networkd3 node id
# networkD3 requires d3 node id rather then neo4j node id 
# D3 node id starts from 0, neo4j node id starts from 1
d3_nodeid <- 1:nrow(G$nodes)-1
G$nodes <- G$nodes %>%  data.frame(., d3_nodeid) %>%
  dplyr::rename(., neo4j_nodeid = "id") %>%
  as_tibble()

# Create neo4j node id and network d3 reference table 
neo4j2d3 <- data.frame(G$nodes$neo4j_nodeid,G$nodes$d3_nodeid, stringsAsFactors=FALSE) %>% 
  dplyr::rename(., neo4j_nid="G.nodes.neo4j_nodeid", d3_nid="G.nodes.d3_nodeid") %>% 
  as_tibble()


# Match the neo4j node id to networkD3 nodeid
match <- G$relationships %>%
  as_tibble() %>%
  subset(., select=-c(type, neo4j_relid)) %>%
  dplyr::rename(., neo4j_nid="startNode") %>%
  left_join(., neo4j2d3, by = "neo4j_nid") %>%
  subset(., select=-c(neo4j_nid)) %>%
  dplyr::rename(., source="d3_nid", neo4j_nid="endNode") %>%
  left_join(., neo4j2d3, by = "neo4j_nid") %>%
  subset(., select=-c(neo4j_nid)) %>%
  dplyr::rename(., target="d3_nid")



# Count words
word_freq <- function(x, omit = c("GO", "OF", "VS", "TO", "AND", "IN", "BY", "AT", "AS", "UP", "DN", "NEGATIVE", "POSITIVE"), n = 4) {
  x$gs_name %>%  # Get all the geneset name
    .[!duplicated(.)] %>% # Remove duplicates if there is
    gsub("\\_"," ", .) %>% # Replace underscore by space 
    strsplit(.," ") %>% # Split each term as multiple words
    unlist() %>% 
    table () %>% # Put them in table
    enframe(name = "Word", value = "Count") %>% 
    arrange(desc(Count)) %>%    # Order by Count from large to small
    dplyr::filter(!Word %in% omit) %>% # Remove the omit words
    dplyr::slice(seq_len(n))
  # dplyr::filter(value > n) # Keep the Count is > 3 term, ignore all the term is number
}


# Community list

community_list <- G$nodes %>%
  split(f = .$louvain_01) %>% 
  lapply(word_freq, n = 4)



word2omit <- c(
  "A", "ACQUIRES", "ACTIVITY", "ALL", "AN", "AND", "ANY", "ARE", 
  "AS", "AT", "BE", "BY", "CELL", "CELLS", "COMBINING", "CONSIDERED", 
  "DN", "ETC", "EXTENT", "FREQUENCY", "FROM", "GO", "GROUP", "IN", 
  "IS", "NEGATIVE", "OF", "OR", "PATHWAYS", "POSITIVE", "PROCESS", 
  "RATE", "RESULT", "RESULTING", "RESULTS", "SOME", "TERMS", "THAT", 
  "THE", "TO", "UP", "VS", "WHICH", "WITH", "WITHIN", "BETWEEN", "INVOLVED",
  "ITS", "BODY"
)
# Count words
definition_freq <- function(x, omit = word2omit, n = 5) {
  Definition(x$go_id) %>%  # Get all the geneset name
    str_to_upper() %>% 
    .[!duplicated(.)] %>% # Remove duplicates if there is
    gsub("\\_"," ", .) %>% # Replace underscore by space 
    strsplit(.," ") %>% # Split each term as multiple words
    unlist() %>% 
    str_remove_all("[,\\.\\(\\)]") %>%
    table () %>% # Put them in table
    enframe(name = "Word", value = "Count") %>% 
    arrange(desc(Count)) %>%    # Order by Count from large to small
    dplyr::filter(
      !Word %in% omit,
      Count > 1
    ) %>% # Remove the omit words
    dplyr::slice(seq_len(n))
  # dplyr::filter(value > n) # Keep the Count is > 3 term, ignore all the term is number
}


# Community list

definition_list <- G$nodes %>%
  split(f = .$louvain_01) %>% 
  lapply(definition_freq, n = 5)





## Community Jaccard index ##
######################################################################################################################

# Return belongingship of each gene set and community
set2com <- "
CALL algo.louvain.stream('GODE','jaccard', {iterations:20, write:true, weightProperty: 'jaccard_weight', direction:'BOTH'})
YIELD nodeId,community
RETURN algo.asNode(nodeId).gs_name as name,community
" %>% 
  call_neo4j(con) 

set2com  <- data.frame(set2com$name, set2com$community) %>%
  as_tibble()

# Match community with the gene name
gene_in_com <- genes2de_go %>% 
  pivot_longer(cols = contains("GO"), names_to = "geneset", values_to = "is_in") %>% 
  dplyr::filter(is_in) %>% 
  dplyr::select(gene_id, geneset) %>%
  left_join(set2com, by = c( "geneset" = "value")) %>% 
  dplyr::rename(community = "value.1") %>% 
  arrange(gene_id)

# Seperate genes by community
comList <- gene_in_com %>% 
  split(f = .$community)

# Get a logical vector for each community
genes2com <- comList %>%
  lapply(function(x){
    topTable$ID %in% x$gene_id
  }) %>%
  as_tibble() %>%
  mutate(gene_id = topTable$ID) %>%
  dplyr::select(gene_id, everything())

# Define community pairs for calculation of jaccard indices
com_Pairs <- combn(names(comList), 2) %>% 
  t()


# Calculate jaccard index between community
com_ji <- com_Pairs %>%
  nrow() %>% #Count how many pairs we have
  seq_len() %>% # Generate a sequence from 1:nPairs
  mclapply( # Apply this process to every value in thea %>% gr sequence
    FUN = function(x){
      s1 <- com_Pairs[x,1] # Set variable s1 as com_Pairs column 1
      s2 <- com_Pairs[x,2] # Set variable s1 as com_Pairs column 2
      n <- length( intersect(comList[[s1]][[1]], comList[[s2]][[1]]) )
      w <- n / length( union(comList[[s1]][[1]], comList[[s2]][[1]]) )
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
  arrange(desc(weight))

# Top 5 gene in each community
top_5_gene <- comList %>% 
  lapply(function(x){sort(table(x$gene_id), decreasing = TRUE)[1:5]})


# Generate data for neo4j
write_csv(x = com_ji, path = "/home/a1696632/neo4j-community-3.5.14/import/com_ji.csv")

com_node <- 
  comList %>% 
  vapply(function(x){c(genes = length(unique(x$gene_id)), genesets = length(unique(x$geneset)))}, integer(2)) %>% 
  cbind (.) %>%
  t() %>% 
  data.frame()

community <- rownames(com_node)
rownames(com_node) <- NULL
com_node <-cbind(community, com_node) %>% 
  write_csv(path = "/home/a1696632/neo4j-community-3.5.14/import/com.csv")





"
USING PERIODIC COMMIT 300
LOAD CSV WITH HEADERS FROM 'file:///com.csv' AS group
MERGE (com:COM{
		comm: group.community,
		genes: group.genes,
		genesets: group.genesets
        })
WITH com
MATCH (com)
RETURN (com)
" %>% 
  call_neo4j(con)


"
USING PERIODIC COMMIT 1000
LOAD CSV WITH HEADERS FROM 'file:///com_ji.csv' AS com_jaccard_index
WITH com_jaccard_index
MATCH (com_set1:COM{comm:com_jaccard_index.Set1}), (com_set2:COM{comm:com_jaccard_index.Set2})
MERGE (com_set1) -[com_ji:com_jaccard {com_jaccard_weight: com_jaccard_index.weight}]-> (com_set2)
SET com_ji.jaccard_weight=toFloat(com_ji.jaccard_weight)
" %>% 
  call_neo4j(con) 



# Load neo4j network to R

G <- "MATCH p=()-[r:com_jaccard]->() RETURN p" %>%
  call_neo4j(con, type = "graph")

# Create dataframe for nodes

G$nodes <- G$nodes %>%
  unnest_nodes(what = "properties")

# Create dataframe for relationships

G$relationships <- G$relationships %>% 
  unnest_relationships() %>%
  dplyr::rename(., neo4j_relid = "id")


# Create the networkd3 node id
# networkD3 requires d3 node id rather then neo4j node id 
# D3 node id starts from 0, neo4j node id starts from 1
d3_nodeid <- 1:nrow(G$nodes)-1
G$nodes <- G$nodes %>%  data.frame(., d3_nodeid) %>%
  dplyr::rename(., neo4j_nodeid = "id") %>%
  as_tibble()

neo4j2d3 <- data.frame(G$nodes$neo4j_nodeid,G$nodes$d3_nodeid, stringsAsFactors=FALSE) %>% 
  dplyr::rename(., neo4j_nid="G.nodes.neo4j_nodeid", d3_nid="G.nodes.d3_nodeid") %>% 
  as_tibble()


match <- G$relationships %>%
  left_join(neo4j2d3, by = c("startNode"= "neo4j_nid")) %>%
  dplyr::rename(source="d3_nid") %>%
  left_join(neo4j2d3, by=c("endNode" = "neo4j_nid")) %>%
  dplyr::rename(target="d3_nid") %>% 
  subset(select=-c(startNode, endNode, type, neo4j_relid)) %>%
  dplyr::rename(value="com_jaccard_weight") %>%
  as_tibble()

match$value <- as.numeric(match$value)
# Return singletons geneset

# This retruns all the genesets are isolatd (= no relationship w/ other nodes)
# S <-"MATCH(n:GOCC) WHERE NOT (n)-[:jaccard]-(:GOCC) RETURN n"  %>%
#   call_neo4j(con)



## Visualization ##
######################################################################################################################

# Create tooltip for network
# Tooltip is for node content preveiw (aka pop-up window)

popup <- "
d3.selectAll('.xtooltip').remove();
d3.select('body').append('div')
.attr('class', 'xtooltip')
.style('position', 'fixed')
.style('border-radius', '0px')
.style('padding', '5px')
.style('opacity', '0.85')
.style('background-color', '#161823')
.style('box-shadow', '2px 2px 6px #161823')

.html(
  '======== Community ' + d.name +' information ========'
  + '<br>' +
  'Community name : ' + d.name
  + '<br>' +
  'Number of genes are included :  ' + d.nodesize
)
.style('right', '50px')
.style('bottom', '50px')
.style('color', d3.select(this).style('fill'))
;"


# Rescale the centrality for node size visulization
# Save distribution of centrality
# no value is zero in centrality in this case
# cause we deleted all the edges jaccard index <0 ???
#centrality_dis <- G$nodes$centrality %>%
#  quantile(probs=seq(0,1,by=0.2)) %>%
#  as_tibble()
#
# Recode
# 0:0.3=1 means from 0 to 0.3 rescale as 1
#G$nodes$size <- car::recode(G$nodes$centrality, "0:3 = 0.7; 3:7 = 1.1; 7:10 = 2; 10:15 = 3; 15:30 = 5; else = 7")
# Construct newtork by network D3



gode_nw <- forceNetwork(Links = match,
                        Nodes = G$nodes,
                        NodeID = "comm", # D3 graph node name
                        Source = "source", # Start node
                        Target = "target", # End node
                        Group = "comm",
                        Nodesize = "genes",
                        radiusCalculation = JS("d.nodesize*0.01"),
                        linkDistance = JS('function(d) {', 'return d.value*1000;', '}'), # Show linkdistance as jaccard index*1000
                        width = 500,
                        height = 500,
                        opacity = 0.9,
                        fontSize = 16,
                        zoom = TRUE,
                        legend = TRUE,
                        bounded = TRUE, # Create a bounded graph
                        clickAction = popup,
                        colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);")
)
gode_nw$x$links$value <- G$relationships$com_jaccard_weight

# Network search box by htmltool
# Searching only for GO id
gode_nw <- htmlwidgets::onRender(
  gode_nw,
  '
function(el,x){
debugger;
  var optArray = [];
  for (var i = 0; i < x.nodes.name.length - 1; i++) {
    optArray.push(x.nodes.name[i]);
  }

  optArray = optArray.sort();

  $(function () {
    $("#search").autocomplete({
      source: optArray
    });
  });

  d3.select(".ui-widget button").node().onclick=searchNode;

  function searchNode() {
    debugger;
    //find the node

    var selectedVal = document.getElementById("search").value;
    var svg = d3.select(el).select("svg");
    var node = d3.select(el).selectAll(".node");

    if (selectedVal == "none") {
      node.style("stroke", "white").style("stroke-width", "1");
    } else {
      var selected = node.filter(function (d, i) {
        return d.name != selectedVal;
      });
      selected.style("opacity", "0");
      var link = svg.selectAll(".link")
      link.style("opacity", "0");
      d3.selectAll(".node, .link").transition()
        .duration(5000)
        .style("opacity", 1);
    }
  }
}
  '
)

browsable(
  attachDependencies(
    tagList(
      tags$head(
        tags$link(
          href="http://code.jquery.com/ui/1.11.0/themes/smoothness/jquery-ui.css",
          rel="stylesheet"
        )
      ),
      HTML(
        '
  <div class="ui-widget">
      <input id="search">
      <button type="button">Search</button>
  </div>
  '
      ),
      gode_nw
    ),
    list(
      rmarkdown::html_dependency_jquery(),
      rmarkdown::html_dependency_jqueryui()
    )
  )
)



gode_nw %>%
  htmlwidgets::prependContent(htmltools::tags$h1("go_cc_networkD3")) %>% # Add title to plot
  saveNetwork(file = 'go_cc_networkD3.html') # Save the plot as html

## Shiny ##
######################################################################################################################




######################################################################################################################
sessionInfo()

