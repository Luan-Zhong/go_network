library(tidyverse)
library(shinyjs) #For load_data()

load_data <- function() {
  Sys.sleep(2)
  hide("loading_page")
  show("main_content")
}
