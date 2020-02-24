

library(shiny)
library(ggplot2)
library(Matrix)
library(reshape2)
library(dplyr)
library(plotly)
library(shinyjs)
library(shinyWidgets)
library(shinyBS)
library(jpeg)
library(scales)
library(viridisLite)
library(ggrepel)
source("R/HeatImage.R")
source("R/VisCelloImport.R")
source("R/strip_module.R")
source("R/uiFunction.R")

cBy_title <- c("Lineages"= "Grouping",
              "sex"= "Gender",
              "res.1"="Cluster",
              "Ontology_ID"="Ontology ID",
              "orig.ident" = "Zonation") 

pmeta <- readRDS("data/pmeta.rds")
kid_proj<- readRDS("data/kid_proj.rds")
neph_proj<- readRDS("data/neph_proj.rds")
uret_proj <- readRDS("data/uret_proj.rds")

gene_tbl <- readRDS("data/gene_tbl.rds")
expr_data <- readRDS("data/expr_data.rds")
final_expr <- readRDS("data/final_expr.rds")
final_frac <- readRDS("data/final_frac.rds")

otg_tbl <- readRDS("data/otg_tbl.rds")






