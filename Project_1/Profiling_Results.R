
library(lineprof)
library(shiny)

benchmark <- readRDS("Profiling_Benchmarkfirst.rds")

shine(benchmark)