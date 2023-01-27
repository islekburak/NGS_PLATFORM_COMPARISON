library(dplyr)
library(networkD3)
library(htmlwidgets)
library(htmltools)
library(manipulateWidget)

nodes = read.csv("/Users/islekbro/Desktop/poster/clingen2seq/nodes.csv", sep=";" , header=TRUE)
links = read.csv("/Users/islekbro/Desktop/poster/clingen2seq/links.csv", sep=";" , header=TRUE)
my_color <- 'd3.scaleOrdinal() .domain(["a","a1","a2","a3","a4","b","b1","b2","b3","b4"]) .range(["#FF7959","#FFB3A3","#6FA0D7","#A0E1B6","#70C775","#ffb6a5","#ffd6cd","#b1cbe9","#cbefd8","#b1e1b4"])'
sankey <- sankeyNetwork(Links = links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              fontSize= 12, nodeWidth = 30, colourScale=my_color, NodeGroup = "group",LinkGroup="group", iterations=0, sinksRight=FALSE)
#save the widget
#library(htmlwidgets)
#saveWidget(p, file=paste0( getwd(), "/HtmlWidget/sankeyColor3.html"))

#NodeGroup yanına LinkGroup="group" eklersen flowlar FZD nodelarının rengine dönüyor.

# create left-right label
leftTx = tags$div(
  style="max-width: 100vw; height: 100%;
                 display: flex; align-items: center;
                 justify-content: center;",
  tags$p("CLINGEN"))

rightTx = tags$p("GENOMIZE-SEQ",
                 style="max-width:100vw; height: 100%;
                 display: flex; align-items: center;
                 justify-content: center")

# final output with label. combineWidget
# can be use by manipulateWidget.
sankey2 <- combineWidgets(sankey,
                         title = tags$h1("Comparison of Pathogenicity Prediction",
                                         style="text-align:center;color:#000000;"),
                         leftCol = leftTx,
                         rightCol = rightTx,
                         nrow = 1)
# view graph
sankey2