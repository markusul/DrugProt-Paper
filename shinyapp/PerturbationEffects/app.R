## app.R ##
library(shinydashboard)
library(plotly)
library(networkD3)
library(HiveR)
require("grid")

ui <- dashboardPage(
  dashboardHeader(title = "Drug Perturbations"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Settings", tabName = "Settings", icon = icon("cog")),
      menuItem("Drug Effects", tabName = "DrugEffects", icon = icon("dashboard")),
      menuItem("Protein Network", tabName = "ProteinNetwork", icon = icon("th"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "DrugEffects", 
              h2("Drug Effects"),
              fluidRow(
                box(title = "Drug Effect Heatmap", status = "primary", solidHeader = TRUE,
                    plotlyOutput("plot1", height = 800), width = 8),
                box(title = "Significant Single Drugs", status = "primary", solidHeader = TRUE,
                    tableOutput("table1"), width = 2),
                box(title = "Significant Drug Interactions", status = "primary", solidHeader = TRUE,
                    tableOutput("table2"), width = 2)
              )
      ),
      tabItem(tabName = "ProteinNetwork", 
              h2("Protein Network"),
              fluidRow(
                box(title = "Summary Graph", status = "primary", solidHeader = TRUE,
                    forceNetworkOutput("SummaryGraph", height = 800), width = 12)
              ),
              fluidRow(
                box(title = "Temporal Graph", status = "primary", solidHeader = TRUE,
                    forceNetworkOutput("TemporalGraph", height = 800), width = 12)
              ),
              fluidRow(
                box(title = "Protein Interaction Network", status = "primary", solidHeader = TRUE,
                    plotOutput("HivePlot", height = 400), width = 8), 
                box(title = "Relevant Proteins", status = "primary", solidHeader = TRUE,
                    fluidRow(column(tableOutput("numProteinEffects"), width = 7), 
                             column(tableOutput("ProteinEffects"), width = 5)), width = 4)
              )
      )
    )

  )
)

server <- function(input, output) {
  print("Start Server")
  
  # load ordering of drugs (sort drugs with experiments together)
  load('../../data/order.RData')
  nDrugs <- length(drugOrder)

  # load data for Drug Effects
  load("../../results/DrugEffects.RData")

  P_selection <- 1:15
  print(prot_names_short[P_selection])

  # collect min p value of drug effect over proteins
  pvec <- apply(matrix(allPvecs[P_selection, ], nrow = length(P_selection)), 2, min)
  pvec <- pmin(pvec * length(P_selection), 1)
  names(pvec) <- colnames(allPvecs)
  
  # Protein Network tab content
  print("Load Protein Network Data")
  load("../../results/proteinNetwork.RData")
  
  # significance level for protein network
  alpha <- 0.05
  alpha <- alpha / length(P_selection) / 3
  
  # transform p value to links using alpha
  Links_all <- lapply(Net, function(net){
    res <- apply(net[[1]], 2, function(pval) which(pval < alpha))
    links <- NULL
    if(length(res) == 0) return(NULL)
    for(i in 1:ncol(net[[1]])){
      if(length(res[[i]]) != 0)
        links <- rbind(links, data.frame(res[[i]], which(net[[2]][i] == prot_names_short)))
    }
    rownames(links) <- 1:nrow(links)
    colnames(links) <- c('source', 'target')
    links
  })

  print("Prepare Summary Graph")
  
  Links_sum <- do.call(rbind, Links_all)
  Links_sum$source <- Links_sum$source - 1
  Links_sum$target <- Links_sum$target - 1
  Links_sum$value <- 1
  Nodes_sum <- data.frame(name = prot_names_short, group = 1, size = 1)
  P.ind <- P_selection - 1

  #add link to itself
  #Links_self <- data.frame(source = P.ind, target = P.ind, value = 1)
  #Links_sum <- rbind(Links_sum, Links_self)

  rel.Links <- Links_sum$source %in% P.ind | Links_sum$target %in% P.ind
  print(sum(rel.Links))

  if(sum(rel.Links) == 0){
    print("No relevant links found, using dummy links")
    Links_sum <- data.frame(source = P.ind, target = P.ind, value = 1)
  }else {
    Links_sum <- Links_sum[rel.Links, ]
  }
  
  rel.Nodes <- sort(unique(unlist(Links_sum[, c('source', 'target')])))
  print(rel.Nodes)
  Nodes_sum$group[P_selection] <- 2
  Nodes_sum <- Nodes_sum[rel.Nodes+1, ]
  
  #reorganize link index
  for(i in 1:length(rel.Nodes)){
    Links_sum[Links_sum == rel.Nodes[i]] <- i - 1
  }
  Links_sum$value <- 1
  
  # temporal graph
  print("Prepare Temporal Graph")
  expTimes <- c(6, 24, 48)
  nodenames <- c(paste(prot_names_short, expTimes[1], sep = '_'), 
                 paste(prot_names_short, expTimes[2], sep = '_'), 
                 paste(prot_names_short, expTimes[3], sep = '_'))
  nodegroups <- rep(paste0(expTimes, 'h'), each = length(prot_names_short))
  
  Links_temp <- lapply(1:2, function(t){
    links <- Links_all[[t]]
    links[, 1] <- links[, 1] + (t-1) * length(prot_names_short)
    links[, 2] <- links[, 2] + (t) * length(prot_names_short)
    links
  })
  Links_temp <- do.call(rbind, Links_temp)
  used_nodes <- sort(unique(unname(unlist(Links_temp))))
  
  Links_temp$source <- Links_temp$source - 1
  Links_temp$target <- Links_temp$target - 1
  Links_temp$value <- 1
  Nodes_temp <- data.frame(name = nodenames, group = nodegroups, size = 1)
  Nodes_temp$size[used_nodes] <- 100
  Nodes_temp$radius <- as.numeric(rep(1:length(prot_names_short), 3))
  
  P.ind <- unlist(lapply(P_selection, function(p) p + 0:2 * length(prot_names_short))) - 1
  #add link to itself
  #Links_self <- data.frame(source = P.ind, target = P.ind, value = 1)
  #Links_temp <- rbind(Links_temp, Links_self)

  rel.Links <- Links_temp$source %in% P.ind | Links_temp$target %in% P.ind
  Links_temp <- Links_temp[rel.Links, ]
  
  rel.Nodes <- sort(unique(unlist(Links_temp[, c('source', 'target')])))
  Nodes_temp <- Nodes_temp[rel.Nodes+1, ]
  
  #reorganize link index
  for(i in 1:length(rel.Nodes)){
    Links_temp[Links_temp == rel.Nodes[i]] <- i - 1
  }
  Links_temp$value <- 1
  
  print("prepare for HivePlot")

  edges <- Links_temp
  edges[, 1:2] <- edges[, 1:2] + 1

  #edges[(nrow(edges)-length(P_selection)*3):nrow(edges), 2] <- c(2:(length(P_selection)*3), 1)
  names(edges) <- c("id1", "id2", "weight")
  # color
  row.names(edges) <- NULL
  edges$id1 <- as.integer(edges$id1)
  edges$id2 <- as.integer(edges$id2)
  edges$color <- "black"
  edges$weight <- 0.1
  
  #edges$color[(nrow(edges)-length(P_selection)*3):nrow(edges)] <- "white"

  nodes <- Nodes_temp
  nodes$size <- 0.01
  names(nodes) <- c("lab", "axis", "size", "radius")
  #nodes$axis <- as.integer(as.numeric(nodes$axis == "24h") + 1)
  nodes$axis[nodes$axis == "6h"] <- 2
  nodes$axis[nodes$axis == "24h"] <- 1
  nodes$axis[nodes$axis == "48h"] <- 3
  nodes$axis <- as.integer(nodes$axis)
  nodes$id <- 1:nrow(nodes)
  nodes$radius <- nodes$radius * 3
  nodes$color <- "black"
  #nodes$color <- c("black", "red", "blue", "green", "violet", "yellow", "darkgreen")[clusters]
  #rep(c("black", "red", "blue", "green", "violet", "yellow", "darkgreen"), times = as.numeric(table(clusters)))
  
  
  HEC <- list()
  HEC$nodes <- nodes
  HEC$edges <- edges
  HEC$type <- "2D"
  HEC$desc <- "HairEyeColor data set"
  HEC$axis.cols <- c("grey", "grey")
  class(HEC) <- "HivePlotData" 

  print("Data loaded and processed")

  output$plot1 <- renderPlotly({
    pMat <- matrix(NA, nrow = nDrugs, ncol = nDrugs)
    rownames(pMat) <- colnames(pMat) <- names(pvec)[1:nDrugs]
    for(l in names(pvec)){
      drugs <- strsplit(l, ":")[[1]]
      if(length(drugs) == 1) drugs <- c(drugs, drugs)
      pMat[drugs[1], drugs[2]] <- pvec[l]
      pMat[drugs[2], drugs[1]] <- pvec[l]
    }

    pMat <- as.matrix(pMat)
    pMat[is.na(pMat)] <- 2
    pMat <- pMat[drugOrder, drugOrder]

    ht <- plot_ly(z = pMat, x = colnames(pMat), y = colnames(pMat), 
                  type = "heatmap", colors = "Greys") %>%
                  layout(title = prot_names_short[P_selection])
    ht
  })

  output$table1 <- renderUI({
    singleInd <- sapply(names(pvec), function(l) length(strsplit(l, ":")[[1]]) == 1)
    pvec_single <- pvec[singleInd]
    pvec_single <- pvec_single[order(pvec_single)]
    if(length(pvec_single) == 0) return(NULL)
    df <- data.frame(Drug = names(pvec_single), PValue = pvec_single, stringsAsFactors = FALSE)
    htmlTable <- paste0(
      '<table class="table table-bordered"><thead><tr><th>Drug</th><th>PValue</th></tr></thead><tbody>',
      paste(
        sapply(1:nrow(df), function(i) {
          pv <- df$PValue[i]
            color <- ifelse(pv < 0.05, ' style="background-color: #add8e6;"', '')
          paste0('<tr><td', color, '>', df$Drug[i], '</td><td', color, '>', format(pv, digits = 4), '</td></tr>')
        }),
        collapse = ""
      ),
      '</tbody></table>'
    )
    HTML(htmlTable)
  })

  output$table2 <- renderUI({
    interactionInd <- sapply(names(pvec), function(l) length(strsplit(l, ":")[[1]]) > 1)
    pvec_interaction <- pvec[interactionInd]
    pvec_interaction <- pvec_interaction[order(pvec_interaction)]
    if(length(pvec_interaction) == 0) return(NULL)
    df <- data.frame(DrugInteraction = names(pvec_interaction), PValue = pvec_interaction, stringsAsFactors = FALSE)
    htmlTable <- paste0(
      '<table class="table table-bordered"><thead><tr><th>Drug Interaction</th><th>PValue</th></tr></thead><tbody>',
      paste(
        sapply(1:nrow(df), function(i) {
          pv <- df$PValue[i]
            color <- ifelse(pv < 0.05, ' style="background-color: #add8e6;"', '')
          paste0('<tr><td', color, '>', df$DrugInteraction[i], '</td><td', color, '>', format(pv, digits = 4), '</td></tr>')
        }),
        collapse = ""
      ),
      '</tbody></table>'
    )
    HTML(htmlTable)
  })

  output$numProteinEffects <- renderTable({
    P_sel_inNodes <- which(Nodes_sum$name %in% prot_names_short[P_selection])
    df <- sapply(P_sel_inNodes - 1, function(i) {
      Protein <- Nodes_sum$name[i + 1]
      Children <- sum(Links_sum$source == i)#-1
      Parents <- sum(Links_sum$target == i)#-1
      df <- data.frame(Protein = Protein, numParents = Parents, num.Children = Children)
      return(df)
    })
    t(df)
  })

  output$ProteinEffects <- renderTable({
    P_sel_inNodes <- which(Nodes_sum$name %in% prot_names_short[P_selection])
    fam <- lapply(P_sel_inNodes - 1, function(i) {
      Children <- Nodes_sum$name[Links_sum$target[Links_sum$source == i & Links_sum$target != i] + 1]
      Parents <- Nodes_sum$name[Links_sum$source[Links_sum$target == i & Links_sum$source != i] + 1]
      return(list(Children = Children, Parents = Parents))
    })

    # vector of all children and parents
    allChildren <- unique(unlist(sapply(fam, "[[", 1)))
    allParents <- unique(unlist(sapply(fam, "[[", 2)))

    # add dummies to smaller vectors
    maxLength <- max(length(allChildren), length(allParents))
    allChildren <- c(allChildren, rep("", maxLength - length(allChildren)))
    allParents <- c(allParents, rep("", maxLength - length(allParents)))
    df <- data.frame(Parents = allParents,
                     Children = allChildren,  
                     stringsAsFactors = FALSE)
    df
  })

  output$SummaryGraph <- renderForceNetwork({
    fN <- forceNetwork(Links = Links_sum, Nodes = Nodes_sum,
                 Source = "source", Target = "target",
                 Value = "value", NodeID = "name",
                 Group = "group", opacity = 0.99, 
                 arrows = T, zoom = T, charge = -5,
                 opacityNoHover = TRUE,
                 colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
    fN
  })

  output$HivePlot <- renderPlot({
    plotHive(HEC, ch = 0.001, bkgnd = "white", 
             axLabs = c("24h", "6h", "48h"), 
             axLab.gpar = gpar(col = "black", fontsize = 24))
  })

  output$TemporalGraph <- renderForceNetwork({
    fN <- forceNetwork(Links = Links_temp, Nodes = Nodes_temp,
             Source = "source", Target = "target",
             Value = "value", NodeID = "name",
             Group = "group", opacity = 0.99,# Nodesize = 3,
             arrows = T, zoom = T, legend=T, charge = -5,
             opacityNoHover = TRUE,
             colourScale = JS("d3.scaleOrdinal(d3.schemeCategory10);"))
    fN
  })

}

shinyApp(ui, server)