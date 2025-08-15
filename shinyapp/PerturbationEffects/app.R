## app.R ##
library(shinydashboard)
library(plotly)
library(networkD3)

ui <- dashboardPage(
  dashboardHeader(title = "Basic dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Drug Effects", tabName = "DrugEffects", icon = icon("dashboard")),
      menuItem("Protein Network", tabName = "Protein Network", icon = icon("th"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "DrugEffects", 
              h2("Drug Effects"),
              fluidRow(
                box(title = "Drug Effect Heatmap", status = "primary", solidHeader = TRUE,
                    plotlyOutput("plot1", height = 600)),
                fluidRow(
                  width = 12,
                  box(title = "Significant Single Drugs", status = "primary", solidHeader = TRUE,
                      tableOutput("table1")),
                  box(title = "Significant Drug Interactions", status = "primary", solidHeader = TRUE,
                      tableOutput("table2"))
                )
              )
      ),
      tabItem(tabName = "Protein Network", 
              h2("Widgets tab content"))
    )

  )
)

server <- function(input, output) {
  # load ordering of drugs (sort drugs with experiments together)
  load('../../data/order.RData')
  nDrugs <- length(drugOrder)

  # load data for Drug Effects
  load("../../results/DrugEffects.RData")

  P <- 1:10

  
  output$plot1 <- renderPlotly({
    # collect min p value of drug effect over proteins
    pvec <- apply(matrix(allPvecs[P, ], nrow = length(P)), 2, min)
    pvec <- pmin(pvec * length(P), 1)
    names(pvec) <- colnames(allPvecs)

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
                  layout(title = prot_names_short[P])
    ht
  })
}

shinyApp(ui, server)