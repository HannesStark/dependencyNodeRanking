# Title     : TODO
# Objective : TODO
# Created by: Hannes
# Created on: 05.06.2019
library(igraph)
library(shiny)
library(jsonlite)
library(tools)
library(DT)
library(visNetwork)

ui <- fluidPage(
    titlePanel("Network Centrality Calculator"),
    sidebarLayout(
        sidebarPanel(
            fileInput(
                "graph_file", "Choose Edgelist or Sharphound Json",
                accept = c(
                    "application/json",
                    "text/comma-separated-values,text/plain",
                    ".json",
                    ".csv",
                    ".txt",
                    "text/json"
                )
            ),
            selectInput("centrality_measure", "Select Your Centrality Measure",
                list(
                    "Degree Centrality",
                    "Closeness Centrality",
                    "Betweenness Centrality",
                    "Eigenvector Centrality",
                    "Dependency Centrality"
                )
            ),
            conditionalPanel(
                condition = "input.centrality_measure == 'Dependency Centrality'",
                numericInput("skSteps", "Sinkhorn-Knopp Steps", 100, min = 1),
                selectInput("perturbationType", "Select Perturbation",
                    list(
                        "No Perturbation",
                        "Diagonal Perturbation",
                        "Full Perturbation"
                    )
                ),
                conditionalPanel(
                    condition = "input.perturbationType != 'No Perturbation'",
                    numericInput("perturbation", "Perturbation by", 0.01)
                )
            )
        ),
        mainPanel(
            conditionalPanel(
                condition = "input.centrality_measure == 'Dependency Centrality'",
                textOutput("matrixType")
            ),
            div(style="display:inline-block; vertical-align: top",
                fluidRow(
                    column(2,
                        #table with centralities
                        tableOutput("table")
                    )
                )
            ),
            div(style="display:inline-block; vertical-align: top",
                visNetworkOutput("network")
            )
        )
    ),
    hr(),
    downloadLink("downloadData", "Documentation")
)


server <- function(input, output) {
    options(shiny.maxRequestSize = 5 * 1024 ^ 3)

    reciprocal <- function(x){
        if (x == 0) return(0)
        else return(1 / x)
    }

    invVec <- function(vec){
        return(sapply(vec, reciprocal))
    }

    countZeros <- function(vec){
        x <- 0
        for (i in vec) {
            if (i == 0)x <- x + 1
        }
        return(x)
    }

    weightedToAdjacency <- function(A){
        for (i in 1 : nrow(A)) {
            for (j in 1 : ncol(A)) {
                if (A[i, j] != 0) {
                    A[i, j] = 1
                }
            }
        }
        return(A)
    }

    #definition from Matrices of O’s and l’s with Total Support RICHARD A. BRUALDI*
    diagonalPerturbation <- function(A, a=1){
        return(A + diag(a, nrow(A), ncol(A)))
    }
    #definition from Matrices of O’s and l’s with Total Support RICHARD A. BRUALDI*
    fullPerturbation <- function(A, a=1){
        return(A + matrix(a, nrow(A), ncol(A)))
    }

    #Kbnig-Egtrvary theorem [12, pp. 55-56]
    termRank <- function(A){
        minNumRows <- NCOL(A) - countZeros(rowSums(A))
        minNumCols <- NROW(A) - countZeros(colSums(A))
        return(min(minNumCols, minNumRows))
    }
    reciprocal <- function(x){
        if (x == 0) return(0)
        else return(1 / x)
    }
    invVec <- function(vec){
        return(sapply(vec, reciprocal))
    }
    #definition from Matrices of O’s and l’s with Total Support RICHARD A. BRUALDI*
    totalSupport <- function(A){
        for (i in 1 : nrow(A)) {
            for (j in 1 : ncol(A)) {
                if (A[i, j] != 0 && termRank(A[- i, - j]) != nrow(A) - 1) {
                    return(FALSE)
                }
            }
        }
        return(TRUE)
    }
    #an adjacency matrix is indecomposable if and only if its graph is a single strongly connected component. Dijkstra, Edsger (1976), A Discipline of Programming, NJ: Prentice Hall, Ch. 25.
    indecomposable <- function(A){
        return(is.connected(graph_from_adjacency_matrix(weightedToAdjacency(A), mode = c("directed")), mode = c("strong")))
    }
    #definition from Matrices of O’s and l’s with Total Support RICHARD A. BRUALDI*
    fullyIndecomposable <- function(A){
        for (i in 1 : nrow(A)) {
            for (j in 1 : ncol(A)) {
                if (termRank(A[- i, - j]) != nrow(A) - 1) {
                    return(FALSE)
                }
            }
        }
        return(TRUE)
    }


    #retruns a matrix with the rownames of A and a column for each x vector of each step of the sinkhornKnopp Method
    sinkhornKnoppSteps <- function(A, steps){
        x = rep(1, nrow(A))
        skVectors <- matrix(x, dimnames = list(rownames(A)))
        for (i in 1 : steps) {
            x <- A %*% invVec(x)
            skVectors <- cbind(skVectors, x)
        }
        return(skVectors)
    }

    getGraph <- reactive({
        req(input$graph_file)
        graph = NULL
        #takes json as Sharphound produces it or Bloodhound exports it and turns it into an Adjacency Matrix that has the Nodes IDs as column and row names
        if (file_ext(input$graph_file) == "json") {
            graph <- read_json(path = input$graph_file$datapath, simplifyVector = TRUE)
            nEdges <- NROW(graph$edges)
            edge_list <- matrix(nrow = 0, ncol = 2)
            for (i in 1 : nEdges) {
                edge_list = rbind(edge_list, c(as.character(graph$edges$source[i]), as.character(graph$edges$target[i])))
            }
            graph = graph_from_edgelist(edge_list)
        } else {
            df = read.table(input$graph_file$datapath)
            graph = graph_from_data_frame(df, directed = TRUE, vertices = NULL)
        }
        graph
    })

    #renders a visualisation of the uploaded graph
    output$network <- renderVisNetwork({
        graph = getGraph()
        visNetGraph = toVisNetworkData(graph)
        nodes = visNetGraph$nodes
        edges = visNetGraph$edges
        visNetwork(nodes, edges)
    })

    perturbBySettings <- function(A){
        perturbationType = input$perturbationType
        perturbBy = input$perturbation
        if (perturbationType == "Diagonal Perturbation") {
            return(diagonalPerturbation(A, perturbBy))
        } else if (perturbationType == "Full Perturbation") {
            return(fullPerturbation(A, perturbBy))
        }
        return(A)
    }

    #renders the information about total support indecomposability and full indecomposability if dependency centrality is selected
    output$matrixType <- renderText({
        graph = getGraph()
        unperturbed = as_adjacency_matrix(graph, sparse=FALSE)
        A = perturbBySettings(unperturbed)
        if(fullyIndecomposable(A)) {
            return("The given network is fully indecomposable. This means there is exactly one solution and the SK-method will converge")
        } else if (totalSupport(A)){
            return("The given network has total support but is not indecomposable. This means there is at least one solution and the SK-method will converge. Consider full perturbation for full indecomposability")
        } else if (indecomposable(A)){
            return("The given network is indecomposable but has no total support. There is no solution. Consider diagonal Perturbation for full indecomposability")
        } else {
            return("The given network is decomposable and has no total support. The SK-method did not converge. Consider full perturbation for full indecomposability or diagonal perturbation for total support")
        }
    })

    # calls the right functions regarding the chosen centrality_measure and returns a table sorted by the weights with the node names as row names
    output$table <- renderTable({
        graph = getGraph()
        measure = input$centrality_measure
        if (measure == "Degree Centrality") {
            table = degree(graph)
        } else if (measure == "Closeness Centrality") {
            table = closeness(graph)
        } else if (measure == "Betweenness Centrality") {
            table = betweenness(graph)
        } else if (measure == "Eigenvector Centrality") {
            table = eigen_centrality(graph)$vector
        } else if (measure == "Dependency Centrality") {
            mat = as_adjacency_matrix(graph, sparse = FALSE)
            numberOfSKSteps = input$skSteps
            allSteps = sinkhornKnoppSteps(perturbBySettings(mat), numberOfSKSteps)
            table = allSteps[, numberOfSKSteps + 1]
        }
        sort(table, decreasing = TRUE)

    },
    rownames = TRUE)

    output$downloadData <- downloadHandler(
    filename = "NetworkCentralityCalculator.pdf",
    content = function(file) {
      file.copy("NetworkCentralityCalculator.pdf", file)
    }
  )
}

shinyApp(ui, server)
