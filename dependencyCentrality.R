# Title     : Dependency centrality
# Objective : Read bloodhound json files and use the sinkhornKnopp method on it to determine a graphs most powerful node according to dependency centrality
# Created by: Hannes Stärk
# Created on: 27.01.2019

library(jsonlite)
library(igraph)

kite = matrix(
c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
0, 1, 0, 1, 1, 0, 0, 0, 0, 0,
0, 0, 1, 0, 1, 1, 1, 0, 1, 0,
0, 0, 1, 1, 0, 0, 1, 1, 0, 1,
0, 0, 0, 1, 0, 0, 1, 0, 1, 0,
0, 0, 0, 1, 1, 1, 0, 1, 1, 1,
0, 0, 0, 0, 1, 0, 1, 0, 0, 1,
0, 0, 0, 1, 0, 1, 1, 0, 0, 1,
0, 0, 0, 0, 1, 0, 1, 1, 1, 0), nrow = 10, ncol = 10, byrow = TRUE)

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
countZeros <- function(vec){
    x <- 0
    for (i in vec) {
        if (i == 0)x <- x + 1
    }
    return(x)
}
#definition from Matrices of O’s and l’s with Total Support RICHARD A. BRUALDI*
doublyStochastic <- function(A){
    return(sum(colSums(A)) == NROW(A) && sum(rowSums(A)) == NCOL(A))
}
#if an entry of the matrix is not 0 it will be replaced with 1 and mirrored
symmetricMatrix <- function(A){
    for (i in 1 : nrow(A)) {
        for (j in 1 : ncol(A)) {
            entry <- A[i, j]
            if (entry != 0) {
                A[j, i] = entry
            }
        }
    }
    return(A)
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
            if (A[i, j] == 1 && termRank(A[- i, - j]) != nrow(A) - 1) {
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
#definition from Matrices of O’s and l’s with Total Support RICHARD A. BRUALDI*
diagonalPerturbation <- function(A, a=1){
    return(A + diag(a, nrow(A), ncol(A)))
}
#definition from Matrices of O’s and l’s with Total Support RICHARD A. BRUALDI*
fullPerturbation <- function(A, a=1){
    return(A + matrix(a, nrow(A), ncol(A)))
}
#takes json as Sharphound produces it or Bloodhound exports it and turns it into an Adjacency Matrix that has the Nodes IDs as column and row names
jsonToAdjacencyMatrix <- function(path){
    graph <- read_json(path = path, simplifyVector = TRUE)
    nNodes <- NROW(graph$nodes$id)
    resultMatrix <- matrix(0, nrow = nNodes, ncol = nNodes, dimnames = list(graph$nodes$id, graph$nodes$id))
    nEdges <- NROW(graph$edges)
    for (i in 1 : nEdges)
    {
        resultMatrix[toString(graph$edges$source[i]), toString(graph$edges$target[i])] <- 1
        #resultMatrix[toString(graph$edges$target[i]), toString(graph$edges$source[i])] <- 1 #if we assume that the graph is undirected a directed edge in bloodhound will represent both directions
    }
    return(resultMatrix)
}
#retruns a matrix with the rownames of A and a column for each x vector of each step of the sinkhornKnopp Method
sinkhornKnoppSteps <- function(A, steps){
    x=rep(1, nrow(A))
    skVectors <- matrix(x, dimnames = list(rownames(A)))
    for (i in 1 : steps) {
        x <- A %*% invVec(x)
        skVectors <- cbind(skVectors, x)
    }
    return(skVectors)
}

#THE FOLDER "convergences" HAS TO EXIST FOR THIS FUNCTION TO RUN. The plots of the convergences are placed in that folder.
#takes an Adjacency Matrix and an Integer. It does numberOfSKSteps iterations of the Sinkhorn-Knopp method and generates:
#1: a plot of the graph
#2: a plot of the values for each stop for every node
#3: plots of the last even and the last odd vector
#4: prints if the matrix is indecomposable fully indecomposable or has total Support
#5: prints every vector for each step as a matrix
#6: prints the most powerful node
analyzeMatrix <- function(A, numberOfSKSteps=100, sideEffects=TRUE){
    if (is.null(rownames(A))) {
        rownames(A) <- 1 : nrow(A)
    }
    erg <- sinkhornKnoppSteps(A, numberOfSKSteps)
    bestNodeName <- rownames(erg)[which.max(erg[, numberOfSKSteps + 1])]
    if (sideEffects) {
        pdf("graph.pdf")
        asGraph = graph_from_adjacency_matrix(A, mode = c("directed"))
        plot(asGraph)
        print(erg)
        print(A)
        print(paste("is indecomposable: ", indecomposable(A)))
        print(paste("has total support: ", totalSupport(A)))
        print(paste("is fully indecomposable: ", fullyIndecomposable(A)))
        if (totalSupport(A)) {
            print(paste("most powerful node is: ", bestNodeName))
        }else {
            print("Sinkhorn-Knopp method did not converge")
        }
        evenResultVector <- erg[, numberOfSKSteps + 1]
        oddResultVector <- erg[, numberOfSKSteps]
        if (numberOfSKSteps %% 2 != 0) {
            evenResultVector <- erg[, numberOfSKSteps]
            oddResultVector <- erg[, numberOfSKSteps + 1]
        }
        for (rowname in rownames(A)) {
            pdf(paste("convergences\\convergence", rowname, ".pdf"))
            plot(1 : (numberOfSKSteps + 1), erg[rowname,], main = paste("Convergence of Node: ", rowname), xlab = "Number of Steps", ylab = "Value of Node")
        }
        pdf("evenResultVector.pdf")
        plot(rownames(A), evenResultVector,main = "Last even Vector" , xlab = "Name of Node", ylab = "Value of Node")
        pdf("oddResultVector.pdf")
        plot(rownames(A), oddResultVector, main = "Last odd Vector" , xlab = "Name of Node", ylab = "Value of Node")
    }
    return(list("skVectors" = erg, "bestNodeName" = bestNodeName))
}

#df = read.table("C:\\Users\\Hannes\\projects\\graph\\test.txt")
#resGraph = graph_from_data_frame(df, directed = TRUE, vertices = NULL)


#jsonToAdjacencyMatrix nimmt einen pfad entgegen an dem die json Ausgabedatei liegt und gibt die Adjazenzmatrix zu dem graphen zurück
#exMatrix = jsonToAdjacencyMatrix("F:\\projects\\graph\\testGraphs\\testGraph.json")

#damit diese Funktion funktioniert muss der ordner "convergences" im selben Verzeichnis wie das Script existieren.
analyzeMatrix(diagonalPerturbation(kite, 0.1), 100)


