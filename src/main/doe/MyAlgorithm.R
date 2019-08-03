#help: MyAlgorithm method for ...
#author: ...
#ref: ...
#tags: ...
#options: ytarget='0.0';ytol='3.e-8';xtol='1.e-8';xsample_size='10';max_iterations='100'
#input: x=list(min=0,max=1)
#output: y=0.01

MyAlgorithm <- function(options) {
    myAlgorithm = new.env()

    myAlgorithm$ytol <- as.numeric(options$ytol)
    myAlgorithm$xtol <- as.numeric(options$xtol)
    myAlgorithm$ytarget <- as.numeric(options$ytarget)
    myAlgorithm$max_iterations <- as.integer(options$max_iterations)
    myAlgorithm$xsample_size <- as.integer(options$xsample_size)
    myAlgorithm$i = NA

    return(myAlgorithm)
}

#' first design building.
#' @param input variables description (min/max, properties, ...)
#' @param output values of interest description
getInitialDesign <- function(myAlgorithm, input, output) {
    myAlgorithm$i <- 0
    myAlgorithm$input <- input
    x = matrix(runif(myAlgorithm$xsample_size * length(input)),ncol=length(input))
    names(x) <- names(input)
    return(from01(x,myAlgorithm$input))
}

## iterated design building.
## @param X data frame of current doe variables
## @param Y data frame of current results
## @return data frame or matrix of next doe step
getNextDesign <- function(myAlgorithm, X, Y) {
    names(X) = names(myAlgorithm$input)
    X = to01(X,myAlgorithm$input)
    Y = matrix(Y,ncol=1) - myAlgorithm$ytarget

    if (myAlgorithm$i >= myAlgorithm$max_iterations) {
        return(NULL)
    }

    Xnext = matrix(runif(myAlgorithm$xsample_size * length(input)),ncol=length(input))
    names(Xnext) <- names(myAlgorithm$input)
    return(from01(Xnext,myAlgorithm$input))
}

## final analysis. Return HTML string
## @param X data frame of doe variables
## @param Y data frame of  results
## @return HTML string of analysis
displayResults <- function(myAlgorithm, X, Y) {
    myAlgorithm$files <- paste("result", myAlgorithm$i, ".png", sep = "")
    height <- 500
    width <- 500

    png(file = myAlgorithm$files,        height = height,        width = width)
    pairs(cbind(X,Y)
    dev.off()

    html <- paste0(' <HTML name="SomeInfo">Here is some information in the end','<br/>',
            '<img src="',  myAlgorithm$files,  '" width="', width, '" height="', height, '"/></HTML>',collapse=';')

    info <- paste0('<info>',mean(Y),'</info>')

    return(paste0(html,info))
}

displayResultsTmp <- displayResults

from01 = function(X, inp) {
    for (i in 1:ncol(X)) {
        namei = names(X)[i]
        X[,i] = X[,i] * (inp[[ namei ]]$max-inp[[ namei ]]$min) + inp[[ namei ]]$min
    }
    return(X)
}

to01 = function(X, inp) {
    for (i in 1:ncol(X)) {
        namei = names(X)[i]
        X[,i] = (X[,i] - inp[[ namei ]]$min) / (inp[[ namei ]]$max-inp[[ namei ]]$min)
    }
    return(X)
}
