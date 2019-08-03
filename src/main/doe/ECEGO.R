#title: ECEGO
#help: Efficient Global Optimization (EGO) algorithm with equality constraints.
#tags: optimization; sparse; contraints
#author: yann.richet@irsn.fr; DiceKriging authors
#require: DiceDesign; DiceKriging; DiceView; pso; jsonlite; randtoolbox
#options: search_ymin='true'; initBatchSize='4'; batchSize='4'; iterations='10'; initBatchBounds='true'; trend='y~1'; covtype='matern3_2'; knots='0'; liar='upper95'; trend_constr='y~1'; covtype_constr='matern3_2'; liar_constr='upper95';  seed='1'
#options.help: search_ymin=minimization or maximisation; initBatchSize=Initial batch size; batchSize=iterations batch size; iterations=number of iterations; initBatchBounds=add input variables bounding values (2^d combinations); trend=(Universal) kriging trend; covtype=Kriging covariance kernel; knots=number of non-stationary points for each Xi; liar=liar value for in-batch loop (when batchsize>1); seed=random seed
#input: x=list(min=0,max=1)
#output: y=0.99

ECEGO <- function(options) {

    library(DiceDesign)
    library(DiceKriging)
    library(DiceView)
    library(pso)
    library(jsonlite)

    ego = new.env()
    ego$i = 0

    ega$search_ymin <- as.logical(options$search_ymin)
    ego$initBatchSize <- as.integer(options$initBatchSize)
    ego$batchSize <- as.integer(options$batchSize)
    ego$iterations <- as.integer(options$iterations)
    ego$initBatchBounds <- as.logical(options$initBatchBounds)
    ego$trend <- as.formula(options$trend)
    ego$covtype <- as.character(options$covtype)
    ego$liar <- as.character(options$liar)
    ego$knots <- as.integer(unlist(strsplit(options$knots,",")))

    ego$seed <- as.integer(options$seed)

    ego$trend_constr <- as.formula(options$trend_constr)
    ego$covtype_constr <- as.character(options$covtype_constr)
    ego$liar_constr <- as.character(options$liar_constr)

    return(ego)
}

getInitialDesign <- function(algorithm, input, output) {
    algorithm$input <- input
    algorithm$output <- output

    set.seed(algorithm$seed)

    d = length(input)
    lhs <- lhsDesign(n = algorithm$initBatchSize, dimension = d,seed=algorithm$seed)$design
    if (isTRUE(algorithm$initBatchBounds)) {
        e = c(0, 1)
        id = 1
        while (id < d) {
            e = rbind(cbind(e, 0), cbind(e, 1))
            id = id + 1
        }
        Xinit = rbind(as.matrix(e), as.matrix(lhs))
    } else {
        Xinit = as.matrix(lhs)
    }

    for (i in 1:d)
        Xinit[,i] = Xinit[,i] * (input[[i]]$max-input[[i]]$min) + input[[i]]$min
    colnames(Xinit) <- names(input)

    return(Xinit)
}

getNextDesign <- function(algorithm, X, Y) {
    if (algorithm$i >= algorithm$iterations) return()

    set.seed(algorithm$seed)

    d = dim(X)[2]
    if (dim(Y)[2] > 2) {
        algorithm$noise.var <- as.array(Y[,2])^2
        y_constr = Y[,3]
    } else {
        algorithm$noise.var <- NULL
        y_constr = matrix(Y[,2],ncol=1)
    }

    if (isTRUE(algorithm$search_ymin)) {
        y = Y[, 1]
    } else {
        y = -Y[, 1]
    }
    y = matrix(y,ncol=1)

    # heurisitc for lower bound of theta : max(1e-6, 0.1 * dX[which.max(dy/rowSums(dX))])
    dX = apply(FUN = dist, X, MARGIN = 2)
    dy = apply(FUN = dist, y, MARGIN = 2)
    dy_constr = apply(FUN = dist, y_constr, MARGIN = 2)

    # define stantionary-changing points
    all_knots <- generate_knots(knots.number = algorithm$knots, d = d, lower = sapply(algorithm$input, "[[", "min"), upper = sapply(algorithm$input, "[[", "max"))

    algorithm$model <- km(algorithm$trend, optim.method = "BFGS",
                        covtype = algorithm$covtype,
                        design = X, response = y, noise.var = algorithm$noise.var,
                        lower = pmax(1e-06, 0.1 * dX[which.max(dy/rowSums(dX)),]),
                        control = list(trace = FALSE),
                        scaling = is.list(all_knots), knots = all_knots)

	  algorithm$model_constr <- km(algorithm$trend_constr, optim.method = "BFGS",
                        covtype = algorithm$covtype_constr,
                        design = X, response = y_constr, noise.var = algorithm$noise.var_constr,
                        lower = pmax(1e-06, 0.1 * dX[which.max(dy_constr/rowSums(dX)),]),
                        control = list(trace = FALSE),
                        scaling = is.list(all_knots), knots = all_knots)

    oEGO <- max_qEITSEE(model = algorithm$model, model.constr=algorithm$model_constr, npoints = algorithm$batchSize,
                    L = algorithm$liar,L.constr = algorithm$liar_constr,
                    lower = sapply(algorithm$input, "[[", "min"),
                    upper = sapply(algorithm$input, "[[", "max"),
                    control = list(trace = FALSE, init = algorithm$i==0))

    if (is.null(oEGO))
        return()

    Xnext <- oEGO$par
    algorithm$i <- algorithm$i + 1

    Xnext = as.matrix(Xnext)
    colnames(Xnext) <- names(algorithm$input)
    return(Xnext)
}


displayResults <- function(algorithm, X, Y) {
    algorithm$files <- paste("view_", algorithm$i,".png", sep = "")
    resolution <- 600

    if (dim(Y)[2] == 2) {
        noise.var <- as.array(Y[, 2])^2
        yname = paste0("N(", colnames(Y)[1], ",", colnames(Y)[2],")")
    } else {
        noise.var <- NULL
        yname = colnames(Y)
    }

    if (isTRUE(algorithm$search_ymin)) {
        y = Y[, 1]
        m = min(Y[, 1])
        x = as.matrix(X)[which(Y[, 1] == m), ]
        html = paste0(sep = "<br/>",
                     paste0("<HTML>minimum is ", m),
                     paste0(sep = "",
                           "found at <br/>",
                           paste0(collapse = "<br/>",paste(sep = "= ", names(X), x)),
                           "<br/><img src='", algorithm$files,
                           "' width='", resolution, "' height='", resolution, "'/></HTML>"))
        html = paste0(html,"<min>",m,"</min><argmin>",toJSON(x),"</argmin>")
    } else {
        y = -Y[, 1]
        m = max(Y[, 1])
        x = as.matrix(X)[which(Y[, 1] == m), ]
        html = paste0(sep = "<br/>",
                     paste0("<HTML>maximum is ", m),
                     paste0(sep = "",
                           "found at <br/>",
                           paste0(collapse = "<br/>", paste(sep = "=", names(X), x)),
                           "<br/><img src='",  algorithm$files,
                           "' width='", resolution, "' height='",  resolution, "'/></HTML>"))
        html = paste0(html,"<max>",m,"</max><argmax>",toJSON(x),"</argmax>")
    }

    if (!exists("model",envir = algorithm)) {
        png(file = algorithm$files, bg = "transparent", height = resolution, width = resolution)
        try(pairs(cbind(X,Y)))
        dev.off()
        return(html)
    }
    print(paste0(html,collapse=';'))

    png(file = algorithm$files, bg = "transparent", height = resolution, width = resolution)
    try(sectionview.km(algorithm$model, center = x, Xname = colnames(X), yname = yname, yscale = ifelse(algorithm$search_ymin,1,-1)))
    dev.off()

    #if (algorithm$i == algorithm$iterations) {
        html = paste0(html,"<data_json>",toJSON(as.data.frame(cbind(X,Y)),dataframe = "columns"),"</data_json>")

        lower = sapply(algorithm$input, "[[", "min")
        upper = sapply(algorithm$input, "[[", "max")
        n = 1000
        set.seed(123) # to get the same points for evaluating model
        Xm = matrix(lower,nrow=n,ncol=length(lower),byrow = T) + matrix(upper-lower,nrow=n,ncol=length(lower),byrow = T) * matrix(runif(n*length(lower)),nrow=n,ncol=length(lower))
        colnames(Xm) <- colnames(X)
        Ym = predict(algorithm$model,newdata = Xm,type = "UK",cov.compute = F, low.memory = T)
        Ym = cbind(ifelse(algorithm$search_ymin,1,-1)*Ym$mean,Ym$sd)
        colnames(Ym) <- c(colnames(Y),paste0("sd_",colnames(Y)))[1:2]

        html = paste0(html,"<model_json>",toJSON(as.data.frame(cbind(Xm,Ym)),dataframe = "columns"),"</model_json>")
    #}

    return(paste0(html,collapse=';'))
}

displayResultsTmp <- displayResults

################### Algorithm dependencies ###################

distXmin <- function (x, Xmin){
    return(min(sqrt(rowSums((Xmin - matrix(x, nrow = nrow(Xmin), ncol = ncol(Xmin), byrow = TRUE))^2))))
}

#' @test X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); EI(runif(100),kmi)
#' @test X=matrix(runif(10),ncol=1); y=-sin(pi*X); kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)EI(x,kmi),dim=1)
#' @test X=matrix(runif(10),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)EI(x,kmi),dim=2)
EI <- function (x, model, plugin = NULL){
    if (is.null(plugin)) {
        if (model@noise.flag)
            plugin <- min(model@y - 2 * sqrt(model@noise.var))
        else plugin <- min(model@y)
    }
    m <- plugin
    if (!is.matrix(x))
        x <- matrix(x, ncol = model@d)
    d <- ncol(x)
    if (d != model@d)
        stop("x does not have the right number of columns (", d, " instead of ", model@d, ")")

    newdata <- x
    colnames(newdata) = colnames(model@X)
    predx <- predict.km(object = model, newdata = newdata, type = "UK", checkNames = FALSE)
    kriging.mean <- predx$mean
    kriging.sd <- predx$sd
    xcr <- (m - kriging.mean)/kriging.sd
    xcr.prob <- pnorm(xcr)
    xcr.dens <- dnorm(xcr)
    res <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens
    too.close = which(kriging.sd/sqrt(model@covariance@sd2) < 1e-06)
    res[too.close] <- max(0, m - kriging.mean)
    return(res)
}

generate_knots <- function (knots.number = NULL, d, lower = NULL, upper = NULL){
    if (is.null(lower))
        lower <- rep(0, times = d)
    if (is.null(upper))
        upper <- rep(1, times = d)
    if (is.null(knots.number))
        return(NULL)
    if (length(knots.number) == 1) {
        if (knots.number > 1) {
            knots.number <- rep(knots.number, times = d)
        } else {
            return(NULL)
        }
    }
    if (length(knots.number) != d) {
        print("Error in function generate_knots. The size of the vector knots.number needs to be equal to d")
        return(NULL)
    }
    knots.number <- pmax(1, knots.number)
    thelist <- NULL
    for (i in 1:d) {
        thelist[[i]] <- seq(from = lower[i], to = upper[i], length = knots.number[i])
    }
    return(thelist)
}

#' @test X=matrix(runif(10),ncol=1); y=-sin(2*pi*X); kmi <- km(design=X,response=y); TSEE(runif(100),kmi)
#' @test X=matrix(runif(10),ncol=1); y=-sin(2*pi*X); kmi <- km(design=X,response=y); DiceView::sectionview.fun(function(x)TSEE(x,kmi),dim=1)
#' @test X=matrix(runif(10),ncol=2); y=branin(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)TSEE(x,kmi),dim=2)
#' @require KrigInv
TSEE = function (x, model, T=0) {
  if (!is.matrix(x)) x <- matrix(x,ncol= model@d)
  y <- t(x)
    if ((nrow(y) == 1) && (ncol(y) == model@d)) {
        z <- y
    }
    else {
        if (ncol(x) == model@d) 
            z <- x
        if (ncol(x) != model@d) 
            z <- y
    }
    newdata <- x
    colnames(newdata) = colnames(model@X)
    krig <- predict.km(object = model, newdata = newdata, type = "UK", se.compute = TRUE)
    mk <- krig$mean
    sk <- krig$sd
    t <- (T - mk)/sk
    ski_dnorm_t <- sk * dnorm(t)
    C <- ((T - mk) * pnorm(t) + ski_dnorm_t) * ((mk - T) * pnorm(-t) + ski_dnorm_t)
    C[is.nan(C)] <- 0
    return(C)
}

#' @test X=matrix(runif(20),ncol=2); y=branin(X); kmi <- km(design=X,response=y); y.constr = rowSums(X)-1; kmi.constr=km(design=X,response=y.constr); DiceView::contourview.fun(function(x)EITSEE(x,kmi,kmi.constr),dim=2)
EITSEE = function(x,model, model.constr){
  return(EI(x,model)*TSEE(x,model.constr))
}

#' @test set.seed(1); X=matrix(runif(20),ncol=2); y=branin(X); kmi <- km(design=X,response=y); y.constr = rowSums(X)-1; kmi.constr=km(design=X,response=y.constr); DiceView::contourview.fun(function(x)EITSEE(x,kmi,kmi.constr),dim=2); points(max_EITSEE(kmi,kmi.constr,lower=c(0,0),upper=c(1,1))$par)
max_EITSEE <-function(model, model.constr, lower, upper, control=NULL) {
  
  d <- ncol(model@X)
  
  if (is.null(control$print.level)) control$print.level <- 1
  if (is.null(control$max.parinit.iter)) control$max.parinit.iter <- 10^d
  if(d<=6) N <- 10*2^d else N <- 100*d 
  if (is.null(control$pop.size))  control$pop.size <- N
  if (is.null(control$solution.tolerance))  control$solution.tolerance <- 1e-15  
  
  pars <- sobol_inside(n = N,dim=d,lower=lower,upper=upper,init=control$init)
  #t=Sys.time()
  eitsee <- EITSEE(pars,model, model.constr)
  #print(capture.output(Sys.time()-t))
  print(cbind(pars,eitsee))
  
  good_start = which(eitsee==max(eitsee,na.rm=T))
  par0=matrix(pars[good_start[sample(1:length(good_start),1)],],nrow=1)
  
  o <- psoptim(par=par0,fn=function(x){
    EITSEE(x,model, model.constr)
  },lower=lower,upper=upper,
  #control=list(vectorize=TRUE, fnscale=-1, trace=control$print.level, hybrid=FALSE, s=control$pop.size, abstol=control$solution.tolerance,maxit=10*d))
  control=list( fnscale=-1, trace=control$print.level,maxit=10*d))
  
  o$par <- t(as.matrix(o$par))
  colnames(o$par) <- colnames(model@X)
  o$value <- as.matrix(o$value)
  colnames(o$value) <- "EITSEE"
  
  return(list(par=o$par, value=o$value, counts=o$counts,par.all=o$par.all))
}

#' @test sobol_inside(10,2,0,1)
#' @test sobol_inside(10,1,0,1)
#' @require randtoolbox
sobol_inside <- function(n,dim,lower=0,upper=1,...) {
  s = sobol(n=n,dim=dim,...)
  #s * (upper-lower) - lower
  matrix(s * (upper-lower) - lower,ncol=dim)
}

#' @test set.seed(1); X=matrix(runif(20),ncol=2); y=branin(X); kmi <- km(design=X,response=y); y.constr = rowSums(X)-1; kmi.constr=km(design=X,response=y.constr); DiceView::contourview.fun(function(x)EITSEE(x,kmi,kmi.constr),dim=2); points(max_qEITSEE(kmi,kmi.constr,npoints=5,L="upper95",L.constr="upper95",lower=c(0,0),upper=c(1,1))$par)
max_qEITSEE <- function(model, model.constr, npoints, L, L.constr, lower, upper,  control=NULL, ...) {
  n1 <- nrow(model@X)
  for (s in 1:npoints) {
    oEGO <- max_EITSEE(model=model, model.constr=model.constr, lower=lower, upper=upper, control, ...)
    
    if (distXmin(oEGO$par,model@X)<=prod(upper-lower)*1E-10) {warning("Proposed a point already in design !");npoints=s-1;break;}
    
    model@X <- rbind(model@X, oEGO$par)
    
    if (L=="min")
      l = min(y)
    else if (L=="max")
      l = max(y)
    else if (L=="upper95") 
      l = predict.km(object = model,newdata = oEGO$par,type="UK",light.return = TRUE)$upper95
    else if (L=="lower95")
      l = predict.km(object = model,newdata = oEGO$par,type="UK",light.return = TRUE)$lower95
    else l = L
    
    model@y <- rbind(model@y, l, deparse.level=0)
    
    model@F <- trendMatrix.update(model, Xnew=data.frame(oEGO$par))
    if (model@noise.flag) { 
      model@noise.var = c(model@noise.var, 0) # here is the fix!
    }
    newmodel = NULL
    try(newmodel <- computeAuxVariables(model))
    if (is.null(newmodel)) {warning("Unable to update model !");npoints=s-1;break;}
    model = newmodel
    
    model.constr@X <- rbind(model.constr@X, oEGO$par)
    
    if (L.constr=="min")
      l.constr = min(y)
    else if (L.constr=="max")
      l.constr = max(y)
    else if (L.constr=="upper95") 
      l.constr = predict.km(object = model.constr,newdata = oEGO$par,type="UK",light.return = TRUE)$upper95
    else if (L.constr=="lower95")
      l.constr = predict.km(object = model.constr,newdata = oEGO$par,type="UK",light.return = TRUE)$lower95
    else l.constr = L.constr
    
    model.constr@y <- rbind(model.constr@y, l.constr, deparse.level=0)
    
    model.constr@F <- trendMatrix.update(model.constr, Xnew=data.frame(oEGO$par))
    if (model.constr@noise.flag) { 
      model.constr@noise.var = c(model.constr@noise.var, 0) # here is the fix!
    }
    newmodel.constr = NULL
    try(newmodel.constr <- computeAuxVariables(model.constr))
    if (is.null(newmodel.constr)) {warning("Unable to update model.constr !");npoints=s-1;break;}
    model.constr = newmodel.constr
    
  }
  #cat("  /max_qEI\n")
  return(list(par = model@X[(n1+1):(n1+npoints),, drop=FALSE], value = model@y[(n1+1):(n1+npoints),, drop=FALSE])) 
}  
