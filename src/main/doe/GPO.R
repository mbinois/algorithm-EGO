#title: GPO
#help: GPareto wrapper for multi-objective Efficient Global Optimization (EGO)
#tags: multi-objective optimization; sparse
#author: yann.richet@irsn.fr; mickael.binois@inria.fr; DiceKriging authors; GPareto authors
#require: DiceDesign; DiceKriging; DiceView; pso; jsonlite; foreach; doParallel; GPareto; xtable
#options: initBatchSize='4'; iterations='10'; initBatchBounds='true'; trend='y~1'; covtype='matern3_2'; knots='0'; nugget='true'; seed='1'; refPoint="NULL"; noise.var="NULL"; ncb='1000'
#options.help: initBatchSize=Initial batch size; batchSize=iterations batch size; iterations=number of iterations; initBatchBounds=add input variables bounding values (2^d combinations); trend=(Universal) kriging trend; covtype=Kriging covariance kernel; knots=number of non-stationary points for each Xi; seed=random seed; refPoint= Vector defining hypervolume reference; noise.var=Vector of known constant noise variance of each objective; ncb=integer number of candidate batches to generate
#input: x=list(min=0,max=1)
#output: y=0.99

GPO <- function(options) {
  
  library(DiceDesign)
  library(DiceKriging)
  library(GPareto)
  library(DiceView)
  library(pso)
  library(jsonlite)
  library(doParallel)
  if (!foreach::getDoParRegistered()) foreach::registerDoSEQ()
  
  gpo = new.env()
  gpo$i <- 0
  
  gpo$initBatchSize <- as.integer(options$initBatchSize)
  gpo$batchSize <- as.integer(options$batchSize)
  gpo$iterations <- as.integer(options$iterations)
  gpo$initBatchBounds <- as.logical(options$initBatchBounds)
  
  gpo$trend <- as.formula(options$trend)
  gpo$covtype <- as.character(options$covtype)
  gpo$liar <- as.character(options$liar)
  gpo$knots <- as.integer(unlist(strsplit(as.character(options$knots),",")))
  gpo$nuggetestim <- isTRUE(as.logical(options$nugget))
  if (gpo$nuggetestim) {
    gpo$nugget <- NULL
  } else {
    gpo$nugget <- as.numeric(options$nugget)
    if (!is.numeric(gpo$nugget) | is.na(gpo$nugget)) gpo$nugget <- NULL
  }
  gpo$noise.var <- options$noise.var
  if(!is.null(gpo$noise.var)) gpo$nugget <- NULL
  gpo$seed <- as.integer(options$seed)
  gpo$refPoint <- options$refPoint
  gpo$ncb <- 1e3
  
  return(gpo)
}

getInitialDesign <- function(algorithm, input, output) {
  algorithm$input <- input
  algorithm$output <- output
  
  d = length(input)
  if (!is.numeric(algorithm$initBatchSize))
    algorithm$initBatchSize = floor(eval(parse(text = algorithm$initBatchSize)))
  
  set.seed(algorithm$seed)
  d = length(input)
  lhs <- maximinESE_LHS(lhsDesign(n = algorithm$initBatchSize, dimension = d,seed=algorithm$seed)$design)$design
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
  nobj <- dim(Y)[2]
  # if (dim(Y)[2] == 2) {
  #   algorithm$noise.var <- as.array(Y[,2])^2
  #   algorithm$nuggetestim = FALSE
  #   algorithm$nugget = NULL
  # } else {
  # algorithm$noise.var <- NULL
  # }
  
  # y = Y[,1,drop=FALSE]
  
  # heuristic for lower bound of theta : max(1e-6, 0.1 * dX[which.max(dy/rowSums(dX))])
  dX = apply(FUN = dist, X, MARGIN = 2)
  
  # dy = apply(FUN = dist, y, MARGIN = 2)
  
  # define stationary-changing points
  all_knots <- generate_knots(knots.number = algorithm$knots, d = d, lower = sapply(algorithm$input, "[[", "min"), upper = sapply(algorithm$input, "[[", "max"))
  models <- NULL
  
  for(i in 1:nobj){
    
    dy = apply(FUN = dist, Y[,i,drop = FALSE], MARGIN = 2)
    if (isTRUE(algorithm$knots>1))
      try(model <- km(algorithm$trend, optim.method = "BFGS",
                      covtype = algorithm$covtype,
                      design = X, response = Y[,i], noise.var = rep(algorithm$noise.var[i],nrow(X)),
                      lower = rep(pmax(1e-06, dX[which.max(dy/rowSums(dX)),]),each=algorithm$knots),
                      control = list(trace = FALSE),
                      nugget.estim = algorithm$nuggetestim,nugget = algorithm$nugget,
                      scaling = TRUE, knots = all_knots))
    else
      try(model <- km(algorithm$trend, optim.method = "BFGS",
                      covtype = algorithm$covtype,
                      design = X, response = Y[,i], noise.var = rep(algorithm$noise.var[i],nrow(X)),
                      lower = pmax(1e-06, dX[which.max(dy/rowSums(dX)),]),
                      control = list(trace = FALSE),
                      nugget.estim = algorithm$nuggetestim,nugget = algorithm$nugget,
                      scaling = FALSE, knots = NULL))
    if (is.null(model)) stop(paste("Error for construction of model ", i)) else models <- c(models, list(model))
  }
  
  observations <- Reduce(cbind, lapply(models, slot, "y"))
  if(is.null(algorithm$noise.var)){
    paretoFront <- unique(t(nondominated_points(t(observations))))
  }else{
    preds <- predict_kms(models, newdata=models[[1]]@X, type="UK", checkNames = FALSE, light.return = TRUE, cov.compute = FALSE)
    observations.denoised <- t(preds$mean)
    paretoFront <- unique(t(nondominated_points(t(observations.denoised))))
  }
  
  if(is.null(algorithm$refPoint)) refPoint <- NULL else refPoint <- algorithm$refPoint
  
  lower = sapply(algorithm$input, "[[", "min")
  upper = sapply(algorithm$input, "[[", "max")
  
  oGPO <- NULL
  if(algorithm$batchSize == 1){
    try(oGPO <- crit_optimizer(crit = "EHI", model = models, lower = lower, upper = upper,
                               optimcontrol = list(trace = 0), paretoFront = paretoFront,
                               critcontrol = list(refPoint = refPoint)))
  }else{
    ncb <- algorithm$ncb
    q <- algorithm$batchSize
    
    Xunif <- matrix(runif(ncb*d), ncb) %*% diag(upper - lower) + matrix(lower,ncb, d, byrow = TRUE)
    EHI_grid <- crit_EHI(x = Xunif, model = models, critcontrol = list(refPoint = refPoint))
    
    Xbcands <- array(NA, dim = c(ncb, q, d))
    for(i in 1:ncb) Xbcands[i,,] <- Xunif[sample(1:ncb, q, prob = pmax(0, EHI_grid)),]
    
    qEHI_grid <- apply(Xbcands, 1, crit_qEHI, model = models,
       critcontrol = list(refPoint = refPoint))
      
    Xq <- Xbcands[which.max(qEHI_grid),,]
    
    # More robust version
    crit_qEHI2 <- function(x, ...){
      tmp <- try(GPareto::crit_qEHI(matrix(x, q), ...), silent = TRUE)
      if(class(tmp) == "try-error") return(-1)
      return(tmp)
    } 
    
    oGPO <- optim(as.vector(Xq), crit_qEHI2,
                 model = models, lower = lower, upper = upper, method = "L-BFGS-B",
                 control = list(fnscale = -1), critcontrol = list(refPoint = refPoint, nb.samp = 100))
    oGPO$par <- matrix(oGPO$par, q)
  }

  
  if (is.null(oGPO))
    return()
  
  Xnext <- oGPO$par
  algorithm$i <- algorithm$i + 1
  
  Xnext = as.matrix(Xnext)
  colnames(Xnext) <- names(algorithm$input)
  return(Xnext)
}

displayResults <- function(algorithm, X, Y) {
  algorithm$files <- paste("GPO_view_", algorithm$i,".png", sep = "")
  resolution <- 600
  
  # if (dim(Y)[2] == 2) {
  #   noise.var <- as.array(Y[, 2])^2
  #   yname = paste0("N(", colnames(Y)[1], ",", colnames(Y)[2],")")
  # } else {
  # noise.var <- NULL
  # yname = colnames(Y)
  # }
  
  m <- t(nondominated_points(t(Y)))
  colnames(m) <- colnames(Y)
  d2 = dist2(Y,m)
  xi = apply(d2,2,function(x)which(x==0))
  x = X[xi,]
  colnames(x) <- colnames(X)
  
  html = paste0("<HTML>minimum is ", print(xtable(m,digits=6),type="html"),
                " found at <br/>",
                print(xtable(x,digits=3),type="html"),
                "<br/><img src='", algorithm$files,
                "' width='", resolution, "' height='", resolution, "'/></HTML>")
  html = paste0(html,"<min>",toJSON(m),"</min><argmin>",toJSON(x),"</argmin>")
  
  png(file = algorithm$files, bg = "transparent", height = resolution, width = resolution)
  reds=rep(0,nrow(X))
  reds[xi]=1
  try(pairs(cbind(X,Y),col=rgb(reds,0,1-reds)))
  dev.off()
  
  #if (algorithm$i == algorithm$iterations) {
  # html = paste0(html,"<data_json>",toJSON(as.data.frame(cbind(X,Y)),dataframe = "columns"),"</data_json>")
  # 
  # lower = sapply(algorithm$input, "[[", "min")
  # upper = sapply(algorithm$input, "[[", "max")
  # n = 1000
  # set.seed(123) # to get the same points for evaluating model
  # Xm = matrix(lower,nrow=n,ncol=length(lower),byrow = T) + matrix(upper-lower,nrow=n,ncol=length(lower),byrow = T) * matrix(runif(n*length(lower)),nrow=n,ncol=length(lower))
  # colnames(Xm) <- colnames(X)
  # Ym = list(mean=rep(NA,n),sd=rep(NA,n))
  # try(Ym <- predict(algorithm$model,newdata = Xm,type = "UK",cov.compute = F, low.memory = T,checkNames=F))
  # Ym = cbind(Ym$mean,Ym$sd)
  # colnames(Ym) <- c(colnames(Y),paste0("sd_",colnames(Y)))[1:2]
  # 
  # html = paste0(html,"<model_json>",toJSON(as.data.frame(cbind(Xm,Ym)),dataframe = "columns"),"</model_json>")
  # 
  # html = paste0(html,"<kriging_json>",toJSON(algorithm$model,force=TRUE,auto_unbox=TRUE,pretty=TRUE,dataframe = "columns"),"</kriging_json>")
  #}
  
  return(paste0(html,collapse=';'))
}

displayResultsTmp <- displayResults

################### Algorithm dependencies ###################

generate_knots <- function(knots.number=NULL,d,lower=NULL,upper=NULL){
  
  if(is.null(lower)) lower <- rep(0,times=d)
  if(is.null(upper)) upper <- rep(1,times=d)
  
  if(is.null(knots.number)) return(NULL)
  
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
  
  knots.number <- pmax(1,knots.number) # 2 knots at least per dimension
  
  thelist <- list()
  for (i in 1:d) {
    thelist[[i]] <- seq(from = lower[i], to = upper[i], length = knots.number[i])
  }
  return(thelist)
}

dist2 = function(X,Y) 
  outer(1:nrow(X),1:nrow(Y),Vectorize(function(ix,iy) sum((X[ix,]-Y[iy,])^2)))

########################################## EXAMPLE ##############################################
# f <- function(X) t(apply(X,1, function (x){
#   if(is.null(dim(x))){
#     x <- matrix(x, nrow = 1)
#   }
#   b1<-15*x[,1]-5
#   b2<-15*x[,2]
#   return(cbind((b2-5.1*(b1/(2*pi))^2+5/pi*b1-6)^2 +10*((1-1/(8*pi))*cos(b1)+1),
#                -sqrt((10.5-b1)*(b1+5.5)*(b2+0.5)) - 1/30*(b2 -5.1*(b1/(2*pi))^2-6)^2 - 1/3*((1-1/(8*pi))*cos(b1)+1)
#   ))
# }))
# 
# Xgrid <- as.matrix(expand.grid(seq(0,1,,101), seq(0,1,,101)))
# Ygrid <- f(Xgrid)
# PFref <- t(nondominated_points(t(Ygrid)))

## Default version: 2objs, noiseless, one design per turn 

# options = list(initBatchSize='8', batchSize='1', iterations='10', initBatchBounds='true', trend='y~1', covtype='matern3_2', knots='2', liar='upper95', nugget='true', seed='1')
# algorithm = GPO(options)
# 
# X0 = getInitialDesign(algorithm, input=list(x1=list(min=0,max=1),x2=list(min=0,max=1)), NULL)
# Y0 = f(X0)
# # X0 = getInitialDesign(gd, input=list(x2=list(min=0,max=1)), NULL)
# # Y0 = f1(X0)
# Xi = X0
# Yi = Y0
# 
# finished = FALSE
# while (!finished) {
#   # print(displayResultsTmp(algorithm,Xi,Yi))
#   Xj = getNextDesign(algorithm,Xi,Yi)
#   if (is.null(Xj) | length(Xj) == 0) {
#     finished = TRUE
#   } else {
#     Yj = f(Xj)
#     Xi = rbind(Xi,Xj)
#     Yi = rbind(Yi,Yj)
#   }
# }
# 
# print(displayResults(algorithm,Xi,Yi))
# 
# cols <- rep(4, nrow(Xi))
# cols[1:options$initBatchSize] <- 1
# plot(Yi, pch = 20, col = cols, xlab = "f1", ylab = "f2")
# PF <- t(nondominated_points(t(Yi)))
# plotParetoEmp(PF, col = "red")

## Noisy version: 2objs, small noise, one design per turn 

# options = list(initBatchSize='8', batchSize='1', iterations='10', initBatchBounds='true', trend='y~1', covtype='matern3_2', knots='2', liar='upper95', nugget='false', seed='1', noise.var = c(1,1))
# algorithm = GPO(options)
# 
# X0 = getInitialDesign(algorithm, input=list(x1=list(min=0,max=1),x2=list(min=0,max=1)), NULL)
# Y0 = f(X0)
# # X0 = getInitialDesign(gd, input=list(x2=list(min=0,max=1)), NULL)
# # Y0 = f1(X0)
# Xi = X0
# Yi = Y0 + matrix(rnorm(nrow(Xi) * 2), ncol = 2)
# 
# finished = FALSE
# while (!finished) {
#   # print(displayResultsTmp(algorithm,Xi,Yi))
#   Xj = getNextDesign(algorithm,Xi,Yi)
#   if (is.null(Xj) | length(Xj) == 0) {
#     finished = TRUE
#   } else {
#     Yj = f(Xj) + rnorm(2)
#     Xi = rbind(Xi,Xj)
#     Yi = rbind(Yi,Yj)
#   }
# }
# 
# print(displayResults(algorithm,Xi,Yi))
# 
# cols <- rep(4, nrow(Xi))
# cols[1:options$initBatchSize] <- 1
# plot(Yi, pch = 20, col = cols, xlab = "f1", ylab = "f2")
# PF <- t(nondominated_points(t(Yi)))
# plotParetoEmp(PF, col = "red")

## Batch version: 2objs, noiseless, fours designs per turn 

# options = list(initBatchSize='8', batchSize='4', iterations='6', initBatchBounds='true', trend='y~1', covtype='matern3_2', knots='2', liar='upper95', nugget='true', seed='1', refPoint = c(300,0))
# algorithm = GPO(options)
# 
# X0 = getInitialDesign(algorithm, input=list(x1=list(min=0,max=1),x2=list(min=0,max=1)), NULL)
# Y0 = f(X0)
# # X0 = getInitialDesign(gd, input=list(x2=list(min=0,max=1)), NULL)
# # Y0 = f1(X0)
# Xi = X0
# Yi = Y0
# PF0 <- t(nondominated_points(t(Yi)))
# 
# finished = FALSE
# while (!finished) {
#   # print(displayResultsTmp(algorithm,Xi,Yi))
#   Xj = getNextDesign(algorithm,Xi,Yi)
#   if (is.null(Xj) | length(Xj) == 0) {
#     finished = TRUE
#   } else {
#     Yj = f(Xj)
#     Xi = rbind(Xi,Xj)
#     Yi = rbind(Yi,Yj)
#   }
# }
# 
# print(displayResults(algorithm,Xi,Yi))
# 
# cols <- rep(4, nrow(Xi))
# cols[1:options$initBatchSize] <- 1
# plot(Yi, pch = 20, col = cols, xlab = "f1", ylab = "f2")
# PF <- t(nondominated_points(t(Yi)))
# plotParetoEmp(PF, col = "red")
# plotParetoEmp(PFref, lty = 3)
# plotParetoEmp(PF0, lty = 3, col = "grey")

# ## Test with 3objs  
# f <- function(X) t(apply(X,1, function (x){
#   nobj <- 3
#   if(is.null(dim(x))){
#     x <- matrix(x, 1) 
#   }
#   n <- ncol(x)
#   
#   y <- matrix(x[,1:(nobj-1)], nrow(x))
#   z <- matrix(x[,nobj:n], nrow(x))
#   
#   g <- rowSums((z-0.5)^2)
#   
#   #   tmp <- c(rev(cumprod(cos(y * pi/2))), 1)
#   #   tmp2 <- c(1, rev(sin(y * pi/2)))
#   tmp <- t(apply(cos(y * pi/2), 1, cumprod))
#   tmp <- cbind(t(apply(tmp, 1, rev)), 1)
#   
#   tmp2 <- cbind(1, t(apply(sin(y * pi/2), 1, rev)))
#   
#   f <- tmp * tmp2 * (1 + g)
#   
# }))
# ngrid <- 5e3
# Xgrid <- matrix(runif(4 * ngrid), ngrid)
# Ygrid <- f(Xgrid)
# PFref <- t(nondominated_points(t(Ygrid)))
# 
# ## Default version: 2objs, noiseless, one design per turn 
# 
# options = list(initBatchSize='40', batchSize='4', iterations='5', initBatchBounds='true', trend='y~1', covtype='matern3_2', knots='2', liar='upper95', nugget='true', seed='1', refPoint=c(2,2,2))
# algorithm = GPO(options)
# 
# X0 = getInitialDesign(algorithm, input=list(x1=list(min=0,max=1),x2=list(min=0,max=1),x3=list(min=0,max=1),x4=list(min=0,max=1)), NULL)
# Y0 = f(X0)
# 
# Xi = X0
# Yi = Y0
# 
# finished = FALSE
# while (!finished) {
#   # print(displayResultsTmp(algorithm,Xi,Yi))
#   Xj = getNextDesign(algorithm,Xi,Yi)
#   if (is.null(Xj) | length(Xj) == 0) {
#     finished = TRUE
#   } else {
#     Yj = f(Xj)
#     Xi = rbind(Xi,Xj)
#     Yi = rbind(Yi,Yj)
#   }
# }
# 
# idsPF <- !is_dominated(t(Yi))
# cols <- rep(3, length(idsPF))
# cols[!idsPF] <- 2
# library(rgl)
# plot3d(Yi, col = cols, size = 5)
# points3d(PFref)
# 
# pairs(Xi, col = cols, pch = 20)

