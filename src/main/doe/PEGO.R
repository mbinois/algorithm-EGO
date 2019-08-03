#title: PEGO
#help: Profile Efficient Global Optimization (PEGO)
#tags: optimization; sparse
#author: Nicolas Garland; DiceKriging authors
#require: DiceDesign; DiceKriging; DiceView; pso; jsonlite
#options: search_ymin='true'; initBatchSize='4'; batchSize='4'; iterations='10'; initBatchBounds='true'; trend='y~1'; covtype='matern3_2'; knots='0'; liar='upper95'; seed='1';opt.index='1',profile.index='2'
#options.help: search_ymin=minimization or maximisation; initBatchSize=Initial batch size; batchSize=iterations batch size; iterations=number of iterations; initBatchBounds=add input variables bounding values (2^d combinations); trend=(Universal) kriging trend; covtype=Kriging covariance kernel; knots=number of non-stationary points for each Xi; liar=liar value for in-batch loop (when batchsize>1); seed=random seed
#input: list(x_opt=list(min=0,max=1),x_profile=list(min=0,max=1))
#output: y=0.99

#' constructor and initializer of R session
PEGO <- function(options) {

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

    ego$opt.index <- as.integer(unlist(strsplit(options$opt.index, ",")))
    ego$profile.index <- as.integer(unlist(strsplit(options$profile.index, ",")))

    ego$seed <- as.integer(options$seed)

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

getNextDesign <- function(pego,X,Y) {
    if (pego$i > pego$iterations) return();

    set.seed(pego$seed)

    d = dim(X)[2]
    if (dim(Y)[2] == 2) {
        algorithm$noise.var <- as.array(Y[,2])^2
    } else {
        algorithm$noise.var <- NULL
    }

    if (pego$search_min) {y=Y[,1]} else {y=-Y[,1]}

    pego$kmi <- km(control=list(trace=FALSE),pego$trend,optim.method='BFGS',covtype=pego$covtype, noise.var = pego$noise.var,design=X,response=y)

    PEGOi <- max_qPEI(model=pego$kmi,npoints=pego$batchSize,L=pego$liar,lower=rep(0,d),upper=rep(1,d),opt.index=pego$opt.index,control=list(trace=FALSE))

    if (is.null(PEGOi)) return()

    Xnext <- PEGOi$par
    pego$i <- pego$i + 1

    Xnext = as.matrix(Xnext)
    colnames(Xnext) <- names(pego$input)
    return(Xnext)
}

displayResults <- function(pego,X,Y) {
    pego$files <- paste("view_",pego$i-1,".o",length(pego$opt.index),".p",length(pego$profile.index),".png",sep="")
    resolution <- 600

    if (dim(Y)[2] == 2) {
        noise.var<- as.array(Y[,2])^2
        yname = paste0("N(", colnames(Y)[1], ",", colnames(Y)[2],")")
    } else {
        noise.var<- NULL
        yname = colnames(Y)
    }

    y_scale = ifelse(pego$search_min,1,-1)
    obj = ifelse(pego$search_min,"min","max")

    low_o = apply(X[,pego$opt.index,drop=F],2,min)
    up_o = apply(X[,pego$opt.index,drop=F],2,max)
    low_p = apply(X[,pego$profile.index,drop=F],2,min)
    up_p = apply(X[,pego$profile.index,drop=F],2,max)

    to01_p = function(xp){
        if (!is.matrix(xp)) xp = as.matrix(xp)
        (xp-matrix(low_p,ncol=ncol(xp),nrow=nrow(xp),byrow = T))/(matrix(up_p,ncol=ncol(xp),nrow=nrow(xp),byrow = T)-matrix(low_p,ncol=ncol(xp),nrow=nrow(xp),byrow = T))
    }

    to01_o = function(xo){
        if (!is.matrix(xo)) xo = as.matrix(xo)
        (xo-matrix(low_o,ncol=ncol(xo),nrow=nrow(xo),byrow = T))/(matrix(up_o,ncol=ncol(xo),nrow=nrow(xo),byrow = T)-matrix(low_o,ncol=ncol(xo),nrow=nrow(xo),byrow = T))
    }


    f = function(x) predict.km(pego$kmi,newdata=x,light.return = T,type = "UK")
    f_p = function(x_p,stat="mean") {
        f_o = function(x_o) {if (!is.matrix(x_o)) x_o=matrix(x_o,ncol=length(pego$opt.index),byrow = T);x=matrix(NaN,nrow(x_o),ncol(X));x[,pego$opt.index]=x_o;x[,pego$profile.index]=x_p;f(x)[[stat]]}

        ## n * BFGS
        #value=array(NA,10)
        #for (i in 1:length(value))
        #    value=min(value,optim(runif(length(pego$opt.index)),f_o,method="L-BFGS-B")$value,na.rm = T)

        ## PSO
        # o = psoptim(rep(0,length(pego$opt.index)),f_o,control = list(vectorize=T))

        ## rand100 > BFGS
        xo = matrix(runif(length(pego$opt.index)*10^length(pego$opt.index)),nrow=10^length(pego$opt.index))
        xo = rbind(xo,to01_o(X[,pego$opt.index]))
        m_p_xo = f_o(xo)
        value = optim(xo[which(m_p_xo==min(m_p_xo))[1],],f_o,method="L-BFGS-B",lower = rep(0,length(pego$opt.index)),upper=rep(1,length(pego$opt.index)))$value

        return(value*y_scale)
    }

    if (length(pego$profile.index)==1) {

        x_p=seq(f=low_p,t=up_p,l=100)

        png(file=pego$files,bg="transparent",height=resolution,width = resolution)
        y_m = y_scale*Vectorize(function(x) f_p(x,stat = "mean"))(to01_p(x_p))
        try(plot(x=x_p,y=y_m,xlab=names(X)[pego$profile.index],ylab=yname, main=paste0(obj,"(",yname,")"),type='lp'))
        try(polygon(x=c(x_p,rev(x_p)),y=c(y_scale*Vectorize(function(x) f_p(x,stat="lower95"))(to01_p(x_p)),rev(y_m)),col='gray',border = F))
        try(points(x=X[,pego$profile.index],y=y_scale*Y[,1],pch=20,col='red'))
        dev.off()

    } else if (length(pego$profile.index) == 2) {

        x_p = expand.grid(seq(f=low_p[1],to=up_p[1],l=20),seq(f=low_p[2],to=up_p[2],l=20))

        png(file=pego$files,bg="transparent",height=resolution,width = resolution)
        y_m = matrix(y_scale*apply(to01_p(x_p),1,function(x) f_p(x,stat = "mean")),ncol=20,nrow=20,byrow = T)
        y_l95 = matrix(y_scale*apply(to01_p(x_p),1,function(x) f_p(x,stat = "lower95")),ncol=20,nrow=20,byrow = T)

        l_m = contourLines(x=seq(f=low_p[1],to=up_p[1],l=20),y=seq(f=low_p[2],to=up_p[2],l=20),z=y_m)
        lev=unique(unlist(lapply(l_m,function(i)i$level)))

        try(contour(x=seq(f=low_p[1],to=up_p[1],l=20),y=seq(f=low_p[2],to=up_p[2],l=20),z=y_m,xlab=names(X)[pego$profile.index],ylab=yname, main=paste0(obj,"(",yname,")"),levels=lev))

        try(contour(x=seq(f=low_p[1],to=up_p[1],l=20),y=seq(f=low_p[2],to=up_p[2],l=20),z=y_l95,xlab=names(X)[pego$profile.index],ylab=yname, main=paste0(obj,"(",yname,")"),levels=lev,add=T,col='gray'))

        try(points(x=X[,pego$profile.index[1]],y=X[,pego$profile.index[2]],pch=20,col='red'))
        dev.off()

    } else {

        png(file=pego$files,bg="transparent",height=resolution,width = resolution)
        try(pairs(X,pch=20,col='red'))
        dev.off()

    }

    return(paste(sep="<br/>","<HTML><img src='",pego$files,"' width='",resolution,"' height='",resolution,"'/></HTML>"))
}

################### Algorithm dependencies ###################

distXmin <- function (x, Xmin){
    return(min(sqrt(rowSums((Xmin - matrix(x, nrow = nrow(Xmin), ncol = ncol(Xmin), byrow = TRUE))^2))))
}

#' @test X=matrix(runif(20),ncol=2); y=-sin(pi*rowSums(X)); kmi <- km(design=X,response=y); PEI(runif(100),kmi, 1, c(0,0), c(1,1))
#' @test X=matrix(runif(30),ncol=2); y=branin_mod(X); kmi <- km(design=X,response=y); DiceView::contourview.fun(function(x)PEI(x,kmi,1,c(0,0),c(1,1)),dim=2)
PEI <- function (x, model, opt.index, lower, upper, plugin=NULL) {
    ########################################################################################
    # Convert x in proper format(s)
    if (!is.matrix(x)) x <- matrix(x,ncol= model@d)
    d <- ncol(x)
    if (d != model@d){ stop("x does not have the right number of columns (",d," instead of ",model@d,")") }
    newdata <- x
    colnames(newdata) = colnames(model@X)

    # commentaire Nicolas : l'essentiel de la difference avec l'EI reside dans le plugin
    # au lieu de se comparer avec le min global, on va chercher a se comparer au min "pseudo-global" de la dimension a optimiser
    # attention : x peut contenir plusieurs points (data.frame ou matrix ?)
    # d'ou => plugin doit etre ici un vecteur de meme taille que le nombre de lignes de x
    if (is.null(plugin)){
        if(length(opt.index)<=6) N <- 10*2^length(opt.index) else N <- 100*length(opt.index)

        plugin <- apply(x, 1, FUN = function(x1) {
            # pour chaque x, on calcule le min de la moyenne de krigeage, en faisant varier uniquement les dimensions opt.index
            xopt=NULL
            for (i in opt.index) xopt=cbind(xopt,matrix(runif(N,lower[i],upper[i]),ncol=1))

            nd1 <- matrix(rep(x1,N), byrow=TRUE, nrow=N, ncol=d)
            nd1[,opt.index] <- xopt
            colnames(nd1) <- colnames(model@X)
            krig.mean <- predict.km(object=model, newdata=nd1, type="UK", checkNames = FALSE, se.compute = FALSE, light.return = TRUE)$mean

            good_start = which(krig.mean==min(krig.mean,na.rm=T))
            xopt0=matrix(xopt[good_start[sample(1:length(good_start),1)],],nrow=1)

            o <- psoptim(par=xopt0,fn=function(x.opt){
                x.tmp <- x1
                x.tmp[opt.index] <- x.opt
                if(!is.matrix(x.tmp)) x.tmp <- matrix(x.tmp,ncol= model@d)
                colnames(x.tmp) <- colnames(model@X)
                krig.mean <- predict.km(object=model, newdata=x.tmp, type="UK", checkNames = FALSE, se.compute = FALSE, light.return = TRUE)$mean
                return(krig.mean)
            },lower=lower[opt.index],upper=upper[opt.index],
            control=list(trace=0, maxit=10*length(opt.index)))

            return(o$value)
        })
        plugin <- pmax(min(model@y), plugin)
        if (model@noise.flag) plugin <- plugin-2*sqrt(model@noise.var)
    }
    m <- plugin

    ########################################################################################
    #cat("predict...")
    predx <- predict.km(object=model, newdata=newdata, type="UK", checkNames = FALSE)
    #cat(" done.\n")
    kriging.mean <- predx$mean
    kriging.sd   <- predx$sd

    xcr <- (m - kriging.mean)/kriging.sd

    xcr.prob <- pnorm(xcr)
    xcr.dens <- dnorm(xcr)
    res <- (m - kriging.mean) * xcr.prob + kriging.sd * xcr.dens

    too.close = which(kriging.sd/sqrt(model@covariance@sd2) < 1e-06)
    res[too.close] <- max(0,m - kriging.mean)

    return(res)
}

#' @test set.seed(1); X=matrix(runif(20),ncol=2); y=branin(X); kmi <- km(design=X,response=y); kmi=km(design=X,response=y); DiceView::contourview.fun(function(x)EI(x,kmi),dim=2); points(max_EI(kmi,lower=c(0,0),upper=c(1,1))$par)
max_PEI <-function(model,  lower, upper, opt.index, control=NULL) {

    d <- ncol(model@X)

    if (is.null(control$print.level)) control$print.level <- 0
    if (is.null(control$max.parinit.iter)) control$max.parinit.iter <- 10^d
    if(d<=6) N <- 10*2^d else N <- 100*d
    if (is.null(control$pop.size))  control$pop.size <- N
    if (is.null(control$solution.tolerance))  control$solution.tolerance <- 1e-15

    pars=NULL
    for (i in 1:d) pars=cbind(pars,matrix(runif(N,lower[i],upper[i]),ncol=1))

    #t=Sys.time()
    ei <- PEI(pars,model,opt.index=opt.index,lower=lower,upper=upper)
    #print(capture.output(Sys.time()-t))
    print(cbind(pars,ei))

    good_start = which(ei==max(ei,na.rm=T))
    par0=matrix(pars[good_start[sample(1:length(good_start),1)],],nrow=1)

    o <- psoptim(par=par0,fn=function(x){
        PEI(x,model,opt.index=opt.index,lower=lower,upper=upper)
    },lower=lower,upper=upper,
    control=list( fnscale=-1, trace=control$print.level,maxit=10*d))

    o$par <- t(as.matrix(o$par))
    colnames(o$par) <- colnames(model@X)
    o$value <- as.matrix(o$value)
    colnames(o$value) <- "EI"

    return(list(par=o$par, value=o$value, counts=o$counts,par.all=o$par.all))
}

#' @test set.seed(1); X=matrix(runif(20),ncol=2); y=apply(FUN=branin,X,1); kmi <- km(design=X,response=y);  kmi=km(design=X,response=y); DiceView::contourview.fun(function(x)EI(x,kmi),dim=2); points(max_qEI(kmi,npoints=5,L="upper95",lower=c(0,0),upper=c(1,1))$par)
max_qPEI <- function(model, npoints, L,  lower, upper, opt.index, control=NULL, ...) {
    n1 <- nrow(model@X)
    for (s in 1:npoints) {
        oEGO <- max_PEI(model=model, lower=lower, upper=upper, opt.index=opt.index, control, ...)

        if (distXmin(oEGO$par,model@X)<=prod(upper-lower)*1E-10) {
            warning("Proposed a point already in design !");
            npoints=s-1;
            break;
        }

        model@X <- rbind(model@X, oEGO$par)

        if (L=="min")
            l = min(model@y)
        else if (L=="max")
            l = max(model@y)
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
    }

    if (npoints==0) return()
    return(list(par = model@X[(n1+1):(n1+npoints),, drop=FALSE], value = model@y[(n1+1):(n1+npoints),, drop=FALSE]))
}
