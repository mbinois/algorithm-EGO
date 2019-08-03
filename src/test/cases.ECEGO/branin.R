## This file should provide following objects, when loaded:
# f : function
# input.f : list of input dimensions, contains list of properties like lower & upper bounds of each dimensions
# output.f : list of output dimensions
# *.f : list of math properties. To be compared with algorithm results
# [print.f] : method to print/plot the function for information

f <- function(x) {
	x1 <- x[,1]*15-5   
	x2 <- x[,2]*15     
	y = matrix((x2 - 5/(4*pi^2)*(x1^2) + 5/pi*x1 - 6)^2 + 10*(1 - 1/(8*pi))*cos(x1) + 10,ncol=1)
    y_constr = matrix((x[,1] - 0.5427730) * (x[,2] - 0.15), ncol=1)
    cbind(y,y_constr)
}
input.f = list(
    x1=list(min=0,max=1),
    x2=list(min=0,max=1)
)
output.f = c("branin","constr")
#argmin1.f = c(0.9616520, 0.15)
#argmin2.f = c(0.1238946, 0.8166644)
#argmin3.f = c(0.5427730, 0.15)
argmin.f = c(0.5427730, 0.15)
min.f = 0.3978874

library(testthat)
test_that("f(armgin.f) == f.min",{expect_equal(f(matrix(argmin.f,nrow=1))[1,1],min.f,tolerance = .0001)})
test_that("f_constr(armgin.f) == 0",{expect_equal(f(matrix(argmin.f,nrow=1))[1,2],0,tolerance = .0001)})

test = function(algorithm_file) {
    results = run.algorithm(algorithm_file, options=NULL,fun=list(input=input.f,output=output.f))
    test_that("branin constr min",{expect_equal(as.numeric(results$min),min.f,tolerance = .1)})
    test_that("branin constr argmin",{expect_equal(sum((jsonlite::fromJSON(results$argmin)-argmin.f)^2),0,tolerance = .01)})
}

