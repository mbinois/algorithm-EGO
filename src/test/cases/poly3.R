## This file should provide following objects, when loaded:
# f : function
# input.f : list of input dimensions, contains list of properties like lower & upper bounds of each dimensions
# output.f : list of output dimensions
# *.f : list of math properties. To be compared with algorithm results
# [print.f] : method to print/plot the function for information

f = function(X) {
    matrix(Vectorize(function(x) {((x-.75)/3)^3})(X),ncol=1)
}
input.f = list(
    x=list(min=0,max=1)
)
output.f = "poly3"
info.f = 0.75

test = function(algorithm_file) {
    results = run.algorithm(algorithm_file, options=NULL,fun=list(input=input.f,output=output.f))
    library(testthat)
    test_that("poly3 info",{expect_equal(as.numeric(results$info),info.f,tolerance = .0001)})
}

