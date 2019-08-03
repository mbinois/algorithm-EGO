[![Build Status](https://travis-ci.org/Funz/algorithm-EGO.png)](https://travis-ci.org/Funz/algorithm-EGO)

# Funz algorithm: EGO

* Efficient Global Optimization (EGO)
* tags: optimization; sparse
* author: yann.richet@irsn.fr; DiceKriging authors
* require: DiceDesign; DiceKriging; DiceView; pso; jsonlite
* options: search_ymin='true'; initBatchSize='4'; batchSize='4'; iterations='10'; initBatchBounds='true'; trend='y~1'; covtype='matern3_2'; knots='0'; liar='upper95'; seed='1'
* options.help: search_ymin=minimization or maximisation; initBatchSize=Initial batch size; batchSize=iterations batch size; iterations=number of iterations; initBatchBounds=add input variables bounding values (2^d combinations); trend=(Universal) kriging trend; covtype=Kriging covariance kernel; knots=number of non-stationary points for each Xi; liar=liar value for in-batch loop (when batchsize>1); seed=random seed
* input: x=list(min=0,max=1)
* output: y=0.99


![Analytics](https://ga-beacon.appspot.com/UA-109580-20/algorithm-EGO)
