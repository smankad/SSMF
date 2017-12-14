library("tm")
data("crude")
x <- DocumentTermMatrix(crude,
	control = list(weighting =
        	function(x)
                	weightTfIdf(x, normalize = FALSE))
	)

n <- nrow(x) ## number of observations/documents
p <- ncol(x) ## number of terms

y <- sample(x=0:2, size=n, replace=T)
apps <- sample(x=c("App 1", "App 2"), size=n, replace=T) ## two apps 
time <- rep("time 1", n)  ## 1 time point


source("ssmf.R")
result <- ssnmf(y,x, numComp = 2)

result$intercept
result$beta

source("ssmf_glmnet.R")
result <- ssnmf_wglmnet(y,x, numComp = 2)

result$intercept
result$beta
result$Lambda

source("ssmf_continuation_ratio.R")
result <- ssmf_continuation_ratio(y,x, apps, time, numComp = 2)

result

