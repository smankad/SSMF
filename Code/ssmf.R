require(Matrix)
require(slam)

ssnmf = function(Y,
                    X,
                    numComp = 1,
                    maxiter = 500,
                    tries = 1) {

  # single stage factorization
  #
  # Args:
  #   Y: a vector of numeric ordinal ratings
  #   X: a document-term matrix from consumer reviews (the input data must be of class simple_triplet_matrix; used by default in the TM package)
  #   numComp: number of topics
  #   maxiter: maximum number of iteration
  #   tries: number of random seeds to try
  #
  # Returns:
  #   Lambda: estimated Lambda matrix (nterm-by-numTopics)
  #   beta: estimated beta coefficients for each topic
  #   intercept: estimated intercepts for each topic
  #   normbest: the objective function value
  #   iter: number of iterations

  nlssubprob = function(Y,
                        intercept,
                        B,
                        Linit,
                        maxiter,
                        X, alpha) {

    # a sub-function for ssmf ()


    gamma = 0.9
    
    L = Linit
    iter = 1
    grad = deriv_lambda(Y, X, L, B, XtX)
    obj2 <- objectiveFn(Y = Y, X = X, Lambda = L, intercept = intercept, beta = B)
    
    while (iter < maxiter) {
      Ln = L - alpha * grad
      Ln[Ln < 0] = 0
      d = Ln - L
      gradd = sum(grad * d)
      obj1 <- objectiveFn(Y = Y, X = X, Lambda = Ln, intercept = intercept, beta = B)
      cond_bert <- obj1 - obj2
      suff_decr <- cond_bert <= 0.01 * gradd
      if (suff_decr) {
        L <- Ln
        break
        alpha = alpha / gamma
        grad = deriv_lambda(Y, X, L, B, XtX)
        obj2 <- objectiveFn(Y = Y, X = X, Lambda = L, intercept = intercept, beta = B)
      }
      else {
        alpha = alpha * gamma
      }
      iter <- iter + 1
    }
    return(list(
      L = L,
      gradL = grad,
      iterL = iter,
      alpha = alpha
    ))
  }
  
  deriv_lambda <- function(Y, X, Lambda, beta, XtX) { 
    # a sub-function for ssmf ()
    deriv <- array(0, c(nrow(Lambda), ncol(Lambda)))

    BBt = tcrossprod(beta, beta) 
    XY <- crossprod_simple_triplet_matrix(X, Y)
    XYBt = tcrossprod(XY, beta) 
    deriv <- deriv + XtX %*% Lambda %*% BBt - XYBt

    return(deriv)
  }
  
  objectiveFn <- function(Y, X, Lambda, beta, intercept) { 
    # a sub-function for ssmf ()
    obj <- sum(as.matrix(Y - intercept - matprod_simple_triplet_matrix(X, Lambda %*% beta)) ^ 2)
    return(obj)
  }
  
  solveBeta = function(Y, Lambda, X) {
    # a sub-function for ssmf ()
    XL <- as.array(X) %*% Lambda
    colnames(XL) <- paste("Topic",1:numComp)
    
    fit = lm(Y ~ XL)
    beta <- coef(fit)[-1]
    intercept <- coef(fit)[1]

    return(list(beta=beta, intercept=intercept))
  }
  
  ssnmf_Helper = function(Y,
                          Linit,
                          X,
                          maxiter,
                          XtX) {
    counter <- 1
    delta <- 999
    Lambda <- Linit
    intercept <- 0
    old <- 1e8
    alpha = 0.1
    while (counter < maxiter) {

      fit <- solveBeta(Y=Y, Lambda=Lambda, X = X)
      B <- fit$beta
      I <- fit$intercept
      
      subsoln = nlssubprob(Y, I, B, Lambda, 10, X, alpha)
      L = subsoln$L
      
      new <- objectiveFn(Y, X, L, B, I)
      
      delta = (old - new) / old
      # print(old)
      # print(new)
      #      write(old, "ssnmf_objectivevalues.txt", append=T)
      if (delta < 0) {
        warning("Monotonicity is broken")
        counter <- counter - 1
        break
      }
      alpha <- max(1e-4,subsoln$alpha)
      old <- new
      Lambda <- L
      beta <- B
      intercept <- I
      counter <-  counter + 1
      if (delta < 1e-4) {break}
    }

    return(
      list(
        Lambda = Lambda,
        beta = beta,
        intercept = intercept,
        counter = counter,
        objectiveFn = old
      )
    )
  }
  normbest <- Inf
  XtX <- crossprod_simple_triplet_matrix(X, X)
  for (j in 1:tries) {
    
    # Get starting value
    Lambda = array(abs(rnorm(ncol(X)*numComp)), c(ncol(X),numComp))
    
    tmp <-
      ssnmf_Helper(Y, Lambda, X, maxiter, XtX) # Perform a factorization
    # Save if this is the best so far
    if (tmp$objectiveFn < normbest) {
      Lambdabest <-  tmp$Lambda
      Betabest <-  tmp$beta
      intercept <-  tmp$intercept
      normbest <-  tmp$objectiveFn
      iters <-  tmp$counter
    }
  }
  if (normbest == Inf) {
    stop("ssnmf: Algorithm could not converge to a finite solution.")
  }
  return(
    list(
      Lambda = Lambdabest,
      beta = Betabest,
      intercept = intercept,
      iters = iters,
      objectiveValue = normbest
    )
  )
}
