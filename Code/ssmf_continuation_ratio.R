ssmf_continuation_ratio <- function(Y, X, Apps, Time, numComp, maxiter = 50, type="single"){
  # single stage estimation based on continuation ratio model
  #
  # Args:
  #   Y: a vector of numeric ordinal ratings
  #   X: a document-term matrix from consumer reviews (the input data must be of class simple_triplet_matrix; used by default in the TM package)
  #   Apps: a vector of mobile apps for which the reviews are
  #   Time: a vector of time stamps
  #   numComp: number of topics
  #   maxiter: maximum number of iteration
  #   type: multiple (different coefficients for each rating category) or single (same coefficients for each rating category)
  #
  # Returns:
  #   Lambda: estimated Lambda matrix (nterm-by-numTopics)
  #   beta: estimated beta coefficients for each topic
  #   intercept: estimated intercepts for each topic
  #   ll: log likelihood
  #   delta: percentage decrease in log likelihood in the last iteration
  #   iter: number of iterations

recursive = function(func, x, y=NULL, z=NULL){
  # Recursively computes the function func taking arguments x, y and z.
  #
  # Args:
  #   func: the function applied to the atomic substructure of the following arguments
  #   x: a nested list of argument x0
  #   y: a nested list of argument y0, the same hierarchical structure as x
  #   z: a nested list of argument z0, the same hierarchical structure as x
  #
  # Returns:
  #   a nested list of func(x0,y0,z0) where x0,y0,z0 are atomic substructures of nested lists x,y,z
  if(is.null(y) && is.null(z)){
    if(is.atomic(x)|is.simple_triplet_matrix(x)){
      match.fun(func)(x)
    } else{
      Map(recursive, x, MoreArgs=list(fun=func))
    }
  } else if(!is.null(y) && is.null(z)){
    if((is.atomic(x)|is.simple_triplet_matrix(x)) && (is.atomic(y)|is.simple_triplet_matrix(y))){
      match.fun(func)(x,y)
    } else{
      Map(recursive, x,y, MoreArgs=list(fun=func))
    }
  }else{
    if((is.atomic(y)|is.simple_triplet_matrix(y)) && (is.simple_triplet_matrix(x)|is.atomic(x)) && (is.atomic(z)|is.simple_triplet_matrix(z))){
      #cat("B")
      match.fun(func)(x,y,z)
    } else{
      #cat("C")
      Map(recursive, x, y, z, MoreArgs=list(fun=func))
    }
  }
}

uneven_recursive = function(func, list1, list2, list3, list4=NULL, atom, apps, timepts){
  # Recursively computes the function func taking nested list arguments list1, list2, list3 and atomic arguments atom, apps, timepts.
  #
  # Args:
  #   func: the function applied to the atomic substructure of the following arguments
  #   list1: a nested list of argument l1
  #   list2: a nested list of argument l2, of the same hierarchical structure as list1
  #   list3: a nested list of argument l3, of the same hierarchical structure as list1
  #   list4: a nested list of argument l4, of the same hierarchical structure as list1
  #   atom: an argument that stays constant throughout recursion
  #   apps: a vector of apps --- one dimension along which recursion is applied
  #   timepts: a vector of time points --- another dimension along which recursion is applied
  # Returns:
  #   numeric sum of func(l1,l2,l3,l4) across all cells in nested lists of list1, list2, list3, list4
  if(is.null(list4)){
    y = 0
    for(a in 1:length(apps)){
      for(t in 1:length(timepts)){
        y = y + match.fun(func)(list1[[a]][[t]], list2[[a]][[t]], list3[[a]][[t]], atom)
        if(sum(is.na(y))>0){stop("uneven_recursive - a")}
      }
    }
  } else{
    y = 0
    for(a in 1:length(apps)){
      for(t in 1:length(timepts)){
        y = y + match.fun(func)(list1[[a]][[t]], list2[[a]][[t]], list3[[a]][[t]], list4[[a]][[t]], atom)
        if(sum(is.na(y))>0){stop("uneven_recursive - b")}
      }
    }
  }
  return(y)
}

  nlssubprob = function(Y_l, intercept, B, Linit, maxiter, X, Apps, Time, alpha = 0.1) {
    # a sub-function for ssmf_continuation_ratio ()
    #
    # Args:
    #   Y_l: a nested list of numeric consumer ratings
    #   X: a nested list of document-term matrices from consumer reviews
    #   Apps: a nested list of mobile apps for which the reviews are
    #   Time: a nested list of time stamps
    #   intercept: a nested list of estimated intercepts
    #   B: a nested list of estimated beta coefficients
    #   Linit: initial Lambda matrix
    #   maxiter: maximum number of iteration
    #
    # Returns:
    #   L: estimated Lambda matrix (nterm-by-numTopics)
    #   gradL: estimated gradient
    #   iterL: number of iterations
    #   alpha: learning rate

    gamma = 0.5
    apps = unique(Apps)
    time = sort(unique(Time))

    L = Linit
    iter = 1
    
    ll2 = uneven_recursive(func = log_likelihood, list1 = Y_l, list2 = X, list3 = B, list4 = intercept, atom = L, apps = apps, time = time)

    grad = uneven_recursive(func = deriv_lambda, list1 = Y_l, list2 = X, list3 = B, list4 = intercept, atom = L, apps = apps, time = time)
    if (all(grad == 0)) { return(list(L=L,gradL=grad,iterL=iter, alpha=alpha)) }
    while(iter < maxiter) {
      Ln = L + alpha*grad 
      Ln[Ln < 0] = 0
      d = Ln-L
      gradd = sum(grad*d)
      ll1 = uneven_recursive(func = log_likelihood, list1 = Y_l, list2 = X, list3 = B, list4 = intercept, atom = Ln, apps = apps, time = time)
      cond_bert <- ll1 - ll2
      suff_decr <- cond_bert <= 0.01 * gradd
      if (suff_decr) {
        L <- Ln
        break
        alpha = alpha/gamma
        ll2 = uneven_recursive(func = log_likelihood, list1 = Y_l, list2 = X, list3 = B, list4 = intercept, atom = L, apps = apps, time = time)
        grad = uneven_recursive(func = deriv_lambda, list1 = Y_l, list2 = X, list3 = B, list4 = intercept, atom = L, apps = apps, time = time)
      }else {
        alpha = alpha * gamma
      }
      iter <- iter + 1
    }
    return(list(L=L,gradL=grad,iterL=iter, alpha=alpha))
  }
  
  deriv_lambda <- function(y_l, x_l, beta0, intercept0, Lambda0) {
    # calculate the partial derivatives of Lambda
    #
    # Args:
    #   y_l: a vector of numeric consumer ratings
    #   x_l: a vector of document-term matrices from consumer reviews
    #   beta0: a vector of estimated beta coefficients
    #   Lambda0: initial Lambda matrix
    #
    # Returns:
    #   deriv: derivative estimates of Lambda

    #list1[[a]][[t]], list2[[a]][[t]], list3[[a]][[t]], atom
    deriv <- array(0, c(nrow(Lambda0), ncol(Lambda0)))
    noBeta <- which(is.na(unlist(beta0)))
    if (length(noBeta) == 0) {	  
	    XLB <- list()
	    sum1_y <- list()
	    for (m in 1:(length(y_l)-1)) {
	      LB <- Lambda0 %*% beta0[[m]]
	      XLB[[m]] <- intercept0[[m]] + as.array(x_l) %*% LB
	      
	      sum1_y[[m]] <- 0
	      for (j in 1:m) {
		sum1_y[[m]] <- sum1_y[[m]] + y_l[[j]]
	      }
	      
	      XLB[[m]][XLB[[m]] > 50] <- 50 ## happens about one every 10,000 elements due to integer overflow
	      c1 <- y_l[[m]] /(1+exp(XLB[[m]]))
	      c2 <- (1 - sum1_y[[m]]) * -1*exp(XLB[[m]])/(1+exp(XLB[[m]]))
	      c_final <- as.numeric((c1 + c2))
	      for (i in 1:length(y_l[[m]])) {

	        deriv <- deriv +  c_final[i] * (t(matrix(x_l[i,], nrow=1)) %*% t(beta0[[m]]))#xt_bt[[m]]
	      }
	    }
     }
     return(deriv)
  }

  log_likelihood <- function(Yi, Xi, beta, intercept, Lambda, neg=T) {
    # estimate the log likelihood
    #
    # Args:
    #   Yi: a vector of numeric consumer ratings
    #   Xi: a vector of document-term matrices from consumer reviews
    #   beta: a vector of estimated beta coefficients
    #   intercept: estimated intercept
    #   Lambda: estimated Lambda matrix
    #   neg: log likelihood being negative or not
    # Returns:
    #   ll: estimated log likelihood value

    ll <- 0
    XLB <- list()
    sum1_y <- list()
    noBeta <- which(is.na(unlist(beta)))
    if (length(noBeta) == 0) {
	    for (m in 1:(length(Yi)-1)) {     
		      LB <- Lambda %*% beta[[m]]
		      XLB[[m]] <- intercept[[m]] + as.array(Xi) %*% LB

		      sum1_y[[m]] <- 0
		      for (j in 1:m) {
      			sum1_y[[m]] <- sum1_y[[m]] + Yi[[j]] 
		      }
	    }
	    for (i in 1:length(Yi[[1]])) {
	      for (m in 1:(length(Yi)-1)) {
			XiLB <- as.numeric(XLB[[m]][i,])
			sum1_yi <- sum1_y[[m]][i]
			sum2_y <- sum1_yi
			if (XiLB > 25) {
			  c1 <- Yi[[m]][i] #sum1_yi
			} else {
			  c1 <- Yi[[m]][i] * (XiLB - log(1+exp(XiLB)))
			}
			ll <- ll+c1
			if (XiLB > 25) {
			    c2 <- (1 - sum2_y) * - XiLB
			} else {
			    c2 <- (1 - sum2_y) * -log(1+exp(XiLB))
			}
			ll <- ll+c2
	      }
	    }
    }
    if (neg) {return(as.numeric(-ll))}
    return(as.numeric(ll))
  }
  
  solveBeta = function(Yi, Lambda, Xi) {
    # estimate the beta coefficients (beta0 as intercept), heterogeneity assumed
    #
    # Args:
    #   Yi: a nested list of numeric consumer ratings
    #   Xi: a nested list of document-term matrices from consumer reviews
    #   Lambda: updated Lambda matrix
    # Returns:
    #   beta: a nested list of estimated beta coefficients
    #   intercept: a nested list of estimated intercept (beta0)

    noY <- which(is.na(unlist(Yi)))
    beta <- NA
    intercept <- NA
    if (length(noY) == 0) {
	    XiL <- as.array(Xi) %*% Lambda
	    colnames(XiL) <- paste("Topic",1:numComp)
	    beta <- list()
	    intercept <- list()
	    for (m in 1:(length(Yi)-1)) {
	      Ytmp <- rep(0, length(Yi[[m]]))
	      for (j in m:length(Yi)) {
	        Ytmp <- Ytmp + Yi[[j]]
	      }
	      Yi_logit <- Yi[[m]][Ytmp > 0]
	      XiLR <- XiL[Ytmp > 0,]
	      fit <- glm(Yi_logit  ~ XiLR, family="binomial")
	      
		beta[[m]] = coef(fit)[-1]
		intercept[[m]] <- coef(fit)[1]
	    }
    }
    return(list(beta=beta, intercept=intercept))
  }
  
  solveBeta_single = function(Yi, Lambda, Xi) {
    # estimate the beta coefficients (beta0 as intercept), homogeneity assumed
    #
    # Args:
    #   Yi: a nested list of numeric consumer ratings
    #   Xi: a nested list of document-term matrices from consumer reviews
    #   Lambda: updated Lambda matrix
    # Returns:
    #   beta: a nested list of estimated beta coefficients
    #   intercept: a nested list of estimated intercept (beta0)
    noY <- which(is.na(unlist(Yi)))
    beta <- NA
    intercept <- NA

    if (length(noY) == 0) {
	    XiL <- as.array(Xi) %*% Lambda
	    colnames(XiL) <- paste("Topic",1:numComp)  
	    
	    beta <- list()
	    intercept <- list()
	    
	    Yi_logit <- NULL
	    XiLFinal <- NULL
	    
	    for (m in 1:(length(Yi)-1)) {
	      Yi_tmp <- rep(0, length(Yi[[m]]))
	      for (j in m:length(Yi)) {
		      Yi_tmp <- Yi_tmp + Yi[[j]]
	      }
	      Yi_logit <- c(Yi_logit, Yi[[m]][Yi_tmp>0])
	      XiLFinal <- rbind(XiLFinal, cbind(rep(m, sum(Yi_tmp>0)) , XiL[Yi_tmp>0,]))
	    }
	    XiLFinal <- as.data.frame(XiLFinal)
	    XiLFinal[,1] <- as.factor(XiLFinal[,1])
	    XiLFinal <- model.matrix(~0+., XiLFinal)
	    
	    fit = glm(Yi_logit~ 0+ XiLFinal, family="binomial")
    
	    for (m in 1:(length(Yi)-1)) {
	      beta[[m]] <- coef(fit)[-1 * 1:(length(Yi)-1)]
	      intercept[[m]] <- coef(fit)[m]
	    }
    }    
    return(list(beta=beta, intercept=intercept))
  }
  
  p = ncol(X) # vocabulary size
  
  Lambda = array(abs(rnorm(p*numComp)), c(p,numComp))
  
  Y_l <- list()
  X_list = list()
  beta <- list()
  intercept = list()
  categories <- sort(unique(Y))
  apps  = (unique(Apps))
  timepts = sort(unique(Time))
  for(a in 1:length(apps)){
    
    Y_l[[a]] <- list()
    X_list[[a]] = list()
    intercept[[a]] = list()
    beta[[a]] = list()
    
    
    
    for(t in 1:length(timepts)){
      Y_l[[a]][[t]] = list()
#      if (sum(Apps==apps[a] & Time == timepts[t]) < 15) {
#	      X_list[[a]][[t]] = NA      
#	      Y_l[[a]][[t]] = NA      
#	      beta[[a]][[t]] = NA      
#	      intercept[[a]][[t]] = NA      
#      } else {
	      X_list[[a]][[t]] = X[Apps==apps[a] & Time == timepts[t],]
	      intercept[[a]][[t]] = list()
	      beta[[a]][[t]] = list()
	      
	      for (i in 1:length(categories)) {		
		Y_l[[a]][[t]][[i]] = sign(Y == categories[i])[Apps==apps[a] & Time == timepts[t]]
		beta[[a]][[t]][[i]] <- array(runif(numComp, min=-1, max=1))
		intercept[[a]][[t]][[i]] <- 0
		
	      }
#      }
    }
  }
  
  old <- uneven_recursive(func = log_likelihood, list1 = Y_l, list2 = X_list, list3 = beta, list4 = intercept, atom = Lambda, apps = apps, time = timepts)
  alpha <- 0.1
  for (i in 1:maxiter) {
    subsoln <- nlssubprob(Y_l, intercept, beta, Lambda, 20, X_list, Apps, Time, alpha)
    L = subsoln$L 

    if (type == "multiple") {
      fit = list()
      for(a in 1:length(apps)){
        fit[[a]] = list()
        for(t in 1:length(timepts)){
          fit[[a]][[t]] <- solveBeta(Y_l[[a]][[t]], L, X_list[[a]][[t]])
        }
      }     
    } else if (type == "single") {
      fit = list()
      for(a in 1:length(apps)){
        fit[[a]] = list()
        for(t in 1:length(timepts)){
          fit[[a]][[t]] <- solveBeta_single(Y_l[[a]][[t]], L, X_list[[a]][[t]]) 
        }
      }
    }
    B <- Map(function(a) Map(function(x) x$beta, a), fit)
    I <- Map(function(a) Map(function(x) x$intercept, a), fit)


    new = uneven_recursive(func = log_likelihood, list1 = Y_l, list2 = X_list, list3 = B, list4 = I, atom = L, apps = apps, time = timepts)
    delta = (old - new) / old

    if (new == 0) {
    	warning("Achieved exact fit to the data. You are likely overfitting!")
	beta <- B
	intercept <- I
	Lambda <- L
    	i <- i - 1
    	break #alpha <- min(alpha*0.01, subsoln$alpha/10)
    }

    if (delta < 0) {
    	warning("Monotonicity is broken. Algorithm ended early.")
    	i <- i - 1
    	break #alpha <- min(alpha*0.01, subsoln$alpha/10)
    }
    alpha <- max(1e-4,subsoln$alpha)
    old <- new
    beta <- B
    intercept <- I
    Lambda <- L
    if (delta < 1e-4 | new == 0) {break}
  }   
  rownames(Lambda) <- colnames(X)
  names(beta) <- apps
  names(intercept) <- apps
  for (i in 1:length(apps)) { 
	names(beta[[i]]) <- timepts
	names(intercept[[i]]) <- timepts
	for (j in 1:length(timepts)) {
		names(beta[[i]][[j]]) <- categories[-1 * length(categories)]
		names(intercept[[i]][[j]]) <- categories[-1 * length(categories)]
	}
  }
  return(list(Lambda = Lambda, beta = beta, intercept = intercept, ll = old, delta=delta, iters=i))
  
}
