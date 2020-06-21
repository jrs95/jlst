#' jlssc
#'
#' jlssc performs the joint location-and-scale score test.
#' @param y vector of outcome values.
#' @param x vector of exposure values.
#' @param covar a data.frame of covariates.
#' @param type type of test (default: 1 [Breusch-Pagan variance test]; options: 1 [Breusch-Pagan variance test], 2 [Brown-Forsythe variance test], 3 [Method of moments version of test 1], 4 [Method of moments version of test 2]).
#' @return a data.frame of results. Q is the test statistic, DF is the degrees of freedom and P is the p-value.
#' @examples
#' x <- rbinom(1000, 1, 0.5)
#' y <- 0.5 + 0.025*x + rnorm(1000, 0, sqrt(0.025^2*x)) + rnorm(1000, 0, 0.1)
#' jlssc(y, x)
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
jlssc <- function(y, x, covar=NULL, type=1){
  
  # Errors
  if(!(is.numeric(y) | is.integer(y))) stop("y has to be a numeric variable")
  if(!(is.numeric(x) | is.integer(x) | is.factor(x))) stop("x has to be either a numeric variable or a factor variable")
  if(!is.null(covar)){if(!is.data.frame(covar)) stop("covar has to be a data.frame")}
  if(length(y)!=length(x)) stop("y is not the same size as x")
  if(!is.null(covar)){if(length(y)!=nrow(covar)) stop("y is not the same size as covar")}
  if(!(type %in% 1:4)) stop("type has to be set to either 1, 2, 3 or 4")
  
  # Missing values
  data <- cbind(y, x); if(!is.null(covar)){data <- cbind(data, covar)}
  keep <- complete.cases(data)
  y <- y[keep]; x <- x[keep]; if(!is.null(covar)){covar <- covar[keep,,drop=F]}
  
  # Covariates
  if(!is.null(covar)){
    covar <- model.matrix(as.formula(~ .), data=covar)[,-1,drop=F]
    if(any(is.na(covar))) stop("there are missing values in the covariates")
  }
  
  # Exposure
  if(is.factor(x)){
    x <- model.matrix(as.formula(~ as.factor(x)))[,-1,drop=F]
    if(any(is.na(x))) stop("there are missing values in the exposure")              
  }
  
  # Regress out covariates
  if(!is.null(covar)){
    y <- lm(y~covar)$resid
    x <- lm(x~covar)$resid
  }
  
  # Location + scale test
  if(type==1){
    n <- length(y)
    yt <- y - mean(y)
    dt <- yt^2 - sum(yt^2)/n
    xt <- scale(x, scale=F)
    beta <- lm(yt~xt-1)$coef
    delta <- lm(dt~xt-1)$coef
    theta <- c(beta, delta)
    sig <- matrix(c(mean(yt^2), mean(yt*dt), mean(yt*dt), mean(dt^2)), nrow=2, ncol=2)
    Q <- ((t(theta)%*%kronecker(solve(sig), crossprod(xt)))%*%theta)
    DF <- length(theta)
    P <- pchisq(Q, df=DF, lower.tail=F)
  }
  if(type==2){
    n <- length(y)
    yt <- y - mean(y)
    ym <- abs(y - median(y))
    dt <- ym - mean(ym)
    xt <- scale(x, scale=F)
    beta <- lm(yt~xt-1)$coef
    delta <- lm(dt~xt-1)$coef
    theta <- c(beta, delta)
    sig <- matrix(c(mean(yt^2), mean(yt*dt), mean(yt*dt), mean(dt^2)), nrow=2, ncol=2)
    Q <- ((t(theta)%*%kronecker(solve(sig), crossprod(xt)))%*%theta)
    DF <- length(theta)
    P <- pchisq(Q, df=DF, lower.tail=F)
  }
  if(type==3){
    n <- length(y)
    yt <- y - mean(y)
    dt<- yt^2 - sum(yt^2)/n
    xt <- scale(x, scale=F)
    g <- cbind(xt*yt, xt*dt); colnames(g) <- NULL
    gbar <- colMeans(g)
    v <- t(g)%*%g/n
    Q <- n*((t(gbar)%*%solve(v))%*%gbar)
    DF <- ncol(g)
    P <- pchisq(Q, df=DF, lower.tail=F)
  }
  if(type==4){
    n <- length(y)
    yt <- y - mean(y)
    ym <- abs(y - median(y))
    dt <- ym - mean(ym)
    xt <- scale(x, scale=F)
    g <- cbind(xt*yt, xt*dt); colnames(g) <- NULL
    gbar <- colMeans(g)
    v <- t(g)%*%g/n
    Q <- n*((t(gbar)%*%solve(v))%*%gbar)
    DF <- ncol(g)
    P <- pchisq(Q, df=DF, lower.tail=F)
  }
  
  # Results
  results <- data.frame(Q=Q, DF=DF, P=P)
  return(results)
  
}

#' jlsp
#'
#' jlsp performs the joint location-and-scale test using Fisher's method.
#' @param y vector of outcome values.
#' @param x vector of exposure values.
#' @param covar a data.frame of covariates.
#' @param covar.var adjust the second stage (variance component) of the approach by the covariates.
#' @param type type of test (default: 1 [Breusch-Pagan variance test]; options: 1 [Breusch-Pagan variance test], 2 [Brown-Forsythe variance test]).
#' @return a list of results. Q/F is the test statistic, DF is the degrees of freedom and P is the p-value. The model coefficients from each part of the model are given in the coef objects.
#' @examples
#' x <- rbinom(1000, 1, 0.5)
#' y <- 0.5 + 0.025*x + rnorm(1000, 0, sqrt(0.025^2*x)) + rnorm(1000, 0, 0.1)
#' jlsp(y, x, type=2)
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
jlsp <- function(y, x, covar=NULL, covar.var=FALSE, var.type=1){
  
  # Errors
  if(!(is.numeric(y) | is.integer(y))) stop("y has to be a numeric variable")
  if(!(is.numeric(x) | is.integer(x) | is.factor(x))) stop("x has to be either a numeric variable or a factor variable")
  if(!is.null(covar)){if(!is.data.frame(covar)) stop("covar has to be a data.frame")}
  if(length(y)!=length(x)) stop("y is not the same size as x")
  if(!is.null(covar)){if(length(y)!=nrow(covar)) stop("y is not the same size as covar")}
  
  # Missing values
  data <- cbind(y, x); if(!is.null(covar)){data <- cbind(data, covar)}
  keep <- complete.cases(data)
  y <- y[keep]; x <- x[keep]; if(!is.null(covar)){covar <- covar[keep,,drop=F]}
  
  # Covariates
  if(!is.null(covar)){
    covar <- model.matrix(as.formula(~ .), data=covar)[,-1,drop=F]
    if(any(is.na(covar))) stop("there are missing values in the covariates")
  }
  
  # Location test
  if(!is.null(covar)){ols <- lm(y~x+covar); ols0 <- lm(y~covar)}else{ols <- lm(y~x); ols0 <- lm(y~1)}
  coef <- summary(ols)$coefficients; rownames(coef) <- sub("covar", "", rownames(coef))
  test <- anova(ols0,ols)[2,c(3,5,6)]; names(test) <- c("DF", "F", "P"); rownames(test) <- 1:nrow(test)
  location_test <- list(coef=coef, test=test)
  p_location <- location_test$test$P
  
  # Scale test
  if(!is.null(covar)){
    if(var.type==1){d <- (resid(lm(y~x+covar)))^2}else{d <- suppressWarnings(abs(resid(rq(y~x+covar, tau=0.5))))}
  }else{
    if(var.type==1){d <- (resid(lm(y~x)))^2}else{d <- suppressWarnings(abs(resid(rq(y~x, tau=0.5))))}
  }
  if(!is.null(covar) & covar.var){mod <- lm(d~x+covar); mod0 <- lm(d~covar)}else{mod <- lm(d~x); mod0 <- lm(d~1)}
  coef <- summary(mod)$coefficients; rownames(coef) <- sub("covar", "", rownames(coef))
  test <- anova(mod0, mod)[2,c(3,5,6)]; names(test) <- c("DF", "F", "P"); rownames(test) <- 1:nrow(test)
  scale_test <- list(coef=coef, test=test)
  p_scale <- scale_test$test$P
  
  # Location + scale test
  Q <- -2*(log(p_location) + log(p_scale))
  DF <- 4
  P <- pchisq(Q, df=DF, lower.tail=F)
  location_scale_test <- data.frame(DF=DF, Q=Q, P=P)
  
  # Results
  results <- list(location_test=location_test, scale_test=scale_test, location_scale_test=location_scale_test)
  return(results)
  
}

#' vartest
#'
#' vartest performs variability tests by either the Breusch-Pagan or Brown-Forsythe methods.
#' @param y vector of outcome values.
#' @param x vector of exposure values.
#' @param covar a data.frame of covariates.
#' @param covar.var adjust the second stage (variance component) of the approach by the covariates.
#' @param type type of test (default: 1 [Breusch-Pagan variance test]; options: 1 [Breusch-Pagan variance test], 2 [Brown-Forsythe variance test]).
#' @return a list of results. F is the test statistic, DF is the degrees of freedom and P is the p-value. The model coefficients from variance part of the model are given in the coef object.
#' @examples
#' x <- rbinom(1000, 1, 0.5)
#' y <- 0.5 + 0.025*x + rnorm(1000, 0, sqrt(0.025^2*x)) + rnorm(1000, 0, 0.1)
#' vartest(y, x, type=2)
#' @author James R Staley <james.staley@bristol.ac.uk>
#' @export
vartest <- function(y, x, covar=NULL, covar.var=FALSE, type=1){
  
  # Errors
  if(!(is.numeric(y) | is.integer(y))) stop("y has to be a numeric variable")
  if(!(is.numeric(x) | is.integer(x) | is.factor(x))) stop("x has to be either a numeric variable or a factor variable")
  if(!is.null(covar)){if(!is.data.frame(covar)) stop("covar has to be a data.frame")}
  if(is.null(covar) & covar.var) stop("covar.var cannot be TRUE if there are no covariates")
  if(length(y)!=length(x)) stop("y is not the same size as x")
  if(!is.null(covar)){if(length(y)!=nrow(covar)) stop("y is not the same size as covar")}
  if(!(type %in% 1:2)) stop("type has to be set to either 1 or 2")
  
  # Missing values
  data <- cbind(y, x); if(!is.null(covar)){data <- cbind(data, covar)}
  keep <- complete.cases(data)
  y <- y[keep]; x <- x[keep]; if(!is.null(covar)){covar <- covar[keep,,drop=F]}
  
  # Covariates
  if(!is.null(covar)){
    covar <- model.matrix(as.formula(~ .), data=covar)[,-1,drop=F]
    if(any(is.na(covar))) stop("there are missing values in the covariates")
  }
  
  # Variance test
  if(!is.null(covar)){
    if(type==1){d <- (resid(lm(y~x+covar)))^2}else{d <- suppressWarnings(abs(resid(rq(y~x+covar, tau=0.5))))}
  }else{
    if(type==1){d <- (resid(lm(y~x)))^2}else{d <- suppressWarnings(abs(resid(rq(y~x, tau=0.5))))}
  }
  if(!is.null(covar) & covar.var){mod <- lm(d~x+covar); mod0 <- lm(d~covar)}else{mod <- lm(d~x); mod0 <- lm(d~1)}
  coef <- summary(mod)$coefficients; rownames(coef) <- sub("covar", "", rownames(coef))
  test <- anova(mod0, mod)[2,c(3,5,6)]; names(test) <- c("DF", "F", "P"); rownames(test) <- 1:nrow(test)
  results <- list(coef=coef, test=test)
  
  # Results  
  return(results)
  
}
