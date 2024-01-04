#' @title Joint location-and-scale test
#'
#' @description `jlssc` performs the joint location-and-scale score test.
#'
#' @name jlssc
#'
#' @param y vector of outcome values
#'
#' @param x vector of exposure values
#'
#' @param covar `data.frame` of covariates
#'
#' @param type type of test, where
#'  `1` = Breusch-Pagan variance test,
#'  `2` = Brown-Forsythe variance test,
#'  `3` = Method of moments version of test 1, and
#'  `4` = Method of moments version of test 2]
#'  (default: `1`)
#'
#' @param x.sq include x-squared in the model
#'
#' @param x.reg regress out the covariates from the exposure terms
#'
#' @return `jlst` returns a `data.frame` of results:
#'
#' @return \item{Q}{the test statistic}
#'
#' @return \item{DF}{the degrees of freedom}
#'
#' @return \item{P}{the p-value}
#'
#' @examples
#' x <- rbinom(1000, 1, 0.5)
#' y <- 0.5 + 0.025 * x + rnorm(1000, 0, sqrt(0.005 * x)) + rnorm(1000, 0, 0.1)
#' jlssc(y, x)
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
jlssc <- function(y, x, covar = NULL, type = 1, x.sq = FALSE, x.reg = TRUE) {

  # Errors
  if (!(is.numeric(y) || is.integer(y))) {
    stop("y has to be a numeric variable")
  }
  if (!(is.numeric(x) || is.integer(x) || is.factor(x))) {
    stop("x has to be either a numeric variable or a factor variable")
  }
  if (!is.null(covar) && !is.data.frame(covar)) {
    stop("covar has to be a data.frame")
  }
  if (length(y) != length(x)) {
    stop("y is not the same size as x")
  }
  if (!is.null(covar) && length(y) != nrow(covar)) {
    stop("y is not the same size as covar")
  }
  if (!(type %in% 1:4)) {
    stop("type has to be set to either 1, 2, 3 or 4")
  }
  if (is.logical(x.sq) == FALSE) {
    stop("x.sq has to logical")
  }
  if (x.sq == TRUE && is.factor(x) == TRUE) {
    stop("x.sq cannot be set to true if x is a factor variable")
  }

  # Missing values
  data <- cbind(y, x)
  if (!is.null(covar)) {
    data <- cbind(data, covar)
  }
  keep <- complete.cases(data)
  y <- y[keep]
  x <- x[keep]
  if (!is.null(covar)) {
    covar <- covar[keep, , drop = FALSE]
  }

  # Covariates
  if (!is.null(covar)) {
    covar <- model.matrix(as.formula(~.), data = covar)[, -1, drop = FALSE]
    if (any(is.na(covar))) {
      stop("there are missing values in the covariates")
    }
  }

  # Exposure
  if (is.factor(x)) {
    x <- model.matrix(as.formula(~ as.factor(x)))[, -1, drop = FALSE]
    if (any(is.na(x))) {
      stop("there are missing values in the exposure")
    }
  }
  if (x.sq == TRUE) {
    x2 <- x^2
    x <- model.matrix(as.formula(~ x + x2))[, -1, drop = FALSE]
    if (any(is.na(x))) stop("there are missing values in the exposure")
  }

  # Regress out covariates
  if (!is.null(covar)) {
    y <- lm(y ~ covar)$resid
    if (x.reg == TRUE) {
      x <- lm(x ~ covar)$resid
    }
  }

  # Location + scale test
  if (type == 1) {
    n <- length(y)
    yt <- y - mean(y)
    dt <- yt^2 - sum(yt^2) / n
    xt <- scale(x, scale = FALSE)
    beta <- lm(yt ~ xt - 1)$coef
    delta <- lm(dt ~ xt - 1)$coef
    theta <- c(beta, delta)
    sig <- matrix(
      c(mean(yt^2), mean(yt * dt), mean(yt * dt), mean(dt^2)),
      nrow = 2,
      ncol = 2
    )
    Q <- ((t(theta) %*% kronecker(solve(sig), crossprod(xt))) %*% theta)
    DF <- length(theta)
    P <- pchisq(Q, df = DF, lower.tail = FALSE)
  }
  if (type == 2) {
    n <- length(y)
    yt <- y - mean(y)
    ym <- abs(y - median(y))
    dt <- ym - mean(ym)
    xt <- scale(x, scale = FALSE)
    beta <- lm(yt ~ xt - 1)$coef
    delta <- lm(dt ~ xt - 1)$coef
    theta <- c(beta, delta)
    sig <- matrix(
      c(mean(yt^2), mean(yt * dt), mean(yt * dt), mean(dt^2)),
      nrow = 2,
      ncol = 2
    )
    Q <- ((t(theta) %*% kronecker(solve(sig), crossprod(xt))) %*% theta)
    DF <- length(theta)
    P <- pchisq(Q, df = DF, lower.tail = FALSE)
  }
  if (type == 3) {
    n <- length(y)
    yt <- y - mean(y)
    dt <- yt^2 - sum(yt^2) / n
    xt <- scale(x, scale = FALSE)
    g <- cbind(xt * yt, xt * dt)
    colnames(g) <- NULL
    gbar <- colMeans(g)
    v <- t(g) %*% g / n
    Q <- n * ((t(gbar) %*% solve(v)) %*% gbar)
    DF <- ncol(g)
    P <- pchisq(Q, df = DF, lower.tail = FALSE)
  }
  if (type == 4) {
    n <- length(y)
    yt <- y - mean(y)
    ym <- abs(y - median(y))
    dt <- ym - mean(ym)
    xt <- scale(x, scale = FALSE)
    g <- cbind(xt * yt, xt * dt)
    colnames(g) <- NULL
    gbar <- colMeans(g)
    v <- t(g) %*% g / n
    Q <- n * ((t(gbar) %*% solve(v)) %*% gbar)
    DF <- ncol(g)
    P <- pchisq(Q, df = DF, lower.tail = FALSE)
  }

  # Output
  output <- data.frame(Q = Q, DF = DF, P = P)

  # Return
  return(output)

}

#' @title Joint location-and-scale test using Fisher's method
#'
#' @description `jlsp` performs the joint location-and-scale test using
#'   Fisher's method.
#'
#' @name jlsp
#'
#' @importFrom quantreg rq
#'
#' @param y vector of outcome values
#'
#' @param x vector of exposure values
#'
#' @param covar `data.frame` of covariates
#'
#' @param covar.var adjust the second stage (variance component) of the approach
#'   by the covariates
#'
#' @param type type of test, where
#'  `1` = Breusch-Pagan variance test, and
#'  `2` = Brown-Forsythe variance test
#'  (default: `1`)
#'
#' @param x.sq include x-squared in the variance part of the model
#'
#' @return `jlsp` returns a `list` of results:
#'
#' @return \item{Q / F}{the test statistic}
#'
#' @return \item{DF}{the degrees of freedom}
#'
#' @return \item{P}{the p-value}
#'
#' @examples
#' x <- rbinom(1000, 1, 0.5)
#' y <- 0.5 + 0.025 * x + rnorm(1000, 0, sqrt(0.005 * x)) + rnorm(1000, 0, 0.1)
#' jlsp(y, x, var.type = 2)
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
jlsp <- function(y, x, covar = NULL, covar.var = FALSE,
  var.type = 1, x.sq = FALSE) {

  # Errors
  if (!(is.numeric(y) || is.integer(y))) {
    stop("y has to be a numeric variable")
  }
  if (!(is.numeric(x) || is.integer(x) || is.factor(x))) {
    stop("x has to be either a numeric variable or a factor variable")
  }
  if (!is.null(covar) && !is.data.frame(covar)) {
    stop("covar has to be a data.frame")
  }
  if (length(y) != length(x)) {
    stop("y is not the same size as x")
  }
  if (!is.null(covar) && length(y) != nrow(covar)) {
    stop("y is not the same size as covar")
  }
  if (is.logical(x.sq) == FALSE) {
    stop("x.sq has to logical")
  }
  if (x.sq == TRUE && is.factor(x) == TRUE) {
    stop("x.sq cannot be set to true if x is a factor variable")
  }

  # Missing values
  data <- cbind(y, x)
  if (!is.null(covar)) {
    data <- cbind(data, covar)
  }
  keep <- complete.cases(data)
  y <- y[keep]
  x <- x[keep]
  if (!is.null(covar)) {
    covar <- covar[keep, , drop = FALSE]
  }

  # Covariates
  if (!is.null(covar)) {
    covar <- model.matrix(as.formula(~.), data = covar)[, -1, drop = FALSE]
    if (any(is.na(covar))) {
      stop("there are missing values in the covariates")
    }
  }

  # Exposure squared
  if (x.sq == TRUE) {
    x2 <- x^2
  }

  # Location test
  if (!is.null(covar)) {
    ols <- lm(y ~ x + covar)
    ols0 <- lm(y ~ covar)
  } else {
    ols <- lm(y ~ x)
    ols0 <- lm(y ~ 1)
  }
  coef <- summary(ols)$coefficients
  rownames(coef) <- sub("covar", "", rownames(coef))
  test <- anova(ols0, ols)[2, c(3, 5, 6)]
  test <- as.data.frame(test)
  names(test) <- c("DF", "F", "P")
  rownames(test) <- seq_len(nrow(test))
  location_test <- list(coef = coef, test = test)
  p_location <- location_test$test$P

  # Scale test
  if (!is.null(covar)) {
    if (var.type == 1) {
      d <- (resid(lm(y ~ x + covar)))^2
    } else {
      d <- suppressWarnings(abs(resid(rq(y ~ x + covar, tau = 0.5))))
    }
  } else {
    if (var.type == 1) {
      d <- (resid(lm(y ~ x)))^2
    } else {
      d <- suppressWarnings(abs(resid(rq(y ~ x, tau = 0.5))))
    }
  }
  if (!is.null(covar) && covar.var) {
    if (x.sq == TRUE) {
      mod <- lm(d ~ x + x2 + covar)
      mod0 <- lm(d ~ covar)
    } else {
      mod <- lm(d ~ x + covar)
      mod0 <- lm(d ~ covar)
    }
  } else {
    if (x.sq == TRUE) {
      mod <- lm(d ~ x + x2)
      mod0 <- lm(d ~ 1)
    } else {
      mod <- lm(d ~ x)
      mod0 <- lm(d ~ 1)
    }
  }
  coef <- summary(mod)$coefficients
  rownames(coef) <- sub("covar", "", rownames(coef))
  test <- anova(mod0, mod)[2, c(3, 5, 6)]
  test <- as.data.frame(test)
  names(test) <- c("DF", "F", "P")
  rownames(test) <- seq_len(nrow(test))
  scale_test <- list(coef = coef, test = test)
  p_scale <- scale_test$test$P

  # Location + scale test
  Q <- -2 * (log(p_location) + log(p_scale))
  DF <- 4
  P <- pchisq(Q, df = DF, lower.tail = FALSE)
  location_scale_test <- data.frame(DF = DF, Q = Q, P = P)

  # Output
  output <- list(
    location_test = location_test,
    scale_test = scale_test,
    location_scale_test = location_scale_test
  )

  # Return
  return(output)

}

#' @title Variability tests
#'
#' @description vartest performs variability tests by either the Breusch-Pagan
#'   or Brown-Forsythe methods.
#'
#' @name vartest
#'
#' @importFrom quantreg rq
#'
#' @param y vector of outcome values
#'
#' @param x vector of exposure values
#'
#' @param covar `data.frame` of covariates
#'
#' @param covar.var adjust the second stage (variance component) of the approach
#'   by the covariates
#'
#' @param type type of test, where
#'  `1` = Breusch-Pagan variance test, and
#'  `2` = Brown-Forsythe variance test
#'  (default: `1`)
#'
#' @param x.sq include x-squared in the variance part of the model
#'
#' @return `vartest` returns a `list` of results:
#'
#' @return \item{coef}{model coefficients from variance part of the model}
#'
#' @return \item{test}{`data.frame` of test results}
#'
#' * `F`: the test statistic
#' * `DF`: the degrees of freedom
#' * `P`: the p-value
#'
#' @examples
#' x <- rbinom(1000, 1, 0.5)
#' y <- 0.5 + 0.025 * x + rnorm(1000, 0, sqrt(0.005 * x)) + rnorm(1000, 0, 0.1)
#' vartest(y, x, type = 2)
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
vartest <- function(y, x, covar = NULL, covar.var = FALSE, type = 1,
  x.sq = FALSE) {

  # Errors
  if (!(is.numeric(y) || is.integer(y))) {
    stop("y has to be a numeric variable")
  }
  if (!(is.numeric(x) || is.integer(x) || is.factor(x))) {
    stop("x has to be either a numeric variable or a factor variable")
  }
  if (!is.null(covar) && !is.data.frame(covar)) {
    stop("covar has to be a data.frame")
  }
  if (is.null(covar) && covar.var) {
    stop("covar.var cannot be TRUE if there are no covariates")
  }
  if (length(y) != length(x)) {
    stop("y is not the same size as x")
  }
  if (!is.null(covar) && length(y) != nrow(covar)) {
    stop("y is not the same size as covar")
  }
  if (!(type %in% 1:2)) {
    stop("type has to be set to either 1 or 2")
  }
  if (is.logical(x.sq) == FALSE) {
    stop("x.sq has to logical")
  }
  if (x.sq == TRUE && is.factor(x) == TRUE) {
    stop("x.sq cannot be set to true if x is a factor variable")
  }

  # Missing values
  data <- cbind(y, x)
  if (!is.null(covar)) {
    data <- cbind(data, covar)
  }
  keep <- complete.cases(data)
  y <- y[keep]
  x <- x[keep]
  if (!is.null(covar)) {
    covar <- covar[keep, , drop = FALSE]
  }

  # Covariates
  if (!is.null(covar)) {
    covar <- model.matrix(as.formula(~.), data = covar)[, -1, drop = FALSE]
    if (any(is.na(covar))) stop("there are missing values in the covariates")
  }

  # Exposure squared
  if (x.sq == TRUE) {
    x2 <- x^2
  }

  # Variance test
  if (!is.null(covar)) {
    if (type == 1) {
      d <- (resid(lm(y ~ x + covar)))^2
    } else {
      d <- suppressWarnings(abs(resid(rq(y ~ x + covar, tau = 0.5))))
    }
  } else {
    if (type == 1) {
      d <- (resid(lm(y ~ x)))^2
    } else {
      d <- suppressWarnings(abs(resid(rq(y ~ x, tau = 0.5))))
    }
  }
  if (!is.null(covar) && covar.var) {
    if (x.sq == TRUE) {
      mod <- lm(d ~ x + x2 + covar)
      mod0 <- lm(d ~ covar)
    } else {
      mod <- lm(d ~ x + covar)
      mod0 <- lm(d ~ covar)
    }
  } else {
    if (x.sq == TRUE) {
      mod <- lm(d ~ x + x2)
      mod0 <- lm(d ~ 1)
    } else {
      mod <- lm(d ~ x)
      mod0 <- lm(d ~ 1)
    }
  }
  coef <- summary(mod)$coefficients
  rownames(coef) <- sub("covar", "", rownames(coef))
  test <- anova(mod0, mod)[2, c(5, 3, 6)]
  test <- as.data.frame(test)
  names(test) <- c("F", "DF", "P")
  rownames(test) <- seq_len(nrow(test))
  output <- list(coef = coef, test = test)

  # Return
  return(output)

}
