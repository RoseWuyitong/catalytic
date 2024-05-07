# lmmlasso code from github.  reference: https://rdrr.io/cran/lmmlasso/
library(emulator)
library(penalized)
ArmijoRule <- function(xGroup, yGroup, LGroup, b, j, cut, HkOldJ, HkJ, JinNonpen, lambda, weights, nonpen, ll1, ll2, converged, control) {
  b.new <- b
  bJ <- b[j]
  grad <- -cut + bJ * HkOldJ
  if (JinNonpen) {
    dk <- -grad / HkJ
  } else {
    dk <- MedianValue(grad, HkJ, lambda / weights[j], bJ)
  }

  if (dk != 0) {
    # calculate delta_k
    if (JinNonpen) {
      deltak <- dk * grad + control$gamma * dk^2 * HkJ
    } else {
      deltak <- dk * grad + control$gamma * dk^2 * HkJ + lambda / weights[j] * (abs(bJ + dk) - abs(bJ))
    }

    fctOld <- ObjFunction(
      xGroup = xGroup, yGroup = yGroup, LGroup = LGroup, b = b, weights = weights,
      lambda = lambda, nonpen = nonpen, ll1 = ll1, ll2 = ll2
    )
    for (l in 0:control$max.armijo)
    {
      b.new[j] <- bJ + control$a_init * control$delta^l * dk

      fctNew <- ObjFunction(
        xGroup = xGroup, yGroup = yGroup, LGroup = LGroup, b = b.new, weights = weights,
        lambda = lambda, nonpen = nonpen, ll1 = ll1, ll2 = ll2
      )
      addDelta <- control$a_init * control$delta^l * control$rho * deltak
      if (fctNew <= fctOld + addDelta) {
        b[j] <- bJ + control$a_init * control$delta^l * dk
        fct <- fctNew
        break
      }
      if (l == control$max.armijo) {
        if (trace > 2) cat("Armijo for b_", j, " not successful", "\n")
        converged <- converged + 2
        fct <- fctOld
      }
    }
  }
  return(list(b = b, fct = fct, converged = converged))
}

HessianMatrix <- function(xGroup, LGroup, activeSet, N, hessian, mat) {
  for (i in 1:N)
  {
    mat[i, ] <- diag(t(xGroup[[i]][, activeSet, drop = FALSE]) %*% LGroup[[i]] %*% xGroup[[i]][, activeSet, drop = FALSE])
  }
  hessian[activeSet] <- apply(mat, 2, sum)
  return(hessian)
}
inverse <- function(X, eps = 1e-12) {
  eig.X <- eigen(X, symmetric = TRUE)
  P <- eig.X[[2]]
  lambda <- eig.X[[1]]

  ind <- lambda > eps

  lambda[ind] <- 1 / lambda[ind]
  lambda[!ind] <- 0

  P %*% diag(lambda) %*% t(P)
}
LambdaInv <- function(Z, ZId, Psi, sigma) {
  # cat("dim of Psi is ",dim(Psi), "dim of Z is ",dim(Z),'\n')
  # print(eigen(quad.tform(Psi,Z))$value)
  tau_square <- Psi[1, 1]
  LambdaInverse <- 1 / sigma^2 * ZId - 1 / sigma^4 * tau_square * Z %*% inverse(diag(ncol(Z)) + t(Z) %*% Z * tau_square / sigma^2) %*% t(Z)
  # inverse(sigma^2*ZId + quad.tform(Psi,Z))
  return(list(LambdaInverse))
}

MLloglik <- function(xGroup, yGroup, LGroup, b, ntot, N, activeSet) {
  l1 <- l2 <- numeric(N)
  ll2b <- 0
  for (i in 1:N)
  {
    l1[i] <- -determinant(LGroup[[i]])$modulus
    l2[i] <- quad.form(LGroup[[i]], yGroup[[i]] - xGroup[[i]] %*% b)
  }

  ll <- -1 / 2 * (sum(l1) + sum(l2) + ntot * log(2 * pi) + ll2b)
  return(loglik = ll)
}

MLpdDiag <- function(zGroup, resGroup, zIdGroup, sigma, pars, q, thres = 10^(-4), ll1, ll4, trace, CovOpt, VarInt) {
  LPsiIter <- diag(pars, nrow = length(pars))

  for (s in 1:q)
  {
    # optimization with optimize
    if (CovOpt == "optimize") {
      optRes <- optimize(MLpdSymFct,
        interval = VarInt, zGroup = zGroup, zIdGroup = zIdGroup, resGroup = resGroup,
        sigma = sigma, a = s, b = s, LPsi = LPsiIter
      )
      if (optRes$minimum < thres) {
        LPsiIter[s, s] <- 0
        if (trace > 2) cat("thres was used", "\n")
      } else {
        LPsiIter[s, s] <- optRes$minimum
      }
      # optimization with nlminb
    } else if (CovOpt == "nlminb") {
      optRes <- nlminb(pars[s], MLpdSymFct,
        zGroup = zGroup, zIdGroup = zIdGroup, resGroup = resGroup,
        sigma = sigma, a = s, b = s, LPsi = LPsiIter
      )
      if (optRes$par < thres) {
        LPsiIter[s, s] <- 0
        if (trace > 2) cat("thres was used", "\n")
      } else {
        LPsiIter[s, s] <- optRes$par
      }
    }

    ll <- ll1 + optRes$objective + ll4
    if (trace > 3) print(ll)
  }

  pars <- diag(LPsiIter)

  return(list(pars = pars, fct = ll))
}

MLpdIdent <- function(zGroup, resGroup, zIdGroup, sigma, pars, q, thres = 10^(-4), ll1, ll4, trace, CovOpt, VarInt) {
  LPsiIter <- diag(q)

  # optimization with optimize (no starting value)
  if (CovOpt == "optimize") {
    optRes <- optimize(MLpdIdentFct,
      interval = VarInt, zGroup = zGroup, zIdGroup = zIdGroup, resGroup = resGroup,
      sigma = sigma, LPsi = LPsiIter
    )
    if (exp(optRes$minimum) < thres) {
      pars <- 0
      if (trace > 2) cat("thres was used", "\n")
    } else {
      pars <- exp(optRes$minimum)
    }

    # optimization with nlminb (starting value provided)
  } else if (CovOpt == "nlminb") {
    optRes <- nlminb(log(pars), MLpdIdentFct,
      zGroup = zGroup, zIdGroup = zIdGroup, resGroup = resGroup,
      sigma = sigma, LPsi = LPsiIter
    )
    if (exp(optRes$par) < thres) {
      pars <- 0
      if (trace > 2) cat("thres was used", "\n")
    } else {
      pars <- exp(optRes$par)
    }
  }

  ll <- ll1 + optRes$objective + ll4
  if (trace > 3) print(ll)

  return(list(pars = pars, fct = ll))
}

MLpdIdentFct <- function(thetak, zGroup, zIdGroup, resGroup, sigma, LPsi) {
  Psi <- exp(thetak)^2 * LPsi

  LambdaGroup <- mapply(MLpdSymLambda, Z = zGroup, ZId = zIdGroup, MoreArgs = list(sigma = sigma, Psi = Psi))

  ll2 <- mapply(MLpdSymObj, LambdaGroup, resGroup)
  1 / 2 * sum(ll2)
}

MLpdSym <- function(zGroup, resGroup, zIdGroup, sigma, pars, q, thres = 10^(-4), ll1, ll4, trace, CovOpt, CovInt, VarInt) {
  LPsiIter <- t(triang(pars, q))

  # optimize over the variance parameters
  for (s in 1:q)
  {
    # optimization with optimize
    if (CovOpt == "optimize") {
      optRes <- optimize(MLpdSymFct,
        interval = VarInt, zGroup = zGroup, zIdGroup = zIdGroup, resGroup = resGroup,
        sigma = sigma, a = s, b = s, LPsi = LPsiIter
      )
      if (optRes$minimum < thres) {
        LPsiIter[s, s] <- 0
        if (trace > 2) cat("thres was used", "\n")
      } else {
        LPsiIter[s, s] <- optRes$minimum
      }
      # optimization with nlminb
    } else if (CovOpt == "nlminb") {
      optRes <- nlminb(LPsiIter[s, s], MLpdSymFct,
        zGroup = zGroup, zIdGroup = zIdGroup, resGroup = resGroup,
        sigma = sigma, a = s, b = s, LPsi = LPsiIter
      )
      if (optRes$par < thres) {
        LPsiIter[s, s] <- 0
        if (trace > 2) cat("thres was used", "\n")
      } else {
        LPsiIter[s, s] <- optRes$par
      }
    }

    ll <- ll1 + optRes$objective + ll4
    if (trace > 3) print(ll)
  }

  # optimize over the covariance parameters

  if (q > 1) {
    for (l in 2:q)
    {
      for (r in 1:(l - 1))
      {
        # optimization with optimize
        if (CovOpt == "optimize") {
          optRes <- optimize(MLpdSymFct,
            interval = CovInt, zGroup = zGroup, zIdGroup = zIdGroup, resGroup = resGroup,
            sigma = sigma, a = l, b = r, LPsi = LPsiIter
          )
          if (optRes$minimum < thres) {
            LPsiIter[r, l] <- 0
            if (trace > 2) cat("thres was used", "\n")
          } else {
            LPsiIter[l, r] <- optRes$minimum
          }
          # optimization with nlminb
        } else if (CovOpt == "nlminb") {
          optRes <- nlminb(LPsiIter[l, r], MLpdSymFct,
            zGroup = zGroup, zIdGroup = zIdGroup, resGroup = resGroup,
            sigma = sigma, a = l, b = r, LPsi = LPsiIter
          )
          if (optRes$par < thres) {
            LPsiIter[r, l] <- 0
            if (trace > 2) cat("thres was used", "\n")
          } else {
            LPsiIter[l, r] <- optRes$par
          }
        }
        ll <- ll1 + optRes$objective + ll4
        if (trace > 3) print(ll)
      }
    }
  }

  # return the parameters
  pars <- vecli(t(LPsiIter))

  return(list(pars = pars, fct = ll))
}

MLpdSymFct <- function(thetak, zGroup, zIdGroup, resGroup, sigma, a, b, LPsi) {
  LPsi[a, b] <- thetak
  Psi <- tcrossprod(LPsi)
  LambdaGroup <- mapply(MLpdSymLambda, Z = zGroup, ZId = zIdGroup, MoreArgs = list(sigma = sigma, Psi = Psi))
  cat("LambdaGroup no use", "\n")
  ll2 <- mapply(MLpdSymObj, LambdaGroup, resGroup)
  1 / 2 * sum(ll2)
}
MLpdSymLambda <- function(Z, ZId, Psi, sigma) {
  Lambda <- sigma^2 * ZId + quad.tform(Psi, Z)
  # print(Lambda)
  return(list(Lambda))
}

MLpdSymObj <- function(Lambda, res) {
  # print(Lambda)
  determinant(Lambda)$modulus + MyQuadFormInv(Lambda, res)
}
MLsigma <- function(zGroup, zIdGroup, resGroup, q, ll1, ll4, true.sigma, Psi, trace, CovOpt, VarInt) {
  # cat("true.sigma is ",true.sigma,"\n")
  ZPZtGroup <- lapply(zGroup, MLsigmaZPsiZt, Psi = Psi)

  if (CovOpt == "optimize") {
    optRes <- optimize(MLsigmaFct, interval = VarInt, ZPZtGroup = ZPZtGroup, zIdGroup = zIdGroup, resGroup = resGroup)
    sigma <- exp(optRes$minimum) # ; cat("after optimize the sigma is ",sigma,"\n")
  } else if (CovOpt == "nlminb") {
    optRes <- nlminb(log(true.sigma), MLsigmaFct, ZPZtGroup = ZPZtGroup, zIdGroup = zIdGroup, resGroup = resGroup, lower = -5, upper = 10)
    sigma <- exp(optRes$par) # ; cat("after nlminb the sigma is ",sigma,"\n")
  }

  ll <- ll1 + optRes$objective + ll4
  if (trace > 3) print(ll)

  return(list(sigma = sigma, fct = ll))
}

MLsigmaLambda <- function(ZId, ZPZt, sigma) {
  Lambda <- sigma^2 * ZId + ZPZt
  return(list(Lambda))
}

MLsigmaFct <- function(sigma, ZPZtGroup, zIdGroup, resGroup) {
  sigma <- exp(sigma)
  LambdaGroup <- mapply(MLsigmaLambda, zIdGroup, ZPZtGroup, MoreArgs = list(sigma = sigma))
  # cat("sigma is ",sigma,"\n")
  ll2 <- mapply(MLpdSymObj, LambdaGroup, resGroup)
  # cat("whether reach here\n")
  1 / 2 * sum(ll2)
}

MLsigmaZPsiZt <- function(Z, Psi) ZPsiZt <- quad.tform(Psi, Z)

MedianValue <- function(grad, hessian, lambda, bj) {
  median(c((lambda - grad) / hessian, -bj, (-lambda - grad) / hessian))
}

MyQuadFormInv <- function(M, x) {
  # print(M)
  crossprod(x, solve(M, x))
}

ObjFunction <- function(xGroup, yGroup, LGroup, b, weights, lambda, nonpen, resGroup = NULL, ll1 = ll1, ll2 = NULL) {
  ResAs <- function(x, y, b, activeSet) y - x[, activeSet, drop = FALSE] %*% b[activeSet, drop = FALSE]

  tResLRes <- function(L, res) 1 / 2 * quad.form(L, res)

  if (missing(resGroup)) {
    activeSet <- which(b != 0)
    resGroup <- mapply(ResAs, x = xGroup, y = yGroup, MoreArgs = list(b = b, activeSet = activeSet), SIMPLIFY = FALSE)
  }

  if (missing(ll2)) ll2 <- nlogdet(LGroup = LGroup)
  # cat("b is ",b,"\n","b[-nonpen] is ",b[-nonpen],"\n","weight is ",weights,"weights[-nonpen] is ",weights[-nonpen],"\n")
  ll <- ll1 + ll2 + sum(mapply(tResLRes, LGroup, resGroup)) + lambda * sum(abs(b[-nonpen]) / weights[-nonpen])

  return(Fc = ll)
}

ResAsSplit <- function(x, y, b, f, activeSet) {
  r <- y - x[, activeSet, drop = FALSE] %*% b[activeSet, drop = FALSE]
  resGroup <- split(r, f)
  return(resGroup)
}

SoftThreshold <- function(z, g) {
  sign(z) * max(abs(z) - g, 0)
}

ZIdentity <- function(Z) {
  ZId <- diag(dim(Z)[[1]])
  return(list(ZId))
}

as1 <- function(xGroup, LGroup, activeSet, N) {
  fs <- function(x, l, a) {
    l %*% x[, a]
  }

  SGroup <- mapply(fs, xGroup, LGroup, MoreArgs = list(a = activeSet), SIMPLIFY = FALSE)

  return(SGroup)
}

as2 <- function(x, y, b, j, activeSet, group, sGroup) {
  r <- y - x[, -c(j), drop = FALSE] %*% b[-j]
  rGroup <- split(r, group)

  as3 <- function(s, r, j) crossprod(r, s[, j])

  ma <- mapply(as3, sGroup, rGroup, MoreArgs = list(j = match(j, activeSet)))

  sumMa <- sum(ma)

  return(sumMa)
}

covStartingValues <- function(xGroup, yGroup, zGroup, zIdGroup, b, ntot, N, lower = -10, upper = 10) {
  optimize1 <- function(x, y, b) y - x %*% b

  optimize2 <- function(gamma, zId, ZtZ) {
    H <- zId + exp(2 * gamma) * ZtZ
    return(list(H = H))
  }

  optimize3 <- function(res, zId, ZtZ, gamma) {
    lambda <- optimize2(gamma, zId, ZtZ)
    logdetH <- determinant(lambda$H)$modulus
    quadH <- quad.form.inv(lambda$H, res)
    return(c(logdetH, quadH))
  }

  optimize4 <- function(gamma) {
    optH <- mapply(optimize3, resGroup, zIdGroup, ZtZ = ztzGroup, MoreArgs = list(gamma = gamma))
    H1 <- optH[1, ]
    H2 <- optH[2, ]

    fn <- ntot * log(sum(H2)) + sum(H1)
    fn
  }

  optimize5 <- function(z) tcrossprod(z)

  resGroup <- mapply(optimize1, x = xGroup, y = yGroup, MoreArgs = list(b = b), SIMPLIFY = FALSE)
  ztzGroup <- mapply(optimize5, z = zGroup, SIMPLIFY = FALSE)

  optRes <- optimize(f = optimize4, interval = c(lower, upper))

  gamma <- optRes$minimum

  quadH <- mapply(optimize3, resGroup, zIdGroup, ztzGroup, MoreArgs = list(gamma = gamma))[2, ]

  sig <- sqrt(1 / ntot * sum(quadH))
  tau <- exp(gamma) * sig
  objfct <- 1 / 2 * (optRes$objective + ntot * (1 - log(ntot)))

  return(list(tau = tau, sigma = sig, opt = objfct))
}

nlogdet <- function(LGroup) {
  nlogdetfun <- function(L) {
    -1 / 2 * determinant(L)$modulus[1]
  }

  sum(mapply(nlogdetfun, LGroup))
}


lmmlasso <- function(x, ...) {
  UseMethod("lmmlasso")
}

lmmlasso.default <- function(x, y, z = x, grp, weights = NULL, coefInit = NULL, lambda, startValue = 1, nonpen = 1, pdMat = c("pdIdent", "pdDiag", "pdSym"),
                             method = "ML", CovOpt = c("nlminb", "optimize"), stopSat = TRUE, standardize = TRUE, control = lmmlassoControl(),
                             ranInd = 1:dim(z)[[2]], ...) {
  # --- Introductory checks and reformulations ---
  # ----------------------------------------------

  CovOpt <- match.arg(CovOpt)
  pdMat <- match.arg(pdMat)

  # do some checks
  if (!is.matrix(x)) stop("x has to be a matrix")
  if (any(is.na(x))) stop("Missing values in x not allowed")
  if (any(is.na(y))) stop("Missing values in y not allowed")
  if (!is.matrix(z)) stop("z has to be a matrix")
  if (any(is.na(z))) stop("Missing values in z not allowed")
  if (!is.numeric(y)) stop("y has to be of type 'numeric'")
  if (nrow(x) != length(y)) stop("x and y have not correct dimensions")
  if (any(x[, 1] != rep(1, dim(x)[[1]]))) stop("first column is not the intercept")
  if (length(levels(grp)) == 1) stop("Only one group. No covariance parameters!")
  if (!all(nonpen %in% 1:dim(x)[[2]])) stop("Error with the argument nonpen")
  if (length(lambda) != 1) stop("lambda can only be one number")
  if (lambda < 0) stop("regularization parameter must be positive")

  # calculate weights (NA: not penalized, 0: drop)
  if (missing(weights)) {
    weights <- rep(1, dim(x)[[2]])
  } else {
    if (!length(weights) == dim(x)[[2]]) stop("Weights vector has not length p")

    nonpen <- which(is.na(weights))
    if (any(weights[-nonpen] < 0)) stop("Weights must be positive")
    remove <- weights == 0
    remove[nonpen] <- FALSE

    x <- x[, !remove, drop = FALSE]

    rem <- which(weights == 0)
    z <- z[, !(ranInd %in% rem), drop = FALSE]

    weights <- weights[!remove]
    nonpen <- which(is.na(weights))

    ran1 <- logical(dim(x)[[2]])
    ran1[ranInd] <- TRUE
    ran2 <- ran1[!remove]
    ranInd <- which(ran2)
  }

  if (standardize) {
    xOr <- x
    meanx <- apply(x[, -1], 2, mean)
    sdx <- apply(x[, -1], 2, sd)
    x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
  }

  # crucial allocations
  grp <- factor(grp)
  N <- length(levels(grp)) # N is the number of groups
  p <- dim(x)[[2]] # p is the number of covariates
  q <- dim(z)[[2]] # q is the number of random effects variables
  ntot <- length(y) # ntot is the total number of observations
  Q <- q * (q + 1) / 2 # maximum number of variance components parameters

  # save the grouped information as components of a list
  yGrp <- split(y, grp)
  xGrp <- split.data.frame(x, grp)
  zGrp <- split.data.frame(z, grp)
  zIdGrp <- mapply(ZIdentity, zGrp)

  # --- Calculating or adopting the starting values ---
  # ---------------------------------------------------

  if (missing(coefInit)) {
    # --- starting value for beta ---
    set.seed(control$seed)
    if (startValue == 1) # 10-fold cv Lasso
      {
        init <- optL1(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else if (startValue == 2) # 10-fold cv Ridge Regression
      {
        init <- optL2(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else # Null start
    {
      betaStart <- c(mean(y), rep(0, p - 1))
    }

    # --- starting values for the variance components parameters ---
    covStart <- covStartingValues(xGrp, yGrp, zGrp, zIdGrp, betaStart, ntot, N)
    sStart <- covStart$sigma
    tStart <- covStart$tau

    if (pdMat == "pdSym") parsStart <- vecli(tStart * diag(q))
    if (pdMat == "pdDiag") parsStart <- rep(tStart, q)
    if (pdMat == "pdIdent") parsStart <- tStart
  } else {
    betaStart <- coefInit[[1]]
    parsStart <- coefInit[[2]]
    sStart <- coefInit[[3]]
  }

  if (pdMat == "pdSym") PsiStart <- crossprod(triang(parsStart, q))
  if (pdMat == "pdDiag") PsiStart <- diag(parsStart^2, nrow = length(parsStart))
  if (pdMat == "pdIdent") PsiStart <- parsStart^2 * diag(q)

  # --- Calculate objective function for the starting values ---
  # ------------------------------------------------------------

  lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiStart, sigma = sStart))
  ll1 <- 1 / 2 * ntot * log(2 * pi)
  fctStart <- ObjFunction(xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaStart, weights = weights, lambda = lambda, nonpen = nonpen, ll1 = ll1)
  if (control$trace > 3) cat("fctStart:", fctStart, "\n")


  # --- Coordinate Gradient Descent-iteration ---
  # ---------------------------------------------

  # some necessary allocations:
  betaIter <- betaStart
  sIter <- sStart
  parsIter <- parsStart
  PsiIter <- PsiStart
  convPar <- crossprod(betaIter)
  convCov <- crossprod(c(sStart, parsStart))
  fctIter <- convFct <- fctStart
  hessian0 <- rep(0, p)
  mat0 <- matrix(0, ncol = p, nrow = N)

  stopped <- FALSE
  doAll <- FALSE
  converged <- 0
  counterIn <- 0
  counter <- 0 # counts the number of outer iterations

  while ((control$maxIter > counter) & (convPar > control$tol | convFct > control$tol | convCov > control$tol | !doAll)) {
    # print(sIter)
    if (control$maxIter == counter + 1) {
      if (control$trace > 2) cat("maxIter reached", "\n")
      converged <- converged + 1
    }
    # cat("PsiIter is ",dim(PsiIter),"\n")
    # arguments that change at each iteration
    counter <- counter + 1
    if (control$trace == 2) cat(counter, "...", "\n")
    betaIterOld <- betaIter
    fctIterOld <- fctIter
    covIterOld <- c(sIter, parsIter)

    # cat("current iteration is ",counter,"\n")
    # --- optimization w.r.t the fixed effects vector beta ---
    # --------------------------------------------------------

    activeSet <- 1:length(betaIter)
    if ((length(activeSet) > min(p, ntot)) & (lambda > 0) & (stopSat) & (counter > 2)) {
      stopped <- TRUE
      break
    }
    # here will have some problem if we are in high dimension-Li Weihao
    if (counterIn == 0 | counterIn > control$number) {
      doAll <- TRUE
      activeSet <- 1:p
      counterIn <- 1
    } else {
      doAll <- FALSE
      counterIn <- counterIn + 1
    }

    # calculate the hessian matrices for j in the activeSet
    HessIter <- HessIterTrunc <- HessianMatrix(xGroup = xGrp, LGroup = lambdaInvGrp, activeSet = activeSet, N = N, hessian = hessian0, mat = mat0[, activeSet, drop = FALSE])
    HessIter[activeSet] <- pmin(pmax(HessIter[activeSet], control$lower), control$upper)

    LxGrp <- as1(xGrp, lambdaInvGrp, activeSet, N = N)

    ll2 <- nlogdet(LGroup = lambdaInvGrp)

    for (j in activeSet)
    {
      cut1 <- as2(x = x, y = y, b = betaIter, j = j, activeSet = activeSet, group = grp, sGroup = LxGrp)
      JinNonpen <- j %in% nonpen

      # optimum can be calculated analytically
      if (HessIterTrunc[j] == HessIter[j]) {
        if (JinNonpen) {
          betaIter[j] <- cut1 / HessIter[j]
        } else {
          betaIter[j] <- SoftThreshold(cut1, lambda / weights[j]) / HessIter[j]
        }
      } else

      # optimimum is determined by the armijo rule
      {
        cat(" optimimum is determined by the armijo rule \n")
        armijo <- ArmijoRule(
          xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, j = j, cut = cut1, HkOldJ = HessIterTrunc[j],
          HkJ = HessIter[j], JinNonpen = JinNonpen, lambda = lambda, weights = weights, nonpen = nonpen,
          ll1 = ll1, ll2 = ll2, converged = converged, control = control
        )
        betaIter <- armijo$b
        converged <- armijo$converged
        fctIter <- armijo$fct
      }
      if (control$trace > 3) {
        fctIter <- ObjFunction(
          xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, weights = weights,
          lambda = lambda, nonpen = nonpen, ll1 = ll1, ll2 = ll2
        )
        print(fctIter)
      }
    }

    if (control$trace > 3) cat("-----------", "\n")

    # --- optimization w.r.t the variance components parameters ---
    # -------------------------------------------------------------

    # calculations before the covariance optimization
    activeSet <- 1:length(betaIter)
    resGrp <- ResAsSplit(x = x, y = y, b = betaIter, f = grp, activeSet = activeSet)
    ll4 <- lambda * sum(abs(betaIter[-nonpen]) / weights[-nonpen])

    # optimization of the free parameters in \Psi
    if (pdMat == "pdSym") {
      covParOpt <- MLpdSym(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt, CovInt = control$CovInt
      )
      PsiIter <- crossprod(triang(covParOpt$pars, q))
    } else if (pdMat == "pdDiag") {
      covParOpt <- MLpdDiag(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- diag(covParOpt$pars^2, nrow = length(covParOpt$pars))
    } else if (pdMat == "pdIdent") {
      covParOpt <- MLpdIdent(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- covParOpt$par^2 * diag(q)
    }

    parsIter <- covParOpt$pars
    fctIter <- covParOpt$fct
    if (control$trace > 3) cat("-----------", "\n")

    # optimization of the error variance \sigma^2
    covParOpt <- MLsigma(
      zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, q = q, ll1 = ll1, ll4 = ll4, true.sigma = sIter, Psi = PsiIter,
      trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
    )
    sIter <- covParOpt$sigma # ; if(sIter<0.0001) sIter=covIterOld[1];print(sIter)
    fctIter <- covParOpt$fct
    # cat("sigma is ",sIter,"\n")
    if (control$trace > 3) cat("-----------", "\n") # cat("obj.fct III:",fctIter,"\n")

    lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiIter, sigma = sIter))
    covIter <- c(sIter, parsIter)

    # --- check convergence ---
    convPar <- sqrt(crossprod(betaIter - betaIterOld)) / (1 + sqrt(crossprod(betaIter)))
    convFct <- abs((fctIterOld - fctIter) / (1 + abs(fctIter)))
    convCov <- sqrt(crossprod(covIter - covIterOld)) / (1 + sqrt(crossprod(covIter)))

    if ((convPar <= control$tol) & (convFct <= control$tol) & (convCov <= control$tol)) counterIn <- 0
  }

  if (standardize) {
    betaIter[-1] <- betaIter[-1] / sdx
    betaIter[1] <- betaIter[1] - sum(meanx * betaIter[-1])
    x <- xOr
    xGrp <- split.data.frame(xOr, grp)
  }


  # --- prediction of the random effects ---
  # ----------------------------------------

  biGroup <- u <- list()
  length(biGroup) <- length(u) <- N
  Psi <- PsiIter
  # cat("psi is ",diag(Psi),"\n")
  corPsi <- cov2cor(Psi)

  if (pdMat == "pdSym") {
    cholPsi <- t(triang(parsIter, q))
  } else if (pdMat == "pdDiag") {
    cholPsi <- diag(parsIter^2, nrow = length(parsIter))
  } else if (pdMat == "pdIdent") cholPsi <- parsIter^2 * diag(q)

  for (i in 1:N)
  {
    u[[i]] <- sIter * ginv(t(zGrp[[i]] %*% cholPsi) %*% zGrp[[i]] %*% cholPsi + sIter^2 * diag(q)) %*% t(zGrp[[i]] %*% cholPsi) %*% resGrp[[i]]
    biGroup[[i]] <- 1 / sIter * cholPsi %*% u[[i]]
  }


  #  --- final calculations ---
  # ---------------------------

  # fitted values and residuals
  residGrp <- fittedGrp <- list()
  length(residGrp) <- length(fittedGrp) <- N
  for (i in 1:N)
  {
    fittedGrp[[i]] <- xGrp[[i]][, activeSet, drop = FALSE] %*% betaIter[activeSet, drop = FALSE] + zGrp[[i]] %*% biGroup[[i]]
    residGrp[[i]] <- yGrp[[i]] - fittedGrp[[i]]
  }
  residuals <- unlist(residGrp)
  fitted <- unlist(fittedGrp)

  # random effects, sorted per subject
  u <- unlist(u) # corresponds to lmer@u
  bi <- unlist(biGroup) # unsorted random effects

  # random effects, sorted per effect
  ranef <- bi[order(rep(1:q, N))] # corresponds to lmer@ranef

  # fixed effects without names
  fixef <- betaIter
  names(fixef) <- NULL

  # --- summary information ---
  # ---------------------------
  npar <- sum(betaIter != 0) + length(c(sIter, parsIter))
  logLik <- MLloglik(xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, ntot = ntot, N = N, activeSet = which(betaIter != 0))
  deviance <- -2 * logLik
  aic <- -2 * logLik + 2 * npar
  bic <- -2 * logLik + log(ntot) * npar

  if (any(parsIter == 0)) cat("Redundant covariance parameters.", "\n")
  if (converged > 0) cat("Algorithm does not properly converge.", "\n")
  if (stopped) {
    cat("|activeSet|>=min(p,ntot): Increase lambda or set stopSat=FALSE.", "\n")
    sIter <- parsIter <- nlogLik <- aic <- bic <- NA
    betaIter <- rep(NA, p)
    bi <- fitted <- residuals <- NULL
  }

  out <- list(
    data = list(x = x, y = y, z = z, grp = grp), weights = weights, coefInit = list(betaStart = betaStart, parsStart = parsStart, sStart = sStart),
    lambda = lambda, sigma = abs(sIter), pars = parsIter, coefficients = betaIter, random = bi, u = u, ranef = ranef, fixef = fixef, fitted.values = fitted,
    residuals = residuals, Psi = Psi, corPsi = corPsi, converged = converged, logLik = logLik, npar = npar, deviance = deviance,
    aic = aic, bic = bic, nonpen = nonpen, counter = counter, pdMat = pdMat, method = method, CovOpt = CovOpt, control = control, call = match.call(),
    stopped = stopped, ranInd = ranInd, objective = fctIter
  )
  out
  structure(out, class = "lmmlasso")
}

lmmlassoControl <- function(tol = 10^(-5), trace = 1, maxIter = 800, maxArmijo = 20, number = 5, a_init = 1,
                            delta = 0.1, rho = 0.001, gamma = 0, lower = 10^(-4), upper = 10^8, seed = 532, VarInt = c(-5, 5), CovInt = c(-5, 5), thres = 10^(-6)) {
  list(
    tol = tol, trace = trace, maxIter = maxIter, max.armijo = maxArmijo, number = number,
    a_init = a_init, delta = delta, rho = rho, gamma = gamma, lower = lower, upper = upper, seed = seed, VarInt = VarInt, CovInt = CovInt, thres = thres
  )
}

# lmmridge modification from lmmlasso
ObjFunction_ridge <- function(xGroup, yGroup, LGroup, b, weights, lambda, nonpen, resGroup = NULL, ll1 = ll1, ll2 = NULL) {
  ResAs <- function(x, y, b, activeSet) y - x[, activeSet, drop = FALSE] %*% b[activeSet, drop = FALSE]

  tResLRes <- function(L, res) 1 / 2 * quad.form(L, res)

  if (missing(resGroup)) {
    activeSet <- which(b != 0)
    resGroup <- mapply(ResAs, x = xGroup, y = yGroup, MoreArgs = list(b = b, activeSet = activeSet), SIMPLIFY = FALSE)
  }

  if (missing(ll2)) ll2 <- nlogdet(LGroup = LGroup)

  ll <- ll1 + ll2 + sum(mapply(tResLRes, LGroup, resGroup)) + lambda * sum((b[-nonpen])^2 / weights[-nonpen]) # only this change--weihao

  return(Fc = ll)
}

HessianMatrix_ridge <- function(xGroup, LGroup, activeSet, N, hessian, mat, lambda) {
  for (i in 1:N)
  {
    mat[i, ] <- diag(t(xGroup[[i]][, activeSet, drop = FALSE]) %*% LGroup[[i]] %*% xGroup[[i]][, activeSet, drop = FALSE]) #+2*lambda*diag(rep(1,ncol(xGroup[[i]]))))
  }
  hessian[activeSet] <- apply(mat, 2, sum)
  return(hessian)
}
ArmijoRule_ridge <- function(xGroup, yGroup, LGroup, b, j, cut, HkOldJ, HkJ, JinNonpen, lambda, weights, nonpen, ll1, ll2, converged, control) {
  b.new <- b
  bJ <- b[j]
  grad <- -cut + bJ * HkOldJ
  if (JinNonpen) {
    dk <- -grad / HkJ
  } else {
    dk <- (-2 * lambda / weights[j] * bJ - grad) / (HkJ + 2 * lambda / weights[j])
  }

  if (dk != 0) {
    # calculate delta_k
    if (JinNonpen) {
      deltak <- dk * grad + control$gamma * dk^2 * HkJ
    } else {
      deltak <- dk * grad + control$gamma * dk^2 * HkJ + lambda / weights[j] * ((bJ + dk)^2 - (bJ)^2)
    }

    fctOld <- ObjFunction_ridge(
      xGroup = xGroup, yGroup = yGroup, LGroup = LGroup, b = b, weights = weights,
      lambda = lambda, nonpen = nonpen, ll1 = ll1, ll2 = ll2
    )
    for (l in 0:control$max.armijo)
    {
      # b.new[j] <- bJ + control$a_init*control$delta^l*dk
      b.new[j] <- bJ + control$a_init * dk

      fctNew <- ObjFunction_ridge(
        xGroup = xGroup, yGroup = yGroup, LGroup = LGroup, b = b.new, weights = weights,
        lambda = lambda, nonpen = nonpen, ll1 = ll1, ll2 = ll2
      )
      # addDelta <- control$a_init*control$delta^l*control$rho*deltak
      addDelta <- control$a_init * control$rho * deltak
      if (fctNew <= fctOld + addDelta) {
        b[j] <- bJ + control$a_init * control$delta^l * dk
        fct <- fctNew
        break
      }
      if (l == control$max.armijo) {
        if (trace > 2) cat("Armijo for b_", j, " not successful", "\n")
        converged <- converged + 2
        fct <- fctOld
      }
    }
  }
  return(list(b = b, fct = fct, converged = converged))
}

lmmridge <- function(x, y, z = x, grp, weights = NULL, coefInit = NULL, lambda, startValue = 2, nonpen = 1, pdMat = c("pdIdent", "pdDiag", "pdSym"),
                     method = "ML", CovOpt = c("nlminb", "optimize"), stopSat = TRUE, standardize = FALSE, control = lmmlassoControl(),
                     ranInd = 1:dim(z)[[2]], ...) {
  # --- Introductory checks and reformulations ---
  # ----------------------------------------------

  CovOpt <- match.arg(CovOpt)
  pdMat <- match.arg(pdMat)

  # do some checks
  if (!is.matrix(x)) stop("x has to be a matrix")
  if (any(is.na(x))) stop("Missing values in x not allowed")
  if (any(is.na(y))) stop("Missing values in y not allowed")
  if (!is.matrix(z)) stop("z has to be a matrix")
  if (any(is.na(z))) stop("Missing values in z not allowed")
  if (!is.numeric(y)) stop("y has to be of type 'numeric'")
  if (nrow(x) != length(y)) stop("x and y have not correct dimensions")
  if (any(x[, 1] != rep(1, dim(x)[[1]]))) stop("first column is not the intercept")
  if (length(levels(grp)) == 1) stop("Only one group. No covariance parameters!")
  if (!all(nonpen %in% 1:dim(x)[[2]])) stop("Error with the argument nonpen")
  if (length(lambda) != 1) stop("lambda can only be one number")
  if (lambda < 0) stop("regularization parameter must be positive")

  # calculate weights (NA: not penalized, 0: drop)
  if (missing(weights)) {
    weights <- rep(1, dim(x)[[2]])
  } else {
    if (!length(weights) == dim(x)[[2]]) stop("Weights vector has not length p")

    nonpen <- which(is.na(weights))
    if (any(weights[-nonpen] < 0)) stop("Weights must be positive")
    remove <- weights == 0
    remove[nonpen] <- FALSE

    x <- x[, !remove, drop = FALSE]

    rem <- which(weights == 0)
    z <- z[, !(ranInd %in% rem), drop = FALSE]

    weights <- weights[!remove]
    nonpen <- which(is.na(weights))

    ran1 <- logical(dim(x)[[2]])
    ran1[ranInd] <- TRUE
    ran2 <- ran1[!remove]
    ranInd <- which(ran2)
  }

  if (standardize) {
    xOr <- x
    meanx <- apply(x[, -1], 2, mean)
    sdx <- apply(x[, -1], 2, sd)
    x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
  }

  # crucial allocations
  grp <- factor(grp)
  N <- length(levels(grp)) # N is the number of groups
  p <- dim(x)[[2]] # p is the number of covariates
  q <- dim(z)[[2]] # q is the number of random effects variables
  ntot <- length(y) # ntot is the total number of observations
  Q <- q * (q + 1) / 2 # maximum number of variance components parameters

  # save the grouped information as components of a list
  yGrp <- split(y, grp)
  xGrp <- split.data.frame(x, grp)
  zGrp <- split.data.frame(z, grp)
  zIdGrp <- mapply(ZIdentity, zGrp)

  # --- Calculating or adopting the starting values ---
  # ---------------------------------------------------

  if (missing(coefInit)) {
    # --- starting value for beta ---
    set.seed(control$seed)
    if (startValue == 1) # 5-fold cv Lasso
      {
        init <- optL1(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else if (startValue == 2) # 5-fold cv Ridge Regression
      {
        init <- optL2(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else # Null start
    {
      betaStart <- c(mean(y), rep(0, p - 1))
    }

    # --- starting values for the variance components parameters ---
    covStart <- covStartingValues(xGrp, yGrp, zGrp, zIdGrp, betaStart, ntot, N)
    sStart <- covStart$sigma
    tStart <- covStart$tau

    if (pdMat == "pdSym") parsStart <- vecli(tStart * diag(q))
    if (pdMat == "pdDiag") parsStart <- rep(tStart, q)
    if (pdMat == "pdIdent") parsStart <- tStart
  } else {
    betaStart <- coefInit[[1]]
    parsStart <- coefInit[[2]]
    sStart <- coefInit[[3]]
  }

  if (pdMat == "pdSym") PsiStart <- crossprod(triang(parsStart, q))
  if (pdMat == "pdDiag") PsiStart <- diag(parsStart^2, nrow = length(parsStart))
  if (pdMat == "pdIdent") PsiStart <- parsStart^2 * diag(q)

  # --- Calculate objective function for the starting values ---
  # ------------------------------------------------------------

  lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiStart, sigma = sStart))
  ll1 <- 1 / 2 * ntot * log(2 * pi)
  fctStart <- ObjFunction_ridge(xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaStart, weights = weights, lambda = lambda, nonpen = nonpen, ll1 = ll1)



  # --- Coordinate Gradient Descent-iteration ---
  # ---------------------------------------------

  # some necessary allocations:
  betaIter <- betaStart
  sIter <- sStart
  parsIter <- parsStart
  PsiIter <- PsiStart
  convPar <- crossprod(betaIter)
  convCov <- crossprod(c(sStart, parsStart))
  fctIter <- convFct <- fctStart
  hessian0 <- rep(0, p)
  mat0 <- matrix(0, ncol = p, nrow = N)

  stopped <- FALSE
  doAll <- FALSE
  converged <- 0
  counterIn <- 0
  counter <- 0 # counts the number of outer iterations

  while ((control$maxIter > counter) & (convPar > control$tol | convFct > control$tol | convCov > control$tol | !doAll)) {
    if (control$maxIter == counter + 1) {
      if (control$trace > 2) cat("maxIter reached", "\n")
      converged <- converged + 1
    }

    # arguments that change at each iteration
    counter <- counter + 1
    if (control$trace == 2) cat(counter, "...", "\n")
    betaIterOld <- betaIter
    fctIterOld <- fctIter
    covIterOld <- c(sIter, parsIter)

    # --- optimization w.r.t the fixed effects vector beta ---
    # --------------------------------------------------------
    # weihao: analytical solution for betaIter, but we need to figure out Psi=PsiIter,sigma=sIter
    "AA=diag(c(0,rep(lambda,p-1)))
    n_i=nrow(X_list[[1]])
    BB=matrix(0,nc=1,nr=p)
    for(i in 1:m){
      Xi=X_list[[i]];Yi=Y_list[[i]];Zi=Z_list[[i]]
      Vi_inv=ginv(sIter^2*diag(n_i)+quad.tform(PsiIter,Zi))
      AA=AA+ quad.form(Vi_inv,Xi)

      BB=BB+ t(Xi)%*%Vi_inv%*%Yi
    }
    betaIter=as.numeric( ginv(AA)%*%BB)"
    # above code does not converge, we need to use original code

    activeSet <- 1:length(betaIter)
    if ((length(activeSet) > min(p, ntot)) & (lambda > 0) & (stopSat) & (counter > 2)) {
      stopped <- TRUE
      break
    }
    # here will have some problem if we are in high dimension-Li Weihao
    if (counterIn == 0 | counterIn > control$number) {
      doAll <- TRUE
      activeSet <- 1:p
      counterIn <- 1
    } else {
      doAll <- FALSE
      counterIn <- counterIn + 1
    }

    # calculate the hessian matrices for j in the activeSet
    HessIter <- HessIterTrunc <- HessianMatrix_ridge(xGroup = xGrp, LGroup = lambdaInvGrp, activeSet = activeSet, N = N, hessian = hessian0, mat = mat0[, activeSet, drop = FALSE], lambda)
    HessIter[activeSet] <- pmin(pmax(HessIter[activeSet], control$lower), control$upper)

    LxGrp <- as1(xGrp, lambdaInvGrp, activeSet, N = N)

    ll2 <- nlogdet(LGroup = lambdaInvGrp)

    for (j in activeSet)
    {
      cut1 <- as2(x = x, y = y, b = betaIter, j = j, activeSet = activeSet, group = grp, sGroup = LxGrp)
      JinNonpen <- j %in% nonpen

      # optimum can be calculated analytically
      if (HessIterTrunc[j] == HessIter[j]) {
        if (JinNonpen) {
          betaIter[j] <- cut1 / HessIter[j]
        } else {
          betaIter[j] <- cut1 / (HessIter[j] + 2 * lambda / weights[j])
        }
      } else

      # optimimum is determined by the armijo rule
      {
        print("optimimum is determined by the armijo rule")
        armijo <- ArmijoRule_ridge(
          xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, j = j, cut = cut1, HkOldJ = HessIterTrunc[j],
          HkJ = HessIter[j], JinNonpen = JinNonpen, lambda = lambda, weights = weights, nonpen = nonpen,
          ll1 = ll1, ll2 = ll2, converged = converged, control = control
        )
        betaIter <- armijo$b
        converged <- armijo$converged
        fctIter <- armijo$fct
      }
      if (control$trace > 3) {
        fctIter <- ObjFunction_ridge(
          xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, weights = weights,
          lambda = lambda, nonpen = nonpen, ll1 = ll1, ll2 = ll2
        )
        print(fctIter)
      }
    }

    if (control$trace > 3) cat("-----------", "\n")


    # --- optimization w.r.t the variance components parameters ---
    # -------------------------------------------------------------

    # calculations before the covariance optimization
    activeSet <- 1:length(betaIter)
    resGrp <- ResAsSplit(x = x, y = y, b = betaIter, f = grp, activeSet = activeSet)
    ll4 <- lambda * sum((betaIter[-nonpen])^2 / weights[-nonpen])

    # optimization of the free parameters in \Psi
    if (pdMat == "pdSym") {
      covParOpt <- MLpdSym(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt, CovInt = control$CovInt
      )
      PsiIter <- crossprod(triang(covParOpt$pars, q))
    } else if (pdMat == "pdDiag") {
      covParOpt <- MLpdDiag(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- diag(covParOpt$pars^2, nrow = length(covParOpt$pars))
    } else if (pdMat == "pdIdent") {
      covParOpt <- MLpdIdent(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- covParOpt$par^2 * diag(q)
    }

    parsIter <- covParOpt$pars
    fctIter <- covParOpt$fct
    if (control$trace > 3) cat("-----------", "\n")

    # optimization of the error variance \sigma^2
    covParOpt <- MLsigma(
      zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, q = q, ll1 = ll1, ll4 = ll4, true.sigma = sIter, Psi = PsiIter,
      trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
    )
    sIter <- covParOpt$sigma # ; if(sIter<0.0001) sIter=covIterOld[1]
    fctIter <- covParOpt$fct



    lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiIter, sigma = sIter))
    covIter <- c(sIter, parsIter)

    # --- check convergence ---
    convPar <- sqrt(crossprod(betaIter - betaIterOld)) / (1 + sqrt(crossprod(betaIter)))
    convFct <- abs((fctIterOld - fctIter) / (1 + abs(fctIter)))
    convCov <- sqrt(crossprod(covIter - covIterOld)) / (1 + sqrt(crossprod(covIter)))

    if ((convPar <= control$tol) & (convFct <= control$tol) & (convCov <= control$tol)) counterIn <- 0
  }

  if (standardize) {
    betaIter[-1] <- betaIter[-1] / sdx
    betaIter[1] <- betaIter[1] - sum(meanx * betaIter[-1])
    x <- xOr
    xGrp <- split.data.frame(xOr, grp)
  }


  # --- prediction of the random effects ---
  # ----------------------------------------

  biGroup <- u <- list()
  length(biGroup) <- length(u) <- N
  Psi <- PsiIter
  corPsi <- cov2cor(Psi)

  if (pdMat == "pdSym") {
    cholPsi <- t(triang(parsIter, q))
  } else if (pdMat == "pdDiag") {
    cholPsi <- diag(parsIter, nrow = length(parsIter))
  } else if (pdMat == "pdIdent") cholPsi <- parsIter * diag(q)

  for (i in 1:N)
  {
    u[[i]] <- sIter * ginv(t(zGrp[[i]] %*% cholPsi) %*% zGrp[[i]] %*% cholPsi + sIter^2 * diag(q)) %*% t(zGrp[[i]] %*% cholPsi) %*% resGrp[[i]]
    biGroup[[i]] <- 1 / sIter * cholPsi %*% u[[i]]
  }


  #  --- final calculations ---
  # ---------------------------

  # fitted values and residuals
  residGrp <- fittedGrp <- list()
  length(residGrp) <- length(fittedGrp) <- N
  for (i in 1:N)
  {
    fittedGrp[[i]] <- xGrp[[i]][, activeSet, drop = FALSE] %*% betaIter[activeSet, drop = FALSE] + zGrp[[i]] %*% biGroup[[i]]
    residGrp[[i]] <- yGrp[[i]] - fittedGrp[[i]]
  }
  residuals <- unlist(residGrp)
  fitted <- unlist(fittedGrp)

  # random effects, sorted per subject
  u <- unlist(u) # corresponds to lmer@u
  bi <- unlist(biGroup) # unsorted random effects

  # random effects, sorted per effect
  ranef <- bi[order(rep(1:q, N))] # corresponds to lmer@ranef

  # fixed effects without names
  fixef <- betaIter
  names(fixef) <- NULL

  # --- summary information ---
  # ---------------------------
  npar <- sum(betaIter != 0) + length(c(sIter, parsIter))
  logLik <- MLloglik(xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, ntot = ntot, N = N, activeSet = which(betaIter != 0))
  deviance <- -2 * logLik
  aic <- -2 * logLik + 2 * npar
  bic <- -2 * logLik + log(ntot) * npar

  if (any(parsIter == 0)) cat("Redundant covariance parameters.", "\n")
  if (converged > 0) cat("Algorithm does not properly converge.", "\n")
  if (stopped) {
    cat("|activeSet|>=min(p,ntot): Increase lambda or set stopSat=FALSE.", "\n")
    sIter <- parsIter <- nlogLik <- aic <- bic <- NA
    betaIter <- rep(NA, p)
    bi <- fitted <- residuals <- NULL
  }

  out <- list(
    data = list(x = x, y = y, z = z, grp = grp), weights = weights, coefInit = list(betaStart = betaStart, parsStart = parsStart, sStart = sStart),
    lambda = lambda, sigma = abs(sIter), pars = parsIter, coefficients = betaIter, random = bi, u = u, ranef = ranef, fixef = fixef, fitted.values = fitted,
    residuals = residuals, Psi = Psi, corPsi = corPsi, converged = converged, logLik = logLik, npar = npar, deviance = deviance,
    aic = aic, bic = bic, nonpen = nonpen, counter = counter, pdMat = pdMat, method = method, CovOpt = CovOpt, control = control, call = match.call(),
    stopped = stopped, ranInd = ranInd, objective = fctIter
  )

  out
  structure(out, class = "lmmlasso")
}






# here I delete the Armijo rule
lmmridge2 <- function(x, y, z = x, grp, weights = NULL, coefInit = NULL, lambda, startValue = 2, nonpen = 1, pdMat = c("pdIdent", "pdDiag", "pdSym"),
                      method = "ML", CovOpt = c("nlminb", "optimize"), stopSat = TRUE, standardize = FALSE, control = lmmlassoControl(),
                      ranInd = 1:dim(z)[[2]], ...) {
  # --- Introductory checks and reformulations ---
  # ----------------------------------------------

  CovOpt <- match.arg(CovOpt)
  pdMat <- match.arg(pdMat)

  # do some checks
  if (!is.matrix(x)) stop("x has to be a matrix")
  if (any(is.na(x))) stop("Missing values in x not allowed")
  if (any(is.na(y))) stop("Missing values in y not allowed")
  if (!is.matrix(z)) stop("z has to be a matrix")
  if (any(is.na(z))) stop("Missing values in z not allowed")
  if (!is.numeric(y)) stop("y has to be of type 'numeric'")
  if (nrow(x) != length(y)) stop("x and y have not correct dimensions")
  if (any(x[, 1] != rep(1, dim(x)[[1]]))) stop("first column is not the intercept")
  if (length(levels(grp)) == 1) stop("Only one group. No covariance parameters!")
  if (!all(nonpen %in% 1:dim(x)[[2]])) stop("Error with the argument nonpen")
  if (length(lambda) != 1) stop("lambda can only be one number")
  if (lambda < 0) stop("regularization parameter must be positive")

  # calculate weights (NA: not penalized, 0: drop)
  if (missing(weights)) {
    weights <- rep(1, dim(x)[[2]])
  } else {
    if (!length(weights) == dim(x)[[2]]) stop("Weights vector has not length p")

    nonpen <- which(is.na(weights))
    if (any(weights[-nonpen] < 0)) stop("Weights must be positive")
    remove <- weights == 0
    remove[nonpen] <- FALSE

    x <- x[, !remove, drop = FALSE]

    rem <- which(weights == 0)
    z <- z[, !(ranInd %in% rem), drop = FALSE]

    weights <- weights[!remove]
    nonpen <- which(is.na(weights))

    ran1 <- logical(dim(x)[[2]])
    ran1[ranInd] <- TRUE
    ran2 <- ran1[!remove]
    ranInd <- which(ran2)
  }

  if (standardize) {
    xOr <- x
    meanx <- apply(x[, -1], 2, mean)
    sdx <- apply(x[, -1], 2, sd)
    x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
  }

  # crucial allocations
  grp <- factor(grp)
  N <- length(levels(grp)) # N is the number of groups
  p <- dim(x)[[2]] # p is the number of covariates
  q <- dim(z)[[2]] # q is the number of random effects variables
  ntot <- length(y) # ntot is the total number of observations
  Q <- q * (q + 1) / 2 # maximum number of variance components parameters

  # save the grouped information as components of a list
  yGrp <- split(y, grp)
  xGrp <- split.data.frame(x, grp)
  zGrp <- split.data.frame(z, grp)
  zIdGrp <- mapply(ZIdentity, zGrp)

  # --- Calculating or adopting the starting values ---
  # ---------------------------------------------------

  if (missing(coefInit)) {
    # --- starting value for beta ---
    set.seed(control$seed)
    if (startValue == 1) # 5-fold cv Lasso
      {
        init <- optL1(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else if (startValue == 2) # 5-fold cv Ridge Regression
      {
        init <- optL2(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else # Null start
    {
      betaStart <- c(mean(y), rep(0, p - 1))
    }

    # --- starting values for the variance components parameters ---
    covStart <- covStartingValues(xGrp, yGrp, zGrp, zIdGrp, betaStart, ntot, N)
    sStart <- covStart$sigma
    tStart <- covStart$tau

    if (pdMat == "pdSym") parsStart <- vecli(tStart * diag(q))
    if (pdMat == "pdDiag") parsStart <- rep(tStart, q)
    if (pdMat == "pdIdent") parsStart <- tStart
  } else {
    betaStart <- coefInit[[1]]
    parsStart <- coefInit[[2]]
    sStart <- coefInit[[3]]
  }

  if (pdMat == "pdSym") PsiStart <- crossprod(triang(parsStart, q))
  if (pdMat == "pdDiag") PsiStart <- diag(parsStart^2, nrow = length(parsStart))
  if (pdMat == "pdIdent") PsiStart <- parsStart^2 * diag(q)

  # --- Calculate objective function for the starting values ---
  # ------------------------------------------------------------

  lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiStart, sigma = sStart))
  ll1 <- 1 / 2 * ntot * log(2 * pi)
  fctStart <- ObjFunction_ridge(xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaStart, weights = weights, lambda = lambda, nonpen = nonpen, ll1 = ll1)



  # --- Coordinate Gradient Descent-iteration ---
  # ---------------------------------------------

  # some necessary allocations:
  betaIter <- betaStart
  sIter <- sStart
  parsIter <- parsStart
  PsiIter <- PsiStart
  convPar <- crossprod(betaIter)
  convCov <- crossprod(c(sStart, parsStart))
  fctIter <- convFct <- fctStart
  hessian0 <- rep(0, p)
  mat0 <- matrix(0, ncol = p, nrow = N)

  stopped <- FALSE
  doAll <- FALSE
  converged <- 0
  counterIn <- 0
  counter <- 0 # counts the number of outer iterations

  while ((control$maxIter > counter) & (convPar > control$tol | convFct > control$tol | convCov > control$tol | !doAll)) {
    if (control$maxIter == counter + 1) {
      if (control$trace > 2) cat("maxIter reached", "\n")
      converged <- converged + 1
    }

    # arguments that change at each iteration
    counter <- counter + 1
    if (control$trace == 2) cat(counter, "...", "\n")
    betaIterOld <- betaIter
    fctIterOld <- fctIter
    covIterOld <- c(sIter, parsIter)

    # --- optimization w.r.t the fixed effects vector beta ---
    # --------------------------------------------------------



    activeSet <- 1:length(betaIter)
    # if ((length(activeSet)>min(p,ntot))&(lambda>0)&(stopSat)&(counter>2)) {stopped <- TRUE ; break}
    # here will have some problem if we are in high dimension-Li Weihao
    if (counterIn == 0 | counterIn > control$number) {
      doAll <- TRUE
      activeSet <- 1:p
      counterIn <- 1
    } else {
      doAll <- FALSE
      counterIn <- counterIn + 1
    }

    # calculate the hessian matrices for j in the activeSet
    HessIter <- HessIterTrunc <- HessianMatrix_ridge(xGroup = xGrp, LGroup = lambdaInvGrp, activeSet = activeSet, N = N, hessian = hessian0, mat = mat0[, activeSet, drop = FALSE], lambda)
    HessIter[activeSet] <- pmin(pmax(HessIter[activeSet], control$lower), control$upper)

    LxGrp <- as1(xGrp, lambdaInvGrp, activeSet, N = N)

    ll2 <- nlogdet(LGroup = lambdaInvGrp)

    for (j in activeSet)
    {
      cut1 <- as2(x = x, y = y, b = betaIter, j = j, activeSet = activeSet, group = grp, sGroup = LxGrp)
      JinNonpen <- j %in% nonpen

      # optimum can be calculated analytically

      if (JinNonpen) {
        betaIter[j] <- cut1 / HessIter[j]
      } else {
        betaIter[j] <- cut1 / (HessIter[j] + 2 * lambda / weights[j])
      }

      if (control$trace > 3) {
        fctIter <- ObjFunction_ridge(
          xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, weights = weights,
          lambda = lambda, nonpen = nonpen, ll1 = ll1, ll2 = ll2
        )
        print(fctIter)
      }
    }

    if (control$trace > 3) cat("-----------", "\n")


    # --- optimization w.r.t the variance components parameters ---
    # -------------------------------------------------------------

    # calculations before the covariance optimization
    activeSet <- 1:length(betaIter)
    resGrp <- ResAsSplit(x = x, y = y, b = betaIter, f = grp, activeSet = activeSet)
    ll4 <- lambda * sum((betaIter[-nonpen])^2 / weights[-nonpen])

    # optimization of the free parameters in \Psi
    if (pdMat == "pdSym") {
      covParOpt <- MLpdSym(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt, CovInt = control$CovInt
      )
      PsiIter <- crossprod(triang(covParOpt$pars, q))
    } else if (pdMat == "pdDiag") {
      covParOpt <- MLpdDiag(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- diag(covParOpt$pars^2, nrow = length(covParOpt$pars))
    } else if (pdMat == "pdIdent") {
      covParOpt <- MLpdIdent(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- covParOpt$par^2 * diag(q)
    }

    parsIter <- covParOpt$pars
    fctIter <- covParOpt$fct
    if (control$trace > 3) cat("-----------", "\n")

    # optimization of the error variance \sigma^2
    covParOpt <- MLsigma(
      zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, q = q, ll1 = ll1, ll4 = ll4, true.sigma = sIter, Psi = PsiIter,
      trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
    )
    sIter <- covParOpt$sigma # ; if(sIter<0.0001) sIter=covIterOld[1]
    fctIter <- covParOpt$fct



    lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiIter, sigma = sIter))
    covIter <- c(sIter, parsIter)

    # --- check convergence ---
    convPar <- sqrt(crossprod(betaIter - betaIterOld)) / (1 + sqrt(crossprod(betaIter)))
    convFct <- abs((fctIterOld - fctIter) / (1 + abs(fctIter)))
    convCov <- sqrt(crossprod(covIter - covIterOld)) / (1 + sqrt(crossprod(covIter)))

    if ((convPar <= control$tol) & (convFct <= control$tol) & (convCov <= control$tol)) counterIn <- 0
  }

  if (standardize) {
    betaIter[-1] <- betaIter[-1] / sdx
    betaIter[1] <- betaIter[1] - sum(meanx * betaIter[-1])
    x <- xOr
    xGrp <- split.data.frame(xOr, grp)
  }


  # --- prediction of the random effects ---
  # ----------------------------------------

  biGroup <- u <- list()
  length(biGroup) <- length(u) <- N
  Psi <- PsiIter
  corPsi <- cov2cor(Psi)

  if (pdMat == "pdSym") {
    cholPsi <- t(triang(parsIter, q))
  } else if (pdMat == "pdDiag") {
    cholPsi <- diag(parsIter^2, nrow = length(parsIter))
  } else if (pdMat == "pdIdent") cholPsi <- parsIter^2 * diag(q)

  for (i in 1:N)
  {
    u[[i]] <- sIter * ginv(t(zGrp[[i]] %*% cholPsi) %*% zGrp[[i]] %*% cholPsi + sIter^2 * diag(q)) %*% t(zGrp[[i]] %*% cholPsi) %*% resGrp[[i]]
    biGroup[[i]] <- 1 / sIter * cholPsi %*% u[[i]]
  }


  #  --- final calculations ---
  # ---------------------------

  # fitted values and residuals
  residGrp <- fittedGrp <- list()
  length(residGrp) <- length(fittedGrp) <- N
  for (i in 1:N)
  {
    fittedGrp[[i]] <- xGrp[[i]][, activeSet, drop = FALSE] %*% betaIter[activeSet, drop = FALSE] + zGrp[[i]] %*% biGroup[[i]]
    residGrp[[i]] <- yGrp[[i]] - fittedGrp[[i]]
  }
  residuals <- unlist(residGrp)
  fitted <- unlist(fittedGrp)

  # random effects, sorted per subject
  u <- unlist(u) # corresponds to lmer@u
  bi <- unlist(biGroup) # unsorted random effects

  # random effects, sorted per effect
  ranef <- bi[order(rep(1:q, N))] # corresponds to lmer@ranef

  # fixed effects without names
  fixef <- betaIter
  names(fixef) <- NULL

  # --- summary information ---
  # ---------------------------
  npar <- sum(betaIter != 0) + length(c(sIter, parsIter))
  logLik <- MLloglik(xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, ntot = ntot, N = N, activeSet = which(betaIter != 0))
  deviance <- -2 * logLik
  aic <- -2 * logLik + 2 * npar
  bic <- -2 * logLik + log(ntot) * npar

  if (any(parsIter == 0)) cat("Redundant covariance parameters.", "\n")
  if (converged > 0) cat("Algorithm does not properly converge.", "\n")
  if (stopped) {
    cat("|activeSet|>=min(p,ntot): Increase lambda or set stopSat=FALSE.", "\n")
    sIter <- parsIter <- nlogLik <- aic <- bic <- NA
    betaIter <- rep(NA, p)
    bi <- fitted <- residuals <- NULL
  }

  out <- list(
    data = list(x = x, y = y, z = z, grp = grp), weights = weights, coefInit = list(betaStart = betaStart, parsStart = parsStart, sStart = sStart),
    lambda = lambda, sigma = abs(sIter), pars = parsIter, coefficients = betaIter, random = bi, u = u, ranef = ranef, fixef = fixef, fitted.values = fitted,
    residuals = residuals, Psi = Psi, corPsi = corPsi, converged = converged, logLik = logLik, npar = npar, deviance = deviance,
    aic = aic, bic = bic, nonpen = nonpen, counter = counter, pdMat = pdMat, method = method, CovOpt = CovOpt, control = control, call = match.call(),
    stopped = stopped, ranInd = ranInd, objective = fctIter
  )

  out
  structure(out, class = "lmmlasso")
}




































# this function is designed for cv, give 4/5 data and 1/5 validation data, do prediction on test data and compute residual on validation set
lmmridge.pred_testdata <- function(x, y, z = x, grp, weights = NULL, coefInit = NULL, lambda, startValue = 2, nonpen = 1, pdMat = c("pdIdent", "pdDiag", "pdSym"),
                                   method = "ML", CovOpt = c("nlminb", "optimize"), stopSat = TRUE, standardize = FALSE, control = lmmlassoControl(),
                                   ranInd = 1:dim(z)[[2]], test_x, test_y, test_z, test_grp, ...) {
  # --- Introductory checks and reformulations ---
  # ----------------------------------------------

  CovOpt <- match.arg(CovOpt)
  pdMat <- match.arg(pdMat)

  # do some checks
  if (!is.matrix(x)) stop("x has to be a matrix")
  if (any(is.na(x))) stop("Missing values in x not allowed")
  if (any(is.na(y))) stop("Missing values in y not allowed")
  if (!is.matrix(z)) stop("z has to be a matrix")
  if (any(is.na(z))) stop("Missing values in z not allowed")
  if (!is.numeric(y)) stop("y has to be of type 'numeric'")
  if (nrow(x) != length(y)) stop("x and y have not correct dimensions")
  if (any(x[, 1] != rep(1, dim(x)[[1]]))) stop("first column is not the intercept")
  if (length(levels(grp)) == 1) stop("Only one group. No covariance parameters!")
  if (!all(nonpen %in% 1:dim(x)[[2]])) stop("Error with the argument nonpen")
  if (length(lambda) != 1) stop("lambda can only be one number")
  if (lambda < 0) stop("regularization parameter must be positive")

  # calculate weights (NA: not penalized, 0: drop)
  if (missing(weights)) {
    weights <- rep(1, dim(x)[[2]])
  } else {
    if (!length(weights) == dim(x)[[2]]) stop("Weights vector has not length p")

    nonpen <- which(is.na(weights))
    if (any(weights[-nonpen] < 0)) stop("Weights must be positive")
    remove <- weights == 0
    remove[nonpen] <- FALSE

    x <- x[, !remove, drop = FALSE]

    rem <- which(weights == 0)
    z <- z[, !(ranInd %in% rem), drop = FALSE]

    weights <- weights[!remove]
    nonpen <- which(is.na(weights))

    ran1 <- logical(dim(x)[[2]])
    ran1[ranInd] <- TRUE
    ran2 <- ran1[!remove]
    ranInd <- which(ran2)
  }

  if (standardize) {
    xOr <- x
    meanx <- apply(x[, -1], 2, mean)
    sdx <- apply(x[, -1], 2, sd)
    x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
  }

  # crucial allocations
  grp <- factor(grp)
  N <- length(levels(grp)) # N is the number of groups
  p <- dim(x)[[2]] # p is the number of covariates
  q <- dim(z)[[2]] # q is the number of random effects variables
  ntot <- length(y) # ntot is the total number of observations
  Q <- q * (q + 1) / 2 # maximum number of variance components parameters

  # save the grouped information as components of a list
  yGrp <- split(y, grp)
  xGrp <- split.data.frame(x, grp)
  zGrp <- split.data.frame(z, grp)
  zIdGrp <- mapply(ZIdentity, zGrp)

  yGrp_test <- split(test_y, test_grp)
  xGrp_test <- split.data.frame(test_x, test_grp)
  zGrp_test <- split.data.frame(test_z, test_grp)
  # --- Calculating or adopting the starting values ---
  # ---------------------------------------------------

  if (missing(coefInit)) {
    # --- starting value for beta ---
    set.seed(control$seed)
    if (startValue == 1) # 5-fold cv Lasso
      {
        init <- optL1(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else if (startValue == 2) # 5-fold cv Ridge Regression
      {
        init <- optL2(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else # Null start
    {
      betaStart <- c(mean(y), rep(0, p - 1))
    }

    # --- starting values for the variance components parameters ---
    covStart <- covStartingValues(xGrp, yGrp, zGrp, zIdGrp, betaStart, ntot, N)
    sStart <- covStart$sigma
    tStart <- covStart$tau

    if (pdMat == "pdSym") parsStart <- vecli(tStart * diag(q))
    if (pdMat == "pdDiag") parsStart <- rep(tStart, q)
    if (pdMat == "pdIdent") parsStart <- tStart
  } else {
    betaStart <- coefInit[[1]]
    parsStart <- coefInit[[2]]
    sStart <- coefInit[[3]]
  }

  if (pdMat == "pdSym") PsiStart <- crossprod(triang(parsStart, q))
  if (pdMat == "pdDiag") PsiStart <- diag(parsStart^2, nrow = length(parsStart))
  if (pdMat == "pdIdent") PsiStart <- parsStart^2 * diag(q)

  # --- Calculate objective function for the starting values ---
  # ------------------------------------------------------------

  lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiStart, sigma = sStart))
  ll1 <- 1 / 2 * ntot * log(2 * pi)
  fctStart <- ObjFunction_ridge(xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaStart, weights = weights, lambda = lambda, nonpen = nonpen, ll1 = ll1)



  # --- Coordinate Gradient Descent-iteration ---
  # ---------------------------------------------

  # some necessary allocations:
  betaIter <- betaStart
  sIter <- sStart
  parsIter <- parsStart
  PsiIter <- PsiStart
  convPar <- crossprod(betaIter)
  convCov <- crossprod(c(sStart, parsStart))
  fctIter <- convFct <- fctStart
  hessian0 <- rep(0, p)
  mat0 <- matrix(0, ncol = p, nrow = N)

  stopped <- FALSE
  doAll <- FALSE
  converged <- 0
  counterIn <- 0
  counter <- 0 # counts the number of outer iterations

  while ((control$maxIter > counter) & (convPar > control$tol | convFct > control$tol | convCov > control$tol | !doAll)) {
    if (control$maxIter == counter + 1) {
      if (control$trace > 2) cat("maxIter reached", "\n")
      converged <- converged + 1
    }

    # arguments that change at each iteration
    counter <- counter + 1
    if (control$trace == 2) cat(counter, "...", "\n")
    betaIterOld <- betaIter
    fctIterOld <- fctIter
    covIterOld <- c(sIter, parsIter)

    # --- optimization w.r.t the fixed effects vector beta ---
    # --------------------------------------------------------
    # weihao: analytical solution for betaIter, but we need to figure out Psi=PsiIter,sigma=sIter
    "AA=diag(c(0,rep(lambda,p-1)))
    n_i=nrow(X_list[[1]])
    BB=matrix(0,nc=1,nr=p)
    for(i in 1:m){
      Xi=X_list[[i]];Yi=Y_list[[i]];Zi=Z_list[[i]]
      Vi_inv=ginv(sIter^2*diag(n_i)+quad.tform(PsiIter,Zi))
      AA=AA+ quad.form(Vi_inv,Xi)

      BB=BB+ t(Xi)%*%Vi_inv%*%Yi
    }
    betaIter=as.numeric( ginv(AA)%*%BB)"
    # above code does not converge, we need to use original code

    activeSet <- which(betaIter != 0)
    if ((length(activeSet) > min(p, ntot)) & (lambda > 0) & (stopSat) & (counter > 2)) {
      stopped <- TRUE
      break
    }
    # here will have some problem if we are in high dimension-Li Weihao
    if (counterIn == 0 | counterIn > control$number) {
      doAll <- TRUE
      activeSet <- 1:p
      counterIn <- 1
    } else {
      doAll <- FALSE
      counterIn <- counterIn + 1
    }

    # calculate the hessian matrices for j in the activeSet
    HessIter <- HessIterTrunc <- HessianMatrix_ridge(xGroup = xGrp, LGroup = lambdaInvGrp, activeSet = activeSet, N = N, hessian = hessian0, mat = mat0[, activeSet, drop = FALSE], lambda)
    HessIter[activeSet] <- pmin(pmax(HessIter[activeSet], control$lower), control$upper)

    LxGrp <- as1(xGrp, lambdaInvGrp, activeSet, N = N)

    ll2 <- nlogdet(LGroup = lambdaInvGrp)

    for (j in activeSet)
    {
      cut1 <- as2(x = x, y = y, b = betaIter, j = j, activeSet = activeSet, group = grp, sGroup = LxGrp)
      JinNonpen <- j %in% nonpen

      # optimum can be calculated analytically
      if (HessIterTrunc[j] == HessIter[j]) {
        if (JinNonpen) {
          betaIter[j] <- cut1 / HessIter[j]
        } else {
          betaIter[j] <- cut1 / (HessIter[j] + 2 * lambda / weights[j])
        }
      } else

      # optimimum is determined by the armijo rule
      {
        print("optimimum is determined by the armijo rule")
        armijo <- ArmijoRule_ridge(
          xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, j = j, cut = cut1, HkOldJ = HessIterTrunc[j],
          HkJ = HessIter[j], JinNonpen = JinNonpen, lambda = lambda, weights = weights, nonpen = nonpen,
          ll1 = ll1, ll2 = ll2, converged = converged, control = control
        )
        betaIter <- armijo$b
        converged <- armijo$converged
        fctIter <- armijo$fct
      }
      if (control$trace > 3) {
        fctIter <- ObjFunction_ridge(
          xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, weights = weights,
          lambda = lambda, nonpen = nonpen, ll1 = ll1, ll2 = ll2
        )
        print(fctIter)
      }
    }

    if (control$trace > 3) cat("-----------", "\n")


    # --- optimization w.r.t the variance components parameters ---
    # -------------------------------------------------------------

    # calculations before the covariance optimization
    activeSet <- which(betaIter != 0)
    resGrp <- ResAsSplit(x = x, y = y, b = betaIter, f = grp, activeSet = activeSet)
    ll4 <- lambda * sum((betaIter[-nonpen])^2 / weights[-nonpen])

    # optimization of the free parameters in \Psi
    if (pdMat == "pdSym") {
      covParOpt <- MLpdSym(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt, CovInt = control$CovInt
      )
      PsiIter <- crossprod(triang(covParOpt$pars, q))
    } else if (pdMat == "pdDiag") {
      covParOpt <- MLpdDiag(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- diag(covParOpt$pars^2, nrow = length(covParOpt$pars))
    } else if (pdMat == "pdIdent") {
      covParOpt <- MLpdIdent(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- covParOpt$par^2 * diag(q)
    }

    parsIter <- covParOpt$pars
    fctIter <- covParOpt$fct
    if (control$trace > 3) cat("-----------", "\n")

    # optimization of the error variance \sigma^2
    covParOpt <- MLsigma(
      zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, q = q, ll1 = ll1, ll4 = ll4, true.sigma = sIter, Psi = PsiIter,
      trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
    )
    sIter <- covParOpt$sigma
    fctIter <- covParOpt$fct



    lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiIter, sigma = sIter))
    covIter <- c(sIter, parsIter)

    # --- check convergence ---
    convPar <- sqrt(crossprod(betaIter - betaIterOld)) / (1 + sqrt(crossprod(betaIter)))
    convFct <- abs((fctIterOld - fctIter) / (1 + abs(fctIter)))
    convCov <- sqrt(crossprod(covIter - covIterOld)) / (1 + sqrt(crossprod(covIter)))

    if ((convPar <= control$tol) & (convFct <= control$tol) & (convCov <= control$tol)) counterIn <- 0
  }

  if (standardize) {
    betaIter[-1] <- betaIter[-1] / sdx
    betaIter[1] <- betaIter[1] - sum(meanx * betaIter[-1])
    x <- xOr
    xGrp <- split.data.frame(xOr, grp)
  }


  # --- prediction of the random effects ---
  # ----------------------------------------

  biGroup <- u <- list()
  length(biGroup) <- length(u) <- N
  Psi <- PsiIter
  corPsi <- cov2cor(Psi)

  if (pdMat == "pdSym") {
    cholPsi <- t(triang(parsIter, q))
  } else if (pdMat == "pdDiag") {
    cholPsi <- diag(parsIter, nrow = length(parsIter))
  } else if (pdMat == "pdIdent") cholPsi <- parsIter * diag(q)

  for (i in 1:N)
  {
    u[[i]] <- sIter * solve(t(zGrp[[i]] %*% cholPsi) %*% zGrp[[i]] %*% cholPsi + sIter^2 * diag(q)) %*% t(zGrp[[i]] %*% cholPsi) %*% resGrp[[i]]
    biGroup[[i]] <- 1 / sIter * cholPsi %*% u[[i]]
  }


  #  --- final calculations ---
  # ---------------------------

  # fitted values and residuals
  residGrp <- fittedGrp <- list()
  length(residGrp) <- length(fittedGrp) <- N
  # cat("betaIter is",dim(zGrp_test[[1]]),length(biGroup[[i]]))
  for (i in 1:N)
  {
    fittedGrp[[i]] <- xGrp_test[[i]][, activeSet, drop = FALSE] %*% betaIter[activeSet, drop = FALSE] + zGrp_test[[i]] %*% biGroup[[i]]
    residGrp[[i]] <- yGrp_test[[i]] - fittedGrp[[i]]
  }
  residuals <- unlist(residGrp)
  fitted <- unlist(fittedGrp)

  # random effects, sorted per subject
  u <- unlist(u) # corresponds to lmer@u
  bi <- unlist(biGroup) # unsorted random effects

  # random effects, sorted per effect
  ranef <- bi[order(rep(1:q, N))] # corresponds to lmer@ranef

  # fixed effects without names
  fixef <- betaIter
  names(fixef) <- NULL

  # --- summary information ---
  # ---------------------------
  npar <- sum(betaIter != 0) + length(c(sIter, parsIter))
  logLik <- MLloglik(xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, ntot = ntot, N = N, activeSet = which(betaIter != 0))
  deviance <- -2 * logLik
  aic <- -2 * logLik + 2 * npar
  bic <- -2 * logLik + log(ntot) * npar

  if (any(parsIter == 0)) cat("Redundant covariance parameters.", "\n")
  if (converged > 0) cat("Algorithm does not properly converge.", "\n")
  if (stopped) {
    cat("|activeSet|>=min(p,ntot): Increase lambda or set stopSat=FALSE.", "\n")
    sIter <- parsIter <- nlogLik <- aic <- bic <- NA
    betaIter <- rep(NA, p)
    bi <- fitted <- residuals <- NULL
  }

  out <- list(
    data = list(x = x, y = y, z = z, grp = grp), weights = weights, coefInit = list(betaStart = betaStart, parsStart = parsStart, sStart = sStart),
    lambda = lambda, sigma = abs(sIter), pars = parsIter, coefficients = betaIter, random = bi, u = u, ranef = ranef, fixef = fixef, fitted.values = fitted,
    residuals = residuals, Psi = Psi, corPsi = corPsi, converged = converged, logLik = logLik, npar = npar, deviance = deviance,
    aic = aic, bic = bic, nonpen = nonpen, counter = counter, pdMat = pdMat, method = method, CovOpt = CovOpt, control = control, call = match.call(),
    stopped = stopped, ranInd = ranInd, objective = fctIter
  )

  out
  structure(out, class = "lmmlasso")
}


# #
ObjFunction_syn <- function(xGroup, yGroup, LGroup, b, weights, lambda, nonpen, resGroup = NULL, ll1 = ll1, ll2 = NULL, xGroup.syn, yGroup.syn, LGroup.syn) {
  ResAs <- function(x, y, b, activeSet) y - x[, activeSet, drop = FALSE] %*% b[activeSet, drop = FALSE]

  tResLRes <- function(L, res) 1 / 2 * quad.form(L, res)

  if (missing(resGroup)) {
    activeSet <- which(b != 0)
    resGroup <- mapply(ResAs, x = xGroup, y = yGroup, MoreArgs = list(b = b, activeSet = activeSet), SIMPLIFY = FALSE)
    resGroup.syn <- mapply(ResAs, x = xGroup.syn, y = yGroup.syn, MoreArgs = list(b = b, activeSet = activeSet), SIMPLIFY = FALSE)
  }

  if (missing(ll2)) ll2 <- nlogdet(LGroup = LGroup)

  ll <- ll1 + ll2 + sum(mapply(tResLRes, LGroup, resGroup)) + sum(mapply(tResLRes, LGroup.syn, resGroup.syn)) # only this change--weihao

  return(Fc = ll)
}

# this function is designed for synthetic version of LMM, with a given tau value
lmmsyn <- function(x, y, z = x, grp, weights = NULL, coefInit = NULL, lambda, startValue = 1, nonpen = 1, pdMat = c("pdIdent", "pdDiag", "pdSym"),
                   method = "ML", CovOpt = c("nlminb", "optimize"), stopSat = TRUE, standardize = TRUE, control = lmmlassoControl(),
                   ranInd = 1:dim(z)[[2]], my_syndata, tau_syn = 0.1, ...) {
  # --- Introductory checks and reformulations ---
  # ----------------------------------------------

  CovOpt <- match.arg(CovOpt)
  pdMat <- match.arg(pdMat)

  # do some checks
  if (!is.matrix(x)) stop("x has to be a matrix")
  if (any(is.na(x))) stop("Missing values in x not allowed")
  if (any(is.na(y))) stop("Missing values in y not allowed")
  if (!is.matrix(z)) stop("z has to be a matrix")
  if (any(is.na(z))) stop("Missing values in z not allowed")
  if (!is.numeric(y)) stop("y has to be of type 'numeric'")
  if (nrow(x) != length(y)) stop("x and y have not correct dimensions")
  if (any(x[, 1] != rep(1, dim(x)[[1]]))) stop("first column is not the intercept")
  if (length(levels(grp)) == 1) stop("Only one group. No covariance parameters!")
  if (!all(nonpen %in% 1:dim(x)[[2]])) stop("Error with the argument nonpen")
  if (length(lambda) != 1) stop("lambda can only be one number")
  if (lambda < 0) stop("regularization parameter must be positive")

  # calculate weights (NA: not penalized, 0: drop)
  if (missing(weights)) {
    weights <- rep(1, dim(x)[[2]])
  } else {
    if (!length(weights) == dim(x)[[2]]) stop("Weights vector has not length p")

    nonpen <- which(is.na(weights))
    if (any(weights[-nonpen] < 0)) stop("Weights must be positive")
    remove <- weights == 0
    remove[nonpen] <- FALSE

    x <- x[, !remove, drop = FALSE]

    rem <- which(weights == 0)
    z <- z[, !(ranInd %in% rem), drop = FALSE]

    weights <- weights[!remove]
    nonpen <- which(is.na(weights))

    ran1 <- logical(dim(x)[[2]])
    ran1[ranInd] <- TRUE
    ran2 <- ran1[!remove]
    ranInd <- which(ran2)
  }

  if (standardize) {
    xOr <- x
    meanx <- apply(x[, -1], 2, mean)
    sdx <- apply(x[, -1], 2, sd)
    x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
  }

  # crucial allocations
  grp <- factor(grp)
  N <- length(levels(grp)) # N is the number of groups
  p <- dim(x)[[2]] # p is the number of covariates
  q <- dim(z)[[2]] # q is the number of random effects variables
  ntot <- length(y) # ntot is the total number of observations
  Q <- q * (q + 1) / 2 # maximum number of variance components parameters

  # turn list to dataframe
  ylist.syn <- my_syndata$Y_list_syn
  xlist.syn <- my_syndata$X_list_syn
  zlist.syn <- my_syndata$Z_list_syn
  n_i <- nrow(xlist.syn[[1]])
  mydata.syn <- matrix(0, nc = p + q + 2, nr = N * n_i)
  for (i in 1:N) {
    mydata.syn[((i - 1) * n_i + 1):(i * n_i), 1] <- ylist.syn[[i]]
    mydata.syn[((i - 1) * n_i + 1):(i * n_i), 2:(p + 1)] <- xlist.syn[[i]]
    mydata.syn[((i - 1) * n_i + 1):(i * n_i), (p + 2):(p + 1 + q)] <- zlist.syn[[i]]
    mydata.syn[((i - 1) * n_i + 1):(i * n_i), (p + q + 2)] <- i
  }

  y.syn <- mydata.syn[, 1]
  x.syn <- mydata.syn[, 2:(p + 1)]
  z.syn <- mydata.syn[, (p + 2):(p + 1 + q)]
  grp.syn <- factor(mydata.syn[, p + q + 2])

  # save the grouped information as components of a list
  yGrp <- split(y, grp)
  xGrp <- split.data.frame(x, grp)
  zGrp <- split.data.frame(z, grp)
  # print(class(z))
  zIdGrp <- mapply(ZIdentity, zGrp)
  # cat("dim of z.syn is ",dim(z.syn),"\n")
  # print(grp.syn)
  yGrp.syn <- split(y.syn, grp.syn)
  xGrp.syn <- split.data.frame(x.syn, grp.syn)
  zGrp.syn <- split.data.frame(z.syn, grp.syn)
  # print(class(z.syn))
  # print(dim(zGrp.syn[[1]]))
  zIdGrp.syn <- mapply(ZIdentity, zGrp.syn)

  M.syn <- nrow(xGrp.syn[[1]])
  tau_syn <- tau_syn * M.syn
  # --- Calculating or adopting the starting values ---
  # ---------------------------------------------------

  if (missing(coefInit)) {
    # --- starting value for beta ---
    set.seed(control$seed)
    if (startValue == 1) # 10-fold cv Lasso
      {
        init <- optL1(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else if (startValue == 2) # 10-fold cv Ridge Regression
      {
        init <- optL2(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else # Null start
    {
      betaStart <- c(mean(y), rep(0, p - 1))
    }

    # --- starting values for the variance components parameters ---
    covStart <- covStartingValues(xGrp, yGrp, zGrp, zIdGrp, betaStart, ntot, N)
    sStart <- covStart$sigma
    tStart <- covStart$tau

    if (pdMat == "pdSym") parsStart <- vecli(tStart * diag(q))
    if (pdMat == "pdDiag") parsStart <- rep(tStart, q)
    if (pdMat == "pdIdent") parsStart <- tStart
  } else {
    betaStart <- coefInit[[1]]
    parsStart <- coefInit[[2]]
    sStart <- coefInit[[3]]
  }

  if (pdMat == "pdSym") PsiStart <- crossprod(triang(parsStart, q))
  if (pdMat == "pdDiag") PsiStart <- diag(parsStart^2, nrow = length(parsStart))
  if (pdMat == "pdIdent") PsiStart <- parsStart^2 * diag(q)

  # --- Calculate objective function for the starting values ---
  # ------------------------------------------------------------
  # cat("dim of zGrp",dim(zGrp[[1]]),'\n')
  lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiStart, sigma = sStart))
  # cat("dim of zGrp.syn",dim(zGrp.syn[[1]]),'\n')
  lambdaInvGrp.syn <- mapply(LambdaInv, Z = zGrp.syn, ZId = zIdGrp.syn, MoreArgs = list(Psi = PsiStart, sigma = sStart))

  ll1 <- 1 / 2 * ntot * log(2 * pi)
  fctStart <- ObjFunction_syn(
    xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaStart, weights = weights, lambda = lambda, nonpen = nonpen, ll1 = ll1,
    xGroup.syn = xGrp.syn, yGroup.syn = yGrp.syn, LGroup.syn = lambdaInvGrp.syn
  )
  if (control$trace > 3) cat("fctStart:", fctStart, "\n")


  # --- Coordinate Gradient Descent-iteration ---
  # ---------------------------------------------

  # some necessary allocations:
  betaIter <- betaStart
  sIter <- sStart
  parsIter <- parsStart
  PsiIter <- PsiStart
  convPar <- crossprod(betaIter)
  convCov <- crossprod(c(sStart, parsStart))
  fctIter <- convFct <- fctStart
  hessian0 <- rep(0, p)
  mat0 <- matrix(0, ncol = p, nrow = N)
  # mat0.syn <- matrix(0,ncol = p,nrow = )
  stopped <- FALSE
  doAll <- FALSE
  converged <- 0
  counterIn <- 0
  counter <- 0 # counts the number of outer iterations

  while ((control$maxIter > counter) & (convPar > control$tol | convFct > control$tol | convCov > control$tol | !doAll)) {
    # print(sIter)
    if (control$maxIter == counter + 1) {
      if (control$trace > 2) cat("maxIter reached", "\n")
      converged <- converged + 1
    }

    # cat("PsiIter is ",dim(PsiIter),"\n")
    # arguments that change at each iteration
    counter <- counter + 1
    if (control$trace == 2) cat(counter, "...", "\n")
    betaIterOld <- betaIter
    fctIterOld <- fctIter
    covIterOld <- c(sIter, parsIter)


    # --- optimization w.r.t the fixed effects vector beta ---
    # --------------------------------------------------------

    activeSet <- 1:length(betaIter)
    if ((length(activeSet) > min(p, ntot)) & (lambda > 0) & (stopSat) & (counter > 2)) {
      stopped <- TRUE
      break
    }
    # here will have some problem if we are in high dimension-Li Weihao
    if (counterIn == 0 | counterIn > control$number) {
      doAll <- TRUE
      activeSet <- 1:p
      counterIn <- 1
    } else {
      doAll <- FALSE
      counterIn <- counterIn + 1
    }

    # calculate the hessian matrices for j in the activeSet
    HessIter <- HessIterTrunc <- HessianMatrix(xGroup = xGrp, LGroup = lambdaInvGrp, activeSet = activeSet, N = N, hessian = hessian0, mat = mat0[, activeSet, drop = FALSE])
    HessIter[activeSet] <- pmin(pmax(HessIter[activeSet], control$lower), control$upper)

    HessIter.syn <- HessIterTrunc.syn <- HessianMatrix(xGroup = xGrp.syn, LGroup = lambdaInvGrp.syn, activeSet = activeSet, N = N, hessian = hessian0, mat = mat0[, activeSet, drop = FALSE])
    HessIter.syn[activeSet] <- pmin(pmax(HessIter.syn[activeSet], control$lower), control$upper)
    # print(lambdaInvGrp.syn)
    # cat("current iteration is ",counter,'\n')

    LxGrp <- as1(xGrp, lambdaInvGrp, activeSet, N = N)
    LxGrp.syn <- as1(xGrp.syn, lambdaInvGrp.syn, activeSet, N = N)

    ll2 <- nlogdet(LGroup = lambdaInvGrp)

    for (j in activeSet)
    {
      cut1 <- as2(x = x, y = y, b = betaIter, j = j, activeSet = activeSet, group = grp, sGroup = LxGrp)
      cut1.syn <- as2(x = x.syn, y = y.syn, b = betaIter, j = j, activeSet = activeSet, group = grp.syn, sGroup = LxGrp.syn)
      JinNonpen <- j %in% nonpen


      if (JinNonpen) {
        betaIter[j] <- cut1 / HessIter[j]
      } else {
        betaIter[j] <- (cut1 + tau_syn / M.syn * cut1.syn) / (HessIter[j] + tau_syn / M.syn * HessIter.syn[j])
      }


      if (control$trace > 3) {
        fctIter <- ObjFunction_syn(
          xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, weights = weights,
          lambda = lambda, nonpen = nonpen, ll1 = ll1, ll2 = ll2, xGroup.syn = xGrp.syn, yGroup.syn = yGrp.syn, LGroup.syn = lambdaInvGrp.syn
        )
        print(fctIter)
      }
    }

    if (control$trace > 3) cat("-----------", "\n")

    # --- optimization w.r.t the variance components parameters ---
    # -------------------------------------------------------------

    # calculations before the covariance optimization
    activeSet <- 1:length(betaIter)
    resGrp <- ResAsSplit(x = x, y = y, b = betaIter, f = grp, activeSet = activeSet)
    ll4 <- lambda * sum(abs(betaIter[-nonpen]) / weights[-nonpen])

    # optimization of the free parameters in \Psi
    if (pdMat == "pdSym") {
      covParOpt <- MLpdSym(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt, CovInt = control$CovInt
      )
      PsiIter <- crossprod(triang(covParOpt$pars, q))
    } else if (pdMat == "pdDiag") {
      covParOpt <- MLpdDiag(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- diag(covParOpt$pars^2, nrow = length(covParOpt$pars))
    } else if (pdMat == "pdIdent") {
      covParOpt <- MLpdIdent(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- covParOpt$par^2 * diag(q)
    }

    parsIter <- covParOpt$pars
    fctIter <- covParOpt$fct
    if (control$trace > 3) cat("-----------", "\n")

    # optimization of the error variance \sigma^2
    covParOpt <- MLsigma(
      zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, q = q, ll1 = ll1, ll4 = ll4, true.sigma = sIter, Psi = PsiIter,
      trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
    )
    sIter <- covParOpt$sigma # ; if(sIter<0.0001) sIter=covIterOld[1];print(sIter)
    fctIter <- covParOpt$fct
    # cat("sigma is ",sIter,"\n")
    if (control$trace > 3) cat("-----------", "\n") # cat("obj.fct III:",fctIter,"\n")

    lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiIter, sigma = sIter))
    lambdaInvGrp.syn <- mapply(LambdaInv, Z = zGrp.syn, ZId = zIdGrp.syn, MoreArgs = list(Psi = PsiIter, sigma = sIter))

    covIter <- c(sIter, parsIter)
    print(covIter)
    # --- check convergence ---
    convPar <- sqrt(crossprod(betaIter - betaIterOld)) / (1 + sqrt(crossprod(betaIter)))
    convFct <- abs((fctIterOld - fctIter) / (1 + abs(fctIter)))
    convCov <- sqrt(crossprod(covIter - covIterOld)) / (1 + sqrt(crossprod(covIter)))

    if ((convPar <= control$tol) & (convFct <= control$tol) & (convCov <= control$tol)) counterIn <- 0
  }

  if (standardize) {
    betaIter[-1] <- betaIter[-1] / sdx
    betaIter[1] <- betaIter[1] - sum(meanx * betaIter[-1])
    x <- xOr
    xGrp <- split.data.frame(xOr, grp)
  }


  # --- prediction of the random effects ---
  # ----------------------------------------

  biGroup <- u <- list()
  length(biGroup) <- length(u) <- N
  Psi <- PsiIter
  # cat("psi is ",diag(Psi),"\n")
  corPsi <- cov2cor(Psi)

  if (pdMat == "pdSym") {
    cholPsi <- t(triang(parsIter, q))
  } else if (pdMat == "pdDiag") {
    cholPsi <- diag(parsIter, nrow = length(parsIter))
  } else if (pdMat == "pdIdent") cholPsi <- parsIter * diag(q)

  for (i in 1:N)
  {
    u[[i]] <- sIter * ginv(t(zGrp[[i]] %*% cholPsi) %*% zGrp[[i]] %*% cholPsi + sIter^2 * diag(q)) %*% t(zGrp[[i]] %*% cholPsi) %*% resGrp[[i]]
    biGroup[[i]] <- 1 / sIter * cholPsi %*% u[[i]]
  }


  #  --- final calculations ---
  # ---------------------------

  # fitted values and residuals
  residGrp <- fittedGrp <- list()
  length(residGrp) <- length(fittedGrp) <- N
  for (i in 1:N)
  {
    fittedGrp[[i]] <- xGrp[[i]][, activeSet, drop = FALSE] %*% betaIter[activeSet, drop = FALSE] + zGrp[[i]] %*% biGroup[[i]]
    residGrp[[i]] <- yGrp[[i]] - fittedGrp[[i]]
  }
  residuals <- unlist(residGrp)
  fitted <- unlist(fittedGrp)

  # random effects, sorted per subject
  u <- unlist(u) # corresponds to lmer@u
  bi <- unlist(biGroup) # unsorted random effects

  # random effects, sorted per effect
  ranef <- bi[order(rep(1:q, N))] # corresponds to lmer@ranef

  # fixed effects without names
  fixef <- betaIter
  names(fixef) <- NULL

  # --- summary information ---
  # ---------------------------
  npar <- sum(betaIter != 0) + length(c(sIter, parsIter))
  logLik <- MLloglik(xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, ntot = ntot, N = N, activeSet = which(betaIter != 0))
  deviance <- -2 * logLik
  aic <- -2 * logLik + 2 * npar
  bic <- -2 * logLik + log(ntot) * npar

  if (any(parsIter == 0)) cat("Redundant covariance parameters.", "\n")
  if (converged > 0) cat("Algorithm does not properly converge.", "\n")
  if (stopped) {
    cat("|activeSet|>=min(p,ntot): Increase lambda or set stopSat=FALSE.", "\n")
    sIter <- parsIter <- nlogLik <- aic <- bic <- NA
    betaIter <- rep(NA, p)
    bi <- fitted <- residuals <- NULL
  }

  out <- list(
    data = list(x = x, y = y, z = z, grp = grp), weights = weights, coefInit = list(betaStart = betaStart, parsStart = parsStart, sStart = sStart),
    lambda = lambda, sigma = abs(sIter), pars = parsIter, coefficients = betaIter, random = bi, u = u, ranef = ranef, fixef = fixef, fitted.values = fitted,
    residuals = residuals, Psi = Psi, corPsi = corPsi, converged = converged, logLik = logLik, npar = npar, deviance = deviance,
    aic = aic, bic = bic, nonpen = nonpen, counter = counter, pdMat = pdMat, method = method, CovOpt = CovOpt, control = control, call = match.call(),
    stopped = stopped, ranInd = ranInd, objective = fctIter
  )
  out
  structure(out, class = "lmmlasso")
}


# this function is designed for cv, give 4/5 data and 1/5 validation data, do prediction on test data and compute residual on validation set
lmmsyn.pred_testdata <- function(x, y, z = x, grp, weights = NULL, coefInit = NULL, lambda, startValue = 1, nonpen = 1, pdMat = c("pdIdent", "pdDiag", "pdSym"),
                                 method = "ML", CovOpt = c("nlminb", "optimize"), stopSat = TRUE, standardize = TRUE, control = lmmlassoControl(),
                                 ranInd = 1:dim(z)[[2]], my_syndata, tau_syn = 0.1, test_x, test_y, test_z, test_grp, ...) {
  # --- Introductory checks and reformulations ---
  # ----------------------------------------------

  CovOpt <- match.arg(CovOpt)
  pdMat <- match.arg(pdMat)

  # do some checks
  if (!is.matrix(x)) stop("x has to be a matrix")
  if (any(is.na(x))) stop("Missing values in x not allowed")
  if (any(is.na(y))) stop("Missing values in y not allowed")
  if (!is.matrix(z)) stop("z has to be a matrix")
  if (any(is.na(z))) stop("Missing values in z not allowed")
  if (!is.numeric(y)) stop("y has to be of type 'numeric'")
  if (nrow(x) != length(y)) stop("x and y have not correct dimensions")
  if (any(x[, 1] != rep(1, dim(x)[[1]]))) stop("first column is not the intercept")
  if (length(levels(grp)) == 1) stop("Only one group. No covariance parameters!")
  if (!all(nonpen %in% 1:dim(x)[[2]])) stop("Error with the argument nonpen")
  if (length(lambda) != 1) stop("lambda can only be one number")
  if (lambda < 0) stop("regularization parameter must be positive")

  # calculate weights (NA: not penalized, 0: drop)
  if (missing(weights)) {
    weights <- rep(1, dim(x)[[2]])
  } else {
    if (!length(weights) == dim(x)[[2]]) stop("Weights vector has not length p")

    nonpen <- which(is.na(weights))
    if (any(weights[-nonpen] < 0)) stop("Weights must be positive")
    remove <- weights == 0
    remove[nonpen] <- FALSE

    x <- x[, !remove, drop = FALSE]

    rem <- which(weights == 0)
    z <- z[, !(ranInd %in% rem), drop = FALSE]

    weights <- weights[!remove]
    nonpen <- which(is.na(weights))

    ran1 <- logical(dim(x)[[2]])
    ran1[ranInd] <- TRUE
    ran2 <- ran1[!remove]
    ranInd <- which(ran2)
  }

  if (standardize) {
    xOr <- x
    meanx <- apply(x[, -1], 2, mean)
    sdx <- apply(x[, -1], 2, sd)
    x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
  }

  # crucial allocations
  grp <- factor(grp)
  N <- length(levels(grp)) # N is the number of groups
  p <- dim(x)[[2]] # p is the number of covariates
  q <- dim(z)[[2]] # q is the number of random effects variables
  ntot <- length(y) # ntot is the total number of observations
  Q <- q * (q + 1) / 2 # maximum number of variance components parameters

  # turn list to dataframe
  ylist.syn <- my_syndata$Y_list_syn
  xlist.syn <- my_syndata$X_list_syn
  zlist.syn <- my_syndata$Z_list_syn
  n_i <- nrow(xlist.syn[[1]])
  mydata.syn <- matrix(0, nc = p + q + 2, nr = N * n_i)
  for (i in 1:N) {
    mydata.syn[((i - 1) * n_i + 1):(i * n_i), 1] <- ylist.syn[[i]]
    mydata.syn[((i - 1) * n_i + 1):(i * n_i), 2:(p + 1)] <- xlist.syn[[i]]
    mydata.syn[((i - 1) * n_i + 1):(i * n_i), (p + 2):(p + 1 + q)] <- zlist.syn[[i]]
    mydata.syn[((i - 1) * n_i + 1):(i * n_i), (p + q + 2)] <- i
  }

  y.syn <- mydata.syn[, 1]
  x.syn <- mydata.syn[, 2:(p + 1)]
  z.syn <- mydata.syn[, (p + 2):(p + 1 + q)]
  grp.syn <- factor(mydata.syn[, p + q + 2])

  # save the grouped information as components of a list
  yGrp <- split(y, grp)
  xGrp <- split.data.frame(x, grp)
  zGrp <- split.data.frame(z, grp)
  # print(class(z))
  zIdGrp <- mapply(ZIdentity, zGrp)
  # cat("dim of z.syn is ",dim(z.syn),"\n")
  # print(grp.syn)
  yGrp.syn <- split(y.syn, grp.syn)
  xGrp.syn <- split.data.frame(x.syn, grp.syn)
  zGrp.syn <- split.data.frame(z.syn, grp.syn)
  # print(class(z.syn))
  # print(dim(zGrp.syn[[1]]))
  zIdGrp.syn <- mapply(ZIdentity, zGrp.syn)

  M.syn <- nrow(xGrp.syn[[1]])
  tau_syn <- tau_syn * M.syn
  # --- Calculating or adopting the starting values ---
  # ---------------------------------------------------

  if (missing(coefInit)) {
    # --- starting value for beta ---
    set.seed(control$seed)
    if (startValue == 1) # 10-fold cv Lasso
      {
        init <- optL1(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else if (startValue == 2) # 10-fold cv Ridge Regression
      {
        init <- optL2(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else # Null start
    {
      betaStart <- c(mean(y), rep(0, p - 1))
    }

    # --- starting values for the variance components parameters ---
    covStart <- covStartingValues(xGrp, yGrp, zGrp, zIdGrp, betaStart, ntot, N)
    sStart <- covStart$sigma
    tStart <- covStart$tau

    if (pdMat == "pdSym") parsStart <- vecli(tStart * diag(q))
    if (pdMat == "pdDiag") parsStart <- rep(tStart, q)
    if (pdMat == "pdIdent") parsStart <- tStart
  } else {
    betaStart <- coefInit[[1]]
    parsStart <- coefInit[[2]]
    sStart <- coefInit[[3]]
  }

  if (pdMat == "pdSym") PsiStart <- crossprod(triang(parsStart, q))
  if (pdMat == "pdDiag") PsiStart <- diag(parsStart^2, nrow = length(parsStart))
  if (pdMat == "pdIdent") PsiStart <- parsStart^2 * diag(q)

  # --- Calculate objective function for the starting values ---
  # ------------------------------------------------------------
  # cat("dim of zGrp",dim(zGrp[[1]]),'\n')
  lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiStart, sigma = sStart))
  # cat("dim of zGrp.syn",dim(zGrp.syn[[1]]),'\n')
  lambdaInvGrp.syn <- mapply(LambdaInv, Z = zGrp.syn, ZId = zIdGrp.syn, MoreArgs = list(Psi = PsiStart, sigma = sStart))

  ll1 <- 1 / 2 * ntot * log(2 * pi)
  fctStart <- ObjFunction_syn(
    xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaStart, weights = weights, lambda = lambda, nonpen = nonpen, ll1 = ll1,
    xGroup.syn = xGrp.syn, yGroup.syn = yGrp.syn, LGroup.syn = lambdaInvGrp.syn
  )
  if (control$trace > 3) cat("fctStart:", fctStart, "\n")


  # --- Coordinate Gradient Descent-iteration ---
  # ---------------------------------------------

  # some necessary allocations:
  betaIter <- betaStart
  sIter <- sStart
  parsIter <- parsStart
  PsiIter <- PsiStart
  convPar <- crossprod(betaIter)
  convCov <- crossprod(c(sStart, parsStart))
  fctIter <- convFct <- fctStart
  hessian0 <- rep(0, p)
  mat0 <- matrix(0, ncol = p, nrow = N)
  # mat0.syn <- matrix(0,ncol = p,nrow = )
  stopped <- FALSE
  doAll <- FALSE
  converged <- 0
  counterIn <- 0
  counter <- 0 # counts the number of outer iterations

  while ((control$maxIter > counter) & (convPar > control$tol | convFct > control$tol | convCov > control$tol | !doAll)) {
    # print(sIter)
    if (control$maxIter == counter + 1) {
      if (control$trace > 2) cat("maxIter reached", "\n")
      converged <- converged + 1
    }
    # cat("PsiIter is ",dim(PsiIter),"\n")
    # arguments that change at each iteration
    counter <- counter + 1
    if (control$trace == 2) cat(counter, "...", "\n")
    betaIterOld <- betaIter
    fctIterOld <- fctIter
    covIterOld <- c(sIter, parsIter)

    # --- optimization w.r.t the fixed effects vector beta ---
    # --------------------------------------------------------

    activeSet <- 1:length(betaIter)
    if ((length(activeSet) > min(p, ntot)) & (lambda > 0) & (stopSat) & (counter > 2)) {
      stopped <- TRUE
      break
    }
    # here will have some problem if we are in high dimension-Li Weihao
    if (counterIn == 0 | counterIn > control$number) {
      doAll <- TRUE
      activeSet <- 1:p
      counterIn <- 1
    } else {
      doAll <- FALSE
      counterIn <- counterIn + 1
    }

    # calculate the hessian matrices for j in the activeSet
    HessIter <- HessIterTrunc <- HessianMatrix(xGroup = xGrp, LGroup = lambdaInvGrp, activeSet = activeSet, N = N, hessian = hessian0, mat = mat0[, activeSet, drop = FALSE])
    HessIter[activeSet] <- pmin(pmax(HessIter[activeSet], control$lower), control$upper)

    HessIter.syn <- HessIterTrunc.syn <- HessianMatrix(xGroup = xGrp.syn, LGroup = lambdaInvGrp.syn, activeSet = activeSet, N = N, hessian = hessian0, mat = mat0[, activeSet, drop = FALSE])
    HessIter.syn[activeSet] <- pmin(pmax(HessIter.syn[activeSet], control$lower), control$upper)



    LxGrp <- as1(xGrp, lambdaInvGrp, activeSet, N = N)
    LxGrp.syn <- as1(xGrp.syn, lambdaInvGrp.syn, activeSet, N = N)

    ll2 <- nlogdet(LGroup = lambdaInvGrp)

    for (j in activeSet)
    {
      cut1 <- as2(x = x, y = y, b = betaIter, j = j, activeSet = activeSet, group = grp, sGroup = LxGrp)
      cut1.syn <- as2(x = x.syn, y = y.syn, b = betaIter, j = j, activeSet = activeSet, group = grp.syn, sGroup = LxGrp.syn)
      JinNonpen <- j %in% nonpen


      if (JinNonpen) {
        betaIter[j] <- cut1 / HessIter[j]
      } else {
        betaIter[j] <- (cut1 + tau_syn / M.syn * cut1.syn) / (HessIter[j] + tau_syn / M.syn * HessIter.syn[j])
      }


      if (control$trace > 3) {
        fctIter <- ObjFunction_syn(
          xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, weights = weights,
          lambda = lambda, nonpen = nonpen, ll1 = ll1, ll2 = ll2, xGroup.syn = xGrp.syn, yGroup.syn = yGrp.syn, LGroup.syn = lambdaInvGrp.syn
        )
        print(fctIter)
      }
    }

    if (control$trace > 3) cat("-----------", "\n")

    # --- optimization w.r.t the variance components parameters ---
    # -------------------------------------------------------------

    # calculations before the covariance optimization
    activeSet <- 1:length(betaIter)
    resGrp <- ResAsSplit(x = x, y = y, b = betaIter, f = grp, activeSet = activeSet)
    ll4 <- lambda * sum(abs(betaIter[-nonpen]) / weights[-nonpen])

    # optimization of the free parameters in \Psi
    if (pdMat == "pdSym") {
      covParOpt <- MLpdSym(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt, CovInt = control$CovInt
      )
      PsiIter <- crossprod(triang(covParOpt$pars, q))
    } else if (pdMat == "pdDiag") {
      covParOpt <- MLpdDiag(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- diag(covParOpt$pars^2, nrow = length(covParOpt$pars))
    } else if (pdMat == "pdIdent") {
      covParOpt <- MLpdIdent(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- covParOpt$par^2 * diag(q)
    }

    parsIter <- covParOpt$pars
    fctIter <- covParOpt$fct
    if (control$trace > 3) cat("-----------", "\n")

    # optimization of the error variance \sigma^2
    covParOpt <- MLsigma(
      zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, q = q, ll1 = ll1, ll4 = ll4, true.sigma = sIter, Psi = PsiIter,
      trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
    )
    sIter <- covParOpt$sigma # ; if(sIter<0.0001) sIter=covIterOld[1];print(sIter)
    fctIter <- covParOpt$fct
    # cat("sigma is ",sIter,"\n")
    if (control$trace > 3) cat("-----------", "\n") # cat("obj.fct III:",fctIter,"\n")

    lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiIter, sigma = sIter))
    lambdaInvGrp.syn <- mapply(LambdaInv, Z = zGrp.syn, ZId = zIdGrp.syn, MoreArgs = list(Psi = PsiIter, sigma = sIter))

    covIter <- c(sIter, parsIter)
    print(covIter)
    # --- check convergence ---
    convPar <- sqrt(crossprod(betaIter - betaIterOld)) / (1 + sqrt(crossprod(betaIter)))
    convFct <- abs((fctIterOld - fctIter) / (1 + abs(fctIter)))
    convCov <- sqrt(crossprod(covIter - covIterOld)) / (1 + sqrt(crossprod(covIter)))

    if ((convPar <= control$tol) & (convFct <= control$tol) & (convCov <= control$tol)) counterIn <- 0
  }

  if (standardize) {
    betaIter[-1] <- betaIter[-1] / sdx
    betaIter[1] <- betaIter[1] - sum(meanx * betaIter[-1])
    x <- xOr
    xGrp <- split.data.frame(xOr, grp)
  }


  # --- prediction of the random effects ---
  # ----------------------------------------

  biGroup <- u <- list()
  length(biGroup) <- length(u) <- N
  Psi <- PsiIter
  # cat("psi is ",diag(Psi),"\n")
  corPsi <- cov2cor(Psi)

  if (pdMat == "pdSym") {
    cholPsi <- t(triang(parsIter, q))
  } else if (pdMat == "pdDiag") {
    cholPsi <- diag(parsIter, nrow = length(parsIter))
  } else if (pdMat == "pdIdent") cholPsi <- parsIter * diag(q)

  for (i in 1:N)
  {
    u[[i]] <- sIter * ginv(t(zGrp[[i]] %*% cholPsi) %*% zGrp[[i]] %*% cholPsi + sIter^2 * diag(q)) %*% t(zGrp[[i]] %*% cholPsi) %*% resGrp[[i]]
    biGroup[[i]] <- 1 / sIter * cholPsi %*% u[[i]]
  }


  #  --- final calculations ---
  # ---------------------------
  yGrp_test <- split(test_y, test_grp)
  xGrp_test <- split.data.frame(test_x, test_grp)
  zGrp_test <- split.data.frame(test_z, test_grp)


  # fitted values and residuals
  residGrp <- fittedGrp <- list()
  length(residGrp) <- length(fittedGrp) <- N
  for (i in 1:N)
  {
    fittedGrp[[i]] <- xGrp_test[[i]][, activeSet, drop = FALSE] %*% betaIter[activeSet, drop = FALSE] + zGrp_test[[i]] %*% biGroup[[i]]
    residGrp[[i]] <- yGrp_test[[i]] - fittedGrp[[i]]
  }
  residuals <- unlist(residGrp)
  fitted <- unlist(fittedGrp)

  # random effects, sorted per subject
  u <- unlist(u) # corresponds to lmer@u
  bi <- unlist(biGroup) # unsorted random effects

  # random effects, sorted per effect
  ranef <- bi[order(rep(1:q, N))] # corresponds to lmer@ranef

  # fixed effects without names
  fixef <- betaIter
  names(fixef) <- NULL

  # --- summary information ---
  # ---------------------------
  npar <- sum(betaIter != 0) + length(c(sIter, parsIter))
  logLik <- MLloglik(xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, ntot = ntot, N = N, activeSet = which(betaIter != 0))
  deviance <- -2 * logLik
  aic <- -2 * logLik + 2 * npar
  bic <- -2 * logLik + log(ntot) * npar

  if (any(parsIter == 0)) cat("Redundant covariance parameters.", "\n")
  if (converged > 0) cat("Algorithm does not properly converge.", "\n")
  if (stopped) {
    cat("|activeSet|>=min(p,ntot): Increase lambda or set stopSat=FALSE.", "\n")
    sIter <- parsIter <- nlogLik <- aic <- bic <- NA
    betaIter <- rep(NA, p)
    bi <- fitted <- residuals <- NULL
  }

  out <- list(
    data = list(x = x, y = y, z = z, grp = grp), weights = weights, coefInit = list(betaStart = betaStart, parsStart = parsStart, sStart = sStart),
    lambda = lambda, sigma = abs(sIter), pars = parsIter, coefficients = betaIter, random = bi, u = u, ranef = ranef, fixef = fixef, fitted.values = fitted,
    residuals = residuals, Psi = Psi, corPsi = corPsi, converged = converged, logLik = logLik, npar = npar, deviance = deviance,
    aic = aic, bic = bic, nonpen = nonpen, counter = counter, pdMat = pdMat, method = method, CovOpt = CovOpt, control = control, call = match.call(),
    stopped = stopped, ranInd = ranInd, objective = fctIter
  )
  out
  structure(out, class = "lmmlasso")
}






lmmsyn.learningrate <- function(x, y, z = x, grp, weights = NULL, coefInit = NULL, lambda, startValue = 1, nonpen = 1, pdMat = c("pdIdent", "pdDiag", "pdSym"),
                                method = "ML", CovOpt = c("nlminb", "optimize"), stopSat = TRUE, standardize = TRUE, control = lmmlassoControl(),
                                ranInd = 1:dim(z)[[2]], my_syndata, tau_syn = 0.1, learningrate = 0.05, ...) {
  # --- Introductory checks and reformulations ---
  # ----------------------------------------------

  CovOpt <- match.arg(CovOpt)
  pdMat <- match.arg(pdMat)

  # do some checks
  if (!is.matrix(x)) stop("x has to be a matrix")
  if (any(is.na(x))) stop("Missing values in x not allowed")
  if (any(is.na(y))) stop("Missing values in y not allowed")
  if (!is.matrix(z)) stop("z has to be a matrix")
  if (any(is.na(z))) stop("Missing values in z not allowed")
  if (!is.numeric(y)) stop("y has to be of type 'numeric'")
  if (nrow(x) != length(y)) stop("x and y have not correct dimensions")
  if (any(x[, 1] != rep(1, dim(x)[[1]]))) stop("first column is not the intercept")
  if (length(levels(grp)) == 1) stop("Only one group. No covariance parameters!")
  if (!all(nonpen %in% 1:dim(x)[[2]])) stop("Error with the argument nonpen")
  if (length(lambda) != 1) stop("lambda can only be one number")
  if (lambda < 0) stop("regularization parameter must be positive")

  # calculate weights (NA: not penalized, 0: drop)
  if (missing(weights)) {
    weights <- rep(1, dim(x)[[2]])
  } else {
    if (!length(weights) == dim(x)[[2]]) stop("Weights vector has not length p")

    nonpen <- which(is.na(weights))
    if (any(weights[-nonpen] < 0)) stop("Weights must be positive")
    remove <- weights == 0
    remove[nonpen] <- FALSE

    x <- x[, !remove, drop = FALSE]

    rem <- which(weights == 0)
    z <- z[, !(ranInd %in% rem), drop = FALSE]

    weights <- weights[!remove]
    nonpen <- which(is.na(weights))

    ran1 <- logical(dim(x)[[2]])
    ran1[ranInd] <- TRUE
    ran2 <- ran1[!remove]
    ranInd <- which(ran2)
  }

  if (standardize) {
    xOr <- x
    meanx <- apply(x[, -1], 2, mean)
    sdx <- apply(x[, -1], 2, sd)
    x <- cbind(1, scale(x[, -1], center = meanx, scale = sdx))
  }

  # crucial allocations
  grp <- factor(grp)
  N <- length(levels(grp)) # N is the number of groups
  p <- dim(x)[[2]] # p is the number of covariates
  q <- dim(z)[[2]] # q is the number of random effects variables
  ntot <- length(y) # ntot is the total number of observations
  Q <- q * (q + 1) / 2 # maximum number of variance components parameters

  # turn list to dataframe
  ylist.syn <- my_syndata$Y_list_syn
  xlist.syn <- my_syndata$X_list_syn
  zlist.syn <- my_syndata$Z_list_syn
  n_i <- nrow(xlist.syn[[1]])
  mydata.syn <- matrix(0, nc = p + q + 2, nr = N * n_i)
  for (i in 1:N) {
    mydata.syn[((i - 1) * n_i + 1):(i * n_i), 1] <- ylist.syn[[i]]
    mydata.syn[((i - 1) * n_i + 1):(i * n_i), 2:(p + 1)] <- xlist.syn[[i]]
    mydata.syn[((i - 1) * n_i + 1):(i * n_i), (p + 2):(p + 1 + q)] <- zlist.syn[[i]]
    mydata.syn[((i - 1) * n_i + 1):(i * n_i), (p + q + 2)] <- i
  }

  y.syn <- mydata.syn[, 1]
  x.syn <- mydata.syn[, 2:(p + 1)]
  z.syn <- mydata.syn[, (p + 2):(p + 1 + q)]
  grp.syn <- factor(mydata.syn[, p + q + 2])

  # save the grouped information as components of a list
  yGrp <- split(y, grp)
  xGrp <- split.data.frame(x, grp)
  zGrp <- split.data.frame(z, grp)
  # print(class(z))
  zIdGrp <- mapply(ZIdentity, zGrp)
  # cat("dim of z.syn is ",dim(z.syn),"\n")
  # print(grp.syn)
  yGrp.syn <- split(y.syn, grp.syn)
  xGrp.syn <- split.data.frame(x.syn, grp.syn)
  zGrp.syn <- split.data.frame(z.syn, grp.syn)
  # print(class(z.syn))
  # print(dim(zGrp.syn[[1]]))
  zIdGrp.syn <- mapply(ZIdentity, zGrp.syn)

  M.syn <- nrow(xGrp.syn[[1]])
  tau_syn <- tau_syn * M.syn
  # --- Calculating or adopting the starting values ---
  # ---------------------------------------------------

  if (missing(coefInit)) {
    # --- starting value for beta ---
    set.seed(control$seed)
    if (startValue == 1) # 10-fold cv Lasso
      {
        init <- optL1(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else if (startValue == 2) # 10-fold cv Ridge Regression
      {
        init <- optL2(y, x[, -1], model = "linear", fold = 5, trace = FALSE)
        betaStart <- c(init$fullfit@unpenalized, init$fullfit@penalized)
      } else # Null start
    {
      betaStart <- c(mean(y), rep(0, p - 1))
    }

    # --- starting values for the variance components parameters ---
    covStart <- covStartingValues(xGrp, yGrp, zGrp, zIdGrp, betaStart, ntot, N)
    sStart <- covStart$sigma
    tStart <- covStart$tau

    if (pdMat == "pdSym") parsStart <- vecli(tStart * diag(q))
    if (pdMat == "pdDiag") parsStart <- rep(tStart, q)
    if (pdMat == "pdIdent") parsStart <- tStart
  } else {
    betaStart <- coefInit[[1]]
    parsStart <- coefInit[[2]]
    sStart <- coefInit[[3]]
  }

  if (pdMat == "pdSym") PsiStart <- crossprod(triang(parsStart, q))
  if (pdMat == "pdDiag") PsiStart <- diag(parsStart^2, nrow = length(parsStart))
  if (pdMat == "pdIdent") PsiStart <- parsStart^2 * diag(q)

  # --- Calculate objective function for the starting values ---
  # ------------------------------------------------------------
  # cat("dim of zGrp",dim(zGrp[[1]]),'\n')
  lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiStart, sigma = sStart))
  # cat("dim of zGrp.syn",dim(zGrp.syn[[1]]),'\n')
  lambdaInvGrp.syn <- mapply(LambdaInv, Z = zGrp.syn, ZId = zIdGrp.syn, MoreArgs = list(Psi = PsiStart, sigma = sStart))

  ll1 <- 1 / 2 * ntot * log(2 * pi)
  fctStart <- ObjFunction_syn(
    xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaStart, weights = weights, lambda = lambda, nonpen = nonpen, ll1 = ll1,
    xGroup.syn = xGrp.syn, yGroup.syn = yGrp.syn, LGroup.syn = lambdaInvGrp.syn
  )
  if (control$trace > 3) cat("fctStart:", fctStart, "\n")


  # --- Coordinate Gradient Descent-iteration ---
  # ---------------------------------------------

  # some necessary allocations:
  betaIter <- betaStart
  sIter <- sStart
  parsIter <- parsStart
  PsiIter <- PsiStart
  convPar <- crossprod(betaIter)
  convCov <- crossprod(c(sStart, parsStart))
  fctIter <- convFct <- fctStart
  hessian0 <- rep(0, p)
  mat0 <- matrix(0, ncol = p, nrow = N)
  # mat0.syn <- matrix(0,ncol = p,nrow = )
  stopped <- FALSE
  doAll <- FALSE
  converged <- 0
  counterIn <- 0
  counter <- 0 # counts the number of outer iterations
  counter_learningrate <- 0
  while ((control$maxIter > counter) & (convPar > control$tol | convFct > control$tol | convCov > control$tol | !doAll)) {
    counter_learningrate <- counter_learningrate + 1
    # print(sIter)
    if (control$maxIter == counter + 1) {
      if (control$trace > 2) cat("maxIter reached", "\n")
      converged <- converged + 1
    }

    # cat("PsiIter is ",dim(PsiIter),"\n")
    # arguments that change at each iteration
    counter <- counter + 1
    if (control$trace == 2) cat(counter, "...", "\n")
    betaIterOld <- betaIter
    fctIterOld <- fctIter
    covIterOld <- c(sIter, parsIter)


    # --- optimization w.r.t the fixed effects vector beta ---
    # --------------------------------------------------------

    activeSet <- 1:length(betaIter)
    if ((length(activeSet) > min(p, ntot)) & (lambda > 0) & (stopSat) & (counter > 2)) {
      stopped <- TRUE
      break
    }
    # here will have some problem if we are in high dimension-Li Weihao
    if (counterIn == 0 | counterIn > control$number) {
      doAll <- TRUE
      activeSet <- 1:p
      counterIn <- 1
    } else {
      doAll <- FALSE
      counterIn <- counterIn + 1
    }

    # calculate the hessian matrices for j in the activeSet
    HessIter <- HessIterTrunc <- HessianMatrix(xGroup = xGrp, LGroup = lambdaInvGrp, activeSet = activeSet, N = N, hessian = hessian0, mat = mat0[, activeSet, drop = FALSE])
    HessIter[activeSet] <- pmin(pmax(HessIter[activeSet], control$lower), control$upper)

    HessIter.syn <- HessIterTrunc.syn <- HessianMatrix(xGroup = xGrp.syn, LGroup = lambdaInvGrp.syn, activeSet = activeSet, N = N, hessian = hessian0, mat = mat0[, activeSet, drop = FALSE])
    HessIter.syn[activeSet] <- pmin(pmax(HessIter.syn[activeSet], control$lower), control$upper)



    LxGrp <- as1(xGrp, lambdaInvGrp, activeSet, N = N)
    LxGrp.syn <- as1(xGrp.syn, lambdaInvGrp.syn, activeSet, N = N)

    ll2 <- nlogdet(LGroup = lambdaInvGrp)

    for (j in activeSet)
    {
      cut1 <- as2(x = x, y = y, b = betaIter, j = j, activeSet = activeSet, group = grp, sGroup = LxGrp)
      cut1.syn <- as2(x = x.syn, y = y.syn, b = betaIter, j = j, activeSet = activeSet, group = grp.syn, sGroup = LxGrp.syn)
      JinNonpen <- j %in% nonpen


      if (JinNonpen) {
        betaIter[j] <- cut1 / HessIter[j]
      } else {
        betaIter[j] <- betaIterOld[j] - learningrate * 2 * (HessIter[j] - cut1 + tau_syn / M.syn * HessIter.syn[j] - tau_syn / M.syn * cut1.syn)
      }


      if (control$trace > 3) {
        fctIter <- ObjFunction_syn(
          xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, weights = weights,
          lambda = lambda, nonpen = nonpen, ll1 = ll1, ll2 = ll2, xGroup.syn = xGrp.syn, yGroup.syn = yGrp.syn, LGroup.syn = lambdaInvGrp.syn
        )
        print(fctIter)
      }
    }

    if (control$trace > 3) cat("-----------", "\n")

    # --- optimization w.r.t the variance components parameters ---
    # -------------------------------------------------------------

    # calculations before the covariance optimization
    activeSet <- 1:length(betaIter)
    resGrp <- ResAsSplit(x = x, y = y, b = betaIter, f = grp, activeSet = activeSet)
    ll4 <- lambda * sum(abs(betaIter[-nonpen]) / weights[-nonpen])

    # optimization of the free parameters in \Psi
    if (pdMat == "pdSym") {
      covParOpt <- MLpdSym(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt, CovInt = control$CovInt
      )
      PsiIter <- crossprod(triang(covParOpt$pars, q))
    } else if (pdMat == "pdDiag") {
      covParOpt <- MLpdDiag(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- diag(covParOpt$pars^2, nrow = length(covParOpt$pars))
    } else if (pdMat == "pdIdent") {
      covParOpt <- MLpdIdent(
        zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, sigma = sIter, pars = parsIter, q = q,
        thres = control$thres, ll1 = ll1, ll4 = ll4, trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
      )
      PsiIter <- covParOpt$par^2 * diag(q)
    }

    parsIter <- covParOpt$pars
    fctIter <- covParOpt$fct
    if (control$trace > 3) cat("-----------", "\n")

    # optimization of the error variance \sigma^2
    covParOpt <- MLsigma(
      zGroup = zGrp, zIdGroup = zIdGrp, resGroup = resGrp, q = q, ll1 = ll1, ll4 = ll4, true.sigma = sIter, Psi = PsiIter,
      trace = control$trace, CovOpt = CovOpt, VarInt = control$VarInt
    )
    sIter <- covParOpt$sigma # ; if(sIter<0.0001) sIter=covIterOld[1];print(sIter)
    fctIter <- covParOpt$fct
    # cat("sigma is ",sIter,"\n")
    if (control$trace > 3) cat("-----------", "\n") # cat("obj.fct III:",fctIter,"\n")

    lambdaInvGrp <- mapply(LambdaInv, Z = zGrp, ZId = zIdGrp, MoreArgs = list(Psi = PsiIter, sigma = sIter))
    lambdaInvGrp.syn <- mapply(LambdaInv, Z = zGrp.syn, ZId = zIdGrp.syn, MoreArgs = list(Psi = PsiIter, sigma = sIter))

    covIter <- c(sIter, parsIter)

    # --- check convergence ---
    convPar <- sqrt(crossprod(betaIter - betaIterOld)) / (1 + sqrt(crossprod(betaIter)))
    convFct <- abs((fctIterOld - fctIter) / (1 + abs(fctIter)))
    convCov <- sqrt(crossprod(covIter - covIterOld)) / (1 + sqrt(crossprod(covIter)))

    if ((convPar <= control$tol) & (convFct <= control$tol) & (convCov <= control$tol)) counterIn <- 0
    if (counter_learningrate > 200) {
      counter_learningrate <- 1
      learningrate <- learningrate / 2
    }
    cat("current iteration is ", counter, "difference with last iteration is ", crossprod(betaIter - betaIterOld), "\n")
  }

  if (standardize) {
    betaIter[-1] <- betaIter[-1] / sdx
    betaIter[1] <- betaIter[1] - sum(meanx * betaIter[-1])
    x <- xOr
    xGrp <- split.data.frame(xOr, grp)
  }


  # --- prediction of the random effects ---
  # ----------------------------------------

  biGroup <- u <- list()
  length(biGroup) <- length(u) <- N
  Psi <- PsiIter
  # cat("psi is ",diag(Psi),"\n")
  corPsi <- cov2cor(Psi)

  if (pdMat == "pdSym") {
    cholPsi <- t(triang(parsIter, q))
  } else if (pdMat == "pdDiag") {
    cholPsi <- diag(parsIter, nrow = length(parsIter))
  } else if (pdMat == "pdIdent") cholPsi <- parsIter * diag(q)

  for (i in 1:N)
  {
    u[[i]] <- sIter * ginv(t(zGrp[[i]] %*% cholPsi) %*% zGrp[[i]] %*% cholPsi + sIter^2 * diag(q)) %*% t(zGrp[[i]] %*% cholPsi) %*% resGrp[[i]]
    biGroup[[i]] <- 1 / sIter * cholPsi %*% u[[i]]
  }


  #  --- final calculations ---
  # ---------------------------

  # fitted values and residuals
  residGrp <- fittedGrp <- list()
  length(residGrp) <- length(fittedGrp) <- N
  for (i in 1:N)
  {
    fittedGrp[[i]] <- xGrp[[i]][, activeSet, drop = FALSE] %*% betaIter[activeSet, drop = FALSE] + zGrp[[i]] %*% biGroup[[i]]
    residGrp[[i]] <- yGrp[[i]] - fittedGrp[[i]]
  }
  residuals <- unlist(residGrp)
  fitted <- unlist(fittedGrp)

  # random effects, sorted per subject
  u <- unlist(u) # corresponds to lmer@u
  bi <- unlist(biGroup) # unsorted random effects

  # random effects, sorted per effect
  ranef <- bi[order(rep(1:q, N))] # corresponds to lmer@ranef

  # fixed effects without names
  fixef <- betaIter
  names(fixef) <- NULL

  # --- summary information ---
  # ---------------------------
  npar <- sum(betaIter != 0) + length(c(sIter, parsIter))
  logLik <- MLloglik(xGroup = xGrp, yGroup = yGrp, LGroup = lambdaInvGrp, b = betaIter, ntot = ntot, N = N, activeSet = which(betaIter != 0))
  deviance <- -2 * logLik
  aic <- -2 * logLik + 2 * npar
  bic <- -2 * logLik + log(ntot) * npar

  if (any(parsIter == 0)) cat("Redundant covariance parameters.", "\n")
  if (converged > 0) cat("Algorithm does not properly converge.", "\n")
  if (stopped) {
    cat("|activeSet|>=min(p,ntot): Increase lambda or set stopSat=FALSE.", "\n")
    sIter <- parsIter <- nlogLik <- aic <- bic <- NA
    betaIter <- rep(NA, p)
    bi <- fitted <- residuals <- NULL
  }

  out <- list(
    data = list(x = x, y = y, z = z, grp = grp), weights = weights, coefInit = list(betaStart = betaStart, parsStart = parsStart, sStart = sStart),
    lambda = lambda, sigma = abs(sIter), pars = parsIter, coefficients = betaIter, random = bi, u = u, ranef = ranef, fixef = fixef, fitted.values = fitted,
    residuals = residuals, Psi = Psi, corPsi = corPsi, converged = converged, logLik = logLik, npar = npar, deviance = deviance,
    aic = aic, bic = bic, nonpen = nonpen, counter = counter, pdMat = pdMat, method = method, CovOpt = CovOpt, control = control, call = match.call(),
    stopped = stopped, ranInd = ranInd, objective = fctIter
  )
  out
  structure(out, class = "lmmlasso")
}
