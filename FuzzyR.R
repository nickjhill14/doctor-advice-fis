fuzzyr.accuracy <- function(f, y, f.ref=0, scale.mase=NULL) {

    e <- y - f

    mae <- mean(abs(e), na.rm=TRUE)

    mse <- mean(e^2, na.rm=TRUE)
    rmse <- sqrt(mse)

    mase <- NULL
    if(!is.null(scale.mase)) {
        mase <- mean(abs(e/scale.mase), na.rm=TRUE)
    }

    e.ref <- y - f.ref
    mrae <- mean(abs(e/e.ref), na.rm=TRUE)  
    gmrae <- exp(mean(log(abs(e/e.ref)), na.rm=TRUE))

    pe <- e / y * 100
    mape <- mean(abs(pe), na.rm=TRUE)

    spe <- e / (abs(y) + abs(f)) * 200
    smape <- mean(abs(spe), na.rm=TRUE)

    out <- c(mae, rmse, mase, mrae, gmrae, mape, smape)

    if(is.null(scale.mase)) {
        names(out) <- c("MAE","RMSE", "MRAE", "GMRAE", "MAPE","sMAPE")
    } else {
        names(out) <- c("MAE","RMSE","MASE", "MRAE", "GMRAE", "MAPE","sMAPE")
    }

    out
}



#unlockBinding("anfis.optimise", e.ns)
#unlockBinding("anfis.optimise", e.pkg)

#e.pkg$anfis.optimise <- e.ns$

anfis.optimise <- function(anfis, data.trn, data.chk=NULL, epoch.total=100, stepsize=0.1, rate.inc=1.1, rate.dec=0.9, method=c("gradient", "lse"), err.log=F, online=0, lambda=1, opt.by="err.opt") {

    LI <- 1
    L1 <- 1 + 1
    L2 <- 1 + 2
    L3 <- 1 + 3
    L4 <- 1 + 4
    L5 <- 1 + 5

    if(method[2] == "lse") {
        for(i in 1:length(anfis$layer[[L4]]$node)) {
            if(anfis$layer[[L4]]$node[[i]]$mf.type != "linearmf") {
                stop("currently, the LSE method can only be applied to 'linearmf' for nodes in Layer 4")
            }
        }
    }

    if(length(online) != 1 || !is.element(online, c(0,1,2))) {
        stop("parameter online: 0 -- batch; 1 -- online; 2 -- semi-online")
    }

    input.num <- length(anfis$layer[[LI]]$node)
    output.num <- length(anfis$layer[[L5]]$node)

    data.num <- nrow(data.trn)
    input.stack <- input.stack.all <- data.trn[, 1:input.num]
    target <- target.all <- as.matrix(data.trn[, -(1:input.num)])
    scale.mase <- mean(abs(data.trn[,input.num] - data.trn[,(input.num+1)]))

    anfis.optimum <- anfis
    err.min <- Inf
    err.min.rec <- NULL
    err.min.rec.tmp <- NULL
    err.dec <- Inf
    epoch.last <- 1

    if(err.log) {
        err.opt <<- matrix(0, nrow=epoch.total, ncol=7)
        err.trn <<- matrix(0, nrow=epoch.total, ncol=7)
        err.chk <<- matrix(0, nrow=epoch.total, ncol=7)
        theta.L1 <<- NULL
        theta.L4 <<- NULL
        colnames(err.opt) <<- colnames(err.trn) <<- colnames(err.chk) <<- c("MAE","RMSE","MASE", "MRAE", "GMRAE", "MAPE","sMAPE")
    }

    epoch <- p <- 1
    if(online != 1) {
        p <- data.num
    } else {
        anfis.output.tmp <- target
        output.L3.tmp <- NULL
    }

    while(epoch <= epoch.total) {
        if(online == 1) {
            input.stack <- matrix(input.stack.all[p,], nrow=1)
            target <- matrix(target.all[p,], nrow=1)
        }

        output.LI <- anfis.LI.eval(anfis, input.stack)
        output.L1 <- anfis.L1.eval(anfis, output.LI, input.stack)
        output.L2 <- anfis.L2.eval(anfis, output.L1)
        output.L4.mf <- anfis.L4.mf.eval(anfis, input.stack)
        output.L2.which <- anfis.L2.which(anfis, output.L2, output.L4.mf)
        output.L3 <- anfis.L3.eval(anfis, output.L2, output.L2.which)

        if(method[2] == "lse") {
            anfis <- anfis.optimise.lse(anfis, output.L3, input.stack, target, online, lambda)
            output.L4.mf <- anfis.L4.mf.eval(anfis, input.stack)
#            if(ncol(output.L3[[1]]) != 1) {
#                output.L2.which <- anfis.L2.which(anfis, output.L2, output.L4.mf)
#                output.L3 <- anfis.L3.eval(anfis, output.L2, output.L2.which)
#            }
        }

        output.L4 <- anfis.L4.eval(output.L3, output.L4.mf)
        output.L5 <- anfis.L5.eval(output.L4)

        if(online == 1) {
            ## option 2
            #anfis.output.tmp[p] <- output.L5

            ## option 3
            output.L3.tmp <- lapply(1:length(output.L3), function(i) rbind(output.L3.tmp[[i]], output.L3[[i]]))
        }

        if(p == data.num) {
            if(online == 1) {
                ## option 1
                #anfis.output <- anfis.eval(anfis, input.stack.all)

                ## option 2
                #anfis.output <- anfis.output.tmp

                ## option 3
                output.L4.mf.tmp <- anfis.L4.mf.eval(anfis, input.stack.all)
                output.L4.tmp <- anfis.L4.eval(output.L3.tmp, output.L4.mf.tmp)
                output.L5.tmp <- anfis.L5.eval(output.L4.tmp)
                anfis.output <- output.L5.tmp
                output.L3.tmp <- NULL
            } else {
                anfis.output <- output.L5
            }

            if(err.log) {
                err.tmp <- fuzzyr.accuracy(anfis.output, data.trn[,(input.num+1)], data.trn[,input.num], scale.mase)
                err.opt[epoch,] <<- err.tmp
                
                err.tmp <- fuzzyr.accuracy(anfis.eval(anfis, data.trn[,1:input.num]), data.trn[,(input.num+1)], data.trn[,input.num], scale.mase)
                err.trn[epoch,] <<- err.tmp

                if(!is.null(data.chk)) {
                    err.tmp <- fuzzyr.accuracy(anfis.eval(anfis, data.chk[,1:input.num]), data.chk[,(input.num+1)], data.chk[,input.num], scale.mase)
                    err.chk[epoch,] <<- err.tmp
                }

                theta.L1.tmp <- NULL
                for( i in 1:length(anfis$layer[[L1]]$node) ) {
                    theta.L1.tmp <- c(theta.L1.tmp, anfis$layer[[L1]]$node[[i]]$mf.params)
                }

                theta.L4.tmp <- NULL
                for( i in 1:length(anfis$layer[[L4]]$node) ) {
                    theta.L4.tmp <- c(theta.L4.tmp, anfis$layer[[L4]]$node[[i]]$mf.params)
                }

                theta.L1 <<- rbind(theta.L1, theta.L1.tmp)
                theta.L4 <<- rbind(theta.L4, theta.L4.tmp)
            }

            err.opt.rmse <- sqrt(mean((anfis.output - target.all)^2))
            if(opt.by == "err.trn") {
                if(err.log) {
                    err.epoch <- err.trn[epoch,2]
                } else {
                    err.tmp <- fuzzyr.accuracy(anfis.eval(anfis, data.trn[,1:input.num]), data.trn[,(input.num+1)], data.trn[,input.num], scale.mase)
                    err.epoch <- err.tmp[2]
                }
            } else if (opt.by == "err.chk") {
                if(err.log) {
                    err.epoch <- err.chk[epoch,2]
                } else {
                    err.tmp <- fuzzyr.accuracy(anfis.eval(anfis, data.chk[,1:input.num]), data.chk[,(input.num+1)], data.chk[,input.num], scale.mase)
                    err.epoch <- err.tmp[2]                
                }
            } else {
                err.epoch <- err.opt.rmse
            }

            #cat(" epoch[", epoch, "], err=", err.epoch, "\n")
            if(err.epoch < err.min) {
                err.min <- err.epoch
                err.min.rec.tmp <- c(err.min.rec.tmp, err.min)
                if(length(err.min.rec.tmp) == 1000) {
                    err.min.rec <- c(err.min.rec, err.min.rec.tmp)
                    err.min.rec.tmp <- NULL
                }
                anfis.optimum <- anfis
            }

            ## to update step_size

            #cat("epoch - epoch.last: ", epoch - epoch.last, "\n")
            #cat("stepsize: ", stepsize, "\n")

            err.dec <- c(err.dec, err.opt.rmse)
            if( epoch - epoch.last >= 4) {
                err.dec <- tail(err.dec, 5)
                decflag <- -sign(diff(err.dec))

                decflag.sum <- sum(decflag)
                if(decflag.sum == 4) {
                    stepsize <- stepsize * rate.inc
                    epoch.last <- epoch
                } else if(decflag.sum <= 0) {
                #} else if( identical(decflag, c(-1,1,-1,1)) ) {
                    stepsize = stepsize * rate.dec
                    epoch.last = epoch
                }
            }
        }


        ## to optimise params in L1 and/or L4 by gradient decent method
        de.do5 <- anfis.dE.dO5(output.L5, target)
        do5.do4 <- anfis.dO5.dO4(output.L4)
        de.do4 <- anfis.dE.dO4(anfis, de.do5, do5.do4)
        de.dp4 <- anfis.dE.dP4(anfis, de.do4, output.L3, input.stack)
        do4.do3 <- anfis.dO4.dO3(output.L4, output.L4.mf)
        de.do3 <- anfis.dE.dO3(de.do4, do4.do3, output.L3)
        do3.do2 <- anfis.dO3.dO2(anfis, output.L2, output.L2.which)
        de.do2 <- anfis.dE.dO2(de.do3, do3.do2)
        do2.do1 <- anfis.dO2.dO1(anfis, output.L2, output.L1)
        de.do1 <- anfis.dE.dO1(anfis, output.L1, de.do2, do2.do1)
        de.dp1 <- anfis.dE.dP1(anfis, de.do1, input.stack)

        ## TODO: to decide whether to update params according to the configuration
        if(method[1] == "gradient") {
            anfis <- anfis.optimise.gradient(anfis, L1, de.dp1, stepsize)
        }

        if(method[2] == "gradient" || method[2] == "both") {
            anfis <- anfis.optimise.gradient(anfis, L4, de.dp4, stepsize)
        }

        if(online == 1 && p < data.num) {
            p <- p + 1
        } else {
            if(online == 1) { p <- 1 }
            epoch <- epoch + 1
        }
    }
    err.min.rec <- c(err.min.rec, err.min.rec.tmp)
    anfis.optimum <- append(anfis.optimum, list(RMSE_MIN_ROUTE=err.min.rec))
    show(err.min.rec)
    anfis.optimum
}

#lockBinding("anfis.optimise", e.ns)
#lockBinding("anfis.optimise", e.pkg)


anfis.optimise.gradient <- function(anfis, L, de.dp, stepsize) {

    de.dp.sum <- lapply(de.dp, apply, 2, sum)
    len = sqrt(sum(unlist(de.dp.sum)^2))

    if(len != 0) {
        eta <- stepsize / len
        for(i in 1:length(de.dp.sum)) {
            mf.params <- anfis$layer[[L]]$node[[i]]$mf.params - eta * de.dp.sum[[i]]
            mf.type <- anfis$layer[[L]]$node[[i]]$mf.type
            if((mf.type == "it2gbellmf" && sign(mf.params[3]) * abs(mf.params[1]) > sign(mf.params[3]) * abs(mf.params[2]))
                || (mf.type == "itlinearmf" && mf.params[1] > mf.params[2])) {
#                    mf.params[1] <- mf.params[2]
                    mf.params[c(1,2)] <- mf.params[c(2,1)]
            }
            anfis$layer[[L]]$node[[i]]$mf.params <- mf.params
        }
    }

    anfis
}


anfis.optimise.lse <- function(anfis, output.L3, input.stack, target, online=0, lambda=1) {

#    if(mean(sapply(output.L3, ncol)) != 1) {
#        stop("lse method here can only be used with type-1 anfis")
#    }

    L <- 4 + 1
    node.num <- length(anfis$layer[[L]]$node)

    input.x <- cbind(1, input.stack)
    data.num <- nrow(input.x)

    if(ncol(output.L3[[1]]) == 1) {
        A <- matrix(sapply(output.L3, function(w) c(w) * input.x), nrow=data.num)
    } else {
        A <- matrix(sapply(output.L3, function(w) apply(w,1,mean) * input.x), nrow=data.num)
    }

    B <- target

    if(!online || is.null(anfis$S)) {
    #if(is.null(anfis$S)) {
        S <- 999999 *  diag(ncol(A))
        if(!online) {
            theta <- matrix(rep(0, node.num * ncol(input.x)), ncol=1)
        } else {
            theta <- matrix(sapply(anfis$layer[[L]]$node, function(x) x$mf.params), ncol=1)
        }
    } else {
        S <- anfis$S
        theta <- matrix(sapply(anfis$layer[[L]]$node, function(x) x$mf.params), ncol=1)
    }

    for( i in 1:data.num ) {
        ai = A[i,]
        bi = B[i,]

        ## (S %*% ai) %*% (t(ai) %*% S) is much more efficient than (S %*% ai %*% t(ai) %*% S)
        S <- (1 / lambda) * (S - (S %*% ai) %*% (t(ai) %*% S) / c(lambda + (t(ai) %*% S %*% ai)))
        theta <- theta + S %*% ai %*% (t(bi) - t(ai) %*% theta)
    }

    theta <- matrix(theta, ncol=node.num)
    for( i in 1:node.num ) {
        anfis$layer[[L]]$node[[i]]$mf.params <- theta[,i]
    }

    if(online) {
        anfis$S <- S
    }

    anfis
}


anfis.dO3.dO2 <- function(anfis, output.L2, output.L2.which) {

    L <- 2 + 1
    node.num <- length(anfis$layer[[L]]$node)
    row.num <- nrow(output.L2[[1]])

    output.L2.sum <- list()
    output.L2.w <- list()
    for( i in 1:length(output.L2.which)) {
        tmp1 <- matrix(unlist(output.L2), nrow=row.num)
        tmp2 <- matrix(unlist(output.L2.which[[i]]), nrow=row.num)
        tmp3 <- t(sapply(1:row.num, function(idx) tmp1[idx,][tmp2[idx,]]))
        output.L2.w <- append(output.L2.w, list(tmp3))
        output.L2.sum <- append(output.L2.sum, list(matrix(apply(tmp3, 1, sum))))
    }

    do3.do2 <- list()
    for ( i in 1:node.num ) {
        #fan.out <- anfis$layer[[L]]$node[[i]]$fan.out
        do3.do2.tmp.k <- list()
        for( k in 1:length(output.L2.which) ) {
            do3.do2.tmp.j <- list()
            for ( j in 1:node.num ) {
                if(i == j) {
                    # output.L2.sum[[k]]^2 may be zero if output.L2.sum[[k]] is too small
                    #tmp <- (output.L2.sum[[k]] - output.L2.w[[k]][,i]) / output.L2.sum[[k]]^2
                    tmp <- (output.L2.sum[[k]] - output.L2.w[[k]][,i]) / output.L2.sum[[k]] / output.L2.sum[[k]]
                } else {
                    #tmp <- -output.L2.w[[k]][,j] / output.L2.sum[[k]]^2
                    tmp <- -output.L2.w[[k]][,j] / output.L2.sum[[k]] / output.L2.sum[[k]]
                }
                do3.do2.tmp.j <- append(do3.do2.tmp.j, list(c(tmp) * output.L2.which[[k]][[i]]))
            }
            do3.do2.tmp.k <- append(do3.do2.tmp.k, list(do3.do2.tmp.j))
        }

        do3.do2 <- append(do3.do2, list(do3.do2.tmp.k))
    }

    do3.do2
}

anfis.dMF.dP.gbellmf <- function(x, mf.params) {

    a = mf.params[1]
    b = mf.params[2]
    c = mf.params[3]

    denom = 0

    if (a == 0) {
        stop("anfis.dE.dP1.gbellmf: a == 0 founded!")
    }

    tmp1 = ((x - c)/a)^2
    tmp2 = tmp1^b
    denom = (1 + tmp2)^2

    dmf.dp = NULL
    # partial mf to partial a
    #tmp = (2*b*tmp2)/(a*denom)
    tmp = (2*b*tmp2)/(1 + tmp2)/a/(1 + tmp2)
    tmp[tmp2==Inf] = 0
    tmp[which(is.nan(tmp))] = 0
    tmp[tmp==Inf] = 0
    dmf.dp = cbind(dmf.dp , tmp)

    # partial mf to partial b
    #tmp = -1*log(tmp1)*tmp2/denom
    tmp = -1*log(tmp1)*tmp2/(1 + tmp2)/(1 + tmp2)
    tmp[tmp1==0] = 0
    tmp[which(is.nan(tmp))] = 0
    tmp[tmp==Inf] = 0
    dmf.dp = cbind(dmf.dp , tmp)

    # partial mf to partial c
    #tmp = (2*b*tmp2)/((x - c)*(denom))
    tmp = (2*b*tmp2)/(x - c)/(1 + tmp2)/(1 + tmp2)
    tmp[x==c] = 0
    tmp[which(is.nan(tmp))] = 0
    tmp[tmp==Inf] = 0
    dmf.dp = cbind(dmf.dp , tmp)

    colnames(dmf.dp) <- NULL
    dmf.dp
}


it2gbellmf <- function(mf.params) {
    
    if(length(mf.params) != 4) {
        stop("improper parameters for gaussian bell membership function")
    }

    a.lower <- mf.params[1]
    a.upper <- mf.params[2]
    b <- mf.params[3]
    c <- mf.params[4]

    gbellmf.lower <- function(x) {
        1 / ( 1 + (((x - c)/a.lower)^2)^b)
    }

    gbellmf.upper <- function(x) {
        1 / ( 1 + (((x - c)/a.upper)^2)^b)
    }

#    it2gbellmf <- function(x) {
#        u <- c(gbellmf.lower(x), gbellmf.upper(x))
#        u <- matrix(u, ncol=2)
#        u
#    }

    it2gbellmf <- c(gbellmf.lower, gbellmf.upper)
}



if(!require(inline)) {
    install.packages('inline', repos='http://cran.uk.r-project.org')
    require(inline)
}

inc <- '
/* This is taken from envir.c in the R 2.15.1 source 
   https://github.com/SurajGupta/r-source/blob/master/src/main/envir.c
*/
#define FRAME_LOCK_MASK (1<<14)
#define FRAME_IS_LOCKED(e) (ENVFLAGS(e) & FRAME_LOCK_MASK)
#define UNLOCK_FRAME(e) SET_ENVFLAGS(e, ENVFLAGS(e) & (~ FRAME_LOCK_MASK))
'

src <- '
  if (TYPEOF(env) == NILSXP)
    error("use of NULL environment is defunct");
  if (TYPEOF(env) != ENVSXP)
    error("not an environment");

  UNLOCK_FRAME(env);

  // Return TRUE if unlocked; FALSE otherwise
  SEXP result = PROTECT( Rf_allocVector(LGLSXP, 1) );
  LOGICAL(result)[0] = FRAME_IS_LOCKED(env) == 0;
  UNPROTECT(1);

  return result;
'

unlockEnvironment <- cfunction(signature(env = "environment"),
                        includes = inc,
                        body = src)

e.ns <- .getNamespace("FuzzyR")
e.pkg <- as.environment('package:FuzzyR')
e.exp <- e.ns$.__NAMESPACE__.$exports

tmp <- function (FUN, descend = TRUE) 
{
    if (is.function(FUN)) 
        return(FUN)
    if (!(is.character(FUN) && length(FUN) == 1L || is.symbol(FUN))) {
        FUN <- eval.parent(substitute(substitute(FUN)))
        if (!is.symbol(FUN)) 
            stop(gettextf("'%s' is not a function, character or symbol", 
                deparse(FUN)), domain = NA)
    }
    envir <- parent.env(environment())
    if (descend) 
        FUN <- get(as.character(FUN), mode = "function", envir = envir)
    else {
        FUN <- get(as.character(FUN), mode = "any", envir = envir)
        if (!is.function(FUN)) 
            stop(gettextf("found non-function '%s'", FUN), domain = NA)
    }
    return(FUN)
}

environment(tmp) <- e.ns

unlockEnvironment(e.ns)
e.ns$match.fun <- tmp
lockEnvironment(e.ns)

unlockEnvironment(e.pkg)
e.pkg$fuzzyr.match.fun <- e.ns$match.fun
lockEnvironment(e.pkg)

e.exp$fuzzyr.match.fun <- c(fuzzyr.match.fun="fuzzyr.match.fun")

unlockBinding("gensurf", e.ns)
unlockBinding("gensurf", e.pkg)

e.pkg$gensurf <- e.ns$gensurf <- function(fis, ix1=1, ix2=2, ox1=1) {
  i1= fis$input[[ix1]]
  i2= fis$input[[ix2]]
  o1= fis$output[[ox1]]
  
  x= seq(i1$range[1], i1$range[2], length=15)
  y= seq(i2$range[1], i2$range[2], length=15)
  m= as.matrix(expand.grid(x, y))

  o= evalfis(m, fis)
  z= matrix(o[,ox1], 15, 15, byrow=F)

  h= (z[-15,-15] + z[-1,-15] + z[-15,-1] + z[-1,-1]) / 4
  h= floor((h-min(h))/(max(h)-min(h))*14+.5)+1

  persp(x, y, z,
    xlab=i1$name, ylab=i2$name, zlab=o1$name,
    theta=-30, phi=30, col=rainbow(15)[16-h], ticktype='detailed')
}

lockBinding("gensurf", e.ns)
lockBinding("gensurf", e.pkg)

rm(tmp, e.ns, e.pkg, e.exp)


