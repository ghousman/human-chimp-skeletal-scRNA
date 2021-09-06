

my.limmafit <- 
  function (exprs, groupid, compid, norm.factor.method = "TMM", voom.normalize.method = "cyclicloess") 
  {
    compnum <- nrow(compid)
    genenum <- nrow(exprs)
    limmat <- matrix(0, genenum, compnum)
    
    colnames(limmat) <- rownames(compid)
    rownames(limmat) <- rownames(exprs)
    
    limmas2 <- rep(0, compnum)
    limmadf <- rep(0, compnum)
    limmav0 <- rep(0, compnum)
    limmag1num <- rep(0, compnum)
    limmag2num <- rep(0, compnum)
    for (i in 1:compnum) {
      message(paste("Running limma for comparision",i,"/",compnum))
      selid1 <- which(groupid == compid[i, 1])
      selid2 <- which(groupid == compid[i, 2])
      # make a new count data frame
      counts <- cbind(exprs[, selid1], exprs[, selid2])
      
      # remove NAs 
      not.nas <- which(apply(counts, 1, function(x) !any(is.na(x))) == TRUE)
      
      # runn voom/limma
      d <- DGEList(counts[not.nas,])
      d <- calcNormFactors(d, method = norm.factor.method)
      g1num <- length(selid1)
      g2num <- length(selid2)
      designmat <- cbind(base = rep(1, (g1num + g2num)), delta = c(rep(0, 
                                                                       g1num), rep(1, g2num)))
      
      y <- voom(d, designmat, normalize.method = voom.normalize.method)
      fit <- lmFit(y, designmat)
      fit <- eBayes(fit)
       
      limmat[not.nas, i] <- fit$t[, 2]
      limmas2[i] <- fit$s2.prior
      limmadf[i] <- fit$df.prior
      limmav0[i] <- fit$var.prior[2]
      limmag1num[i] <- g1num
      limmag2num[i] <- g2num
    }
    limmacompnum <- nrow(compid)
    result <- list(t = limmat, 
                   v0 = limmav0, 
                   df0 = limmadf, 
                   s20 = limmas2,
                   g1num = limmag1num, 
                   g2num = limmag2num, 
                   compnum = limmacompnum)
  }

my.generatetype <- 
  function (limfitted) 
  {
    jtype <- list()
    df <- limfitted$g1num + limfitted$g2num - 2 + limfitted$df0
    for (j in 1:limfitted$compnum) {
      jtype[[j]] <- list(f0 = Cormotif:::modt.f0.loglike, f0.param = df[j], 
                         f1 = Cormotif:::modt.f1.loglike, f1.param = c(df[j], limfitted$g1num[j], 
                                                                       limfitted$g2num[j], limfitted$v0[j]))
    }
    jtype
  }



my.cmfit <- 
  function (x, type, K = 1, tol = 0.001, max.iter = 100) 
  {
    xrow <- nrow(x)
    xcol <- ncol(x)
    loglike0 <- list()
    loglike1 <- list()
    p <- rep(1, K)/K
    q <- matrix(runif(K * xcol), K, xcol)
    q[1, ] <- rep(0.01, xcol)
    for (i in 1:xcol) {
      f0 <- type[[i]][[1]]
      f0param <- type[[i]][[2]]
      f1 <- type[[i]][[3]]
      f1param <- type[[i]][[4]]
      loglike0[[i]] <- f0(x[, i], f0param)
      loglike1[[i]] <- f1(x[, i], f1param)
    }
    condlike <- list()
    for (i in 1:xcol) {
      condlike[[i]] <- matrix(0, xrow, K)
    }
    loglike.old <- -1e+10
    for (i.iter in 1:max.iter) {
      if ((i.iter%%50) == 0) {
        print(paste("We have run the first ", i.iter, " iterations for K=", 
                    K, sep = ""))
      }
      err <- tol + 1
      clustlike <- matrix(0, xrow, K)
      templike <- matrix(0, xrow, 2)
      for (j in 1:K) {
        for (i in 1:xcol) {
          templike[, 1] <- log(q[j, i]) + loglike1[[i]]
          templike[, 2] <- log(1 - q[j, i]) + loglike0[[i]]
          tempmax <- pmax(templike[, 1], templike[, 2])
          for (z in 1:2) {
            templike[, z] <- exp(templike[, z] - tempmax)
          }
          tempsum <- templike[, 1] + templike[, 2]
          clustlike[, j] <- clustlike[, j] + tempmax + 
            log(tempsum)
          condlike[[i]][, j] <- templike[, 1]/tempsum
        }
        clustlike[, j] <- clustlike[, j] + log(p[j])
      }
      tempmax <- apply(clustlike, 1, max)
      for (j in 1:K) {
        clustlike[, j] <- exp(clustlike[, j] - tempmax)
      }
      tempsum <- apply(clustlike, 1, sum)
      for (j in 1:K) {
        clustlike[, j] <- clustlike[, j]/tempsum
      }
      p.new <- (apply(clustlike, 2, sum) + 1)/(xrow + K)
      q.new <- matrix(0, K, xcol)
      for (j in 1:K) {
        clustpsum <- sum(clustlike[, j])
        for (i in 1:xcol) {
          q.new[j, i] <- (sum(clustlike[, j] * condlike[[i]][, 
                                                             j]) + 1)/(clustpsum + 2)
        }
      }
      err.p <- max(abs(p.new - p)/p)
      err.q <- max(abs(q.new - q)/q)
      err <- max(err.p, err.q)
      loglike.new <- (sum(tempmax + log(tempsum)) + sum(log(p.new)) + 
                        sum(log(q.new) + log(1 - q.new)))/xrow
      p <- p.new
      q <- q.new
      loglike.old <- loglike.new
      if (err < tol) {
        break
      }
    }
    clustlike <- matrix(0, xrow, K)
    for (j in 1:K) {
      for (i in 1:xcol) {
        templike[, 1] <- log(q[j, i]) + loglike1[[i]]
        templike[, 2] <- log(1 - q[j, i]) + loglike0[[i]]
        tempmax <- pmax(templike[, 1], templike[, 2])
        for (z in 1:2) {
          templike[, z] <- exp(templike[, z] - tempmax)
        }
        tempsum <- templike[, 1] + templike[, 2]
        clustlike[, j] <- clustlike[, j] + tempmax + log(tempsum)
        condlike[[i]][, j] <- templike[, 1]/tempsum
      }
      clustlike[, j] <- clustlike[, j] + log(p[j])
    }
    tempmax <- apply(clustlike, 1, max)
    for (j in 1:K) {
      clustlike[, j] <- exp(clustlike[, j] - tempmax)
    }
    tempsum <- apply(clustlike, 1, sum)
    for (j in 1:K) {
      clustlike[, j] <- clustlike[, j]/tempsum
    }
    p.post <- matrix(0, xrow, xcol)
    for (j in 1:K) {
      for (i in 1:xcol) {
        # clustlike is the likelihood that a gene belongs to a cluster
        # condlike is the likelihood of DE in tissue i given membership in a cluster
        p.post[, i] <- p.post[, i] + clustlike[, j] * condlike[[i]][, 
                                                                    j]
      }
    }
    loglike.old <- loglike.old - (sum(log(p)) + sum(log(q) + 
                                                      log(1 - q)))/xrow
    loglike.old <- loglike.old * xrow
    result <- list(p.post = p.post, motif.prior = p, motif.q = q, 
                   loglike = loglike.old, clustlike=clustlike, condlike = condlike)
  }




cormotiffit.counts <- 
  function (exprs, groupid, compid, K = 1, tol = 0.001, max.iter = 100, 
            BIC = TRUE, norm.factor.method = "TMM", voom.normalize.method = "cyclicloess") 
  {
    limfitted <- my.limmafit(exprs, groupid, compid, norm.factor.method = norm.factor.method, voom.normalize.method = voom.normalize.method) 
    jtype <- my.generatetype(limfitted)
    fitresult <- list()
    for (i in 1:length(K)) {
      message(paste("Fitting cormotif model for K =", K[i]))
      fitresult[[i]] <- my.cmfit(limfitted$t, 
                                 type = jtype, K = K[i], max.iter = max.iter, tol = tol)
    }
    bic <- rep(0, length(K))
    aic <- rep(0, length(K))
    loglike <- rep(0, length(K))
    for (i in 1:length(K)) loglike[i] <- fitresult[[i]]$loglike
    for (i in 1:length(K)) bic[i] <- -2 * fitresult[[i]]$loglike + 
      (K[i] - 1 + K[i] * limfitted$compnum) * log(dim(exprs)[1])
    for (i in 1:length(K)) aic[i] <- -2 * fitresult[[i]]$loglike + 
      2 * (K[i] - 1 + K[i] * limfitted$compnum)
    if (BIC == TRUE) {
      bestflag = which(bic == min(bic))
    } else {
      bestflag = which(aic == min(aic))
    }
    result <- list(bestmotif = fitresult[[bestflag]], bic = cbind(K, 
                                                                  bic), aic = cbind(K, aic), loglike = cbind(K, loglike), allmotifs = fitresult)
  }

