works_with_R("3.0.2", rankSVMcompare="2013.9.3", ggplot2="0.9.3.1",
             microbenchmark="1.3.0")

source("svmlight.R")

funs <- list(l2=function(x)sum(x*x),
             l1=function(x)sum(abs(x))^2,
             linf=function(x)max(abs(x))^2)
deltas <- list(l2=function()runif(2,-1,1/2),
               l1=function()runif(2,-1/2,1/2),
               linf=function()runif(2,-1,1))
set.seed(1)
pair.sets <- list()
for(norm in names(funs)){
  delta.fun <- deltas[[norm]]
  f <- funs[[norm]]
  Xi <- c()
  Xip <- c()
  yi <- c()
  for(i in 1:1000){
    x <- runif(2,-2,2)
    delta <- delta.fun()
    xp <- x+delta
    noise <- 0
    noise <- rnorm(1,sd=1/4)#comment for noiseless simulation.
    fxdiff <- f(xp)-f(x)+noise
    y <- ifelse(fxdiff < -1, -1L,
                ifelse(fxdiff > 1, 1L, 0L))
    Xi <- rbind(Xi, x)
    Xip <- rbind(Xip, xp)
    rownames(Xi) <- rownames(Xip) <- NULL
    yi <- c(yi, y)
  }
  pair.sets[[norm]] <- list(Xi=Xi, Xip=Xip, yi=yi)
}
N <- 100
lapply(pair.sets, with, table(yi))
## Pick exactly N equality and inequality pairs.
pairs.picked <- list()
unused.sets <- list()
for(norm in names(pair.sets)){
  Pairs <- pair.sets[[norm]]
  is.zero <- Pairs$yi == 0
  equal <- which(is.zero)
  not.equal <- which(!is.zero)
  i <- c(sample(equal, N), sample(not.equal, N))
  is.unused <- !(seq_along(Pairs$yi) %in% i)
  pairs.picked[[norm]] <-
    list(Xi=Pairs$Xi[i,],
         Xip=Pairs$Xip[i,],
         yi=Pairs$yi[i])
  unused.sets[[norm]] <-
    list(Xi=Pairs$Xi[is.unused,],
         Xip=Pairs$Xip[is.unused,],
         yi=Pairs$yi[is.unused])
}
lapply(pairs.picked, with, table(yi))
lapply(unused.sets, with, table(yi))
## Plot the points.
point.df <- data.frame()
seg.df <- data.frame()
arrow.df <- data.frame()
for(norm in names(pairs.picked)){
  Pairs <- pairs.picked[[norm]]
  m <- with(Pairs, rbind(Xi, Xip))
  point.df <- rbind(point.df, data.frame(m, norm))
  yi <- Pairs$yi
  segs <- with(Pairs, data.frame(Xi, Xip))[yi == 0,]
  seg.df <- rbind(seg.df, data.frame(norm, segs))
  arrow.df <- with(Pairs,{
    rbind(arrow.df,
          data.frame(norm, Xip, Xi)[yi == -1,],
          data.frame(norm, Xi, Xip)[yi == 1,])
  })
}
library(grid)
basePlot <- ggplot(,aes(X1, X2))+
  facet_grid(.~norm)+
  theme_bw()+
  theme(panel.margin=unit(0,"cm"))+
  coord_equal()
pointPlot <- basePlot+
  geom_point(data=point.df)
print(pointPlot)
segPlot <- basePlot+
  aes(xend=X1.1,yend=X2.1)+
  geom_segment(data=seg.df)+
  geom_segment(data=arrow.df, arrow=arrow(type="closed",length=unit(0.05,"in")),
               color="red")
print(segPlot)
## Looks fine.
levs <- seq(-3,3,l=41)
X.grid <- as.matrix(expand.grid(x1=levs,x2=levs))
all.ranks <- data.frame()
## fit SVM.
unused.err <- data.frame()
train.pairs <- list()
for(norm in names(pairs.picked)){
  Pairs <- pairs.picked[[norm]]
  is.zero <- Pairs$yi == 0
  equal <- which(is.zero)
  not.equal <- which(!is.zero)
  i <- c(sample(equal, N/2), sample(not.equal, N/2))
  is.train <- seq_along(Pairs$yi) %in% i
  Pairs.train <- list(Xi=Pairs$Xi[is.train,],
                      Xip=Pairs$Xip[is.train,],
                      yi=Pairs$yi[is.train])
  train.pairs[[norm]] <- Pairs.train
  err.df <- data.frame()
  Cvals <- 4^seq(-4,5,by=1)
  models <- list()
  kvals <- 2^seq(-7, 1, by=1)
  model.df <- expand.grid(C=Cvals, k.width=kvals)
  for(model.i in 1:nrow(model.df)){
  ##for(model.i in 1:2){
    model <- model.df[model.i,]
    Cval <- model$C
    k.width <- model$k.width
    ker <- rbfdot(k.width)
    ##ker <- laplacedot(k.width)
    cat(sprintf("%4d / %4d C=%5.2f k.width=%5.2f\n",
                model.i, nrow(model.df), Cval, k.width))
    ## softCompareQP(Pairs.train, ker, C=Cval,cache=4000)
    ## microbenchmark(cache40=softCompareQP(Pairs.train, ker, C=Cval,cache=40),
    ##                cache400=softCompareQP(Pairs.train, ker, C=Cval,cache=4000),
    ##                times=10)
    fits <- list(compare=softCompareQP(Pairs.train, ker, C=Cval),
                 rank=svmlight(Pairs.train, Cval, k.width),
           rank2=svmlight(Pairs.train, Cval, k.width, equality="bothpairs"))
    models[[model.i]] <- fits
    ##X.grid$rank <- fit$rank(as.matrix(X.grid))
    ## pretty good train err:
    for(fit.name in names(fits)){
      fit <- fits[[fit.name]]
      yhat <- with(Pairs, fit$predict(Xi, Xip))
      yi <- Pairs$yi
      table(yi, yhat)
      sets <- list(train=is.train, test=!is.train)
      for(set in names(sets)){
        s <- sets[[set]]
        pred <- yhat[s]
        true <- yi[s]
        err.df <- rbind(err.df, {
          data.frame(FpFnInv(true, pred), Cval, k.width, set, norm, fit.name)
        })
      }
    }
  }
  ## train/test error curves.
  test.big <- subset(err.df, set=="test")
  test.dfs <- split(test.big, test.big$fit.name)
  f <- funs[[norm]]
  rank.df <- data.frame(X.grid, rank=apply(X.grid, 1, f), what="latent")
  chosen.df <- data.frame()
  unused <- unused.sets[[norm]]
  for(fit.name in names(test.dfs)){
    test.err <- test.dfs[[fit.name]]
    chosen <- which.min(test.err$error)
    chosen.df <- rbind(chosen.df, test.err[chosen,])
    fit <- models[[chosen]][[fit.name]]
    r <- fit$rank(X.grid)
    rank.df <- rbind(rank.df, {
      data.frame(X.grid, rank=r-min(r), what=fit.name)
    })
    yhat <- with(unused, fit$predict(Xi, Xip))
    unused.err <- rbind(unused.err, {
      data.frame(norm, fit.name, FpFnInv(unused$yi, yhat))
    })
  }
  normContour <- ggplot(rank.df, aes(x1, x2, z=rank))+
    geom_contour(colour="black")+
    coord_equal()+
    theme_bw()+
    theme(panel.margin=unit(0,"cm"))+
    facet_grid(.~what)
  print(normContour)
  overfitPlot <- ggplot(err.df, aes(log2(Cval), error, colour=fit.name))+
    geom_line(aes(group=interaction(set, fit.name), linetype=set))+
    facet_wrap("k.width")+
    theme_bw()+
    theme(panel.margin=unit(0,"cm"))+
    geom_point(data=chosen.df)
  print(overfitPlot)
  ## if the optimal model occurs on the min/max of the tested
  ## hyperparameters, then this is probably sub-optimal and we need to
  ## define a larger grid.
  stopifnot(! test.err[chosen,"Cval"] %in% range(Cvals))
  stopifnot(! test.err[chosen,"k.width"] %in% range(kvals))
    
  all.ranks <- rbind(all.ranks, data.frame(rank.df, norm))
}

simulation <- list(rank=all.ranks, error=unused.err, train.pairs=train.pairs)

save(simulation, file="simulation.RData")

