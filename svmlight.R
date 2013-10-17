pick.best.index <- function(err){
  nparam <- length(err)
  candidates <- which(err==min(err))
  if(length(err)==1)return(candidates)
  st <- abs(median(candidates)-candidates)
  middle <- candidates[which.min(st)]
  if(all(diff(err)==0))return(middle)
  if(nparam %in% candidates && 1 %in% candidates){
    cat("Warning: strange error profile, picking something near the center\n")
    print(as.numeric(err))
    d <- diff(candidates)>1
    if(any(d)){
      which(d)[1]
    }else{
      middle
    }
  }else if(1 %in% candidates){
    max(candidates)
  }else if(nparam %in% candidates){
    min(candidates)
  }else {
    middle
  }
### Integer index of the minimal error.
}

mat2lines <- function(x){
  stopifnot(is.matrix(x))
  m <- matrix(paste0(col(x),":",x),nrow(x))
  apply(m, 1, paste0, collapse=" ")
}

pairs2lines <- function(equality, L){
  y <- L$yi
  flip <- y==-1
  tmp <- L$Xi[flip,]
  L$Xi[flip,] <- L$Xip[flip,]
  L$Xip[flip,] <- tmp
  y[flip] <- -y[flip]
  stopifnot(all(y %in% c(0, 1)))
  Xi <- mat2lines(L$Xi)
  Xip <- mat2lines(L$Xip)
  plines <- c()
  plines <- makeLines(Xi[y==1], 0, Xip[y==1], 1, plines)
  if(equality == "bothpairs"){
    plines <- makeLines(Xi[y==1], 0, Xip[y==1], 1, plines)
    plines <- makeLines(Xi[y==0], 0, Xip[y==0], 1, plines)
    plines <- makeLines(Xi[y==0], 1, Xip[y==0], 0, plines)
  }else if(equality=="nopairs"){
    ## this is OK, do nothing.
  }else{
    stop("undefined pairs for equality ",equality)
  }
  plines
}

makeLine <- function(z, qid, x){
  sprintf("%d qid:%d %s", z, qid, x)
}

makeLines <- function(Xi, zi, Xip, zip, plines){
  first.qid <- length(plines) + 1
  last.qid <- first.qid + length(Xi) - 1
  qid <- first.qid:last.qid
  m <- cbind(makeLine(zi, qid, Xi),
             makeLine(zip, qid, Xip))
  newlines <- apply(m, 1, paste0, collapse="\n")
  c(plines,newlines)
}
  
svmlight <- function
### Fit an SVMlight comparison model: min_{w,xi} w'w/2 + C sum xi_{i}
### such that w'x_i > w'x_ip + 1 - xi_i for all pairs of layouts x_i,
### x_ip with inequality constraints. Assumes that svm_learn and
### svm_classify programs are in the current directory.
(pairs,
### List of pairs of training layouts, with labels.
 Cvalue=1,
### positive numeric scalar: slack penalization coefficient in the
### objective function.
 kernel.width=1,
### positive numeric scalar: kernel width parameter.
 equality="nopairs",
### When y_i=0 we can either ignore it (equality="nopairs") or treat
### it as two pairs, one in each direction (equality="bothpairs").
 filebase=tempfile(),
### Base file name of svmlight interface files.
 verbose=0
 ){
  for(scalar in list(Cvalue, kernel.width)){
    stopifnot(is.numeric(scalar))
    stopifnot(scalar > 0)
    stopifnot(length(scalar) == 1)
  }
  ## Reality checks about Xi, Xip, yi.
  stopifnot(is.matrix(pairs$Xi))
  stopifnot(is.matrix(pairs$Xip))
  stopifnot(dim(pairs$Xi)==dim(pairs$Xip))
  stopifnot(is.numeric(pairs$Xi))
  stopifnot(is.numeric(pairs$Xip))
  stopifnot(nrow(pairs$Xi)==length(pairs$yi))
  stopifnot(is.numeric(pairs$yi))
  stopifnot(pairs$yi %in% c(-1,0,1))

  ## Scale input data before throwing it to SVMlight!
  X.scaled <- scale(with(pairs, rbind(Xi, Xip)))
  smat <- function(one.row, target.matrix){
    stopifnot(is.numeric(one.row))
    stopifnot(is.matrix(target.matrix))
    stopifnot(length(one.row) == ncol(target.matrix))
    matrix(one.row, nrow(target.matrix), ncol(target.matrix),
           byrow=TRUE)
  }
  mu.mat <- smat(attr(X.scaled,"scaled:center"), pairs$Xi)
  sigma.mat <- smat(attr(X.scaled,"scaled:scale"), pairs$Xi)
  pairs$Xi <- (pairs$Xi - mu.mat)/sigma.mat
  pairs$Xip <- (pairs$Xip - mu.mat)/sigma.mat

  ## temporary files we will use.
  exfile <- sprintf("%s.examples",filebase)
  trainfile <- sprintf("%s.train",filebase)
  mfile <- sprintf("%s.model",filebase)
  pfile <- sprintf("%s.predictions",filebase)

  ## Vectorized calculation of lines to export to svmlight.
  plines <- pairs2lines(equality, pairs)

  cat(plines, file=trainfile, sep="\n")
  cmd <- sprintf("./svm_learn -z p -c %f -t 2 -g %f %s %s",
                 Cvalue, kernel.width, trainfile, mfile)
  ##cmd <- sprintf("./svm_learn -z p -c %f %s %s", Cvalue, trainfile, mfile)

  ## From svm_rank web page: You can in principle use kernels in
  ## SVMrank using the '-t' option just like in SVMlight, but it is
  ## painfully slow and you are probably better off using SVMlight.
  ##cmd <- sprintf("./svm_rank_learn -c %f %s %s", Cvalue, trainfile, mfile)
  if(verbose>0){
    print(cmd)
  }
  system(cmd,ignore.stdout=verbose<2)

  ### This is the learned function.
  f <- function(X){
    if(is.data.frame(X)){
      X <- as.matrix(X)
    }
    stopifnot(is.matrix(X))
    stopifnot(is.numeric(X))
    ## verify column names!
    stopifnot(colnames(X) == colnames(pairs$Xi))
    ## before applying the learned function, scale the test data using
    ## the same constants we found during training.
    mu.test <- smat(attr(X.scaled,"scaled:center"), X)
    sigma.test <- smat(attr(X.scaled,"scaled:scale"), X)
    Xtest <- (X-mu.test)/sigma.test
    test.lines <- makeLine(0, 0, mat2lines(Xtest))
    cat(test.lines, file=exfile, sep="\n")
    ##testcmd <- sprintf("./svm_rank_classify %s %s %s", exfile, mfile, pfile)
    testcmd <- sprintf("./svm_classify %s %s %s", exfile, mfile, pfile)
    system(testcmd,ignore.stdout=TRUE)
    scan(pfile, what=numeric(), quiet=TRUE)
  }
  ## For visualization, it is useful to have all the level curves on
  ## the same scale so we don't have to update the legend and
  ## colours. So here we apply a linear transformation to the output
  ## function values and we will use that instead of f.
  f.values <- with(pairs, c(f(pairs$Xip),f(pairs$Xi)))
  mu <- mean(f.values)
  sigma <- sd(f.values)
  fnorm <- function(X){
    (f(X) - mu)/sigma
  }

  ## Read the model file for interpretation of support vectors.
  mlines <- readLines(mfile)
  svec.lines <- sub(" #","",gsub("[0-9]+:","",mlines[12:length(mlines)]))
  stopifnot(length(svec.lines)>1)
  cvecs <- strsplit(svec.lines,split=" ")
  model.mat <- sapply(cvecs, as.numeric)
  X.sv <- t(model.mat[-1,])
  colnames(X.sv) <- colnames(pairs$Xi)
  yalpha <- model.mat[1,]
  sv.type <- factor(ifelse(Cvalue != abs(yalpha),"notC (margin)",
                           ifelse(sign(yalpha)>0,
                                  "C (wrongside)","-C (wrongside)")))
  sv.df <- data.frame(X.sv, yalpha, sv.type)

  ## Calculate the optimal threshold using grid search.
  fxdiff <- with(pairs, fnorm(pairs$Xip)-fnorm(pairs$Xi))
  adiff <- abs(fxdiff)
  min.abs <- min(adiff)
  if(min.abs == 0) min.abs <- 1e-6
  ## Do grid search using these 100 possible thresholds...
  tgrid <- 10^seq(log10(min.abs),log10(max(adiff)),l=100)
  ## ...and the differences found in the data.
  thresholds <- sort(unique(c(tgrid, adiff)))
  rug.df <- data.frame(threshold=adiff,
                       label=ifelse(pairs$yi==0,"same","different"))
  err.df <- data.frame()
  for(threshold in thresholds){
    yhat <- ifelse(fxdiff < -threshold, -1,
                   ifelse(fxdiff > threshold, 1, 0))
    fp <- sum(pairs$yi == 0 & yhat != 0)
    inv <- sum(pairs$yi %in% c(-1,1) & yhat != pairs$yi & yhat %in% c(-1,1))
    fn <- sum(pairs$yi!=0 & yhat==0)
    err <- sum(pairs$yi != yhat)
    stopifnot(err == fp+inv+fn)
    df <- data.frame(threshold,
                     false.positive=fp/length(yhat),
                     false.negative=fn/length(yhat),
                     inversion=inv/length(yhat),
                     error=err/length(yhat))
    err.df <- rbind(err.df, df)
  }
  thresh <- err.df[pick.best.index(err.df$error),"threshold"]
  list(rank=fnorm,predict=function(Xi, Xip){
    fxdiff <- fnorm(Xip)-fnorm(Xi)
    ifelse(fxdiff < -thresh, -1L,
           ifelse(fxdiff > thresh, 1L, 0L))
  },threshold=thresh,threshold.error=err.df,sv.df=sv.df)
### List of model fit results. You can do fit$rank(X) with a feature
### matrix X to get a vector of learned model predictions (numeric
### vector) or fit$compare(Xi, Xip) to get comparison values (integer
### vector with elements -1,0,1).
}
