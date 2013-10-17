works_with_R("3.0.1",rankSVMcompare="2013.9.3",quadmod="2013.8.23",
             proto="0.3.10")

data(separable)
source("tikz.R")

only <- 1:50
one <- with(separable, list(yi=yi[only],Xi=Xi[only,],Xip=Xip[only,]))
is.one <- one$yi==1
is.zero <- one$yi == 0
Di <- with(one, Xip - Xi)
Di[,"angle"] <- Di[,"angle"]*2
Di.other <- Di
Di.other[is.zero,] <- -Di.other[is.zero,]
setdf <- function(yi, m,set){
  data.frame(yi, m, set)
}
Di.both <- rbind(Di,Di.other[is.zero,])
yi.both <- c(one$yi, rep(0, sum(is.zero)))
neg <- yi.both == -1
yi.both[neg] <- 1
Di.both[neg,] <- -Di.both[neg,]
both.scaled <- Di.both
bsd <- both.scaled[,"distance"]
both.scaled[,"distance"] <- bsd/100
point.df <- rbind(setdf(one$yi,Di,"one"),
                  setdf(one$yi,Di.other,"other"),
                  setdf(yi.both,Di.both,"both"),
                  setdf(yi.both,both.scaled,"scaled"))
svm.color <- "black"
lp.color <- "grey50"
yi.colors <- c("0"="#f8756e", #orange
               "1"="#00ba38", #green
               "-1"="#619cff",
               "2"=svm.color) #blue
library(grid)
p <- ggplot()+
  geom_point(aes(distance, angle, colour=factor(yi)), data=point.df)+
  facet_grid(.~set)+
  scale_colour_manual(values=yi.colors)+
  theme_bw()+
  theme(panel.margin=unit(0,"cm"))
print(p)

sets <- split(point.df, point.df$set)

model.points <- data.frame()
model.segs <- data.frame()
model.lines <- data.frame()
model.sv <- data.frame()
blank.df <- data.frame()
for(set.name in c("one","both","scaled")){
  pairs <- sets[[set.name]]

  diffs <- pairs[,-4]

  vars <- make.ids(margin=1, weight=2)
  constraints <- list()
  for(i in 1:nrow(diffs)){
    if(diffs$yi[i] == 0){
      right.side <- -1
      yi.vec <- c(-1,1)
    }else{
      right.side <- 1
      yi.vec <- diffs$yi[i]
    }
    di <- unlist(diffs[i,-1])
    for(yi in yi.vec){
      const <- with(vars,{
        weight*di*yi + margin*-1 >= right.side
      })
      constraints <- c(constraints,list(const))
    }
  }
  n.vars <- length(unlist(vars))
  d <- rep(0, n.vars)
  d[vars$margin] <- 1
  sol <- run.lpSolveAPI(vars, d, constraints)

  fxdiff <- as.matrix(diffs[,-1]) %*% sol$weight
  thresh <- function(x)ifelse(x>1,1,ifelse(abs(x)<1,0,-1))
  ## check to make sure we have perfect prediction.
  stopifnot(thresh(fxdiff) == pairs$yi)
  margin <- ifelse(pairs$yi==0,{
    1-abs(fxdiff)
  },{
    -1 + pairs$yi * fxdiff
  })
  on.margin <- abs(margin - sol$margin)<1e-4
  diffs$constraint <- ifelse(on.margin, "active", "inactive")

  slope <- with(sol, -weight[1]/weight[2])
  linedf <- function(line, intercept){
    data.frame(slope, intercept, line)
  }
  line.df <-
    rbind(linedf("decision",c(-1,1)/sol$weight[2]),
          linedf("margin",(c(-1,1,-1,1)*sol$margin+c(1,1,-1,-1))/sol$weight[2]))

  dots <- ggplot()+
    geom_point(aes(distance,angle,colour=factor(yi),
                   size=constraint),
               data=diffs)+
    scale_size_manual(values=c(active=1,inactive=1), guide="none")+
    scale_colour_manual(values=yi.colors)
  arange <- range(diffs$angle)
  seg <- function(v, line){
    d <- with(sol, (v-weight[2]*arange)/weight[1])
    data.frame(t(c(distance=d, angle=arange)), line)
  }
  seg.df <- rbind(seg(1-sol$margin,"margin"),
                  seg(1+sol$margin,"margin"),
                  seg(-1-sol$margin,"margin"),
                  seg(-1+sol$margin,"margin"),
                  seg(1,"decision"),
                  seg(-1,"decision"))
  p <- dots+
    scale_size_manual(values=c(active=2,inactive=1))+
    geom_segment(aes(distance1,angle1,xend=distance2,yend=angle2,
                     linetype=line),data=seg.df)+
    scale_linetype_manual(values=c(decision="solid",margin="dotted"))
  print(p)

  model.points <- rbind(model.points, {
    data.frame(diffs, set.name)
  })
  blank.df <- rbind(blank.df, {
    data.frame(rbind(diffs[,2:3], -diffs[,2:3]), set.name, model="compare")
  })
  ## model.segs <- rbind(model.segs, {
  ##   data.frame(seg.df, set.name, model="compare")
  ## })
  model.lines <- rbind(model.lines, {
    data.frame(line.df, set.name, model="compare")
  })
}

for(set.name in c("both","scaled")){
  pairs <- sets[[set.name]]
  ## Hard margin SVM.
  stopifnot(pairs$yi %in% c(0, 1))
  svm.y <- ifelse(pairs$yi==0, -1, 1)
  X <- as.matrix(pairs[,2:3])
  K <- X %*% t(X)
  N <- nrow(X)
  ## Construct the QP constraints using the quadmod modeling language.
  vars <- make.ids(alpha=N)
  constraints <-
    c(vars$alpha[]*svm.y >= 0,
      list(sum(vars$alpha) == 0))
  diag(K) <- diag(K)+1e-9
  sol <- run.quadprog(vars, K, svm.y, constraints)
  is.sv <- abs(sol$alpha)>1e-4
  a.sv <- sol$alpha[is.sv]
  X.sv <- X[is.sv,]
  qplot(distance, angle, data=pairs, colour=factor(yi))+
    geom_point(aes(colour=NULL),data=data.frame(X.sv), pch=1)
  y.sv <- svm.y[is.sv]
  sv.evals <- (a.sv * X.sv) %*% t(X.sv)
  bias.values <- y.sv - colSums(sv.evals)
  bias <- mean(bias.values)
  ## calc svm decision boundary and margin.
  weight <- -colSums(X.sv * a.sv / bias)
  mu <- -1/bias
  arange <- range(pairs$angle)
  seg <- function(v, line){
    d <- (v-weight[2]*arange)/weight[1]
    data.frame(t(c(distance=d, angle=arange)), line)
  }
  seg.df <- rbind(seg(1-mu,"margin"),
                  seg(1+mu,"margin"),
                  ## seg(-1-mu,"margin"),
                  ## seg(-1+mu,"margin"),
                  ## seg(-1,"decision"),
                  seg(1,"decision"))
  model.segs <- rbind(model.segs, {
    data.frame(seg.df, set.name)
  })
  model.sv <- rbind(model.sv, {
    data.frame(yi=2, X.sv, constraint="sv", set.name)
  })
}

plotted.points <- rbind(model.points, model.sv)
qp <- data.frame(x=200,y=3,label="QP",set.name="both")
lp <- data.frame(x=-200,y=0,label="LP",set.name="both")
lab.df <-
  data.frame(set.name=c("one","one","one","scaled","scaled"),
             yi=c(-1,0,1,0,1),
             distance=c(-190, -50, 180, -0.3, 1.8),
             angle=c(-1.5, 1.9, 2.4, -1.8, 3),
             label=c("$y_i=-1$","$y_i=0$","$y_i=1$",
               "$\\tilde y_i=-1$", "$\\tilde y_i=1$"))
mplot <- ggplot()+
  geom_abline(aes(slope=slope,intercept=intercept,linetype=line),
              data=model.lines, size=1, color=lp.color)+
  geom_segment(aes(distance1,angle1,xend=distance2,yend=angle2,
                    linetype=line),data=model.segs,
               colour=svm.color, size=0.8)+
  geom_blank(aes(distance, angle),data=blank.df)+
  geom_text(aes(x,y,label=label),data=qp,colour=svm.color,size=3)+
  geom_text(aes(x,y,label=label),data=lp,colour=lp.color,size=3)+
  geom_point(aes(distance,angle,colour=factor(yi),shape=constraint),
             data=plotted.points, size=1.8)+
  scale_colour_manual("label",values=yi.colors,
                      labels=c("$y_i=-1$",
                        "$y_i=0$, $\\tilde y_i=-1$",
                        "$y_i=1$, $\\tilde y_i=1$",
                        "QP support vector"))+
  ##scale_size_manual(values=c(active=2,inactive=1))+
  scale_shape_manual("point",values=c(active=8,inactive=1, sv=20),
                     labels=c("LP constraint active",
                       "LP constraint inactive",
                       "QP support vector"))+
  scale_linetype_manual(values=c(decision="solid",margin="dotted"),
                        labels=c("decision $r(x)=\\pm 1$",
                          "margin $r(x)=\\pm 1\\pm\\mu$"))+
  geom_text(aes(distance, angle, colour=factor(yi), label=label),
            data=lab.df, size=3)+
  facet_grid(.~set.name,scales="free",labeller=function(var,val){
    c(one="original $x'-x$",
      both="flipped $\\tilde x'-\\tilde x$",
      scaled="scaled and flipped $\\tilde x'-\\tilde x$")[val]
  })+
  theme_bw()+
  ##geom_point(aes(distance, angle), data=model.sv, size=1, pch=1)+
  theme(panel.margin=unit(0,"cm"))+
  xlab("difference feature 1")+
  ylab("difference feature 2")+
  guides(colour="none")
print(mplot)
tikz("figure-hard-margin.tex",h=2.8)
print(mplot)
dev.off()
