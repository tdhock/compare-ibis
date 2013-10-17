works_with_R("3.0.2", plyr="1.8", reshape2="1.2.2")

source("tikz.R")

load("simulation.samples.RData")

## matrix versions of the norm.
funs <- list(l2=function(x)rowSums(x*x),
             l1=function(x)rowSums(abs(x))^2,
             linf=function(x)apply(abs(x), 1, max)^2)

size.list <- simulation.samples$data
err <- simulation.samples$err
err$percent <- err$error / err$count * 100
## sets of training data and bayes error on test data.
sets <- dcast(err, N + seed + norm ~ fit.name, value.var="percent")
sets$diff <- sets$compare-sets$rank
sets$set.id <- 1:nrow(sets)
diff.df <- ddply(sets, .(N, norm), summarize,
                   N=N[1], norm=norm[1],
                   mean=mean(diff), sd=sd(diff))
train.df <- data.frame()
bayes.df <- data.frame()
for(set.id in sets$set.id){
  e <- sets[set.id,]
  N <- as.character(e$N)
  norm <- as.character(e$norm)
  seed <- as.character(e$seed)
  err$set.id[err$norm == norm & err$N == N & err$seed == seed] <- set.id
  set.list <- size.list[[N]][[seed]][[norm]]
  info <- data.frame(N=as.integer(as.character(N)), norm, seed, set.id)
  ## The Bayes error on the test data set.
  test <- set.list$test
  fun <- funs[[norm]]
  fxdiff <- with(test, fun(Xip)-fun(Xi))
  yhat <- ifelse(fxdiff > 1, 1L,
                 ifelse(fxdiff < -1, -1L, 0))
  table(yhat, test$yi)
  percent <- mean(yhat != test$yi) * 100
  bayes.df <- rbind(bayes.df, data.frame(info, percent))
  ## Train pairs, oriented in the same way:
  pair.df <- with(set.list$train,{
    rbind(data.frame(Xt=Xi[yi==1,],Xtp=Xip[yi==1,],yt=1),
          data.frame(Xt=Xip[yi==-1,],Xtp=Xi[yi==-1,],yt=1),
          data.frame(Xt=Xi[yi==0,],Xtp=Xip[yi==0,],yt=-1))
  })
  train.df <- rbind(train.df, data.frame(pair.df, info))
}
bayes.df$fit.name <- "latent"
combined <- rbind(err[,names(bayes.df)],
                  bayes.df)
percents <-
  ddply(combined, .(N, fit.name, norm), summarize,
        mean=mean(percent),
        sd=sd(percent),
        se=sd(percent)/sqrt(length(percent)))

model.colors <- c(
                  rank="#f8766d",
                  compare="#00bfc4", #bluish
                  latent="grey")
model.labels <- c("rank (\\ref{eq:svmrank})",
                  "compare (\\ref{eq:svm-dual})","latent $r$")
model.labels <- names(model.colors)
library(grid)
percents$fit.name <- factor(percents$fit.name, names(model.colors))
labels <- c(l1="||x||_1^2",
            l2="||x||_2^2",
            linf="||x||_\\infty^2")
makelabel <- function(x)sprintf("$r(x) = %s$", labels[as.character(x)])
percents$label <- makelabel(percents$norm)
err$label <- makelabel(err$norm)
leg <- "learned\nfunction"
library(directlabels)
boring <- ggplot(percents, aes(N, mean, group=fit.name))+
  geom_vline(xintercept=100,size=2)+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),alpha=1/3)+
  geom_dl(aes(label=fit.name),method="last.qp")+
  geom_line(lwd=1.5)+
  ## Plot actual data:
  ##geom_point(aes(N, error/count*100, colour=fit.name), data=err)+
  facet_grid(.~label)+
  theme_bw()+
  theme(panel.margin=unit(0,"cm"))+
  ##scale_colour_manual(leg,values=model.colors,labels=model.labels)+
  ##scale_fill_manual(leg,values=model.colors,labels=model.labels)+
  ylab("percent incorrectly\npredicted test pairs")+
  scale_x_continuous("$n=$ number of labeled pairs in the training set",
                     breaks=c(50, 200, 400, 800))+
  coord_cartesian(xlim=c(50, 1000))

tikz("figure-simulation-samples.tex",h=2)
print(boring)
dev.off()
