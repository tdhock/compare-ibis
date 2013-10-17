load("simulation.RData")

source("tikz.R")

fun.levs <- c("training data", "rank", "rank2", "compare")
model.colors <- c()
eq.lab <- "equality\npair\n$y_i=0$"
ineq.lab <- "inequality\npair\n$y_i\\in\\{-1,1\\}$"
model.colors[[eq.lab]] <- "black"
model.colors[[ineq.lab]] <- "red"
model.colors <- c(model.colors, {
  c(latent="grey",
    rank="#f8766d",
    rank2="green",
    compare="#00bfc4") #bluish
})
what.levs <- names(model.colors)

norm.list <- simulation$train
seg.df <- data.frame()
arrow.df <- data.frame()
for(norm in names(norm.list)){
  Pairs <- norm.list[[norm]]
  m <- with(Pairs, rbind(Xi, Xip))
  yi <- Pairs$yi
  segs <- with(Pairs, data.frame(Xi, Xip))[yi == 0,]
  seg.df <- rbind(seg.df, data.frame(norm, segs))
  arrow.df <- with(Pairs,{
    rbind(arrow.df,
          data.frame(norm, Xip, Xi)[yi == -1,],
          data.frame(norm, Xi, Xip)[yi == 1,])
  })
}
seg.df$what <- factor(eq.lab, what.levs)
arrow.df$what <- factor(ineq.lab, what.levs)
train.fun <- factor("training data", fun.levs)
seg.df$fun <- train.fun
arrow.df$fun <- train.fun
library(grid)
basePlot <- ggplot(,aes(X1, X2))+
  facet_grid(.~norm)+
  theme_bw()+
  theme(panel.margin=unit(0,"cm"))+
  coord_equal()
segPlot <- basePlot+
  aes(xend=X1.1,yend=X2.1)+
  geom_segment(data=seg.df)+
  geom_segment(data=arrow.df, arrow=arrow(type="closed",length=unit(0.05,"in")),
               color="red")
print(segPlot)


all.ranks <- simulation$rank
labels <- c(l1="||x||_1^2",
            l2="||x||_2^2",
            linf="||x||_\\infty^2")
all.ranks$label <- sprintf("$r(x) = %s$", labels[as.character(all.ranks$norm)])
seg.df$label <- sprintf("$r(x) = %s$", labels[as.character(seg.df$norm)])
arrow.df$label <- sprintf("$r(x) = %s$", labels[as.character(arrow.df$norm)])
toplot <- data.frame()
plot.funs <- c("rank",
               ##"rank2",
               "compare")
for(fun in plot.funs){
  these <- subset(all.ranks, what %in% c("latent", fun))
  fun <- factor(fun, fun.levs)
  toplot <- rbind(toplot, data.frame(these, fun))
}
toplot$what <- factor(toplot$what, what.levs)
p <- ggplot()+
  geom_segment(aes(X1, X2, xend=X1.1, yend=X2.1, color=what), data=seg.df)+
  geom_segment(aes(X1, X2, xend=X1.1, yend=X2.1, color=what), data=arrow.df,
               arrow=arrow(type="closed",length=unit(0.025,"in")))+
  geom_contour(aes(x1, x2, z=rank, colour=what, group=what),
               data=toplot, size=1)+
  facet_grid(fun~label)+
  theme_bw()+
  theme(panel.margin=unit(0,"cm"))+
  coord_equal()+
  scale_colour_manual("lines",values=model.colors, breaks=what.levs,
                      labels=c(eq.lab, ineq.lab, "latent $r$",
                        "SVMrank\nignore $y_i=0$",
                        "SVMrank\ndouble $y_i=0$",
                        "SVMcompare\nmodel"))+
  xlab("feature 1")+
  ylab("feature 2")+
  guides(colour=guide_legend(keyheight=3))
print(p)

tikz("figure-norm-level-curves.tex", h=5.7)
print(p)
dev.off()

