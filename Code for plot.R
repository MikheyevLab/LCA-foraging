t1 <- read.nexus("treenew.nexus") %>% mutate(tip.label=gsub("_"," ",tip.label))

plotdat1 <- read.csv("dataforplot.csv")
plotdat2 <- data.frame(label=plotdat1$Species, Foraging_preference=plotdat1$foraged, Palatability_coefficient=plotdat1$Species_mean)
row.names(plotdat2) <- plotdat2[,1]

p1 <- ggtree(t1, layout = 'fan') %<+% plotdat2 + geom_tiplab2(size=2,allign=TRUE, offset=0.5,aes(color=as.factor(Foraging_preference))) + scale_color_manual(values = c("black", "red"))
p2 <- gheatmap(p1, plotdat2[,"Palatability_coefficient", drop=FALSE], offset = -1, width=0.05, colnames=FALSE, colnames_position = "bottom") + scale_fill_gradient2(low = "black", mid = "white", high = "red", midpoint = 0) + theme(legend.position = "none")

