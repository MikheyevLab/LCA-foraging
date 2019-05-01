t1 <- read.nexus("treenew.nexus") %>% mutate(tip.label=gsub("_"," ",tip.label))

df6 <- read.csv("dataforplot.csv")
df12 <- data.frame(label=df6$Species, Foraging_preference=df6$foraged, Palatability_coefficient=df6$Species_mean)
row.names(df12) <- df12[,1]
view(df12)

p3 <- ggtree(t1a, layout = 'fan') %<+% df12 + geom_tiplab2(size=2,allign=TRUE, offset=0.5,aes(color=as.factor(Foraging_preference))) + scale_color_manual(values = c("black", "red"))
plot(p3)
p4 <- gheatmap(p3, df12[,"Palatability_coefficient", drop=FALSE], offset = -1, width=0.05, colnames=FALSE, colnames_position = "bottom") + scale_fill_gradient2(low = "black", mid = "white", high = "red", midpoint = 0) + theme(legend.position = "none")
plot(p4)