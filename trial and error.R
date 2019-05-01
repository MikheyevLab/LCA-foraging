library(ape)
library(phylobase)
library(tidyverse)
require(dplyr)
library(tidybayes)
library(brms)
library(phytools)
library(phyloseq)
library(ggtree)
library(ggplot2)


speciesdat2 <- read_csv("/Users/Manasee/OneDrive/Documents/Phylocom/redo/completedataR.csv") %>% transmute(Species = gsub(" ", "_", Species), DBH =  DBH, Colony = Colony, Foraged = Foraged, Distance = Distance, Succession = Succession) %>% dplyr::distinct()
head(speciesdat2)
tre2 <- read.nexus("/Users/Manasee/OneDrive/Documents/Phylocom/redo/tree2.nexus")

inv.phylo <- MCMCglmm::inverseA(tre2, nodes = "TIPS", scale = FALSE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)
priorscomplete <- get_prior(Foraged ~ DBH + Distance + (1 | Colony) + (1 | Species), data = speciesdat2, family = bernoulli())
mixedmodelcomplete <- brm(Foraged ~ DBH + Distance + (1 | Colony) + (1 | Species), data = speciesdat2, family = bernoulli, cov_ranef = list(phylo = A), prior = priorscomplete, save_all_pars = TRUE)

mixedmodelcomp <- readRDS("/Users/Manasee/OneDrive/Documents/Phylocom/redo/mixedmodelcompletef.brmsfit")
mixedmodelcomp %>% spread_draws(b_Intercept, r_Species[Species ]) %>% group_by(Species) %>% median_qi(Species_mean =  r_Species)  %>%  ggplot(aes(x = reorder(Species, - Species_mean), y = Species_mean, ymin = .lower, ymax = .upper)) + geom_pointrange() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab("Species") + ylab("species mean") +geom_hline(yintercept = 0, color = "blue") + scale_x_discrete(breaks=NULL)

View(mixedmodelcomp %>% spread_draws(b_Intercept, r_Species[Species ]) %>% group_by(Species) %>% median_qi(Species_mean = b_Intercept + r_Species))
palatcoef <- mixedmodelcomp %>% spread_draws(b_Intercept, r_Species[Species ]) %>% group_by(Species) %>% median_qi(Species_mean =  r_Species)
loo(mixedmodelcomp)
bayes_R2(mixedmodelcomp)



mixedmodelPhyloOnly <- readRDS("/Users/Manasee/OneDrive/Documents/Phylocom/redo/mixedmodelPhyloOnlyf.brmsfit")
summary(mixedmodelPhyloOnly)
plot(mixedmodelPhyloOnly)
hyp <- paste(
  "sd_Species__Intercept^2 /", 
  "(sd_Species__Intercept^2 + 3.141593^(2/3)) = 0"
)

(hyp <- hypothesis(mixedmodelPhyloOnly, hyp, class = NULL))
plot(hyp)
tre3 <- read_tree("/Users/Manasee/OneDrive/Documents/Phylocom/redo/tree.tre")
tiplab <- tre3$tip.label
write.csv(tiplab, file="/Users/Manasee/OneDrive/Documents/Phylocom/redo/tiplab.csv")
tiplab2 <- read.csv("/Users/Manasee/OneDrive/Documents/Phylocom/redo/tiplab.csv")
dat <- read.csv("/Users/Manasee/OneDrive/Documents/Phylocom/redo/dataforphylosig.csv")
head(dat)

meanForaged <- dat %>% select(Species, Foraged) %>% na.omit() %>% group_by(Species)   %>% mutate(f = mean(Foraged))
head(meanForaged)
dat2 <- merge(meanForaged, tiplab2, by = 'Species')
head(dat2)
dat3 <- dat2[,-2]
head(dat3)
dat4 <- distinct(dat3)
head(dat4)
head(meanForaged)
foraged <- as.vector(dat4$f)
names(foraged) <- dat4$Species
phylosigmodel <- phylosig(tre3, foraged, method = "lambda", test = T, se = NULL, start = NULL)
summary(phylosigmodel)
phylosigmodel

foraged1 <- as.vector(meanForaged$f)
names(foraged1) <- meanForaged$Species
phylosigmodel1 <- phylosig(tre3, foraged1, method = "lambda", test = T, se = NULL, start = NULL)
summary(phylosigmodel)
phylosigmodel1

library(ggtree)
library(ggplot2)
tre2 <- read.nexus("/Users/Manasee/OneDrive/Documents/Phylocom/redo/treenew.nexus")
tre2$tip.label<-gsub("_"," ",tre2$tip.label)
wasForaged <- speciesdat2 %>% select(Species, Foraged) %>% na.omit() %>% filter(Species %in% tre2$tip.label) %>% group_by(Species) %>% transmute(foraged = max(Foraged), taxa = Species) %>% unique()
wasForaged2 <- unique(wasForaged)
data <- merge(wasForaged2,palatcoef, by = "Species") %>% unique()
write.csv(data,file = "/Users/Manasee/OneDrive/Documents/Phylocom/redo/datafortree.csv")
data <- read.csv("/Users/Manasee/OneDrive/Documents/Phylocom/redo/datafortree.csv")
data2 <- data  %>% mutate(Species=gsub("_"," ",data$Species))
head(data2)
data
data3 <- data2 %>% transmute(Foraged=as.factor(foraged))
head(data3)
data4 <- cbind(data2,data3)
head(data4)
tr <- ggtree(tre2, layout = "fan")
print(tr)
tr %<+% data2 + text(gsub("_"," ",tre2$tip.label),font = 3)
tree2 <- tr %<+% data4 + geom_tiplab2(aes(color=Foraged))
print(tree2)
ggplo
x <- plotTree(tre2, size = 2, layout = "fan")
plot.ne
print(x)
ggsave("/Users/Manasee/OneDrive/Documents/Phylocom/redo/tree.jpg", width = 50, height = 50, units = "cm", limitsize = FALSE)

dev.new(width=200, height=200, unit="cm", limitsize = FALSE)
tree3 <- ggtree(tre2, layout = "fan") 
tree4 <- tree3 %<+% data4 + geom_tiplab2(size=2, offset=0.25)  + geom_tippoint(size=2.5,aes(color=Species_mean, shape=Foraged, pch=21)) + scale_color_gradient2(high = "red", mid = "white", low = "black", midpoint = 0) + theme(legend.position = "bottom")
tree5 <- tree4 %<+% data2 + geom_tiplab2(aes(color=foraged))
print(tree4)
print(tree4) + geom_treescale(x=20, y=1,color = "black", offset = NULL, fontsize = 0.2, linesize = 0.2)
ggsave("/Users/Manasee/OneDrive/Documents/Phylocom/redo/scaledtree.jpg")
data5 <- data4 %>% select(Species,Species_mean)
head(data5)
as.matrix(data5)
tree5 <- gheatmap(tre2,as.factor(data4$Species_mean),offset = 0, width = 1, color = colorRampPalette(c("red", "white", "blue"))(256), colnames = TRUE, colnames_position = "bottom", colnames_level = NULL, font.size = 4)
tree5 <- heatmap(as.matrix(data5))
print(tree5)
colorRampPalette(c("red", "white", "blue"))(256)
heatmap(as.matrix(data5), scale = "none", col =  col, 
        RowSideColors = rep(c("blue", "pink"), each = 16),
        ColSideColors = c(rep("purple", 5), rep("orange", 6)))


ggsave("/Users/Manasee/OneDrive/Documents/Phylocom/redo/treedot.pdf", width = 200, height = , units = "cm", limitsize = FALSE)
speciesdat3 <- speciesdat2 %>% select(Species,Foraged,Colony) %>% na.omit() %>%filter(Species %in% tre2$tip.label) %>% transmute(Species = Species, Colony = as.factor(Colony), Foraged=Foraged) %>% unique()
head(speciesdat3,20)
speciesdat4 <- speciesdat3 %>% group_by(Species) %>% sample_n(1) %>% unique()
view(speciesdat4)
tre <- read.nexus("/Users/Manasee/OneDrive/Documents/Phylocom/redo/treenew.nexus")
library(phylolm)
fit <- phyloglm(Foraged ~ Colony, phy = tre, data = speciesdat4, boot=100)



set.seed(123456)
tre = rtree(50)
x = rTrait(n=1,phy=tre)
x
X = cbind(rep(1,50),x)
x
y = rbinTrait(n=1,phy=tre, beta=c(-1,0.5), alpha=1 ,X=X)
y
dat = data.frame(trait01 = y, predictor = x)
dat
fit = phyloglm(trait01~predictor,phy=tre,data=dat,boot=100)
summary(fit)
coef(fit)
vcov(fit)