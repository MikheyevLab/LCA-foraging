library(tidyverse)
require(dplyr)
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

tx1 <- read.csv("taxa_unique.csv")
unique(tx1)
tx2 <- transmute(tx1,Species=as.character(Species_f), Family=as.character(Family_f))
tx3 <- cbind(tx1,tx2)
tx4 <- tx3[,-(1:2)]
taxa_list <- c()
for (i in 1:nrow(tx4)) {
  family_name <- tolower(tx4$Family[i]) # converts Family name to lowercase
  split_binomial_name <- strsplit(tx4$Species[i], split=' ')
  genus_name <- tolower(split_binomial_name[[1]][1])
  species_name <- paste(split_binomial_name[[1]][1], split_binomial_name[[1]][2], sep='_')
  phylomatic_syntax <- paste(family_name, genus_name, species_name, sep='/')
  taxa_list <- c(taxa_list, phylomatic_syntax)
}
taxa <- unique(taxa_list)
taxalist <- list(taxa)
head(taxalist)
lapply(taxalist, write, "taxalistnew.txt", append=TRUE) #Input for Phylomatic

#Input for Sango
speciesdat2 <- read_csv("completedataR.csv") %>% transmute(Species = gsub(" ", "_", Species), DBH =  DBH, Colony = Colony, Foraged = Foraged, Distance = Distance, Succession = Succession) %>% dplyr::distinct()
tre2 <- read.nexus("treenew.nexus") #Output from Phylomatic and FigTree

inv.phylo <- MCMCglmm::inverseA(tre2, nodes = "TIPS", scale = FALSE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)
priorscomplete <- get_prior(Foraged ~ DBH + Distance + (1 | Colony) + (1 | Species), data = speciesdat2, family = bernoulli())
mixedmodelcompletef <- brm(Foraged ~ DBH + Distance + (1 | Colony) + (1 | Species), data = speciesdat2, family = bernoulli, cov_ranef = list(phylo = A), prior = priorscomplete, save_all_pars = TRUE)
saveRDS("mixedmodelcompletef", file = "mixedmodelcompletef.brmsfit")
priors.noPhylo <- get_prior(Foraged ~ DBH + Distance + (1 | Colony), data = speciesdat2, family = bernoulli())
mixedmodel.noPhylo <- brm(Foraged ~ DBH + Distance + (1 | Colony), data = speciesdat2, family = bernoulli, cov_ranef = list(phylo = A), prior = priors.noPhylo, save_all_pars = TRUE)
saveRDS("mixedmodel.noPhylo", file = "mixedmodel_noPhylof.brmsfit")
priorsPhyloOnly <- get_prior(Foraged ~ (1 | Species), data = speciesdat2, family = bernoulli())
mixedmodelPhyloOnly <- brm(Foraged ~ (1 | Species), data = speciesdat2, family = bernoulli, cov_ranef = list(phylo = A), prior = priorsPhyloOnly, save_all_pars = TRUE)
saveRDS("mixedmodelPhyloOnly", file = "mixedmodelPhyloOnlyf.brmsfit")

#Outputfrom Sango
mixedmodelcomp <- readRDS("mixedmodelcompletef.brmsfit")
summary(mixedmodelcomp)
plot(mixedmodeldcomp)
mixedmodelcomp %>% spread_draws(b_Intercept, r_Species[Species ]) %>% group_by(Species) %>% median_qi(Species_mean =  r_Species)  %>%  ggplot(aes(x = reorder(Species, - Species_mean), y = Species_mean, ymin = .lower, ymax = .upper)) + geom_pointrange() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab("Species") + ylab("species mean") +geom_hline(yintercept = 0, color = "blue") + scale_x_discrete(breaks=NULL)

View(mixedmodelcomp %>% spread_draws(b_Intercept, r_Species[Species ]) %>% group_by(Species) %>% median_qi(Species_mean = b_Intercept + r_Species))
palatcoef <- mixedmodelcomp %>% spread_draws(b_Intercept, r_Species[Species ]) %>% group_by(Species) %>% median_qi(Species_mean =  r_Species)
mixedmodelcompLOO <- loo(mixedmodelcomp)
bayes_R2(mixedmodelcomp)


mixedmodel.noPhylo <- readRDS("mixedmodel_noPhylof.brmsfit")
summary(mixedmoded.noPhylo)
noPhyloLoo <- loo(mixedmodel.noPhylo)
bayes_R2(mixedmodel.noPhylo)



View(mixedmodelcomp %>% spread_draws(b_Intercept, r_Species[Species ]) %>% group_by(Species) %>% median_qi(Species_mean = b_Intercept + r_Species))
palatcoef <- mixedmodelcomp %>% spread_draws(b_Intercept, r_Species[Species ]) %>% group_by(Species) %>% median_qi(Species_mean =  r_Species)

mixedmodelPhyloOnly <- readRDS("mixedmodelPhyloOnlyf.brmsfit")
summary(mixedmodelPhyloOnly)
plot(mixedmodelPhyloOnly)
hyp <- paste("sd_Species__Intercept^2 /", "(sd_Species__Intercept^2 + 3.141593^(2/3)) = 0")
(hyp <- hypothesis(mixedmodelPhyloOnly, hyp, class = NULL))
plot(hyp)

tre3 <- read_tree("tree.tre")
tiplab <- tre3$tip.label
write.csv(tiplab, file="tiplab.csv")
tiplab2 <- read.csv("tiplab.csv")
dat <- read.csv("dataforphylosig.csv")
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

tre2$tip.label<-gsub("_"," ",tre2$tip.label)
wasForaged <- speciesdat2 %>% select(Species, Foraged) %>% na.omit() %>% filter(Species %in% tre2$tip.label) %>% group_by(Species) %>% transmute(foraged = max(Foraged), taxa = Species) %>% unique()
wasForaged2 <- unique(wasForaged)
data <- merge(wasForaged2,palatcoef, by = "Species") %>% unique()
write.csv(data,file = "datafortree.csv")
data <- read.csv("datafortree.csv")
data2 <- data  %>% mutate(Species=gsub("_"," ",data$Species))
data3 <- data2 %>% transmute(Foraged=as.factor(foraged))
data4 <- cbind(data2,data3)
tree3 <- ggtree(tre2, layout = "fan") 
tree4 <- tree3 %<+% data4 + geom_tiplab2(size=2, offset=0.25)  + geom_tippoint(size=2.5,aes(color=Species_mean, shape=Foraged, pch=21)) + scale_color_gradient2(high = "red", mid = "white", low = "black", midpoint = 0) + theme(legend.position = "bottom")
tree5 <- tree4 %<+% data2 + geom_tiplab2(aes(color=foraged))
print(tree4)
print(tree4) + geom_treescale(x=20, y=1,color = "black", offset = NULL, fontsize = 0.2, linesize = 0.2)
ggsave("scaledtree.jpg")



