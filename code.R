library(ape)
library(phytools)
library(tidyverse)
library(brms)

speciesdat2 <- read_csv("speciesdist.csv", col_types = "-nnfcccn") %>% transmute(Species = gsub(" ", "_", Species), DBH =  DBH, Colony = Colony, Foraged = Foraged, Distance = Distance) %>% dplyr::distinct()
tre2 <- read.nexus("tre.nex")
inv.phylo <- MCMCglmm::inverseA(tre2, nodes = "TIPS", scale = FALSE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)
priors <- get_prior(Foraged ~ DBH + Distance + (1 | Colony) + (1 | Species), data = speciesdat2, family = bernoulli())
mixedmodel2 <- brm(Foraged ~ DBH + Distance + (1 | Colony) + (1 | Species), data = speciesdat2, family = bernoulli, cov_ranef = list(phylo = A), prior = priors, save_all_pars = TRUE, seed = 123, cores = 4)
summary(mixedmodel2)
mixedmodelL <- loo(mixedmodel2)
plot(mixedmodel2)
bayes_R2(mixedmodel2)


priors2 <- get_prior(Foraged ~ DBH + Distance + (1 | Colony) , data = speciesdat2, family = bernoulli())
mixedmodel.noPhylo <- brm(Foraged ~ DBH + Distance + (1 | Colony) , data = speciesdat2, family = bernoulli,  prior = priors2, save_all_pars = TRUE, seed = 123, cores = 4)
mixedmodel.noPhylo.LOO <- loo(mixedmodel.noPhylo)
loo_compare(mixedmodelL, mixedmodel.noPhylo.LOO)
bayes_R2(mixedmodel.noPhylo)

#Since higher LOOIC values indicate better fit, the phylogenetic model fits better
