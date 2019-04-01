library(ape)
library(phytools)
library(tidyverse)
library(brms)
library(tidybayes)

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
plot(marginal_effects(mixedmodel2), points = TRUE)
bayes_R2(mixedmodel2)

# this doesn't work because no sigma in bernoulli
hyp <- "sd_Species__Intercept^2 / (sd_Species__Intercept^2 + sd_Colony__Intercept^2 + Distance__Intercept^2 ) = 0"
(hyp <- hypothesis(mixedmodel2, hyp, class = NULL))
plot(hyp)
# https://groups.google.com/forum/#!msg/brms-users/8ADbWh7v9kc/-_wHkmLMCAAJ
priors2 <- get_prior(Foraged ~ DBH + Distance + (1 | Colony) , data = speciesdat2, family = bernoulli())
mixedmodel.noPhylo <- brm(Foraged ~ DBH + Distance + (1 | Colony) , data = speciesdat2, family = bernoulli,  prior = priors2, save_all_pars = TRUE, seed = 123, cores = 4)
mixedmodel.noPhylo.LOO <- loo(mixedmodel.noPhylo)
loo_compare(mixedmodelL, mixedmodel.noPhylo.LOO)
bayes_R2(mixedmodel.noPhylo)


mixedmodel2 %>% spread_draws(b_Intercept, r_Species[Species ]) %>% group_by(Species) %>% median_qi(Species_mean =  r_Species)  %>%  ggplot(aes(x = reorder(Species, - Species_mean), y = Species_mean, ymin = .lower, ymax = .upper)) + geom_pointrange() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlab("Species") + ylab("species mean") +geom_hline(yintercept = 0, color = "blue") + scale_x_discrete(breaks=NULL)

View(mixedmodel2 %>% spread_draws(b_Intercept, r_Species[Species ]) %>% group_by(Species) %>% median_qi(Species_mean = b_Intercept + r_Species))
