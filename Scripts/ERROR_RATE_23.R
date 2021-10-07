library(OUwie)
library(phytools)
ak_wt <- function(x){
   dak <- x - min(x)
   rl <- sapply(dak, function(y)(exp(-0.5*y)))
   weights <- rl/sum(rl)
   return(weights)
}

#### ERROR RATE ----
data <- data.frame(sim = numeric(), b=numeric(), BM_AICc=numeric(), OU_AICc=numeric(),
  BS_AICc = numeric(), BM_wt = numeric(), OU_wt = numeric(),  BS_wt = numeric())

args <- commandArgs(trailingOnly = TRUE)
r <- as.numeric(args[1])

lb = 1e-100
for (b in c(25, 100)){
   for (i in 1:100){
      print(paste("branch", b, "for simulation", i))
      v = 0
      sigma = 1
      tree <- sim.bdtree(stop = "taxa", n = b)
      Xn <- fastBM(tree, v, sig2=sigma)
      bs <- sample(c(1:b), 1)
      Genus_species = tree$tip.label
      Reg = c(rep(1, b))
      Reg[bs] = 2
      names(Reg) <- Genus_species
      tree <- make.simmap(tree, Reg)
      trait_df <- data.frame(Genus_species, Reg, Xn)
      fitBM<-OUwie(tree,trait_df,model="BM1", simmap.tree = TRUE, quiet = FALSE, lb = lb)
      fitOU<-OUwie(tree,trait_df,model="OU1", simmap.tree = TRUE, quiet = TRUE, lb = lb)
      fitBS <- OUwie(tree, trait_df, model = "OUM", simmap.tree = TRUE, quiet = TRUE, lb = lb)
   
      weights <- ak_wt(c(fitBM$AICc, fitOU$AICc, fitBS$AICc))
      temp_data <- data.frame(sim = i, b=b, BM_AICc=fitBM$AICc, OU_AICc=fitOU$AICc,
                              BS_AICc = fitBS$AICc, BM_wt = weights[1], OU_wt = weights[2], BS_wt = weights[3])
      data = rbind(data, temp_data)
  }
}


write.csv(data, paste(r, "error_rate.csv", sep = "_"), row.names = F )
