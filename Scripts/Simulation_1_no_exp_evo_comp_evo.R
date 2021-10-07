install.packages("geiger")
install.packages("phylolm")

library(geiger)
library(phylolm)
library(dplyr)
library(MuMIn)
library(phytools)
fc = 1

ak_wt <- function(x){
  dak <- x - min(x)
  rl <- sapply(dak, function(y)(exp(-0.5*y)))
  weights <- rl/sum(rl)
  return(weights)
}
data <- data.frame(fc = as.numeric(), b = as.numeric(), sim= numeric(),
                   BM_AICc_exp = numeric(), OU_AICc_exp = numeric(),
                   WN_exp = numeric(), Static_exp = numeric(), 
                   BM_wt = numeric(), OU_wt = numeric(), WN_wt = numeric())

c=0
for (b in c(25, 100)){
  for (r in 1:1000){
    print(paste("running simulation ", c, "for tree of ", b, "branches"))
    c=c+1
    tree <- sim.bdtree(stop = "taxa", n = b)
    x<-fastBM(tree, a = 0, sig2=1)
    y = (x-min(x))/(max(x) - min(x))
    X1 <- rep(1, b)
    X2 <- X1 + fc
    Xn <- (X1*y) + (X2*(1-y))
    Xn <- data.frame(x = Xn)
    
    BM <- (phylolm(x ~ 1, Xn, tree, model="BM"))
    
    OU <- (phylolm(x ~ 1, Xn, tree, model="OUfixedRoot"))

    fix <- data.frame(lambda=0)
    WN <- (phylolm(x ~ 1, Xn, tree, model="lambda", starting.value = fix, lower.bound = fix$lambda, upper.bound = fix$lambda))
    
    # to test whether static model is rejected...
    # use bootstraps to estimate 95% CIs on sigma2 estimate - if does not overlap 0 then static model rejected
    Static <- summary(phylolm(x ~ 1, Xn, tree, model="BM", boot = 1000, full.matrix = TRUE))
    if (range(Static$bootconfint95[,2])[1]<=0){
      static_conclusion = TRUE
    } else {
      static_conclusion = FALSE
    }
    
    weights = ak_wt(c(AICc(BM), AICc(OU)))
    
    temp_data <- data.frame(fc = fc, b = b, sim = c, BM_AICc_exp = AICc(BM), 
                            OU_AICc_exp = AICc(OU), WN_exp = AICc(WN), 
                            Static_exp = static_conclusion, BM_wt = weights[1], 
                            OU_wt = weights[2], WN_wt = weights[3])
    data = rbind(data, temp_data)

  }
}


###write.csv(data, file = "data/final_figure_data/Simulation_1_static_expression_comp_evol.csv", row.names = F)
##
