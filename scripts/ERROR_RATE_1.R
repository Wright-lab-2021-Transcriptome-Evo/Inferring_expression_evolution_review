library(geiger)
library(phylolm)
library(dplyr)
library(MuMIn)
library(phytools)

ak_wt <- function(x){
   dak <- x - min(x)
   rl <- sapply(dak, function(y)(exp(-0.5*y)))
   weights <- rl/sum(rl)
   return(weights)
}
data <- data.frame(exp_value = as.numeric(), b = as.numeric(), sim= numeric(),
                   Static_exp = numeric())

args <- commandArgs(trailingOnly = TRUE)
r <- as.numeric(args[1])

c=0
for (b in c(25, 100)){
   for (r in 1:50){
      tree <- sim.bdtree(stop = "taxa", n = b)
      c=c+1
      for (e in c(1,2)){
         print(paste("running simulation ", c, "for tree of ", b, "branches"))
         Xn <- (rep(e, b))
         names(Xn) <- tree$tip.label 
         Xn <- data.frame(x = Xn)
         
         # to test whether static model is rejected...
         # use bootstraps to estimate 95% CIs on sigma2 estimate - if does not overlap 0 then static model rejected
         Static <- summary(phylolm(x ~ 1, Xn, tree, model="BM", boot = 1000))
         if (range(Static$bootconfint95[,2])[1]<=0){
            static_conclusion = TRUE
         } else {
            static_conclusion = FALSE
         }
         
         temp_data <- data.frame(exp_value = e, b = b, sim = c, Static_exp = static_conclusion)
         
         data = rbind(data, temp_data)
      }
   }
}


write.csv(data, paste("data/", r, "_S1_error_rate.csv", sep = ""), row.names = F )



