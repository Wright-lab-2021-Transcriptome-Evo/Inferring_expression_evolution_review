library(phytools)
library(OUwie)
library(dplyr)
library(geiger)

args <- commandArgs(trailingOnly = TRUE)
b <- as.numeric(args[1])
r <- as.numeric(args[2])
fc = as.numeric(args[3])
mtd = (args[4])
print(mtd)

ak_wt <- function(x){
  dak <- x - min(x)
  rl <- sapply(dak, function(y)(exp(-0.5*y)))
  weights <- rl/sum(rl)
  return(weights)
}

data <- data.frame(run = numeric(), b=numeric(), sim = numeric(), BM_AICc=numeric(),
                   OU_AICc=numeric(),  extra_AICc = numeric(), pr=numeric(),
                   BM_wt = numeric(), OU_wt = numeric(),  extra_wt = numeric(),
                   fc = numeric(),
                   lambda = numeric())

c=0
for (i in 1:25){
  c=c+1
  si = (25 *(r-1)) + c
  # ----- #
  v = 0
  sigma = 1
  tree <- sim.bdtree(stop = "taxa", n = b)
  # ----- #
  
  x1 <- fastBM(tree, v, sig2=sigma)
  x2 <- fastBM(tree, v+fc, sig2=sigma)
  Xn <-  (x1 +x2)/2
  bs <- sample(c(1:b), 1)
  Genus_species = tree$tip.label
  Reg = c(rep(1, b))
  Reg[bs] = 2
  names(Reg) <- Genus_species
  tree <- make.simmap(tree, Reg)  
  
  for (pr in seq(0,1,0.1)){
    print(paste("running sim", c, "for pr", pr, "for tree of", b, "branches"))
    Xn[bs] <- (x1[bs]*pr) + (x2[bs]*(1-pr))
    names(Xn) <-  names(x1)
    fitLAM <- fitContinuous(tree, Xn, model = "lambda")
    lb = 1e-999
    if (mtd == "ouwie"){
      print(mtd)
      trait_df <- data.frame(Genus_species, Reg, Xn)
      
      fitBM<-OUwie(tree,trait_df,model="BM1", simmap.tree = TRUE, quiet = FALSE, lb = lb)
      fitOU<-OUwie(tree,trait_df,model="OU1", simmap.tree = TRUE, quiet = TRUE, lb = lb)
      fitBS <- OUwie(tree, trait_df, model = "OUM", simmap.tree = TRUE, quiet = TRUE, lb = lb)
      
      BMAICC <- fitBM$AICc
      OUAICC <- fitOU$AICc
      extraAICC <- fitBS$AICc
      
    } else if (mtd == "geiger"){
      print(mtd)
      fitBM <- fitContinuous(tree, Xn, model = "BM")
      fitOU <- fitContinuous(tree, Xn, model = "OU")
      fitWN <- fitContinuous(tree, Xn, model = "white")
      
      OUAICC <- fitOU$opt$aicc
      BMAICC <- fitBM$opt$aicc
      extraAICC <- fitWN$opt$aicc
    }
    
    
    
    weights <- ak_wt(c(BMAICC, OUAICC, extraAICC))
    
    temp_data <- data.frame(run = r, b=b, sim = c, BM_AICc=BMAICC,
                            OU_AICc=OUAICC, extra_AICc = extraAICC,
                            pr=pr, BM_wt = weights[1], OU_wt = weights[2], extra_wt = weights[3],
                            fc = fc,
                            lambda = fitLAM$opt$lambda)
    
    
    data = rbind(data, temp_data)
  }
}
if (mtd == "geiger"){
  colnames(data)[c(6,10)] <- c("WN_AICc", "WN_wt")
} else {
  colnames(data)[c(6,10)] <- c("BS_AICc", "BS_wt")
}
write.csv(data, file = paste(mtd, "_bs_", b, "_branches_run_", r, "_fc_", fc, ".csv", sep = "" ))








