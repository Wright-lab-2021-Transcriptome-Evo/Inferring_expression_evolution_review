library(phytools)
library(OUwie)

ak_wt <- function(x){
  dak <- x - min(x)
  rl <- sapply(dak, function(y)(exp(-0.5*y)))
  weights <- rl/sum(rl)
  return(weights)
}

args <- commandArgs(trailingOnly = TRUE)
b <- as.numeric(args[1])
r <- as.numeric(args[2])
comp_evol = args[3]
mtd = args[4]

print(mtd)
print(paste("composition is ", comp_evol, " evolving"))

fc = 2
sigma = 1.0001



data <- data.frame(run = numeric(), b=numeric(), sim = numeric(), BM_AICc=numeric(),
                   OU_AICc=numeric(),  extra_AICc = numeric(), cv=numeric(),
                   BM_wt = numeric(), OU_wt = numeric(),  extra_wt = numeric(),
                   fc = numeric(),
                   lambda = numeric())


c=0
for (cv in rep(seq(-1,1,0.05), 20)){
  c=c+1
  si = (820 *(r-1 )) + c
  # ----- #
  tree <- sim.bdtree(stop = "taxa", n = b)
  print(paste("running cv ", cv, " for repeat of ", c, sep = ""))
  xt<-fastBM(tree, a = 0.5, sig2=sigma,internal=FALSE, nsim = 1)
  y = (xt-min(xt))/(max(xt) - min(xt))
  
  xe<-sim.corrs(tree, vcv=matrix(c(sigma, cv, cv, sigma),2,2), anc = c(4,(4+fc)))
  
  X1 <- xe[,1]
  X2 <- xe[,2]
  
  if (comp_evol == TRUE){
    Xn <- (X1*y) + (X2*(1-y))
  } else {
    Xn <- (X1 + X2)/2
  }
  bs <- sample(c(1:b), 1)
  Genus_species = tree$tip.label
  Reg = c(rep(1, b))
  Reg[bs] = 2
  names(Reg) <- Genus_species
  tree <- make.simmap(tree, Reg)  
  
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
                          cv=cv, BM_wt = weights[1], OU_wt = weights[2], extra_wt = weights[3],
                          fc = fc,
                          lambda = fitLAM$opt$lambda)
  
  data <- rbind(data, temp_data)
}

if (mtd == "geiger"){
  colnames(data)[c(6,10)] <- c("WN_AICc", "WN_wt")
} else {
  colnames(data)[c(6,10)] <- c("BS_AICc", "BS_wt")
}

sav_name3 <- paste("Two_traits_cv_evo_comp_evo_branches_", b, "_rep_", r, 
                   "_fc_", fc, "_compevol_", comp_evol, ".csv", sep = "")
write.csv(data, sav_name3, row.names = FALSE)

