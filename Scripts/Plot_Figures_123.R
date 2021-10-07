library(tidyr)
library(dplyr)
library(gridExtra)
library(ggpmisc)
library(ggpubr)
library(ggplot2)
library(ggExtra)
library(plotrix)
library(geiger)
library(phytools)
library(OUwie)
library(rlist)

ak_wt <- function(x){
  dak <- x - min(x)
  rl <- sapply(dak, function(y)(exp(-0.5*y)))
  weights <- rl/sum(rl)
  return(weights)
}


#### ERROR RATE ----
error_data <- do.call(rbind,lapply(Sys.glob("data/final_figure_data/error_rate_S23/*"),read.csv))
error_data$OUtop <- error_data$OU_wt > error_data$BM_wt | error_data$BS_wt > error_data$BM_wt
error_data$OUsig <- (error_data$OU_AICc - error_data$BM_AICc) < -2 | ((error_data$BS_AICc - error_data$BM_AICc) < -2)
test1_error_rate_25 <- mean(error_data$OUtop[error_data$b == 25])
test2_error_rate_25 <- mean(error_data$OUsig[error_data$b == 25])
test1_error_rate_100 <- mean(error_data$OUtop[error_data$b == 100])
test2_error_rate_100 <- mean(error_data$OUsig[error_data$b == 100])

#### Sim 1  ----
s1_error <- read.csv("data/final_figure_data/S1_error.csv")
s1e_25 <- mean(s1_error$Static_exp[data$b == 25])
s1e_100 <- mean(s1_error$Static_exp[data$b == 100])

data_1 <- read.csv("data/final_figure_data/Simulation_1_static_expression_comp_evol.csv")
tidy_data_1 = data_1 %>% 
  pivot_longer(c("BM_wt", "OU_wt", "WN_wt"), names_to = "weight_type",
               values_to = "weight_value")

sum_data_1 = tidy_data_1 %>%
  group_by(b, weight_type) %>%
  summarise(mean = mean(weight_value),
            sd = sd(weight_value)) 

sum_data_1[sum_data_1 == "BM_wt"] <- "BM"
sum_data_1[sum_data_1 == "OU_wt"] <- "OU"
sum_data_1[sum_data_1 == "WN_wt"] <- "WN"
Bar_1 <- ggplot(sum_data_1, aes(x = weight_type, y = mean, fill = as.factor(b))) +
  geom_bar(stat = "identity", position=position_dodge(), color = "black") +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  theme_classic() + xlab("") + ylab("AICc Weight") + 
  scale_fill_brewer(name = "tree size (tips)", palette = "Pastel1") +ylim(0,1)
Bar_1

tiff("Plots/Final_plots/Sim1_bar.tiff", units="in", width=2, height=2, res=400)
Bar_1 + theme(legend.position = "none")
dev.off()

Tree_size_legend <- get_legend(Bar_1)
tiff("Plots/Final_plots/TreeSize_legend.tiff", units="in", width=5, height=5, res=300)
as_ggplot(Tree_size_legend)
dev.off()



#### Sim 2 ----
files_2 = Sys.glob("data/final_figure_data/Simulation_2_branch_shifts/*")
total_data_2 <- do.call(rbind,lapply(files_B,read.csv))
total_data_2$OUtop <- total_data_2$OU_wt > total_data_2$BM_wt & total_data_2$BS_wt > total_data_2$BM_wt
total_data_2$OUsig <- (total_data_2$OU_AICc - total_data_2$BM_AICc) < -2 | (total_data_2$BS_AICc - total_data_2$BM_AICc < -2)

sum_data_2 <- total_data_2 %>% 
  group_by(b, pr) %>% 
  summarise(bm_wt = mean(BM_wt), 
            ou_wt = mean(OU_wt),
            bs_wt = mean(BS_wt),
            test1 = mean(OUtop), 
            test2 = mean(OUsig),
            lambda = mean(lambda),
            n = n())
sum_data_2$test2[sum_data_2$b == 25] <- 
  sum_data_2$test2[sum_data_2$b == 25]-test2_error_rate_25
sum_data_2$test2[sum_data_2$b == 100] <- 
  sum_data_2$test2[sum_data_2$b == 100]-test2_error_rate_100
  
fig_2_line <- ggplot(sum_data_2, aes(x = pr, y = test2, group = b, colour = as.factor(b))) + 
  geom_line() + xlab("Branch specific pr value") + ylab("Relative type 1 error rate") +
  scale_color_discrete(name = "tree size (tips)") +
  theme_classic() #+ 
  #geom_hline(aes(yintercept = test2_error_rate_25, linetype = "Single cell type T1 error rate (25 tips)"), alpha = 0.4) +
  #geom_hline(aes(yintercept = test2_error_rate_100, linetype = "Single cell type T1 error rate (100 tips)"), colour = "black", alpha = 0.5) +
  #scale_linetype_manual(name ="", values = c('dotted','dashed'))

tiff("Plots/Final_plots/Sim2_line.tiff", units="in", width=5, height=5, res=300)
fig_2_line + theme(legend.position = "none")
dev.off()

fig_2_line_legend <- get_legend(fig_2_line)
tiff("Plots/Final_plots/Sim2_line_legend.tiff", units="in", width=5, height=5, res=300)
as_ggplot(fig_2_line_legend)
dev.off()

tidy_data_2 = total_data_2 %>% 
  pivot_longer(c("BM_wt", "OU_wt", "BS_wt"), names_to = "weight_type", 
               values_to = "weight_value")

sum_data_2.2 = tidy_data_2 %>% 
  group_by(b, pr, weight_type) %>% 
  summarise(mean = mean(weight_value), 
            sd = sd(weight_value))
sum_data_2.2[sum_data_2.2 == "BM_wt"] <- "BM"
sum_data_2.2[sum_data_2.2 == "OU_wt"] <- "OU"
sum_data_2.2[sum_data_2.2 == "BS_wt"] <- "BS"
level_order = c("BM", "OU", "BS")

bar_2 <- ggplot(filter(sum_data_2.2, (pr == 0)), aes(x = factor(weight_type, level = level_order), y = mean, fill = as.factor(b))) + 
  geom_bar(stat = "identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  theme_classic() + xlab("") + ylab("AICc Weight")  + 
  scale_fill_brewer(palette = "Pastel1") +ylim(0,1)

tiff("Plots/Final_plots/Sim2_bar.tiff", units="in", width=2, height=2, res=400)
bar_2 + theme(legend.position = "none")
dev.off()

Sim_2_25_T1 <- mean(filter(sum_data_2, (b == 25 & (pr == 0 | pr ==1)))$T1_error) - test2_error_rate_25
Sim_2_100_T1 <- mean(filter(sum_data_2, (b == 100 & (pr == 0 | pr ==1)))$T1_error) - test2_error_rate_100



#### Sim 3 FIGURE ----
files_3T = Sys.glob("data/final_figure_data/Simulation_3_Two_traits_cv/*TRUE*")
total_data3T <- do.call(rbind,lapply(files_CT,read.csv))
total_data3T$test = ", composition evolving"
files3F = Sys.glob("data/final_figure_data/Simulation_3_Two_traits_cv/*FALSE*")
total_data3F <- do.call(rbind,lapply(files3F,read.csv))
total_data3F$test = ", composition not evolving"

total_data3 <- rbind(total_data3F, total_data3T)

total_data3$OUtop <- total_data3$OU_wt > total_data3$BM_wt & total_data3$BS_wt > total_data3$BM_wt
total_data3$OUsig <- (total_data3$OU_AICc - total_data3$BM_AICc) < -2 | (total_data3$BS_AICc - total_data3$BM_AICc < -2)


sum_data3 <- total_data3 %>% 
  group_by(b, cv, test) %>% 
  dplyr::summarise(bm_wt = mean(BM_wt), 
            ou_wt = mean(OU_wt), 
            test1 = mean(OUtop), 
            test2 = mean(OUsig),
            lambda = mean(lambda),
            n = n(), T1_error = mean(OUsig))


sum_data3$test2[sum_data3$b == 25] <- 
  sum_data3$test2[sum_data3$b == 25]-test2_error_rate_25
sum_data3$test2[sum_data3$b == 100] <- 
  sum_data3$test2[sum_data3$b == 100]-test2_error_rate_100



fig_3_line <- ggplot(sum_data3, aes(x = cv, y = test2, colour = paste(b,test, sep = ""), 
                      group = paste(b,test, sep = ""))) + geom_line() +
  theme_classic() + scale_color_manual(values=c("turquoise", "blue", "green", "red"), 
                     name = "tree size (tips) & cellular composition") +
  xlab("Gene expression covariance") + ylab("Relative type 1 error rate") #+
  #geom_hline(aes(yintercept = test2_error_rate_25, linetype = "Single cell type T1 error rate (25 tips)"), alpha = 0.4) +
  #geom_hline(aes(yintercept = test2_error_rate_100, linetype = "Single cell type T1 error rate (100 tips)"), colour = "black", alpha = 0.5) +
  #scale_linetype_manual(name ="", values = c('dotted','dashed')) 

tiff("Plots/Final_plots/Sim3_line.tiff", units="in", width=5, height=5, res=300)
fig_3_line + theme(legend.position = "none")
dev.off()

fig_3_line_legend <- get_legend(fig_3_line)
tiff("Plots/Final_plots/Sim3_line_legend.tiff", units="in", width=5, height=5, res=300)
as_ggplot(fig_3_line_legend)
dev.off()




tidy_data_3 = filter(total_data3, test == ', composition evolving') %>% 
  pivot_longer(c("BM_wt", "OU_wt", "BS_wt"), names_to = "weight_type", 
               values_to = "weight_value")

sum_data_3.2 = tidy_data_3 %>% 
  group_by(b, cv, weight_type) %>% 
  summarise(mean = mean(weight_value), 
            sd = sd(weight_value))
sum_data_3.2[sum_data_3.2 == "BM_wt"] <- "BM"
sum_data_3.2[sum_data_3.2 == "OU_wt"] <- "OU"
sum_data_3.2[sum_data_3.2 == "BS_wt"] <- "BS"

bar_3 <- ggplot(filter(sum_data_3.2, (cv == 0)), aes(x = factor(weight_type, level = level_order), y = mean, fill = as.factor(b))) + 
  geom_bar(stat = "identity", position=position_dodge(), color = "black") + 
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  theme_classic() + xlab("Model") + ylab("AICc Weight") + 
  scale_fill_brewer(palette = "Pastel1")+ylim(0,1)

tiff("Plots/Final_plots/Sim3_bar.tiff", units="in", width=2, height=2, res=400)
bar_3 + theme(legend.position = "none")
dev.off()

Sim_3_25_T1 <- mean(filter(sum_data3, (b == 25 & cv == 0  & test == ", composition evolving"))$T1_error) - test2_error_rate_25
Sim_3_100_T1 <- mean(filter(sum_data3, (b == 100 & cv == 0  & test == ", composition evolving"))$T1_error) - test2_error_rate_100





