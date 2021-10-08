## Description of scripts for simulations and plots 

 **[Simulation_1_no_exp_evo_comp_evo.R](https://github.com/Wright-lab-2021-Transcriptome-Evo/Inferring_expression_evolution_review/blob/main/scripts/Simulation_1_no_exp_evo_comp_evo.R "Simulation_1_no_exp_evo_comp_evo.R")**

- Simulates scenario i where expression is static, however abundance evolves under Brownian Motion


**[Simulation_2_bs_exp_evo.R](https://github.com/Wright-lab-2021-Transcriptome-Evo/Inferring_expression_evolution_review/blob/main/scripts/Simulation_2_bs_exp_evo.R "Simulation_2_bs_exp_evo.R")**

- Simulates scenario ii where expression evolves under Brownian motion and abundance is 0.5 bar one tip, where it fluctuates (0 to 1)

**[Simulation_3_expandcomp_evo.R](https://github.com/Wright-lab-2021-Transcriptome-Evo/Inferring_expression_evolution_review/blob/main/scripts/Simulation_3_expandcomp_evo.R "Simulation_3_expandcomp_evo.R")**

- Simulates scenario iii where expression evoles under BM with varying covariances, and abundance is either 0.5 or is evolving under BM. 

**[ERROR_RATE_1.R](https://github.com/Wright-lab-2021-Transcriptome-Evo/Inferring_expression_evolution_review/blob/main/scripts/ERROR_RATE_1.R "ERROR_RATE_1.R")**

- Script to calculate the type 1 error rate when running the same set of models as in scenario i but on a non-composite dataset (i.e. a single cell type). Type 1 error rate in this case is the percentage of cases where sigma<sup>2</sup> is not identified as 0 in a boostrapped BM model 

**[ERROR_RATE_23.R](https://github.com/Wright-lab-2021-Transcriptome-Evo/Inferring_expression_evolution_review/blob/main/scripts/ERROR_RATE_23.R "ERROR_RATE_23.R")**

- Script to calcualte the type 1 error rate when running the same set of models as in scenario ii and iii but on a non-composite dataset (i.e a single cell type). Type 1 error rate in this case is the percentage of cases when the AICc value for an alternate model (OU or OU with a shift) is at least 2 less than that for a Brownian motion model. 

**[Plot_Figures_123.R](https://github.com/Wright-lab-2021-Transcriptome-Evo/Inferring_expression_evolution_review/blob/main/scripts/Plot_Figures_123.R "Plot_Figures_123.R")**

- Script to calculate relative error rates in each of the scenarios and reproduce figures used in the manuscript. 

