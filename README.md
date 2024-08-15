# Inferring Cell-Type-Specific Co-Expressed Genes from Single Cell Data

Shan, X., Zhao, H.: Inferring Cell-Type-Specific Co-Expressed Genes from Single Cell Data.

## Benchmark of simulation
* Prepare the 10 dataset used for simulation: Benchmark_simu/prepare_data.R
* Simulate IND data using scSimu(NB), scDesign2, SPsimSeq, ESCO: Benchmark_simu/simulate_ind.R
* Simulate IND data using scDesign3: Benchmark_simu/simulate_ind_scDesign3.R
* Simulate IND data using permutation: Benchmark_simu/simulate_ind_permu.R
* Simulate COR data using scSimu(NB), muscat, POWSC, powsimR, scDesign2, SCRIP, Splat, ESCO, SymSim, ZINB-WaVE: Benchmark_simu/simulate_cor.R
* Simulate COR data using hierarchicell: Benchmark_simu/simulate_cor_hierarchicell.R
* Simulate COR data using scDesign3: Benchmark_simu/simulate_cor_scDesign3.R
* Simulate IND/COR data using SPARSim: Benchmark_simu/simulate_SPARSim.R

* Estimate the data metrics for each simulation: Benchmark_simu/evaluation.R and Benchmark_simu/evaluation_scDesign3.R
* Estimate the KDE test statistics and KS test statistics: Benchmark_simu/evaluate_part2.R and Benchmark_simu/evaluate_part2_scDesign3.R
* Estimate the data metrics for each simulation for IND task: Benchmark_simu/evaluation_ind.R
* Estimate the KDE test statistics and KS test statistics for IND task: Benchmark_simu/evaluate_part2_ind.R

* Summarize the result under correlated cases (the heatmap for benchmark): Benchmark_simu/compare_cor.R
* Summarize the result under independent cases (the heatmap for benchmark, heatmap for cor mat, density): Benchmark_simu/compare_ind.R
* Computational time: Benchmark_simu/computational_time.R

## QQ-plot for the empirical p-values
* Estimate gene marginal parameters in PNAS NC Oli: Simulation/marginal_fit_PNAS_NC_Oli.R
* Estimate gene marginal parameters in ROSMAP NC Oli: Simulation/marginal_fit_ROSMAP_NC_Oli.R
* Ks fit between mu and alpha for ROSMAP/PNAS NC Oli: Simulation/KS_fit_PNAS_and_ROSMAP_NC_Oli.R
* Prepare for simulating real data (marginal para, sct): Simulation/prepare_NB_simu_PNAS_and_ROSMAP_NC_Oli_PD.R
* Simulate real data: Simulation/NB_simu_PNAS_and_ROSMAP_NC_Oli_PD.R
* Estimate the cor of simulated real data: Simulation/estimate_cor_1000_abs_thresh_PNAS_and_ROSMAP_NC_Oli_PD.R
* Estimate marginal para of simulated data: Simulation/marginal_fit_simulated_PNAS_ROSMAP_NC_Oli_PD.R
* Simulate IND data and estimate the correlation: Simulation/NB_simu_IND_PNAS_and_ROSMAP_NC_Oli_PD.R
* Estimate the empirical p-values (qqplot and the comparison with t-dist p): Simulation/norm_p_value_PNAS_and_ROSMAP_NC_Oli_PD.R

## Spurious better performance 
* Prepare of STRING: STRING/Preprocess_STRING.R
* Prepare of Reactome: Reactome/Explore_Reactome.R
* Reproducibility: Simulation/reproducibility_ROSMAP_and_PNAS_1000_abs_thresh_simu_norm_p_PD_v2.R
* Overlap with STRING: Simulation/overlap_with_biologNet_STRING_norm_p_PD.R
* Overlap with Reactome: Simulation/overlap_with_biologNet_Reactome_norm_p_PD.R

## Real data (Oli)
* Gene name from ROSMAP_NC_Oli_cscore_cor1000.rds: Real_data/prepare_NB_simu_PNAS_and_ROSMAP_NC_Oli.R
* Est cor in real data: Real_data/estimate_cor_ROSMAP_NC_Oli.R
* Simulate IND data: Real_data/NB_simu_IND_ROSMAP_NC_Oli.R
* Norm p: Real_data/norm_p_value_ROSMAP_NC_Oli.R
* Share GO terms: Real_data/identify_cor_v2.R

## Real data (Ex) 
* Ex marginal fit: Real_data/Ex/prepare_dat.R
* Est cor in real data: Real_data/Ex/estimate_cor_ROSMAP_NC_ex.R
* Simulate IND data: Real_data/Ex/NB_simu_IND_ROSMAP_NC_Ex.R
* Norm p: Real_data/Ex/norm_p_value_ROSMAP_NC_Ex.R
* Share GO terms: Real_data/Ex/identify_cor_ex_v2.R

## Randomness and bias 
* Randomness_bias/expression_compared_with_background_v2.R
* Randomness_bias/expression_compared_with_background_lung_v2.R
* Randomness_bias/expression_compared_with_background_reactome_v2.R
* Randomness_bias/expression_compared_with_background_reactome_lung_v2.R

