# A Unified Framework for Selecting and Evaluating Cell-Type-Specific Gene Co-expressions in Single-Cell Data
Shan, X., Zhao, H.: A Unified Framework for Selecting and Evaluating Cell-Type-Specific Gene Co-expressions in Single-Cell Data.

## Benchmark of simulation

1. Prepare the 10 datasets used for simulation: Benchmark_simu/prepare_data.R
2. Independent data comparison:
    * Simulate data:   
        * Simulate IND data using scSimu(NB), scDesign2, SPsimSeq, ESCO: Benchmark_simu/simulate_ind.R
        * Simulate IND data using scDesign3: Benchmark_simu/simulate_ind_scDesign3.R
        * Simulate IND data using permutation: Benchmark_simu/simulate_ind_permu.R
    * Estimate data metrics, for example, gene mean, gene var, ...
        * Benchmark_simu/evaluation_ind.R
    * Estimate the KDE test statistics and KS test statistics:
        * Benchmark_simu/evaluate_part2_ind.R
    * Summarize results (the heatmap for benchmark, heatmap for cor mat, density):
        * Benchmark_simu/compare_ind.R

3. Correlated data comparison:
    * Simulate data:
        * Simulate COR data using scSimu(NB), muscat, POWSC, powsimR, scDesign2, SCRIP, Splat, ESCO, SymSim, ZINB-WaVE: Benchmark_simu/simulate_cor.R
        * Simulate COR data using hierarchicell: Benchmark_simu/simulate_cor_hierarchicell.R
        * Simulate COR data using scDesign3: Benchmark_simu/simulate_cor_scDesign3.R
        * Simulate IND/COR data using SPARSim: Benchmark_simu/simulate_SPARSim.R
    * Estimate data metrics, for example, gene mean, gene var, ...
        * Benchmark_simu/evaluation.R
        * Benchmark_simu/evaluation_scDesign3.R
    * Estimate the KDE test statistics and KS test statistics:
        * Benchmark_simu/evaluate_part2.R
        * Benchmark_simu/evaluate_part2_scDesign3.R
    * Summarize results (the heatmap for benchmark):
        * Benchmark_simu/compare_cor.R
    * Computational time for both cor and ind:
        * Benchmark_simu/computational_time.R

## Generate simulated oligodendrocytes from the ROSMAP and PNAS controls
1. Estimate gene margnal parameters
    * PNAS NC Oli: Simulation/marginal_fit_PNAS_NC_Oli.R
    * ROSMAP NC Oli: Simulation/marginal_fit_ROSMAP_NC_Oli.R
    * Ks fit between mu and alpha for ROSMAP/PNAS NC Oli: Simulation/KS_fit_PNAS_and_ROSMAP_NC_Oli.R
2. Estimate gene correlation matrix:
    * Simulation/prepare_NB_simu_PNAS_and_ROSMAP_NC_Oli_PD.R
3. Simulate data that mimic real data:
    * Simulation/NB_simu_PNAS_and_ROSMAP_NC_Oli_PD.R
4. Estimate the marginal parameters of the simulated data for the analysis below:
    * Simulation/marginal_fit_simulated_PNAS_ROSMAP_NC_Oli_PD.R
5. Estimate the empirical p-values of simulated data
    * Estimate the correlation sturctures of simulated real data using seven estimation methods: Simulation/estimate_cor_1000_abs_thresh_PNAS_and_ROSMAP_NC_Oli_PD.R
    * Simulate IND data of the simulated real data and estimate the correlation: Simulation/NB_simu_IND_PNAS_and_ROSMAP_NC_Oli_PD.R
    * Estimate the empirical p-values: Simulation/norm_p_value_PNAS_and_ROSMAP_NC_Oli_PD.R
        *  QQ-plots for the empirical p-values and the comparison with t-dist p-value is generated in this script.

## Comparison between p-value-based approach and correlation-strength-based approach using simulation
1. Prepare known biological networks
    * Prepare STRING: Simulation/Preprocess_STRING.R
    * Prepare Reactome: Simulation/Explore_Reactome.R
2. Compare p-value-based approach and correlation-strength-based approach using CS-CORE
    * Simulation/plot_cscore.R
3. Compare p-value-based approach and correlation-strength-based approach using other estimation methods
    * Reproducibility: Simulation/reproducibility_ROSMAP_and_PNAS.R
    * Overlap with STRING: Simulation/overlap_with_biologNet_STRING.R
    * Overlap with Reactome: Simulation/overlap_with_biologNet_Reactome.R

## Comparison between p-value-based approach and correlation-strength-based approach using real data
1. Excitatory neurons
    * Estimate empirical p-values:
        * Ex marginal fit: Real_data/Ex/prepare_dat.R
        * Est cor in real data: Real_data/Ex/estimate_cor_ROSMAP_NC_ex.R
        * Simulate IND data: Real_data/Ex/NB_simu_IND_ROSMAP_NC_Ex.R
        * Empirical p: Real_data/Ex/norm_p_value_ROSMAP_NC_Ex.R
    * Estimate shared GO terms: Real_data/Ex/identify_cor_ex_v2.R

2. Oligodendrocytes
     * Estimate empirical p-values:
          * The marginal parameters is from the first point of "Generate simulated oligodendrocytes from the ROSMAP and PNAS controls" section
          * Correlated gene name is from the second point of "Generate simulated oligodendrocytes from the ROSMAP and PNAS controls" section
          * Est cor in real data: Real_data/estimate_cor_ROSMAP_NC_Oli.R
          * Simulate IND data: Real_data/NB_simu_IND_ROSMAP_NC_Oli.R
          * Empirical p: Real_data/norm_p_value_ROSMAP_NC_Oli.R
      * Estimate shared GO terms: Real_data/identify_cor_v2.R


## Randomness and bias 
1. Monocytes in the lung data:
    * Estimate empirical p-values:
        * Marginal fit: Randomness_bias/lung_monocyte.R
        * Est cor in real data: Randomness_bias/estimate_cor_lung_NC_monocyte.R
        * Simulate IND data: Randomness_bias/NB_simu_IND_lung_NC_monocyte.R
        * Norm p: Randomness_bias/norm_p_value_lung_NC_monocyte.R
    * Check expression level bias and estimate z-scores:
        * Overall expression level bias and estimate z-scores: Randomness_bias/expression_bias_summarize_lung_NC_monocyte.r
        * Expression level bias on a subset of genes: Randomness_bias/expression_bias_summarize_lung_NC_monocyte_supp.r
2. Oligodendrocytes in the ROSMAP data:
    * Check expression level bias:
        * Expression bias in STRING: Randomness_bias/expression_compared_with_background_v2.R
        * Expression bias in Reactome: Randomness_bias/expression_compared_with_background_reactome_v2.R
    * The correlations and corresponding empirical p-values for Oligodendrocytes in the ROSMAP data have been estimated in "Comparison between p-value-based approach and correlation-strength-based approach using real data".
    * Estimate correlations and corresponding empirical p-values for Oligodendrocytes in the PNAS data (This is used for the reproducibility analysis):
        * The marginal parameters is from the first point of "Generate simulated oligodendrocytes from the ROSMAP and PNAS controls" section
        * Correlated gene name is from the second point of "Generate simulated oligodendrocytes from the ROSMAP and PNAS controls" section
        * Est cor in real data: Randomness_bias/estimate_cor_PNAS_NC_Oli_same_gene.R
        * Simulate IND data: Randomness_bias/NB_simu_IND_PNAS_NC_Oli_same_gene.R
        * Empirical p: Randomness_bias/norm_p_value_PNAS_NC_oli_same_gene.R
    * Estimate z-scores for overlaps with known biological networks and reproducibility:
        * Randomness_bias/expression_bias_summarize_ROSMAP_oli.R

## Others
* Comparing the performance of scSimu and scDesign3 in simulating Smart-Seq data:
    * Simulate using scSimu: Benchmark_simu/smart_seq/scSimu_smart.r
    * Simulate using scDesign3: Benchmark_simu/smart_seq/scdesign3_smart.r
    * Estimate data metrics and test similarity:
        * For simulated correlated data: Benchmark_simu/smart_seq/compare_smart.r
        * For simulated independent data: Benchmark_simu/smart_seq/compare_smart_ind.r
    * Summarize results: Benchmark_simu/smart_seq/plot_comparison.R
* scSimu v.s. simulation in CS-CORE:
    * Benchmark_simu/comparison_scSimu_cscoreSimu.r
* Use simulated data to show the inherent expression bias in estimation methods:
    * Randomness_bias/inherent_methods_bias.R






