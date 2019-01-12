# PASS

PASS is a program to calculate Pathway Activation for Single Sample and use them to classify similar complex diseases.
This program includes Matlab scripts and several datasets for demo of PASS approach:

(1）PASS_main.m is the main script to call PASS.

(2) There are some Matlab scripts for each step of PASS analysis, and called in PASS_main.m
    1. CalAUCPathNode.m
    
    2. CalAUCPath.m
    
    3. Ref_network_construction.m
    
    4. Perturbed_network_construction.m
    
(3) The input datasets include:  
    1. Genes_expression_data.mat
    
    2. The folder named Pathway_data contains 50 pathways for demo of PASS.
    
(4) The analysis results are saved as PASSresult.mat

(5) As a demo, users can directly run PASS_main.m in Matlab. This package has been tested in different computer environments as: Window 7 or above; Matlab 2016 or above.
