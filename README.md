# ImputationBound: 
code for paper xxx

compile: g++ -std=c++11

"Compute_case_I(II)" computes the conditional probability of having completely wrong or ambiguous imputed dosage given reference panel allele count, reference size and population growth model (can be specified by exprected time between coalescent events). 

"Simulation" performs coalescent simulations to generate both the above probabilities and the upper bound of average $r^2$ between imputed dosage and true genotypes. 
