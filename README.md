## A Bayesian approach for optimizing metamodeling under uncertainty: Application to complex biological systems

### Summary
Bayesian metamodeling divides and conquers this task by integrating a collection of input models. It proceeds through three stages: (i) convert input models into probabilistic surrogate models; (ii) couple surrogate models; and (iii) update surrogate models and input models via backpropagation. Comprehensive understanding and quantitative assessment of uncertainty in Bayesian metamodeling are crucial for its proper interpretation, although challenging due to its inherent high-complexity. We analyze the propagation of uncertainty across metamodeling stages using both an analytical example of two Gaussian time-independent input models and a numerical example of three time-dependent input models describing glucose-stimulated insulin secretion in pancreatic $\beta$-cells. We elucidate criteria in the selection of surrogate models and coupling strategies, as well as identification and resolution of potential conflicts among the input models. The optimized metamodel demonstrates better alignment with experimental data compared to the non-optimized one. 

### List of files
- ``metamodel`` contains the scripts for the construction of a glucose-stimulated-insulin-secretion (GSIS) metamodel with three input models, using the BNET package in MATLAB by Kevin Murphy [http:// github.com/bayesnet/bnt]{http:// github.com/bayesnet/bnt}. 
- ``analysis`` contains the scripts for metamodel analysis and plots for the figures in main text.


### Information
Author: Chenxi Wang
