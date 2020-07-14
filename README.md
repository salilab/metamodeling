These scripts demonstrate the use of bayesian meta-modeling in the modeling of pancreatic beta-cells.
## Prerequisites:

- **matlab:** The scripts are built upon and work with Matlab.  
- **bnt:** Bayes Net Toolbox for Matlab.

## To get started:
- Download bnt here: 
[https://github.com/bayesnet/bnt](https://github.com/bayesnet/bnt).  

- Apply two minor bugfixes on bnt:  
	- Disable line 132 of /Your-Path-To/bnt-master/BNT/learning/learn_params_dbn_em.m. 

	```matlabscript
	125- loglik = 0;
	126- for l=1:length(cases)
	127-   evidence = cases{l};
	128-   if ~iscell(evidence)
	129-     error('training data must be a cell array of cell arrays')
	130-   end
	131-   [engine, ll] = enter_evidence(engine, evidence);
	132-   % assert(~isnan(ll))
	133-   loglik = loglik + ll;
	134-   T = size(evidence, 2);
	```

	- Wrap line 85 with the following if statement in /Your-Path-To/bnt-master/BNT/general/mk_bnet.m. 

	```matlabscript
	85-  if length(mems)>=1
	86-    bnet.rep_of_eclass(e) = mems(1);
	87-  end
	```

- Add the path of the bnt package:  
```matlabscript
1- addpath(genpathKPM('/Your-Path-To/bnt-master'))
```

## List of files and directories: 

- `data`     contains the data of individual models including:  
	- .json files with the input values and parameters of model variables   
	- DGintake_dt1_sigmoid.dat with the observed values for glucose intake after a meal  
	- Gb_k_input_errors21.dat  with the input values for different accuracy and precision of model
	variables k_T and G_B  
- `scripts`     contains all bnet scripts for meta-modeling, please refer to READNE.md in `scripts` for more details
- `bnt-master`     contains Bayes Net Toolbox for Matlab after bugfixes.

## Information

_Author(s)_: Barak Raveh, Liping Sun, Tanmoy Sanyal, Kate White, Jeremy Tempkin, etc.

