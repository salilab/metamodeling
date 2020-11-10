These scripts demonstrate the use of bayesian metamodeling of complex biological systems across varying representations.

Authors: Liping Sun, Barak Raveh

License:

Publications:

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

- `data` contains the data of six input models and the metamodel including:  
	- JSON files with the values of model parameters and variables:   
	- `GI.dat` with the observed values for the glucose intake after a meal  
	- `Gb_kt_input_err101.dat` and `Gb_kt_input_sigma101.dat`  with the input values for different accuracy and precision of model
	variables G_B and kt
	- `072919-INS1e-30min-Enrichment-analysis-cleaned-summary.xlsx` with the data for the metabolism model
- `scripts`     contains all the bnet scripts for metamodeling, please refer to READNE.md in `scripts` for more details
- `bnt-master`     contains Bayes Net Toolbox for Matlab with the bugfixes.
