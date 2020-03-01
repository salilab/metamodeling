These scripts demonstrate the use of meta-modeling in the modeling of pancreatic beta-cells.
## Prerequisites:

**- matlab:** The scripts are built upon and work with Matlab.  
**- bnt:** Bayes Net Toolbox for Matlab.

## To get started:
1. Download bnt here: 
[https://github.com/bayesnet/bnt](https://github.com/bayesnet/bnt).  

2. Apply two minor bugfixes on bnt:  
(a) Disable line 132 of /Your-Path-To/bnt-master/BNT/learning/learn_params_dbn_em.m

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

(b) wrap line 85 with the following if statement in /Your-Path-To/bnt-master/BNT/general/mk_bnet.m,

```matlabscript
85-  if length(mems)>=1
86-    bnet.rep_of_eclass(e) = mems(1);
87-  end
```

3. Add the path of the bnt package:  
```matlabscript
1- addpath(genpathKPM('/Your-Path-To/bnt-master'))
```

## List of files and directories: 

- `dataset`     contains all relevant input data of individual models
- `bnet`     contains all bnet scripts for meta-modeling, please refer to READNE.md in `bnet` for more details

## Information

_Author(s)_: Barak Raveh, Liping Sun, Tanmoy Sanyal, Kate White, Jeremy Tempkin, etc.

