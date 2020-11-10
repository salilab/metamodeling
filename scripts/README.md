Created by: 
Liping Sun, Barak Raveh

Date created: 
1-7-19

Synopsis:
This repository is for developing Bayesian metamodeling implementations. 

## Add the bnet toolbox
`add_bnet.m`

## Factories and Functions:
`CPDFactory.m`  
`DBNFactory.m`  
`get_eclass_from_maps.m`  
`get_reverse_nodes_map.m`  
`get_valid_nodes_graph.m`  
`merge_dbn_factories.m`  

## To make the dbn factories for six input models and the metamodel:
`make_postprandial_dbn.m`  
`make_pancreas_dbn.m`  
`make_signaling_dbn.m`   
`make_exocytosis_dbn.m`   
`make_metabolism_dbn.m`  
`make_screening_dbn.m`   
`make_meta_dbn6.m`  

## To run the model priors:
`postprandial_prior_normal.m`  
`postprandial_prior_t2d.m`   
`pancreas_prior.m`   
`network_prior.m`   
`spt_prior.m`  
`screening_prior.m`  
`metabolism_prior.m`   

## To run the metamodel for the incretin effects (Fig. 3):
`metamodel_normal_GLP1.m`  - metamodel of nomal subjects after a meal,
with different GLP1 concentrations.   
`metamodel_normal_incretin.m` -metamodel of nomal subjects after a
meal, with different incretin concentrations.   
`metamodel_t2d_GLP1.m`  - metamodel of t2d subjects after a meal,
with different GLP1 concentrations.   
`metamodel_t2d_incretin.m` -metamodel of t2d subjects after a meal,
with different incretin concentrations.   

## To run the metamodel for the accuracy and precision of model variables (Fig. 4):
`metamodel_normal_kt_mean51_cov32_Gb_mean51.m` - an example script for
the metamodel of a nomal subject to sample output kt with different
input sigma of Gb.   
