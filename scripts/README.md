Authors in random order: 
Jeremy Tempkin, Liping Sun, Barak Raveh

Date created: 
1-7-19

Synopsis:
This repository is for developing meta-modeling implementations. 

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
`make_GLP1R_dbn.m` 
`make_meta_dbn6.m`

## To run individual model priors:
`postprandial_prior_normal.m`
`meal_model_prior_t2d.m` 
`pancreas_model_prior_normal.m` 
`network_model_prior_normal.m` 
`spt_model_prior_normal.m`
`GLP1R_model_prior_normal.m`
`KEGG_model_prior_normal.m` 

## To run the metamodel for the incretin effects:
`metamodel_normal_GLP1.m`  - metamodel of nomal subjects after a meal, with GLP1 levels. 
`metamodel_normal_incretin.m` -metamodel of nomal subjects after a meal, with incretin levels. 
`metamodel_normal_GLP1.m`  - metamodel of t2d subjects after a meal, with GLP1 levels. 
`metamodel_normal_incretin.m` -metamodel of t2d subjects after a meal, with incretin levels. 
`run.py` -python script to run the individual models and the metamodel above using matlab python engine.

##To run the metamodel for the accuracy and precision of model variables:
`metamodel_normal_k_mean11_cov11_Gb.m` - an example script for the metamodel of nomal subjects to sample output k_T with different input G_B. 
`metamodel_normal_k_mean11_cov11_Gb.py`-python script to run this metamodel above using matlab python engine.
