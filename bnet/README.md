Synopsis:
This repository is for developing meta-modeling implementations. 

## Factories and Functions:
* CPDFactory.m  
* DBNFactory.m   
* get_eclass_from_maps.m    
* get_reverse_nodes_map.m   
* get_valid_nodes_graph.m   
* merge_dbn_factories.m  

## DBN factories for (meta-)models:
* make_meal_dbn_factory_eq.m --- The meal model.
* make_pancreas_dbn_factory.m --- The pancreas model. 
* make_network_dbn_factory.m --- The network model  
* make_spt_dbn_factory.m --- The SPT model. 
* make_KEGG_dbn_factory.m --- The KEGG model.   
* make_glp1r_dbn_factory.m --- The GLP1R model. 
* make_meta_bnet5.m --- The metamodel of meal, network, SPT, KEGG and GP1R models

## Running individual model priors and posteriors:
* meal_model_prior_normal.m --- The meal model for normal subject. 
* meal_model_prior_t2d.m --- The meal model for T2D subject. 
* pancreas_model_prior_normal.m --- The pancreas model. 
* network_model_prior_normal.m --- The network model 
* spt_model_prior_normal.m --- The SPT model 
* GLP1R_model_prior_normal.m --- The GLP1R model  
* KEGG_model_prior_normal.m --- The KEGG model  

## Running the meta-model priors and posteriors:
### Normal subjects
***GLP1 effects on the meal model***
* metamodeling5_normal_meal.m --- To run the metamodel of a nomal subject after a meal. 
* metamodeling5_normal_meal_glp1_medium.m ---  To run the metamodel after of a nomal subject a meal, with medium GLP1. 
* metamodeling5_normal_meal_glp1_high.m --- To run the metamodel of a nomal subject after a meal, with high GLP1. 

***GLP1 effects on the model variables under basal condition***
* metamodeling5_normal_variables.m --- To run the metamodel of a nomal subject at basal condition (fasting). 
* metamodeling5_normal_variables_glp1_medium.m --- To run the metamodel of a nomal subject at basal condition (fasting), with medium GLP1. 
* metamodeling5_normal_variables_glp1_high.m --- To run the metamodel ofa nomal subject at basal condition (fasting), with high GLP1. 

***Effects of different basal glucose concentrations on model variables under basal condition***
* metamodeling5_normal_variables_Gb.m --- To run the metamodel of a nomal subject at basal condition (fasting), with different basal glucoes Gb. 
* metamodeling5_normal_variables_highGb.m --- To run the metamodel of a nomal subject at basal condition (fasting), with high basal glucoes Gb. 

### T2D subjects
***GLP1 effects on the meal model***
* metamodeling5_t2d_meal.m --- To run the metamodel of a t2d subject after a meal. 
* metamodeling5_t2d_meal_glp1_medium.m ---  To run the metamodel after of a t2d subject a meal, with medium GLP1. 
* metamodeling5_t2d_meal_glp1_high.m --- To run the metamodel of a t2d subject after a meal, with high GLP1. 

***GLP1 effects on the model variables under basal condition***
* metamodeling5_t2d_variables.m --- To run the metamodel of a t2d subject at basal condition (fasting). 
* metamodeling5_t2d_variables_glp1_medium.m --- To run the metamodel of a t2d subject at basal condition (fasting), with medium GLP1. 
* metamodeling5_t2d_variables_glp1_high.m --- To run the metamodel of a t2d subject at basal condition (fasting), with high GLP1. 
