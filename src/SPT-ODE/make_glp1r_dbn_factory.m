% Make a DBN for the GLP1R model with the following variables
%
% Time-dependent variables
%  -> GLP1R.GLP1(t)  ->  GLP1R.GLP1(t+1) ->
%  -> GLP1R.GLP1R(t)  ->  GLP1R.GLP1R(t+1) ->
%  -> GLP1R.Galpha(t)  ->  GLP1R.Galpha(t+1) ->
%  -> GLP1R.cAMP(t)  ->  GLP1R.cAMP(t+1) ->
%
% Reference variables
% 
%
% Observed variables
% 
%
% Time-invariant variables
% 
%
% Parameters
%
% To generate a conditional gaussian model
function [dbn_factory]= make_glp1r_dbn_factory(Im, time)
    node_names=  {'GLP1R.GLP1','GLP1R.GLP1R', 'GLP1R.Galpha','GLP1R.cAMP'}; 
    n= length(node_names);
    % Intra - in one time slice
    edges_intra= {'GLP1R.GLP1','GLP1R.GLP1R';'GLP1R.GLP1R','GLP1R.Galpha'};
    % Inter - between time slices
    edges_inter= {'GLP1R.GLP1', 'GLP1R.GLP1'; 'GLP1R.GLP1R','GLP1R.GLP1R';'GLP1R.Galpha','GLP1R.Galpha';'GLP1R.Galpha','GLP1R.cAMP';'GLP1R.cAMP','GLP1R.cAMP'}; % BARAK comment: switched G->I and h->I to (G-h)->I
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('GLP1R.GLP1')= 'GLP1R.GLP1.inter';
    eclass2_map('GLP1R.GLP1R')= 'GLP1R.GLP1R.inter';
    eclass2_map('GLP1R.Galpha')= 'GLP1R.Galpha.inter';   
    eclass2_map('GLP1R.cAMP')= 'GLP1R.cAMP.inter';       
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.GLP1', 0, ...
        {'mean', 6.1, 'cov', 0.1} ); 
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.GLP1R', 0, ...
        {'mean', 10, 'cov', 0.1, 'weights',1});

    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.Galpha', 0, ...
        {'mean', 2, 'cov', 0.01,'weights',1});
    % works only when there is ATP !!!!!!!!TODO!!!!!!!!
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.cAMP', 0, ...
        {'mean', 0.1, 'cov', 0.1});
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    
    % CPD for G(t+1)
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'GLP1R.GLP1', 1, ...
        {'mean', 0.0, 'cov', 0.1,  'weights', 1.0} );  

    weights_GLP1R_GLP1R_T0= containers.Map(); % parents in slice t
    weights_GLP1R_GLP1R_T1= containers.Map(); % parents in slice t+1
    weights_GLP1R_GLP1R_T1('GLP1R.GLP1')= 0.4;
    weights_GLP1R_GLP1R_T0('GLP1R.GLP1R')= 0.4;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.GLP1R', 1, ...
        {'mean', 0.1, 'cov', 0.01}, weights_GLP1R_GLP1R_T0, weights_GLP1R_GLP1R_T1);  % NETWORK_ATP ~ 0.4* NETWORK.Gpm + 0.4*NETWORK.PFK1 + 0.2*NETWORK.GCK; 
    
    weights_GLP1R_Galpha_T0= containers.Map(); % parents in slice t
    weights_GLP1R_Galpha_T1= containers.Map(); % parents in slice t+1
    weights_GLP1R_Galpha_T1('GLP1R.GLP1R')= 0.4;
    weights_GLP1R_Galpha_T0('GLP1R.Galpha')= 0.4;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.Galpha', 1, ...
        {'mean', 0.1, 'cov', 0.01}, weights_GLP1R_Galpha_T0, weights_GLP1R_Galpha_T1);

    weights_GLP1R_cAMP_T0= containers.Map(); % parents in slice t
    weights_GLP1R_cAMP_T1= containers.Map(); % parents in slice t+1
    weights_GLP1R_cAMP_T0('GLP1R.Galpha')= 0.4;
    weights_GLP1R_cAMP_T0('GLP1R.cAMP')= 0.4;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.cAMP', 1, ...
        {'mean', 0.1, 'cov', 0.01}, weights_GLP1R_cAMP_T0, weights_GLP1R_cAMP_T1);

    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

