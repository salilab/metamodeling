% Make a DBN for the GLP1R model with the following variables
%
% Time-dependent variables
%  -> GLP1.GLP1R(t)  ->  GLP1.GLP1R(t+1) ->
%  -> GLP1R.GLP1R(t)  ->  GLP1R.GLP1R(t+1) ->
%  -> Galpha.GLP1R(t)  ->  Galpha.GLP1R(t+1) ->
%  -> cAMP.GLP1R(t)  ->  cAMP.GLP1R(t+1) ->
%
% Reference variables
% cAMP.ref(t), cAMP.ref(t+1) 
%
% Observed variables
% cAMP.obs(t), cAMP.obs(t+1)
%
% Time-invariant variables
% 
%
% Parameters
%
% To generate a conditional gaussian model
function [dbn_factory]= make_glp1r_dbn_factory(Gex_GLP1R,GLP1_GLP1R, cAMP_GLP1R)
    node_names=  {'Gex.GLP1R','Gex.GLP1R.ref','Gex.GLP1R.obs','GLP1.GLP1R','GLP1.ref','GLP1.obs','cAMP.GLP1R','cAMP.ref'}; 
    n= length(node_names);
    % Intra - in one time slice
    edges_intra= {'Gex.GLP1R','Gex.GLP1R.ref';'Gex.GLP1R.ref','Gex.GLP1R.obs';'cAMP.GLP1R','cAMP.ref';'GLP1.GLP1R','GLP1.ref';'GLP1.ref','GLP1.obs'};
    % Inter - between time slices                             
    edges_inter= {'Gex.GLP1R','Gex.GLP1R';'GLP1.GLP1R','GLP1.GLP1R';'GLP1.GLP1R','cAMP.GLP1R';'cAMP.GLP1R','cAMP.GLP1R'}; 
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('Gex.GLP1R')= 'Gex.GLP1R.inter';  
    eclass2_map('GLP1.GLP1R')= 'GLP1.GLP1R.inter';  
    eclass2_map('cAMP.GLP1R')= 'cAMP.GLP1R.inter';       
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Gex.GLP1R', 0, ...
        {'mean', Gex_GLP1R, 'cov', 1E-8} ); 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Gex.GLP1R.ref', 0, ...
        {'mean', 0.0, 'cov', 1E-8,  'weights', 1.0} ); 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Gex.GLP1R.obs', 0, ...
        {'mean', 0.0, 'cov', 1E-8,  'weights', 1.0} ); 
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1.GLP1R', 0, ...
        {'mean', GLP1_GLP1R, 'cov', 1E-8} ); 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1.ref', 0, ...
        {'mean', 0.0, 'cov', 1E-8,  'weights', 1.0} ); 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1.obs', 0, ...
        {'mean', 0.0, 'cov', 1E-8,  'weights', 1.0} ); 
    
    % works only when there is ATP !!!!!!!!TODO!!!!!!!!
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cAMP.GLP1R', 0, ...
        {'mean', cAMP_GLP1R, 'cov', 1E-12});
    
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'cAMP.ref', 0, ...
        {'mean', 0.0, 'cov',1E-12,  'weights', 1.0} );  

    % eclass2 (time-slice t+1 with parents in the previous time slice)
    
    % CPD for G(t+1)
    weights_GLP11_T0= containers.Map(); % parents in slice t
    weights_GLP11_T1= containers.Map(); % parents in slice t+1
    weights_GLP11_T0('Gex.GLP1R')= 1.0;
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'Gex.GLP1R', 1, ...
        {'mean', 0.0, 'cov', 1E-8}, weights_GLP11_T0, weights_GLP11_T1);  
    
    weights_GLP11_T0= containers.Map(); % parents in slice t
    weights_GLP11_T1= containers.Map(); % parents in slice t+1
    weights_GLP11_T0('GLP1.GLP1R')= 1.0;
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'GLP1.GLP1R', 1, ...
        {'mean', 0.0, 'cov', 1E-8}, weights_GLP11_T0, weights_GLP11_T1);  

   
    weights_cAMP1_T0= containers.Map(); % parents in slice t
    weights_cAMP1_T1= containers.Map(); % parents in slice t+1
    weights_cAMP1_T0('GLP1.GLP1R')= 80;
    weights_cAMP1_T0('cAMP.GLP1R')= 0.5;
    weights_cAMP1_T0('Gex.GLP1R')= 1.0;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cAMP.GLP1R', 1, ...
        {'mean', 0.0, 'cov', 1E-12}, weights_cAMP1_T0, weights_cAMP1_T1);

    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

