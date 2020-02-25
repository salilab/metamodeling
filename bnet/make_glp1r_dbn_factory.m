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
function [dbn_factory]= make_glp1r_dbn_factory(GLP1a_GLP1R, cons_GLP1R,GLP1R_GLP1R)
    node_names=  {'GLP1a.GLP1R','GLP1a.ref','GLP1a.obs','cons.GLP1R','cons.ref','cons.obs','GLP1R.GLP1R','GLP1R.ref'}; 
    n= length(node_names);
    % Intra - in one time slice
    edges_intra= {'GLP1a.GLP1R','GLP1a.ref';'GLP1a.ref','GLP1a.obs';'cons.GLP1R','cons.ref';'cons.ref','cons.obs';...
        'GLP1R.GLP1R','GLP1R.ref'};
    % Inter - between time slices                             
    edges_inter= {'GLP1a.GLP1R','GLP1a.GLP1R';'cons.GLP1R','cons.GLP1R';'GLP1R.GLP1R','GLP1R.GLP1R';...
        'GLP1a.GLP1R','GLP1R.GLP1R';'cons.GLP1R','GLP1R.GLP1R'}; 
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('GLP1a.GLP1R')= 'GLP1a.GLP1R.inter';  
    eclass2_map('cons.GLP1R')= 'cons.GLP1R.inter';  
    eclass2_map('GLP1R.GLP1R')= 'GLP1R.GLP1R.inter';        
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1a.GLP1R', 0, ...
        {'mean', GLP1a_GLP1R, 'cov', 1E-18} ); % GLP1a - GLP1R analogues
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1a.ref', 0, ...
        {'mean', 0.0, 'cov', 1E-18,  'weights', 1.0} ); % GLP1a.ref = 1.0 * GLP1a
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1a.obs', 0, ...
        {'mean', 0.0, 'cov', 1E-18,  'weights', 1.0} ); % GLP1a.obs = 1.0 * GLP1a.reef
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cons.GLP1R', 0, ...
        {'mean', cons_GLP1R, 'cov', 1E-18} ); % cons - compound concentration
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cons.ref', 0, ...
        {'mean', 0.0, 'cov', 1E-18,  'weights', 1.0} ); % cons.ref = 1.0 * cons
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cons.obs', 0, ...
        {'mean', 0.0, 'cov', 1E-18,  'weights', 1.0} ); % cons.obs = 1.0 * cons.ref
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.GLP1R', 0, ...
        {'mean', GLP1R_GLP1R, 'cov', 1E-18} ); % GLP1R 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.ref', 0, ...
        {'mean', 0.0, 'cov', 1E-18,  'weights', 1.0}); % GLP1R.ref = 1.0 * GLP1R
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)
 
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'GLP1a.GLP1R', 1, ...
        {'mean', 0.0, 'cov', 1E-12, 'weights', 1.0});  % GLP1a(t+1) = 1.0 * GLP1a(t)
    
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'cons.GLP1R', 1, ...
        {'mean', 0.0, 'cov', 1E-12, 'weights', 1.0});  % cons(t+1) = 1.0 * cons(t)
    
    weights_GLP1R1_T0= containers.Map(); % parents in slice t
    weights_GLP1R1_T1= containers.Map(); % parents in slice t+1
    weights_GLP1R1_T0('GLP1a.GLP1R')= 1.2;
    weights_GLP1R1_T0('cons.GLP1R')= 1.2;
    weights_GLP1R1_T0('GLP1R.GLP1R')= 0.0;
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'GLP1R.GLP1R', 1, ...
        {'mean', 0.0, 'cov', 1E-18}, weights_GLP1R1_T0, weights_GLP1R1_T1);  % GLP1R(t+1) = 0.0 * GLP1R(t) + 1.0 * GLP1a(t) + 1.0 * cons(t)
   
    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end