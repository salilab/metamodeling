% Make a DBN for the exocytosis model with the following variables
%
% Time-dependent variables
%
% Reference variables
% Gex.ref, Gin.ref, Scell.ref
%
% Observed variables
%
% Time-invariant variables
%
% Parameters
%
% To generate a conditional gaussian model
%
% All variables display gaussian distributions.

% To generate a conditional gaussian model
function [dbn_factory]= make_GLP1R_dbn(GLP1a_mean_GLP1R, GLP1a_cov_GLP1R, conc_mean_GLP1R,...
                                        conc_conv_GLP1R, GLP1R_mean_GLP1R, GLP1R_cov_GLP1R,...
                                        min_cov_GLP1R, GLP1a_w_GLP1R_GLP1R, conc_w_GLP1R_GLP1R);
    
    node_names=  {'GLP1a.GLP1R','GLP1a.ref','GLP1a.obs','conc.GLP1R','conc.ref','conc.obs','GLP1R.GLP1R','GLP1R.ref'}; 
    n= length(node_names);
    
    % Intra - in one time slice
    edges_intra= {'GLP1a.GLP1R','GLP1a.ref';'GLP1a.ref','GLP1a.obs';'conc.GLP1R','conc.ref';...
                  'conc.ref','conc.obs';'GLP1R.GLP1R','GLP1R.ref'};
    
    % Inter - between time slices                             
    edges_inter= {'GLP1a.GLP1R','GLP1a.GLP1R';'conc.GLP1R','conc.GLP1R';'GLP1R.GLP1R','GLP1R.GLP1R';...
                  'GLP1a.GLP1R','GLP1R.GLP1R';'conc.GLP1R','GLP1R.GLP1R'}; 
    
    % 'Equivalence classes' specify how the template is initiated and rolled
    % Specify which CPD is associates with each node in either time
    % slice 1 (eclass1) or in slice 2 onwards (eclass2)
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('GLP1a.GLP1R')= 'GLP1a.GLP1R.inter';  
    eclass2_map('conc.GLP1R')= 'conc.GLP1R.inter';  
    eclass2_map('GLP1R.GLP1R')= 'GLP1R.GLP1R.inter';        
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1a.GLP1R', 0, ...
        {'mean', GLP1a_mean_GLP1R, 'cov', GLP1a_cov_GLP1R} ); % GLP1a - GLP1R analogues
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1a.ref', 0, ...
        {'mean', 0.0, 'cov', GLP1a_cov_GLP1R*min_cov_GLP1R,  'weights', 1.0} ); % GLP1a.ref = 1.0 * GLP1a
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1a.obs', 0, ...
        {'mean', 0.0, 'cov', GLP1a_cov_GLP1R*min_cov_GLP1R,  'weights', 1.0} ); % GLP1a.obs = 1.0 * GLP1a.ref
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'conc.GLP1R', 0, ...
        {'mean', conc_mean_GLP1R, 'cov', conc_conv_GLP1R} ); % conc - compound concentration
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'conc.ref', 0, ...
        {'mean', 0.0, 'cov', conc_conv_GLP1R*min_cov_GLP1R,  'weights', 1.0} ); % conc.ref = 1.0 * conc
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'conc.obs', 0, ...
        {'mean', 0.0, 'cov', conc_conv_GLP1R*min_cov_GLP1R,  'weights', 1.0} ); % conc.obs = 1.0 * conc.ref
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.GLP1R', 0, ...
        {'mean', GLP1R_mean_GLP1R, 'cov', GLP1R_cov_GLP1R} ); % GLP1R 
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.ref', 0, ...
        {'mean', 0.0, 'cov', GLP1R_cov_GLP1R*min_cov_GLP1R,  'weights', 1.0}); % GLP1R.ref = 1.0 * GLP1R
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'GLP1a.GLP1R', 1, ...
        {'mean', 0.0, 'cov', GLP1a_cov_GLP1R*min_cov_GLP1R, 'weights', 1.0});  % GLP1a(t+1) = 1.0 * GLP1a(t)
    
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'conc.GLP1R', 1, ...
        {'mean', 0.0, 'cov', conc_conv_GLP1R*min_cov_GLP1R, 'weights', 1.0});  % conc(t+1) = 1.0 * conc(t)
    
    weights_GLP1R1_T0= containers.Map(); % parents in slice t
    weights_GLP1R1_T1= containers.Map(); % parents in slice t+1
    weights_GLP1R1_T0('GLP1a.GLP1R')= GLP1a_w_GLP1R_GLP1R;
    weights_GLP1R1_T0('conc.GLP1R')= conc_w_GLP1R_GLP1R;
    weights_GLP1R1_T0('GLP1R.GLP1R')= 0.0;
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'GLP1R.GLP1R', 1, ...
        {'mean', 0.0, 'cov', GLP1R_cov_GLP1R*min_cov_GLP1R}, weights_GLP1R1_T0, weights_GLP1R1_T1);  % GLP1R(t+1) = 0.0 * GLP1R(t) + GLP1a_w_GLP1R_GLP1R * GLP1a(t) + conc_w_GLP1R_GLP1R * conc(t)
   
    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end