% Make a DBN for the the virtual screening model with the following variables
%
% Time-dependent variables
%
% Coupling variables
% Gex.C, Gin.C, Scell.C
%
% Observed variables
% The node 'A.obs' corresponds to the 'A.GL' in the GL data model
%
% Time-invariant variables
%
% Parameters
%
% To generate a conditional gaussian model
%
% All variables display gaussian distributions.

% To generate a conditional gaussian model
function [dbn_factory]= make_screening_dbn(A_mean_screening, A_cov_screening, conc_mean_screening,...
                                        conc_conv_screening, GLP1R_mean_screening, GLP1R_cov_screening,...
                                        cov_scale_screening, k1_screening, k2_screening);
    
    node_names=  {'A.screening','A.C','A.obs','conc.screening','conc.C','conc.obs','GLP1R.screening','GLP1R.C'}; 
    n= length(node_names);
    
    % Intra - in one time slice
    edges_intra= {'A.screening','A.C';'A.C','A.obs';'conc.screening','conc.C';...
                  'conc.C','conc.obs';'GLP1R.screening','GLP1R.C'};
    
    % Inter - between time slices                             
    edges_inter= {'A.screening','A.screening';'conc.screening','conc.screening';'GLP1R.screening','GLP1R.screening';...
                  'A.screening','GLP1R.screening';'conc.screening','GLP1R.screening'}; 
    
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
    eclass2_map('A.screening')= 'A.screening.inter';  
    eclass2_map('conc.screening')= 'conc.screening.inter';  
    eclass2_map('GLP1R.screening')= 'GLP1R.screening.inter';        
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'A.screening', 0, ...
        {'mean', A_mean_screening, 'cov', A_cov_screening} ); % A - GLP1 agonist
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'A.C', 0, ...
        {'mean', 0.0, 'cov', A_cov_screening*cov_scale_screening,  'weights', 1.0} ); % A.C = 1.0 * A
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'A.obs', 0, ...
        {'mean', 0.0, 'cov', A_cov_screening*cov_scale_screening,  'weights', 1.0} ); % A.obs = 1.0 * A.C
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'conc.screening', 0, ...
        {'mean', conc_mean_screening, 'cov', conc_conv_screening} ); % conc - compound concentration
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'conc.C', 0, ...
        {'mean', 0.0, 'cov', conc_conv_screening*cov_scale_screening,  'weights', 1.0} ); % conc.C = 1.0 * conc
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'conc.obs', 0, ...
        {'mean', 0.0, 'cov', conc_conv_screening*cov_scale_screening,  'weights', 1.0} ); % conc.obs = 1.0 * conc.C
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.screening', 0, ...
        {'mean', GLP1R_mean_screening, 'cov', GLP1R_cov_screening} ); % GLP1R 
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.C', 0, ...
        {'mean', 0.0, 'cov', GLP1R_cov_screening*cov_scale_screening,  'weights', 1.0}); % GLP1R.C = 1.0 * GLP1R
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'A.screening', 1, ...
        {'mean', 0.0, 'cov', A_cov_screening*cov_scale_screening, 'weights', 1.0});  % A(t+1) = 1.0 * A(t)
    
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'conc.screening', 1, ...
        {'mean', 0.0, 'cov', conc_conv_screening*cov_scale_screening, 'weights', 1.0});  % conc(t+1) = 1.0 * conc(t)
    
    weights_screening1_T0= containers.Map(); % parents in slice t
    weights_screening1_T1= containers.Map(); % parents in slice t+1
    weights_screening1_T0('A.screening')= k1_screening;
    weights_screening1_T0('conc.screening')= k2_screening;
    weights_screening1_T0('GLP1R.screening')= 0.0;
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'GLP1R.screening', 1, ...
        {'mean', 0.0, 'cov', GLP1R_cov_screening*cov_scale_screening}, weights_screening1_T0, weights_screening1_T1);  % GLP1R(t+1) = 0.0 * GLP1R(t) + k1_screening * A(t) + k2_screening * conc(t)
   
    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end