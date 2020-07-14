% Make a DBN for the pancreas model with the following variables
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
function [dbn_factory]= make_pancreas_dbn(Scell_mean_pancreas, Scell_cov_pancreas, Sislet_mean_pancreas, ...
                                          Sislet_cov_pancreas, Spancreas_mean_pancreas, Spancreas_cov_pancreas, ...
                                          min_cov_pancreas, Scell_w_Sislet_pancreas, Sislet_w_Spancreas_pancreas);
    
    node_names=  {'Scell.ref','Scell.pancreas','Sislet.pancreas','Spancreas.pancreas','Spancreas.ref'}; 
    n= length(node_names);
    
    % Intra - in one time slice
    edges_intra= {'Scell.ref','Scell.pancreas';'Scell.pancreas','Sislet.pancreas';'Sislet.pancreas','Spancreas.pancreas';...
                  'Spancreas.pancreas','Spancreas.ref'};
    
    % Inter - between time slices
    edges_inter= {}; 
    
    % 'Equivalence classes' specify how the template is initiated and rolled
    % Specify which CPD is associates with each node in either time
    % slice 1 (eclass1) or in slice 2 onwards (eclass2)
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; 
    end
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Scell.ref', 0, ...
        {'mean', Scell_mean_pancreas, 'cov',  Scell_cov_pancreas} ); % Scell.ref
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Scell.pancreas', 0, ...
        {'mean', 0.0, 'cov', Scell_cov_pancreas*min_cov_pancreas, 'weights', 1.0} ); % Scell = Scell.ref
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Sislet.pancreas', 0, ...
        {'mean', 0.0, 'cov', Sislet_cov_pancreas*min_cov_pancreas, 'weights', Scell_w_Sislet_pancreas} ); % Sislet = Scell_w_Sislet_pancreas * Scell
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Spancreas.pancreas', 0, ...
        {'mean', 0.0, 'cov', Spancreas_cov_pancreas*min_cov_pancreas, 'weights', Sislet_w_Spancreas_pancreas} ); % Spancreas = Sislet_w_Spancreas_pancreas * Sislet
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Spancreas.ref', 0, ...
        {'mean', 0.0, 'cov', Spancreas_cov_pancreas*min_cov_pancreas, 'weights', 1.0} ); % Spancreas.ref = Spancreas
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)

    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

