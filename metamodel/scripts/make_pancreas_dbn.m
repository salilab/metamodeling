% Make a DBN for the pancreas model with the following variables
%
% Time-dependent variables
%
% Coupling variables
% Gex.C, Gin.C, Scell.C
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
                                          cov_scale_pancreas, Nc_pancreas, Ni_pancreas, ...
                                          S_pa_weight, S_pa_obs_weight, S_pa_obs_mean, S_pa_obs_cov);
    
    node_names=  {'Scell.C','Scell.pancreas','Sis.pancreas','Spa.pancreas','Spa.obs','Spa.C'}; 
    n= length(node_names);
    
    % Intra - in one time slice
    edges_intra= {'Scell.C','Scell.pancreas';'Scell.pancreas','Sis.pancreas';'Sis.pancreas','Spa.pancreas';...
                  'Spa.pancreas','Spa.C';'Spa.obs','Spa.C'};
              
%     node_names=  {'Scell.C','Scell.pancreas','Sis.pancreas','Spa.pancreas','Spa.C'}; 
%     n= length(node_names);
%     
%     % Intra - in one time slice
%     edges_intra= {'Scell.C','Scell.pancreas';'Scell.pancreas','Sis.pancreas';'Sis.pancreas','Spa.pancreas';...
%                   'Spa.pancreas','Spa.C'};
    
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
        CPDFactory('Gaussian_CPD', 'Scell.C', 0, ...
        {'mean', Scell_mean_pancreas, 'cov',  Scell_cov_pancreas} ); % Scell.C
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Scell.pancreas', 0, ...
        {'mean', 0.0, 'cov', Scell_cov_pancreas*cov_scale_pancreas, 'weights', 1.0} ); % Scell = Scell.C
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Sis.pancreas', 0, ...
        {'mean', 0.0, 'cov', Sislet_cov_pancreas*cov_scale_pancreas, 'weights', Nc_pancreas} ); % Sislet = Nc_pancreas * Scell
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Spa.pancreas', 0, ...
        {'mean', 0.0, 'cov', Spancreas_cov_pancreas*cov_scale_pancreas, 'weights', Ni_pancreas} ); % Spancreas = Ni_pancreas * Sislet
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Spa.C', 0, ...
        {'mean', 0.0, 'cov', Spancreas_cov_pancreas*cov_scale_pancreas, 'weights', [S_pa_weight, S_pa_obs_weight]} ); % Spancreas.C = Spancreas
 
%     CPDFactories{end+1}=  ...
%         CPDFactory('Gaussian_CPD', 'Spa.C', 0, ...
%         {'mean', 0.0, 'cov', Spancreas_cov_pancreas*cov_scale_pancreas, 'weights', 1.0} ); % Spancreas.C = Spancreas

        CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Spa.obs', 0, ...
        {'mean',  S_pa_obs_mean, 'cov', S_pa_obs_cov} ); % Spancreas_cov_pancreas*cov_scale_pancreas
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)

    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
end


