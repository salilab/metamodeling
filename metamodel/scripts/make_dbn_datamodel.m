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
function [dbn_factory]= make_dbn_datamodel(Gcell_mean_data, Gcell_cov_data, Gcell_mean_obs, Gcell_cov_obs, G_cov_data, cov_scale_data);
    
    node_names=  {'GcellD.obs','Gcell.data','Gcell.C'}; 
    n= length(node_names);
    
    % Intra - in one time slice
    edges_intra= {'GcellD.obs','Gcell.data';'Gcell.data','Gcell.C'};
    
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
        CPDFactory('Gaussian_CPD', 'Gcell.data', 0, ...
            {'mean',  0.0, 'cov', G_cov_data*cov_scale_data, 'weights', 1.0} ); %
    
    CPDFactories{end+1}=  ...
    CPDFactory('Gaussian_CPD', 'GcellD.obs', 0, ...
    {'mean',  Gcell_mean_obs, 'cov', Gcell_cov_obs} ); % 

    CPDFactories{end+1}=  ...
    CPDFactory('Gaussian_CPD', 'Gcell.C', 0, ...
    {'mean',  0.0, 'cov', G_cov_data*cov_scale_data, 'weights', 1.0} ); % 
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)

    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
end


