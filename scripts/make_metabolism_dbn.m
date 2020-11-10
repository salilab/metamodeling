% Make a DBN for the exocytosis model with the following variables
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
function [dbn_factory]= make_metabolism_dbn(Gex_mean_metabolism, Gex_cov_metabolism, Ex4_mean_metabolism, ...
                                            Ex4_cov_metabolism, Pathway_mean_metabolism, Pathway_cov_metabolism, ATP_mean_metabolism, ...
                                            ATP_cov_metabolism, cov_scale_metabolism,...
                                            k1_metabolism, k2_metabolism, ...
                                            k3_metabolism);
    
    node_names=  {'Gex.metabolism','Gex.C','Gex.obs','Ex4.metabolism','Ex4.C',...
                  'Ex4.obs','Pathway.metabolism','ATP.metabolism','ATP.C','ATP.obs'}; 
    n= length(node_names);
    
    % Intra - in one time slice
    edges_intra= {'Gex.metabolism','Gex.C';'Gex.C','Gex.obs';'Gex.metabolism','Pathway.metabolism';...
                  'Ex4.metabolism','Pathway.metabolism';'Ex4.metabolism','Ex4.C';'Ex4.C','Ex4.obs';...
                  'Pathway.metabolism','ATP.metabolism';'ATP.metabolism','ATP.C';...
                  'ATP.C','ATP.obs'};
    
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
        CPDFactory('Gaussian_CPD', 'Gex.metabolism', 0, ...
        {'mean', Gex_mean_metabolism, 'cov', Gex_cov_metabolism} ); % Gex.C
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Gex.C', 0, ...
        {'mean', 0.0, 'cov', Gex_cov_metabolism*cov_scale_metabolism, 'weights', 1.0} ); % Gex = Gex.C
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Gex.obs', 0, ...
        {'mean', 0.0, 'cov', Gex_cov_metabolism*cov_scale_metabolism, 'weights', 1.0} ); % Gex.obs = Gex.C
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Ex4.metabolism', 0, ...
        {'mean', Ex4_mean_metabolism, 'cov', Ex4_cov_metabolism} ); % Ex4
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Ex4.C', 0, ...
        {'mean', 0.0, 'cov', Ex4_cov_metabolism*cov_scale_metabolism, 'weights', 1.0} ); % Ex4.C = Ex4
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Ex4.obs', 0, ...
        {'mean', 0.0, 'cov', Ex4_cov_metabolism*cov_scale_metabolism, 'weights', 1.0} ); % Ex4.obs = Ex4.C
    
    weights_Pathway_map_T0= containers.Map();
    weights_Pathway_map_T1= containers.Map();
    weights_Pathway_map_T0('Gex.metabolism')= k1_metabolism;
    weights_Pathway_map_T0('Ex4.metabolism')= k2_metabolism;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Pathway.metabolism', 0, ...
        {'mean', Pathway_mean_metabolism, 'cov', Pathway_cov_metabolism*cov_scale_metabolism}, ...
        weights_Pathway_map_T0, weights_Pathway_map_T1); % Pathway = k1_metabolism * Gex + k2_metabolism * Ex4
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ATP.metabolism', 0, ...
        {'mean', ATP_mean_metabolism, 'cov', ATP_cov_metabolism*cov_scale_metabolism, 'weights', k3_metabolism} ); % ATP =  k3_metabolism * Pathwhay
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ATP.C', 0, ...
        {'mean', 0.0, 'cov', ATP_cov_metabolism*cov_scale_metabolism, 'weights', 1.0} ); % ATP.C = ATP
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ATP.obs', 0, ...
        {'mean', 0.0, 'cov', ATP_cov_metabolism*cov_scale_metabolism, 'weights', 1.0} ); % ATP.obs = ATP.C
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)

    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

