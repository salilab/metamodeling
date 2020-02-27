% Make a DBN for the pancreas model with the following variables
%
% Time-dependent variables
%
% Reference variables
% Scell.ref
% Spancreas.ref
%
% Observed variables
%
% Time-invariant variables
% Scell
% Sislet
% Spancreas
%
% Parameters
%
% To generate a conditional gaussian model
function [dbn_factory]= make_pancreas_dbn_factory(Scell_Pancreas);
    node_names=  {'Scell.Pancreas','Scell.ref','Scell.obs','Sislet.Pancreas','Spancreas.Pancreas',...
        'Spancreas.ref', 'Spancreas.obs'}; 
    n= length(node_names);
    % Intra - in one time slice
    edges_intra= {'Scell.ref','Scell.obs';'Scell.ref','Scell.Pancreas';'Scell.Pancreas', 'Sislet.Pancreas'; ...
       'Sislet.Pancreas','Spancreas.Pancreas';'Spancreas.Pancreas', 'Spancreas.ref';'Spancreas.ref','Spancreas.obs' };
    % Inter - between time slices
    edges_inter= {}; 
    % edges_inter= {'Gex.Pancreas', 'Gex.Pancreas'; 'NES.Pancreas','NES.Pancreas';'GLP1.Pancreas','GLP1.Pancreas';...
        % 'NES.Pancreas','NES.Pancreas'}; 
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    % eclass2_map('Gex.Pancreas')= 'Gex.Pancreas.inter';
    % eclass2_map('NES.Pancreas')= 'NES.Pancreas.inter';
    % eclass2_map('GLP1.Pancreas')= 'GLP1.Pancreas.inter';   
    % eclass2_map('GLP1.Pancreas')= 'GLP1.Pancreas.inter';   
    % eclass2_map('NES.Pancreas')= 'NES.Pancreas.inter';   
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Scell.ref', 0, ...
        {'mean', Scell_Pancreas, 'cov', 1E-20, 'clamp_cov', 1E-20} ); % 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Scell.Pancreas', 0, ...
        {'mean', 0.0, 'cov', 1E-20, 'clamp_cov', 1E-20, 'weights', 1.0, 'clamp_weights', 1} ); % 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Scell.obs', 0, ...
        {'mean', 0.0, 'cov', 1E-20, 'clamp_cov', 1E-20, 'weights', 1.0, 'clamp_weights', 1} ); % 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Sislet.Pancreas', 0, ...
        {'mean', 0.0, 'cov', 1E-20, 'clamp_cov', 1E-20, 'weights', 300, 'clamp_weights', 300} ); % Approx. 300 beta-cells in one islet, Rorsman and Braun 2013
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Spancreas.Pancreas', 0, ...
        {'mean', 0.0, 'cov', 1E-15, 'clamp_cov', 1E-15, 'weights', 1E+6, 'clamp_weights', 1E+6} ); % Approx. 1 million islets in one pancreas, Hellman 1959; Rorsman and Braun 2013
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Spancreas.ref', 0, ...
        {'mean', 0.0, 'cov', 1E-15, 'clamp_cov', 1E-15, 'weights', 1.0, 'clamp_weights', 1} ); % 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Spancreas.obs', 0, ...
        {'mean', 0.0, 'cov', 1E-15, 'clamp_cov', 1E-15, 'weights', 1.0, 'clamp_weights', 1} ); % 
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    % CPD for G(t+1)
    
    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

