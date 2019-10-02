% Make a DBN for the network model with the following variables
%
% Time-dependent variables
%  -> Gin.Network(t)  ->  Gin.Network(t+1) ->
%  -> ATP.Network(t)  ->  ATP.Network(t+1) ->
%  -> cAMP.Network(t)  ->  cAMP.Network(t+1) ->
%  -> Ca.Network(t)  ->  Ca.Network(t+1) ->
%  -> S.Network(t)  ->  S.Network(t+1) ->
%
% Reference variables
% Gin.ref(t), Gin.ref(t+1)
% cAMP.ref(t), cAMP.ref(t+1)
% S.ref(t), S.ref(t+1)
%
% Observed variables
% Gin.obs(t), Gin.obs(t+1)
% cAMP.obs(t), cAMP.obs(t+1)
% S.obs(t), S.obs(t+1)
%
% Time-invariant variables
%
% Parameters
% PFK.Network
%
% To generate a conditional gaussian model
function [dbn_factory]= make_network_dbn_factory(Gin_Network, PFK_Network, ATP_Network, GLP1_Network, GLP1R_Network, cAMP_Network, Ca_Network, S_Network, I_Network)
    node_names=  {'Gin.Network','Gin.ref','ATP.Network','GLP1.Network','GLP1R.Network','GLP1R.ref','cAMP.Network','Ca.Network',...
        'S.Network','I.Network','Scell.ref'}; 
    n= length(node_names);
    % Intra - in one time slice
    edges_intra= {'Gin.ref','Gin.Network';'GLP1R.ref','GLP1R.Network';...
       'S.Network','Scell.ref'};
    % Inter - between time slices
    edges_inter= {'Gin.Network', 'Gin.Network'; 'Gin.Network','ATP.Network';'ATP.Network','ATP.Network';...
        'ATP.Network','cAMP.Network';'GLP1.Network','GLP1.Network';'GLP1.Network','GLP1R.Network'; 'GLP1R.Network', 'GLP1R.Network';...
        'GLP1R.Network','cAMP.Network';'cAMP.Network', 'cAMP.Network'; 'cAMP.Network', 'Ca.Network';...
        'Ca.Network','Ca.Network'; 'Ca.Network','S.Network';'S.Network','S.Network';'S.Network','I.Network';'I.Network','I.Network'}; 
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('Gin.Network')= 'Gin.Network.inter';
    eclass2_map('ATP.Network')= 'ATP.Network.inter';
    eclass2_map('GLP1.Network')= 'GLP1.Network.inter';   
    eclass2_map('GLP1R.Network')= 'GLP1R.Network.inter';   
    eclass2_map('cAMP.Network')= 'cAMP.Network.inter';   
    eclass2_map('Ca.Network')= 'Ca.Network.inter';    
    eclass2_map('S.Network')= 'S.Network.inter';   
    eclass2_map('I.Network')= 'I.Network.inter'; 
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Gin.ref', 0, ...
        {'mean', Gin_Network, 'cov', 1E-12} ); % Gin.ref = 1.0 * Gin
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Gin.Network', 0, ...
        {'mean', 0.0, 'cov', 1E-12,   'weights', 1.0} ); % Gin
    %CPDFactories{end+1}=  ...
    %    CPDFactory('Gaussian_CPD', 'Gin.obs', 0, ...
    %    {'mean', 0.0, 'cov', 1E-12,   'weights', 1.0} ); % Gin.obs = 1.0 * Gin.ref
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ATP.Network', 0, ...
        {'mean', ATP_Network, 'cov', 1E-12}); % ATP
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1.Network', 0, ...
        {'mean', GLP1_Network, 'cov', 1E-18}); % GLP1
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.Network', 0, ...
        {'mean', GLP1R_Network, 'cov', 1E-12,   'weights', 0.0}); % GLP1R = 1.0 * GLP1R.ref
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.ref', 0, ...
        {'mean', GLP1R_Network, 'cov', 1E-12}); % GLP1R.ref
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cAMP.Network', 0, ...
        {'mean', cAMP_Network, 'cov', 1E-12}); % cAMP
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Ca.Network', 0, ...
        {'mean', Ca_Network, 'cov', 1E-12}); % Ca2+
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'S.Network', 0, ...
        {'mean', S_Network, 'cov', 1E-4}); % S
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'I.Network', 0, ...
        {'mean', I_Network, 'cov', 1E-4}); % I
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Scell.ref', 0, ...
        { 'mean', 0.0,   'cov', 1E-4,   'weights', 1.0} ); % Scell.ref = 1.0 * I

    % eclass2 (time-slice t+1 with parents in the previous time slice)
    % CPD for G(t+1)
    weights_Gin1_T0= containers.Map(); % parents in slice t
    weights_Gin1_T1= containers.Map(); % parents in slice t+1
    weights_Gin1_T0('Gin.Network')= 0.0;
    weights_Gin1_T1('Gin.ref')= 1.0;
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'Gin.Network', 1, ...
        {'mean', 0.0, 'cov', 1E-12}, weights_Gin1_T0, weights_Gin1_T1);  % Gin(t+1) = 1.0 * Gin(t+1)

    weights_ATP1_T0= containers.Map(); % parents in slice t
    weights_ATP1_T1= containers.Map(); % parents in slice t+1
    INITIAL_PFK= PFK_Network;
    weights_ATP1_T0('Gin.Network')= INITIAL_PFK;
    weights_ATP1_T0('ATP.Network')= 1.65/3.3;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ATP.Network', 1, ...
        {'mean', 0.0, 'cov', 1E-12}, weights_ATP1_T0, weights_ATP1_T1);  % ATP(t+1) = 1.65/3.3 * ATP(t) + PFK_Network * Gin (t)
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1.Network', 1, ...
        {'mean', 0.0, 'cov', 1E-18, 'weights', 1.0}); % GLP1(t+1) = 1.0 * GLP1(t)
    
    weights_GLP1R1_T0= containers.Map(); % parents in slice t
    weights_GLP1R1_T1= containers.Map(); % parents in slice t+1
    weights_GLP1R1_T0('GLP1.Network')= 1/GLP1_Network;
    weights_GLP1R1_T0('GLP1R.Network')= 0.0;
    weights_GLP1R1_T1('GLP1R.ref')= 0.0;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.Network', 1, ...
        {'mean', 0.0, 'cov', 1E-12}, weights_GLP1R1_T0, weights_GLP1R1_T1); % GLP1R(t+1) = 0.0 * GLP1R(t) + 0.5 * GLP1R.ref(t) + 0.5/GLP1_Network * GLP1(t)
    
    weights_cAMP1_T0= containers.Map(); % parents in slice t
    weights_cAMP1_T1= containers.Map(); % parents in slice t+1
    weights_cAMP1_T0('ATP.Network')= 0.00013;
    weights_cAMP1_T0('GLP1R.Network')= 1E-3/3;
    weights_cAMP1_T0('cAMP.Network')= 1/3;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cAMP.Network', 1, ...
        {'mean', 0.0, 'cov', 1E-12}, weights_cAMP1_T0, weights_cAMP1_T1);  % cAMP(t+1) = 1/3 * cAMP(t) + 1E-3/3*GLP1R(t) + 0.00013* ATP(t)
    
    weights_Ca1_T0= containers.Map(); % parents in slice t
    weights_Ca1_T1= containers.Map(); % parents in slice t+1
    weights_Ca1_T0('cAMP.Network')= 0.05/1.3;
    weights_Ca1_T0('Ca.Network')= 0.5;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Ca.Network', 1, ...
        {'mean', 0.0, 'cov', 1E-12}, weights_Ca1_T0, weights_Ca1_T1); % Ca(t+1) = 0.5 * Ca(t) + 0.05/1.3 * cAMP(t)

    weights_S1_T0= containers.Map(); % parents in slice t
    weights_S1_T1= containers.Map(); % parents in slice t+1
    weights_S1_T0('Ca.Network')= 34/1E-4;
    weights_S1_T0('S.Network')= 0.0;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'S.Network', 1, ...
        {'mean', 0.0, 'cov', 1E-4}, weights_S1_T0, weights_S1_T1); % S(t+1) = 0.0 * S(t) + 34/1E-4 * Ca(t) 

    weights_I1_T0= containers.Map(); % parents in slice t
    weights_I1_T1= containers.Map(); % parents in slice t+1
    weights_I1_T0('S.Network')= 25/2/34;
    weights_I1_T0('I.Network')= 0.5;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'I.Network', 1, ...
        {'mean', 0.0, 'cov', 1E-4}, weights_I1_T0, weights_I1_T1); % I(t+1) = 0.5 * I(t) + 25/2/34 * S(t)
    
    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

