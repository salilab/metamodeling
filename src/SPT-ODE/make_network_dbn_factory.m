% Make a DBN for the network model with the following variables
%
% Time-dependent variables
%  -> G.Network(t)  ->  G.Network(t+1) ->
%  -> ATP.Network(t)  ->  ATP.Network(t+1) ->
%  -> cAMP.Network(t)  ->  cAMP.Network(t+1) ->
%  -> Ca.Network(t)  ->  Ca.Network(t+1) ->
%  -> S.Network(t)  ->  S.Network(t+1) ->
%
% Reference variables
% G.ref(t), G.ref(t+1)
% cAMP.ref(t), cAMP.ref(t+1)
% S.ref(t), S.ref(t+1)
%
% Observed variables
% G.obs(t), G.obs(t+1)
% cAMP.obs(t), cAMP.obs(t+1)
% S.obs(t), S.obs(t+1)
%
% Time-invariant variables
%
% Parameters
% PFK.Network
%
% To generate a conditional gaussian model
function [dbn_factory]= make_network_dbn_factory(G_Network, PFK_Network,S_Network)
    node_names=  {'G.Network','ATP.Network','cAMP.Network','cAMP.ref','cAMP.obs','Ca.Network',...
        'S.Network','S.ref','S.obs'}; 
    n= length(node_names);
    % Intra - in one time slice
    edges_intra= {'cAMP.Network','cAMP.ref';'cAMP.ref','cAMP.obs';'S.Network','S.ref';'S.ref','S.obs'};
    % Inter - between time slices
    edges_inter= {'G.Network', 'G.Network'; 'G.Network','ATP.Network';'ATP.Network','ATP.Network';...
        'ATP.Network','cAMP.Network';'cAMP.Network', 'cAMP.Network'; 'cAMP.Network', 'Ca.Network';...
        'Ca.Network','Ca.Network'; 'Ca.Network','S.Network';'S.Network','S.Network'}; 
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('G.Network')= 'G.Network.inter';
    eclass2_map('ATP.Network')= 'ATP.Network.inter';
    eclass2_map('cAMP.Network')= 'cAMP.Network.inter';   
    eclass2_map('Ca.Network')= 'Ca.Network.inter';    
    eclass2_map('S.Network')= 'S.Network.inter';   
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'G.Network', 0, ...
        {'mean', G_Network, 'cov', 0.0,'clamp_mean', 1, 'clamp_cov', 1, 'clamp_weights', 1} ); 
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ATP.Network', 0, ...
        {'mean', 3.3, 'cov', 1E-12});
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cAMP.Network', 0, ...
        {'mean', 1.3E-3, 'cov', 1E-12});
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cAMP.ref', 0, ...
        {'mean', 0.0, 'cov', 1E-12,   'weights', 1.0});
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cAMP.obs', 0, ...
        {'mean', 0.0, 'cov', 1E-12,   'weights', 1.0});

    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Ca.Network', 0, ...
        {'mean', 1E-4, 'cov', 1E-12,'clamp_mean', 1, 'clamp_cov', 1, 'clamp_weights', 1});

    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'S.Network', 0, ...
        {'mean', 30, 'cov', 1E-6});
    
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'S.ref', 0, ...
        { 'mean', 0.0,   'cov', 1E-6,   'weights', 1.0} ); % E= [1.0*I(t)] 
    
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'S.obs', 0, ...   
        {'mean', 0.0,   'cov', 1E-6,   'weights', 1.0} ); % E= [1.0*Iref(t)]
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    % CPD for G(t+1)
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'G.Network', 1, ...
        {'mean', 0.0, 'cov', 0.0,'clamp_mean', 1, 'clamp_cov', 1, 'clamp_weights', 1,'weights', 1.0});  

    weights_ATP1_T0= containers.Map(); % parents in slice t
    weights_ATP1_T1= containers.Map(); % parents in slice t+1
    INITIAL_PFK= PFK_Network;
    weights_ATP1_T0('G.Network')= INITIAL_PFK;
    weights_ATP1_T0('ATP.Network')= -1.0;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ATP.Network', 1, ...
        {'mean', 0.0, 'cov', 1E-12}, weights_ATP1_T0, weights_ATP1_T1);  % NETWORK_ATP ~ 0.4* NETWORK.Gpm + 0.4*NETWORK.PFK1 + 0.2*NETWORK.GCK; 
    
    weights_cAMP1_T0= containers.Map(); % parents in slice t
    weights_cAMP1_T1= containers.Map(); % parents in slice t+1
    weights_cAMP1_T0('ATP.Network')= -0.01/3;
    weights_cAMP1_T0('cAMP.Network')= 1;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cAMP.Network', 1, ...
        {'mean', 0.0105, 'cov', 1E-12}, weights_cAMP1_T0, weights_cAMP1_T1);  % NETWORK_ATP ~ 0.4* NETWORK.Gpm + 0.4*NETWORK.PFK1 + 0.2*NETWORK.GCK; 
    
    weights_Ca1_T0= containers.Map(); % parents in slice t
    weights_Ca1_T1= containers.Map(); % parents in slice t+1
    weights_Ca1_T0('cAMP.Network')= 0.6/0.25;
    weights_Ca1_T0('Ca.Network')= 1.0;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Ca.Network', 1, ...
        {'mean', -0.6*1.05*1E-3/0.25, 'cov', 1E-12}, weights_Ca1_T0, weights_Ca1_T1);

    weights_S1_T0= containers.Map(); % parents in slice t
    weights_S1_T1= containers.Map(); % parents in slice t+1
    weights_S1_T0('Ca.Network')= -1E+5/3;
    weights_S1_T0('S.Network')= 1.0;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'S.Network', 1, ...
        {'mean', 40/3, 'cov', 1E-6}, weights_S1_T0, weights_S1_T1);
   
    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

