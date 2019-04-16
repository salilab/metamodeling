% Make a DBN for the network model with the following variables
%
% Time-dependent variables
%  -> NETWORK.Gpm(t)  ->  NETWORK.Gpm(t+1) ->
%  -> NETWORK.ATP(t)  ->  NETWORK.ATP(t+1) ->
%  -> NETWORK.cAMP(t)  ->  NETWORK.cAMP(t+1) ->
%  -> NETWORK.Calcium(t)  ->  NETWORK.Calcium(t+1) ->
%  -> NETWORK.Ipm(t)  ->  NETWORK.Ipm(t+1) ->
%
% Reference variables
% R.Ipm(t), R.Ipm(t+1)
%
% Observed variables
% E.Ipm(t), E.Ipm(t+1)
%
% Time-invariant variables
% NETWORK.PFK1, NETWORK.GCK
%
% Parameters
%
% To generate a conditional gaussian model
function [dbn_factory]= make_network_dbn_factory(Im, time)
    node_names=  {'NETWORK.Gpm','NETWORK.PFK1', 'NETWORK.GCK','NETWORK.ATP','NETWORK.cAMP','NETWORK.Calcium','NETWORK.Ipm','R.Ipm','E.Ipm'}; 
    n= length(node_names);
    % Intra - in one time slice
    edges_intra= {'NETWORK.PFK1','NETWORK.ATP';'NETWORK.GCK','NETWORK.ATP';'NETWORK.Ipm','R.Ipm';'R.Ipm','E.Ipm'};
    % Inter - between time slices
    edges_inter= {'NETWORK.Gpm', 'NETWORK.Gpm'; 'NETWORK.Gpm','NETWORK.ATP';'NETWORK.ATP','NETWORK.ATP';'NETWORK.ATP','NETWORK.cAMP'; 'NETWORK.cAMP', 'NETWORK.cAMP'; 'NETWORK.cAMP', 'NETWORK.Calcium';'NETWORK.Calcium','NETWORK.Calcium'; 'NETWORK.Calcium','NETWORK.Ipm';'NETWORK.Ipm','NETWORK.Ipm'}; % BARAK comment: switched G->I and h->I to (G-h)->I
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('NETWORK.Gpm')= 'NETWORK.Gpm.inter';
    eclass2_map('NETWORK.ATP')= 'NETWORK.ATP.inter';
    eclass2_map('NETWORK.cAMP')= 'NETWORK.cAMP.inter';   
    eclass2_map('NETWORK.Calcium')= 'NETWORK.Calcium.inter';    
    eclass2_map('NETWORK.Ipm')= 'NETWORK.Ipm.inter';   
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'NETWORK.Gpm', 0, ...
        {'mean', 6.1, 'cov', 0.1} ); 
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'NETWORK.PFK1', 0, ...
        {'mean', 10, 'cov', 0.001});

    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'NETWORK.GCK', 0, ...
        {'mean', 2, 'cov', 0.001});

    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'NETWORK.ATP', 0, ...
        {'mean', 0.1, 'cov', 0.001, 'weights', [1.0 1.0]});
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'NETWORK.cAMP', 0, ...
        {'mean', 0.1, 'cov', 0.001});

    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'NETWORK.Calcium', 0, ...
        {'mean', 0.1, 'cov', 0.001});

    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'NETWORK.Ipm', 0, ...
        {'mean', 0.1, 'cov', 0.001});
    
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'R.Ipm', 0, ...
        { 'mean', 0.0,   'cov', 0.5,   'weights', 1.0} ); % E= [1.0*I(t)] 
    
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'E.Ipm', 0, ...   
        {'mean', 0.0,   'cov', 0.5,   'weights', 1.0} ); % E= [1.0*Iref(t)]
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    
    % CPD for G(t+1)
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'NETWORK.Gpm', 1, ...
        {'mean', 0.0, 'cov', 1e-2,  'weights', 1.0} );  

    weights_NETWORK_ATP_T0= containers.Map(); % parents in slice t
    weights_NETWORK_ATP_T1= containers.Map(); % parents in slice t+1
    weights_NETWORK_ATP_T0('NETWORK.Gpm')= 0.4;
    weights_NETWORK_ATP_T0('NETWORK.ATP')= 0.4;
    weights_NETWORK_ATP_T1('NETWORK.PFK1')= 0.4;
    weights_NETWORK_ATP_T1('NETWORK.GCK')= 0.2;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'NETWORK.ATP', 1, ...
        {'mean', 0.1, 'cov', 0.001}, weights_NETWORK_ATP_T0, weights_NETWORK_ATP_T1);  % NETWORK_ATP ~ 0.4* NETWORK.Gpm + 0.4*NETWORK.PFK1 + 0.2*NETWORK.GCK; 
    
    weights_NETWORK_cAMP_T0= containers.Map(); % parents in slice t
    weights_NETWORK_cAMP_T1= containers.Map(); % parents in slice t+1
    weights_NETWORK_cAMP_T0('NETWORK.ATP')= 0.4;
    weights_NETWORK_cAMP_T0('NETWORK.cAMP')= 0.4;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'NETWORK.cAMP', 1, ...
        {'mean', 0.1, 'cov', 0.001}, weights_NETWORK_cAMP_T0, weights_NETWORK_cAMP_T1);

    weights_NETWORK_Calcium_T0= containers.Map(); % parents in slice t
    weights_NETWORK_Calcium_T1= containers.Map(); % parents in slice t+1
    weights_NETWORK_Calcium_T0('NETWORK.cAMP')= 0.4;
    weights_NETWORK_Calcium_T0('NETWORK.Calcium')= 0.4;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'NETWORK.Calcium', 1, ...
        {'mean', 0.1, 'cov', 0.001}, weights_NETWORK_Calcium_T0, weights_NETWORK_Calcium_T1);

    weights_NETWORK_Ipm_T0= containers.Map(); % parents in slice t
    weights_NETWORK_Ipm_T1= containers.Map(); % parents in slice t+1
    weights_NETWORK_Ipm_T0('NETWORK.Calcium')= 0.4;
    weights_NETWORK_Ipm_T0('NETWORK.Ipm')= 0.4;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'NETWORK.Ipm', 1, ...
        {'mean', 0.1, 'cov', 0.001}, weights_NETWORK_Ipm_T0, weights_NETWORK_Ipm_T1);
   
    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

