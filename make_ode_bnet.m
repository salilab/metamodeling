% Make a DBN for the ODE model with the following variables
%
% Time-dependent variables
%  -> G(t)  ->  G(t+1) ->
%  -> I(t)  ->  I(t+1) ->
%
% Reference variables
% Gref(t), Gref(t+1)
% Iref(t), Iref(t+1)
%
% Observed variables
% Gexp(t), Gexp(t+1)
% Iexp(t), Iexp(t+1)
%
% Time-invariant variables
% alpha beta h 
%

% To generate a conditional gaussian model
function [bnet, node_names, edges_intra, edges_inter, eclass1_map, eclass2_map, CPDFactories]= make_ode_bnet(Gm, Im, time)
    node_names=  horzcat(strcat('ODE.', {'h', 'G', 'G_minus_h', 'I', 'Gref', 'Gexp', 'Iexp'}), 'Reference.I'); % BARAK comment: removed alpha, beta (turned to weights) + added intermediate (G-h)
    n= length(node_names);
    ns = ones(1, n);% all cts nodes are scalar
    % Intra - in one time slice
    edges_intra1= strcat('ODE.', {'h', 'G_minus_h'; 'G', 'G_minus_h'; 'G', 'Gref'; 'Gref', 'Gexp'}); % BARAK comment: changed in accordance with change in nodes list
    edges_intra2= {'ODE.I','Reference.I'; 'Reference.I','ODE.Iexp'};
    edges_intra= [edges_intra1; edges_intra2];
    % Inter - between time slices
    edges_inter= strcat('ODE.', { 'G', 'G'; 'G_minus_h', 'I'; 'I', 'I' }); % BARAK comment: switched G->I and h->I to (G-h)->I
    [intra, inter, nodes_map, ~]= get_valid_nodes_graph(node_names, edges_intra, edges_inter);
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('ODE.G')= 'ODE.G.inter';
    eclass2_map('ODE.I')= 'ODE.I.inter';   

    %%
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    %%
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ODE.h', 0, ...
        {'mean', 6.1,   'cov', 0.1} );
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ODE.G', 0, ...
        {'mean', Gm(1), 'cov', 0.1} ); 
    weights_G_minus_h_map_T0= containers.Map(); % parents in slice t
    weights_G_minus_h_map_T1= containers.Map(); % parents in slice t+1
    weights_G_minus_h_map_T0('ODE.G')= 1.0;
    weights_G_minus_h_map_T0('ODE.h')= -1.0;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ODE.G_minus_h', 0, ...
        {'mean', 0.0, 'cov', 0.00001, 'clamp_mean', 1, 'clamp_cov', 1, 'clamp_weights', 1}, ...
        weights_G_minus_h_map_T0, weights_G_minus_h_map_T1);  % G_minus_h ~ Norm(E(G)-E(h))
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.Gref', 0,   ...
        {'mean', 0.0,   'cov', 0.2, 'weights', 1.0} ); % E= [1.0*G(t+1)] 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.Gexp', 0, ...
        {'mean', 0.0,   'cov', 0.1, 'weights', 1.0}); % E= [1.0*Gref(t+1)]
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.I', 0, ...
        {'mean', Im(1), 'cov', 5.0} ); 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Reference.I', 0, ...
        { 'mean', 0.0,   'cov', 10.0,   'weights', 1.0} ); 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.Iexp', 0, ...   
        {'mean', 0.0,   'cov', 5.0,   'weights', 1.0} ); 
    
    %%
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    %%
    % CPD for G(t+1)
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.G', 1, ...
        {'mean', 0.0, 'cov', 0.1,  'weights', 1.0} ); % E[ G(t+1) ] = 1.0*G(t)
    % CPD for I(t+1), assume for now all parents are continuous
    % I(t+1) := Normal dist. E = (1-alpha) * I(t) + alpha * beta * (G(t)-h)
    weights_I1_map_T0= containers.Map(); % parents in slice t
    weights_I1_map_T1= containers.Map(); % parents in slice t+1
    INITIAL_ALPHA= 0.05;
    INITIAL_BETA= 0.11;
    weights_I1_map_T0('ODE.I')= 1.0 - INITIAL_ALPHA;
    weights_I1_map_T0('ODE.G_minus_h')= INITIAL_ALPHA * INITIAL_BETA;
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.I', 1, ...
        {'mean', 0.0, 'cov', 5.0}, ...
        weights_I1_map_T0, weights_I1_map_T1);
    bnet= get_dynamic_bnet_from_maps(node_names, edges_intra, edges_inter, eclass1_map, eclass2_map, CPDFactories);
 end

