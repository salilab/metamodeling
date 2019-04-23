% Make a DBN for the ODE model with the following variables
%
% Time-dependent variables
%  -> ODE.Gpm(t)  ->  ODE.Gpm(t+1) ->
%  -> ODE.Ipm(t)  ->  ODE.Ipm(t+1) ->
%
% Reference variables
% R.Gpm(t), R.Gpm(t+1)
% R.Ipm(t), R.Ipm(t+1)
%
% Observed variables
% E.Gpm(t), E.Gpm(t+1)
% E.Ipm(t), E.Ipm(t+1)
%
% Time-invariant variables
% ODE.h 
%
% Parameters
% ALPHA BETA
%
% To generate a conditional gaussian model
function [dbn_factory]= make_ode_dbn_factory(Gm, Im, time)
    node_names=  {'ODE.Gpm','ODE.Ipm','R.Gpm','R.Ipm','E.Gpm','E.Ipm','ODE.h', 'ODE.alpha', 'ODE.beta'}; % BARAK comment: removed alpha, beta (turned to weights) + added intermediate (G-h)
    n= length(node_names);
    % Intra - in one time slice
    edges_intra= {'ODE.Gpm', 'R.Gpm';'R.Gpm','E.Gpm';'ODE.Ipm','R.Ipm';'R.Ipm','E.Ipm';'ODE.h','ODE.Ipm';'ODE.alpha','ODE.Ipm';'ODE.beta','ODE.Ipm'};
    % Inter - between time slices
    edges_inter= { 'ODE.Gpm', 'ODE.Gpm'; 'ODE.Gpm','ODE.Ipm';'ODE.Ipm' 'ODE.Ipm'; 'ODE.h','ODE.h'}; % BARAK comment: switched G->I and h->I to (G-h)->I
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('ODE.h')= 'ODE.h.inter';
    eclass2_map('ODE.Gpm')= 'ODE.Gpm.inter';
    eclass2_map('ODE.Ipm')= 'ODE.Ipm.inter';   
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ODE.h', 0, ...
        {'mean', 6.1,   'cov', 0.01} );
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ODE.alpha', 0, ...
        {'mean', 0.05,   'cov', 0.01} );
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ODE.beta', 0, ...
        {'mean', 0.11,   'cov', 0.01} );
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ODE.Gpm', 0, ...
        {'mean', Gm(1), 'cov', 0.1} ); 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'R.Gpm', 0,   ...
        {'mean', 0.0,   'cov', 0.1, 'weights', 1.0} ); % E= [1.0*G(t)] 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'E.Gpm', 0, ...
        {'mean', 0.0,   'cov', 0.1, 'weights', 1.0}); % E= [1.0*Gref(t)]
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.Ipm', 0, ...
        {'mean', Im(1), 'cov', 0.5,'weights', [1.0 1.0 1.0]} ); 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'R.Ipm', 0, ...
        { 'mean', 0.0,   'cov', 0.5,   'weights', 1.0} ); % E= [1.0*I(t)] 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'E.Ipm', 0, ...   
        {'mean', 0.0,   'cov', 0.5,   'weights', 1.0} ); % E= [1.0*Iref(t)]
    
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    
    % CPD for h(t+1), forcing h to be nearly constant
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'ODE.h', 1, ...
        {'mean', 0.0, 'cov', 0.01,  'weights', 1.0} ); % E[ h(t+1) ] = 1.0*h(t) +-  sqrt(1E-8)    

    % CPD for G(t+1)
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.Gpm', 1, ...
        {'mean', 0.0, 'cov', 0.1,  'weights', 1.0} ); % E[ G(t+1) ] = 1.0*G(t)
    % I(t+1) := Normal dist. E = (1-alpha) * I(t) + alpha * beta * (G(t)-h)
    weights_I1_map_T0= containers.Map(); % parents in slice t
    weights_I1_map_T1= containers.Map(); % parents in slice t+1
    weights_I1_map_T0('ODE.Ipm')= 1.0 ;
    weights_I1_map_T1('ODE.alpha')= 1.0;
    weights_I1_map_T1('ODE.beta')= 1.0;
    weights_I1_map_T0('ODE.Gpm')= 1.0;
    weights_I1_map_T1('ODE.h')= 1.0;
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.Ipm', 1, ...
        {'mean', 0.0, 'cov', 0.5}, ...
        weights_I1_map_T0, weights_I1_map_T1);
    
    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

