% Make a DBN for the meal model with the following variables
%
% Time-dependent variables
%  -> G.Meal(t)  ->  G.Meal(t+1) ->
%  -> I.Meal(t)  ->  I.Meal(t+1) ->
%
% Reference variables
% G.ref(t), G.ref(t+1)
% I.ref(t), I.ref(t+1)
%
% Observed variables
% G.obs(t), G.obs(t+1)
% I.obs(t), I.obs(t+1)
%
% Time-invariant variables
% Gb.Meal
%
% Parameters
% ALPHA BETA
%
% To generate a conditional gaussian model
function [dbn_factory]= make_meal_dbn_factory(G_Model, I_Model, alpha_Model,beta_Model, Gb_Model)
    node_names=  {'G.Meal','I.Meal','G.ref','I.ref','G.obs','I.obs','Gb.Meal','G_minus_Gb.Meal'}; % BARAK comment: removed alpha, beta (turned to weights) + added intermediate (G-h)
    n= length(node_names);
    % Intra - in one time slice
    edges_intra= {'Gb.Meal','G_minus_Gb.Meal'; 'G.Meal','G_minus_Gb.Meal'; 'G.Meal', 'G.ref';'G.ref','G.obs';'I.Meal','I.ref';'I.ref','I.obs'};
    % Inter - between time slices
    edges_inter= { 'G.Meal', 'G.Meal'; 'Gb.Meal', 'Gb.Meal';'G_minus_Gb.Meal','G_minus_Gb.Meal'; 'G_minus_Gb.Meal','I.Meal';'I.Meal', 'I.Meal' }; % BARAK comment: switched G->I and h->I to (G-h)->I
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('Gb.Meal')= 'Gb.Meal.inter';
    eclass2_map('G.Meal')= 'G.Meal.inter';
    eclass2_map('G_minus_Gb.Meal')= 'G.minus.Gb.Meal.inter';
    eclass2_map('I.Meal')= 'I.Meal.inter';   
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Gb.Meal', 0, ...
        {'mean', Gb_Model,   'cov', 0.001} );
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'G.Meal', 0, ...
        {'mean', G_Model, 'cov', 0.1} ); 
    weights_G_minus_h_map_T0= containers.Map(); % parents in slice t
    weights_G_minus_h_map_T1= containers.Map(); % parents in slice t+1
    weights_G_minus_h_map_T0('G.Meal')= 1.0;
    weights_G_minus_h_map_T0('Gb.Meal')= -1.0;
    % When using clamp, the root node is clamped to the N(0,I) distribution, so that we will not update these parameters during learning. 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'G_minus_Gb.Meal', 0, ...
        {'mean', 0.0, 'cov', 0.1}, ...
        weights_G_minus_h_map_T0, weights_G_minus_h_map_T1); % G_minus_h ~ Norm(E(G)-E(h))
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'G.ref', 0,   ...
        {'mean', 0.0,'cov', 0.1, 'weights', 1.0} ); % E= [1.0*G(t)] 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'G.obs', 0, ...
        {'mean', 0.0,'cov', 0.1, 'weights', 1.0}); % E= [1.0*Gref(t)]
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'I.Meal', 0, ...
        {'mean', I_Model, 'cov', 0.5} ); 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'I.ref', 0, ...
        { 'mean', 0.0,'cov', 0.5,   'weights', 1.0} ); % E= [1.0*I(t)] 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'I.obs', 0, ...   
        {'mean', 0.0,'cov', 0.5,   'weights', 1.0} ); % E= [1.0*Iref(t)]
    
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    
    % CPD for h(t+1), forcing h to be nearly constant
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'Gb.Meal', 1, ...
        {'mean',0.0,'cov', 0.5, 'weights', 1.0} ); % E[Gb(t+1) ] = 1.0*Gb(t) +- 0.2    

    % CPD for G(t+1)
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'G.Meal', 1, ...
        {'mean',0.0,'cov', 0.5, 'weights', 1.0} ); % E[ G(t+1) ] = 1.0*G(t) +- 0.1
    % I(t+1) := Normal dist. E = (1-alpha) * I(t) + alpha * beta * (G(t)-h)
   
    weights_G_minus_h_map_T0= containers.Map(); % parents in slice t
    weights_G_minus_h_map_T1= containers.Map(); % parents in slice t+1
    weights_G_minus_h_map_T0('G_minus_Gb.Meal')= 0.0;
    weights_G_minus_h_map_T1('G.Meal')= 1.0;
    weights_G_minus_h_map_T1('Gb.Meal')= -1.0;
    % When using clamp, the root node is clamped to the N(0,I) distribution, so that we will not update these parameters during learning. 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'G_minus_Gb.Meal', 1, ...
        {'mean', 0.0, 'cov', 0.1}, ...
        weights_G_minus_h_map_T0, weights_G_minus_h_map_T1); % G_minus_h ~ Norm(E(G)-E(h))
    
    weights_I1_map_T0= containers.Map(); % parents in slice t
    weights_I1_map_T1= containers.Map(); % parents in slice t+1
    INITIAL_ALPHA= alpha_Model;
    INITIAL_BETA= beta_Model;
    weights_I1_map_T0('I.Meal')= 1.0 - 2 * INITIAL_ALPHA;% assume dt = 2 min
    weights_I1_map_T0('G_minus_Gb.Meal')= 2* INITIAL_ALPHA * INITIAL_BETA * 1000;
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'I.Meal', 1, ...
        {'mean',0.0,'cov', 0.5}, ...
        weights_I1_map_T0, weights_I1_map_T1);
    
    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

