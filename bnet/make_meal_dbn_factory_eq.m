
% Make a DBN for the meal model with the following variables
%
% Time-dependent variables
%  -> Gex.Meal(t)  ->  Gex.Meal(t+1) ->
%  -> Y.Meal(t)  ->  Y.Meal(t+1) ->
%
% Reference variables
% Gex.ref(t), Gex.ref(t+1)
% I.ref(t), I.ref(t+1)
%
% Observed variables
% Gex.obs(t), Gex.obs(t+1)
% I.obs(t), I.obs(t+1)
%
% Time-invariant variables
% Gb.Meal
%
% Parameters
% ALPHA BETA
%
% To generate a conditional gaussian model
function [dbn_factory]= make_meal_dbn_factory_eq(Gex_Meal, Y_Meal, alpha_Meal,beta_Meal, gamma_Meal, k1_Meal, k2_Meal, K_Meal, dt_Meal_min, Gb_Meal, DGintake_Meal, Sb_Meal, I_Meal)
    node_names = {'DGintake.Meal','DGintake.ref','DGintake.obs','Gex.Meal','Gex.ref','Gin.ref','Gb.Meal','G_minus_Gb.Meal','Y.Meal','Scell.ref',...
        'S.Meal','I.Meal','I.ref','I.obs'}; % BARAK comment: removed alpha, beta (turned to weights) + added intermediate (G-h)
    n= length(node_names);
    % Intra - in one time slice
    edges_intra= {'DGintake.Meal','DGintake.ref';'DGintake.ref','DGintake.obs';'Gb.Meal','G_minus_Gb.Meal'; 'Gex.Meal','G_minus_Gb.Meal'; 'Gex.Meal', 'Gex.ref';...
        'Gex.ref','Gin.ref'; 'DGintake.Meal','S.Meal';'Y.Meal','S.Meal';'Scell.ref','S.Meal';'I.Meal','I.ref';'I.ref','I.obs'};
    % Inter - between time slices
    edges_inter= { 'DGintake.Meal','DGintake.Meal';'DGintake.Meal','Gex.Meal';'I.Meal','Gex.Meal'; 'Gex.Meal', 'Gex.Meal'; 'Gb.Meal', 'Gb.Meal';'G_minus_Gb.Meal','G_minus_Gb.Meal'; ...
        'G_minus_Gb.Meal','Y.Meal';'Y.Meal', 'Y.Meal'; 'S.Meal','S.Meal';'S.Meal','I.Meal';'I.Meal','I.Meal' }; % BARAK comment: switched G->I and h->I to (G-h)->I
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('DGintake.Meal')= 'DGintake.Meal.inter';
    eclass2_map('Gb.Meal')= 'Gb.Meal.inter';
    eclass2_map('Gex.Meal')= 'Gex.Meal.inter';
    eclass2_map('G_minus_Gb.Meal')= 'Gex.minus.Gb.Meal.inter';
    eclass2_map('Y.Meal')= 'Y.Meal.inter';   
    eclass2_map('S.Meal')= 'S.Meal.inter';   
    eclass2_map('I.Meal')= 'I.Meal.inter';  
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    % When using clamp, the root node is clamped to the N(0,I) distribution, so that we will not update these parameters during learninGex. 
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'DGintake.Meal', 0, ...
        {'mean', DGintake_Meal,   'cov', 1E-9','clamp_cov', 10e-10} ); % DGintake
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'DGintake.ref', 0, ...
        {'mean', 0.0,   'cov', 1E-10, 'clamp_cov', 10e-10,'weights', 1.0} ); % DGintake.ref = 1.0 * DGintake
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'DGintake.obs', 0, ...
        {'mean', 0.0,   'cov', 1E-10,'clamp_cov', 10e-10, 'weights', 1.0} ); % DGintake.obs = 1.0 * DGintake.ref
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Gb.Meal', 0, ...
        {'mean', Gb_Meal,   'cov', 0.001} ); % Gb
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Gex.Meal', 0, ...
        {'mean', Gex_Meal, 'cov', 0.001} ); % Gex 
    weights_G_minus_h0_map_T0= containers.Map(); % parents in slice t
    weights_G_minus_h0_map_T1= containers.Map(); % parents in slice t+1
    weights_G_minus_h0_map_T0('Gex.Meal')= 1.0;
    weights_G_minus_h0_map_T0('Gb.Meal')= -1.0;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'G_minus_Gb.Meal', 0, ...
        {'mean', 0.0, 'cov', 0.001}, ...
        weights_G_minus_h0_map_T0, weights_G_minus_h0_map_T1); % G_minus_Gb = 1.0 * Gex - 1.0 * Gb
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Gex.ref', 0,   ...
        {'mean', 0.0,'cov', 0.001, 'weights', 1.0} ); % Gex.ref = 1.0 * Gex
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Gin.ref', 0, ...
        {'mean', 0.0,'cov', 0.001, 'weights', 0.5}); % Gex.obs = 0.5 * Gex.ref
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Y.Meal', 0, ...
        {'mean', Y_Meal, 'cov', 0.001} ); % Y
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Scell.ref', 0, ...
        {'mean', Sb_Meal, 'cov', 0.005} ); % Sb
    
    weights_S0_map_T0= containers.Map(); % parents in slice t
    weights_S0_map_T1= containers.Map(); % parents in slice t+1
    INITIAL_K = K_Meal;
    weights_S0_map_T0('DGintake.Meal')= 0.0;
    weights_S0_map_T0('Scell.ref')= 0.0;
    weights_S0_map_T0('Y.Meal')= 0.0;
    %weights_S0_map_T0('Span.ref')= 2/3;
    % When using clamp, the root node is clamped to the N(0,I) distribution, so that we will not update these parameters during learninGex. 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'S.Meal', 0, ...
        {'mean', Sb_Meal, 'cov', 0.001}, ...
        weights_S0_map_T0, weights_S0_map_T1); % S = 1.0 * Sb
    %,'clamp_mean', 1, 'clamp_cov', 1, 'clamp_weights', 1); % G_minus_h ~ Norm(E(G)-E(h))
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'I.Meal', 0, ...
        { 'mean', I_Meal,'cov', 0.001} ); % I
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'I.ref', 0, ...
        { 'mean', 0.0,'cov', 0.001,   'weights', 1.0} ); % I.ref = 1.0 * I
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'I.obs', 0, ...
        { 'mean', 0.0,'cov', 0.001,   'weights', 1.0} ); % I.obs = 1.0 * I.ref 
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'DGintake.Meal', 1, ...
        {'mean',DGintake_Meal,'cov', 10e-10, 'clamp_cov', 10e-10,'weights', 0.0} ); % Gintake(t+1) = 0.0 * Gintake(t+1)
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'Gb.Meal', 1, ...
        {'mean',0.0,'cov', 0.001, 'weights', 1.0} ); % Gb(t+1) = 1.0 * Gb(t)    
    % CPD for G(t+1)
    weights_Gex1_map_T0= containers.Map(); % parents in slice t
    weights_Gex1_map_T1= containers.Map(); % parents in slice t+1
    INITIAL_k1= k1_Meal;
    INITIAL_k2= k2_Meal;
    weights_Gex1_map_T0('DGintake.Meal')= dt_Meal_min;
    weights_Gex1_map_T0('I.Meal')= -INITIAL_k1*dt_Meal_min; % parents in slice t
    weights_Gex1_map_T0('Gex.Meal')= 1.0-INITIAL_k2*dt_Meal_min; % parents in slice t+1
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Gex.Meal', 1, ...
        {'mean',0.0,'cov', 0.001}, ...
        weights_Gex1_map_T0, weights_Gex1_map_T1); % Gex (t+1) = 1.0 * Gex(t) + 1.0 * DGintake(t+1) - INITIAL_k1 * dt_Meal_min * I(t) - INITIAL_k2 * dt_Meal_min * Gex(t)
    % I(t+1) := Normal dist. E = (1-alpha) * I(t) + alpha * beta * (G(t)-h
    weights_G_minus_h1_map_T0= containers.Map(); % parents in slice t
    weights_G_minus_h1_map_T1= containers.Map(); % parents in slice t+1
    weights_G_minus_h1_map_T0('G_minus_Gb.Meal')= 0.0;
    weights_G_minus_h1_map_T1('Gex.Meal')= 1.0;
    weights_G_minus_h1_map_T1('Gb.Meal')= -1.0;
    % When using clamp, the root node is clamped to the N(0,I) distribution, so that we will not update these parameters during learninGex. 
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'G_minus_Gb.Meal', 1, ...
        {'mean', 0.0, 'cov', 0.001}, ...
        weights_G_minus_h1_map_T0, weights_G_minus_h1_map_T1); % G_minus_Gb(t+1) = 1.0 * Gex(t+1) - 1.0 * Gb(t+1)
    
    weights_Y1_map_T0= containers.Map(); % parents in slice t
    weights_Y1_map_T1= containers.Map(); % parents in slice t+1
    INITIAL_ALPHA= alpha_Meal;
    INITIAL_BETA= beta_Meal;
    weights_Y1_map_T0('Y.Meal')= 1.0 - dt_Meal_min * INITIAL_ALPHA;% assume dt = 2 min
    weights_Y1_map_T0('G_minus_Gb.Meal')= dt_Meal_min* INITIAL_ALPHA * INITIAL_BETA;
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Y.Meal', 1, ...
        {'mean',0.0,'cov', 0.001}, ...
        weights_Y1_map_T0, weights_Y1_map_T1); % Y(t+1) = (1.0 - dt_Meal_min * INITIAL_ALPHA) * Y(t) + (dt_Meal_min* INITIAL_ALPHA * INITIAL_BETA) * G_minus_Gb
    
    weights_S1_map_T0= containers.Map(); % parents in slice t
    weights_S1_map_T1= containers.Map(); % parents in slice t+1
    INITIAL_K = K_Meal;
    weights_S1_map_T1('DGintake.Meal')= INITIAL_K;
    weights_S1_map_T1('Scell.ref')= 0.0;
    weights_S1_map_T1('Y.Meal')= 1.0;
    weights_S1_map_T0('S.Meal')= 0.0;
    %weights_S1_map_T1('Span.Meal')= 2/3;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'S.Meal', 1, ...
        {'mean', Sb_Meal, 'cov', 0.001}, ...
        weights_S1_map_T0, weights_S1_map_T1); % S(t+1) = Sb_Meal + 0.0 * S(t) + INITIAL_K * DGex(t+1) + 1.0 * Sb(t+1) + 1.0 * Y(t+1)
    
    weights_I1_map_T0= containers.Map(); % parents in slice t
    weights_I1_map_T1= containers.Map(); % parents in slice t+1
    INITIAL_GAMMA= gamma_Meal;
    weights_I1_map_T0('I.Meal')= 1 - INITIAL_GAMMA*dt_Meal_min;
    weights_I1_map_T0('S.Meal')= dt_Meal_min; % fast
    CPDFactories{end+1}= ...
        CPDFactory('Gaussian_CPD', 'I.Meal', 1, ...
        {'mean', 0.0, 'cov', 0.001}, ...
        weights_I1_map_T0, weights_I1_map_T1); % I(t+1) = (1 - INITIAL_GAMMA*dt_Meal_min) * I(t) + dt_Meal_min * S(t)
    
    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

