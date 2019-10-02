% Make a DBN for the spatiotemporal BD model with the following variables
%
% Time-dependent variables
%  -> Gin.SPT(t)  ->  Gin.SPT(t+1) ->
%  -> S.SPT(t)  ->  S.SPT(t+1) ->
%
% Reference variables
% Gin.ref(t), Gin.ref(t+1)
% S.ref(t), S.ref(t+1)
%
% Observed variables
% Gin.obs(t), Gin.obs(t+1)
% S.obs(t), S.obs(t+1)
%
% Time-invariant variables
% lambda.SPT k.SPT Npatch.SPT Nisg.SPT Disg.SPT Ninsulin.SPT Rpbc.SPT
%
% All variables display gaussian distributions.

% To generate a conditional gaussian model
function [dbn_factory]= make_spt_dbn_factory(Gin_SPT, S_SPT, k_SPT, Nisg_SPT, Npatch_SPT, Ninsulin_SPT, I_SPT)
    % Define nodes and intra-slice and inter-slice edges between them 
    node_names=  { 'lambda.SPT','Gin.ref','Gin.SPT', 'k.SPT', 'Npatch.SPT' ,...
        'Nisg.SPT', 'Disg.SPT', 'Ninsulin.SPT','Rpbc.SPT', 'S.SPT', 'I.SPT', 'Scell.ref'};
    % Intra - in one time slice
    edges_intra= { 'lambda.SPT','Gin.SPT';'Gin.ref','Gin.SPT';'k.SPT', 'S.SPT';...
        'Npatch.SPT', 'S.SPT'; 'Nisg.SPT', 'S.SPT'; 'Disg.SPT', 'S.SPT'; 'Ninsulin.SPT', 'S.SPT';...
        'Rpbc.SPT', 'S.SPT';'S.SPT', 'Scell.ref'};
    % Inter - between time slices
    edges_inter= { 'Gin.SPT', 'Gin.SPT'; 'Gin.SPT', 'S.SPT'; 'S.SPT', 'S.SPT';'k.SPT', 'k.SPT'; ...
        'Nisg.SPT', 'Nisg.SPT'; 'Npatch.SPT', 'Npatch.SPT'; 'Ninsulin.SPT', 'Ninsulin.SPT'; 'S.SPT', 'I.SPT'; 'I.SPT', 'I.SPT' };
    % "Equivalence classes" specify how the template is initiated and rolled
    % Specify which CPD is associates with each node in either time
    % slice 1 (eclass1) or in slice 2 onwards (eclass2)
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('Gin.SPT')= 'Gin.SPT.inter';
    eclass2_map('S.SPT')= 'S.SPT.inter';
    eclass2_map('I.SPT')= 'I.SPT.inter';
    eclass2_map('k.SPT')= 'k.SPT.inter'; 
    eclass2_map('Nisg.SPT')= 'Nisg.SPT.inter';   
    eclass2_map('Npatch.SPT')= 'Npatch.SPT.inter';   
    eclass2_map('Ninsulin.SPT')= 'Ninsulin.SPT.inter';   
    
    % Specify distributions for CPDs, mu is mean, Sigma is cov,  W is
    % weights
    % - no parents: Y ~ N(mu, Sigma)
    % - cts parents : Y|X=x ~ N(mu + W x, Sigma)
    % - discrete parents: Y|Q=i ~ N(mu(:,i), Sigma(:,:,i))
    % - cts and discrete parents: Y|X=x,Q=i ~ N(mu(:,i) + W(:,:,i) * x, Sigma(:,:,i))
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    
    CPDFactories= {};
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'lambda.SPT', 0, ...
        {'mean', 1.0, 'cov', 1E-4});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Gin.ref', 0, ...
        {'mean', Gin_SPT, 'cov', 1E-4});
    weights_Gin0_map_T0= containers.Map(); % parents in slice t
    weights_Gin0_map_T1= containers.Map(); % parents in slice t+1
    weights_Gin0_map_T0('lambda.SPT')= 0.0;
    weights_Gin0_map_T0('Gin.ref')= 1.0;
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Gin.SPT', 0, ...
        {'mean', 0.0, 'cov', 1E-4}, ...
        weights_Gin0_map_T0, weights_Gin0_map_T1);
    %CPDFactories{end+1} = ...         
    %    CPDFactory('Gaussian_CPD', 'Gin.obs', 0, ...
    %    {'mean', 0.0, 'cov', 1E-4,   'weights', 1.0});
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'k.SPT', 0, ...
        {'mean', k_SPT, 'cov', 0.005});
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Npatch.SPT', 0, ...
        {'mean', Npatch_SPT, 'cov', 0.1});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Nisg.SPT', 0, ...
        {'mean', Nisg_SPT, 'cov', 1});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Ninsulin.SPT', 0, ...
        {'mean', Ninsulin_SPT, 'cov', 1E-8});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Rpbc.SPT', 0, ...
        {'mean', 6.0, 'cov', 0.1});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Disg.SPT', 0, ...
        {'mean', 3.2e-3, 'cov', 1E-5});
    weights_S0_map_T0= containers.Map(); % parents in slice t
    weights_S0_map_T1= containers.Map(); % parents in slice t+1
    weights_S0_map_T0('k.SPT')= 0.0;
    weights_S0_map_T0('Npatch.SPT')= 0.0;
    weights_S0_map_T0('Ninsulin.SPT')= 0.0;
    weights_S0_map_T0('Disg.SPT')= 0.0;
    weights_S0_map_T0('Nisg.SPT')= 0.0;
    weights_S0_map_T0('Rpbc.SPT')= 0.0;
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'S.SPT', 0, ...
        {'mean', S_SPT, 'cov', 1E-4}, ...
        weights_S0_map_T0, weights_S0_map_T1);
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'I.SPT', 0, ...
        { 'mean', I_SPT,'cov', 1E-4} ); % I
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Scell.ref', 0, ...
        { 'mean', 0.0,'cov', 1E-4,   'weights', 1.0} ); % E= [1.0*I(t)] 
  
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    
    % CPD for G(t+1)
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'k.SPT', 1, ...
        {'mean', 0.0, 'cov', 0.005,  'weights', 1.0} ); % E[ k(t+1) ] = 1.0*k(t) +-  sqrt(1E-8)    
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'Npatch.SPT', 1, ...
        {'mean', 0.0, 'cov', 0.1,  'weights', 1.0} ); % E[ k(t+1) ] = 1.0*k(t) +-  sqrt(1E-8)   
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'Nisg.SPT', 1, ...
        {'mean', 0.0, 'cov', 1,  'weights', 1.0} ); % E[ k(t+1) ] = 1.0*k(t) +-  sqrt(1E-8)   
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'Ninsulin.SPT', 1, ...
        {'mean', 0.0, 'cov', 1e-8,  'weights', 1.0} ); % E[ k(t+1) ] = 1.0*k(t) +-  sqrt(1E-8)   
    weights_Gin1_map_T0= containers.Map(); % parents in slice t
    weights_Gin1_map_T1= containers.Map(); % parents in slice t+1
    weights_Gin1_map_T0('Gin.SPT')= 0.0; % weightA
    weights_Gin1_map_T1('lambda.SPT')= 0.0; % weightA
    weights_Gin1_map_T1('Gin.ref')= 1.0; % weightB
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Gin.SPT', 1, ...
        {'mean', 0.0, 'cov', 1E-4}, ...
        weights_Gin1_map_T0, weights_Gin1_map_T1);
    
    % CPD for I(t+1), assume for now all parents are continuous
    % TODO: The following simplified structure will need to be updated according to the SPT model.
    % I(t+1) := Normal dist. E = weightA * I(t) + weightB * G(t) + weightC
    % * k + weightD * Npatch + weightE  * Nisg + weightF  * Rpbc
    
    weights_S1_map_T0= containers.Map(); % parents in slice t
    weights_S1_map_T1= containers.Map(); % parents in slice t+1
    weights_S1_map_T0('S.SPT')= 0.1; % weightA
    weights_S1_map_T0('Gin.SPT')= 18; % weightB
    weights_S1_map_T1('k.SPT')= 0.5; % weightC
    weights_S1_map_T1('Npatch.SPT')= 0.5; % weightD
    weights_S1_map_T1('Nisg.SPT')= 0.02; % weightE
    weights_S1_map_T1('Ninsulin.SPT')= 0.5; % weightE
    weights_S1_map_T1('Disg.SPT')= 1E-5; % weightF
    weights_S1_map_T1('Rpbc.SPT')= -5; % weightF
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'S.SPT', 1, ...
        {'mean', 0.0, 'cov', 1E-4}, ...
        weights_S1_map_T0, weights_S1_map_T1);

    weights_I1_T0= containers.Map(); % parents in slice t
    weights_I1_T1= containers.Map(); % parents in slice t+1
    weights_I1_T0('S.SPT')= 25/2/34;
    weights_I1_T0('I.SPT')= 0.5;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'I.SPT', 1, ...
        {'mean', 0.0, 'cov', 1E-4}, weights_I1_T0, weights_I1_T1); % I(t+1) = 0.5 * I(t) + 25/2/34 * S(t)
    
    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
end    

