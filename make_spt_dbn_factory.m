% Make a DBN for the spatiotemporal BD model with the following variables
%
% Time-dependent variables
%  -> SPT.Gisg(t)  ->  SPT.Gisg(t+1) ->
%  -> SPT.Ipm(t)  ->  SPT.Ipm(t+1) ->
%
% Reference variables
% R.Gisg(t), R.Gisg(t+1)
% R.Ipm(t), R.Ipm(t+1)
%
% Observed variables
% E.Gisg(t), E.Gisg(t+1)
% E.Ipm(t), E.Ipm(t+1)
%
% Time-invariant variables
% SPT.lambda SPT.k SPT.Npatch SPT.Nisg SPT.Rpbc
%
% All variables display gaussian distributions.

% To generate a conditional gaussian model
function [dbn_factory]= make_spt_dbn_factory(Go, Io, time)
    % Define nodes and intra-slice and inter-slice edges between them 
    node_names=  { 'SPT.lambda', 'SPT.k', 'SPT.Npatch' ,'SPT.Nisg', 'SPT.Rpbc', 'SPT.Gisg', 'SPT.Ipm', 'R.Gisg', 'E.Gisg','R.Ipm', 'E.Ipm'};
    % Intra - in one time slice
    edges_intra= {'SPT.Gisg', 'R.Gisg'; 'R.Gisg', 'E.Gisg'; 'SPT.k', 'SPT.Ipm'; 'SPT.Npatch', 'SPT.Ipm';'SPT.Nisg', 'SPT.Ipm'; 'SPT.Rpbc', 'SPT.Ipm'; 'SPT.Ipm', 'R.Ipm'; 'R.Ipm', 'E.Ipm'};
    % Inter - between time slices
    edges_inter= { 'SPT.Gisg', 'SPT.Gisg'; 'SPT.Gisg', 'SPT.Ipm'; 'SPT.Ipm', 'SPT.Ipm'; 'SPT.k', 'SPT.k' };
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
    eclass2_map('SPT.k')= 'SPT.k.inter';
    eclass2_map('SPT.Gisg')= 'SPT.Gisg.inter';
    eclass2_map('SPT.Ipm')= 'SPT.Ipm.inter';   
    
    % Specify distributions for CPDs, mu is mean, Sigma is cov,  W is
    % weights
    % - no parents: Y ~ N(mu, Sigma)
    % - cts parents : Y|X=x ~ N(mu + W x, Sigma)
    % - discrete parents: Y|Q=i ~ N(mu(:,i), Sigma(:,:,i))
    % - cts and discrete parents: Y|X=x,Q=i ~ N(mu(:,i) + W(:,:,i) * x, Sigma(:,:,i))
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    
    CPDFactories= {};
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'SPT.lambda', 0, ...
        {'mean', 0.1, 'cov', 1e-4});
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'SPT.k', 0, ...
        {'mean', 10, 'cov', 1e-2});
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'SPT.Npatch', 0, ...
        {'mean', 6, 'cov', 1e-2});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.Nisg', 0, ...
        {'mean', 300, 'cov', 1e-2});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.Rpbc', 0, ...
        {'mean', 4, 'cov', 1e-2});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.Gisg', 0, ...
        {'mean', Go(1), 'cov', 0.1});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'R.Gisg', 0, ...
        {'mean', 0.0, 'cov', 0.1,'weights', 1.0}); % E= [1.0*G(t)] 
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'E.Gisg', 0, ...
        {'mean', 0.0, 'cov', 0.1, 'weights', 1.0}); % E= [1.0*Gref(t)]
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.Ipm', 0, ...
        {'mean', Io(1), 'cov', 0.5, 'weights', [1.0 1.0 1.0 1.0]});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'R.Ipm', 0, ...
        {'mean', 0.0, 'cov', 0.5, 'weights', 1.0}); % E= [1.0*I(t)] 
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'E.Ipm', 0, ...
        {'mean', 0.0, 'cov', 0.5, 'weights', 1.0}); % E= [1.0*Iref(t)] 
  
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    
    % CPD for G(t+1)
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'SPT.k', 1, ...
        {'mean', 0.0, 'cov', 1e-2,  'weights', 1.0} ); % E[ k(t+1) ] = 1.0*k(t) +-  sqrt(1E-8)    

    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.Gisg', 1, ...
        {'mean', 0.0, 'cov', 0.1, 'weights', 1.0});
    
    % CPD for I(t+1), assume for now all parents are continuous
    % TODO: The following simplified structure will need to be updated according to the SPT model.
    % I(t+1) := Normal dist. E = weightA * I(t) + weightB * G(t) + weightC
    % * k + weightD * Npatch + weightE  * Nisg + weightF  * Rpbc
    
    weights_I1_map_T0= containers.Map(); % parents in slice t
    weights_I1_map_T1= containers.Map(); % parents in slice t+1
    weights_I1_map_T0('SPT.Ipm')= 0.1; % weightA
    weights_I1_map_T0('SPT.Gisg')= 0.1; % weightB
    weights_I1_map_T1('SPT.k')= 0.8; % weightC
    weights_I1_map_T1('SPT.Npatch')= 0.1; % weightD
    weights_I1_map_T1('SPT.Nisg')= 0.2; % weightE
    weights_I1_map_T1('SPT.Rpbc')= 0.2; % weightF
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.Ipm', 1, ...
        {'mean', 0.0, 'cov', 0.5}, ...
        weights_I1_map_T0, weights_I1_map_T1);
    
    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
end    

