% Make a DBN for the spatiotemporal BD model with the following variables
%
% Time-dependent variables
%  -> SPT.G(t)  ->  SPT.G(t+1) ->
%  -> SPT.I(t)  ->  SPT.I(t+1) ->
%
% Reference variables
% SPT.Gref(t), SPT.Gref(t+1)
% SPT.Iref(t), SPT.Iref(t+1)
%
% Observed variables
% SPT.Gobs(t), SPT.Gobs(t+1)
% SPT.Iobs(t), SPT.Iobs(t+1)
%
% Time-invariant variables
% SPT.lambda SPT.k SPT.Npatch SPT.Nisg SPT.Rpbc
%
% All variables display gaussian distributions.

% To generate a conditional gaussian model
function [dbn_factory]= make_spt_dbn_factory(Go, Io, time)
    % Define nodes and intra-slice and inter-slice edges between them 
    node_names=  { 'SPT.lambda', 'SPT.k', 'SPT.Npatch' ,'SPT.Nisg', 'SPT.Rpbc', 'SPT.G', 'SPT.I', 'SPT.Gref', 'SPT.Gobs','Reference.I', 'SPT.Iobs'};
    % Intra - in one time slice
    edges_intra= {'SPT.G', 'SPT.Gref'; 'SPT.Gref', 'SPT.Gobs'; 'SPT.k', 'SPT.I'; 'SPT.Npatch', 'SPT.I';'SPT.Nisg', 'SPT.I'; 'SPT.Rpbc', 'SPT.I'; 'SPT.I', 'Reference.I'; 'Reference.I', 'SPT.Iobs'};
    % Inter - between time slices
    edges_inter= { 'SPT.G', 'SPT.G'; 'SPT.G', 'SPT.I'; 'SPT.I', 'SPT.I' };
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
    eclass2_map('SPT.G')= 'SPT.G.inter';
    eclass2_map('SPT.I')= 'SPT.I.inter';   
    
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
        {'mean', 0.1, 'cov', 0.01});
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'SPT.k', 0, ...
        {'mean', 10, 'cov', 1});
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'SPT.Npatch', 0, ...
        {'mean', 6, 'cov', 0.5});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.Nisg', 0, ...
        {'mean', 300, 'cov', 20});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.Rpbc', 0, ...
        {'mean', 4, 'cov', 0.1});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.G', 0, ...
        {'mean', Go(1), 'cov', 0.00001, 'clamp_mean', 1, 'clamp_cov', 1, 'clamp_weights', 1});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.Gref', 0, ...
        {'mean', Go(1), 'cov', 0.00001,'weights', 1}); % E= [1.0*G(t)] 
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.Gobs', 0, ...
        {'mean', Go(1), 'cov', 0.00001, 'weights', 1}); % E= [1.0*Gref(t)]
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.I', 0, ...
        {'mean', Io(1), 'cov', 5.0, 'weights', zeros(1,4)});
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Reference.I', 0, ...
        {'mean', Io(1), 'cov', 5.0, 'weights', 1}); % E= [1.0*I(t)] 
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.Iobs', 0, ...
        {'mean', Io(1), 'cov', 5.0, 'weights', 1}); % E= [1.0*Iref(t)] 
  
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    
    % CPD for G(t+1)
    
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.G', 1, ...
        {'mean', Go(time), 'cov', 2, 'weights', 1.0});
    
    % CPD for I(t+1), assume for now all parents are continuous
    % TODO: The following simplified structure will need to be updated according to the SPT model.
    % I(t+1) := Normal dist. E = weightA * I(t) + weightB * G(t) + weightC
    % * k + weightD * Npatch + weightE  * Nisg + weightF  * Rpbc
    
    weights_I1_map_T0= containers.Map(); % parents in slice t
    weights_I1_map_T1= containers.Map(); % parents in slice t+1
    weights_I1_map_T0('SPT.I')= 0.6; % weightA
    weights_I1_map_T0('SPT.G')= 0.05; % weightB
    weights_I1_map_T1('SPT.k')= 1; % weightC
    weights_I1_map_T1('SPT.Npatch')= 0; % weightD
    weights_I1_map_T1('SPT.Nisg')= 0; % weightE
    weights_I1_map_T1('SPT.Rpbc')= 0; % weightF
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'SPT.I', 1, ...
        {'mean', Io(time), 'cov', 5.0}, ...
        weights_I1_map_T0, weights_I1_map_T1);
    
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
end    

