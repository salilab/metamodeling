% Make a DBN for the spatiotemporal BD model with the following variables
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
% Gobs(t), Gobs(t+1)
% Iobs(t), Iobs(t+1)
%
% Time-invariant variables
% lamda k Npatch Nisg Rpbc
%
% All variables display gaussian distributions.

% To generate a conditional gaussian model
function [bnet, nodes_map, node_names, edges_intra, edges_inter, ns, eclass1, eclass2]= make_spt_bnet(Gm, Im, time)
    % Define nodes and intra-slice and inter-slice edges between them 
    node_names=  { 'SPT.lambda', 'SPT.k', 'SPT.Npatch' ,'SPT.Nisg', 'SPT.Rpbc', 'SPT.G', 'SPT.I', 'SPT.Gref', 'SPT.Gobs','Reference.I', 'SPT.Iobs'};
    % Intra - in one time slice
    edges_intra= {'SPT.lambda', 'SPT.G'; 'SPT.G', 'SPT.Gref'; 'SPT.Gref', 'SPT.Gobs'; 'SPT.k', 'SPT.I'; 'SPT.Npatch', 'SPT.I';'SPT.Nisg', 'SPT.I'; 'SPT.Rpbc', 'SPT.I'; 'SPT.I', 'Reference.I'; 'Reference.I', 'SPT.Iobs'};
    % Inter - between time slices
    edges_inter= { 'SPT.G', 'SPT.G'; 'SPT.G', 'SPT.I'; 'SPT.I', 'SPT.I' };
    [intra, inter, nodes_map, reverse_nodes_map]= get_valid_nodes_graph(node_names, edges_intra, edges_inter);
    n= length(node_names);
    %disp("size")
    %disp(n)
    ns = ones(1, n);% all cts nodes are scalar
    %dnodes= []; % descrete nodes
    %cnodes = mysetdiff(1:n, dnodes); % all are continuous nodes except for dnodes
    %onodes= cnodes; % observed nodes
    %ns(dnodes) = 1; % descrete nodes have one value, one simulation for each setup
    % "Equivalence classes" specify how the template is initiated and rolled
    % Specify which CPD is associates with each node in either time
    % slice 1 (eclass1) or in slice 2 onwards (eclass2)

    eclass1_map= containers.Map();
    eclass1_map('SPT.lambda')= 1;
    eclass1_map('SPT.k')= 2;
    eclass1_map('SPT.Npatch')= 3;
    eclass1_map('SPT.Nisg')= 4;
    eclass1_map('SPT.Rpbc')= 5;
    eclass1_map('SPT.G')= 6;
    eclass1_map('SPT.Gref')= 7;
    eclass1_map('SPT.Gobs')= 8;
    eclass1_map('SPT.I')= 9;
    eclass1_map('Reference.I')= 10;
    eclass1_map('SPT.Iobs')= 11;
    
    eclass2_map= containers.Map();
    eclass2_map('SPT.lambda')= 1;
    eclass2_map('SPT.k')= 2;
    eclass2_map('SPT.Npatch')= 3;
    eclass2_map('SPT.Nisg')= 4;
    eclass2_map('SPT.Rpbc')= 5;
    eclass2_map('SPT.G')= 12;
    eclass2_map('SPT.Gref')= 7;
    eclass2_map('SPT.Gobs')= 8;
    eclass2_map('SPT.I')= 13;
    eclass2_map('Reference.I')= 10;
    eclass2_map('SPT.Iobs')= 11;
    
    eclass1= get_eclass_from_maps(eclass1_map, nodes_map);
    eclass2= get_eclass_from_maps(eclass2_map, nodes_map);  
    
    % make the dbn
    bnet = mk_dbn(intra, inter, ns, ...
        'discrete', [ ], ...
        'eclass1', eclass1, ...
        'eclass2', eclass2);
    % Specify distributions for CPDs, mu is mean, Sigma is cov,  W is
    % weights
    % - no parents: Y ~ N(mu, Sigma)
    % - cts parents : Y|X=x ~ N(mu + W x, Sigma)
    % - discrete parents: Y|Q=i ~ N(mu(:,i), Sigma(:,:,i))
    % - cts and discrete parents: Y|X=x,Q=i ~ N(mu(:,i) + W(:,:,i) * x, Sigma(:,:,i))
    % Create gaussian CPDs for alpha, beta, and h, all with no parents.
    % elcass1
    bnet.CPD{1} = gaussian_CPD(bnet, nodes_map('SPT.lambda'), 'mean', 0.1, 'cov', 0.01);
    bnet.CPD{2} = gaussian_CPD(bnet, nodes_map('SPT.k'), 'mean', 10, 'cov', 1);
    bnet.CPD{3} = gaussian_CPD(bnet, nodes_map('SPT.Npatch'), 'mean', 6, 'cov', 0.5);
    bnet.CPD{4} = gaussian_CPD(bnet, nodes_map('SPT.Nisg'), 'mean', 300, 'cov', 20);
    bnet.CPD{5} = gaussian_CPD(bnet, nodes_map('SPT.Rpbc'), 'mean', 4, 'cov', 0.1);
    bnet.CPD{6} = gaussian_CPD(bnet, nodes_map('SPT.G'),   'mean', Gm(1), 'cov', 2);
    bnet.CPD{7} = gaussian_CPD(bnet, nodes_map('SPT.Gref'),   'mean', Gm(1), 'cov', 2,'weights', 1);
    bnet.CPD{8} = gaussian_CPD(bnet, nodes_map('SPT.Gobs'), 'mean', Gm(1), 'cov', 2, 'weights', 1);
    bnet.CPD{9} = gaussian_CPD(bnet, nodes_map('SPT.I'),   'mean', Im(1), 'cov', 5);
    bnet.CPD{10} = gaussian_CPD(bnet, nodes_map('Reference.I'),   'mean', Im(1), 'cov', 5, 'weights', 1);
    bnet.CPD{11} = gaussian_CPD(bnet, nodes_map('SPT.Iobs'), 'mean', Im(1), 'cov', 5, 'weights', 1);
  
    %eclass2
    weights_G= [0.5 0.5];
    bnet.CPD{12} = gaussian_CPD(bnet, nodes_map('SPT.G')+n, 'mean', Gm(time), 'cov', 2, 'weights', weights_G);
    
    % CPD for I(t+1), assume for now all parents are continuous
    parents_I1= parents(bnet.dag, nodes_map('SPT.I')+n); % parents of I(t+1)
    weights_I1_map_T0= containers.Map(); % parents in slice t
    weights_I1_map_T1= containers.Map(); % parents in slice t+1
    weights_I1_map_T0('SPT.G')= 0.05;
    weights_I1_map_T0('SPT.I')= 0.6;
    weights_I1_map_T1('SPT.k')= 0.6;
    weights_I1_map_T1('SPT.Npatch')= 0.6;
    weights_I1_map_T1('SPT.Nisg')= 0.6;
    weights_I1_map_T1('SPT.Rpbc')= 0.05;
    weights_I1= zeros(length(parents_I1),1);
    for i=1:length(parents_I1) 
        parent_index= parents_I1(i);
        if (parent_index <= n) % parent in slice t
            parent_name= reverse_nodes_map(parent_index);
            weights_I1(i)= weights_I1_map_T0(parent_name);
        else % parent in slice t+1
            parent_name= reverse_nodes_map(parent_index-n); % we only have a map for 1..9
            weights_I1(i)= weights_I1_map_T1(parent_name);
        end
    end
    bnet.CPD{13} = gaussian_CPD(bnet, nodes_map('SPT.I')+n, 'mean', Im(time), 'cov', 5, 'weights', weights_I1);
    %bnet.CPD{13} = gaussian_CPD(bnet, nodes_map('I')+n, 'mean', 70, 'cov', 5);
end    

