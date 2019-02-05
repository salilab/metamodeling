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
function [bnet, nodes_map, node_names, edges_intra, edges_inter, ns, eclass1, eclass2]= make_ode_bnet(Gm, Im, time)
    % Define nodes and intra-slice and inter-slice edges between them 
    node_names=  horzcat(strcat('ODE.', {'h', 'G', 'G_minus_h', 'I', 'Gref', 'Gexp', 'Iexp'}), 'Reference.I'); % BARAK comment: removed alpha, beta (turned to weights) + added intermediate (G-h)
    % Intra - in one time slice
    edges_intra1= strcat('ODE.', {'h', 'G_minus_h'; 'G', 'G_minus_h'; 'G', 'Gref'; 'Gref', 'Gexp'}); % BARAK comment: changed in accordance with change in nodes list
    edges_intra2= {'ODE.I','Reference.I'; 'Reference.I','ODE.Iexp'};
    edges_intra= [edges_intra1; edges_intra2];
    % Inter - between time slices
    edges_inter= strcat('ODE.', { 'G', 'G'; 'G_minus_h', 'I'; 'I', 'I' }); % BARAK comment: switched G->I and h->I to (G-h)->I
    [intra, inter, nodes_map, reverse_nodes_map]= get_valid_nodes_graph(node_names, edges_intra, edges_inter);
    n= length(node_names);
    ns = ones(1, n);% all cts nodes are scalar
    %dnodes= [ ]; % descrete nodes
    %cnodes = mysetdiff(1:n, dnodes); % all are continuous nodes except for dnodes
    %onodes= cnodes; % observed nodes
    % "Equivalence classes" specify how the template is initiated and rolled
    % Specify which CPD is associates with each node in either time
    % slice 1 (eclass1) or in slice 2 onwards (eclass2)
    %disp("check")
    %disp(nodes_map('ODE.alpha'));
    
    eclass1_map= containers.Map();
    eclass1_map('ODE.h')= 1;
    eclass1_map('ODE.G')= 2;
    eclass1_map('ODE.G_minus_h')= 3;
    eclass1_map('ODE.Gref')= 4;
    eclass1_map('ODE.Gexp')= 5;
    eclass1_map('ODE.I')= 6;
    eclass1_map('Reference.I')= 7;
    eclass1_map('ODE.Iexp')= 8;
    
    eclass2_map= containers.Map();
    eclass2_map('ODE.h')= 1;
    eclass2_map('ODE.G')= 9;
    eclass2_map('ODE.G_minus_h')= 3;
    eclass2_map('ODE.Gref')= 4;
    eclass2_map('ODE.Gexp')= 5;
    eclass2_map('ODE.I')= 10;
    eclass2_map('Reference.I')= 7;
    eclass2_map('ODE.Iexp')= 8;
    
    eclass1= get_eclass_from_maps(eclass1_map, nodes_map);
    eclass2= get_eclass_from_maps(eclass2_map, nodes_map);  

    
    % make the dbn
    bnet = mk_dbn(intra, inter, ns, ...
        'discrete', [], ...
        'eclass1', eclass1, ...
        'eclass2', eclass2);
    % Specify distributions for CPDs, mu is mean, Sigma is cov,  W is
    % weights
    % - no parents: Y ~ N(mu, Sigma)
    % - cts parents : Y|X=x ~ N(mu + W x, Sigma)
    % - discrete parents: Y|Q=i ~ N(mu(:,i), Sigma(:,:,i))
    % - cts and discrete parents: Y|X=x,Q=i ~ N(mu(:,i) + W(:,:,i) * x, Sigma(:,:,i))
    % Create gaussian CPDs for alpha, beta, and h, all with no parents.

    %%
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    %%
    bnet.CPD{1} = gaussian_CPD(bnet, nodes_map('ODE.h'),      'mean', 6.1,   'cov', 0.1);
    bnet.CPD{2} = gaussian_CPD(bnet, nodes_map('ODE.G'),      'mean', Gm(1), 'cov', 0.1); 
    parents_G_minus_h= parents(bnet.dag, nodes_map('ODE.G_minus_h')); % parents of G_minus_h(t)
    weights_G_minus_h_map_T0= containers.Map(); % parents in slice t
    weights_G_minus_h_map_T0('ODE.G')= 1.0;
    weights_G_minus_h_map_T0('ODE.h')= -1.0;
    weights_G_minus_h= zeros(length(parents_G_minus_h),1);
   for i=1:length(parents_G_minus_h) 
        parent_index= parents_G_minus_h(i);
        if (parent_index <= n) % parent in slice t
            parent_name= reverse_nodes_map(parent_index);
            weights_G_minus_h(i)= weights_G_minus_h_map_T0(parent_name);
        else % parent in slice t+1
            assert(boolean(0))
        end
    end

    bnet.CPD{3} = gaussian_CPD(bnet, nodes_map('ODE.G_minus_h'), 'mean', 0.0, 'cov', 0.00001, 'weights', weights_G_minus_h, 'clamp_mean', 1, 'clamp_cov', 1, 'clamp_weights', 1); 
    bnet.CPD{4} = gaussian_CPD(bnet, nodes_map('ODE.Gref'),   'mean', 0.0,   'cov', 0.2, 'weights', 1.0); % COMMENT BR: - changed mean from [Gm(1)+1.0*G(t+1)]    to [1.0*G(t+1)] ; increased variance to reflect the intuition that ODE is a less accurate model of Gref than Gexp (this is quite arbitraty)
    bnet.CPD{5} = gaussian_CPD(bnet, nodes_map('ODE.Gexp'),   'mean', 0.0,   'cov', 0.1, 'weights', 1.0);  % COMMENT BR: Changed mean from  [Gm(1)+1.0*Gref(t+1)] to [1.0*Gref(t+1)]
    bnet.CPD{6} = gaussian_CPD(bnet, nodes_map('ODE.I'),      'mean', Im(1), 'cov', 5.0); 
    bnet.CPD{7} = gaussian_CPD(bnet, nodes_map('Reference.I'),'mean', 0.0,   'cov', 10.0,   'weights', 1.0); % COMMENT BR: see comment above for Gref, Gexp
    bnet.CPD{8} = gaussian_CPD(bnet, nodes_map('ODE.Iexp'),   'mean', 0.0,   'cov', 5.0,   'weights', 1.0); % COMMENT BR: see comment above for Gref, Gexp
    
    %%
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    %%
    % CPD for G(t+1)
    bnet.CPD{9} = gaussian_CPD(bnet, nodes_map('ODE.G')+n, 'mean', 0.0, 'cov', 0.1,  'weights', 1.0); % COMMENT BR: changed mean of G(t+1) from 0.5*G(t)+Gm(time) to 1.0*G(t)
    % CPD for I(t+1), assume for now all parents are continuous
    parents_I1= parents(bnet.dag, nodes_map('ODE.I')+n); % parents of I(t+1)
    weights_I1_map_T0= containers.Map(); % parents in slice t
    weights_I1_map_T1= containers.Map(); % parents in slice t+1
    INITIAL_ALPHA= 0.05;
    INITIAL_BETA= 0.11;
    weights_I1_map_T0('ODE.I')= 1.0 - INITIAL_ALPHA;
    weights_I1_map_T0('ODE.G_minus_h')= INITIAL_ALPHA * INITIAL_BETA;
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
    bnet.CPD{10} = gaussian_CPD(bnet, nodes_map('ODE.I')+n,   'mean', 0.0, 'cov', 5.0, 'weights', weights_I1);
 end

