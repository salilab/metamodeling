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
    node_names=  horzcat(strcat('ODE.', {'alpha', 'beta', 'h' ,'G',  'I', 'Gref', 'Gexp', 'Iexp'}), 'Reference.I');
    % Intra - in one time slice
    edges_intra1= strcat('ODE.', {'alpha', 'I'; 'beta', 'I'; 'h', 'I'; 'G', 'Gref'; 'Gref', 'Gexp'});
    edges_intra2= {'ODE.I','Reference.I'; 'Reference.I','ODE.Iexp'};
    edges_intra= [edges_intra1; edges_intra2];
    % Inter - between time slices
    edges_inter= strcat('ODE.', { 'G', 'G'; 'G', 'I'; 'I', 'I' });
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
    eclass1_map('ODE.alpha')=1;
    eclass1_map('ODE.beta')=2;
    eclass1_map('ODE.h')= 3;
    eclass1_map('ODE.G')= 4;
    eclass1_map('ODE.Gref')= 5;
    eclass1_map('ODE.Gexp')= 6;
    eclass1_map('ODE.I')= 7;
    eclass1_map('Reference.I')= 8;
    eclass1_map('ODE.Iexp')= 9;
    
    eclass2_map= containers.Map();
    eclass2_map('ODE.alpha')= 1;
    eclass2_map('ODE.beta')= 2;
    eclass2_map('ODE.h')= 3;
    eclass2_map('ODE.G')= 10;
    eclass2_map('ODE.Gref')= 5;
    eclass2_map('ODE.Gexp')= 6;
    eclass2_map('ODE.I')= 11;
    eclass2_map('Reference.I')= 8;
    eclass2_map('ODE.Iexp')= 9;
    
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
    % elcass1
    
    bnet.CPD{1} = gaussian_CPD(bnet, nodes_map('ODE.alpha'),   'mean', 0.05, 'cov', 0.001);
    bnet.CPD{2} = gaussian_CPD(bnet, nodes_map('ODE.beta'),   'mean', 0.11, 'cov', 0.002);
    bnet.CPD{3} = gaussian_CPD(bnet, nodes_map('ODE.h'),   'mean', 6.1, 'cov', 0.1);
    bnet.CPD{4} = gaussian_CPD(bnet, nodes_map('ODE.G'),   'mean', Gm(1), 'cov', 0.1);
    bnet.CPD{5} = gaussian_CPD(bnet, nodes_map('ODE.Gref'),   'mean', Gm(1), 'cov', 0.1,'weights', 1);
    bnet.CPD{6} = gaussian_CPD(bnet, nodes_map('ODE.Gexp'), 'mean', Gm(1), 'cov', 0.1, 'weights', 1);
    bnet.CPD{7} = gaussian_CPD(bnet, nodes_map('ODE.I'),   'mean', Im(1), 'cov', 5);
    bnet.CPD{8} = gaussian_CPD(bnet, nodes_map('Reference.I'),   'mean', Im(1), 'cov', 5, 'weights', 1);
    bnet.CPD{9} = gaussian_CPD(bnet, nodes_map('ODE.Iexp'), 'mean', Im(1), 'cov', 5, 'weights', 1);
    
    % eclass2
    bnet.CPD{10} = gaussian_CPD(bnet, nodes_map('ODE.G')+n, 'mean', Gm(time), 'cov', 0.1,  'weights', 0.5);
    
    % CPD for I(t+1), assume for now all parents are continuous
    parents_I1= parents(bnet.dag, nodes_map('ODE.I')+n); % parents of I(t+1)
    weights_I1_map_T0= containers.Map(); % parents in slice t
    weights_I1_map_T1= containers.Map(); % parents in slice t+1
    weights_I1_map_T0('ODE.G')= 0.05;
    weights_I1_map_T0('ODE.I')= 0.6;
    weights_I1_map_T1('ODE.alpha')= 0.6;
    weights_I1_map_T1('ODE.beta')= 0.6;
    weights_I1_map_T1('ODE.h')= 0.6;
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
    bnet.CPD{11} = gaussian_CPD(bnet, nodes_map('ODE.I')+n,   'mean', Im(time), 'cov', 5, 'weights', weights_I1);
    %bnet.CPD{16} = gaussian_CPD(bnet, nodes_map('I')+n,   'mean', 72, 'cov', 5);
    
    %disp('here');
    %disp(max(eclass2));
    %disp('here');
end

