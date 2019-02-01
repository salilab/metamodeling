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
function [bnet, nodes_map, node_names, edges_intra, edges_inter, ns, eclass1, eclass2]= make_ode_bnet(Gm, Im)
    % Define nodes and intra-slice and inter-slice edges between them 
    node_names=  {'alpha', 'beta', 'h' ,'G',  'I', 'Gref', 'Gexp','Iref', 'Iexp'};
    % Intra - in one time slice
    edges_intra= {'alpha', 'I'; 'beta', 'I'; 'h', 'I'; 'G', 'Gref'; 'I', 'Iref'; 'Gref', 'Gexp'; 'Iref', 'Iexp'};
    % Inter - between time slices
    edges_inter= { 'G', 'G'; 'G', 'I'; 'I', 'I' };
    [intra, inter, nodes_map]= get_valid_nodes_graph(node_names, edges_intra, edges_inter);
    n= length(node_names);
    ns = ones(1, n);% all cts nodes are scalar
    %dnodes= [ ]; % descrete nodes
    %cnodes = mysetdiff(1:n, dnodes); % all are continuous nodes except for dnodes
    %onodes= cnodes; % observed nodes
    % "Equivalence classes" specify how the template is initiated and rolled
    % Specify which CPD is associates with each node in either time
    % slice 1 (eclass1) or in slice 2 onwards (eclass2)
    disp("check")
    disp(nodes_map('alpha'));
    
    eclass1_map= containers.Map();
    eclass1_map('alpha')=1;
    eclass1_map('beta')=2;
    eclass1_map('h')= 3;
    eclass1_map('G')= 4;
    eclass1_map('Gref')= 5;
    eclass1_map('Gexp')= 6;
    eclass1_map('I')= 7;
    eclass1_map('Iref')= 8;
    eclass1_map('Iexp')= 9;
    
    eclass2_map= containers.Map();
    eclass2_map('alpha')= 1;
    eclass2_map('beta')= 2;
    eclass2_map('h')= 3;
    eclass2_map('G')= 10;
    eclass2_map('Gref')= 11;
    eclass2_map('Gexp')= 12;
    eclass2_map('I')= 13;
    eclass2_map('Iref')= 14;
    eclass2_map('Iexp')= 15;
    
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
    
    bnet.CPD{1} = gaussian_CPD(bnet, nodes_map('alpha'),   'mean', 0.05, 'cov', 0.001);
    bnet.CPD{2} = gaussian_CPD(bnet, nodes_map('beta'),   'mean', 0.11, 'cov', 0.002);
    bnet.CPD{3} = gaussian_CPD(bnet, nodes_map('h'),   'mean', 6.1, 'cov', 0.1);
    
    bnet.CPD{4} = gaussian_CPD(bnet, nodes_map('G'),   'mean', Gm(1), 'cov', 0.1);
    bnet.CPD{5} = gaussian_CPD(bnet, nodes_map('Gref'),   'mean', Gm(1), 'cov', 0.1, 'weights', 1);
    bnet.CPD{6} = gaussian_CPD(bnet, nodes_map('Gexp'), 'mean', Gm(1), 'cov', 0.1, 'weights', 1);
    
    bnet.CPD{7} = gaussian_CPD(bnet, nodes_map('I'),   'mean', Im(1), 'cov', 5);
    bnet.CPD{8} = gaussian_CPD(bnet, nodes_map('Iref'),   'mean', Im(1), 'cov', 5, 'weights', 1);
    bnet.CPD{9} = gaussian_CPD(bnet, nodes_map('Iexp'), 'mean', Im(1), 'cov', 5, 'weights', 1);
    
    % eclass2
    weights_G= 0.5;
    bnet.CPD{10} = gaussian_CPD(bnet, nodes_map('G')+n, 'mean', Gm(n), 'cov', 0.1,  'weights', weights_G);
    
    weights_Gref= 1;n_CPD(bnet, nodes_map('Gref')+n, 'mean', Gm(n), 'cov', 0.1,  'weights', weights_Gref);
    
    weights_Gexp= 1;
    bnet.CPD{12} = gaussian_CPD(bnet, nodes_map('Gexp')+n, 'mean', Gm(n), 'cov', 0.1,  'weights', weights_Gexp);
    
    weights_I=zeros(5,1);
    bnet.CPD{11} = gaussia
    disp(nodes_map('alpha'));
    disp(nodes_map('beta'));
    disp(nodes_map('h'));
    disp(nodes_map('G'));
    disp(nodes_map('I'));
    weights_I(nodes_map('alpha'))= 0.6;
    weights_I(nodes_map('beta'))= 0.6;
    weights_I(nodes_map('h'))= 0.6;
    weights_I(nodes_map('G'))= 0.05;
    weights_I(nodes_map('I'))= 0.6;
    bnet.CPD{13} = gaussian_CPD(bnet, nodes_map('I')+n,   'mean', Im(n), 'cov', 5, 'weights', weights_I);
    %bnet.CPD{16} = gaussian_CPD(bnet, nodes_map('I')+n,   'mean', 72, 'cov', 5);
    
    weights_Iref= 1;
    bnet.CPD{14} = gaussian_CPD(bnet, nodes_map('Iref')+n, 'mean', Im(n), 'cov', 5,  'weights', weights_Iref);
    
    weights_Iexp= 1;
    bnet.CPD{15} = gaussian_CPD(bnet, nodes_map('Iexp')+n, 'mean', Im(n), 'cov', 5,  'weights', weights_Iexp);
end

