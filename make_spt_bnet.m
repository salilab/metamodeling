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
function [bnet, nodes_map, node_names, edges_intra, edges_inter, ns, eclass1, eclass2]= make_spt_bnet(Gm, Im)
    % Define nodes and intra-slice and inter-slice edges between them 
    node_names=  { 'lamda', 'k', 'Npatch' ,'Nisg', 'Rpbc', 'G', 'I', 'Gref', 'Gobs','Iref', 'Iobs'};
    % Intra - in one time slice
    edges_intra= {'lamda', 'G'; 'G', 'Gref'; 'Gref', 'Gobs'; 'k', 'I'; 'Npatch', 'I';'Nisg', 'I'; 'Rpbc', 'I'; 'I', 'Iref'; 'Iref', 'Iobs'};
    % Inter - between time slices
    edges_inter= { 'G', 'G'; 'G', 'I'; 'I', 'I' };
    [intra, inter, nodes_map]= get_valid_nodes_graph(node_names, edges_intra, edges_inter);
    n= length(node_names);
    disp("size")
    disp(n)
    ns = ones(1, n);% all cts nodes are scalar
    %dnodes= []; % descrete nodes
    %cnodes = mysetdiff(1:n, dnodes); % all are continuous nodes except for dnodes
    %onodes= cnodes; % observed nodes
    %ns(dnodes) = 1; % descrete nodes have one value, one simulation for each setup
    % "Equivalence classes" specify how the template is initiated and rolled
    % Specify which CPD is associates with each node in either time
    % slice 1 (eclass1) or in slice 2 onwards (eclass2)

    eclass1 = zeros(n,1); 
    eclass1(nodes_map('lamda'))= 1;
    eclass1(nodes_map('k'))= 2;
    eclass1(nodes_map('Npatch'))= 3;
    eclass1(nodes_map('Nisg'))= 4;
    eclass1(nodes_map('Rpbc'))= 5;
    eclass1(nodes_map('G'))= 6;
    eclass1(nodes_map('Gref'))= 7;
    eclass1(nodes_map('Gobs'))= 8;
    eclass1(nodes_map('I'))= 9;
    eclass1(nodes_map('Iref'))= 10;
    eclass1(nodes_map('Iobs'))= 11;
    
    eclass2 = zeros(n,1); 
    eclass2(nodes_map('lamda'))= 1;
    eclass2(nodes_map('k'))= 2;
    eclass2(nodes_map('Npatch'))= 3;
    eclass2(nodes_map('Nisg'))= 4;
    eclass2(nodes_map('Rpbc'))= 5;
    eclass2(nodes_map('G'))= 12;
    eclass2(nodes_map('Gref'))= 13;
    eclass2(nodes_map('Gobs'))= 14;
    eclass2(nodes_map('I'))= 15;
    eclass2(nodes_map('Iref'))= 16;
    eclass2(nodes_map('Iobs'))= 17;
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
    bnet.CPD{1} = gaussian_CPD(bnet, nodes_map('lamda'), 'mean', 0.1, 'cov', 0.01);
    bnet.CPD{2} = gaussian_CPD(bnet, nodes_map('k'), 'mean', 10, 'cov', 1);
    bnet.CPD{3} = gaussian_CPD(bnet, nodes_map('Npatch'), 'mean', 6, 'cov', 0.5);
    bnet.CPD{4} = gaussian_CPD(bnet, nodes_map('Nisg'), 'mean', 300, 'cov', 20);
    bnet.CPD{5} = gaussian_CPD(bnet, nodes_map('Rpbc'), 'mean', 4, 'cov', 0.1);
    bnet.CPD{6} = gaussian_CPD(bnet, nodes_map('G'),   'mean', Gm(1), 'cov', 2);
    bnet.CPD{7} = gaussian_CPD(bnet, nodes_map('Gref'),   'mean', Gm(1), 'cov', 2, 'weights', 1);
    bnet.CPD{8} = gaussian_CPD(bnet, nodes_map('Gobs'), 'mean', Gm(1), 'cov', 2, 'weights', 1);
    bnet.CPD{9} = gaussian_CPD(bnet, nodes_map('I'),   'mean', Im(1), 'cov', 5);
    bnet.CPD{10} = gaussian_CPD(bnet, nodes_map('Iref'),   'mean', Im(1), 'cov', 5, 'weights', 1);
    bnet.CPD{11} = gaussian_CPD(bnet, nodes_map('Iobs'), 'mean', Im(1), 'cov', 5,'weights', 1);
  
    %eclass2
    weights_G= [0.5 0.5];
    bnet.CPD{12} = gaussian_CPD(bnet, nodes_map('G')+n, 'mean', Gm(n), 'cov', 2, 'weights', weights_G);
    
    weights_Gref= 1;
    bnet.CPD{13} = gaussian_CPD(bnet, nodes_map('Gref')+n, 'mean', Gm(n), 'cov', 2,  'weights', weights_Gref);
    
    weights_Gobs= 1;
    bnet.CPD{14} = gaussian_CPD(bnet, nodes_map('Gobs')+n, 'mean', Gm(n), 'cov', 2,  'weights', weights_Gobs);
    
    weights_I=zeros(6,1);
    disp(nodes_map('k'));
    disp(nodes_map('Npatch'));
    disp(nodes_map('Nisg'));
    disp(nodes_map('Rpbc'));
    disp(nodes_map('G'));
    disp(nodes_map('I'));
    weights_I(nodes_map('k')-1)= 0.6;
    weights_I(nodes_map('Npatch')-1)= 0.6;
    weights_I(nodes_map('Nisg')-1)= 0.6;
    weights_I(nodes_map('Rpbc')-1)= 0.05;
    weights_I(nodes_map('G')-1)= 0.6;
    weights_I(nodes_map('I')-1)= 0.6;
    bnet.CPD{15} = gaussian_CPD(bnet, nodes_map('I')+n, 'mean', Im(n), 'cov', 5, 'weights', weights_I);
    %bnet.CPD{15} = gaussian_CPD(bnet, nodes_map('I')+n, 'mean', 70, 'cov', 5);
    
    weights_Iref= 1;
    bnet.CPD{16} = gaussian_CPD(bnet, nodes_map('Iref')+n, 'mean', Im(n), 'cov', 5,  'weights', weights_Iref);
    
    weights_Iobs= 1;
    bnet.CPD{17} = gaussian_CPD(bnet, nodes_map('Iobs')+n, 'mean', Im(n), 'cov', 5,  'weights', weights_Iobs);
end

