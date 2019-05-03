% Make a DBN with the following topology
%   X1    X1
%   ^     ^ 
%   |     |
%   X2 -> X2
%   |     |
%   v     v
%   X3    X3

p1= 0.25; % p1 as in manuscript - X2 swithching probability
p2= 0.8;  % p2 as in manuscript - X1/X3|X2 switching probability
[bnet, nodes_map]= make_main_bnet(p1, p2);
n= bnet.nnodes_per_slice;
T = 10; % lengthhs of sequences to explore

% Sample from the posterior:
%evidence=cell(n,T);
%evidence{2,1}=1;
seq=sample_dbn(bnet, 'length', T);%,'evidence', evidence);
order= cell2mat(values(nodes_map));
fprintf("Sampled time-series of length %d", T);
%sample_data= cell2mat(seq(order,:));

% Infer posterior probability given evidence:
fprintf("Infering posterioer probability giving evidence...\n")
%engine = bk_inf_engine(bnet, 'clusters', 'ff');
%engine = bk_inf_engine(bnet, 'clusters', 'exact');
engine = jtree_dbn_inf_engine(bnet); 
% I.
evidence= cell(n, T);
evidence{nodes_map('X2'), 1} = 2; 
marg= marginal_nodes(enter_evidence(engine, evidence), ...
                     nodes_map('X1'), ...
                     3);
fprintf("Posterior probability distribution of X1(time=3) given X2(time=1)==2 is:\n"); 
disp(marg.T)

% Parameter estimation from submodels and correct overall network topology
ncases = 200;
fprintf(strcat("Learning CPD parameters from %d instances of %d frames time-series of model 1\n", ...
    " and %d instances of model 2, with correct network topology\n"), ...
    ncases/2, T, ncases/2);
cases = cell(1, ncases);
onodes = [nodes_map('X1')]; 
%A = [1 2; 3 1; 2 3; 1 1; 1 2 ; 4 5; 6 3; 4 5; 1 1 ;1 2]
A = [1 2 2 1 1 2 1 1 2 2]
%A = 1:10

for i=1:ncases
  ev = sample_dbn(bnet, T);
  cases{i} = cell(n, T);
  %cases{i}(onodes,:) = ev(onodes, :);
  cases{i}(onodes,:) = num2cell(A(1:T));
end

disp(cases{1}(onodes,:))
% LLtrace is the learning curve: the vector of log-likelihood scores at each iteration.
[bnet2, ll] = learn_params_dbn_em(engine, cases, 'max_iter', 2);
disp(bnet.CPD)
for cpd_id=1:length(bnet.CPD)
    disp(cpd_id);
    fprintf("CPD %d - original parameters:\n", cpd_id);
    disp(CPD_to_CPT(bnet.CPD{cpd_id}));
    fprintf("CPD %d - trained parameters:\n", cpd_id);
    disp(CPD_to_CPT(bnet2.CPD{cpd_id,2}));
    disp(CPD_to_CPT(bnet2.CPD{cpd_id,3}));
    disp(CPD_to_CPT(bnet2.CPD{nodes_map('X2'),3}));
end

function [intra, inter, nodes_map]= get_valid_nodes_graph(node_names, edges_intra, edges_inter)
% Returns inter and intra adjacency matrices based on:
% node_names - a cell array of the nodes names
% edges_intra - intra-edges between nodes in the same time slice
% edges_inter - inter-edges between nodes in consecutive time slices
%
% The computed nodes indexes in inter and intra are topologically sorted within
% the same time slice, such that parents in intra always appear before
% their children
%
% Return:
% intra - intra-edge adjacency matrix such that parents always come before
%         children
% inter - inter-edge adjacency matrix in the same order as intra
% nodes_map - a map from node names to node indexes in intra/inter

    % Create intra and inter adjacency matrixes
    n= length(node_names);
    nodes_map= containers.Map(node_names, 1:n);
    intra= zeros(n);
    for row=1:size(edges_intra,1)
        from= edges_intra{row, 1};
        to= edges_intra{row, 2};
        intra(nodes_map(from), nodes_map(to))= 1;
    end
    inter= zeros(n);
    for row=1:size(edges_inter, 1)
        from= edges_inter{row, 1};
        to= edges_inter{row, 2};
        inter(nodes_map(from), nodes_map(to))= 1;
    end
    % Resort all in topological order
    order= topological_sort(intra);
    intra= intra(order, order);
    inter= inter(order, order);
    nodes_map= containers.Map(node_names(order), 1:3);
end


function [bnet, nodes_map]= make_main_bnet(p1, p2)
    % Define nodes and intra-slice and inter-slice edges between them 
    node_names=  { 'X1', 'X2', 'X3' };
    edges_intra= { 'X2', 'X1'; 
                   'X2', 'X3' };
    edges_inter= { 'X2', 'X2'; };
    [intra, inter, nodes_map]= get_valid_nodes_graph(node_names, edges_intra, edges_inter);
    n= length(node_names);
    node_sizes= 2*ones( 1, n ); % all binary
    dnodes= 1:n; % all are discrete nodes
    % "Equivalence classes" specify how the template is initiated and rolled
    % Specify which CPD is associates with each node in either time
    % slice 1 (eclass1) or in slice 2 onwards (eclass2)
    eclass1= zeros(n,1); 
    eclass1(nodes_map('X1'))= 3;
    eclass1(nodes_map('X2'))= 1;
    eclass1(nodes_map('X3'))= 3;
    eclass2= zeros(n,1);
    eclass2(nodes_map('X1'))= 3;
    eclass2(nodes_map('X2'))= 2;
    eclass2(nodes_map('X3'))= 3;
    onodes = dnodes; %[nodes_map('X3')]; 
    bnet = mk_dbn(intra, inter, node_sizes, ...
            'discrete', dnodes, ...
            'eclass1', eclass1, ...
            'eclass2', eclass2, ...
            'observed', onodes);
    % Assign CPDs to match equivalence classes
    % Specify distributions for CPDs
    prior_X2 = [0.5, 0.5];
    transmat_inter_X2_X2 = [ 1-p1,   p1; 
                             p1,     1-p1 ];
    transmat_intra_X2_X1 = [ 1-p2,   p2; 
                             p2,     1-p2 ];
    bnet.CPD{1} = tabular_CPD(bnet, nodes_map('X2'),   'CPT', prior_X2);
    bnet.CPD{2} = tabular_CPD(bnet, nodes_map('X2')+n, 'CPT', transmat_inter_X2_X2);
    bnet.CPD{3} = tabular_CPD(bnet, nodes_map('X1'),   'CPT', transmat_intra_X2_X1);
end   

function [bnet, nodes_map]= make_other_bnet()
    % Define nodes and intra-slice and inter-slice edges between them 
    node_names=  { 'X1', 'X2', 'X3' };
    edges_intra= { 'X1', 'X2'; 
                   'X1', 'X3' };
    edges_inter= { 'X1', 'X1'; };
    [intra, inter, nodes_map]= get_valid_nodes_graph(node_names, edges_intra, edges_inter);
    n= length(node_names);
    node_sizes= 2*ones( 1, n ); % all binary
    dnodes= 1:n; % all are discrete nodes
    eclass1= zeros(n,1); 
    eclass1(nodes_map('X1'))= 1;
    eclass1(nodes_map('X2'))= 3;
    eclass1(nodes_map('X3'))= 3;
    eclass2= zeros(n,1);
    eclass2(nodes_map('X1'))= 2;
    eclass2(nodes_map('X2'))= 3;
    eclass2(nodes_map('X3'))= 3;
    onodes = dnodes; %[nodes_map('X3')]; 
    bnet = mk_dbn(intra, inter, node_sizes, ...
            'discrete', dnodes, ...
            'eclass1', eclass1, ...
            'eclass2', eclass2, ...
            'observed', onodes);
    bnet.CPD{1} = tabular_CPD(bnet, nodes_map('X1'));
    bnet.CPD{2} = tabular_CPD(bnet, nodes_map('X1')+n);
    bnet.CPD{3} = tabular_CPD(bnet, nodes_map('X2'));
end   


