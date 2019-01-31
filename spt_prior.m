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

% Read in the experimental measurements
sptm1 = importdata('spt_obs1_avr.dat');
Go1 = sptm1(:,2); % Gexp in measurement number 1, vector along time
Io1 = sptm1(:,3); % Iexp in measurement number 1, vector along time
disp(Go1(1));

[bnet, nodes_map, intra, inter]= make_spt_bnet(Go1, Io1);
n = 11;

% parameter learning
npers= bnet.nnodes_per_slice;
T = 100; % lengthhs of sequences to explore
disp(npers);
% Sample from the posterior:
%disp(length(bnet.intra))
%disp(length(bnet.inter))
%evidence=cell(n,T);
%evidence{2,1}=1;
%sample_seq=  cell2mat(sample_dbn(bnet, 'length', T,'evidence', evidence));
sample_seq=  cell2mat(sample_dbn(bnet, 'length', T));
nodes_order= cell2mat(values(nodes_map));

% Plot
disp('plot');
fprintf("Sampled time-series of length %d", T);
y = sample_seq(nodes_map('Npatch'),:);
nbins = 10;

[hts,ctrs] = hist(y,nbins);
h = bar(ctrs,hts,'hist');
set(h,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0 0 0]);
area = sum(hts) * (ctrs(2)-ctrs(1));
xx = linspace(3,9);
hold on; 
plot(xx,area*normpdf(xx,mean(y),std(y)),'k-','LineWidth',2);
fprintf("Normal probability density function of Npatch");
disp(mean(y));
disp(std(y));
%f = ksdensity(y,xx);
%plot(xx,area*f,'g-')
legend('Npatch, prior');
hold off;

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
    %order= topological_sort(intra);
    %intra= intra(order, order);
    %inter= inter(order, order);
    %nodes_map= containers.Map(node_names(order), 1:n);
end

% To generate a conditional gaussian model
function [bnet, nodes_map, intra, inter]= make_spt_bnet(Gm, Im)
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
    eclass2(nodes_map('Gref'))= 7;
    eclass2(nodes_map('Gobs'))= 8;
    eclass2(nodes_map('I'))= 13;
    eclass2(nodes_map('Iref'))= 10;
    eclass2(nodes_map('Iobs'))= 11;
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
    bnet.CPD{1} = gaussian_CPD(bnet, nodes_map('lamda'), 'mean', 0.1, 'cov', 0.01);
    bnet.CPD{2} = gaussian_CPD(bnet, nodes_map('k'), 'mean', 10, 'cov', 1);
    bnet.CPD{3} = gaussian_CPD(bnet, nodes_map('Npatch'), 'mean', 6, 'cov', 0.5);
    bnet.CPD{4} = gaussian_CPD(bnet, nodes_map('Nisg'), 'mean', 300, 'cov', 20);
    bnet.CPD{5} = gaussian_CPD(bnet, nodes_map('Rpbc'), 'mean', 4, 'cov', 0.1);
    bnet.CPD{6} = gaussian_CPD(bnet, nodes_map('G'),   'mean', Gm(1), 'cov', 2);
    bnet.CPD{7} = gaussian_CPD(bnet, nodes_map('Gref'),   'mean', Gm(1), 'cov', 2,'weights', 1);
    bnet.CPD{8} = gaussian_CPD(bnet, nodes_map('Gobs'), 'mean', Gm(1), 'cov', 2, 'weights', 1);
    bnet.CPD{9} = gaussian_CPD(bnet, nodes_map('I'),   'mean', Im(1), 'cov', 5);
    bnet.CPD{10} = gaussian_CPD(bnet, nodes_map('Iref'),   'mean', Im(1), 'cov', 5, 'weights', 1);
    bnet.CPD{11} = gaussian_CPD(bnet, nodes_map('Iobs'), 'mean', Im(1), 'cov', 5, 'weights', 1);
  
    %eclass2
    weights_G= [0.5 0.5];
    bnet.CPD{12} = gaussian_CPD(bnet, nodes_map('G')+n, 'mean', Gm(n), 'cov', 2, 'weights', weights_G);
    
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
    bnet.CPD{13} = gaussian_CPD(bnet, nodes_map('I')+n, 'mean', Im(n), 'cov', 5, 'weights', weights_I);
    %bnet.CPD{13} = gaussian_CPD(bnet, nodes_map('I')+n, 'mean', 70, 'cov', 5);
end    

