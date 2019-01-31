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

% check the time step of the equation, and see the G and I.

% Read in the experimental measurements
odem1 = importdata('ode_exp1_avr.dat');
Gm1 = odem1(:,2); % Gexp in measurement number 1, vector along time
Im1 = odem1(:,3); % Iexp in measurement number 1, vector along time
disp(Gm1(1));

[bnet, nodes_map, intra, inter]= make_ode_bnet(Gm1, Im1);
n = 9;

% parameter learning
npers= bnet.nnodes_per_slice;
T = 200; % lengthhs of sequences to explore
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
y = sample_seq(nodes_map('beta'),:);
nbins = 10;

[hts,ctrs] = hist(y, nbins);
h = bar(ctrs,hts,'hist');
set(h,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0 0 0]);
area = sum(hts) * (ctrs(2)-ctrs(1));
xx = linspace(-0.1,0.3);
hold on; 
plot(xx,area*normpdf(xx,mean(y),std(y)),'k-','LineWidth',2);
fprintf("Normal probability density function of beta");
disp(mean(y));
disp(std(y));
%f = ksdensity(y,xx);
%plot(xx,area*f,'g-')
legend('beta, prior');
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
function [bnet, nodes_map, intra, inter]= make_ode_bnet(Gm, Im)
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
    
    eclass1 = zeros(n,1); 
    eclass1(nodes_map('alpha'))= 1;
    eclass1(nodes_map('beta'))= 2;
    eclass1(nodes_map('h'))= 3;
    eclass1(nodes_map('G'))= 4;
    eclass1(nodes_map('Gref'))= 5;
    eclass1(nodes_map('Gexp'))= 6;
    eclass1(nodes_map('I'))= 7;
    eclass1(nodes_map('Iref'))= 8;
    eclass1(nodes_map('Iexp'))= 9;
    
    eclass2 = zeros(n,1); 
    eclass2(nodes_map('alpha'))= 1;
    eclass2(nodes_map('beta'))= 2;
    eclass2(nodes_map('h'))= 3;
    eclass2(nodes_map('G'))= 10;
    eclass2(nodes_map('Gref'))= 5;
    eclass2(nodes_map('Gexp'))= 6;
    eclass2(nodes_map('I'))= 11;
    eclass2(nodes_map('Iref'))= 8;
    eclass2(nodes_map('Iexp'))= 9;
    
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
    bnet.CPD{5} = gaussian_CPD(bnet, nodes_map('Gref'),   'mean', Gm(1), 'cov', 0.1,'weights', 1);
    bnet.CPD{6} = gaussian_CPD(bnet, nodes_map('Gexp'), 'mean', Gm(1), 'cov', 0.1, 'weights', 1);
    bnet.CPD{7} = gaussian_CPD(bnet, nodes_map('I'),   'mean', Im(1), 'cov', 5);
    bnet.CPD{8} = gaussian_CPD(bnet, nodes_map('Iref'),   'mean', Im(1), 'cov', 5, 'weights', 1);
    bnet.CPD{9} = gaussian_CPD(bnet, nodes_map('Iexp'), 'mean', Im(1), 'cov', 5, 'weights', 1);
    
    % eclass2
    bnet.CPD{10} = gaussian_CPD(bnet, nodes_map('G')+n, 'mean', Gm(n), 'cov', 0.1,  'weights', 0.5);
    
    weights_I=zeros(5,1);
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
    bnet.CPD{11} = gaussian_CPD(bnet, nodes_map('I')+n,   'mean', Im(n), 'cov', 5, 'weights', weights_I);
    %bnet.CPD{16} = gaussian_CPD(bnet, nodes_map('I')+n,   'mean', 72, 'cov', 5);
    
    disp('here');
    disp(max(eclass2));
    disp('here');
end

