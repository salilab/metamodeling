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

for time  = 1:length(Gm1)
    [bnet, nodes_map, intra, inter]= make_ode_bnet(Gm1, Im1,time);
end
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

% Plot ODE.G and ODE.I
%figure()
%yyaxis left;
%plot(1:T, sample_seq(nodes_map('ODE.G'),:));
%yyaxis right;
%plot(1:T, sample_seq(nodes_map('ODE.I'),:));
%legend('ODE.G','ODE.I');  

% Plot the distribution of ODE.beta
disp('plot');
fprintf("Sampled time-series of length %d", T);
y = sample_seq(nodes_map('ODE.beta'),:);
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
legend('ODE.beta, prior');
hold off;

% To generate a conditional gaussian model
function [bnet, nodes_map, intra, inter]= make_ode_bnet(Gm, Im, time)
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
    disp(eclass2);
    
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
    %disp(bnet)
end

