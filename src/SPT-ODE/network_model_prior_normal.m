% Make a DBN for the network model with the following variables
%
% Time-dependent variables
%  -> G.Network(t)  ->  G.Network(t+1) ->
%  -> ATP.Network(t)  ->  ATP.Network(t+1) ->
%  -> cAMP.Network(t)  ->  cAMP.Network(t+1) ->
%  -> Ca.Network(t)  ->  Ca.Network(t+1) ->
%  -> S.Network(t)  ->  S.Network(t+1) ->
%
% Reference variables
% G.ref(t), G.ref(t+1)
% S.ref(t), S.ref(t+1)
%
% Observed variables
% G.obs(t), G.obs(t+1)
% S.obs(t), S.obs(t+1)
%
% Time-invariant variables
%
% Parameters
% PFK.Network
%
% To generate a conditional gaussian model

warning('off','MATLAB:singularMatrix');

%Read in the experimental measurements
clear;
% G_Model, PFK_Model,S_Model
[network_dbn_factory]= make_network_dbn_factory(0.1, 63, 30);
[dbn, ~, ~, nodes_map] = create_dbn(network_dbn_factory);
% parameter learning
npers= dbn.nnodes_per_slice;
T = 20; % lengthhs of sequences to explore
%disp(npers);
%%%%%%% try smaller timesteps %%%%%%%

%Sample from the posterior:
%disp(length(bnet.intra))
%disp(length(bnet.inter))
%evidence=cell(n,T);
%evidence{2,1}=1;
%sample_seq=  cell2mat(sample_dbn(bnet, 'length', T,'evidence', evidence));
sample_seq=  cell2mat(sample_dbn(dbn, 'length', T));
nodes_order= cell2mat(values(nodes_map));

% Plot ODE.G and ODE.I
figure()
yyaxis left;
Gvalues = sample_seq(nodes_map('cAMP.Network'),:);
plot(1:T, Gvalues);
yyaxis right;
Ivalues = sample_seq(nodes_map('S.Network'),:);
plot(1:T, Ivalues);
legend('cAMP.Network','S.Network');  
legend('boxoff');

%compute the unconditional marginals of h(20)
dbn_engine = jtree_dbn_inf_engine(dbn);
T = 100;
evidence= cell(npers, T);
[dbn_engine, ll] = enter_evidence(dbn_engine, evidence); % ll is the log marginal likelihood
marg = marginal_nodes(dbn_engine, nodes_map('cAMP.Network'),10);

tol = 1; % tolerance
% evaluates EXPRESSION and, if it is false, displays the error message 'Assertion Failed'
fprintf("Unconditional probability distribution of k(10) is:\n"); 
fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)) % mean +- stddev

%assert(approxeq(marg.mu, 6.1, tol)); % marg.mu equals to 4 +- tol
%assert(approxeq(marg.Sigma, 3.0, tol)); % marg.Sigma equals to 4 +- tol

%Posterior marginal of h(20) given Iexp
T = 400;
evidence{nodes_map('S.ref'),2} = 25.0; 
i=10
disp(i)
marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
                     nodes_map('cAMP.Network')+npers, ...
                     i);
fprintf("Posterior probability distribution of h(10) given Iexp(10) is:\n"); 
fprintf("%f sigma %f +- %f", marg.mu, marg.Sigma, sqrt(marg.Sigma)) % mean +- stddev
xx = linspace(8,12,100);
plot(xx,area*normpdf(xx,marg.mu,marg.Sigma),'k-','LineWidth',2);
legend('network.k, posterior');

% Create a table with the data and variable names
yy=normpdf(xx,marg.mu,marg.Sigma)/sum(normpdf(xx,marg.mu,marg.Sigma))
variable = [xx(:) yy(:)];
size(variable);
dlmwrite('cAMP_posterior.txt',variable);

%compute other posteror marginals with evidence
%T = 400;
%n=nodes_map.Count;
%evidence= cell(npers, T);
%evidence{nodes_map('network.Ipm'),10} = 1.0; 
%evidence{nodes_map('network.k'),10} = 70.0; 

%i=12
%disp(i)
%marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
%                     nodes_map('network.k')+npers, ...
%                     i);
%fprintf("Posterior probability distribution of k is:\n"); 
%fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)) % mean +- stddev