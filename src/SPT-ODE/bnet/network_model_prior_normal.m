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

clear;

% Initial values of variables
Gin_Network =  5.111/2 % Basal intracellular glucose concentration at the first time slice, mM
PFK_Network = 0.65 % Enzyme activily, no unit
ATP_Network =  3.3 % intracellular ATP concentration at the first time slice, mM
GLP1_Network = 12.5 % intracellular GLP1 concentration at the first time slice, pM 
GLP1R_Network = 1.0 % GLP1R activity at the first time slice, no unit
cAMP_Network =  1.3E-3 % intracellular cAMP concentration at the first time slice, mM
Ca_Network =  1E-4 % intracellular Ca concentration at the first time slice, mM
S_Network = 34 % Basal insulin secretion at the first time slice, normalized for pancreas, pM/min 
I_Network = 25; % Basal plasma insulin level at the first time slice, normalized for pancreas, pM

% G_Model, PFK_Model,S_Model
[network_dbn_factory]= make_network_dbn_factory(Gin_Network, PFK_Network, ATP_Network, GLP1_Network, GLP1R_Network, cAMP_Network, Ca_Network, S_Network, I_Network);
[dbn, ~, ~, nodes_map] = create_dbn(network_dbn_factory);
% parameter learning
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);
T = 420; % lengthhs of sequences to explore, which is 420 min
%disp(npers);
% Sample from the posterior:
%disp(length(bnet.intra))
%disp(length(bnet.inter))
evidence=cell(npers,T);
%evidence{2,1}=1;
%sample_seq=  cell2mat(sample_dbn(bnet, 'length', T,'evidence', evidence));
sample_seq=  cell2mat(sample_dbn(dbn, 'length', T));
nodes_order= cell2mat(values(nodes_map));

% Plot ODE.G and ODE.I
%figure()
yyaxis left;
Gvalues = sample_seq(nodes_map('Gin.Network'),:);
plot(1:T, Gvalues);
yyaxis right;
Ivalues = sample_seq(nodes_map('I.Network'),:);
plot(1:T, Ivalues);
legend('Gin.Network','I.Network');  
legend('boxoff');

%fprintf("prior distribution of G.Meal, mean %f, error, %f \n", mean(Gvalues), std(Gvalues));
%fprintf("prior distribution of I.Meal, mean %f, error, %f ", mean(Ivalues), std(Ivalues));


% testing
%CPD = struct(dbn.CPD{nodes_map('G.ref')});
%disp('here');
%disp(CPD.mean);
%disp(CPD.cov);

%T=10
evidence= cell(npers, T);
%evidence{nodes_map('G.obs'),1} = Gexp(1)

%[engine, ll] = enter_evidence(dbn_engine, evidence);
%disp(ll);
%margG= marginal_nodes(engine,nodes_map('G.Meal'),1);

%disp(Gexp(1))
%disp(margG.mu)
%disp(sqrt(margG.Sigma))
% Parameter estimation from submodels and correct overall network topology
Gin={}
GLP1R={}
cAMP={}
Ca={}
S={}
I={}
for measure = 1:420
    evidence{nodes_map('Gin.ref'),measure} = Gin_Network; % evidence at time slice 2
end
%evidence{nodes_map('I.obs'),measure} = Iexp(measure); 
[engine, ll] = enter_evidence(dbn_engine, evidence);

for measure = 1:420
    margGin= marginal_nodes(engine,nodes_map('Gin.Network'),measure);
    margGLP1R= marginal_nodes(engine,nodes_map('GLP1R.Network'),measure);
    margcAMP= marginal_nodes(engine,nodes_map('cAMP.Network'),measure);
    margCa= marginal_nodes(engine,nodes_map('Ca.Network'),measure);
    margS= marginal_nodes(engine,nodes_map('S.Network'),measure);
    margI= marginal_nodes(engine,nodes_map('I.Network'),measure);

    %For tabular nodes, we display marg.T(index of node)
    Gin(end+1,:) = {margGin.mu, margGin.Sigma, sqrt(margGin.Sigma)};
    GLP1R(end+1,:) = {margGLP1R.mu, margGLP1R.Sigma, sqrt(margGLP1R.Sigma)};
    cAMP(end+1,:) = {margcAMP.mu, margcAMP.Sigma, sqrt(margcAMP.Sigma)};
    Ca(end+1,:) = {margCa.mu, margCa.Sigma, sqrt(margCa.Sigma)};
    I(end+1,:) = {margI.mu, margI.Sigma, sqrt(margI.Sigma)};
    S(end+1,:) = {margS.mu, margS.Sigma, sqrt(margS.Sigma)};
    %fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)); % mean +- stddev
end
%disp("separate");
%disp(Gintake);
%disp("separate");
%disp(G);
%disp("separate");
%disp(DG);
%disp("separate");
%disp(I);
disp("separate");
disp(S);

% Create a table with the data and variable names
T = table(Gin, GLP1R, cAMP, Ca, S, I);
% Write data to text file
writetable(T, 'Network_prior.txt');
