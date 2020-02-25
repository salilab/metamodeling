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
Gin_KEGG = 0 % Extracellular glucose concentration, mM 
Ex4_KEGG = 0 % GLP-1 treatmeent, nM
NES_KEGG = 1 % Normalized enrichment score, no unit
Pathway_KEGG = 1 % Pathway number, no unit
ATP_KEGG = 0 % intracellular ATP concentration, mM

% G_Model, PFK_Model,S_Model
[network_dbn_factory]= make_KEGG_dbn_factory(Gin_KEGG, Ex4_KEGG, NES_KEGG, Pathway_KEGG, ATP_KEGG);
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
Gvalues = sample_seq(nodes_map('Gin.KEGG'),:);
plot(1:T, Gvalues);
yyaxis right;
Ivalues = sample_seq(nodes_map('ATP.KEGG'),:);
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
%evidence{nodes_map('G.obs'),1} = Ginp(1)

%[engine, ll] = enter_evidence(dbn_engine, evidence);
%disp(ll);
%margG= marginal_nodes(engine,nodes_map('G.Meal'),1);

%disp(Ginp(1))
%disp(margG.mu)
%disp(sqrt(margG.Sigma))
% Parameter estimation from submodels and correct overall network topology

for measure = 1:420
    evidence{nodes_map('Gin.KEGG'),measure} = 0; % evidence at time slice 2
    evidence{nodes_map('Ex4.ref'),measure} = 0; % evidence at time slice 2
end
%evidence{nodes_map('I.obs'),measure} = Iexp(measure); 
[engine, ll] = enter_evidence(dbn_engine, evidence);

Gin={};
Ex4={};
NES={}; 
Pathway={};
ATP={};

for measure = 1:420
    margGin= marginal_nodes(engine,nodes_map('Gin.KEGG'),measure);
    marGin4= marginal_nodes(engine,nodes_map('Ex4.KEGG'),measure);
    margNES= marginal_nodes(engine,nodes_map('NES.KEGG'),measure);
    margPathway= marginal_nodes(engine,nodes_map('Pathway.KEGG'),measure);
    margATP= marginal_nodes(engine,nodes_map('ATP.KEGG'),measure);

    %For tabular nodes, we display marg.T(index of node)
    Gin(end+1,:) = {margGin.mu, margGin.Sigma, sqrt(margGin.Sigma)};
    Ex4(end+1,:) = {marGin4.mu, marGin4.Sigma, sqrt(marGin4.Sigma)};
    NES(end+1,:) = {margNES.mu, margNES.Sigma, sqrt(margNES.Sigma)};
    Pathway(end+1,:) = {margPathway.mu, margPathway.Sigma, sqrt(margPathway.Sigma)};
    ATP(end+1,:) = {margATP.mu, margATP.Sigma, sqrt(margATP.Sigma)};
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
disp(ATP);

% Create a table with the data and variable names
T = table(Gin, Ex4, NES, Pathway, ATP);
% Write data to text file
writetable(T, 'KEGG_prior.txt');
