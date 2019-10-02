% Make a DBN for the meal model with the following variables
%
% Time-dependent variables
% 
% Reference variables
%
% Observed variables
%
% Time-invariant variables
%
% Parameters
%
% The ODE equation is only when G > h
% TODO: check the time step of the equation, and see the G and I.
warning('off','MATLAB:singularMatrix');

%Read in the experimental measurements
clear;
%Gex_Meal, Y_Meal, alpha_Meal,beta_Meal, K_Meal, dt_Meal, DGex_Meal, Gb_Meal, Sb_Meal
%Gin_Network, PFK_Network,S_Network
[dbn, nodes_map] = make_meta_bnet4(5.111, 0, 0.05, 39.6,828, 3, 0.056, 5.111, 34, ...
    5.111/2, 34, 10, 300, 6, 1.8e-6, ...
    5.111/2, 3, 34, ...
    0.0, 0.0);

[meal_dbn_factory]= make_meal_dbn_factory_eq2(5.111, 0, 0.05, 39.6,828, 3, 0.056, 5.111, 34);
[spt_dbn_factory]= make_spt_dbn_factory(0.1, 34, 10, 300, 6, 1.8e-6);
[network_dbn_factory]= make_network_dbn_factory(0.1, 3, 34);
[glp1r_dbn_factory]= make_glp1r_dbn_factory(0.0,0.0);

% Sample
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);

[meal_dbn, ~, ~, meal_nodes_map] = create_dbn(meal_dbn_factory);
%meal_dbn_engine = jtree_dbn_inf_engine(meal_dbn);

[spt_dbn, ~, ~, spt_nodes_map] = create_dbn(spt_dbn_factory);
%spt_dbn_engine = jtree_dbn_inf_engine(spt_dbn);

[network_dbn, ~, ~, network_nodes_map] = create_dbn(network_dbn_factory);
%network_dbn_engine = jtree_dbn_inf_engine(network_dbn);

[glp1r_dbn, ~, ~, glp1r_nodes_map] = create_dbn(glp1r_dbn_factory);
%glp1r_dbn_engine = jtree_dbn_inf_engine(glp1r_dbn);

meal_npers= meal_dbn.nnodes_per_slice;
spt_npers= spt_dbn.nnodes_per_slice;
network_npers= network_dbn.nnodes_per_slice;
glp1r_npers= glp1r_dbn.nnodes_per_slice;

% sampling
T = 140;

sample_seq=  cell2mat(sample_dbn(dbn, 'length', T));
nodes_order= cell2mat(values(nodes_map));

% Plot ODE.G and ODE.I
%figure()
%yyaxis left;
%Gvalues = sample_seq(nodes_map('cAMP.Network'),:);
%plot(1:T, Gvalues);
%yyaxis right;
%Ivalues = sample_seq(nodes_map('S.ref'),:);
%plot(1:T, Ivalues);
%legend('cAMP.Network','S.Network');  
%legend('boxoff');

% test different measurement
Experiment1 = importdata('/Users/lipingsun/SynologyDrive/research/Meta-modeling/bnt/Meta-modelingS/dataset/meal/140points/data_normal_G-S.dat');
Gexp = Experiment1(:,1); % Gexp in measurement number 1, vector along time
Sexp = Experiment1(:,2); % Gexp in measurement number 1, vector along time
Experiment2 = importdata('/Users/lipingsun/SynologyDrive/research/Meta-modeling/bnt/Meta-modelingS/dataset/meal/140points/meal_exp1_normal2_DG.dat');
DGexp = Experiment2(:,2); % Gexp in measurement number 1, vector along time

%compute the unconditional marginals of h(20)
meal_dbn_engine = jtree_dbn_inf_engine(meal_dbn);
spt_dbn_engine = jtree_dbn_inf_engine(spt_dbn);
network_dbn_engine = jtree_dbn_inf_engine(network_dbn);
glp1r_dbn_engine = jtree_dbn_inf_engine(glp1r_dbn);

T = 140;
spt_k={};
meta_k={};
spt_Nisg={};
meta_Nisg={};
spt_Ninsulin={};
meta_Ninsulin={};
spt_Npatch={};
meta_Npatch={};
meal_evidence= cell(meal_npers, T);
spt_evidence= cell(spt_npers, T);
network_evidence= cell(network_npers, T);
glp1r_evidence= cell(glp1r_npers, T);
evidence= cell(npers, T);

for measure = 1:140
    new = measure
    disp(measure);
    evidence{nodes_map('Gex.obs'),measure} = Gexp(new); % evidence at time slice 2
    evidence{nodes_map('S.obs'),measure} = Sexp(new); % evidence at time slice 2
    evidence{nodes_map('DGex.obs'),measure} = DGexp(new); % evidence at time slice 2
    evidence{nodes_map('Gin.SPT.obs'),measure} = Gexp(new)/2; % evidence at time slice 2
    evidence{nodes_map('Gin.Network.obs'),measure} = Gexp(new)/2; % evidence at time slice 2
    evidence{nodes_map('GLP1.obs'),measure} = 0.0; % evidence at time slice 2
    spt_evidence{spt_nodes_map('Gin.SPT.obs'),measure} = Gexp(new)/2; % evidence at time slice 2
    %spt_evidence{spt_nodes_map('S.obs'),measure} = Sexp(new); % evidence at time slice 2
end

[spt_engine, ll] = enter_evidence(spt_dbn_engine, spt_evidence);
[engine, ll] = enter_evidence(dbn_engine, evidence);

for measure = 1:140
    marg_spt_k= marginal_nodes(spt_engine,spt_nodes_map('k.SPT'),measure);
    marg_meta_k= marginal_nodes(engine,nodes_map('k.SPT'),measure);
    spt_k(end+1,:) = {marg_spt_k.mu, marg_spt_k.Sigma, sqrt(marg_spt_k.Sigma)};
    meta_k(end+1,:) = {marg_meta_k.mu, marg_meta_k.Sigma, sqrt(marg_meta_k.Sigma)};
    marg_spt_Nisg= marginal_nodes(spt_engine,spt_nodes_map('Nisg.SPT'),measure);
    marg_meta_Nisg= marginal_nodes(engine,nodes_map('Nisg.SPT'),measure);
    spt_Nisg(end+1,:) = {marg_spt_Nisg.mu, marg_spt_Nisg.Sigma, sqrt(marg_spt_Nisg.Sigma)};
    meta_Nisg(end+1,:) = {marg_meta_Nisg.mu, marg_meta_Nisg.Sigma, sqrt(marg_meta_Nisg.Sigma)};
    marg_spt_Ninsulin= marginal_nodes(spt_engine,spt_nodes_map('Ninsulin.SPT'),measure);
    marg_meta_Ninsulin= marginal_nodes(engine,nodes_map('Ninsulin.SPT'),measure);
    spt_Ninsulin(end+1,:) = {marg_spt_Ninsulin.mu, marg_spt_Ninsulin.Sigma, sqrt(marg_spt_Ninsulin.Sigma)};
    meta_Ninsulin(end+1,:) = {marg_meta_Ninsulin.mu, marg_meta_Ninsulin.Sigma, sqrt(marg_meta_Ninsulin.Sigma)};
end
disp("separate");
disp(spt_k);
disp("separate");
disp(meta_k);
disp("separate");
disp(spt_Nisg);
disp("separate");
disp(meta_Nisg);
disp("separate");
disp(spt_Ninsulin);
disp("separate");
disp(meta_Ninsulin);
% GLP-1 effect
meal_evidence= cell(meal_npers, T);
spt_evidence= cell(spt_npers, T);
network_evidence= cell(network_npers, T);
glp1r_evidence= cell(glp1r_npers, T);
evidence= cell(npers, T);

glp1_k={};
prior_k={};
glp1_Nisg={};
prior_Nisg={};
glp1_Ninsulin={};
prior_Ninsulin={};
for measure = 1:140
    new = measure
    disp(measure);
    evidence{nodes_map('Gex.obs'),measure} = Gexp(new); % evidence at time slice 2
    evidence{nodes_map('S.obs'),measure} = Sexp(new); % evidence at time slice 2
    evidence{nodes_map('DGex.obs'),measure} = DGexp(new); % evidence at time slice 2
    evidence{nodes_map('Gin.SPT.obs'),measure} = Gexp(new)/2; % evidence at time slice 2
    evidence{nodes_map('Gin.Network.obs'),measure} = Gexp(new)/2; % evidence at time slice 2
    evidence{nodes_map('GLP1.obs'),measure} = 5E-5
    ; % evidence at time slice 2
end

[spt_engine, ll] = enter_evidence(spt_dbn_engine, spt_evidence);
[engine, ll] = enter_evidence(dbn_engine, evidence);

for measure = 1:140
    marg_prior_k= marginal_nodes(spt_engine,spt_nodes_map('k.SPT'),measure);
    marg_glp1_k= marginal_nodes(engine,nodes_map('k.SPT'),measure);
    prior_k(end+1,:) = {marg_prior_k.mu, marg_prior_k.Sigma, sqrt(marg_prior_k.Sigma)};
    glp1_k(end+1,:) = {marg_glp1_k.mu, marg_glp1_k.Sigma, sqrt(marg_glp1_k.Sigma)};
    marg_prior_Nisg= marginal_nodes(spt_engine,spt_nodes_map('Nisg.SPT'),measure);
    marg_glp1_Nisg= marginal_nodes(engine,nodes_map('Nisg.SPT'),measure);
    prior_Nisg(end+1,:) = {marg_prior_Nisg.mu, marg_prior_Nisg.Sigma, sqrt(marg_prior_Nisg.Sigma)};
    glp1_Nisg(end+1,:) = {marg_glp1_Nisg.mu, marg_glp1_Nisg.Sigma, sqrt(marg_glp1_Nisg.Sigma)};
    marg_prior_Ninsulin= marginal_nodes(spt_engine,spt_nodes_map('Ninsulin.SPT'),measure);
    marg_glp1_Ninsulin= marginal_nodes(engine,nodes_map('Ninsulin.SPT'),measure);
    prior_Ninsulin(end+1,:) = {marg_prior_Ninsulin.mu, marg_prior_Ninsulin.Sigma, sqrt(marg_prior_Ninsulin.Sigma)};
    glp1_Ninsulin(end+1,:) = {marg_glp1_Ninsulin.mu, marg_glp1_Ninsulin.Sigma, sqrt(marg_glp1_Ninsulin.Sigma)};
end
disp("separate");
disp(glp1_k);
disp("separate");
disp(prior_k);
disp("separate");
disp(glp1_Nisg);
disp("separate");
disp(prior_Nisg);
disp("separate");
disp(glp1_Ninsulin);
disp("separate");
disp(prior_Ninsulin);

% Create a table with the data and variable names
T = table(spt_k, meta_k, glp1_k, prior_k);
% Write data to text file
writetable(T, 'meta_k_normal.txt');

% Create a table with the data and variable names
T = table(spt_Nisg, meta_Nisg, glp1_Nisg, prior_Nisg);
% Write data to text file
writetable(T, 'meta_Nisg_normal.txt');

% Create a table with the data and variable names
T = table(spt_Ninsulin, meta_Ninsulin, glp1_Ninsulin, prior_Ninsulin);
% Write data to text file
writetable(T, 'meta_Ninsulin_normal.txt');
