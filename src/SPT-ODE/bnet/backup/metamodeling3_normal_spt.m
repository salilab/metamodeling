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

[dbn, nodes_map] = make_meta_bnet3(5.111, 0, 0.05, 39.6,828, 3, 0.056, 5.111, 34, ...
    0.1, 34, 10, 300, 6, 1.8e-6, ...
    0.1, 63, 34);

[meal_dbn_factory]= make_meal_dbn_factory_eq2(5.111, 0, 0.05, 39.6,828, 3, 0.056, 5.111, 34);
[spt_dbn_factory]= make_spt_dbn_factory(0.1, 34, 10, 300, 6, 1.8e-6);
[network_dbn_factory]= make_network_dbn_factory(0.1, 63, 34);

% Sample
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);

[meal_dbn, ~, ~, meal_nodes_map] = create_dbn(meal_dbn_factory);
%meal_dbn_engine = jtree_dbn_inf_engine(meal_dbn);

[spt_dbn, ~, ~, spt_nodes_map] = create_dbn(spt_dbn_factory);
%spt_dbn_engine = jtree_dbn_inf_engine(spt_dbn);

[network_dbn, ~, ~, network_nodes_map] = create_dbn(network_dbn_factory);
%network_dbn_engine = jtree_dbn_inf_engine(network_dbn);

meal_npers= meal_dbn.nnodes_per_slice;
spt_npers= spt_dbn.nnodes_per_slice;
network_npers= network_dbn.nnodes_per_slice;

% sampling
T = 400;

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

%compute the posterior marginals of h(20)
spt_dbn_engine = jtree_dbn_inf_engine(spt_dbn);
T = 400;
spt_k={};
meta_k={};
spt_evidence= cell(spt_npers, T);
evidence= cell(npers, T);

for measure = 1:400
    disp(measure);
    spt_evidence{spt_nodes_map('S.obs'),measure} = 34; % evidence at time slice 2
    evidence{nodes_map('S.obs'),measure} = 34; % evidence at time slice 2
    %evidence{nodes_map('cAMP.obs'),measure} = 0.008; % evidence at time slice 2
end
%spt_evidence{spt_nodes_map('S.obs'),2} = 57; % evidence at time slice 2
%evidence{nodes_map('S.obs'),2} = 57; % evidence at time slice 2
%evidence{nodes_map('cAMP.obs'),2} = 1.3E-3; % evidence at time slice 2
[spt_engine, ll] = enter_evidence(spt_dbn_engine, spt_evidence);
[engine, ll] = enter_evidence(dbn_engine, evidence);

for measure = 1:400
    marg_spt_k= marginal_nodes(spt_engine,spt_nodes_map('k.SPT'),measure);
    marg_meta_k= marginal_nodes(engine,nodes_map('k.SPT'),measure);
   
    %For tabular nodes, we display marg.T(index of node)
    spt_k(end+1,:) = {marg_spt_k.mu, marg_spt_k.Sigma, sqrt(marg_spt_k.Sigma)};
    meta_k(end+1,:) = {marg_meta_k.mu, marg_meta_k.Sigma, sqrt(marg_meta_k.Sigma)};
    %fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)); % mean +- stddev
end
disp("separate");
disp(spt_k);
disp("separate");
disp(meta_k);


%compute the posterior -2  marginals of h(20)
spt_dbn_engine = jtree_dbn_inf_engine(spt_dbn);
T = 400;
spt_k_2={};
meta_k_2={};
spt_evidence= cell(spt_npers, T);
evidence= cell(npers, T);

for measure = 1:400
    disp(measure);
    spt_evidence{spt_nodes_map('S.obs'),measure} = 34; % evidence at time slice 2
    evidence{nodes_map('S.obs'),measure} = 34; % evidence at time slice 2
    evidence{nodes_map('cAMP.obs'),measure} = 1.3E-3; % evidence at time slice 2
end
%spt_evidence{spt_nodes_map('S.obs'),2} = 57; % evidence at time slice 2
%evidence{nodes_map('S.obs'),2} = 57; % evidence at time slice 2
%evidence{nodes_map('cAMP.obs'),2} = 1.3E-3; % evidence at time slice 2
[spt_engine, ll] = enter_evidence(spt_dbn_engine, spt_evidence);
[engine, ll] = enter_evidence(dbn_engine, evidence);

for measure = 1:400
    marg_spt_k= marginal_nodes(spt_engine,spt_nodes_map('k.SPT'),measure);
    marg_meta_k= marginal_nodes(engine,nodes_map('k.SPT'),measure);
   
    %For tabular nodes, we display marg.T(index of node)
    spt_k_2(end+1,:) = {marg_spt_k.mu, marg_spt_k.Sigma, sqrt(marg_spt_k.Sigma)};
    meta_k_2(end+1,:) = {marg_meta_k.mu, marg_meta_k.Sigma, sqrt(marg_meta_k.Sigma)};
    %fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)); % mean +- stddev
end
disp("separate");
disp(spt_k_2);
disp("separate");
disp(meta_k_2);

% Create a table with the data and variable names
T = table(meta_k, meta_k_2 )
% Write data to text file
writetable(T, 'k_meta3_posterior1-2.txt')
