% Make a DBN for the spatiotemporal BD model with the following variables
%
% Time-dependent variables
%  -> G.SPT(t)  ->  G.SPT(t+1) ->
%  -> S.SPT(t)  ->  S.SPT(t+1) ->
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
% lambda.SPT k.SPT Npatch.SPT Nisg.SPT Disg.SPT Ninsulin.SPT Rpbc.SPT
%
% All variables display gaussian distributions.
warning('off','MATLAB:singularMatrix');

%Read in the experimental measurements
clear;

% Initial values of variables
Gin_SPT =  5.111/2 % Basal intracellular glucose concentration at the first time slice, mM
S_SPT =  34 % Basal insulin secretion at the first time slice, normalized for pancreas, pM/min 
k_SPT =  10 % Force coefficient, m/s
Nisg_SPT = 300 %  Number of ISGs
Npatch_SPT = 6 % Number of activation patches per ISG
Ninsulin_SPT = 1.8e-6 % Amount of insulin molecules in one ISG, pmol
I_SPT = 25; % Basal plasma insulin level at the first time slice, normalized for pancreas, pM

% G_Model, S_Model, k_Model, Nisg_Model, Npatch_Model, Ninsulin_Model
[spt_dbn_factory]= make_spt_dbn_factory(Gin_SPT, S_SPT, k_SPT, Nisg_SPT, Npatch_SPT, Ninsulin_SPT, I_SPT);
[dbn, ~, ~, nodes_map] = create_dbn(spt_dbn_factory);
% parameter learning
npers= dbn.nnodes_per_slice;
T = 420; % lengthhs of sequences to explore
%disp(npers);

%Sample from the prior:
sample_seq=  cell2mat(sample_dbn(dbn, 'length', T));
nodes_order= cell2mat(values(nodes_map));

% Plot the distribution of SPT.beta
%disp('plot');
%fprintf("Sampled time-series of length %d", T);
y = sample_seq(nodes_map('k.SPT'),:);
nbins = 10;
[hts,ctrs] = hist(y, nbins);
h = bar(ctrs,hts,'hist');
set(h,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0 0 0]);
area = sum(hts) * (ctrs(2)-ctrs(1));
xx = linspace(8,12,100);
hold on; 
plot(xx,area*normpdf(xx,mean(y),std(y)),'k-','LineWidth',2);
fprintf("Normal probability density function of SPT.k");
disp(mean(y));
disp(std(y)*std(y));
disp(std(y));
%f = ksdensity(y,xx);
%plot(xx,area*f,'g-')
%legend('SPT.k, prior');
%hold off;

% test different measurement
Experiment1 = importdata('/Users/lipingsun/SynologyDrive/research/Meta-modeling/bnt/Meta-modelingS/dataset/meal/140points/data_normal_G-S.dat');
Gexp = Experiment1(:,1); % Gexp in measurement number 1, vector along time
Sexp = Experiment1(:,2); % Gexp in measurement number 1, vector along time
Experiment2 = importdata('/Users/lipingsun/SynologyDrive/research/Meta-modeling/bnt/Meta-modelingS/dataset/meal/140points/meal_exp1_normal2_DG.dat');
DGexp = Experiment2(:,2); % Gexp in measurement number 1, vector along time

%compute the prior marginals of h(20)
spt_dbn_engine = jtree_dbn_inf_engine(dbn);
T = 420;
spt_k_prior={};
spt_Npatch={}
spt_S={};
spt_evidence= cell(npers, T);

for measure = 1:420
    disp(measure);
    %spt_evidence{nodes_map('Npatch.SPT'),measure} = 6; % evidence at time slice 2
    %spt_evidence{spt_nodes_map('Gin.SPT.obs'),measure} = Gexp(measure)/2; % evidence at time slice 2
    %evidence{nodes_map('DG.obs'),measure} = DGexp(measure); % evidence at time slice 2
end
%spt_evidence{nodes_map('S.obs'),2} = 55; % evidence at time slice 2
[spt_engine, ll] = enter_evidence(spt_dbn_engine, spt_evidence);

for measure = 1:420
    marg_spt_k= marginal_nodes(spt_engine,nodes_map('k.SPT'),measure);
    marg_spt_Npatch= marginal_nodes(spt_engine,nodes_map('Npatch.SPT'),measure);
    marg_spt_S= marginal_nodes(spt_engine,nodes_map('S.SPT'),measure);
    %For tabular nodes, we display marg.T(index of node)
    spt_k_prior(end+1,:) = {marg_spt_k.mu, marg_spt_k.Sigma, sqrt(marg_spt_k.Sigma)};
    spt_Npatch(end+1,:) = {marg_spt_Npatch.mu, marg_spt_Npatch.Sigma, sqrt(marg_spt_Npatch.Sigma)};
    spt_S(end+1,:) = {marg_spt_S.mu, marg_spt_S.Sigma, sqrt(marg_spt_S.Sigma)};
    %fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)); % mean +- stddev
end
disp("separate");
disp(spt_k_prior);
disp("separate");
disp(spt_Npatch);
disp("separate");
disp(spt_S);

%compute the posterior marginals of h(20)
spt_dbn_engine = jtree_dbn_inf_engine(dbn);
T = 420;
spt_k={};
spt_Npatch={}
spt_S={};
spt_evidence= cell(npers, T);

for measure = 1:420
    disp(measure);
    spt_evidence{nodes_map('Scell.ref'),measure} = 34; % evidence at time slice 2
end
%spt_evidence{nodes_map('S.obs'),2} = 55; % evidence at time slice 2
[spt_engine, ll] = enter_evidence(spt_dbn_engine, spt_evidence);

for measure = 1:420
    marg_spt_k= marginal_nodes(spt_engine,nodes_map('k.SPT'),measure);
    marg_spt_Npatch= marginal_nodes(spt_engine,nodes_map('Npatch.SPT'),measure);
    marg_spt_S= marginal_nodes(spt_engine,nodes_map('S.SPT'),measure);
    %For tabular nodes, we display marg.T(index of node)
    spt_k(end+1,:) = {marg_spt_k.mu, marg_spt_k.Sigma, sqrt(marg_spt_k.Sigma)};
    spt_Npatch(end+1,:) = {marg_spt_Npatch.mu, marg_spt_Npatch.Sigma, sqrt(marg_spt_Npatch.Sigma)};
    spt_S(end+1,:) = {marg_spt_S.mu, marg_spt_S.Sigma, sqrt(marg_spt_S.Sigma)};
    %fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)); % mean +- stddev
end
disp("separate");
disp(spt_k);
disp("separate");
%disp(spt_Npatch);
disp("separate");
%disp(spt_S);

% Create a table with the data and variable names
T = table(spt_k_prior, spt_k )
% Write data to text file
writetable(T, 'k_spt_prior-posterior.txt')

