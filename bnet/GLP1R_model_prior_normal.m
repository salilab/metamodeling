% Make a DBN for the GLP1R model with the following variables
%
% Time-dependent variables
%  -> GLP1.GLP1R(t)  ->  GLP1.GLP1R(t+1) ->
%  -> GLP1R.GLP1R(t)  ->  GLP1R.GLP1R(t+1) ->
%  -> Galpha.GLP1R(t)  ->  Galpha.GLP1R(t+1) ->
%  -> cAMP.GLP1R(t)  ->  cAMP.GLP1R(t+1) ->
%
% Reference variables
% cAMP.ref(t), cAMP.ref(t+1) 
%
% Observed variables
% cAMP.obs(t), cAMP.obs(t+1)
%
% Time-invariant variables
% 
%
% Parameters
%
% To generate a conditional gaussian model

warning('off','MATLAB:singularMatrix');

%Read in the experimental measurements
clear;

% Initial values of variables
GLP1a_GLP1R =  0.0 % Type of the GLP1 analogues, no unit - analogue 1-10 with increasing affinity to GLP1R
cons_GLP1R = 0.0 % Compound concentration of the GLP1 analogues , nM - 1-10
GLP1R_GLP1R =  0.0 % GLP1R activity at the first time slice, no unit

% GLP1_Model,GLP1R_Model, Galpha_Model, cAMP_Model
[glp1r_dbn_factory]= make_glp1r_dbn_factory(GLP1a_GLP1R, cons_GLP1R, GLP1R_GLP1R);
[dbn, ~, ~, nodes_map] = create_dbn(glp1r_dbn_factory);
% parameter learning
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);
T = 420; % lengthhs of sequences to explore
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
Gvalues = sample_seq(nodes_map('GLP1a.GLP1R'),:);
plot(1:T, Gvalues);
yyaxis right;
Ivalues = sample_seq(nodes_map('GLP1R.GLP1R'),:);
plot(1:T, Ivalues);
legend('GLP1a.GLP1R','GLP1R.GLP1R');  
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
GLP1a={}
cons={}
GLP1R={}

for measure = 1:420
    evidence{nodes_map('GLP1a.obs'),measure} = 10; % evidence at time slice 2
    evidence{nodes_map('cons.obs'),measure} = 10; % evidence at time slice 2
end
%evidence{nodes_map('I.obs'),measure} = Iexp(measure); 
[engine, ll] = enter_evidence(dbn_engine, evidence);

for measure = 1:420
    margGLP1a= marginal_nodes(engine,nodes_map('GLP1a.GLP1R'),measure);
    margcons= marginal_nodes(engine,nodes_map('cons.GLP1R'),measure);
    margGLP1R= marginal_nodes(engine,nodes_map('GLP1R.GLP1R'),measure);

    %For tabular nodes, we display marg.T(index of node)
    GLP1a(end+1,:) = {margGLP1a.mu, margGLP1a.Sigma, sqrt(margGLP1a.Sigma)};
    cons(end+1,:) = {margcons.mu, margcons.Sigma, sqrt(margcons.Sigma)};
    GLP1R(end+1,:) = {margGLP1R.mu, margGLP1R.Sigma, sqrt(margGLP1R.Sigma)};

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
disp(GLP1R);

% Create a table with the data and variable names
T = table(GLP1a,cons, GLP1R);
% Write data to text file
writetable(T, 'GLP1R_prior.txt');
