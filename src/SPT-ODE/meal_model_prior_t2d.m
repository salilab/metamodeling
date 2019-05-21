% Make a DBN for the meal model with the following variables
%
% Time-dependent variables
%  -> G.Meal(t)  ->  G.Meal(t+1) ->
%  -> I.Meal(t)  ->  I.Meal(t+1) ->
%
% Reference variables
% G.ref(t), G.ref(t+1)
% I.ref(t), I.ref(t+1)
%
% Observed variables
% G.obs(t), G.obs(t+1)
% I.obs(t), I.obs(t+1)
%
% Time-invariant variables
% Gb.Meal
%
% Parameters
% ALPHA BETA
%
% The ODE equation is only when G > h
% TODO: check the time step of the equation, and see the G and I.
warning('off','MATLAB:singularMatrix');

%Read in the experimental measurements
clear;
% G_Model, I_Model, alpha_Model,beta_Model, Gb_Model
%[meal_dbn_factory]= make_meal_dbn_factory_eq1(5.1, 34, 0.05, 39.6, 2, 5.1);
%G_Model, I_Model, alpha_Model,beta_Model, K_Model, dt_Model, DG_Model, Gb_Model, Sb_Model
[meal_dbn_factory]= make_meal_dbn_factory_eq2(9.167, 0, 0.013, 22.5,445.5, 3, 0.03, 9.167, 102.5);
[dbn, ~, ~, nodes_map] = create_dbn(meal_dbn_factory);
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);
T = 140; % lengthhs of sequences to explore
%disp(npers);
% Sample from the posterior:
%disp(length(bnet.intra))
%disp(length(bnet.inter))
%evidence=cell(n,T);
%evidence{2,1}=1;
%sample_seq=  cell2mat(sample_dbn(bnet, 'length', T,'evidence', evidence));
sample_seq=  cell2mat(sample_dbn(dbn, 'length', T));
nodes_order= cell2mat(values(nodes_map));

% Plot ODE.G and ODE.I
%figure()
%yyaxis left;
%Gvalues = sample_seq(nodes_map('G.Meal'),:);
%plot(1:T, Gvalues);
%yyaxis right;
%Ivalues = sample_seq(nodes_map('I.Meal'),:);
%plot(1:T, Ivalues);
%legend('G.Meal','I.Meal');  
%legend('boxoff');

%fprintf("prior distribution of G.Meal, mean %f, error, %f \n", mean(Gvalues), std(Gvalues));
%fprintf("prior distribution of I.Meal, mean %f, error, %f ", mean(Ivalues), std(Ivalues));

% test different measurement
Experiment1 = importdata('/Users/lipingsun/research/Meta-modeling/bnt/Meta-modeling2/dataset/meal/140points/meal_exp1_t2d2.dat');
Gexp = Experiment1(:,2); % Gexp in measurement number 1, vector along time
Experiment2 = importdata('/Users/lipingsun/research/Meta-modeling/bnt/Meta-modeling2/dataset/meal/140points/meal_exp1_t2d2_DG.dat');
DGexp = Experiment2(:,2); % Gexp in measurement number 1, vector along time

% Parameter estimation from submodels and correct overall network topology
G={}
Gb={}
DG={}
S={}

% testing
%CPD = struct(dbn.CPD{nodes_map('G.ref')});
%disp('here');
%disp(CPD.mean);
%disp(CPD.cov);

%T=10
%evidence= cell(npers, T);
%evidence{nodes_map('G.obs'),1} = Gexp(1)

%[engine, ll] = enter_evidence(dbn_engine, evidence);
%disp(ll);
%margG= marginal_nodes(engine,nodes_map('G.Meal'),1);

%disp(Gexp(1))
%disp(margG.mu)
%disp(sqrt(margG.Sigma))

for measure = 1:140
    %disp(measure);
    evidence{nodes_map('G.obs'),measure} = Gexp(measure); % evidence at time slice 2
    evidence{nodes_map('DG.obs'),measure} = DGexp(measure); % evidence at time slice 2
end
%evidence{nodes_map('I.obs'),measure} = Iexp(measure); 
[engine, ll] = enter_evidence(dbn_engine, evidence);

for measure = 1:140
    margG= marginal_nodes(engine,nodes_map('G.Meal'),measure);
    margGb= marginal_nodes(engine,nodes_map('Gb.Meal'),measure);
    margDG= marginal_nodes(engine,nodes_map('DG.Meal'),measure);
    margS= marginal_nodes(engine,nodes_map('S.Meal'),measure);
   
    %For tabular nodes, we display marg.T(index of node)
    G(end+1,:) = {margG.mu, sqrt(margG.Sigma)};
    Gb(end+1,:) = {margGb.mu, sqrt(margGb.Sigma)};
    DG(end+1,:) = {margDG.mu, sqrt(margDG.Sigma)};
    S(end+1,:) = {margS.mu, sqrt(margS.Sigma)};
    %fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)); % mean +- stddev
end
disp(G);
disp("separate");
disp(Gb);
disp("separate");
disp(DG);
disp("separate");
disp(S);

% Create a table with the data and variable names
T = table(G, S);
% Write data to text file
writetable(T, 'meal_G-S_t2d.txt');


%evidence= cell(npers, T);
%evidence{nodes_map('G.obs'),2} = 90; % evidence at time slice 2
%evidence{nodes_map('I.obs'),measure} = Iexp(measure); 
%[dbn_engine, ll] = enter_evidence(dbn_engine, evidence);
%disp(ll)
%post = zeros(1, T);
%for cpd_id=1:length(dbn.CPD)
%  m = marginal_nodes(dbn_engine, i);
%  post(i) = m.T(2);
%end

% For tabular nodes, we display marg.T(index of node)
%fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)) % mean +- stddev
%Imeal(end+1,:) = {marg.mu, sqrt(marg.Sigma)};

% Parameter estimation from submodels and correct overall network topology
%ncases = 20;
%cases = cell(1, ncases);
%onodes = [nodes_map('G.obs')]; 
%evidence = cell(1, N);

%for i=1:ncases
%  ev = sample_dbn(dbn, 'length', T);
%  cases{i} = cell(npers, T);
  %cases{i}(onodes,:) = ev(onodes, :);
%  cases{i}(onodes,:) = num2cell(Gexp(1:T));
%end

%sample_seq=  cell2mat(sample_dbn(bnet, 'length', T,'evidence', evidence));
% LLtrace is the learning curve: the vector of log-likelihood scores at each iteration.
%disp("checking");
%[dbn2, ll] = learn_params_dbn_em(dbn_engine, cases, "max_iter", 2);
%disp(ll);

%sample_seq=  cell2mat(sample_dbn(dbn2, 'length', T));

%CPD = struct(dbn2.CPD{nodes_map('G.Meal')});
%disp(CPD.mean);
%disp(CPD.cov);
    
%for i=1:T
%    disp(sample_seq(nodes_map('G.Meal'),i));
%    disp(sample_seq(nodes_map('I.Meal'),i));
    %fprintf("CPD %d - trained parameters:\n", i);
    %CPD = struct(dbn2.CPD{nodes_map('G.Meal')});
    %disp(CPD.mean);
    %disp(CPD.cov);
%end
