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
% Tricks - DG is 
warning('off','MATLAB:singularMatrix');

%Read in the experimental measurements
clear;
% G_Model, I_Model, alpha_Model,beta_Model, Gb_Model
%[meal_dbn_factory]= make_meal_dbn_factory_eq1(5.1, 34, 0.05, 39.6, 2, 5.1);
%Gex_Meal, Y_Meal, alpha_Meal,beta_Meal, K_Meal, dt_Meal, DGex_Meal, Gb_Meal, Sb_Meal
[meal_dbn_factory]= make_meal_dbn_factory_eq2(5.111, 0, 0.05, 39.6,828, 3, 0.056, 5.111,5.111, 34, 25);
[dbn, ~, ~, nodes_map] = create_dbn(meal_dbn_factory);
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);
T = 140; % lengthhs of sequences to explore
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
Gvalues = sample_seq(nodes_map('Gin.Meal'),:);
plot(1:T, Gvalues);
yyaxis right;
Ivalues = sample_seq(nodes_map('I.Meal'),:);
plot(1:T, Ivalues);
legend('Gex.Meal','I.Meal');  
legend('boxoff');

%fprintf("prior distribution of G.Meal, mean %f, error, %f \n", mean(Gvalues), std(Gvalues));
%fprintf("prior distribution of I.Meal, mean %f, error, %f ", mean(Ivalues), std(Ivalues));

% test different measurement
Experiment1 = importdata('/Users/lipingsun/SynologyDrive/research/Meta-modeling/bnt/Meta-modelingS_Gmeal/dataset/meal/140points/data_normal_G-S.dat');
Gexp = Experiment1(:,1); % Gexp in measurement number 1, vector along time
Sexp = Experiment1(:,2); 
Experiment2 = importdata('/Users/lipingsun/SynologyDrive/research/Meta-modeling/bnt/Meta-modelingS_Gmeal/dataset/meal/140points/meal_exp1_normal2.dat');
DGexp = Experiment2(:,2);
Iexp = Experiment2(:,3); 
Experiment3 = importdata('/Users/lipingsun/SynologyDrive/research/Meta-modeling/bnt/Meta-modelingS_Gmeal/dataset/meal/140points/meal_exp1_normal2_Gin.dat');
DGexp = Experiment3(:,2);
Ginexp = Experiment3(:,3); 

% testing
%CPD = struct(dbn.CPD{nodes_map('G.ref')});
%disp('here');
%disp(CPD.mean);
%disp(CPD.cov);

%T=10
evidence= cell(npers, T);
%evidence{nodes_map('G.obs'),1} = Gexp(1)

% Parameter estimation from submodels and correct overall network topology
ncases = 2;
cases = cell(1, ncases);
onodes1 = [nodes_map('S.obs')]; 
onodes2 = [nodes_map('I.obs')];
onodes3 = [nodes_map('Gex.obs')];
onodes4 = [nodes_map('DGex.obs')];
onodes5 = [nodes_map('Gin.obs')];

%onodes = 1:npers
evidence = cell(1, npers);

for measure = 1:140
    %disp(measure);
    %disp(Ginexp(measure));
    evidence{nodes_map('S.obs'),measure} = Sexp(measure); % evidence at time slice 2
    evidence{nodes_map('I.obs'),measure} = Iexp(measure); % evidence at time slice 2
    evidence{nodes_map('Gex.obs'),measure} = Gexp(measure); % evidence at time slice 2
    evidence{nodes_map('DGex.obs'),measure} = DGexp(measure); % evidence at time slice 2
    evidence{nodes_map('Gin.obs'),measure} = Ginexp(measure); % evidence at time slice 2
end

for i=1:ncases
  %ev = sample_dbn(dbn, 'length', T);
  cases{i} = cell(npers, T);
  %cases{i}(onodes,:) = ev(onodes, :);
  cases{i}(onodes1,:) = num2cell(Sexp(1:T));
  cases{i}(onodes2,:) = num2cell(Iexp(1:T));
  cases{i}(onodes3,:) = num2cell(Gexp(1:T));
  cases{i}(onodes4,:) = num2cell(DGexp(1:T));
  cases{i}(onodes5,:) = num2cell(Ginexp(1:T));
end

sample_seq=  cell2mat(sample_dbn(dbn, 'length', T,'evidence', evidence));
% LLtrace is the learning curve: the vector of log-likelihood scores at each iteration.
disp("checking");
[dbn2, ll] = learn_params_dbn_em(dbn_engine, cases, "max_iter", 2);
disp(ll);

fprintf("CPD %d - trained parameters:\n", i);
Gex_CPD = struct(dbn2.CPD{nodes_map('Gex.Meal')});
disp(Gex_CPD.mean);
disp(Gex_CPD.cov);
I_CPD = struct(dbn2.CPD{nodes_map('I.Meal')});
disp(I_CPD.mean);
disp(I_CPD.cov);

sample_seq=  cell2mat(sample_dbn(dbn2, 'length', T));
    
for i=T-1:T
    disp(sample_seq(nodes_map('Gex.Meal'),i));
    disp(sample_seq(nodes_map('I.Meal'),i));
end


