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
[meal_dbn_factory]= make_meal_dbn_factory(165, 52, 0.013, 0.05, 165);
[dbn, ~, ~, nodes_map] = create_dbn(meal_dbn_factory);
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);
T = 200; % lengthhs of sequences to explore
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
figure()
yyaxis left;
Gvalues = sample_seq(nodes_map('G.Meal'),:);
plot(1:T, Gvalues);
yyaxis right;
Ivalues = sample_seq(nodes_map('I.Meal'),:);
plot(1:T, Ivalues);
legend('G.Meal','I.Meal');  
legend('boxoff');

fprintf("prior distribution of G.Meal, mean %f, error, %f \n", mean(Gvalues), std(Gvalues));
fprintf("prior distribution of I.Meal, mean %f, error, %f ", mean(Ivalues), std(Ivalues));

% test different measurement
Experiments = importdata('dataset/meal/meal_exp1_t2b.dat');
Gexp = Experiments(:,2); % Gexp in measurement number 1, vector along time
Iexp = Experiments(:,3); % Gexp in measurement number 1, vector along time
% Parameter estimation from submodels and correct overall network topology
G={}
Gb={}
G_minus_Gb={}
I={}

for measure = 1:198
    %disp(measure);
    evidence= cell(npers, T);
    evidence{nodes_map('G.obs'),measure} = Gexp(measure); % evidence at time slice 2
    %evidence{nodes_map('I.obs'),measure} = Iexp(measure); 
    [engine, ll] = enter_evidence(dbn_engine, evidence);
    %disp(ll);
    margG= marginal_nodes(engine,nodes_map('G.Meal')+npers,measure+1);
    margGb= marginal_nodes(engine,nodes_map('Gb.Meal')+npers,measure+1);
    margG_minus_Gb= marginal_nodes(engine,nodes_map('G_minus_Gb.Meal')+npers,measure+1);
    margI= marginal_nodes(engine,nodes_map('I.Meal')+npers,measure+1);
   
    %For tabular nodes, we display marg.T(index of node)
    G(end+1,:) = {margG.mu, sqrt(margG.Sigma)};
    Gb(end+1,:) = {margGb.mu, sqrt(margGb.Sigma)};
    G_minus_Gb(end+1,:) = {margG_minus_Gb.mu, sqrt(margG_minus_Gb.Sigma)};
    I(end+1,:) = {margI.mu, sqrt(margI.Sigma)};
    %fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)); % mean +- stddev
end
disp(G);
disp("separate");
disp(Gb);
disp("separate");
disp(G_minus_Gb);
disp("separate");
disp(I);


% Create a table with the data and variable names
T = table(G, I);
% Write data to text file
writetable(T, 'meal_G-I_t2b.txt');
  
