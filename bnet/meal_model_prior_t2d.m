% Make a DBN for the meal model with the following variables
%
% Time-dependent variables
%  -> Gex.Meal(t)  ->  Gex.Meal(t+1) ->
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
warning('off','MATLAB:singularMatrix');

% Read in the experimental measurements
Experiment3 = importdata('../dataset/meal/140points/meal_exp1_normal2_Gin_dt1_sigmoid.dat');
DGexp = Experiment3(:,2); % Rate of change in glucose intake, mM
Gintakeexp = Experiment3(:,3); % Glucose intake, mM

% Initial values of variables
Gex_Meal = 9.167; % Basal plasma glucose level, mM
Y_Meal = 0; % provision of new insulin to the beta-cells at t=0, pM/min
alpha_Meal = 0.013; % Delay between the glucose signal and insulin secretion, min^{-1}
beta_Meal = 22.5; % Pancreatic responsivity to glucose, pM/min per mM
gamma_Meal = 0.5; % transfer rate constant between portal vein and liver, min^{-1}
k1_Meal = 0.0001; % Coefficient for I reducing glucose level, significantly affect the periodicity of the curves.
k2_Meal = 0.00001; % Coefficient for Gex reducing glucose level
K_Meal = 10; % Pancreatic responsivity to the glucose rate of change, pmol/L per mM
dt_Meal_min = 1; % time scale, min
Gb_Meal = 9.167; % Basal plasma glucose level, mM
DGintake_Meal = DGexp(1); % Rate of glucose intake from food at t = 0, mM/min
Sb_Meal = 102.5; %Basal insulin secretion, pM/min
I_Meal = 52; % Basal plasma insulin level, pM

[meal_dbn_factory]= make_meal_dbn_factory_eq(Gex_Meal, Y_Meal, alpha_Meal,beta_Meal, gamma_Meal, k1_Meal, k2_Meal, K_Meal, dt_Meal_min, Gb_Meal, DGintake_Meal, Sb_Meal, I_Meal);
[dbn, ~, ~, nodes_map] = create_dbn(meal_dbn_factory);
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
Gvalues = sample_seq(nodes_map('Gex.Meal'),:);
plot(1:T, Gvalues);
yyaxis right;
Ivalues = sample_seq(nodes_map('I.Meal'),:);
plot(1:T, Ivalues);
legend('Gex.Meal','I.Meal');  
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

G={}
I={}
S={}
for measure = 1:420
    %disp(measure);
    %disp(Gintakeexp(measure));
    evidence{nodes_map('DGintake.obs'),measure} = DGexp(measure); % evidence at time slice 2
end
%evidence{nodes_map('I.obs'),measure} = Iexp(measure); 
[engine, ll] = enter_evidence(dbn_engine, evidence);

for measure = 1:420
    margG= marginal_nodes(engine,nodes_map('Gex.Meal'),measure);
    margS= marginal_nodes(engine,nodes_map('S.Meal'),measure);
    margI= marginal_nodes(engine,nodes_map('I.Meal'),measure);

    %For tabular nodes, we display marg.T(index of node)
    G(end+1,:) = {margG.mu, margG.Sigma, sqrt(margG.Sigma)};
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
T = table(G, S, I);
% Write data to text file
writetable(T, 'meal_t2d.txt');

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
%nodes = [nodes_map('G.obs')]; 
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
