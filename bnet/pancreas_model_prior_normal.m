% Make a DBN for the pancreas model with the following variables
%
% Time-dependent variables
%
% Reference variables
% Scell.ref
% Spancreas.ref
%
% Observed variables
%
% Time-invariant variables
% Scell
% Sislet
% Spancreas
%
% Parameters
%
% To generate a conditional gaussian model
warning('off','MATLAB:singularMatrix');

% Read in the experimental measurements
Experiment3 = importdata('../dataset/meal/140points/meal_exp1_normal2_Gin_dt1_sigmoid.dat');
DGexp = Experiment3(:,2); % Rate of change in glucose intake, mM
Gintakeexp = Experiment3(:,3); % Glucose intake, mM

% Initial values of variables
Scell_Pancreas = 34/300/1e+6 % insulin secretion rate of the pancreas, pM/min

[pancreas_dbn_factory]= make_pancreas_dbn_factory(Scell_Pancreas);
[dbn, ~, ~, nodes_map] = create_dbn(pancreas_dbn_factory);
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
Gvalues = sample_seq(nodes_map('Sislet.Pancreas'),:);
plot(1:T, Gvalues);
yyaxis right;
Ivalues = sample_seq(nodes_map('Spancreas.Pancreas'),:);
plot(1:T, Ivalues);
legend('Sislet.Pancreas','Spancreas.Pancreas');  
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

Scell={}
Sislet={}
Spancreas={}
for measure = 1:100
    %disp(measure);
    %disp(Gintakeexp(measure));
    %evidence{nodes_map('DGintake.obs'),measure} = DGexp(measure); % evidence at time slice 2
end
%evidence{nodes_map('I.obs'),measure} = Iexp(measure); 
[engine, ll] = enter_evidence(dbn_engine, evidence);

for measure = 1:100
    margScell= marginal_nodes(engine,nodes_map('Scell.Pancreas'),measure);
    margSislet= marginal_nodes(engine,nodes_map('Sislet.Pancreas'),measure);
    margSpancreas= marginal_nodes(engine,nodes_map('Spancreas.Pancreas'),measure);

    %For tabular nodes, we display marg.T(index of node)
    Scell(end+1,:) = {margScell.mu, margScell.Sigma, sqrt(margScell.Sigma)};
    Sislet(end+1,:) = {margSislet.mu, margSislet.Sigma, sqrt(margSislet.Sigma)};
    Spancreas(end+1,:) = {margSpancreas.mu, margSpancreas.Sigma, sqrt(margSpancreas.Sigma)};
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
disp(Spancreas);

% Create a table with the data and variable names
T = table(Scell, Sislet, Spancreas);
% Write data to text file
writetable(T, 'pancreas_normal.txt');

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
