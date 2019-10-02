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

[meal_dbn_factory,n]= make_meal_dbn_factory(90, 25, 0.05, 0.11, 6.1);
[dbn, ~, ~, nodes_map] = create_dbn(meal_dbn_factory);
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);
T = 50; % lengthhs of sequences to explore
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
Experiments = importdata('meal_exp1_avr.dat');
Gexp = Experiments(:,2); % Gexp in measurement number 1, vector along time

% Parameter estimation from submodels and correct overall network topology
%ncases = 200;
%fprintf(strcat("Learning CPD parameters from %d instances of %d frames time-series of model 1\n", ...
%    " and %d instances of model 2, with correct network topology\n"), ...
%    ncases/2, T, ncases/2);
%cases = cell(1, ncases);
%onodes = [nodes_map('G.obs')]; 
%for i=1:ncases
%  ev = sample_dbn(dbn, T);
%  cases{i} = cell(n, T);
%  cases{i}(onodes,:) = ev(onodes, :);
%  marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
%                     nodes_map('I.Meal')+npers, ...
%                     i);
%  fprintf("%f sigma %f +- %f", marg.mu, marg.Sigma, sqrt(marg.Sigma)) % mean +- stddev
%  Imeal(end+1,:) = {marg.mu, sqrt(marg.Sigma)};
%end

for measure = 1:200  
    
Imeal={}    
evidence= cell(npers, T);
evidence{nodes_map('G.obs'),2} = Gexp(measure); 
i=1;
disp(i);
marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
                     nodes_map('I.Meal')+npers, ...
                     i);
fprintf("%f sigma %f +- %f", marg.mu, sqrt(marg.Sigma)) % mean +- stddev
Imeal(end+1,:) = {marg.mu, sqrt(marg.Sigma)};
end

% Create a table with the data and variable names
T = table(Gexp, Imeal )
% Write data to text file
writetable(T, 'mean_I.txt')

