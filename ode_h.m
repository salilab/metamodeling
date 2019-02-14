% Make a DBN for the ODE model with the following variables
%
% Time-dependent variables
%  -> G(t)  ->  G(t+1) ->
%  -> I(t)  ->  I(t+1) ->
%
% Reference variables
% Gref(t), Gref(t+1)
% Iref(t), Iref(t+1)
%
% Observed variables
% Gexp(t), Gexp(t+1)
% Iexp(t), Iexp(t+1)
%
% Time-invariant variables
% alpha beta h 
%

% TODO: check the time step of the equation, and see the G and I.
%%
% Read in the experimental measurements
odem1 = importdata('ode_exp1_avr.dat');
Gm1 = odem1(:,2); % Gexp in measurement number 1, vector along time
Im1 = odem1(:,3); % Iexp in measurement number 1, vector along time
disp(Gm1(1));

%for time  = 1:length(Gm1)
%    [bnet, nodes_map, intra, inter]= make_ode_bnet(Gm1, Im1,time);
%end
n = 9;

% TODO: To implement model along time
time = 1
[dbn_factory]= make_ode_dbn_factory(Gm1, Im1, time);
[dbn, ~, ~, nodes_map] = create_dbn(dbn_factory);

% parameter learning
npers= bnet.nnodes_per_slice;
T = 400; % lengthhs of sequences to explore
disp(npers);
% Sample from the posterior:
%disp(length(bnet.intra))
%disp(length(bnet.inter))
%evidence=cell(n,T);
%evidence{2,1}=1;
%sample_seq=  cell2mat(sample_dbn(bnet, 'length', T,'evidence', evidence));
sample_seq=  cell2mat(sample_dbn(bnet, 'length', T));
nodes_order= cell2mat(values(nodes_map));
%%
% Plot ODE.G and ODE.I
%figure()
%yyaxis left;
%plot(1:T, sample_seq(nodes_map('ODE.G'),:));
%yyaxis right;
%plot(1:T, sample_seq(nodes_map('ODE.I'),:));
%legend('ODE.G','ODE.I');  

% Plot the distribution of ODE.beta
disp('plot');
fprintf("Sampled time-series of length %d", T);
y = sample_seq(nodes_map('ODE.h'),:);
nbins = 10;

[hts,ctrs] = hist(y, nbins);
h = bar(ctrs,hts,'hist');
set(h,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0 0 0]);
area = sum(hts) * (ctrs(2)-ctrs(1));
xx = linspace(4,8);
hold on; 
plot(xx,area*normpdf(xx,mean(y),std(y)),'k-','LineWidth',2);
fprintf("Normal probability density function of ODE.h");
disp(mean(y));
disp(std(y));
%f = ksdensity(y,xx);
%plot(xx,area*f,'g-')
legend('ODE.h, prior');
hold off;

%%
% compute the unconditional marginals
engine = jtree_inf_engine(bnet);
T = 2;
evidence= cell(npers, T);
[engine, ll] = enter_evidence(engine, evidence); % ll is the log marginal likelihood
marg = marginal_nodes(engine, nodes_map('ODE.h'));

tol = 1; % tolerance
% evaluates EXPRESSION and, if it is false, displays the error message 'Assertion Failed'
marg.mu
marg.Sigma
assert(approxeq(marg.mu, 6.1, tol)); % marg.mu equals to 4 +- tol
assert(approxeq(sqrt(marg.Sigma), 0.1, tol)); % marg.Sigma equals to 4 +- tol
%%
% compute the posteror marginals with evidence
T = 100;
evidence= cell(npers, T);
evidence{nodes_map('ODE.G_minus_h'),1} = 30; 
marg= marginal_nodes(enter_evidence(engine, evidence), ...
                     nodes_map('ODE.h'), ...
                     1);
fprintf("Posterior probability distribution of I given G==2 is:\n"); 

% evaluates EXPRESSION and, if it is false, displays the error message 'Assertion Failed'
%tol = 1; % tolerance
marg.mu % mean of the marginal
marg.Sigma % variance of the marginal
%assert(approxeq(marg.mu, 6.1, tol)); % marg.mu equals to 4 +- tol
%assert(approxeq(sqrt(marg.Sigma), 0.1, tol)); % marg.Sigma equals to 4 +- tol