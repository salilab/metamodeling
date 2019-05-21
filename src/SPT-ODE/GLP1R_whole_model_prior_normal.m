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
% GLP1_Model,GLP1R_Model, Galpha_Model, cAMP_Model
[glp1r_dbn_factory]= make_glp1r_dbn_factory_whole(0.1, 0.1, 0.1, 0.1);
[dbn, ~, ~, nodes_map] = create_dbn(glp1r_dbn_factory);
% parameter learning
npers= dbn.nnodes_per_slice;
T = 400; % lengthhs of sequences to explore
%disp(npers);

%Sample from the posterior:
%disp(length(bnet.intra))
%disp(length(bnet.inter))
%evidence=cell(n,T);
%evidence{2,1}=1;
%sample_seq=  cell2mat(sample_dbn(bnet, 'length', T,'evidence', evidence));
sample_seq=  cell2mat(sample_dbn(dbn, 'length', T));
nodes_order= cell2mat(values(nodes_map));

% Plot GLP1R.G and GLP1R.I
%figure()
%yyaxis left;
%plot(1:T, sample_seq(nodes_map('GLP1R.G'),:));
%yyaxis right;
%plot(1:T, sample_seq(nodes_map('GLP1R.I'),:));
%legend('GLP1R.G','GLP1R.I');  

% Plot the distribution of GLP1R.beta
disp('plot');
fprintf("Sampled time-series of length %d", T);
y = sample_seq(nodes_map('cAMP.GLP1R'),:);
nbins = 10;

[hts,ctrs] = hist(y, nbins);
h = bar(ctrs,hts,'hist');
set(h,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0 0 0]);
area = sum(hts) * (ctrs(2)-ctrs(1));
xx = linspace(0,0.2,100);
hold on; 
plot(xx,area*normpdf(xx,mean(y),std(y)),'k-','LineWidth',2);
fprintf("Normal probability density function of GLP1R.k");
disp(mean(y));
disp(std(y));
%f = ksdensity(y,xx);
%plot(xx,area*f,'g-')
legend('cAMP.GLP1R, prior');
hold off;
% Create a table with the data and variable names

yy=area*normpdf(xx,mean(y),std(y))/sum(area*normpdf(xx,mean(y),std(y)))
variable = [xx(:) yy(:)];
size(variable);
dlmwrite('GLP1R_cAMP_prior_normal.txt',variable);

%compute the unconditional marginals of h(20)
dbn_engine = jtree_dbn_inf_engine(dbn);
T = 400;
evidence= cell(npers, T);
[dbn_engine, ll] = enter_evidence(dbn_engine, evidence); % ll is the log marginal likelihood
marg = marginal_nodes(dbn_engine, nodes_map('cAMP.GLP1R'),400);

tol = 1; % tolerance
% evaluates EXPRESSION and, if it is false, displays the error message 'Assertion Failed'
fprintf("Unconditional probability distribution of k(10) is:\n"); 
fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)) % mean +- stddev

%assert(approxeq(marg.mu, 6.1, tol)); % marg.mu equals to 4 +- tol
%assert(approxeq(marg.Sigma, 3.0, tol)); % marg.Sigma equals to 4 +- tol

%Posterior marginal of h(20) given Iexp
T = 400;
evidence{nodes_map('cAMP.ref'),2} = 1.0; 
i=10
disp(i)
marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
                     nodes_map('cAMP.GLP1R')+npers, ...
                     i);
fprintf("Posterior probability distribution of h(10) given Iexp(10) is:\n"); 
fprintf("%f sigma %f +- %f", marg.mu, marg.Sigma, sqrt(marg.Sigma)) % mean +- stddev
xx = linspace(8,12,100);
plot(xx,area*normpdf(xx,marg.mu,marg.Sigma),'k-','LineWidth',2);
legend('GLP1R.k, posterior');

% Create a table with the data and variable names
yy=normpdf(xx,marg.mu,marg.Sigma)/sum(normpdf(xx,marg.mu,marg.Sigma))
variable = [xx(:) yy(:)];
size(variable);
dlmwrite('cAMP_posterior.txt',variable);

%compute other posteror marginals with evidence
%T = 400;
%n=nodes_map.Count;
%evidence= cell(npers, T);
%evidence{nodes_map('GLP1R.Ipm'),10} = 1.0; 
%evidence{nodes_map('GLP1R.k'),10} = 70.0; 

%i=12
%disp(i)
%marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
%                     nodes_map('GLP1R.k')+npers, ...
%                     i);
%fprintf("Posterior probability distribution of k is:\n"); 
%fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)) % mean +- stddev