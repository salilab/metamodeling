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
% G_Model, S_Model, k_Model, Nisg_Model, Npatch_Model, Ninsulin_Model
[spt_dbn_factory]= make_spt_dbn_factory(0.1, 34, 10, 300, 6, 1.8e-6);
[dbn, ~, ~, nodes_map] = create_dbn(spt_dbn_factory);
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

% Plot SPT.G and SPT.I
%figure()
%yyaxis left;
%plot(1:T, sample_seq(nodes_map('SPT.G'),:));
%yyaxis right;
%plot(1:T, sample_seq(nodes_map('SPT.I'),:));
%legend('SPT.G','SPT.I');  

% Plot the distribution of SPT.beta
disp('plot');
fprintf("Sampled time-series of length %d", T);
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
disp(std(y));
%f = ksdensity(y,xx);
%plot(xx,area*f,'g-')
legend('SPT.k, prior');
hold off;
% Create a table with the data and variable names

yy=area*normpdf(xx,mean(y),std(y))/sum(area*normpdf(xx,mean(y),std(y)))
variable = [xx(:) yy(:)];
size(variable);
dlmwrite('spt_k_prior_normal.txt',variable);

%compute the unconditional marginals of h(20)
dbn_engine = jtree_dbn_inf_engine(dbn);
T = 400;
evidence= cell(npers, T);
[dbn_engine, ll] = enter_evidence(dbn_engine, evidence); % ll is the log marginal likelihood
marg = marginal_nodes(dbn_engine, nodes_map('k.SPT'),400);

tol = 1; % tolerance
% evaluates EXPRESSION and, if it is false, displays the error message 'Assertion Failed'
fprintf("Unconditional probability distribution of k(10) is:\n"); 
fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)) % mean +- stddev

%assert(approxeq(marg.mu, 6.1, tol)); % marg.mu equals to 4 +- tol
%assert(approxeq(marg.Sigma, 3.0, tol)); % marg.Sigma equals to 4 +- tol

%Posterior marginal of h(20) given Iexp
T = 400;
evidence{nodes_map('E.Ipm'),2} = 25.0; 
i=10
disp(i)
marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
                     nodes_map('SPT.k')+npers, ...
                     i);
fprintf("Posterior probability distribution of h(10) given Iexp(10) is:\n"); 
fprintf("%f sigma %f +- %f", marg.mu, marg.Sigma, sqrt(marg.Sigma)) % mean +- stddev
xx = linspace(8,12,100);
plot(xx,area*normpdf(xx,marg.mu,marg.Sigma),'k-','LineWidth',2);
legend('SPT.k, posterior');

% Create a table with the data and variable names
yy=normpdf(xx,marg.mu,marg.Sigma)/sum(normpdf(xx,marg.mu,marg.Sigma))
variable = [xx(:) yy(:)];
size(variable);
dlmwrite('k_posterior.txt',variable);

%compute other posteror marginals with evidence
%T = 400;
%n=nodes_map.Count;
%evidence= cell(npers, T);
%evidence{nodes_map('SPT.Ipm'),10} = 1.0; 
%evidence{nodes_map('SPT.k'),10} = 70.0; 

%i=12
%disp(i)
%marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
%                     nodes_map('SPT.k')+npers, ...
%                     i);
%fprintf("Posterior probability distribution of k is:\n"); 
%fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)) % mean +- stddev