% Make a DBN for the spatiotemporal BD model with the following variables
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
% Gobs(t), Gobs(t+1)
% Iobs(t), Iobs(t+1)
%
% Time-invariant variables
% lambda k Npatch Nisg Rpbc
%
% All variables display gaussian distributions.

% TODO: check the time step of the equation, and see the G and I.
%%
% Read in the experimental measurements
SPTm1 = importdata('ode_exp1_avr.dat');
Gm1 = SPTm1(:,2); % Gexp in measurement number 1, vector along time
Im1 = SPTm1(:,3); % Iexp in measurement number 1, vector along time
disp(Gm1(1));

%n = 11;

% TODO: To implement model along time
time = 1
[spt_dbn_factory]= make_spt_dbn_factory(Gm1, Im1, time);
[dbn, ~, ~, nodes_map] = create_dbn(spt_dbn_factory);
% parameter learning

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


% Plot ODE.G and ODE.I

figure()
yyaxis left;
plot(1:T, sample_seq(nodes_map('SPT.Gisg'),:));
yyaxis right;
plot(1:T, sample_seq(nodes_map('SPT.Ipm'),:));
legend('SPT.Gisg','SPT.Ipm');  
%%
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
y = sample_seq(nodes_map('SPT.k'),:);
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
dlmwrite('k_prior.txt',variable);

%%
%compute the unconditional marginals of h(20)
dbn_engine = jtree_dbn_inf_engine(dbn);
T = 400;
evidence= cell(npers, T);
[dbn_engine, ll] = enter_evidence(dbn_engine, evidence); % ll is the log marginal likelihood
marg = marginal_nodes(dbn_engine, nodes_map('SPT.k'),10);

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

