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
% The ODE equation is only when G > h
% TODO: check the time step of the equation, and see the G and I.
warning('off','MATLAB:singularMatrix');

%Read in the experimental measurements
clear;
odem1 = importdata('ode_exp1_avr.dat');
Gm1 = odem1(:,2); % Gexp in measurement number 1, vector along time
Im1 = odem1(:,3); % Iexp in measurement number 1, vector along time
%disp(Gm1(1));

% TODO: To implement model along time
time = 1
[ode_dbn_factory]= make_ode_dbn_factory(Gm1, Im1, time);
[dbn, ~, ~, nodes_map] = create_dbn(ode_dbn_factory);
% parameter learning
npers= dbn.nnodes_per_slice;
T = 400; % lengthhs of sequences to explore
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
plot(1:T, sample_seq(nodes_map('ODE.Gpm'),:));
yyaxis right;
plot(1:T, sample_seq(nodes_map('ODE.Ipm'),:));
legend('ODE.Gpm','ODE.Ipm');  


%Plot ODE.h distribution
disp('plot');
fprintf("Sampled time-series of length %d", T);
y = sample_seq(nodes_map('ODE.h'),:);
nbins = 10;

[hts,ctrs] = hist(y, nbins);
%h = bar(ctrs,hts,'hist');
%set(h,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0 0 0]);
area = sum(hts) * (ctrs(2)-ctrs(1));
xx = linspace(2,10);
hold on; 
plot(xx,area*normpdf(xx,mean(y),std(y))/T,'k-','LineWidth',2);
fprintf("Normal probability density function of ODE.h");
disp(mean(y));
disp(std(y));
%f = ksdensity(y,xx);
%plot(xx,area*f,'g-')
set(gca,'linewidth', 2,'fontsize',24,'fontname','Times New Roman') % Sets the width of the axis lines, font size, font
gca.XAxis.MinorTick = 'on';
gca.XAxis.MinorTickValues = xx(1):.5:xx(2);
legend('Meal.h, prior');
legend('boxoff');
hold off;

% Create a table with the data and variable names
%yy=area*normpdf(xx,mean(y),std(y))/sum(area*normpdf(xx,mean(y),std(y)))
%variable = [xx(:) yy(:)];
%size(variable);
%dlmwrite('h_prior.txt',variable);
%type 'h_prior.txt'

%compute the unconditional marginals of h(20)
disp(npers)
dbn_engine = jtree_dbn_inf_engine(dbn);
T = 400;
evidence= cell(npers, T);
[dbn_engine, ll] = enter_evidence(dbn_engine, evidence); % ll is the log marginal likelihood
marg = marginal_nodes(dbn_engine, nodes_map('ODE.h'),10);

tol = 1; % tolerance
% evaluates EXPRESSION and, if it is false, displays the error message 'Assertion Failed'
fprintf("Unconditional probability distribution of h(10) is:\n"); 
fprintf("%f sigma %f +- %f", marg.mu, marg.Sigma, sqrt(marg.Sigma)) % mean +- stddev

%assert(approxeq(marg.mu, 6.1, tol)); % marg.mu equals to 4 +- tol
%assert(approxeq(marg.Sigma, 3.0, tol)); % marg.Sigma equals to 4 +- tol

%Posterior marginal of h(20) given Iexp
T = 400;
evidence= cell(npers, T);
%evidence{nodes_map('ODE.Iexp'),10} = [40.0 50.0]; % we may use soft
%evidence if the observed node has some distributions over its values
%soft_evidence{nodes_map('ODE.Iexp'),10} = [0.6 0.4]; 
%[dbn_engine, ll] = enter_evidence(dbn_engine, evidence, 'soft', soft_evidence);
evidence{nodes_map('E.Ipm'),2} = 25.0; 
i=10
disp(i)
marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
                     nodes_map('ODE.h')+npers, ...
                     i);
fprintf("Posterior probability distribution of h(10) given Iexp(10) is:\n"); 
fprintf("%f sigma %f +- %f", marg.mu, marg.Sigma, sqrt(marg.Sigma)) % mean +- stddev
xx = linspace(4,8);
plot(xx,area*normpdf(xx,marg.mu,marg.Sigma),'k-','LineWidth',2);


% Create a table with the data and variable names
yy=normpdf(xx,marg.mu,marg.Sigma)/sum(normpdf(xx,marg.mu,marg.Sigma))
variable = [xx(:) yy(:)];
size(variable);
dlmwrite('h_posterior.txt',variable);

%assert(approxeq(marg.mu, 6.1, tol)); % marg.mu equals to 4 +- tol
%assert(approxeq(sqrt(marg.Sigma), 0.1, tol)); % marg.Sigma equals to 4 +- tol

%compute other posteror marginals with evidence
%T = 400;
%n=nodes_map.Count;
%evidence= cell(npers, T);
%evidence{nodes_map('ODE.Gpm_minus_h'),10} = 1.0; 
%evidence{nodes_map('ODE.Gpm'),10} = 70.0; 

%i=12
%disp(i)
%marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
%                     nodes_map('ODE.h')+npers, ...
%                     i);
%fprintf("Posterior probability distribution of h is:\n"); 
%fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)) % mean +- stddev

