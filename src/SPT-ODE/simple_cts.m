% Make a simple PGM with two gaussian nodes
%    G -> I
%

node_names=  { 'G', 'I'};
n = length(node_names);

dag = zeros(n);
dag(1,2)=1;

% node sizes - all cts nodes are scalar, all discrete nodes are binary
ns = ones(1, n);
dnodes = [];
cnodes = mysetdiff(1:n, dnodes);
%ns(dnodes) = 2;

bnet = mk_bnet(dag, ns, 'discrete', dnodes);

% - no parents: Y ~ N(mu, Sigma)
% - cts parents : Y|X=x ~ N(mu + W x, Sigma)
bnet.CPD{1} = gaussian_CPD(bnet, 1, 'mean', 3, 'cov', 0.3);
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', 2, 'cov', 0.1,  'weights', 1);

disp(bnet.CPD{1})
CPD = struct(bnet.CPD{1}); % Peek inside CPD{1}
disp(CPD.mean);

%

% compute the unconditional marginals
engine = jtree_inf_engine(bnet);
evidence = cell(1,n);
[engine, ll] = enter_evidence(engine, evidence);
marg = marginal_nodes(engine, 2);

tol = 1e-2; % tolerance
% evaluates EXPRESSION and, if it is false, displays the error message 'Assertion Failed'
marg.mu
marg.Sigma
assert(approxeq(marg.mu, 4, tol)); % marg.mu equals to 4 +- tol
assert(approxeq(sqrt(marg.Sigma), 0.63, tol)); % marg.Sigma equals to 4 +- tol

% test 2
evidence = cell(1,n);
evidence{1} = 10; % industrial
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, 2);
marg.mu
marg.Sigma

assert(approxeq(marg.mu, 11, tol));
assert(approxeq(sqrt(marg.Sigma), 0.32, tol));

% plot
%gaussplot2d(marg.mu, marg.Sigma);