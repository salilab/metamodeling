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
bnet.CPD{1} = gaussian_CPD(bnet, 1, 'mean', 3.0, 'cov', 0.3);
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean', 1.0, 'cov', 0.1,  'weights', 0.8); % E(x2)= 1.0 + 0.8*x1

fprintf('Test 1: unconditional marginals')
engine = jtree_inf_engine(bnet);
evidence = cell(1,n);
[engine, ll] = enter_evidence(engine, evidence);
marg = marginal_nodes(engine, 2);

tol = 1e-3; % tolerance
% evaluates EXPRESSION and, if it is false, displays the error message 'Assertion Failed'
marg
assert(approxeq(marg.mu, 3.4, tol)); 
assert(approxeq(marg.Sigma, 0.3*0.8^2 + 0.1, tol)); 

fprintf('Test 2: conditional probability of x2 given x1=1.0')
evidence = cell(1, n);
evidence{1} = 1; % industrial
[engine, ll] = enter_evidence(engine, evidence);

marg = marginal_nodes(engine, 2);
marg


assert(approxeq(marg.mu, 1.8, tol));
assert(approxeq(marg.Sigma, 0.1, tol));
