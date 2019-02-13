% Make a simple PGM with three gaussian nodes in a v-structure
%    BD -> Iexp <- SPT
%

tol = 1e-2; % tolerance for assertion tests
node_names=  { 'BD', 'SPT', 'Iexp'};
n = length(node_names);

dag = zeros(n);
dag(1,3)=1;
dag(2,3)=1;

% node sizes - all cts nodes are scalar, all discrete nodes are binary
ns = ones(1, n);
dnodes = [];
cnodes = mysetdiff(1:n, dnodes);
%ns(dnodes) = 2;

bnet = mk_bnet(dag, ns, 'discrete', dnodes);

% - no parents: Y ~ N(mu, Sigma)
% - cts parents : Y|X=x ~ N(mu + W x, Sigma)
bnet.CPD{1} = gaussian_CPD(bnet, 1, 'mean',  3.0, 'cov', 1.0); % E(BD)=  3.0 +- sqrt(0.3)
bnet.CPD{2} = gaussian_CPD(bnet, 2, 'mean',  -3.0, 'cov', 10.0);%  E(SPT)= 3.0 +- sqrt(0.3)
bnet.CPD{3} = gaussian_CPD(bnet, 3, 'mean',  1.0, 'cov', 0.1, ...
    'weights', [1.0, 1.0]); % E(x2)= 1.0 + BD + SPT  +-  sqrt(0.1)

fprintf("==\nTest1: compute the unconditional marginals\n")
engine = jtree_inf_engine(bnet);
evidence = cell(1,n);
[engine, ll] = enter_evidence(engine, evidence);
marg = marginal_nodes(engine, 1);

% evaluates EXPRESSION and, if it is false, displays the error message 'Assertion Failed'
disp('marginal for BD (parent 1)')
marg = marginal_nodes(engine, 1);
marg
assert(approxeq(marg.mu, 3.0, tol));
assert(approxeq(marg.Sigma, 1.0, tol)); 

disp('marginal for SPT (parent 2)')
marg = marginal_nodes(engine, 2);
marg
assert(approxeq(marg.mu, -3.0, tol));
assert(approxeq(marg.Sigma, 10.0, tol)); 

fprintf('==\nTest2: compute the conditional marginalss given I-exp (the child node)''s value is 1.0\n');
% Test 2:
evidence = cell(1, n);
evidence{3} = 3.0; 
[engine, ll] = enter_evidence(engine, evidence);

disp('marginal for BD (parent 1)')
marg = marginal_nodes(engine, 1);
marg

disp('marginal for SPT (parent 2)')
marg = marginal_nodes(engine, 2);
marg

assert(approxeq(marg.mu, -1.1982, tol));
assert(approxeq(marg.Sigma, 0.991, tol));
