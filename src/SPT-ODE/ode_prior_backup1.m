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

% check the time step of the equation, and see the G and I.

% Read in the experimental measurements
odem1 = importdata('ode_exp1_avr.dat');
Gm1 = odem1(:,2); % Gexp in measurement number 1, vector along time
Im1 = odem1(:,3); % Iexp in measurement number 1, vector along time
disp(Gm1(1));

%for time  = 1:length(Gm1)
%    [bnet, nodes_map, intra, inter]= make_ode_bnet(Gm1, Im1,time);
%end
n = 9;

[dbn_factory]= make_ode_dbn_factory(Gm1, Im1);
% parameter learning
npers= bnet.nnodes_per_slice;
T = 200; % lengthhs of sequences to explore
disp(npers);
% Sample from the posterior:
%disp(length(bnet.intra))
%disp(length(bnet.inter))
%evidence=cell(n,T);
%evidence{2,1}=1;
%sample_seq=  cell2mat(sample_dbn(bnet, 'length', T,'evidence', evidence));
sample_seq=  cell2mat(sample_dbn(bnet, 'length', T));
nodes_order= cell2mat(values(nodes_map));

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
xx = linspace(-0.1,0.3);
hold on; 
plot(xx,area*normpdf(xx,mean(y),std(y)),'k-','LineWidth',2);
fprintf("Normal probability density function of beta");
disp(mean(y));
disp(std(y));
%f = ksdensity(y,xx);
%plot(xx,area*f,'g-')
legend('ODE.h, prior');
hold off;

% To generate a conditional gaussian model
function [dbn_factory]= make_ode_dbn_factory(Gm, Im)
    node_names=  horzcat(strcat('ODE.', {'h', 'G', 'G_minus_h', 'I', 'Gref', 'Gexp', 'Iexp'}), 'Reference.I'); % BARAK comment: removed alpha, beta (turned to weights) + added intermediate (G-h)
    n= length(node_names);
    % Intra - in one time slice
    edges_intra1= strcat('ODE.', {'h', 'G_minus_h'; 'G', 'G_minus_h'; 'G', 'Gref'; 'Gref', 'Gexp'}); % BARAK comment: changed in accordance with change in nodes list
    edges_intra2= {'ODE.I','Reference.I'; 'Reference.I','ODE.Iexp'};
    edges_intra= [edges_intra1; edges_intra2];
    % Inter - between time slices
    edges_inter= strcat('ODE.', { 'G', 'G'; 'G_minus_h', 'I'; 'I', 'I' }); % BARAK comment: switched G->I and h->I to (G-h)->I
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('ODE.G')= 'ODE.G.inter';
    eclass2_map('ODE.I')= 'ODE.I.inter';   

    %%
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    %%
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ODE.h', 0, ...
        {'mean', 6.1,   'cov', 0.1} );
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ODE.G', 0, ...
        {'mean', Gm(1), 'cov', 0.1} ); 
    weights_G_minus_h_map_T0= containers.Map(); % parents in slice t
    weights_G_minus_h_map_T1= containers.Map(); % parents in slice t+1
    weights_G_minus_h_map_T0('ODE.G')= 1.0;
    weights_G_minus_h_map_T0('ODE.h')= -1.0;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'ODE.G_minus_h', 0, ...
        {'mean', 0.0, 'cov', 0.00001, 'clamp_mean', 1, 'clamp_cov', 1, 'clamp_weights', 1}, ...
        weights_G_minus_h_map_T0, weights_G_minus_h_map_T1);  % G_minus_h ~ Norm(E(G)-E(h))
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.Gref', 0,   ...
        {'mean', 0.0,   'cov', 0.2, 'weights', 1.0} ); % E= [1.0*G(t+1)] 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.Gexp', 0, ...
        {'mean', 0.0,   'cov', 0.1, 'weights', 1.0}); % E= [1.0*Gref(t+1)]
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.I', 0, ...
        {'mean', Im(1), 'cov', 5.0} ); 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Reference.I', 0, ...
        { 'mean', 0.0,   'cov', 10.0,   'weights', 1.0} ); 
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.Iexp', 0, ...   
        {'mean', 0.0,   'cov', 5.0,   'weights', 1.0} ); 
    
    %%
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    %%
    % CPD for G(t+1)
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.G', 1, ...
        {'mean', 0.0, 'cov', 0.1,  'weights', 1.0} ); % E[ G(t+1) ] = 1.0*G(t)
    % CPD for I(t+1), assume for now all parents are continuous
    % I(t+1) := Normal dist. E = (1-alpha) * I(t) + alpha * beta * (G(t)-h)
    weights_I1_map_T0= containers.Map(); % parents in slice t
    weights_I1_map_T1= containers.Map(); % parents in slice t+1
    INITIAL_ALPHA= 0.05;
    INITIAL_BETA= 0.11;
    weights_I1_map_T0('ODE.I')= 1.0 - INITIAL_ALPHA;
    weights_I1_map_T0('ODE.G_minus_h')= INITIAL_ALPHA * INITIAL_BETA;
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'ODE.I', 1, ...
        {'mean', 0.0, 'cov', 5.0}, ...
        weights_I1_map_T0, weights_I1_map_T1);
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

