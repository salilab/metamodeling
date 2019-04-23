% Make a meta-modeling DBN for the ODE and SPT model
warning('off','MATLAB:singularMatrix');

%Read in the experimental measurements
clear;
% Read in the experimental measurements
odem1 = importdata('ode_exp1_avr.dat');
Gm1 = odem1(:,2); % Gexp in measurement number 1, vector along time
Im1 = odem1(:,3); % Iexp in measurement number 1, vector along time

sptm1 = importdata('ode_exp1_avr.dat');
Go1 = sptm1(:,2); % Gexp in measurement number 1, vector along time
Io1 = sptm1(:,3); % Iexp in measurement number 1, vector along time

%for time  = 1:length(Gm1)
%    [bnet, nodes_map] = make_meta_bnet(Gm1, Im1, Go1, Io1, time);
%end
time = 1;
[dbn, nodes_map] = make_meta_bnet(Gm1, Im1, Go1, Io1, time);

% parameter learning
npers= dbn.nnodes_per_slice;
T = 400; % lengthhs of sequences to explore
disp(npers);
% Sample from the posterior:
%disp(length(bnet.intra))
%disp(length(bnet.inter))
%evidence=cell(n,T);
%evidence{2,1}=1;
%sample_seq=  cell2mat(sample_dbn(bnet, 'length', T,'evidence', evidence));
sample_seq=  cell2mat(sample_dbn(dbn, 'length', T));
nodes_order= cell2mat(values(nodes_map));

% Plot the distribution of ODE.beta
disp('plot');
fprintf("Sampled time-series of length %d", T);
y = sample_seq(nodes_map('ODE.h'),:);
nbins = 10;

[hts,ctrs] = hist(y, nbins);
%h = bar(ctrs,hts,'hist');
%set(h,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0 0 0]);
area = sum(hts) * (ctrs(2)-ctrs(1));
%xx = linspace(4,8);
%hold on; 
%plot(xx,area*normpdf(xx,mean(y),std(y)),'k-','LineWidth',2);
%fprintf("Normal probability density function of ODE.h");
%disp(mean(y));
%disp(std(y));
%f = ksdensity(y,xx);
%plot(xx,area*f,'g-')
%legend('ODE.h, posterior');
%hold off;

%compute the unconditional marginals of h(20)
disp(npers)
dbn_engine = jtree_dbn_inf_engine(dbn);
%T = 400;
%evidence= cell(npers, T);
%[dbn_engine, ll] = enter_evidence(dbn_engine, evidence); % ll is the log marginal likelihood
%marg = marginal_nodes(dbn_engine, nodes_map('ODE.h'),10);

%tol = 1; % tolerance
% evaluates EXPRESSION and, if it is false, displays the error message 'Assertion Failed'
%fprintf("Unconditional probability distribution of h(10) is:\n"); 
%fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)) % mean +- stddev

%assert(approxeq(marg.mu, 6.1, tol)); % marg.mu equals to 4 +- tol
%assert(approxeq(marg.Sigma, 3.0, tol)); % marg.Sigma equals to 4 +- tol

%Posterior marginal of h(20) given Iexp

T = 400;
evidence= cell(npers, T);
evidence{nodes_map('E.Ipm'),2} = 85.0; 
i=10
disp(i)
marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
                     nodes_map('ODE.h')+npers, ...
                     i);
fprintf("Posterior probability distribution of h(10) given Iexp(10) is:\n"); 
fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)) % mean +- stddev

% plot h
%xx = linspace(4.5,7.5);

% create some plot with a legend
%hAx(1) = axes();
%hLine(1) = plot(xx,area*normpdf(xx,6.065716,0.299999),'k-','LineWidth',2, 'Parent',hAx(1));
%set(hAx(1), 'Box','off');
%legend(hLine(1), {'ODE.h, prior'},'Location','NorthWest');

% copy the axis
%hAx(2) = copyobj(hAx(1),gcf);
%delete( get(hAx(2),'Children') )            %# delete its children
%hLine(2) = plot(xx,area*normpdf(xx,marg.mu,marg.Sigma),'r-','LineWidth',2,'Parent',hAx(2));
%set(hAx(2), 'Color','none', 'XTick',[], ...
%    'YAxisLocation','right', 'Box','off');   %# make it transparent
%legend(hLine(2), {'ODE.h, posterior'}, 'Color','w');

%Posterior marginal of h(20) given Iexp

T = 400;
evidence= cell(npers, T);
evidence{nodes_map('E.Ipm'),2} = 25.0; 
i=10
disp(i)
marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
                     nodes_map('SPT.k')+npers, ...
                     i);
fprintf("Posterior probability distribution of h(10) given Iexp(10) is:\n"); 
fprintf("%f +- %f", marg.mu, sqrt(marg.Sigma)) % mean +- stddev

% plot k
xx = linspace(8,11);

% create some plot with a legend
hAx(1) = axes();
hLine(1) = plot(xx,area*normpdf(xx,9.023592,0.109810),'k-','LineWidth',2, 'Parent',hAx(1));
set(hAx(1), 'Box','off');
legend(hLine(1), {'SPT.k, prior'},'Location','NorthWest');

% copy the axis
hAx(2) = copyobj(hAx(1),gcf);
delete( get(hAx(2),'Children') )            %# delete its children
hLine(2) = plot(xx,area*normpdf(xx,marg.mu,marg.Sigma),'r-','LineWidth',2,'Parent',hAx(2));
set(hAx(2), 'Color','none', 'XTick',[], ...
    'YAxisLocation','right', 'Box','off');   %# make it transparent
legend(hLine(2), {'SPT.k, posterior'}, 'Color','w');

function [meta_dbn, nodes_map]=make_meta_bnet(Gm, Im, Go, Io, time)

    % make ODE and BD models
    [ode_dbn_factory]= ...
        make_ode_dbn_factory(Gm, Im, time);
    [spt_dbn_factory]= ...
        make_spt_dbn_factory(Go, Io, time);
    meta_dbn_factory= ...
        merge_dbn_factories(ode_dbn_factory, ...
                              spt_dbn_factory);
    weights_Iref_map_T0= containers.Map(); % parents in slice t
    weights_Iref_map_T1= containers.Map(); % parents in slice t
    weights_Iref_map_T0('ODE.Ipm')= 0.5;
    weights_Iref_map_T0('SPT.Ipm')= 0.5;
    CPDFactory_Iref = ...
        CPDFactory('Gaussian_CPD', 'R.Ipm', 0, ...
        {'mean', 0.0, 'cov', 5.0}, ...
        weights_Iref_map_T0, ...
        weights_Iref_map_T1);
    add_CPD_factories(meta_dbn_factory, {CPDFactory_Iref}, false);

    [meta_dbn, ~, ~, nodes_map] = create_dbn(meta_dbn_factory);
    
end
