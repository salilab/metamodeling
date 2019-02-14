% Make a meta-modeling DBN for the ODE and SPT model

% Read in the experimental measurements
odem1 = importdata('ode_exp1_avr.dat');
Gm1 = odem1(:,2); % Gexp in measurement number 1, vector along time
Im1 = odem1(:,3); % Iexp in measurement number 1, vector along time

sptm1 = importdata('spt_obs1_avr.dat');
Go1 = sptm1(:,2); % Gexp in measurement number 1, vector along time
Io1 = sptm1(:,3); % Iexp in measurement number 1, vector along time

%for time  = 1:length(Gm1)
%    [bnet, nodes_map] = make_meta_bnet(Gm1, Im1, Go1, Io1, time);
%end
[bnet, nodes_map] = make_meta_bnet(Gm1, Im1, Go1, Io1, 1);

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

% I.
T = 2
engine = jtree_dbn_inf_engine(bnet); 
evidence= cell(npers, T);
evidence{nodes_map('SPT.k'), 1} = 10; 
marg= marginal_nodes(enter_evidence(engine, evidence), ...
                     nodes_map('ODE.h'), ...
                     2);
fprintf("Posterior probability distribution of ODE.h(time =2) given SPT.k(time=1)== 100is:\n"); 
disp(marg.T)

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
legend('ODE.h, posterior');
hold off;
time =3
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
    weights_Iref_map_T0('ODE.I')= 0.5;
    weights_Iref_map_T0('SPT.I')= 0.5;
    CPDFactory_Iref = ...
        CPDFactory('Gaussian_CPD', 'Reference.I', 0, ...
        {'mean', 0.0, 'cov', 5.0}, ...
        weights_Iref_map_T0, ...
        weights_Iref_map_T1);
    add_CPD_factories(meta_dbn_factory, {CPDFactory_Iref}, false);

    [meta_dbn, ~, ~, nodes_map] = create_dbn(meta_dbn_factory);
    
end
