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

% check the time step of the equation, and see the G and I.

% Read in the experimental measurements
SPTm1 = importdata('SPT_obs1_avr.dat');
Gm1 = SPTm1(:,2); % Gexp in measurement number 1, vector along time
Im1 = SPTm1(:,3); % Iexp in measurement number 1, vector along time
disp(Gm1(1));

%for time  = 1:length(Gm1)
%    [bnet, nodes_map, intra, inter]= make_ode_bnet(Gm1, Im1,time);
%end
n = 9;

time = 3
[dbn_factory]= make_spt_dbn_factory(Gm1, Im1, time);
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
xx = linspace(5,15);
hold on; 
plot(xx,area*normpdf(xx,mean(y),std(y)),'k-','LineWidth',2);
fprintf("Normal probability density function of SPT.k");
disp(mean(y));
disp(std(y));
%f = ksdensity(y,xx);
%plot(xx,area*f,'g-')
legend('SPT.k, prior');
hold off;
