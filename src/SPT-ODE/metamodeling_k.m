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

m1 = importdata('insulin_exp1.dat');
I1 = m1(:,2); % Gexp in measurement number 1, vector along time

% TODO: To implement model along time
time = 1;
[dbn, nodes_map] = make_meta_bnet(Gm1, Im1, Go1, Io1, time);
[ode_dbn_factory]= make_ode_dbn_factory(Gm1, Im1, time);
[spt_dbn_factory]= make_spt_dbn_factory(Go1, Io1, time);
[ode_dbn, ~, ~, ode_nodes_map] = create_dbn(ode_dbn_factory);
[spt_dbn, ~, ~, spt_nodes_map] = create_dbn(spt_dbn_factory);

% Sample
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);
ode_dbn_engine = jtree_dbn_inf_engine(ode_dbn);
spt_dbn_engine = jtree_dbn_inf_engine(spt_dbn);

ode_npers= ode_dbn.nnodes_per_slice;
spt_npers= spt_dbn.nnodes_per_slice;
T = 400;

%Posterior marginal of h(20) given Iexp

spt_evidence= cell(spt_npers, T);
spt_evidence{spt_nodes_map('E.Ipm'),2} = 55.0; 
i=10;
disp(i);
disp(spt_npers);
spt_marg= marginal_nodes(enter_evidence(spt_dbn_engine, spt_evidence), ...
                     spt_nodes_map('SPT.k')+spt_npers, ...
                     i);
fprintf("%f sigma %f +- %f", spt_marg.mu, spt_marg.Sigma, sqrt(spt_marg.Sigma)) % mean +- stddev
 
evidence= cell(npers, T);
evidence{nodes_map('E.Ipm'),2} = 25.0; 
marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
                     nodes_map('SPT.k')+npers, ...
                     i);
fprintf("%f sigma %f +- %f", marg.mu, marg.Sigma, sqrt(marg.Sigma)) % mean +- stddev

% plot k
xx = linspace(8,11);

% create some plot with a legend
hAx(1) = axes();
hLine(1) = plot(xx,normpdf(xx,spt_marg.mu,spt_marg.Sigma),'k-','LineWidth',2, 'Parent',hAx(1));
set(hAx(1), 'Box','off');
legend(hLine(1), {'SPT.k, '},'Location','NorthWest');

% copy the axis
hAx(2) = copyobj(hAx(1),gcf);
delete(get(hAx(2),'Children'));            %# delete its children
hLine(2) = plot(xx,normpdf(xx,marg.mu,marg.Sigma),'r-','LineWidth',2,'Parent',hAx(2));
set(hAx(2), 'Color','none', 'XTick',[], ...
    'YAxisLocation','right', 'Box','off');   %# make it transparent
legend(hLine(2), {'SPT.k, posterior'}, 'Color','w');

% test different measurement
sptk = {};
metak= {};
for measure = 1:50  
    
spt_evidence= cell(spt_npers, T);
spt_evidence{spt_nodes_map('E.Ipm'),2} = I1(measure); 
i=10;
disp(i);
disp(spt_npers);
spt_marg= marginal_nodes(enter_evidence(spt_dbn_engine, spt_evidence), ...
                     spt_nodes_map('SPT.k')+spt_npers, ...
                     i);
fprintf("%f sigma %f +- %f", spt_marg.mu, spt_marg.Sigma, sqrt(marg.Sigma)) % mean +- stddev
sptk(end+1,:) = {spt_marg.mu, sqrt(spt_marg.Sigma)};
evidence= cell(npers, T);
evidence{nodes_map('E.Ipm'),2} = I1(measure); 
marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
                     nodes_map('SPT.k')+npers, ...
                     i);
fprintf("%f sigma %f +- %f", marg.mu, marg.Sigma, sqrt(marg.Sigma)) % mean +- stddev
metak(end+1,:) = {marg.mu, sqrt(marg.Sigma)};
end

% Create a table with the data and variable names
T = table(sptk, metak )
% Write data to text file
writetable(T, 'k.txt')
