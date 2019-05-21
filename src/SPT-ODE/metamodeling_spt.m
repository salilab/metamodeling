% Make a DBN for the meal model with the following variables
%
% Time-dependent variables
% 
% Reference variables
%
% Observed variables
%
% Time-invariant variables
%
% Parameters
%
% The ODE equation is only when G > h
% TODO: check the time step of the equation, and see the G and I.
warning('off','MATLAB:singularMatrix');

%Read in the experimental measurements
clear;

[dbn, nodes_map] = make_meta_bnet(5.111, 0, 0.05, 39.6,828, 3, 0.056, 5.111, 34, ...
    0.1, 34, 10, 300, 6, 1.8e-6, ...
    0.1, 63, 34, ...
    0.05E-3, 0.008);

[meal_dbn_factory]= make_meal_dbn_factory_eq2(5.111, 0, 0.05, 39.6,828, 3, 0.056, 5.111, 34);
[spt_dbn_factory]= make_spt_dbn_factory(0.1, 34, 10, 300, 6, 1.8e-6);
[network_dbn_factory]= make_network_dbn_factory(0.1, 63, 34);
[glp1r_dbn_factory]= make_glp1r_dbn_factory(0.05E-3, 0.008);

% Sample
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);

[meal_dbn, ~, ~, meal_nodes_map] = create_dbn(meal_dbn_factory);
meal_dbn_engine = jtree_dbn_inf_engine(meal_dbn);

[spt_dbn, ~, ~, spt_nodes_map] = create_dbn(spt_dbn_factory);
spt_dbn_engine = jtree_dbn_inf_engine(spt_dbn);

[network_dbn, ~, ~, network_nodes_map] = create_dbn(network_dbn_factory);
network_dbn_engine = jtree_dbn_inf_engine(network_dbn);

[glp1r_dbn, ~, ~, glp1r_nodes_map] = create_dbn(glp1r_dbn_factory);
glp1r_dbn_engine = jtree_dbn_inf_engine(glp1r_dbn);

meal_npers= meal_dbn.nnodes_per_slice;
spt_npers= spt_dbn.nnodes_per_slice;
network_npers= network_dbn.nnodes_per_slice;
glp1r_npers= glp1r_dbn.nnodes_per_slice;

T = 100;

sample_seq=  cell2mat(sample_dbn(dbn, 'length', T));
nodes_order= cell2mat(values(nodes_map));

% Plot ODE.G and ODE.I
figure()
yyaxis left;
Gvalues = sample_seq(nodes_map('cAMP.Network'),:);
plot(1:T, Gvalues);
yyaxis right;
Ivalues = sample_seq(nodes_map('G.obs'),:);
plot(1:T, Ivalues);
legend('cAMP.Network','S.Network');  
legend('boxoff');

%Posterior marginal of h(20) given Iexp

% Plot ODE.G and ODE.I
%sample_seq=  cell2mat(sample_dbn(dbn, 'length', T));
%nodes_order= cell2mat(values(nodes_map));
%figure()
%yyaxis left;
%plot(1:T, sample_seq(nodes_map('G.Meal'),:));
%yyaxis right;
%plot(1:T, sample_seq(nodes_map('cAMP.Network'),:));
%legend('G.Meal','S.Meal');  
T = 10;
spt_evidence= cell(spt_npers, T);
spt_evidence{spt_nodes_map('S.obs'),1} = 30.0; 
i=2;
disp(i);
disp(spt_npers);
spt_marg= marginal_nodes(enter_evidence(spt_dbn_engine, spt_evidence), ...
                     spt_nodes_map('k.SPT')+spt_npers, ...
                     i);
fprintf("%f sigma %f +- %f", spt_marg.mu, spt_marg.Sigma, sqrt(spt_marg.Sigma)) % mean +- stddev
 
evidence= cell(npers, T);
%evidence{nodes_map('S.obs'),1} = 30.0; 
%evidence{nodes_map('cAMP.Network'),1} = 1.3E-3;
marg= marginal_nodes(enter_evidence(dbn_engine, evidence), ...
                     nodes_map('k.SPT')+npers, ...
                     i);
fprintf("%f sigma %f +- %f", marg.mu, marg.Sigma, sqrt(marg.Sigma)) % mean +- stddev

% plot k
xx = linspace(7,13 );

% create some plot with a legend
hAx(1) = axes();
hLine(1) = plot(xx,normpdf(xx,spt_marg.mu,spt_marg.Sigma),'k-','LineWidth',2, 'Parent',hAx(1));
set(hAx(1), 'Box','off');
legend(hLine(1), {'SPT.k, posterior'},'Location','NorthWest');

% copy the axis
hAx(2) = copyobj(hAx(1),gcf);
delete(get(hAx(2),'Children'));            %# delete its children
hLine(2) = plot(xx,normpdf(xx,marg.mu,marg.Sigma),'r-','LineWidth',2,'Parent',hAx(2));
set(hAx(2), 'Color','none', 'XTick',[], ...
    'YAxisLocation','right', 'Box','off');   %# make it transparent
legend(hLine(2), {'Meta.k, posterior'}, 'Color','w');

% Create a table with the data and variable names
yy=normpdf(xx,marg.mu,marg.Sigma)/sum(normpdf(xx,marg.mu,marg.Sigma))
variable = [xx(:) yy(:)];
size(variable);
dlmwrite('I_meta.txt',variable);

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
writetable(T, 'I.txt')
