% Sample the meta-modeling DBN containing following five input models:
% 1. The meal model of T2D subject
% 2. The SPT model
% 3. The network model
% 4. The GLP1R model
% 5. The KEGG model
%
% Sampling parameter:
% Evidence: 
%   1. The DGintake of the meal model is set to the experimentally obtained
%      profile after a meal for 420 time slices.
% 
% Return:
% Model variables in the meal model, the SPT model and the
% netwok model at different time slices before and after meta-modeling.
%
warning('off','MATLAB:singularMatrix');

% Read in the experimental measurements
Experiment3 = importdata('/Users/lipingsun/SynologyDrive/research/Meta-modeling/bnt/Metamodeling-bnet/dataset/meal/140points/meal_exp1_normal2_Gin_dt1_sigmoid.dat');
DGexp = Experiment3(:,2); % Rate of change in glucose intake, mM
Gintakeexp = Experiment3(:,3); % Glucose intake, mM

%----------------------------------------
% Initial input values of model variables
%----------------------------------------
% Meal model normal
Gex_Meal = 9.167; % Basal plasma glucose level, mM
Y_Meal = 0; % provision of new insulin to the beta-cells at t=0, pM/min
alpha_Meal = 0.013; % Delay between the glucose signal and insulin secretion, min^{-1}
beta_Meal = 22.5; % Pancreatic responsivity to glucose, pM/min per mM
gamma_Meal = 0.5; % transfer rate constant between portal vein and liver, min^{-1}
k1_Meal = 0.0001; % Coefficient for I reducing glucose level, significantly affect the periodicity of the curves.
k2_Meal = 0.00001; % Coefficient for Gex reducing glucose level
K_Meal = 10; % Pancreatic responsivity to the glucose rate of change, pmol/L per mM
dt_Meal_min = 1; % time scale, min
Gb_Meal = 9.167; % Basal plasma glucose level, mM
DGintake_Meal = DGexp(1); % Rate of glucose intake from food at t = 0, mM/min
Sb_Meal = 102.5; %Basal insulin secretion, pM/min
I_Meal = 52; % Basal plasma insulin level, pM
% SPT
Gin_SPT =  5.111/2 % Basal intracellular glucose concentration at the first time slice, mM
S_SPT =  34 % Basal insulin secretion at the first time slice, normalized for pancreas, pM/min 
k_SPT =  10 % Force coefficient, m/s
Nisg_SPT = 300 %  Number of ISGs
Npatch_SPT = 6 % Number of activation patches per ISG
Ninsulin_SPT = 1.8e-6 % Amount of insulin molecules in one ISG, pmol
I_SPT = 25; % Basal plasma insulin level at the first time slice, normalized for pancreas, pM
% Network
Gin_Network =  5.111/2 % Basal intracellular glucose concentration at the first time slice, mM
PFK_Network = 0.65 % Enzyme activily, no unit
ATP_Network =  3.3 % intracellular ATP concentration at the first time slice, mM
GLP1_Network = 12.5 % intracellular GLP1 concentration at the first time slice, pM 
GLP1R_Network = 1.0 % GLP1R activity at the first time slice, no unit
cAMP_Network =  1.3E-3 % intracellular cAMP concentration at the first time slice, mM
Ca_Network =  1E-4 % intracellular Ca concentration at the first time slice, mM
S_Network = 34 % Basal insulin secretion at the first time slice, normalized for pancreas, pM/min 
I_Network = 25; % Basal plasma insulin level at the first time slice, normalized for pancreas, pM
% GLP1R
GLP1a_GLP1R =  0.0 % Type of the GLP1 analogues, no unit - analogue 1-10 with increasing affinity to GLP1R
cons_GLP1R = 0.0 % Compound concentration of the GLP1 analogues , nM - 1-10
GLP1R_GLP1R =  0.0 % GLP1R activity at the first time slice, no unit
% KEGG
Gin_KEGG = 0 % Extracellular glucose concentration, mM 
Ex4_KEGG = 0 % GLP-1 treatmeent, nM
NES_KEGG = 1 % Normalized enrichment score, no unit
Pathway_KEGG = 1 % Pathway number, no unit
ATP_KEGG = 0 % intracellular ATP concentration, mM

%----------------------------------------
% Make bnet and DBN
%----------------------------------------
[dbn, nodes_map] = make_meta_bnet5(Gex_Meal, Y_Meal, alpha_Meal,beta_Meal, gamma_Meal, k1_Meal, k2_Meal, K_Meal, dt_Meal_min, Gb_Meal, DGintake_Meal, Sb_Meal, I_Meal, ...
    Gin_SPT, S_SPT, k_SPT, Nisg_SPT, Npatch_SPT, Ninsulin_SPT, I_SPT, ...
    Gin_Network, PFK_Network, ATP_Network, GLP1_Network, GLP1R_Network, cAMP_Network, Ca_Network, S_Network, I_Network, ...
    GLP1a_GLP1R, cons_GLP1R,GLP1R_GLP1R,...
    Gin_KEGG, Ex4_KEGG, NES_KEGG, Pathway_KEGG, ATP_KEGG);

[meal_dbn_factory]= make_meal_dbn_factory_eq(Gex_Meal, Y_Meal, alpha_Meal,beta_Meal, gamma_Meal, k1_Meal, k2_Meal, K_Meal, dt_Meal_min, Gb_Meal, DGintake_Meal, Sb_Meal, I_Meal);
[spt_dbn_factory]= make_spt_dbn_factory(Gin_SPT, S_SPT, k_SPT, Nisg_SPT, Npatch_SPT, Ninsulin_SPT, I_SPT);
[network_dbn_factory]= make_network_dbn_factory(Gin_Network, PFK_Network, ATP_Network, GLP1_Network, GLP1R_Network, cAMP_Network, Ca_Network, S_Network, I_Network);
[glp1r_dbn_factory]= make_glp1r_dbn_factory(GLP1a_GLP1R, cons_GLP1R, GLP1R_GLP1R);
[kegg_dbn_factory]= make_KEGG_dbn_factory(Gculture_KEGG, Ex4_KEGG, NES_KEGG, Pathway_KEGG, ATP_KEGG);

% Sample
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);

[meal_dbn, ~, ~, meal_nodes_map] = create_dbn(meal_dbn_factory);
%meal_dbn_engine = jtree_dbn_inf_engine(meal_dbn);

[spt_dbn, ~, ~, spt_nodes_map] = create_dbn(spt_dbn_factory);
%spt_dbn_engine = jtree_dbn_inf_engine(spt_dbn);

[network_dbn, ~, ~, network_nodes_map] = create_dbn(network_dbn_factory);
%network_dbn_engine = jtree_dbn_inf_engine(network_dbn);

[glp1r_dbn, ~, ~, glp1r_nodes_map] = create_dbn(glp1r_dbn_factory);
%glp1r_dbn_engine = jtree_dbn_inf_engine(glp1r_dbn);

[kegg_dbn, ~, ~, kegg_nodes_map] = create_dbn(kegg_dbn_factory);
%glp1r_dbn_engine = jtree_dbn_inf_engine(glp1r_dbn);

meal_npers= meal_dbn.nnodes_per_slice;
spt_npers= spt_dbn.nnodes_per_slice;
network_npers= network_dbn.nnodes_per_slice;
glp1r_npers= glp1r_dbn.nnodes_per_slice;
kegg_npers= kegg_dbn.nnodes_per_slice;

%----------------------------------------
% Sampling
%----------------------------------------
T = 420;

sample_seq=  cell2mat(sample_dbn(dbn, 'length', T));
nodes_order= cell2mat(values(nodes_map));

%compute the unconditional marginals of h(20)
meal_dbn_engine = jtree_dbn_inf_engine(meal_dbn);
spt_dbn_engine = jtree_dbn_inf_engine(spt_dbn);
network_dbn_engine = jtree_dbn_inf_engine(network_dbn);
glp1r_dbn_engine = jtree_dbn_inf_engine(glp1r_dbn);
kegg_dbn_engine = jtree_dbn_inf_engine(kegg_dbn);

%meta model
T = 420;
meal_Gex={};
meal_S={};
meal_Gb={};
meal_I={};
spt_Gin={};
spt_k={};
spt_Npatch={};
spt_Nisg={};
spt_Ninsulin={};
spt_S={};
spt_I={};
network_Gin={};
network_cAMP={};
network_Ca={};
network_S={};
network_I={};

Gex_ref={};
Gin_ref={};
Scell_ref={};
GLP1R_ref={};
meta_meal_Gex={};
meta_meal_S={};
meta_meal_Gb={};
meta_meal_I={};
meta_spt_Gin={};
meta_spt_k={};
meta_spt_Npatch={};
meta_spt_Nisg={};
meta_spt_Ninsulin={};
meta_spt_S={};
meta_spt_I={};
meta_network_Gin={};
meta_network_cAMP={};
meta_network_Ca={};
meta_network_S={};
meta_network_I={};

% Enter the evidence
meal_evidence= cell(meal_npers, T);
spt_evidence= cell(spt_npers, T);
network_evidence= cell(network_npers, T);
glp1r_evidence= cell(glp1r_npers, T);
kegg_evidence= cell(kegg_npers, T);
evidence= cell(npers, T);

% Sampling along the time slices
for measure = 1:420
    disp(measure);
    evidence{nodes_map('DGintake.obs'),measure} = DGexp(measure);
    evidence{nodes_map('GLP1.obs'),measure} = GLP1_Network; 
    evidence{nodes_map('GLP1a.obs'),measure} = 0; 
    evidence{nodes_map('cons.obs'),measure} = 0; 
    evidence{nodes_map('Gculture.obs'),measure} = 0; 
    evidence{nodes_map('Ex4.obs'),measure} = 0; 
    
    meal_evidence{meal_nodes_map('DGintake.obs'),measure} = DGexp(measure); 
    network_evidence{network_nodes_map('GLP1.obs'),measure} = GLP1_Network;
    glp1r_evidence{glp1r_nodes_map('GLP1a.obs'),measure} = 0;
    glp1r_evidence{glp1r_nodes_map('cons.obs'),measure} = 0;
    kegg_evidence{kegg_nodes_map('Gculture.obs'),measure} = 0;
    kegg_evidence{kegg_nodes_map('Ex4.obs'),measure} = 0; 
end

[meal_engine, ll] = enter_evidence(meal_dbn_engine, meal_evidence);
[spt_engine, ll] = enter_evidence(spt_dbn_engine, spt_evidence);
[network_engine, ll] = enter_evidence(network_dbn_engine, network_evidence);
[glp1r_engine, ll] = enter_evidence(glp1r_dbn_engine, glp1r_evidence);
[kegg_engine, ll] = enter_evidence(kegg_dbn_engine, kegg_evidence);
[engine, ll] = enter_evidence(dbn_engine, evidence);

% Sampling along the time slices
for measure = 1:420
    marg_meal_Gex= marginal_nodes(meal_engine,meal_nodes_map('Gex.Meal'),measure);
    marg_meal_S= marginal_nodes(meal_engine,meal_nodes_map('S.Meal'),measure);
    marg_meal_Gb= marginal_nodes(meal_engine,meal_nodes_map('Gb.Meal'),measure);
    marg_meal_I= marginal_nodes(meal_engine,meal_nodes_map('I.Meal'),measure);
    
    marg_spt_Gin= marginal_nodes(spt_engine,spt_nodes_map('Gin.SPT'),measure);
    marg_spt_k= marginal_nodes(spt_engine,spt_nodes_map('k.SPT'),measure);
    marg_spt_Npatch= marginal_nodes(spt_engine,spt_nodes_map('Npatch.SPT'),measure);
    marg_spt_Nisg= marginal_nodes(spt_engine,spt_nodes_map('Nisg.SPT'),measure);
    marg_spt_Ninsulin= marginal_nodes(spt_engine,spt_nodes_map('Ninsulin.SPT'),measure);
    marg_spt_S= marginal_nodes(spt_engine,spt_nodes_map('S.SPT'),measure);
    marg_spt_I= marginal_nodes(spt_engine,spt_nodes_map('I.SPT'),measure);
    
    marg_network_Gin= marginal_nodes(network_engine,network_nodes_map('Gin.Network'),measure);
    marg_network_cAMP= marginal_nodes(network_engine,network_nodes_map('cAMP.Network'),measure);
    marg_network_Ca= marginal_nodes(network_engine,network_nodes_map('Ca.Network'),measure);
    marg_network_S= marginal_nodes(network_engine,network_nodes_map('S.Network'),measure);
    marg_network_I= marginal_nodes(network_engine,network_nodes_map('I.Network'),measure);
    
    marg_Gex_ref= marginal_nodes(engine,nodes_map('Gex.ref'),measure);
    marg_Gin_ref= marginal_nodes(engine,nodes_map('Gin.ref'),measure);
    marg_Scell_ref= marginal_nodes(engine,nodes_map('Scell.ref'),measure);
    marg_GLP1R_ref= marginal_nodes(engine,nodes_map('GLP1R.ref'),measure);
    marg_meta_meal_Gex= marginal_nodes(engine,nodes_map('Gex.Meal'),measure);
    marg_meta_meal_S= marginal_nodes(engine,nodes_map('S.Meal'),measure);
    marg_meta_meal_Gb= marginal_nodes(engine,nodes_map('Gb.Meal'),measure);
    marg_meta_meal_I= marginal_nodes(engine,nodes_map('I.Meal'),measure);
    
    marg_meta_spt_Gin= marginal_nodes(engine,nodes_map('Gin.SPT'),measure);
    marg_meta_spt_k= marginal_nodes(engine,nodes_map('k.SPT'),measure);
    marg_meta_spt_Npatch= marginal_nodes(engine,nodes_map('Npatch.SPT'),measure);
    marg_meta_spt_Nisg= marginal_nodes(engine,nodes_map('Nisg.SPT'),measure);
    marg_meta_spt_Ninsulin= marginal_nodes(engine,nodes_map('Ninsulin.SPT'),measure);
    marg_meta_spt_S= marginal_nodes(engine,nodes_map('S.SPT'),measure);
    marg_meta_spt_I= marginal_nodes(engine,nodes_map('I.SPT'),measure);
    
    marg_meta_network_Gin= marginal_nodes(engine,nodes_map('Gin.Network'),measure);
    marg_meta_network_cAMP= marginal_nodes(engine,nodes_map('cAMP.Network'),measure);
    marg_meta_network_Ca= marginal_nodes(engine,nodes_map('Ca.Network'),measure);
    marg_meta_network_S= marginal_nodes(engine,nodes_map('S.Network'),measure);
    marg_meta_network_I= marginal_nodes(engine,nodes_map('I.Network'),measure);

    meal_Gex(end+1,:) = {marg_meal_Gex.mu, marg_meal_Gex.Sigma, sqrt(marg_meal_Gex.Sigma)};
    meal_S(end+1,:) = {marg_meal_S.mu, marg_meal_S.Sigma, sqrt(marg_meal_S.Sigma)};
    meal_Gb(end+1,:) = {marg_meal_Gb.mu, marg_meal_Gb.Sigma, sqrt(marg_meal_Gb.Sigma)};
    meal_I(end+1,:) = {marg_meal_I.mu, marg_meal_I.Sigma, sqrt(marg_meal_I.Sigma)};
    
    spt_Gin(end+1,:) = {marg_spt_Gin.mu, marg_spt_Gin.Sigma, sqrt(marg_spt_Gin.Sigma)};
    spt_k(end+1,:) = {marg_spt_k.mu, marg_spt_k.Sigma, sqrt(marg_spt_k.Sigma)};
    spt_Npatch(end+1,:) = {marg_spt_Npatch.mu, marg_spt_Npatch.Sigma, sqrt(marg_spt_Npatch.Sigma)};
    spt_Ninsulin(end+1,:) = {marg_spt_Ninsulin.mu, marg_spt_Ninsulin.Sigma, sqrt(marg_spt_Ninsulin.Sigma)};
    spt_Nisg(end+1,:) = {marg_spt_Nisg.mu, marg_spt_Nisg.Sigma, sqrt(marg_spt_Nisg.Sigma)};
    spt_S(end+1,:) = {marg_spt_S.mu, marg_spt_S.Sigma, sqrt(marg_spt_S.Sigma)};
    spt_I(end+1,:) = {marg_spt_I.mu, marg_spt_I.Sigma, sqrt(marg_spt_I.Sigma)};
   
    network_Gin(end+1,:) = {marg_network_Gin.mu, marg_network_Gin.Sigma, sqrt(marg_network_Gin.Sigma)};
    network_cAMP(end+1,:) = {marg_network_cAMP.mu, marg_network_cAMP.Sigma, sqrt(marg_network_cAMP.Sigma)};
    network_Ca(end+1,:) = {marg_network_Ca.mu, marg_network_Ca.Sigma, sqrt(marg_network_Ca.Sigma)};
    network_S(end+1,:) = {marg_network_S.mu, marg_network_S.Sigma, sqrt(marg_network_S.Sigma)};
    network_I(end+1,:) = {marg_network_I.mu, marg_network_I.Sigma, sqrt(marg_network_I.Sigma)};
    
    Gex_ref(end+1,:) = {marg_Gex_ref.mu, marg_Gex_ref.Sigma, sqrt(marg_Gex_ref.Sigma)};
    Gin_ref(end+1,:) = {marg_Gin_ref.mu, marg_Gin_ref.Sigma, sqrt(marg_Gin_ref.Sigma)};
    Scell_ref(end+1,:) = {marg_Scell_ref.mu, marg_Scell_ref.Sigma, sqrt(marg_Scell_ref.Sigma)};
    GLP1R_ref(end+1,:) = {marg_GLP1R_ref.mu, marg_GLP1R_ref.Sigma, sqrt(marg_GLP1R_ref.Sigma)};
    meta_meal_Gex(end+1,:) = {marg_meta_meal_Gex.mu, marg_meta_meal_Gex.Sigma, sqrt(marg_meta_meal_Gex.Sigma)};
    meta_meal_S(end+1,:) = {marg_meta_meal_S.mu, marg_meta_meal_S.Sigma, sqrt(marg_meta_meal_S.Sigma)};
    meta_meal_Gb(end+1,:) = {marg_meta_meal_Gb.mu, marg_meta_meal_Gb.Sigma, sqrt(marg_meta_meal_Gb.Sigma)};
    meta_meal_I(end+1,:) = {marg_meta_meal_I.mu, marg_meta_meal_I.Sigma, sqrt(marg_meta_meal_I.Sigma)};
    meta_spt_Gin(end+1,:) = {marg_meta_spt_Gin.mu, marg_meta_spt_Gin.Sigma, sqrt(marg_meta_spt_Gin.Sigma)};
    meta_spt_k(end+1,:) = {marg_meta_spt_k.mu, marg_meta_spt_k.Sigma, sqrt(marg_meta_spt_k.Sigma)};
    meta_spt_Npatch(end+1,:) = {marg_meta_spt_Npatch.mu, marg_meta_spt_Npatch.Sigma, sqrt(marg_meta_spt_Npatch.Sigma)};
    meta_spt_Nisg(end+1,:) = {marg_meta_spt_Nisg.mu, marg_meta_spt_Nisg.Sigma, sqrt(marg_meta_spt_Nisg.Sigma)};
    meta_spt_Ninsulin(end+1,:) = {marg_meta_spt_Ninsulin.mu, marg_meta_spt_Ninsulin.Sigma, sqrt(marg_meta_spt_Ninsulin.Sigma)};
    meta_spt_S(end+1,:) = {marg_meta_spt_S.mu, marg_meta_spt_S.Sigma, sqrt(marg_meta_spt_S.Sigma)};
    meta_spt_I(end+1,:) = {marg_meta_spt_I.mu, marg_meta_spt_I.Sigma, sqrt(marg_meta_spt_I.Sigma)};
    meta_network_Gin(end+1,:) = {marg_meta_network_Gin.mu, marg_meta_network_Gin.Sigma, sqrt(marg_meta_network_Gin.Sigma)};
    meta_network_cAMP(end+1,:) = {marg_meta_network_cAMP.mu, marg_meta_network_cAMP.Sigma, sqrt(marg_meta_network_cAMP.Sigma)};
    meta_network_Ca(end+1,:) = {marg_meta_network_Ca.mu, marg_meta_network_Ca.Sigma, sqrt(marg_meta_network_Ca.Sigma)};
    meta_network_S(end+1,:) = {marg_meta_network_S.mu, marg_meta_network_S.Sigma, sqrt(marg_meta_network_S.Sigma)};
    meta_network_I(end+1,:) = {marg_meta_network_I.mu, marg_meta_network_I.Sigma, sqrt(marg_meta_network_I.Sigma)};
    
end
disp("separate");
disp(meal_I);
disp("separate");
disp(meal_Gb);
disp("separate");
disp(spt_I);
disp("separate");
disp(network_Ca);
disp("separate");
disp(network_I);
disp("separate");
disp(meta_meal_I);

%----------------------------------------
% Write out the sampled model variables 
%----------------------------------------
% Create a table with the data and variable names
T = table(meal_Gex, meal_S, meal_Gb, meal_I, meta_meal_Gex, meta_meal_S,meta_meal_Gb, meta_meal_I, Scell_ref, GLP1R_ref, Gex_ref, Gin_ref);
% Write data to text file
writetable(T, 'meta_meal_t2d.txt');

% Create a table with the data and variable names
T1 = table(spt_Gin, spt_k, spt_Npatch, spt_Nisg, spt_Ninsulin, spt_S, spt_I,meta_spt_Gin, meta_spt_k, meta_spt_Npatch, meta_spt_Nisg, meta_spt_Ninsulin, meta_spt_S, meta_spt_I);
% Write data to text file
writetable(T1, 'meta_spt_t2d.txt');

% Create a table with the data and variable names
T2 = table(network_Gin, network_cAMP, network_Ca, network_S, network_I, meta_network_Gin, meta_network_cAMP, meta_network_Ca, meta_network_S, meta_network_I );
% Write data to text file
writetable(T2, 'meta_network_t2d.txt');
