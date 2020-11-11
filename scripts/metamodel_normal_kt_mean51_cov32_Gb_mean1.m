% Sample the meta-modeling DBN containing following five input models:
% 1. The postprandial model
% 2. The exocytosis model
% 3. The signaling model
% 4. The GLP1R model
% 5. The metabolism model

SCRIPTS_PATH=pwd(); % Change this line if not running from scripts path
OUTPUT_PATH= SCRIPTS_PATH+"/Output/accuracy_precision/Gb_kt/";
if ~exist(OUTPUT_PATH, 'dir')
    mkdir(OUTPUT_PATH)
end
% Add bnet
cd ../bnt_master
addpath(genpathKPM('../bnt_master'))
cd(SCRIPTS_PATH) %../bnet_scripts_Gb-kt

warning('off','MATLAB:singularMatrix');

% ---------------------------------
% Read data as input and evidence
% ---------------------------------
% postprandial model
Json_postprandial = jsondecode(fileread('../data/postprandial_normal.json'));
DGexp = importdata(Json_postprandial.EvidenceDG);
EvidenceDG = DGexp(:,2);
% pancreas model
Json_pancreas = jsondecode(fileread(['../data/pancreas.json']));
% exocytosis model
Json_exocytosis = jsondecode(fileread(['../data/exocytosis.json']));
% signaling model
Json_signaling = jsondecode(fileread('../data/signaling.json'));
% screening model
Json_screening = jsondecode(fileread('../data/screening.json'));
% metabolism model
Json_metabolism = jsondecode(fileread('../data/metabolism.json'));

% meta model
Json_meta = jsondecode(fileread('../data/meta_normal.json'));
Gb_k_input = importdata(Json_meta.InputErr);

% meta model
Json_meta = jsondecode(fileread('../data/meta_normal.json'));
Gb_kt_mu_input = importdata(Json_meta.InputErr); % an nx2 array with possible values of Gb and kt to scan
Gb_kt_sigma_input = importdata(Json_meta.InputSigma);
  
% Create a table with the data and variable names
inputGb = {};
inputGb(end+1,:) = {'Gb_mean_error' 'Gb_cov'};
inputkt = {};
inputkt(end+1,:) = {'kt_mean_error' 'kt_cov'};
outputGb = {};
outputGb(end+1,:) = {'Gb.postprandial' 'Gb.postprandial' 'Gb.postprandial'};
outputkt = {};
outputkt(end+1,:) = {'kt.exocytosis' 'kt.exocytosis' 'kt.exocytosis'};

outputGex = {};
outputGex(end+1,:) = {'G.postprandial' 'G.postprandial' 'G.postprandial'};
outputGin = {};
outputGin(end+1,:) = {'G.exocytosis' 'G.exocytosis' 'G.exocytosis'};
outputSp = {};
outputSp(end+1,:) = {'S.postprandial' 'S.postprandial' 'S.postprandial'};
outputSpa = {};
outputSpa(end+1,:) = {'Spa.pancreas' 'Spa.pancreas' 'Spa.pancreas'};
outputSis = {};
outputSis(end+1,:) = {'Sis.pancreas' 'Sis.pancreas' 'Sis.pancreas'};
outputScell = {};
outputScell(end+1,:) = {'Scell.pancreas' 'Scell.pancreas' 'Scell.pancreas'};
outputSe = {};
outputSe(end+1,:) = {'S.exocytosis' 'S.exocytosis' 'S.exocytosis'};

disp(length(Gb_kt_mu_input(:,1)));

kt_mean_min_ndx = 51; % minimal 
kt_mean_max_ndx = 51;
kt_cov_min_ndx = 32;
kt_cov_max_ndx = 32;

Gb_mean_min_ndx = 1;
Gb_mean_max_ndx = 1;
Gb_cov_min_ndx = 1;
Gb_cov_max_ndx = 101;

for i = kt_mean_min_ndx:  kt_mean_max_ndx
    kt_mean_error = Gb_kt_mu_input(i,2)
    for j = kt_cov_min_ndx:  kt_cov_max_ndx
        kt_cov = Gb_kt_sigma_input(j,2)
        for m = Gb_mean_min_ndx: Gb_mean_max_ndx
            Gb_mean_error = Gb_kt_mu_input(m,1)
            for n =Gb_cov_min_ndx: Gb_cov_max_ndx
                Gb_cov = Gb_kt_sigma_input(n,1)
                inputGb(end+1,:) = {Gb_mean_error Gb_cov};
                inputkt(end+1,:) = {kt_mean_error kt_cov};
                [dbn, nodes_map]= make_meta_dbn6(Json_postprandial.DataInput.DGd_mean_postprandial, Json_postprandial.DataInput.DGd_cov_postprandial, Gb_mean_error,...
                                 Gb_cov, Json_postprandial.DataInput.G_mean_postprandial, Json_postprandial.DataInput.G_cov_postprandial,...
                                 Json_postprandial.DataInput.DG_mean_postprandial, Json_postprandial.DataInput.DG_cov_postprandial, Json_postprandial.DataInput.Y_mean_postprandial,...
                                 Json_postprandial.DataInput.Y_cov_postprandial, Json_postprandial.DataInput.S_mean_postprandial, Json_postprandial.DataInput.S_cov_postprandial,...
                                 Json_postprandial.DataInput.I_mean_postprandial, Json_postprandial.DataInput.I_cov_postprandial, Json_postprandial.DataInput.Sb_mean_postprandial,...
                                 Json_postprandial.DataInput.Sb_cov_postprandial, Json_postprandial.DataInput.alpha_postprandial, Json_postprandial.DataInput.beta_postprandial,...
                                 Json_postprandial.DataInput.gamma_postprandial, Json_postprandial.DataInput.k1_postprandial,Json_postprandial.DataInput.k2_postprandial,...
                                 Json_postprandial.DataInput.k3_postprandial, Json_postprandial.DataInput.k4_postprandial, Json_postprandial.DataInput.K_postprandial,...
                                 Json_postprandial.DataInput.dt_postprandial, Json_postprandial.DataInput.cov_scale_postprandial,...
                                 Json_pancreas.DataInput.Scell_mean_pancreas, Json_pancreas.DataInput.Scell_cov_pancreas, Json_pancreas.DataInput.Sislet_mean_pancreas, ...
                                 Json_pancreas.DataInput.Sislet_cov_pancreas, Json_pancreas.DataInput.Spancreas_mean_pancreas, Json_pancreas.DataInput.Spancreas_cov_pancreas, ...
                                 Json_pancreas.DataInput.cov_scale_pancreas, Json_pancreas.DataInput.Nc_pancreas, Json_pancreas.DataInput.Ni_pancreas,...
                                 Json_exocytosis.DataInput.G_mean_exocytosis, kt_mean_error, Json_exocytosis.DataInput.Npatch_mean_exocytosis,...
                                 Json_exocytosis.DataInput.Nvesicle_mean_exocytosis, Json_exocytosis.DataInput.Ninsulin_mean_exocytosis, Json_exocytosis.DataInput.Rcell_mean_exocytosis,...
                                 Json_exocytosis.DataInput.Dvesicle_mean_exocytosis, Json_exocytosis.DataInput.S_mean_exocytosis, Json_exocytosis.DataInput.G_cov_exocytosis,...
                                 kt_cov, Json_exocytosis.DataInput.Npatch_cov_exocytosis, Json_exocytosis.DataInput.Nvesicle_cov_exocytosis, ...
                                 Json_exocytosis.DataInput.Ninsulin_cov_exocytosis, Json_exocytosis.DataInput.Rcell_cov_exocytosis,Json_exocytosis.DataInput.Dvesicle_cov_exocytosis, ...
                                 Json_exocytosis.DataInput.S_cov_exocytosis,...
                                 Json_exocytosis.DataInput.cov_scale_exocytosis,...
                                 Json_exocytosis.DataInput.kG_exocytosis, Json_exocytosis.DataInput.alpha_exocytosis, ...
                                 Json_exocytosis.DataInput.kp_exocytosis,Json_exocytosis.DataInput.beta_exocytosis, Json_exocytosis.DataInput.kinsulin_exocytosis, ...
                                 Json_exocytosis.DataInput.kD_exocytosis, Json_exocytosis.DataInput.kR_exocytosis,...
                                 Json_signaling.DataInput.G_mean_signaling, Json_signaling.DataInput.ATP_mean_signaling, Json_signaling.DataInput.GLP1_mean_signaling,...
                                 Json_signaling.DataInput.GLP1R_mean_signaling, Json_signaling.DataInput.cAMP_mean_signaling, Json_signaling.DataInput.Ca_mean_signaling,...
                                 Json_signaling.DataInput.S_mean_signaling, ...
                                 Json_signaling.DataInput.G_cov_signaling, Json_signaling.DataInput.ATP_cov_signaling, Json_signaling.DataInput.GLP1_cov_signaling, ...
                                 Json_signaling.DataInput.GLP1R_cov_signaling, Json_signaling.DataInput.cAMP_cov_signaling, Json_signaling.DataInput.Ca_cov_signaling, ...
                                 Json_signaling.DataInput.S_cov_signaling, Json_signaling.DataInput.cov_scale_signaling,...
                                 Json_signaling.DataInput.alpha_signaling, Json_signaling.DataInput.kATP_signaling, Json_signaling.DataInput.kGLP1_signaling,...
                                 Json_signaling.DataInput.beta_signaling, Json_signaling.DataInput.gamma_signaling, Json_signaling.DataInput.kcAMP_signaling,...
                                 Json_signaling.DataInput.delta_signaling, Json_signaling.DataInput.kCa_signaling, Json_signaling.DataInput.epsilon_signaling,...
                                 Json_screening.DataInput.A_mean_screening, Json_screening.DataInput.A_cov_screening, Json_screening.DataInput.conc_mean_screening,...
                                 Json_screening.DataInput.conc_conv_screening, Json_screening.DataInput.GLP1R_mean_screening, Json_screening.DataInput.GLP1R_cov_screening,...
                                 Json_screening.DataInput.cov_scale_screening, Json_screening.DataInput.k1_screening, Json_screening.DataInput.k2_screening,...
                                 Json_metabolism.DataInput.Gex_mean_metabolism, Json_metabolism.DataInput.Gex_cov_metabolism, Json_metabolism.DataInput.Ex4_mean_metabolism, ...
                                 Json_metabolism.DataInput.Ex4_cov_metabolism,...
                                 Json_metabolism.DataInput.Pathway_mean_metabolism, Json_metabolism.DataInput.Pathway_cov_metabolism, Json_metabolism.DataInput.ATP_mean_metabolism, ...
                                 Json_metabolism.DataInput.ATP_cov_metabolism, Json_metabolism.DataInput.cov_scale_metabolism,...
                                 Json_metabolism.DataInput.k1_metabolism, Json_metabolism.DataInput.k2_metabolism,...
                                 Json_metabolism.DataInput.k3_metabolism);
                
                npers= dbn.nnodes_per_slice;
                dbn_engine = jtree_dbn_inf_engine(dbn);
                
                % Time slices to sample
                evidence=cell(npers,Json_meta.DataInput.T);

                [engine, ll] = enter_evidence(dbn_engine, evidence);
                %disp(ll);
                marg = marginal_nodes(engine,nodes_map('Gb.postprandial'),Json_meta.DataInput.T);
                outputGb(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                marg = marginal_nodes(engine,nodes_map('kt.exocytosis'),Json_meta.DataInput.T);
                outputkt(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                
                marg = marginal_nodes(engine,nodes_map('G.postprandial'),Json_meta.DataInput.T);
                outputGex(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                marg = marginal_nodes(engine,nodes_map('G.exocytosis'),Json_meta.DataInput.T);
                outputGin(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                
                marg = marginal_nodes(engine,nodes_map('S.postprandial'),Json_meta.DataInput.T);
                outputSp(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                marg = marginal_nodes(engine,nodes_map('Spa.pancreas'),Json_meta.DataInput.T);
                outputSpa(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                marg = marginal_nodes(engine,nodes_map('Sis.pancreas'),Json_meta.DataInput.T);
                outputSis(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                marg = marginal_nodes(engine,nodes_map('Scell.pancreas'),Json_meta.DataInput.T);
                outputScell(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                marg = marginal_nodes(engine,nodes_map('S.exocytosis'),Json_meta.DataInput.T);
                outputSe(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
            end
        end
        T = table(inputkt,inputGb,outputkt,outputGb,outputGex,outputGin,outputSp,outputSpa,outputSis,outputScell,outputSe);
        % Write data to text file
        f1 = sprintf(OUTPUT_PATH+ 'metamodel_kt_mean%d_cov%d_Gb_mean%d_DataError.csv', i,j,m);
        writetable(T, f1)
    end
end


