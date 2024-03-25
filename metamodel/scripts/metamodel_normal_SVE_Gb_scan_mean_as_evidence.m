
function []=metamodel_normal_SVE_Gb_scan_mean_as_evidence()
% Sample the meta-modeling DBN containing following five input models:
% 1. The postprandial model
% 2. The exocytosis model
% 3. The pancreas model

SCRIPTS_PATH=pwd(); % Change this line if not running from scripts path
OUTPUT_PATH= SCRIPTS_PATH+"/Output/";
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
% metamodel
Json_meta = jsondecode(fileread('../data/meta_normal.json'));
Gb_SVE_mu_input = importdata(Json_meta.InputErr); % an nx2 array with possible values of Gb and kt to scan
Gb_SVE_cov_input = importdata(Json_meta.InputCov);
  
% Create a table with the data and variable names
inputS_VE = {};
inputS_VE(end+1,:) = {'INPUT_S_mu' 'INPUT_S_cov'};
outputGb_PR = {};
outputGb_PR(end+1,:) = {'Gb.postprandial_mu' 'Gb.postprandial_sigma2' 'Gb.postprandial_sigma'};
outputGex_PR = {};
outputGex_PR(end+1,:) = {'G.postprandial_mu' 'G.postprandial_sigma2' 'G.postprandial_sigma'};
outputGin_VE = {};
outputGin_VE(end+1,:) = {'G.exocytosis_mu' 'G.exocytosis_sigma2' 'G.exocytosis_sigma'};
outputSp_PR = {};
outputSp_PR(end+1,:) = {'S.postprandial_mu' 'S.postprandial_sigma2' 'S.postprandial_sigma'};
outputSpa_Pa = {};
outputSpa_Pa(end+1,:) = {'Spa.pancreas_mu' 'Spa.pancreas_sigma2' 'Spa.pancreas_sigma'};
outputSis_Pa = {};
outputSis_Pa(end+1,:) = {'Sis.pancreas_mu' 'Sis.pancreas_sigma2' 'Sis.pancreas_sigma'};
outputScell_Pa = {};
outputScell_Pa(end+1,:) = {'Scell.pancreas_mu' 'Scell.pancreas_sigma2' 'Scell.pancreas_sigma'};
outputkt_VE = {};
outputkt_VE(end+1,:) = {'kt.exocytosis_mu' 'kt.exocytosis_sigma2' 'kt.exocytosis_sigma'};

% disp(length(Gb_SVE_mu_input(:,1)));

% S_VE_mean_min_ndx = 1;
% S_VE_mean_max_ndx = length(Gb_SVE_mu_input(:,1));
% S_VE_cov_min_ndx = round(length(Gb_SVE_mu_input(:,1))/2);
% S_VE_cov_max_ndx = round(length(Gb_SVE_mu_input(:,1))/2);

% repeat for both err and sigma
S_VE_mean_min_ndx = 1;
S_VE_mean_max_ndx = length(Gb_SVE_mu_input(:,1));
S_VE_cov_min_ndx = 1;
S_VE_cov_max_ndx = length(Gb_SVE_mu_input(:,1));

Scell_init = readmatrix('../data/S_cell_obs.csv');
for i = S_VE_mean_min_ndx:S_VE_mean_max_ndx
    S_VE_mean_error = Gb_SVE_mu_input(i,1);
    for j = S_VE_cov_min_ndx:S_VE_cov_max_ndx
        S_VE_cov = Gb_SVE_cov_input(j,1);
        inputS_VE(end+1,:) = {S_VE_mean_error S_VE_cov};
        S_cell_obs_mean = Scell_init((i-1)*10+j,1);
        S_cell_obs_cov = Scell_init((i-1)*10+j,2);
%         disp(S_cell_obs_mean);
%         disp(S_cell_obs_cov);
        [dbn, nodes_map]= make_meta_dbn3(Json_postprandial.DataInput.DGd_mean_postprandial, Json_postprandial.DataInput.DGd_cov_postprandial, Json_postprandial.DataInput.Gb_mean_postprandial,...
                         Json_postprandial.DataInput.Gb_cov_postprandial, Json_postprandial.DataInput.G_mean_postprandial, Json_postprandial.DataInput.G_cov_postprandial,...
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
                         Json_exocytosis.DataInput.G_mean_exocytosis, Json_exocytosis.DataInput.kt_mean_exocytosis, Json_exocytosis.DataInput.Npatch_mean_exocytosis,...
                         Json_exocytosis.DataInput.Nvesicle_mean_exocytosis, Json_exocytosis.DataInput.Ninsulin_mean_exocytosis, Json_exocytosis.DataInput.Rcell_mean_exocytosis,...
                         Json_exocytosis.DataInput.Dvesicle_mean_exocytosis, ...
                         S_VE_mean_error, ...
                         Json_exocytosis.DataInput.G_cov_exocytosis,...
                         Json_exocytosis.DataInput.kt_mean_exocytosis, Json_exocytosis.DataInput.Npatch_cov_exocytosis, Json_exocytosis.DataInput.Nvesicle_cov_exocytosis, ...
                         Json_exocytosis.DataInput.Ninsulin_cov_exocytosis, Json_exocytosis.DataInput.Rcell_cov_exocytosis,Json_exocytosis.DataInput.Dvesicle_cov_exocytosis, ...
                         S_VE_cov, ...
                         Json_exocytosis.DataInput.cov_scale_exocytosis,...
                         Json_exocytosis.DataInput.kG_exocytosis, Json_exocytosis.DataInput.alpha_exocytosis, ...
                         Json_exocytosis.DataInput.kp_exocytosis,Json_exocytosis.DataInput.beta_exocytosis, Json_exocytosis.DataInput.kinsulin_exocytosis, ...
                         Json_exocytosis.DataInput.kD_exocytosis, Json_exocytosis.DataInput.kR_exocytosis, ...
                         Json_meta.DataInput.S_ve_weight, Json_meta.DataInput.G_pr_weight, ...
                         Json_meta.DataInput.G_c_weight, Json_meta.DataInput.S_pa_weight, ...
                         Json_meta.DataInput.S_cell_obs_weight, Json_meta.DataInput.G_obs_weight, ...
                         Json_meta.DataInput.G_cell_obs_weight, Json_meta.DataInput.S_pa_obs_weight, ...
                         S_cell_obs_mean, S_cell_obs_cov, ...
                         Json_meta.DataInput.S_pa_obs_mean, Json_meta.DataInput.S_pa_obs_cov, ...
                         Json_meta.DataInput.G_pl_obs_mean, Json_meta.DataInput.G_pl_obs_cov, ...
                         Json_meta.DataInput.G_cell_obs_mean, Json_meta.DataInput.G_cell_obs_cov);

        npers= dbn.nnodes_per_slice;
        dbn_engine = jtree_dbn_inf_engine(dbn);

        % Time slices to sample
        evidence=cell(npers,Json_meta.DataInput.T);
        [engine, ~] = enter_evidence(dbn_engine, evidence);
        %disp(ll);

        marg = marginal_nodes(engine,nodes_map('Gb.postprandial'),Json_meta.DataInput.T);
        outputGb_PR(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
        marg = marginal_nodes(engine,nodes_map('G.postprandial'),Json_meta.DataInput.T);
        outputGex_PR(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
        marg = marginal_nodes(engine,nodes_map('G.exocytosis'),Json_meta.DataInput.T);
        outputGin_VE(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
        marg = marginal_nodes(engine,nodes_map('S.postprandial'),Json_meta.DataInput.T);
        outputSp_PR(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
        marg = marginal_nodes(engine,nodes_map('Spa.pancreas'),Json_meta.DataInput.T);
        outputSpa_Pa(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
        marg = marginal_nodes(engine,nodes_map('Sis.pancreas'),Json_meta.DataInput.T);
        outputSis_Pa(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
        marg = marginal_nodes(engine,nodes_map('Scell.pancreas'),Json_meta.DataInput.T);
        outputScell_Pa(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
        marg = marginal_nodes(engine,nodes_map('kt.exocytosis'),Json_meta.DataInput.T);
        outputkt_VE(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};

        T = table(inputS_VE, outputGb_PR, outputGex_PR, outputGin_VE, ...
            outputSp_PR, outputSpa_Pa, outputSis_Pa, outputScell_Pa, outputkt_VE);
        % Write data to text file
        f1 = sprintf(OUTPUT_PATH+ 'metamodel_scan_input_SVE_mean_-8_10_cov_1e-21_1e-15_obs_0.2.csv', i,j);
        writetable(T, f1, 'WriteVariableNames',false)
    end
end