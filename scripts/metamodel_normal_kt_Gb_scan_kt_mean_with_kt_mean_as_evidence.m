function []=metamodel_normal_kt_Gb_scan_kt_mean_with_kt_mean_as_evidence()
% Sample the meta-modeling DBN containing following five input models:
% 1. The postprandial model
% 2. The exocytosis model
% 3. The signaling model
% 4. The GLP1R model
% 5. The metabolism model

SCRIPTS_PATH=pwd(); % Change this line if not running from scripts path
OUTPUT_PATH= SCRIPTS_PATH+"../Output/accuracy_precision/Gb_kt/";
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
inputkt_VE = {};
inputkt_VE(end+1,:) = {'INPUT_kt_mu' 'INPUT_kt_cov'};
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
outputSe_VE = {};
outputSe_VE(end+1,:) = {'S.exocytosis_mu' 'S.exocytosis_sigma2' 'S.exocytosis_sigma'};

disp(length(Gb_kt_mu_input(:,1)));

kt_VE_mean_min_ndx = 1;
kt_VE_mean_max_ndx = length(Gb_kt_mu_input(:,2));
kt_VE_cov_min_ndx = round(length(Gb_kt_mu_input(:,2))/2);
kt_VE_cov_max_ndx = round(length(Gb_kt_mu_input(:,2))/2);


for i = kt_VE_mean_min_ndx:kt_VE_mean_max_ndx
    kt_VE_mean_error_evidence = Gb_kt_mu_input(i,2);
    for j = kt_VE_cov_min_ndx:kt_VE_cov_max_ndx
        kt_VE_cov = Gb_kt_sigma_input(j,2);
        inputkt_VE(end+1,:) = {kt_VE_mean_error_evidence kt_VE_cov};
        [dbn, nodes_map]= make_meta_dbn6(Json_postprandial.DataInput.DGd_mean_postprandial, Json_postprandial.DataInput.DGd_cov_postprandial, Json_postprandial.DataInput.Gb_mean_postprandial,...
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
                         Json_exocytosis.DataInput.Dvesicle_mean_exocytosis, Json_exocytosis.DataInput.S_mean_exocytosis, Json_exocytosis.DataInput.G_cov_exocytosis,...
                         kt_VE_cov, Json_exocytosis.DataInput.Npatch_cov_exocytosis, Json_exocytosis.DataInput.Nvesicle_cov_exocytosis, ...
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
        kt_VE_i= nodes_map('kt.exocytosis');
        for T_i=1:Json_meta.DataInput.T
            evidence{kt_VE_i, T_i}= kt_VE_mean_error_evidence;
        end

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
        marg = marginal_nodes(engine,nodes_map('S.exocytosis'),Json_meta.DataInput.T);
        outputSe_VE(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
 
        T = table(inputkt_VE, outputGb_PR, outputGex_PR,outputGin_VE, ...
            outputSp_PR, outputSpa_Pa, outputSis_Pa, outputScell_Pa, outputSe_VE);
        % Write data to text file
        f1 = sprintf(OUTPUT_PATH+ 'metamodel_scan_input_kt_mean_output_Gb_with_kt_mean_as_evidence.csv', i,j);
        writetable(T, f1, 'WriteVariableNames',false)
    end
end


