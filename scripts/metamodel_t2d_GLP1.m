% Sample the meta-modeling DBN containing following six input models:
% 1. The postprandial model
% 2. The pancreas model
% 3. The exocytosis model
% 4. The signaling model
% 5. The metabolism model
% 6. The screening model

warning('off','MATLAB:singularMatrix');


% ---------------------------------
% Read data as input and evidence
% ---------------------------------
% postprandial model
Json_postprandial = jsondecode(fileread('../data/postprandial_t2d.json'));
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
                                 Json_exocytosis.DataInput.kt_cov_exocytosis, Json_exocytosis.DataInput.Npatch_cov_exocytosis, Json_exocytosis.DataInput.Nvesicle_cov_exocytosis, ...
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
for GLP1_evidence = [0 Json_signaling.EvidenceInput.GL_GLP1_1 ...
                     Json_signaling.EvidenceInput.GL_GLP1_2 ...
                     Json_signaling.EvidenceInput.GL_GLP1_3 Json_signaling.EvidenceInput.GL_GLP1_4]
    evidence=cell(npers,Json_postprandial.DataInput.T);
    for measure = 1:Json_postprandial.DataInput.T
        evidence{nodes_map('DGd.obs'),measure} = EvidenceDG(measure);
        evidence{nodes_map('GLP1.obs'),measure} = GLP1_evidence; 
        %evidence{nodes_map('GLP1a.obs'),measure} = 0; 
        %evidence{nodes_map('conc.obs'),measure} = 0; 
        %evidence{nodes_map('Gculture.obs'),measure} = 0; 
        %evidence{nodes_map('Ex4.obs'),measure} = 0;
    end

    [engine, ll] = enter_evidence(dbn_engine, evidence);
    %disp(ll);

    % Create a table with the data and variable names
    T = table();

    for node_name = ["G.C" "Gcell.C" "Spa.C" "GLP1R.C" ...
                 "DGd.postprandial" "G.postprandial" "Gb.postprandial" ...
                 "S.postprandial" "Sb.postprandial" "I.postprandial"  ...
                 "Spa.pancreas" "Scell.pancreas" "G.exocytosis" ...
                 "kt.exocytosis" "Npatch.exocytosis" "Nvesicle.exocytosis" ...
                 "Ninsulin.exocytosis" "S.exocytosis" "G.signaling" ...
                 "cAMP.signaling" "Ca.signaling" "S.signaling" "GLP1R.signaling"]
        node_values = {};
        node_values(end+1,:) = {node_name,node_name,node_name}
        for slice = 1:Json_postprandial.DataInput.T
            marg = marginal_nodes(engine,nodes_map(node_name),slice);
            node_values(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
        end
        T = [T node_values];
    end

    % Write data to text file
    fname = sprintf('../results/GLP1_incretin/metamodel_t2d_GLP1_%d.csv', GLP1_evidence);
    writetable(T, fname);
end