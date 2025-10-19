% Sample the meta-modeling DBN containing following five input models:
% 1. The postprandial model
% 2. The exocytosis model
% 3. The signaling model
% 4. The GLP1R model
% 5. The metabolism model

% Add bnet
cd ../bnt_master
addpath(genpathKPM('../bnt_master'))
cd ../scripts


warning('off','MATLAB:singularMatrix');

% ---------------------------------
% Read data as input and evidence
% ---------------------------------
% postprandial model
Json_postprandial = jsondecode(fileread(['../data/surrogate/postprandial_normal_s1.json']));
DGexp = importdata(Json_postprandial.EvidenceDG);
EvidenceDG = DGexp(:,1);

% pancreas model
Json_pancreas = jsondecode(fileread(['../data/surrogate/pancreas.json']));
% exocytosis model
Json_exocytosis = jsondecode(fileread(['../data/surrogate/exocytosis.json']));
% meta model
Json_meta = jsondecode(fileread(['../data/metamodel/meta_normal_fig5.json']));

EvidenceGPRerror = importdata(Json_meta.InputErr); % an nx2 array with possible values of Gb and kt to scan
Gb_SVE_cov_input = importdata(Json_meta.InputCov);

for i = 1:11
     error_scale = 0.5+(i-1)/10;
    for j = 1:11
        G_PR_cov = Gb_SVE_cov_input(j,1);

        
        [dbn, nodes_map]= make_meta_dbn3(Json_postprandial.DataInput.DGd_mean_postprandial, Json_postprandial.DataInput.DGd_cov_postprandial, ...
                                         Json_postprandial.DataInput.Gb_mean_postprandial, Json_postprandial.DataInput.Gb_cov_postprandial,...
                                         Json_postprandial.DataInput.G_mean_postprandial, G_PR_cov, ...
                                         Json_postprandial.DataInput.DG_mean_postprandial, Json_postprandial.DataInput.DG_cov_postprandial, Json_postprandial.DataInput.Y_mean_postprandial,...
                                         Json_postprandial.DataInput.Y_cov_postprandial, Json_postprandial.DataInput.S_mean_postprandial, Json_postprandial.DataInput.S_cov_postprandial,...
                                         Json_postprandial.DataInput.I_mean_postprandial, Json_postprandial.DataInput.I_cov_postprandial, Json_postprandial.DataInput.Sb_mean_postprandial,...
                                         Json_postprandial.DataInput.Sb_cov_postprandial, Json_postprandial.DataInput.alpha_postprandial, Json_postprandial.DataInput.beta_postprandial,...
                                         Json_postprandial.DataInput.gamma_postprandial, Json_postprandial.DataInput.k1_postprandial,Json_postprandial.DataInput.k2_postprandial,...
                                         Json_postprandial.DataInput.k3_postprandial, Json_postprandial.DataInput.k4_postprandial, Json_postprandial.DataInput.K_postprandial,...
                                         Json_postprandial.DataInput.dt_postprandial, Json_postprandial.DataInput.cov_scale_postprandial, ...
                                         Json_pancreas.DataInput.Scell_mean_pancreas, Json_pancreas.DataInput.Scell_cov_pancreas, Json_pancreas.DataInput.Sislet_mean_pancreas, ...
                                         Json_pancreas.DataInput.Sislet_cov_pancreas, Json_pancreas.DataInput.Spancreas_mean_pancreas, Json_pancreas.DataInput.Spancreas_cov_pancreas, ...
                                         Json_pancreas.DataInput.cov_scale_pancreas, Json_pancreas.DataInput.Nc_pancreas, Json_pancreas.DataInput.Ni_pancreas,...
                                         Json_exocytosis.DataInput.G_mean_exocytosis, Json_exocytosis.DataInput.kt_mean_exocytosis, Json_exocytosis.DataInput.Npatch_mean_exocytosis,...
                                         Json_exocytosis.DataInput.Nvesicle_mean_exocytosis, Json_exocytosis.DataInput.Ninsulin_mean_exocytosis, Json_exocytosis.DataInput.Rcell_mean_exocytosis,...
                                         Json_exocytosis.DataInput.Dvesicle_mean_exocytosis, ...
                                         Json_exocytosis.DataInput.S_mean_exocytosis, ...
                                         Json_exocytosis.DataInput.G_cov_exocytosis,Json_exocytosis.DataInput.kt_cov_exocytosis, Json_exocytosis.DataInput.Npatch_cov_exocytosis, Json_exocytosis.DataInput.Nvesicle_cov_exocytosis, ...
                                         Json_exocytosis.DataInput.Ninsulin_cov_exocytosis, Json_exocytosis.DataInput.Rcell_cov_exocytosis,Json_exocytosis.DataInput.Dvesicle_cov_exocytosis, ...
                                         Json_exocytosis.DataInput.S_cov_exocytosis,...
                                         Json_exocytosis.DataInput.cov_scale_exocytosis,...
                                         Json_exocytosis.DataInput.kG_exocytosis, Json_exocytosis.DataInput.alpha_exocytosis, ...
                                         Json_exocytosis.DataInput.kp_exocytosis,Json_exocytosis.DataInput.beta_exocytosis, Json_exocytosis.DataInput.kinsulin_exocytosis, ...
                                         Json_exocytosis.DataInput.kD_exocytosis, Json_exocytosis.DataInput.kR_exocytosis, ...
                                         Json_meta.DataInput.S_ve_weight, Json_meta.DataInput.G_pr_weight, ...
                                         Json_meta.DataInput.G_c_weight, Json_meta.DataInput.S_pa_weight, ...
                                         Json_meta.DataInput.S_cell_obs_weight, Json_meta.DataInput.G_obs_weight, ...
                                         Json_meta.DataInput.G_cell_obs_weight, Json_meta.DataInput.S_pa_obs_weight, ...
                                         Json_meta.DataInput.S_cell_obs_mean, Json_meta.DataInput.S_cell_obs_cov, ...
                                         Json_meta.DataInput.S_pa_obs_mean, Json_meta.DataInput.S_pa_obs_cov, ...
                                         Json_meta.DataInput.G_pl_obs_mean, Json_meta.DataInput.G_pl_obs_cov, ...
                                         Json_meta.DataInput.G_cell_obs_mean, Json_meta.DataInput.G_cell_obs_cov);


        npers= dbn.nnodes_per_slice;
        dbn_engine = jtree_dbn_inf_engine(dbn);     
        
        % Time slices to sample
        evidence=cell(npers,Json_meta.DataInput.T);
        
        G_PR_i= nodes_map('Gpr.obs');
        for T_i=1:Json_meta.DataInput.T
            evidence{G_PR_i, T_i}= EvidenceGPRerror(T_i) * error_scale;
        end
        
        Gpl_C = importdata(Json_meta.EvidenceInput.Evidence_Gpl);
        for measure = 1:Json_meta.DataInput.T
          evidence{nodes_map('G.obs'),measure} = Gpl_C(measure);
        end
        
        % % observations for the coupling variable
        Scell_C = importdata(Json_meta.EvidenceInput.Evidence_Scell);
        for measure = 1:Json_meta.DataInput.T
           evidence{nodes_map('Scell.obs'), measure} = Scell_C(measure);
        end
        % % disp(Scell_C);

        Spa_C = importdata(Json_meta.EvidenceInput.Evidence_Spa);
        for measure = 1:Json_meta.DataInput.T   
           evidence{nodes_map('Spa.obs'),measure} = Spa_C(measure);
        end
        % disp(Spa_C);

        Gcell_C = importdata(Json_meta.EvidenceInput.Evidence_Gcell);
        for measure = 1:Json_meta.DataInput.T
          evidence{nodes_map('Gcell.obs'),measure} = Gcell_C(measure);
        end
        %disp(Gcell_C);

        [engine, ll] = enter_evidence(dbn_engine, evidence);

        % Create a table with the data and variable names
        T = table();

        for node_name = ["DG.postprandial" "DGd.postprandial" "G.postprandial" "Gb.postprandial" ...
                         "S.postprandial" "Sb.postprandial" "I.postprandial" "Y.postprandial" ...
                         "Spa.pancreas" "Scell.pancreas" "Sis.pancreas" ...
                         "Rcell.exocytosis" "G.exocytosis" "kt.exocytosis" "Npatch.exocytosis" ...
                         "Nvesicle.exocytosis" "Dvesicle.exocytosis" "Ninsulin.exocytosis" "S.exocytosis"...
                         "Gcell.C" "G.C" "Scell.C" "Spa.C"]
            node_values = {};
            node_values(end+1,:) = {node_name,node_name,node_name}
            for slice = 1:Json_meta.DataInput.T
                marg = marginal_nodes(engine,nodes_map(node_name),slice);
                node_values(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
            end
            T = [T node_values];
        end

        % Write data to text file
        fn = ['../Output/fig5_test_Gpl/metamodel_scan_input_G_PR_test_mean_scale_' num2str(error_scale) '_cov_' num2str(G_PR_cov) '.csv'];
        fname = sprintf(fn);
        writetable(T, fname);   
    end
end
