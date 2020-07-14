% Sample the meta-modeling DBN containing following five input models:
% 1. The postprandial model
% 2. The exocytosis model
% 3. The signaling model
% 4. The GLP1R model
% 5. The metabolism model

% Add bnet
cd ../bnt-master
addpath(genpathKPM('../bnt-master'))
cd ../bnet_scripts_spt_G-S-k_Gb_k

warning('off','MATLAB:singularMatrix');
% ---------------------------------
% Read data as input and evidence
% ---------------------------------
% postprandial model
Json_postprandial = jsondecode(fileread('../data_highstd/postprandial_prior_normal.json'));
DGexp = importdata(Json_postprandial.EvidenceDG);
EvidenceDG = DGexp(:,2);

% pancreas model
Json_pancreas = jsondecode(fileread(['../data_highstd/pancreas_prior.json']));
% exocytosis model
Json_exocytosis = jsondecode(fileread(['../data_highstd/exocytosis_prior.json']));
% signaling model
Json_signaling = jsondecode(fileread('../data_highstd/signaling_prior_normal.json'));
% GLP1R model
Json_GLP1R = jsondecode(fileread('../data_highstd/GLP1R_prior.json'));
% metabolism model
Json_metabolism = jsondecode(fileread('../data_highstd/metabolism_prior.json'));

% meta model
Json_meta = jsondecode(fileread('../data_highstd/meta_normal.json'));
Gb_k_input = importdata(Json_meta.DataError21);
  
% Create a table with the data and variable names
inputGb = {};
inputGb(end+1,:) = {'Gb_mean_error' 'Gb_cov_error'};
inputk = {};
inputk(end+1,:) = {'k_mean_error' 'k_cov_error'};
outputGb = {};
outputGb(end+1,:) = {'Gb.postprandial' 'Gb.postprandial' 'Gb.postprandial'};
outputk = {};
outputk(end+1,:) = {'k.exocytosis' 'k.exocytosis' 'k.exocytosis'};

outputGex = {};
outputGex(end+1,:) = {'Gex.postprandial' 'Gex.postprandial' 'Gex.postprandial'};
outputGin = {};
outputGin(end+1,:) = {'Gin.exocytosis' 'Gin.exocytosis' 'Gin.exocytosis'};
outputSp = {};
outputSp(end+1,:) = {'S.postprandial' 'S.postprandial' 'S.postprandial'};
outputSe = {};
outputSe(end+1,:) = {'S.exocytosis' 'S.exocytosis' 'S.exocytosis'};

disp(length(Gb_k_input(:,1)));
k_mean_min_ndx = 11;
k_mean_max_ndx = 11;
k_cov_min_ndx = 11;
k_cov_max_ndx = 11;

for i = k_mean_min_ndx:  k_mean_max_ndx
    k_mean_error = Gb_k_input(i,3)
    for j = k_cov_min_ndx:  k_cov_max_ndx
        k_cov_error = Gb_k_input(j,4)
        for m = 1: length(Gb_k_input(:,1))
            Gb_mean_error = Gb_k_input(m,1)
            for n = 2:  length(Gb_k_input(:,1))
                Gb_cov_error = Gb_k_input(n,2)
                inputGb(end+1,:) = {Gb_mean_error Gb_cov_error};
                inputk(end+1,:) = {k_mean_error k_cov_error};
                [dbn, nodes_map]= make_meta_dbn6(Json_postprandial.DataInput.Gex_mean_postprandial, Json_postprandial.DataInput.Y_mean_postprandial, Json_postprandial.DataInput.alpha_mean_postprandial,...
                                   Json_postprandial.DataInput.beta_mean_postprandial, Json_postprandial.DataInput.gamma_mean_postprandial, Json_postprandial.DataInput.k1_mean_postprandial, ...
                                   Json_postprandial.DataInput.k2_mean_postprandial, Json_postprandial.DataInput.K_mean_postprandial, Json_postprandial.DataInput.dt_mean_postprandial_min, ...
                                   Gb_mean_error, Json_postprandial.DataInput.DGintake_mean_postprandial, Json_postprandial.DataInput.Sb_mean_postprandial, ...
                                   Json_postprandial.DataInput.I_mean_postprandial,...
                                   Json_postprandial.DataInput.Gex_cov_postprandial, Json_postprandial.DataInput.Y_cov_postprandial, Gb_cov_error, ...
                                   Json_postprandial.DataInput.DGintake_cov_postprandial, Json_postprandial.DataInput.Sb_cov_postprandial, Json_postprandial.DataInput.I_cov_postprandial,...
                                   Json_postprandial.DataInput.min_cov_postprandial,Json_postprandial.DataInput.G_minus_Gb_w_G_postprandial,Json_postprandial.DataInput.S_w_I_postprandial,...
                                   Json_pancreas.DataInput.Scell_mean_pancreas, Json_pancreas.DataInput.Scell_cov_pancreas, Json_pancreas.DataInput.Sislet_mean_pancreas, ...
                                   Json_pancreas.DataInput.Sislet_cov_pancreas, Json_pancreas.DataInput.Spancreas_mean_pancreas, Json_pancreas.DataInput.Spancreas_cov_pancreas, ...
                                   Json_pancreas.DataInput.min_cov_pancreas, Json_pancreas.DataInput.Scell_w_Sislet_pancreas, Json_pancreas.DataInput.Sislet_w_Spancreas_pancreas,...
                                   Json_exocytosis.DataInput.Gin_mean_exocytosis, k_mean_error, Json_exocytosis.DataInput.Npatch_mean_exocytosis,...
                                   Json_exocytosis.DataInput.Nisg_mean_exocytosis, Json_exocytosis.DataInput.Ninsulin_mean_exocytosis, Json_exocytosis.DataInput.Rpbc_mean_exocytosis,...
                                   Json_exocytosis.DataInput.Disg_mean_exocytosis, Json_exocytosis.DataInput.S_mean_exocytosis, Json_exocytosis.DataInput.Gin_cov_exocytosis,...
                                   k_cov_error, Json_exocytosis.DataInput.Npatch_cov_exocytosis, Json_exocytosis.DataInput.Nisg_cov_exocytosis, ...
                                   Json_exocytosis.DataInput.Ninsulin_cov_exocytosis, Json_exocytosis.DataInput.Rpbc_cov_exocytosis,Json_exocytosis.DataInput.Disg_cov_exocytosis, ...
                                   Json_exocytosis.DataInput.S_cov_exocytosis,...
                                   Json_exocytosis.DataInput.min_cov_exocytosis,...
                                   Json_exocytosis.DataInput.S_w_S_exocytosis, Json_exocytosis.DataInput.Gin_w_S_exocytosis, Json_exocytosis.DataInput.S_w_k_exocytosis, ...
                                   Json_exocytosis.DataInput.Npatch_w_S_exocytosis,Json_exocytosis.DataInput.S_w_Nisg_exocytosis, Json_exocytosis.DataInput.Ninsulin_w_S_exocytosis, ...
                                   Json_exocytosis.DataInput.Disg_w_S_exocytosis,Json_exocytosis.DataInput.Rpbc_w_S_exocytosis,...
                                   Json_signaling.DataInput.Gin_mean_signaling, Json_signaling.DataInput.ATP_mean_signaling, Json_signaling.DataInput.GLP1_mean_signaling,...
                                   Json_signaling.DataInput.GLP1R_mean_signaling, Json_signaling.DataInput.cAMP_mean_signaling, Json_signaling.DataInput.Ca_mean_signaling,...
                                   Json_signaling.DataInput.S_mean_signaling, ...
                                   Json_signaling.DataInput.Gin_cov_signaling, Json_signaling.DataInput.ATP_cov_signaling, Json_signaling.DataInput.GLP1_cov_signaling, ...
                                   Json_signaling.DataInput.GLP1R_cov_signaling, Json_signaling.DataInput.cAMP_cov_signaling, Json_signaling.DataInput.Ca_cov_signaling, ...
                                   Json_signaling.DataInput.S_cov_signaling, Json_signaling.DataInput.min_cov_signaling,...
                                   Json_signaling.DataInput.Gin_w_ATP_signaling, Json_signaling.DataInput.ATP_w_ATP_signaling, Json_signaling.DataInput.GLP1_w_GLP1R_signaling,...
                                   Json_signaling.DataInput.ATP_w_cAMP_signaling, Json_signaling.DataInput.GLP1R_w_cAMP_signaling, Json_signaling.DataInput.cAMP_w_cAMP_signaling,...
                                   Json_signaling.DataInput.cAMP_w_Ca_signaling, Json_signaling.DataInput.Ca_w_Ca_signaling, Json_signaling.DataInput.Ca_w_S_signaling,...
                                   Json_GLP1R.DataInput.GLP1a_mean_GLP1R, Json_GLP1R.DataInput.GLP1a_cov_GLP1R, Json_GLP1R.DataInput.conc_mean_GLP1R,...
                                   Json_GLP1R.DataInput.conc_conv_GLP1R, Json_GLP1R.DataInput.GLP1R_mean_GLP1R, Json_GLP1R.DataInput.GLP1R_cov_GLP1R,...
                                   Json_GLP1R.DataInput.min_cov_GLP1R, Json_GLP1R.DataInput.GLP1a_w_GLP1R_GLP1R, Json_GLP1R.DataInput.conc_w_GLP1R_GLP1R,...
                                   Json_metabolism.DataInput.Gculture_mean_metabolism, Json_metabolism.DataInput.Gculture_cov_metabolism, Json_metabolism.DataInput.Ex4_mean_metabolism, ...
                                   Json_metabolism.DataInput.Ex4_cov_metabolism, Json_metabolism.DataInput.NES_mean_metabolism, Json_metabolism.DataInput.NES_cov_metabolism, ...
                                   Json_metabolism.DataInput.Pathway_mean_metabolism, Json_metabolism.DataInput.Pathway_cov_metabolism, Json_metabolism.DataInput.ATP_mean_metabolism,...
                                   Json_metabolism.DataInput.ATP_cov_metabolism, Json_metabolism.DataInput.min_cov_metabolism,...
                                   Json_metabolism.DataInput.Gculture_w_NES_metabolism, Json_metabolism.DataInput.Ex4_w_NES_metabolism, Json_metabolism.DataInput.NES_w_Pathway_metabolism,...
                                   Json_metabolism.DataInput.Pathway_w_ATP_metabolism);
                
                npers= dbn.nnodes_per_slice;
                dbn_engine = jtree_dbn_inf_engine(dbn);
                
                % Time slices to sample
                evidence=cell(npers,Json_meta.DataInput.T);

                [engine, ll] = enter_evidence(dbn_engine, evidence);
                %disp(ll);
                marg = marginal_nodes(engine,nodes_map('Gb.postprandial'),Json_meta.DataInput.T);
                outputGb(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                marg = marginal_nodes(engine,nodes_map('k.exocytosis'),Json_meta.DataInput.T);
                outputk(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                
                marg = marginal_nodes(engine,nodes_map('Gex.postprandial'),Json_meta.DataInput.T);
                outputGex(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                marg = marginal_nodes(engine,nodes_map('Gin.exocytosis'),Json_meta.DataInput.T);
                outputGin(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                
                marg = marginal_nodes(engine,nodes_map('S.postprandial'),Json_meta.DataInput.T);
                outputSp(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
                marg = marginal_nodes(engine,nodes_map('S.exocytosis'),Json_meta.DataInput.T);
                outputSe(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
            end
        end
        T = table(inputGb,inputk,outputGb,outputk,outputGex,outputGin,outputSp,outputSe);
        % Write data to text file
        f1 = sprintf('../Gb_k_results_highstd/metamodel_k_mean%d_cov%d_Gb_DataError.csv', i,j);
        writetable(T, f1)
    end
end


