% Make a DBN for the postprandial model S3 with enhanced edges
% Add bnet
cd ../bnt_master
addpath(genpathKPM('../bnt_master'))
cd ../scripts
warning('off','MATLAB:singularMatrix');

% ---------------------------------
% Read data as input and evidence
% ---------------------------------
Json_postprandial = jsondecode(fileread('../data/surrogate/postprandial_normal_s3.json'));
DGexp = importdata(Json_postprandial.EvidenceDG);
EvidenceDG = DGexp(:,2);

[postprandial_dbn_factory]= make_postprandial_dbn_s3(Json_postprandial.DataInput.DGd_mean_postprandial, Json_postprandial.DataInput.DGd_cov_postprandial, Json_postprandial.DataInput.Gb_mean_postprandial,...
                                              Json_postprandial.DataInput.Gb_cov_postprandial, Json_postprandial.DataInput.G_mean_postprandial, Json_postprandial.DataInput.G_cov_postprandial,...
                                              Json_postprandial.DataInput.DG_mean_postprandial, Json_postprandial.DataInput.DG_cov_postprandial, Json_postprandial.DataInput.Y_mean_postprandial,...
                                              Json_postprandial.DataInput.Y_cov_postprandial, Json_postprandial.DataInput.S_mean_postprandial, Json_postprandial.DataInput.S_cov_postprandial,...
                                              Json_postprandial.DataInput.I_mean_postprandial, Json_postprandial.DataInput.I_cov_postprandial, Json_postprandial.DataInput.Sb_mean_postprandial,...
                                              Json_postprandial.DataInput.Sb_cov_postprandial, Json_postprandial.DataInput.alpha_postprandial, Json_postprandial.DataInput.beta_postprandial,...
                                              Json_postprandial.DataInput.gamma_postprandial, Json_postprandial.DataInput.k1_postprandial,Json_postprandial.DataInput.k2_postprandial,...
                                              Json_postprandial.DataInput.k3_postprandial, Json_postprandial.DataInput.k4_postprandial, Json_postprandial.DataInput.K_postprandial,...
                                              Json_postprandial.DataInput.dt_postprandial, Json_postprandial.DataInput.cov_scale_postprandial);

[dbn, intra, inter, nodes_map] = create_dbn(postprandial_dbn_factory);
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);

evidence=cell(npers,Json_postprandial.DataInput.T);

for measure = 1:Json_postprandial.DataInput.T
    evidence{nodes_map('DGd.obs'),measure} = EvidenceDG(measure);
end

[engine, ll] = enter_evidence(dbn_engine, evidence);

% Display network statistics
fprintf('============================================\n');
fprintf('Postprandial Normal S3 Model Statistics\n');
fprintf('============================================\n');
fprintf('Number of nodes per time slice: %d\n', npers);
fprintf('Number of time slices: %d\n', Json_postprandial.DataInput.T);

% Count edges in intra
num_edges_intra = size(postprandial_dbn_factory.edges_intra, 1);
fprintf('Number of intra-slice edges: %d\n', num_edges_intra);

% Count edges in inter
num_edges_inter = size(postprandial_dbn_factory.edges_inter, 1);
fprintf('Number of inter-slice edges: %d\n', num_edges_inter);

fprintf('Total edges: %d\n', num_edges_intra + num_edges_inter);
fprintf('============================================\n');

% writing model variables in a .json file for all time slices.
keys =  keys(nodes_map);

% Create a table with the data and variable names
T = table();

for node_ndx = 1:npers
    keydnx = cellfun(@(x)isequal(x,node_ndx),values(nodes_map));
    node_name = keys(keydnx);
    node_values = {};
    node_values(end+1,:) = {node_name,node_name,node_name}
    for slice = 1:Json_postprandial.DataInput.T 
        marg = marginal_nodes(engine,node_ndx,slice);
        node_values(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
    end
    T = [T node_values];
end 

% Write data to text file
writetable(T, ['../Output/postprandial_prior_normal_wDGd_s3.csv']);
fprintf('Output saved to: ../Output/postprandial_prior_normal_wDGd_s3.csv\n');
fprintf('S3 model completed successfully!\n');

