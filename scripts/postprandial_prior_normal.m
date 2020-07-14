% Make a DBN for the postprandial model with the following variables
% Add bnet
cd ../bnt-master
addpath(genpathKPM('../bnt-master'))
cd ../bnet_scripts_spt_G-S-k

warning('off','MATLAB:singularMatrix');

% ---------------------------------
% Read data as input and evidence
% ---------------------------------
Json_postprandial = jsondecode(fileread('../data/postprandial_prior_normal.json'));
DGexp = importdata(Json_postprandial.EvidenceDG);
EvidenceDG = DGexp(:,2);
%disp(EvidenceDG)

%return;
[postprandial_dbn_factory]= make_postprandial_dbn(Json_postprandial.DataInput.Gex_mean_postprandial, Json_postprandial.DataInput.Y_mean_postprandial, Json_postprandial.DataInput.alpha_mean_postprandial,...
                                          Json_postprandial.DataInput.beta_mean_postprandial, Json_postprandial.DataInput.gamma_mean_postprandial, Json_postprandial.DataInput.k1_mean_postprandial, ...
                                          Json_postprandial.DataInput.k2_mean_postprandial, Json_postprandial.DataInput.K_mean_postprandial, Json_postprandial.DataInput.dt_mean_postprandial_min, ...
                                          Json_postprandial.DataInput.Gb_mean_postprandial, Json_postprandial.DataInput.DGintake_mean_postprandial, Json_postprandial.DataInput.Sb_mean_postprandial, ...
                                          Json_postprandial.DataInput.I_mean_postprandial,...
                                          Json_postprandial.DataInput.Gex_cov_postprandial, Json_postprandial.DataInput.Y_cov_postprandial, Json_postprandial.DataInput.Gb_cov_postprandial, ...
                                          Json_postprandial.DataInput.DGintake_cov_postprandial, Json_postprandial.DataInput.Sb_cov_postprandial, Json_postprandial.DataInput.I_cov_postprandial,...
                                          Json_postprandial.DataInput.min_cov_postprandial,Json_postprandial.DataInput.G_minus_Gb_w_G_postprandial,Json_postprandial.DataInput.S_w_I_postprandial);

[dbn, intra, inter, nodes_map] = create_dbn(postprandial_dbn_factory);
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);

evidence=cell(npers,Json_postprandial.DataInput.T);

for measure = 1:Json_postprandial.DataInput.T
    evidence{nodes_map('DGintake.obs'),measure} = EvidenceDG(measure);
end

[engine, ll] = enter_evidence(dbn_engine, evidence);
%disp(ll);

%disp(nodes_map('DGintake.postprandial'));
%keydnx = cellfun(@(x)isequal(x,9),values(nodes_map));
%disp(keys(keys);)

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
writetable(T, '../results/postprandial_prior_normal.csv');
