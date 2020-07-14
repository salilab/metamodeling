% Make a DBN for the exocytosis model with the following variables
%

warning('off','MATLAB:singularMatrix');

% ---------------------------------
% Read data as input and evidence
% ---------------------------------
Json_exocytosis = jsondecode(fileread('../data/exocytosis_prior.json'));

[exocytosis_dbn_factory]= make_exocytosis_dbn(Json_exocytosis.DataInput.Gin_mean_exocytosis, Json_exocytosis.DataInput.k_mean_exocytosis, Json_exocytosis.DataInput.Npatch_mean_exocytosis,...
                                              Json_exocytosis.DataInput.Nisg_mean_exocytosis, Json_exocytosis.DataInput.Ninsulin_mean_exocytosis, Json_exocytosis.DataInput.Rpbc_mean_exocytosis,...
                                              Json_exocytosis.DataInput.Disg_mean_exocytosis, Json_exocytosis.DataInput.S_mean_exocytosis, Json_exocytosis.DataInput.Gin_cov_exocytosis,...
                                              Json_exocytosis.DataInput.k_cov_exocytosis, Json_exocytosis.DataInput.Npatch_cov_exocytosis, Json_exocytosis.DataInput.Nisg_cov_exocytosis, ...
                                              Json_exocytosis.DataInput.Ninsulin_cov_exocytosis, Json_exocytosis.DataInput.Rpbc_cov_exocytosis,Json_exocytosis.DataInput.Disg_cov_exocytosis, ...
                                              Json_exocytosis.DataInput.S_cov_exocytosis,...
                                              Json_exocytosis.DataInput.min_cov_exocytosis,...
                                              Json_exocytosis.DataInput.S_w_S_exocytosis, Json_exocytosis.DataInput.Gin_w_S_exocytosis, Json_exocytosis.DataInput.S_w_k_exocytosis, ...
                                              Json_exocytosis.DataInput.Npatch_w_S_exocytosis,Json_exocytosis.DataInput.S_w_Nisg_exocytosis, Json_exocytosis.DataInput.Ninsulin_w_S_exocytosis, ...
                                              Json_exocytosis.DataInput.Disg_w_S_exocytosis,Json_exocytosis.DataInput.Rpbc_w_S_exocytosis);

[dbn, intra, inter, nodes_map] = create_dbn(exocytosis_dbn_factory);
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);

evidence=cell(npers,Json_exocytosis.DataInput.T);

[engine, ll] = enter_evidence(dbn_engine, evidence);
disp(ll);

%disp(nodes_map('DGintake.exocytosis'));
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
    for slice = 1:Json_exocytosis.DataInput.T
        marg = marginal_nodes(engine,node_ndx,slice);
        node_values(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
    end
    T = [T node_values];
end 

% Write data to text file
writetable(T, '../results/exocytosis_prior_normal.csv');
