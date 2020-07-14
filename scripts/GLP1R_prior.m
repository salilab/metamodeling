% Make a DBN for the GLP1R model with the following variables
%

warning('off','MATLAB:singularMatrix');

% Read data as input and evidence.
Json_GLP1R = jsondecode(fileread('../data/GLP1R_prior.json'));

[GLP1R_dbn_factory]= make_GLP1R_dbn(Json_GLP1R.DataInput.GLP1a_mean_GLP1R, Json_GLP1R.DataInput.GLP1a_cov_GLP1R, Json_GLP1R.DataInput.conc_mean_GLP1R,...
                                    Json_GLP1R.DataInput.conc_conv_GLP1R, Json_GLP1R.DataInput.GLP1R_mean_GLP1R, Json_GLP1R.DataInput.GLP1R_cov_GLP1R,...
                                    Json_GLP1R.DataInput.min_cov_GLP1R, Json_GLP1R.DataInput.GLP1a_w_GLP1R_GLP1R, Json_GLP1R.DataInput.conc_w_GLP1R_GLP1R)

[dbn, intra, inter, nodes_map] = create_dbn(GLP1R_dbn_factory);
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);

evidence=cell(npers,Json_GLP1R.DataInput.T);

[engine, ll] = enter_evidence(dbn_engine, evidence);
disp(ll);

%disp(nodes_map('DGintake.GLP1R'));
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
    for slice = 1:Json_GLP1R.DataInput.T
        marg = marginal_nodes(engine,node_ndx,slice);
        node_values(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
    end
    T = [T node_values];
end 

% Write data to text file
writetable(T, '../results/GLP1R_prior_normal.csv');
