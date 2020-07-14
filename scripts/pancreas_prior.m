% Make a DBN for the pancreas model with the following variables
%

warning('off','MATLAB:singularMatrix');

% Read data as input and evidence.
Json_pancreas = jsondecode(fileread('../data/pancreas_prior.json'));

[pancreas_dbn_factory]= make_pancreas_dbn(Json_pancreas.DataInput.Scell_mean_pancreas, Json_pancreas.DataInput.Scell_cov_pancreas, Json_pancreas.DataInput.Sislet_mean_pancreas, ...
                                          Json_pancreas.DataInput.Sislet_cov_pancreas, Json_pancreas.DataInput.Spancreas_mean_pancreas, Json_pancreas.DataInput.Spancreas_cov_pancreas, ...
                                          Json_pancreas.DataInput.min_cov_pancreas, Json_pancreas.DataInput.Scell_w_Sislet_pancreas, Json_pancreas.DataInput.Sislet_w_Spancreas_pancreas)

[dbn, intra, inter, nodes_map] = create_dbn(pancreas_dbn_factory);
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);

evidence=cell(npers,Json_pancreas.DataInput.T);

[engine, ll] = enter_evidence(dbn_engine, evidence);
disp(ll);

%disp(nodes_map('DGintake.pancreas'));
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
    for slice = 1:Json_pancreas.DataInput.T
        marg = marginal_nodes(engine,node_ndx,slice);
        node_values(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
    end
    T = [T node_values];
end 

% Write data to text file
writetable(T, '../results/pancreas_prior_normal.csv');
