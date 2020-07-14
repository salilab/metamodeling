% Make a DBN for the metabolism model with the following variables
%

warning('off','MATLAB:singularMatrix');

% Read data as input and evidence.
Json_metabolism = jsondecode(fileread('../data/metabolism_prior.json'));

[metabolism_dbn_factory]= make_metabolism_dbn(Json_metabolism.DataInput.Gculture_mean_metabolism, Json_metabolism.DataInput.Gculture_cov_metabolism, Json_metabolism.DataInput.Ex4_mean_metabolism, ...
                                              Json_metabolism.DataInput.Ex4_cov_metabolism, Json_metabolism.DataInput.NES_mean_metabolism, Json_metabolism.DataInput.NES_cov_metabolism, ...
                                              Json_metabolism.DataInput.Pathway_mean_metabolism, Json_metabolism.DataInput.Pathway_cov_metabolism, Json_metabolism.DataInput.ATP_mean_metabolism, ...
                                              Json_metabolism.DataInput.ATP_cov_metabolism, Json_metabolism.DataInput.min_cov_metabolism,...
                                              Json_metabolism.DataInput.Gculture_w_NES_metabolism, Json_metabolism.DataInput.Ex4_w_NES_metabolism, Json_metabolism.DataInput.NES_w_Pathway_metabolism,...
                                              Json_metabolism.DataInput.Pathway_w_ATP_metabolism)

[dbn, intra, inter, nodes_map] = create_dbn(metabolism_dbn_factory);
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);

evidence=cell(npers,Json_metabolism.DataInput.T);

[engine, ll] = enter_evidence(dbn_engine, evidence);
disp(ll);

%disp(nodes_map('DGintake.metabolism'));
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
    for slice = 1:Json_metabolism.DataInput.T
        marg = marginal_nodes(engine,node_ndx,slice);
        node_values(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
    end
    T = [T node_values];
end 

% Write data to text file
writetable(T, '../results/metabolism_prior_normal.csv');
