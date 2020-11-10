% Make a DBN for the virtual screening model with the following variables
%

warning('off','MATLAB:singularMatrix');

% Read data as input and evidence.
Json_screening = jsondecode(fileread('../data/screening.json'));

[screening_dbn_factory]= make_screening_dbn(Json_screening.DataInput.A_mean_screening, Json_screening.DataInput.A_cov_screening, Json_screening.DataInput.conc_mean_screening,...
                                    Json_screening.DataInput.conc_conv_screening, Json_screening.DataInput.GLP1R_mean_screening, Json_screening.DataInput.GLP1R_cov_screening,...
                                    Json_screening.DataInput.cov_scale_screening, Json_screening.DataInput.k1_screening, Json_screening.DataInput.k2_screening)

[dbn, intra, inter, nodes_map] = create_dbn(screening_dbn_factory);
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);

evidence=cell(npers,Json_screening.DataInput.T);

[engine, ll] = enter_evidence(dbn_engine, evidence);
disp(ll);

%disp(nodes_map('DGintake.screening'));
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
    for slice = 1:Json_screening.DataInput.T
        marg = marginal_nodes(engine,node_ndx,slice);
        node_values(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
    end
    T = [T node_values];
end 

% Write data to text file
writetable(T, '../results/models/screening.csv');
