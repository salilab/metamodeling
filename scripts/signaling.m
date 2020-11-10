% Make a DBN for the signaling model with the following variables
%

warning('off','MATLAB:singularMatrix');

% ---------------------------------
% Read data as input and evidence
% ---------------------------------
Json_signaling = jsondecode(fileread('../data/signaling.json'));

[signaling_dbn_factory]= make_signaling_dbn(Json_signaling.DataInput.G_mean_signaling, Json_signaling.DataInput.ATP_mean_signaling, Json_signaling.DataInput.GLP1_mean_signaling,...
                                           Json_signaling.DataInput.GLP1R_mean_signaling, Json_signaling.DataInput.cAMP_mean_signaling, Json_signaling.DataInput.Ca_mean_signaling,...
                                           Json_signaling.DataInput.S_mean_signaling, ...
                                           Json_signaling.DataInput.G_cov_signaling, Json_signaling.DataInput.ATP_cov_signaling, Json_signaling.DataInput.GLP1_cov_signaling, ...
                                           Json_signaling.DataInput.GLP1R_cov_signaling, Json_signaling.DataInput.cAMP_cov_signaling, Json_signaling.DataInput.Ca_cov_signaling, ...
                                           Json_signaling.DataInput.S_cov_signaling, Json_signaling.DataInput.cov_scale_signaling,...
                                           Json_signaling.DataInput.alpha_signaling, Json_signaling.DataInput.kATP_signaling, Json_signaling.DataInput.kGLP1_signaling,...
                                           Json_signaling.DataInput.beta_signaling, Json_signaling.DataInput.gamma_signaling, Json_signaling.DataInput.kcAMP_signaling,...
                                           Json_signaling.DataInput.delta_signaling, Json_signaling.DataInput.kCa_signaling, Json_signaling.DataInput.epsilon_signaling)

[dbn, intra, inter, nodes_map] = create_dbn(signaling_dbn_factory);
npers= dbn.nnodes_per_slice;
dbn_engine = jtree_dbn_inf_engine(dbn);

evidence=cell(npers,Json_signaling.DataInput.T);

[engine, ll] = enter_evidence(dbn_engine, evidence);
disp(ll);

% writing model variables in a .json file for all time slices.
keys =  keys(nodes_map);

% Create a table with the data and variable names
T = table();

for node_ndx = 1:npers
    keydnx = cellfun(@(x)isequal(x,node_ndx),values(nodes_map));
    node_name = keys(keydnx);
    node_values = {};
    node_values(end+1,:) = {node_name,node_name,node_name}
    for slice = 1:Json_signaling.DataInput.T
        marg = marginal_nodes(engine,node_ndx,slice);
        node_values(end+1,:) = {marg.mu, marg.Sigma, sqrt(marg.Sigma)};
    end
    T = [T node_values];
end 

% Write data to text file
writetable(T, '../results/models/signaling.csv');
