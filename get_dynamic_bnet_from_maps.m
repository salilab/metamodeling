function [dbn]= get_dynamic_bnet_from_maps(nodes_names, edges_inter, edges_intra, eclass1_map, eclass2_map, cpd_names)
    node_names=  horzcat(strcat('ODE.', {'h', 'G', 'G_minus_h', 'I', 'Gref', 'Gexp', 'Iexp'}), 'Reference.I'); % BARAK comment: removed alpha, beta (turned to weights) + added intermediate (G-h)
    n= length(node_names);
    ns = ones(1, n);% all cts nodes are scalar
    % Intra - in one time slice
    edges_intra1= strcat('ODE.', {'h', 'G_minus_h'; 'G', 'G_minus_h'; 'G', 'Gref'; 'Gref', 'Gexp'}); % BARAK comment: changed in accordance with change in nodes list
    edges_intra2= {'ODE.I','Reference.I'; 'Reference.I','ODE.Iexp'};
    edges_intra= [edges_intra1; edges_intra2];
    % Inter - between time slices
    edges_inter= strcat('ODE.', { 'G', 'G'; 'G_minus_h', 'I'; 'I', 'I' }); % BARAK comment: switched G->I and h->I to (G-h)->I
    [intra, inter, nodes_map, ~]= get_valid_nodes_graph(node_names, edges_intra, edges_inter);
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('ODE.G')= 'ODE.G.inter';
    eclass2_map('ODE.I')= 'ODE.I.inter';   
    cpd_names=  unique( [values(eclass1_map), values(eclass2_map)] );
    
    eclass1= get_eclass_from_maps(eclass1_map, nodes_map, cpd_names);
    eclass2= get_eclass_from_maps(eclass2_map, nodes_map, cpd_names);  

    % make the dbn
    dbn = mk_dbn(intra, inter, ns, ...
        'discrete', [], ...
        'eclass1', eclass1, ...
        'eclass2', eclass2);

    cpdfs={};
    cpdfs{end+1}= CPDFactory ('Gaussian_CPD','ODE.h', 1,{'mean', 3.0, 'cov', 5.0});
    for i=1:numel(cpdfs)
        cpdf= cpdfs{i};
        create_and_associate_cpd(cpdf, dbn, nodes_map)
    end
end

