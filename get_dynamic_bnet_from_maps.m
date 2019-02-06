function [dbn]= get_dynamic_bnet_from_maps(node_names, edges_intra, edges_inter, eclass1_map, eclass2_map, CPD_factories)
    n= length(node_names);
    ns = ones(1, n);% all cts nodes are scalar
    % Create DAG
    [intra, inter, nodes_map, ~]= get_valid_nodes_graph(node_names, edges_intra, edges_inter);
    cpd_names=  unique( [values(eclass1_map), values(eclass2_map)] );    
    eclass1= get_eclass_from_maps(eclass1_map, nodes_map, cpd_names);
    eclass2= get_eclass_from_maps(eclass2_map, nodes_map, cpd_names);  
    % Make the dbn:
    dbn = mk_dbn(intra, inter, ns, ...
        'discrete', [], ...
        'eclass1', eclass1, ...
        'eclass2', eclass2);
    % Generate the CPDs:
    for i=1:numel(CPD_factories)
        CPD_factory= CPD_factories{i};
        disp(CPD_factory);
        create_and_associate_cpd(CPD_factory, dbn, nodes_map)
    end
end

