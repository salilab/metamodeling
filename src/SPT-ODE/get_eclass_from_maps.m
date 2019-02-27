function [eclass]= get_eclass_from_maps(eclass_map, nodes_map, CPD_names)
% eclass_map - a map from node names to CPD names
% nodes_map - a map from node names to node indexes in DBN
% cpds - a cell array of CPD names
%
% Return an eclass for make_dbn, indexed by node index in DBN graph, with CPD
% equivalence class in entry #node_index
    n= length(nodes_map);
    eclass= zeros(n,1);
    cpd_to_id_map= containers.Map(CPD_names, 1:numel(CPD_names));
    for node_name_cell= keys(nodes_map)
        node_name= node_name_cell{1}; 
        node_index= nodes_map(node_name);
        cpd_name= eclass_map(node_name);
        cpd_id= cpd_to_id_map(cpd_name);
        eclass(node_index)= cpd_id;
    end
    
end