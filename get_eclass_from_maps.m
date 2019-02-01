function [eclass]= get_eclass_from_maps(eclass_map, nodes_map)
% eclass_map - a map from node names to CPD indexes
% nodes_map - a map from node names to node indexes in DBN
%
% Return an eclass for make_dbn, indexed by node index in DBN graph, with CPD
% equivalence class in entry #node_index
    n= length(nodes_map);
    eclass= zeros(n,1);
    % TODO: asserting that the keys of eclass_map and nodes_map are
    % identical
    for node_name_cell= keys(nodes_map)
        node_name= node_name_cell{1}; % TODO: ugly hack - iterate properly
        node_index= nodes_map(node_name);
        eclass(node_index)= eclass_map(node_name);
    end
    
end