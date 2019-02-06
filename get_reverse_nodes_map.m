function reverse_nodes_map= get_reverse_nodes_map(nodes_map)
% Create a reverse map from node index to node name
%
% nodes_map - input map from node name to node index
    node_names= keys(nodes_map);
    reverse_nodes_map= containers.Map('KeyType','uint32','ValueType','any');
    for i= 1:length(node_names)
        node_name= node_names{i};
        node_index= nodes_map(node_name);
        reverse_nodes_map(node_index)= node_name;
    end
end