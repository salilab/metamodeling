function [intra, inter, nodes_map, reverse_nodes_map]= get_valid_nodes_graph(node_names, edges_intra, edges_inter)
% Returns inter and intra adjacency matrices based on:
% node_names - a cell array of the nodes names
% edges_intra - intra-edges between nodes in the same time slice
% edges_inter - inter-edges between nodes in consecutive time slices
%
% The computed nodes indexes in inter and intra are topologically sorted within
% the same time slice, such that parents in intra always appear before
% their children
%
% Return:
% intra - intra-edge adjacency matrix such that parents always come before
%         children
% inter - inter-edge adjacency matrix in the same order as intra
% nodes_map - a map from node names to node indexes in intra/inter
% reverse_nodes_map a map from node indexes in intra/inter to node names

    % Create intra and inter adjacency matrixes
    n= length(node_names);
    nodes_map= containers.Map(node_names, 1:n);
    intra= zeros(n);
    for row=1:size(edges_intra,1)
        from= edges_intra{row, 1};
        to= edges_intra{row, 2};
        intra(nodes_map(from), nodes_map(to))= 1;
    end
    inter= zeros(n);
    for row=1:size(edges_inter, 1)
        from= edges_inter{row, 1};
        to= edges_inter{row, 2};
        inter(nodes_map(from), nodes_map(to))= 1;
    end
    % Resort all in topological order
    order= topological_sort(intra);
    intra= intra(order, order);
    inter= inter(order, order);
    nodes_map= containers.Map(node_names(order), 1:n);
    reverse_nodes_map= get_reverse_nodes_map(nodes_map)
end
