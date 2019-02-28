classdef DBNFactory < matlab.mixin.Copyable
    % Factory for DBNs that uses string labels instead of numeral ids (so
    % more flexible) and supports some (not all yet) CPD types and params
    % Currently - only supporting Gaussin_CPD with continuous parents
    
    properties %(GetAccess=private)
        node_names
        edges_intra
        edges_inter
        eclass1_map
        eclass2_map
        CPD_factories_map
    end
    
    
    methods(Access= public)
       
        
        function obj=DBNFactory( ...
                node_names, ...
                edges_intra, ...
                edges_inter, ...
                eclass1_map, ...
                eclass2_map, ...
                CPD_factories)
            % node_names - cell array of node names in the DAG
            % edges_intra - cell array of edges within the same time
            %     slice in the DBN. Each edge is a 1x2 cell array of node
            %     names
            % edges_inter - cell array of edges from time slice t to time
            %     slice t+1
            % eclass1_map - map from node names to the CPD equivalence
            %     class label of each node at time slice 0. Nodes that map
            %     to the same label use the same CPD
            % eclass2_map - map from node names to the CPD equivalence
            %     class label of each node at time slice 1. Nodes that map
            %     to the same label use the same CPD
            if nargin<1
                node_names= {};
            end
            if nargin<2
                edges_intra= {};
            end
            if nargin<3
                edges_inter= {};
            end
            if nargin<4
                eclass1_map= containers.Map();
                eclass2_map= containers.Map();
                CPD_factories= {};
            end
            obj.node_names= node_names;
            obj.edges_intra= edges_intra;
            obj.edges_inter= edges_inter;
            obj.eclass1_map= eclass1_map;
            obj.eclass2_map= eclass2_map;
            obj.CPD_factories_map= containers.Map();
            add_CPD_factories(obj, CPD_factories);
        end % function
       
        function [dbn, intra, inter, nodes_map]= create_dbn(obj)
            % Create a new dbn according to factory parameters
            n= length(obj.node_names);
            ns = ones(1, n);% all cts nodes are scalar
            % Create DAG
            [intra, inter, nodes_map, ~]= get_valid_nodes_graph(...
                obj.node_names, ...
                obj.edges_intra, ...
                obj.edges_inter);
            cpd_names=  unique(...
                [ values(obj.eclass1_map), ...
                  values(obj.eclass2_map) ] );
            eclass1= get_eclass_from_maps(...
                obj.eclass1_map, nodes_map, cpd_names);
            eclass2= get_eclass_from_maps(...
                obj.eclass2_map, nodes_map, cpd_names);
            % Make the dbn:
            dbn = mk_dbn(intra, inter, ns, ...
                'discrete', [], ...
                'eclass1', eclass1, ...
                'eclass2', eclass2);
            % Generate the CPDs:
            CPD_factories= values(obj.CPD_factories_map);
            for i=1:numel(CPD_factories)
                CPD_factory= CPD_factories{i};
                disp(CPD_factory);
                [~,dbn]=create_and_associate_cpd(CPD_factory, dbn, nodes_map);
            end
        end % function 
        
        
        function []= add_CPD_factories(obj, CPD_factories, error_if_duplicates)
            % Add cpd factories in cell array CPD_factories to the DBN
            % factory. 
            %
            % If error_if_duplicates is true, raise error if a duplicate
            % CPD is found for the same equivalence class of the
            % representative nodes at the same time slice.
            if nargin<3
                error_if_duplicates= true;
            end
            for i=1:numel(CPD_factories)
                CPD_factory= CPD_factories{i};
                if CPD_factory.rep_node_time_slice==0
                    eclass_map= obj.eclass1_map;
                else
                    eclass_map= obj.eclass2_map;
                end
                if ~isKey(eclass_map, CPD_factory.rep_node_name)
                    msg= sprintf( ...
                        "No CPD is associated with the equivalence class of representative node %s at time slice %d", ...
                        CPD_factory.rep_node_name, ...
                        CPD_factory.rep_node_time_slice);
                    error(msg);
                end
                CPD_tag= eclass_map(CPD_factory.rep_node_name);
                if isKey(obj.CPD_factories_map, CPD_tag)
                    msg= sprintf(...
                        ['Duplicate CPD exists for the equivalence class %s', ...
                         ' (representative node %s at time slice %d)'], ...
                         CPD_tag, ...
                         CPD_factory.rep_node_name, ...
                         CPD_factory.rep_node_time_slice);
                    if error_if_duplicates
                        error(msg);
                    else
                        warning(msg)
                    end
                end
                obj.CPD_factories_map(CPD_tag)= CPD_factory;
            end      
        end
        
        
        function []= merge(obj, other, error_if_duplicates)
            % Merge another CPD factory into this factory. The DBN graph is
            % merged based on node labels to identify nodes uniquely.
            
            % error_if_duplicates [default=true]
            %  If true, an error is raised if EClass1 map or EClass2map
            %  point the same node at the same time slice to a different
            %  CPD equivalence class, or if there are duplicates in the
            %  CPD factories of the same representative node at the same
            %  time slice.
            if nargin<3
                error_if_duplicates=1
            end
           if ~strcmp(class(obj), class(other))
               error("Cannot merge objects of class %s to %s", class(other), class(obj));
           end
           obj.node_names= union(obj.node_names, other.node_names);
           obj.edges_intra= [obj.edges_intra; other.edges_intra];
           obj.edges_inter= [obj.edges_inter; other.edges_inter];
           for key_cell=keys(other.eclass1_map)
               key= key_cell{1};
               value= other.eclass1_map(key);
               if isKey(obj.eclass1_map, key) 
                   if(obj.eclass1_map(key) ~= value)
                       msg= sprintf("EClass1 maps contain non-identical values for key %s", key);
                       if error_if_duplicates
                           error(msg);
                       else
                           warn(msg)
                       end
                   end
               else
                   obj.eclass1_map(key)= value;
               end
           end
           for key_cell=keys(other.eclass2_map)
               key= key_cell{1};
               value= other.eclass2_map(key);
               if isKey(obj.eclass2_map, key) 
                   if(obj.eclass2_map(key) ~= value)
                       msg= sprintf("EClass2 maps contain non-identical values for key %s", key);
                       if error_if_duplicates
                           error(msg);
                       else
                           warn(msg)
                       end
                   end
               else
                   obj.eclass2_map(key)= value;
               end
           end   
           obj.add_CPD_factories(...
               values(other.CPD_factories_map), ...
               error_if_duplicates);               
        end
        
    end % public methods
end

