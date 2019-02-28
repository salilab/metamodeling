classdef CPDFactory
    % Factory for CPDs that uses string labels instead of numeral ids (so
    % more flexible) and supports some (not all yet) CPD types and params
    % Currently - only supporting Gaussin_CPD with continuous parents
    
    properties %(GetAccess=private)
        type_string
        rep_node_name
        rep_node_time_slice
        non_weights_params
        weights_map_T0
        weights_map_T1
    end   
    
    
    methods(Access= public)
        
        function obj=CPDFactory(type_string, ...
                rep_node_name, ...
                rep_node_time_slice, ...
                non_weights_params, ...
                weights_map_T0, weights_map_T1)
            % type_string - for now, supports just 'Gaussian_CPD'
            % rep_node_name - a representative of self node name in the
            %   final dbn
            % rep_node_time_slice - 0 or 1, to indicate time slice t or t+1
            % non_weights_params - a cell array with all non-weights params
            %   to CPD contructure in dbn apckage
            % weights_map_T0/1 - map from continuous parent names to their
            %   weights in time slice 0 and 1, respectively
            obj.type_string= type_string;
            obj.rep_node_name= rep_node_name;
            obj.non_weights_params= non_weights_params;
            assert(rep_node_time_slice==0 || rep_node_time_slice==1)
            obj.rep_node_time_slice= rep_node_time_slice;
            if nargin>=5
                if any(strcmp('weights', obj.non_weights_params))
                    error('Cannot have weights specified both in non_weights_params and in weights_map_T0/1')
                end
                obj.weights_map_T0= weights_map_T0;
                obj.weights_map_T1= weights_map_T1;
            end
        end  
       
        function [cpd, dbn]=create_and_associate_cpd(obj, dbn, nodes_map)
            % Create a new cpd according to factory parameters and
            % associate it with dynamic bayesian network dbn using
            % nodes_map to map node names to node indexes in dbn.dag
            
            % Check if any discrete parents
            rep_node_index= nodes_map(obj.rep_node_name);
            n= dbn.nnodes_per_slice;
            if obj.rep_node_time_slice==1
                rep_node_index= rep_node_index + n;
            end
            ps= parents(dbn.dag, rep_node_index);
            dps = myintersect(ps, dbn.dnodes);
            if ~isempty(dps)
                error('CPDFactory does not yet support discrete parents') % if we want to support, will need to reorder mean and cov matrices if specified by order of parents
            end            
            if strcmp(obj.type_string, 'Gaussian_CPD')
                cpd= create_gaussian_cpd(obj, dbn, nodes_map);
            else
                error(printf('CPDFactory does not yet support CPDs of type %s\n', obj.type_string));
            end
            % Associate cpd with dbn
            if obj.rep_node_time_slice==0
                cpd_index= dbn.eclass1(rep_node_index);
            else
                cpd_index= dbn.eclass2(rep_node_index - n);
            end
            dbn.CPD{cpd_index}= cpd;            
        end
        
    end
    
    
    methods (Access= protected)
        
        function cpd=create_gaussian_cpd(obj, dbn, nodes_map)
            % Create a new Gaussian cpd according to factory parameters
            % based on topology of dynamic bayesian network dbn
            assert(strcmp(obj.type_string, 'Gaussian_CPD'))
            % Get weights in right order
            n= dbn.nnodes_per_slice;
            assert(n == nodes_map.Count);
            rep_node_index= nodes_map(obj.rep_node_name);
            if obj.rep_node_time_slice==1
                rep_node_index= rep_node_index + n;
            end
            ps= parents(dbn.dag, rep_node_index);
            reverse_nodes_map= get_reverse_nodes_map(nodes_map);
            if any(strcmp('weights', obj.non_weights_params))
                cpd= gaussian_CPD(dbn, ...
                    rep_node_index, ...
                    obj.non_weights_params{:});
            else
                weights= zeros(length(ps),1);
                for i=1:length(ps)
                    parent_index= ps(i);                    
                    if (parent_index <= n) % parent in slice t
                        parent_name= reverse_nodes_map(parent_index);
                        weights(i)= obj.weights_map_T0(parent_name);
                    else % parent in slice t+1
                        parent_name= reverse_nodes_map(parent_index);
                        weights(i)= obj.weights_map_T1(parent_name);
                    end
                end
                % Create CPD
                cpd= gaussian_CPD(dbn, ...
                    rep_node_index, ...
                    obj.non_weights_params{:}, ...
                    'weights', weights);            
            end
        end
    end
end

