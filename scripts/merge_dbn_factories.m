% Merge a list of DBNFactory classes, outputed by e.g.
% make_ode_dbn_descriptor()
function [merged_dbn_factory]= merge_dbn_factories(varargin)
   if nargin==0
       error('Cannot merge zero dbn descriptors')
   end
   merged_dbn_factory= copy(varargin{1});
   for i=2:nargin
       merged_dbn_factory.merge(varargin{i}, false);
   end
end