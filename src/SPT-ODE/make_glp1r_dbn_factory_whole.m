% Make a DBN for the GLP1R model with the following variables
%
% Time-dependent variables
%  -> GLP1.GLP1R(t)  ->  GLP1.GLP1R(t+1) ->
%  -> GLP1R.GLP1R(t)  ->  GLP1R.GLP1R(t+1) ->
%  -> Galpha.GLP1R(t)  ->  Galpha.GLP1R(t+1) ->
%  -> cAMP.GLP1R(t)  ->  cAMP.GLP1R(t+1) ->
%
% Reference variables
% cAMP.ref(t), cAMP.ref(t+1) 
%
% Observed variables
% cAMP.obs(t), cAMP.obs(t+1)
%
% Time-invariant variables
% 
%
% Parameters
%
% To generate a conditional gaussian model
function [dbn_factory]= make_glp1r_dbn_factory_whole(GLP1_GLP1R,GLP1R_GLP1R, Galpha_GLP1R, cAMP_GLP1R)
    node_names=  {'GLP1.GLP1R','GLP1R.GLP1R', 'Galpha.GLP1R','cAMP.GLP1R','cAMP.ref','cAMP.obs'}; 
    n= length(node_names);
    % Intra - in one time slice
    edges_intra= {'GLP1.GLP1R','GLP1R.GLP1R';'GLP1R.GLP1R','Galpha.GLP1R';'cAMP.GLP1R','cAMP.ref';...
        'cAMP.ref','cAMP.obs'};
    % Inter - between time slices
    edges_inter= {'GLP1.GLP1R', 'GLP1.GLP1R'; 'GLP1R.GLP1R','GLP1R.GLP1R';'Galpha.GLP1R','Galpha.GLP1R';...
        'Galpha.GLP1R','cAMP.GLP1R';'cAMP.GLP1R','cAMP.GLP1R'}; 
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('GLP1.GLP1R')= 'GLP1.GLP1R.inter';
    eclass2_map('GLP1R.GLP1R')= 'GLP1R.GLP1R.inter';
    eclass2_map('Galpha.GLP1R')= 'Galpha.GLP1R.inter';   
    eclass2_map('cAMP.GLP1R')= 'cAMP.GLP1R.inter';       
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    CPDFactories= {};
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1.GLP1R', 0, ...
        {'mean', GLP1_GLP1R, 'cov', 0.001} ); 
    
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.GLP1R', 0, ...
        {'mean', GLP1R_GLP1R, 'cov', 0.001, 'weights',1.0});

    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'Galpha.GLP1R', 0, ...
        {'mean', Galpha_GLP1R, 'cov', 0.001,'weights',1.0});
    % works only when there is ATP !!!!!!!!TODO!!!!!!!!
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cAMP.GLP1R', 0, ...
        {'mean', cAMP_GLP1R, 'cov', 0.001});
    
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'cAMP.ref', 0, ...
        {'mean', 0.0, 'cov', 0.001,  'weights', 1.0} );  
    
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'cAMP.obs', 0, ...
        {'mean', 0.0, 'cov', 0.001,  'weights', 1.0} );  

    % eclass2 (time-slice t+1 with parents in the previous time slice)
    
    % CPD for G(t+1)
    weights_GLP11_T0= containers.Map(); % parents in slice t
    weights_GLP11_T1= containers.Map(); % parents in slice t+1
    weights_GLP11_T0('GLP1.GLP1R')= 1.0;
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'GLP1.GLP1R', 1, ...
        {'mean', 0.0, 'cov', 0.001}, weights_GLP11_T0, weights_GLP11_T1);  

    weights_GLP1R1_T0= containers.Map(); % parents in slice t
    weights_GLP1R1_T1= containers.Map(); % parents in slice t+1
    weights_GLP1R1_T0('GLP1R.GLP1R')= 0.0;
    weights_GLP1R1_T1('GLP1.GLP1R')= 1.0;
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'GLP1R.GLP1R', 1, ...
        {'mean', 0.0, 'cov', 0.001}, weights_GLP1R1_T0, weights_GLP1R1_T1);  

    weights_Galpha1_T0= containers.Map(); % parents in slice t
    weights_Galpha1_T1= containers.Map(); % parents in slice t+1
    weights_Galpha1_T0('Galpha.GLP1R')= 0.0;
    weights_Galpha1_T1('GLP1R.GLP1R')= 1.0;
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'Galpha.GLP1R', 1, ...
        {'mean', 0.0, 'cov', 0.001}, weights_Galpha1_T0, weights_Galpha1_T1);  
    
    weights_cAMP1_T0= containers.Map(); % parents in slice t
    weights_cAMP1_T1= containers.Map(); % parents in slice t+1
    weights_cAMP1_T0('Galpha.GLP1R')= 0.5;
    weights_cAMP1_T0('cAMP.GLP1R')= 0.5;
    CPDFactories{end+1}=  ...
        CPDFactory('Gaussian_CPD', 'cAMP.GLP1R', 1, ...
        {'mean', 0.0, 'cov', 0.001}, weights_cAMP1_T0, weights_cAMP1_T1);

    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
  end

