% Make a DBN for the exocytosis model with the following variables
%
% Time-dependent variables
%
% Coupling variables
% Gex.C, Gcell.C, Scell.C
%
% Observed variables
%
% Time-invariant variables
%
% Parameters
%
% To generate a conditional gaussian model
%
% All variables display gaussian distributions.

% To generate a conditional gaussian model
function [dbn_factory]= make_exocytosis_dbn_nocouple(G_mean_exocytosis, kt_mean_exocytosis, Npatch_mean_exocytosis,...
                                            Nvesicle_mean_exocytosis, Ninsulin_mean_exocytosis, Rcell_mean_exocytosis,...
                                            Dvesicle_mean_exocytosis, S_mean_exocytosis, G_cov_exocytosis,...
                                            kt_cov_exocytosis, Npatch_cov_exocytosis, Nvesicle_cov_exocytosis, ...
                                            Ninsulin_cov_exocytosis, Rcell_cov_exocytosis,Dvesicle_cov_exocytosis, ...
                                            S_cov_exocytosis,...
                                            cov_scale_exocytosis,...
                                            kG_exocytosis, alpha_exocytosis, ...
                                            kp_exocytosis,beta_exocytosis, kinsulin_exocytosis, ...
                                            kD_exocytosis, kR_exocytosis);

    node_names=  {'Gcell.C','G.exocytosis', 'kt.exocytosis', 'Npatch.exocytosis','Nvesicle.exocytosis',...
                  'Dvesicle.exocytosis', 'Ninsulin.exocytosis','Rcell.exocytosis', ...
                  'S.exocytosis', 'Scell.C', 'Scell.obs'};
    
%     node_names=  {'Gcell.C','G.exocytosis', 'kt.exocytosis', 'Npatch.exocytosis','Nvesicle.exocytosis',...
%                   'Dvesicle.exocytosis', 'Ninsulin.exocytosis','Rcell.exocytosis', 'S.exocytosis', 'Scell.C'};
    
              % Intra - in one time slice
    edges_intra= { 'Gcell.C','G.exocytosis';
                    'Scell.C','Scell.obs';
                    'S.exocytosis','kt.exocytosis';
                    'Npatch.exocytosis', 'S.exocytosis';...
                   'S.exocytosis','Nvesicle.exocytosis'; 
                   'Dvesicle.exocytosis', 'S.exocytosis'; 
                   'Ninsulin.exocytosis','S.exocytosis';...
                   'Rcell.exocytosis','S.exocytosis';
                   'S.exocytosis','Scell.C'};
               
%     edges_intra= { 'Gcell.C','G.exocytosis';
%                     'S.exocytosis','kt.exocytosis';
%                     'Npatch.exocytosis', 'S.exocytosis';...
%                    'S.exocytosis','Nvesicle.exocytosis'; 
%                    'Dvesicle.exocytosis', 'S.exocytosis'; 
%                    'Ninsulin.exocytosis','S.exocytosis';...
%                    'Rcell.exocytosis','S.exocytosis';
%                    'S.exocytosis','Scell.C'};
    
    % Inter - between time slices
    edges_inter= { 'G.exocytosis', 'G.exocytosis'; 
                    'G.exocytosis', 'S.exocytosis'; 
                    'S.exocytosis', 'S.exocytosis';...
                   'kt.exocytosis', 'kt.exocytosis';
                   'Nvesicle.exocytosis', 'Nvesicle.exocytosis'; 
                   'Npatch.exocytosis', 'Npatch.exocytosis';...
                   'Ninsulin.exocytosis', 'Ninsulin.exocytosis'};
    
    % 'Equivalence classes' specify how the template is initiated and rolled
    % Specify which CPD is associates with each node in either time
    % slice 1 (eclass1) or in slice 2 onwards (eclass2)
    eclass1_map= containers.Map();
    eclass2_map= containers.Map();
    for i=1:numel(node_names)
        node_name= node_names{i};
        cpd_name= [ node_name '.intra' ];
        eclass1_map(node_name) = cpd_name;
        eclass2_map(node_name) = cpd_name; % default - to be changed for some special cases
    end
    eclass2_map('G.exocytosis')= 'G.exocytosis.inter';
    eclass2_map('S.exocytosis')= 'S.exocytosis.inter';
    eclass2_map('kt.exocytosis')= 'kt.exocytosis.inter'; 
    eclass2_map('Nvesicle.exocytosis')= 'Nvesicle.exocytosis.inter';   
    eclass2_map('Npatch.exocytosis')= 'Npatch.exocytosis.inter';   
    eclass2_map('Ninsulin.exocytosis')= 'Ninsulin.exocytosis.inter';    
    
    
    
    % elcass1 (time-slice 0 or all parents are in the same time slice)
    % When using clamp, the root node is clamped to the N(0,I) distribution, so that we will not update these parameters during learninGex. 
    CPDFactories= {};
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Gcell.C', 0, ...
        {'mean', G_mean_exocytosis, 'cov', G_cov_exocytosis}); % Gcell.C
    
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'G.exocytosis', 0, ...
        {'mean', 0.0, 'cov', G_cov_exocytosis*cov_scale_exocytosis, 'weights', 1.0}); % G = Gcell.C
    
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'kt.exocytosis', 0, ...
        {'mean', kt_mean_exocytosis/2, 'cov', kt_cov_exocytosis, 'weights', alpha_exocytosis/2}); % k
    
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Npatch.exocytosis', 0, ...
        {'mean', Npatch_mean_exocytosis, 'cov', Npatch_cov_exocytosis}); % Npatch
    
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Nvesicle.exocytosis', 0, ...
        {'mean', 0.0, 'cov', Nvesicle_cov_exocytosis*cov_scale_exocytosis, 'weights', beta_exocytosis}); % Nvesicle
    
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Ninsulin.exocytosis', 0, ...
        {'mean', Ninsulin_mean_exocytosis, 'cov', Ninsulin_cov_exocytosis}); % Ninsulin
    
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Rcell.exocytosis', 0, ...
        {'mean', Rcell_mean_exocytosis, 'cov', Rcell_cov_exocytosis}); % Rcell
    
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'Dvesicle.exocytosis', 0, ...
        {'mean', Dvesicle_mean_exocytosis, 'cov', Dvesicle_cov_exocytosis}); % Dvesicle
    
    weights_S0_map_T0= containers.Map(); % parents in slice t
    weights_S0_map_T1= containers.Map(); % parents in slice t+1
    weights_S0_map_T0('Npatch.exocytosis')= 0.0;
    weights_S0_map_T0('Ninsulin.exocytosis')= 0.0;
    weights_S0_map_T0('Dvesicle.exocytosis')= 0.0;
    weights_S0_map_T0('Nvesicle.exocytosis')= 0.0;
    weights_S0_map_T0('Rcell.exocytosis')= 0.0;
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'S.exocytosis', 0, ...
        {'mean', S_mean_exocytosis, 'cov', S_cov_exocytosis}, ...
        weights_S0_map_T0, weights_S0_map_T1); % S
    
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Scell.C', 0, ...
        { 'mean', 0.0,'cov', S_cov_exocytosis*cov_scale_exocytosis, 'weights', 1.0} ); % Scell.C = S
       
    CPDFactories{end+1} = ...
        CPDFactory('Gaussian_CPD', 'Scell.obs', 0, ...
        { 'mean', 0.0,'cov', S_cov_exocytosis*cov_scale_exocytosis, 'weights', 1.0} ); % Scell.obs = Scell.C
    
    
    % eclass2 (time-slice t+1 with parents in the previous time slice)
    weights_k1_map_T0= containers.Map(); 
    weights_k1_map_T1= containers.Map(); 
    weights_k1_map_T0('kt.exocytosis')= 0.0; 
    weights_k1_map_T1('S.exocytosis')= alpha_exocytosis/2;
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'kt.exocytosis', 1, ...
        {'mean', kt_mean_exocytosis/2, 'cov', kt_cov_exocytosis*cov_scale_exocytosis}, ...
        weights_k1_map_T0, weights_k1_map_T1); % k(t+1) = 1.0 * k(t) 
    
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'Npatch.exocytosis', 1, ...
        {'mean', 0.0, 'cov', Npatch_cov_exocytosis*cov_scale_exocytosis,  'weights', 1.0} ); % Npatch(t+1) = 1.0 * Npatch(t) 
    
    weights_Nvesicle1_map_T0= containers.Map(); 
    weights_Nvesicle1_map_T1= containers.Map(); 
    weights_Nvesicle1_map_T0('Nvesicle.exocytosis')= 0.0; 
    weights_Nvesicle1_map_T1('S.exocytosis')= beta_exocytosis;
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'Nvesicle.exocytosis', 1, ...
        {'mean', 0.0, 'cov', Nvesicle_cov_exocytosis*cov_scale_exocytosis}, ...
        weights_Nvesicle1_map_T0, weights_Nvesicle1_map_T1); % k(t+1) = 1.0 * k(t)     
           
    CPDFactories{end+1} = ...
       CPDFactory('Gaussian_CPD', 'Ninsulin.exocytosis', 1, ...
        {'mean', 0.0, 'cov', Ninsulin_cov_exocytosis*cov_scale_exocytosis,  'weights', 1.0} ); % Npatch(t+1) = 1.0 * Npatch(t) 
    
    weights_G1_map_T0= containers.Map(); 
    weights_G1_map_T1= containers.Map(); 
    weights_G1_map_T0('G.exocytosis')= 0.0; 
    weights_G1_map_T1('Gcell.C')= 1.0;
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'G.exocytosis', 1, ...
        {'mean', 0.0, 'cov', G_cov_exocytosis*cov_scale_exocytosis}, ...
        weights_G1_map_T0, weights_G1_map_T1); % G(t+1) = 1.0 * Gcell.C + 0.0 * G(t)
    
    weights_S1_map_T0= containers.Map(); 
    weights_S1_map_T1= containers.Map(); 
    weights_S1_map_T0('S.exocytosis')= 0; 
    weights_S1_map_T0('G.exocytosis')= kG_exocytosis; 
    weights_S1_map_T1('Npatch.exocytosis')= kp_exocytosis; 
    weights_S1_map_T1('Ninsulin.exocytosis')= kinsulin_exocytosis; 
    weights_S1_map_T1('Dvesicle.exocytosis')= kD_exocytosis;
    weights_S1_map_T1('Rcell.exocytosis')= kR_exocytosis;
    CPDFactories{end+1} = ...         
        CPDFactory('Gaussian_CPD', 'S.exocytosis', 1, ...
        {'mean', 0.0, 'cov', S_cov_exocytosis*cov_scale_exocytosis}, ...
        weights_S1_map_T0, weights_S1_map_T1); % S(t+1) = 0*S(t) + kG_exocytosis*G(t) + kt_w_S_exocytosis * k(t+1) 
                                               % kp_exocytosis*Npatch(t+1) + Nvesicle_w_S_exocytosis * Nvesicle(t+1) + kinsulin_exocytosis * Ninsulin(t+1)
                                               % + kD_exocytosis*Dvesicle + kR_exocytosis*Rcell 

    % Final DBN factory
    dbn_factory= DBNFactory( ...
        node_names, edges_intra, edges_inter, ...
        eclass1_map, eclass2_map, CPDFactories);
end     


