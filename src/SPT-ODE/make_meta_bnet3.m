
function [meta_dbn, nodes_map]=make_meta_bnet3(G_Network, PFK_Network,S_Network, ...
    GLP1_GLP1R, cAMP_GLP1R)

    % make ODE and BD models
    [network_dbn_factory]= ...
        make_network_dbn_factory(G_Network, PFK_Network,S_Network);
    [glp1r_dbn_factory]= ...
        make_glp1r_dbn_factory(GLP1_GLP1R, cAMP_GLP1R);
    
    meta_dbn_factory= ...
        merge_dbn_factories(network_dbn_factory, glp1r_dbn_factory);
                          
    
    weights_cAMP0_map_T0= containers.Map(); % parents in slice t
    weights_cAMP0_map_T1= containers.Map(); % parents in slice t+1
    weights_cAMP0_map_T0('cAMP.GLP1R')= 0.5;
    weights_cAMP0_map_T0('cAMP.Network')= 0.5;
    CPDFactory_cAMP0 = ...
        CPDFactory('Gaussian_CPD', 'cAMP.ref', 0, ...
        {'mean', 0.0, 'cov', 0.001}, ...
        weights_cAMP0_map_T0, ...
        weights_cAMP0_map_T1);
    
    add_CPD_factories(meta_dbn_factory, {CPDFactory_cAMP0}, false);
    %add_CPD_factories(meta_dbn_factory, {CPDFactory_cAMP0}, false);

    [meta_dbn, ~, ~, nodes_map] = create_dbn(meta_dbn_factory);
    
end
