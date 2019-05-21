
function [meta_dbn, nodes_map]=make_meta_bnet2(G_Meal, I_Meal, alpha_Meal,beta_Meal, K_Meal, dt_Meal, DG_Meal, Gb_Meal, Sb_Meal, ...
    G_SPT, S_SPT, k_SPT, Nisg_SPT, Npatch_SPT, Ninsulin_SPT)

    % make ODE and BD models
    [meal_dbn_factory]= ...
        make_meal_dbn_factory_eq2(G_Meal, I_Meal, alpha_Meal,beta_Meal, K_Meal, dt_Meal, DG_Meal, Gb_Meal, Sb_Meal);
    [spt_dbn_factory]= ...
        make_spt_dbn_factory(G_SPT, S_SPT, k_SPT, Nisg_SPT, Npatch_SPT, Ninsulin_SPT);
    
    meta_dbn_factory= ...
        merge_dbn_factories(meal_dbn_factory, spt_dbn_factory);
                          
    weights_S0_map_T0= containers.Map(); % parents in slice t
    weights_S0_map_T1= containers.Map(); % parents in slice t+1
    weights_S0_map_T0('S.Meal')= 0.5;
    weights_S0_map_T0('S.SPT')= 0.5;
    CPDFactory_S0 = ...
        CPDFactory('Gaussian_CPD', 'S.ref', 0, ...
        {'mean', 0.0, 'cov', 1.0}, ...
        weights_S0_map_T0, ...
        weights_S0_map_T1);
 
    add_CPD_factories(meta_dbn_factory, {CPDFactory_S0}, false);
    %add_CPD_factories(meta_dbn_factory, {CPDFactory_cAMP0}, false);

    [meta_dbn, ~, ~, nodes_map] = create_dbn(meta_dbn_factory);
    
end
