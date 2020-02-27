
function [meta_dbn, nodes_map]=make_meta_bnet5(Gex_Meal, Y_Meal, alpha_Meal,beta_Meal, gamma_Meal, k1_Meal, k2_Meal, K_Meal, dt_Meal_min, Gb_Meal, DGintake_Meal, Sb_Meal, I_Meal, ...
    Gin_SPT, S_SPT, k_SPT, Nisg_SPT, Npatch_SPT, Ninsulin_SPT, I_SPT, ...
    Gin_Network, PFK_Network, ATP_Network, GLP1_Network, GLP1R_Network, cAMP_Network, Ca_Network, S_Network, I_Network, ...
    GLP1a_GLP1R, cons_GLP1R,GLP1R_GLP1R,...
    Gculture_KEGG, Ex4_KEGG, NES_KEGG, Pathway_KEGG, ATP_KEGG)

    % make ODE and BD models
    [meal_dbn_factory]= ...
        make_meal_dbn_factory_eq(Gex_Meal, Y_Meal, alpha_Meal,beta_Meal, gamma_Meal, k1_Meal, k2_Meal, K_Meal, dt_Meal_min, Gb_Meal, DGintake_Meal, Sb_Meal, I_Meal);
    [spt_dbn_factory]= ...
        make_spt_dbn_factory(Gin_SPT, S_SPT, k_SPT, Nisg_SPT, Npatch_SPT, Ninsulin_SPT, I_SPT);
    [network_dbn_factory]= ...
        make_network_dbn_factory(Gin_Network, PFK_Network, ATP_Network, GLP1_Network, GLP1R_Network, cAMP_Network, Ca_Network, S_Network, I_Network);
    [glp1r_dbn_factory]= ...
        make_glp1r_dbn_factory(GLP1a_GLP1R, cons_GLP1R,GLP1R_GLP1R);
    [KEGG_dbn_factory]= ...
        make_KEGG_dbn_factory(Gculture_KEGG, Ex4_KEGG, NES_KEGG, Pathway_KEGG, ATP_KEGG);
    
    meta_dbn_factory= ...
        merge_dbn_factories(meal_dbn_factory, spt_dbn_factory, network_dbn_factory, glp1r_dbn_factory, KEGG_dbn_factory);
                          
    weights_Scell0_map_T0= containers.Map(); % parents in slice t
    weights_Scell0_map_T1= containers.Map(); % parents in slice t+1
    weights_Scell0_map_T0('S.SPT')= 1/2;
    weights_Scell0_map_T0('S.Network')= 1/2;
    CPDFactory_Scell0 = ...
        CPDFactory('Gaussian_CPD', 'Scell.ref', 0, ...
        {'mean', 0.0, 'cov', 1E-4}, ...
        weights_Scell0_map_T0, ...
        weights_Scell0_map_T1);
    
    weights_Gin0_map_T0= containers.Map(); % parents in slice t
    weights_Gin0_map_T1= containers.Map(); % parents in slice t+1
    weights_Gin0_map_T0('Gex.ref')= 0.5;
    CPDFactory_Gin0 = ...
        CPDFactory('Gaussian_CPD', 'Gin.ref', 0, ...
        {'mean', 0.0, 'cov', 1E-4}, ...
        weights_Gin0_map_T0, ...
        weights_Gin0_map_T1);
    
    weights_S0_map_T0= containers.Map(); % parents in slice t
    weights_S0_map_T1= containers.Map(); % parents in slice t+1
    INITIAL_K = K_Meal;
    weights_S0_map_T0('DGintake.Meal')= 0.0;
    weights_S0_map_T0('Scell.ref')= 2/3;
    weights_S0_map_T0('Y.Meal')= 0.0;
    CPDFactory_S0 =  ...
        CPDFactory('Gaussian_CPD', 'S.Meal', 0, ...
        {'mean', Sb_Meal/3, 'cov', 0.001}, ...
        weights_S0_map_T0, weights_S0_map_T1); % S = 1.0 * Scell.ref
  
    weights_S1_map_T0= containers.Map(); % parents in slice t
    weights_S1_map_T1= containers.Map(); % parents in slice t+1
    INITIAL_K = K_Meal;
    weights_S1_map_T1('DGintake.Meal')= INITIAL_K/3;
    weights_S1_map_T1('Scell.ref')= 2.0/3;
    weights_S1_map_T1('Y.Meal')= 1.0/3;
    weights_S1_map_T0('S.Meal')= 0.0;
    CPDFactory_S1=  ...
        CPDFactory('Gaussian_CPD', 'S.Meal', 1, ...
        {'mean', Sb_Meal, 'cov', 0.001}, ...
        weights_S1_map_T0, weights_S1_map_T1); % S(t+1) = 0.0 * S(t) + INITIAL_K * DGex(t+1) + 1.0 * Scell.ref(t+1) + 1.0 * Y(t+1)
    
    CPDFactory_GLP1R0=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.Network', 0, ...
        {'mean', GLP1R_Network/2, 'cov', 1E-12,   'weights', 0.5}); % GLP1R = 1.0 * GLP1R.ref
    
    weights_GLP1R1_T0= containers.Map(); % parents in slice t
    weights_GLP1R1_T1= containers.Map(); % parents in slice t+1
    weights_GLP1R1_T0('GLP1.Network')= 0.5/GLP1_Network;
    weights_GLP1R1_T0('GLP1R.Network')= 0.0;
    weights_GLP1R1_T1('GLP1R.ref')= 0.5;
    CPDFactory_GLP1R1=  ...
        CPDFactory('Gaussian_CPD', 'GLP1R.Network', 1, ...
        {'mean', 0.0, 'cov', 1E-12}, weights_GLP1R1_T0, weights_GLP1R1_T1); % GLP1R(t+1) = 0.0 * GLP1R(t) + 0.5 * GLP1R.ref(t) + 0.5/GLP1_Network * GLP1(t)
    
    weights_ATP1_T0= containers.Map(); % parents in slice t
    weights_ATP1_T1= containers.Map(); % parents in slice t+1
    INITIAL_PFK= PFK_Network;
    weights_ATP1_T0('Gin.Network')= INITIAL_PFK;
    weights_ATP1_T0('ATP.Network')= 1.65/3.3;
    weights_ATP1_T1('ATP.ref')= 0.5;
    CPDFactory_ATP=  ...
        CPDFactory('Gaussian_CPD', 'ATP.Network', 1, ...
        {'mean', 0.0, 'cov', 1E-12}, weights_ATP1_T0, weights_ATP1_T1);  % ATP(t+1) = 1.65/3.3 * ATP(t) + PFK_Network * Gin (t)
    
    add_CPD_factories(meta_dbn_factory, {CPDFactory_Scell0, CPDFactory_Gin0, CPDFactory_S0, CPDFactory_S1, CPDFactory_GLP1R0, CPDFactory_GLP1R1, CPDFactory_ATP}, false);
    %add_CPD_factories(meta_dbn_factory, {CPDFactory_cAMP0}, false);

    [meta_dbn, ~, ~, nodes_map] = create_dbn(meta_dbn_factory);
    
end