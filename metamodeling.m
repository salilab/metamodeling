% Make a meta-modeling DBN for the ODE and SPT model

% Read in the experimental measurements
odem1 = importdata('ode_exp1_avr.dat');
Gm1 = odem1(:,2); % Gexp in measurement number 1, vector along time
Im1 = odem1(:,3); % Iexp in measurement number 1, vector along time

sptm1 = importdata('spt_obs1_avr.dat');
Go1 = sptm1(:,2); % Gexp in measurement number 1, vector along time
Io1 = sptm1(:,3); % Iexp in measurement number 1, vector along time

[meta_meta_bnet, nodes_map] = make_meta_bnet(Gm1, Im1, Go1, Io1);

function [meta_bnet, nodes_map]=make_meta_bnet(Gm, Im, Go, Io)

% make ODE and BD models
[ode_bnet, ode_nodes_map, ode_node_names, ode_edges_intra, ode_edges_inter, ode_ns, ode_eclass1, ode_eclass2]= make_ode_bnet(Gm, Im);
ode_n = 9;

[spt_bnet, spt_nodes_map, spt_node_names, spt_edges_intra, spt_edges_inter, spt_ns, spt_eclass1, spt_eclass2]= make_spt_bnet(Go, Io);
spt_n = 11;

% change names and get a valid nodes_graph

meta_node_names=unique([strrep(strcat(ode_node_names,'-ODE'),'Iref-ODE','Iref'),strrep(strcat(spt_node_names,'-SPT'),'Iref-SPT','Iref')],'stable');
disp(meta_node_names);

% Total number of nodes minus the shared node
meta_n = ode_n + spt_n -1;

% Total edges
meta_edges_intra = [strrep(strcat(ode_edges_intra,'-ODE'),'Iref-ODE','Iref'); strrep(strcat(spt_edges_intra,'-SPT'),'Iref-SPT','Iref')];
meta_edges_inter = [strrep(strcat(ode_edges_inter,'-ODE'),'Iref-ODE','Iref'); strrep(strcat(spt_edges_inter,'-SPT'),'Iref-SPT','Iref')];
disp(meta_edges_intra);

% get valid nodes graph for the meta-model
[intra, inter, nodes_map]= get_valid_nodes_graph(meta_node_names, meta_edges_intra, meta_edges_inter);
disp('check');
disp(nodes_map('Iref'));

% ???? TO CORRECT ECLASS AS WELL ????
meta_eclass1 = [ode_eclass1; spt_eclass1+max(ode_eclass1)-1];
meta_eclass2 = [ode_eclass2; spt_eclass2+max(ode_eclass2)-1];
%disp(meta_eclass1);

% make the new dbn
meta_bnet = mk_dbn(intra, inter, meta_n, ...
        'discrete', [], ...
        'eclass1', meta_eclass1, ...
        'eclass2', meta_eclass2);
    
% eclass1
weights_Iref= [0.5 0.5];
bnet.CPD{8} = gaussian_CPD(bnet, nodes_map('Iref')+n, 'mean', Im(n), 'cov', 5,  'weights', weights_Iref);

% eclass2
weights_Iref= [0.5 0.5];
bnet.CPD{10} = gaussian_CPD(bnet, nodes_map('Iref')+n, 'mean', Im(n), 'cov', 5,  'weights', weights_Iref);

end
