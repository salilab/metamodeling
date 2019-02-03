% Make a meta-modeling DBN for the ODE and SPT model

% Read in the experimental measurements
odem1 = importdata('ode_exp1_avr.dat');
Gm1 = odem1(:,2); % Gexp in measurement number 1, vector along time
Im1 = odem1(:,3); % Iexp in measurement number 1, vector along time

sptm1 = importdata('spt_obs1_avr.dat');
Go1 = sptm1(:,2); % Gexp in measurement number 1, vector along time
Io1 = sptm1(:,3); % Iexp in measurement number 1, vector along time

%for time  = 1:length(Gm1)
%    [bnet, nodes_map] = make_meta_bnet(Gm1, Im1, Go1, Io1, time);
%end
[bnet, nodes_map] = make_meta_bnet(Gm1, Im1, Go1, Io1, 1);

% parameter learning
npers= bnet.nnodes_per_slice;
T = 20; % lengthhs of sequences to explore
disp(npers);
% Sample from the posterior:
%disp(length(bnet.intra))
%disp(length(bnet.inter))
%evidence=cell(n,T);
%evidence{2,1}=1;
%sample_seq=  cell2mat(sample_dbn(bnet, 'length', T,'evidence', evidence));
sample_seq=  cell2mat(sample_dbn(bnet, 'length', T));
nodes_order= cell2mat(values(nodes_map));

% Plot ODE.G and ODE.I
%figure()
%yyaxis left;
%plot(1:T, sample_seq(nodes_map('ODE.G'),:));
%yyaxis right;
%plot(1:T, sample_seq(nodes_map('ODE.I'),:));
%legend('ODE.G','ODE.I');  

% Plot the distribution of ODE.beta
disp('plot');
fprintf("Sampled time-series of length %d", T);
y = sample_seq(nodes_map('ODE.beta'),:);
nbins = 10;

[hts,ctrs] = hist(y, nbins);
h = bar(ctrs,hts,'hist');
set(h,'FaceColor',[0.8 0.8 0.8],'EdgeColor',[0 0 0]);
area = sum(hts) * (ctrs(2)-ctrs(1));
xx = linspace(-0.1,0.3);
hold on; 
plot(xx,area*normpdf(xx,mean(y),std(y)),'k-','LineWidth',2);
fprintf("Normal probability density function of beta");
disp(mean(y));
disp(std(y));
%f = ksdensity(y,xx);
%plot(xx,area*f,'g-')
legend('ODE.beta, posterior');
hold off;

function [bnet, nodes_map]=make_meta_bnet(Gm, Im, Go, Io, time)

    % make ODE and BD models
    for time  = 1:length(Gm)
        [ode_bnet, ode_nodes_map, ode_node_names, ode_edges_intra, ode_edges_inter, ode_ns, ode_eclass1_map, ode_eclass2_map]= make_ode_bnet(Gm, Im,time);
        [spt_bnet, spt_nodes_map, spt_node_names, spt_edges_intra, spt_edges_inter, spt_ns, spt_eclass1_map, spt_eclass2_map]= make_spt_bnet(Go, Io,time);
    end

    ode_n = 9;
    spt_n = 11;

    % change names and get a valid nodes_graph
    %meta_node_names=unique([strrep(strcat(ode_node_names,'-ODE'),'Iref-ODE','Iref'),strrep(strcat(spt_node_names,'-SPT'),'Iref-SPT','Iref')],'stable');
    meta_node_names=unique([ode_node_names,spt_node_names],'stable');
    %disp(meta_node_names);

    % Total number of nodes minus the shared node
    n = length(meta_node_names);
    ns = ones(1, n);% all cts nodes are scalar
    % Total edges
    %meta_edges_intra = [strrep(strcat(ode_edges_intra,'-ODE'),'Iref-ODE','Iref'); strrep(strcat(spt_edges_intra,'-SPT'),'Iref-SPT','Iref')];
    meta_edges_intra = [ode_edges_intra; spt_edges_intra];
    %meta_edges_inter = [strrep(strcat(ode_edges_inter,'-ODE'),'Iref-ODE','Iref'); strrep(strcat(spt_edges_inter,'-SPT'),'Iref-SPT','Iref')];
    meta_edges_inter = [ode_edges_inter; spt_edges_inter];
    disp(meta_edges_intra);

    % get valid nodes graph for the meta-model
    [meta_intra, meta_inter, nodes_map, reverse_nodes_map]= get_valid_nodes_graph(meta_node_names, meta_edges_intra, meta_edges_inter);

    disp('check');
    %disp(meta_reverse_nodes_map);
    disp(nodes_map('Reference.I'));
    %disp(meta_nodes_map);
    %disp(ode_eclass1);
    %disp(spt_eclass1([1:spt_nodes_map('Reference.I')-1 spt_nodes_map('Reference.I')+1:end]));
    % ???? TO CORRECT ECLASS AS WELL ????
    %meta_eclass1 = [ode_eclass1; spt_eclass1+max(ode_eclass1)-1];
    %meta_eclass2 = [ode_eclass2; spt_eclass2+max(ode_eclass2)-1];
    %disp(meta_eclass1);
    eclass1_map= containers.Map();
    eclass1_map('ODE.alpha')=1;
    eclass1_map('ODE.beta')=2;
    eclass1_map('ODE.h')= 3;
    eclass1_map('ODE.G')= 4;
    eclass1_map('ODE.Gref')= 5;
    eclass1_map('ODE.Gexp')= 6;
    eclass1_map('ODE.I')= 7;
    eclass1_map('ODE.Iexp')= 8;
    eclass1_map('SPT.lambda')= 9;
    eclass1_map('SPT.k')= 10;
    eclass1_map('SPT.Npatch')= 11;
    eclass1_map('SPT.Nisg')= 12;
    eclass1_map('SPT.Rpbc')= 13;
    eclass1_map('SPT.G')= 14;
    eclass1_map('SPT.Gref')= 15;
    eclass1_map('SPT.Gobs')= 16;
    eclass1_map('SPT.I')= 17;
    eclass1_map('Reference.I')= 18;
    eclass1_map('SPT.Iobs')= 19;

    eclass2_map= containers.Map();
    eclass2_map('ODE.alpha')= 1;
    eclass2_map('ODE.beta')= 2;
    eclass2_map('ODE.h')= 3;
    eclass2_map('ODE.G')= 20;
    eclass2_map('ODE.Gref')= 5;
    eclass2_map('ODE.Gexp')= 6;
    eclass2_map('ODE.I')= 21;
    eclass2_map('ODE.Iexp')= 8;
    eclass2_map('SPT.lambda')= 9;
    eclass2_map('SPT.k')= 10;
    eclass2_map('SPT.Npatch')= 11;
    eclass2_map('SPT.Nisg')= 12;
    eclass2_map('SPT.Rpbc')= 13;
    eclass2_map('SPT.G')= 22;
    eclass2_map('SPT.Gref')= 14;
    eclass2_map('SPT.Gobs')= 15;
    eclass2_map('SPT.I')= 23;
    eclass2_map('Reference.I')= 18;
    eclass2_map('SPT.Iobs')= 19;
    
    eclass1= get_eclass_from_maps(eclass1_map, nodes_map);
    eclass2= get_eclass_from_maps(eclass2_map, nodes_map);  
    
    % make the new dbn
    bnet = mk_dbn(meta_intra, meta_inter, ns, ...
        'discrete', [], ...
        'eclass1', eclass1, ...
        'eclass2', eclass2);
    
    % eclass1
    bnet.CPD{1} = gaussian_CPD(bnet, nodes_map('ODE.alpha'),   'mean', 0.05, 'cov', 0.001);
    bnet.CPD{2} = gaussian_CPD(bnet, nodes_map('ODE.beta'),   'mean', 0.11, 'cov', 0.002);
    bnet.CPD{3} = gaussian_CPD(bnet, nodes_map('ODE.h'),   'mean', 6.1, 'cov', 0.1);
    bnet.CPD{4} = gaussian_CPD(bnet, nodes_map('ODE.G'),   'mean', Gm(1), 'cov', 0.1);
    bnet.CPD{5} = gaussian_CPD(bnet, nodes_map('ODE.Gref'),   'mean', Gm(1), 'cov', 0.1,'weights', 1);
    bnet.CPD{6} = gaussian_CPD(bnet, nodes_map('ODE.Gexp'), 'mean', Gm(1), 'cov', 0.1, 'weights', 1);
    bnet.CPD{7} = gaussian_CPD(bnet, nodes_map('ODE.I'),   'mean', Im(1), 'cov', 5);
    bnet.CPD{8} = gaussian_CPD(bnet, nodes_map('ODE.Iexp'), 'mean', Im(1), 'cov', 5, 'weights', 1);
    bnet.CPD{9} = gaussian_CPD(bnet, nodes_map('SPT.lambda'), 'mean', 0.1, 'cov', 0.01);
    bnet.CPD{10} = gaussian_CPD(bnet, nodes_map('SPT.k'), 'mean', 10, 'cov', 1);
    bnet.CPD{11} = gaussian_CPD(bnet, nodes_map('SPT.Npatch'), 'mean', 6, 'cov', 0.5);
    bnet.CPD{12} = gaussian_CPD(bnet, nodes_map('SPT.Nisg'), 'mean', 300, 'cov', 20);
    bnet.CPD{13} = gaussian_CPD(bnet, nodes_map('SPT.Rpbc'), 'mean', 4, 'cov', 0.1);
    bnet.CPD{14} = gaussian_CPD(bnet, nodes_map('SPT.G'),   'mean', Gm(1), 'cov', 2);
    bnet.CPD{15} = gaussian_CPD(bnet, nodes_map('SPT.Gref'),   'mean', Gm(1), 'cov', 2,'weights', 1);
    bnet.CPD{16} = gaussian_CPD(bnet, nodes_map('SPT.Gobs'), 'mean', Gm(1), 'cov', 2, 'weights', 1);
    bnet.CPD{17} = gaussian_CPD(bnet, nodes_map('SPT.I'),   'mean', Im(1), 'cov', 5);
    bnet.CPD{18} = gaussian_CPD(bnet, nodes_map('Reference.I'),   'mean', Im(1), 'cov', 5, 'weights', [0.5 0.5]);
    bnet.CPD{19} = gaussian_CPD(bnet, nodes_map('SPT.Iobs'), 'mean', Im(1), 'cov', 5, 'weights', 1);
  
    % eclass2
    bnet.CPD{20} = gaussian_CPD(bnet, nodes_map('ODE.G')+n, 'mean', Gm(time), 'cov', 0.1,  'weights', 0.5);
    
    % CPD for I(t+1), assume for now all parents are continuous
    parents_I1= parents(bnet.dag, nodes_map('ODE.I')+n); % parents of I(t+1)
    weights_I1_map_T0= containers.Map(); % parents in slice t
    weights_I1_map_T1= containers.Map(); % parents in slice t+1
    weights_I1_map_T0('ODE.G')= 0.05;
    weights_I1_map_T0('ODE.I')= 0.6;
    weights_I1_map_T1('ODE.alpha')= 0.6;
    weights_I1_map_T1('ODE.beta')= 0.6;
    weights_I1_map_T1('ODE.h')= 0.6;
    weights_I1= zeros(length(parents_I1),1);
    for i=1:length(parents_I1) 
        parent_index= parents_I1(i);
        if (parent_index <= n) % parent in slice t
            parent_name= reverse_nodes_map(parent_index);
            weights_I1(i)= weights_I1_map_T0(parent_name);
        else % parent in slice t+1
            parent_name= reverse_nodes_map(parent_index-n); % we only have a map for 1..9
            weights_I1(i)= weights_I1_map_T1(parent_name);
        end
    end
    bnet.CPD{21} = gaussian_CPD(bnet, nodes_map('ODE.I')+n,   'mean', Im(time), 'cov', 5, 'weights', weights_I1);
    %bnet.CPD{16} = gaussian_CPD(bnet, nodes_map('I')+n,   'mean', 72, 'cov', 5);
    weights_G= [0.5 0.5];
    bnet.CPD{22} = gaussian_CPD(bnet, nodes_map('SPT.G')+n, 'mean', Gm(time), 'cov', 2, 'weights', weights_G);
    
    % CPD for I(t+1), assume for now all parents are continuous
    parents_I1= parents(bnet.dag, nodes_map('SPT.I')+n); % parents of I(t+1)
    weights_I1_map_T0= containers.Map(); % parents in slice t
    weights_I1_map_T1= containers.Map(); % parents in slice t+1
    weights_I1_map_T0('SPT.G')= 0.05;
    weights_I1_map_T0('SPT.I')= 0.6;
    weights_I1_map_T1('SPT.k')= 0.6;
    weights_I1_map_T1('SPT.Npatch')= 0.6;
    weights_I1_map_T1('SPT.Nisg')= 0.6;
    weights_I1_map_T1('SPT.Rpbc')= 0.05;
    weights_I1= zeros(length(parents_I1),1);
    for i=1:length(parents_I1) 
        parent_index= parents_I1(i);
        if (parent_index <= n) % parent in slice t
            parent_name= reverse_nodes_map(parent_index);
            weights_I1(i)= weights_I1_map_T0(parent_name);
        else % parent in slice t+1
            parent_name= reverse_nodes_map(parent_index-n); % we only have a map for 1..9
            weights_I1(i)= weights_I1_map_T1(parent_name);
        end
    end
    
    bnet.CPD{23} = gaussian_CPD(bnet, nodes_map('SPT.I')+n, 'mean', Im(time), 'cov', 5, 'weights', weights_I1);
    %bnet.CPD{13} = gaussian_CPD(bnet, nodes_map('I')+n, 'mean', 70, 'cov', 5);
    
end
