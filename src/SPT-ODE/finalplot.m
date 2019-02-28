% plot

% Read in the experimental measurements
fid = fopen('k.txt','rt');
indata = textscan(fid, '%f', 'HeaderLines',2);
fclose(fid);
disp(indata)

lines = regexp(k, '\r\n|\r|\n|\r|\n|\r|\n', 'split');
disp(lines{1}(1));
sptk = k(1); 
disp(sptk);

metak = k(:,3); 

h = importdata('h.txt');
sptk = h(:,1); 
metak = h(:,3); 

disp(sptk)
% plot
figure()
%yyaxis left;
%plot(1:T, sample_seq(nodes_map('SPT.G'),:));
plot(I1,spt_marg.mu);
plot(I1,marg.mu);
%hold off;