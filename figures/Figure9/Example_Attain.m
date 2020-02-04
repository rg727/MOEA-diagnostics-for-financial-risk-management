clc
clear all
 
m= load('Borg_Concatenate_Average_Metrics.txt');
 
%Best Metric values Borg
 
%Best GD
BestGD=1-min(m(:,2));
 
%Best Epsilon normalized
Bestepsilon=1-min(m(:,5));
 
% Best hypervolume
Besthypervol=max(m(:,1));
 
%Attainment
%Normalized metrics
% Hypervolume
RFHV=0.61; %Reference set hypervolume
NormalizedHV= Besthypervol/RFHV(1); % Besthypervolume/reference_set
 
k= 1:-0.01:0.0; % k is the threshold
n=1;
 
%Hypervolume
count=0;
count1=0;
count2=0;
algorithms=1;
percent_attainment=zeros(length(algorithms),length(k));
percent_attainmentgd=zeros(length(algorithms),length(k));
percent_attainmentei=zeros(length(algorithms),length(k));
for l=1:length(algorithms)
% load file
for n=1:length(k)
 
for i= 1:length(m(:,2))
 
if (m(i,1)/RFHV) >= k(n)
count=count+1;
percent_attainment(l,n)=(count/length(m(:,1)));
end
if (1-m(i,2)) >= k(n)
count1=count1+1;
percent_attainmentgd(l,n)=(count1/length(m(:,2)));
end
if (1-m(i,5)) >= k(n)
count2=count2+1;
percent_attainmentei(l,n)=(count2/length(m(:,2)));
end
 
end
count=0;
count1=0;
count2=0;
 
 
end
end
 
percent_attainment_c= {percent_attainment', percent_attainmentgd', percent_attainmentei'};
Best_metrics=[NormalizedHV, BestGD, Bestepsilon];
Names={'Hypervolume' 'Generational Distance' 'Epsilon Indicator'};

fid=fopen('Borg_Attainment.txt','w');
%fprintf(fid, [ header1 ' ' header2 '\n']);
fprintf(fid, '%f \n', percent_attainment);
fclose(fid);

fid=fopen('Borg_Attainment_epsilon.txt','w');
%fprintf(fid, [ header1 ' ' header2 '\n']);
fprintf(fid, '%f \n', percent_attainmentei);
fclose(fid);


fid=fopen('Borg_Attainment_gd.txt','w');
%fprintf(fid, [ header1 ' ' header2 '\n']);
fprintf(fid, '%f \n', percent_attainmentgd);
fclose(fid);

% l= [2,4,6];
% for i=1:length(percent_attainment)
%  
% subplot(1,3,i)
% imagesc(1, k, percent_attainment{i})
% set(gca,'ydir','normal')
% colormap(flipud(bone))
% hold on
% scatter(1, Best_metrics(i), 80, 'k', 'filled')
% set(gca,'XTick')
% xlabel(Names(i))
% ylabel('Probability of Best Metric Value');
% if i==2
% title('Attainment probabilities for Borg, DTLZ2_3')
% legend('Best metric value across all runs')
% end
% end