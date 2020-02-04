clc
clear all
 
m= load('FINAL_TABLE_NSGAIII.txt');


%Best Metric values Borg
 
%Best GD
%BestGD=1-min(m(:,2));
 
%Best Epsilon normalized
%Bestepsilon=1-min(m(:,5));
 
% Best hypervolume
Besthypervol=max(m(:,3));

RFHV=0.5738930505; %Reference set hypervolume
NormalizedHV= Besthypervol/RFHV; % Besthypervolume/reference_set


