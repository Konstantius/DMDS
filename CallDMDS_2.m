global N NR KN K BA delta MT Anew eps ObsTbl InitCorrMtr;
N = 4;
K = 4;
KN = K^N;
NR = N^2;

BA = zeros(2^N, N);
for i =1:2^N 
  BA(i,:) = BasePres(N, 2, i-1);
end;

% Observation table
ObsTbl = zeros(KN, N);
for i = 1:KN 
 ObsTbl(i, :) = BasePres(N, K, i-1)+1;
end;


MT = zeros(KN, KN);
Anew = zeros(1, N);
delta = 0.5;
P = zeros(N, 3);
lb=zeros(1, N^2);
lb(:) = -1;
ub=zeros(1, N^2);
ub(:) = 1;
PopSize = 50;
NofGenerations = 60;
StallLimit = 50;
InitCorrMtr = xlsread('c:\MDS\MATLAB\DMDS\TestTable.xlsx', '@@1', 'B5:E8');

[Xq, Fval, Exitflag, Output, Population, Score]=StartGA(N^2, lb, ub, PopSize, NofGenerations, StallLimit);
% 
% 
% 
% M = [-0.2 -0.8 0.45 0.87 ...
%     0.4 0.6 -0.7 0.6 ...
%     0.3 -0.4 -0.7 0.7 ...
%     -0.8 0.3 0.3 0.7];
% CORPEAR(M);