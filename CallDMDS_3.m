function [Xq, Fval, Exitflag, Output, Population, Score] = CallDMDS_3(PopSize, ...
NofGenerations, StallLimit, NofEps, FName, SheetName, RangeName)
global N NR KN K BA delta MT Anew eps ObsTbl InitCorrMtr SteadyStType;
%  CallDMDS_3(PopSize, NofGenerations, StallLimit, NofEps, FName, SheetName, RangeName)
%          Variables K, N should be declared as global!
%           global K N; K = KValue; N = 5; SteadyStType = 2;
% PopSize           -- size of population
% NofGenerations -- number of generations
% NofEps  -- ratio of mean steady vector component taken as accuracy
% SteadyStType -- type of Steady state vector calculation
% SteadyStType =
%                      1: Power method (iteration, use of eps) 
%                      2: LU Decomposition if Q^T (Q -- infinitismal generator)

KN = K^N;
NR = N^2;
BA = zeros(2^N, N);
for i =1:2^N 
  BA(i,:) = BasePres(N, 2, i-1);
end;

eps = 1/KN/NofEps;

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



InitCorrMtr = xlsread(FName, SheetName, RangeName);

[Xq, Fval, Exitflag, Output, Population, Score]=StartGA(N^2, lb, ub, PopSize, NofGenerations, StallLimit);
