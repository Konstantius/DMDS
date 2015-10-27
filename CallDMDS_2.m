function [Xq, Fval, Exitflag, Output, Population, Score] = CallDMDS_2(PopSize, ...
NofGenerations, StallLimit, NofEps, FName, SheetName, RangeName)
global N NR KN K MT Anew eps ObsTbl InitCorrMtr PR;
% Pearson correlation,  von Liebig approach
%  CallDMDS_2(PopSize, NofGenerations, StallLimit, NofEps, FName, SheetName, RangeName)
% --------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------
%          Variables K, N, Cmatz, Cmat should be declared as global!
%           global K N Cmatz Cmat; K = KValue; N = 5; SteadyStType = 2; Cmat = ones(N,N); Cmatz = ones(N,N)+2;
% PopSize           -- size of population
% NofGenerations -- number of generations
% NofEps  -- ratio of mean steady vector component taken as accuracy (only when SteadyStType=2)
% SteadyStType -- type of Steady state vector calculation
% SteadyStType =
%                      1: Power method (iteration, use of eps) 
%                      2: LU Decomposition if Q^T (Q -- infinitismal generator)

KN = K^N;
NR = N^2;

eps = 1/KN/NofEps;

% Observation table
ObsTbl = zeros(KN, N);
for i = 1:KN 
 ObsTbl(i, :) = BasePres(N, K, i-1)+1;
end;

PR = zeros(3^N, N);
for i =1:3^N
  PR(i, 1:N) = BasePres(N, 3, i-1);
end;

MT = zeros(KN, KN);
Anew = zeros(1, N);

lb=zeros(1, NR);
lb(:) = -1;
ub=ones(1, NR);

InitCorrMtr = xlsread(FName, SheetName, RangeName);

[Xq, Fval, Exitflag, Output, Population, Score]=StartGA2(NR, lb, ub, PopSize, NofGenerations, StallLimit);
