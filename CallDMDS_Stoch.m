global N NR KN K BA delta MT Anew eps ObsTbl InitCorrMtr;
N = 4;
K = 3;
KN = K^N;
NR = N^2;
% Observation table
ObsTbl = zeros(KN, N);
for i = 1:KN 
 ObsTbl(i, :) = BasePres(N, K, i-1)+1;
end;

BA = zeros(2^N, N);
for i =1:2^N 
  BA(i,:) = BasePres(N, 2, i-1);
end;
eps = 1.0/KN/1000;      % Accurancy of steady-state

Anew = zeros(1, N);
delta = 0.5;
P = zeros(N, 3);
InitCorrMtr=[1       0.2    0.3    -0.5; ...
                      0.2     1    -0.48   0.7; ...
                      0.3 -0.48     1     -0.5; ...
                      -0.5  0.7     -0.5    1];   
    

M = [-0.2 -0.8 0.45 0.87 ...
    0.4 0.6 -0.7 -0.6 ...
    0.3 -0.4 -0.7 0.7 ...
    -0.8 0.3 0.3 0.7];
outp =CORPEAR(M);
disp(outp);