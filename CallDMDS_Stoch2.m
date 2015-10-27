global Cmatz Cmat N NR KN K   Anew eps ObsTbl InitCorrMtr PR;
N = 5;
K = 3;
KN = K^N;
NR = N^2;

% Observation table
ObsTbl = zeros(KN, N);

Cmatz = ones(N, N)+2;
Cmat =ones(N, N);

for i = 1:KN 
 ObsTbl(i, :) = BasePres(N, K, i-1)+1;
end;

eps = 1.0/KN/1000;      % Accurancy of steady-state

Anew = zeros(1, N);
% delta = 0.5;
% P = zeros(N, 3);

InitCorrMtr= xlsread('c:\MDS\MATLAB\DMDS\TestTable.xlsx', 'TestSheet', 'B12:F16');
% 
% [1       0.2    0.3    -0.5; ...
%                       0.2     1    -0.48   0.7; ...
%                       0.3 -0.48     1     -0.5; ...
%                       -0.5  0.7     -0.5    1];   
 PR = zeros(3^N, N);
 for i =1:3^N
   PR(i, 1:N) = BasePres(N, 3, i-1);   
 end;
 
% M = [-0.2 -0.8 0.45 0.87 ...
%     0.4 0.6 -0.7 -0.6 ...
%     0.3 -0.4 -0.7 0.7 ...
%     -0.8 0.3 0.3 0.7];


M = xlsread('c:\MDS\MATLAB\DMDS\TestTable.xlsx', 'TestSheet', 'B2:F6');
M1 = [];
for i =1:N
M1 = [M1 M(i, 1:N)];
end;


outp =CORPEAR2(M1);
disp(outp);