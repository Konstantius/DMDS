function [outs] = CORPEAR(x)
% x -- vector of relationships'  matric entries
% x = 
% m_11, m_12, m_1N
% m_21, m_22, m_2N
% ................
% m_N1, m_N2, m_NN
global N NR PearCorrMat KN K BA delta eps ObsTbl InitCorrMtr;

P=zeros(N, 3);
Ps=zeros(1, 3^N);
zls=0;
M =zeros(N, N);
MT=zeros(KN, KN);
CorrMtr = ones(N, N);
for i1 = 1:N
  M(i1, :) = x(1+(i1-1)*N:i1*N);
end;

for i1 = 1:N
%   m1 =min(x(1+(i1-1)*N:i1*N));
%   m2 = max(x(1+(i1-1)*N:i1*N));
if min(x(1+(i1-1)*N:i1*N))>0 || max(x(1+(i1-1)*N:i1*N)) < 0
  outs = N*(N-1);
  return;
end
end;


% Ps := Vector[row](3^N):                                       # # vector of probabilities
for i0 =1:KN                                                        %%      State's number
Astate=  BasePres(N, K, i0-1)+1;

  
 for j =1:N                                       %%                 # #       cycle by components
 PL =0; PC =0; PU =0;

for k =1:2^N                                    %   cycle by combinations k1, k2, kN
  p = 1;
 for j1 =1:N
   p=p*((2*abs(M(j, j1))-1)*BA(k, j1)+1-abs(M(j, j1)));   
 end;
 d =0;
 for j1 = 1:N
    d = d+sign(M(j, j1))*BA(k, j1)*(Astate(j1)-zls); 
 end;
 
if d < -delta
    PL =PL +p;
elseif d > delta 
    PU =PU+p;
elseif (d <= delta) && (d >= -delta) 
    PC =PC+p;
end;                     
end;                                                      %   /cycle by combinations k1, k2, kN
P(j,1)=PL; P(j,2)=PC; P(j,3) = PU;

end;                                 %%                                     /cycle by components

for j = 1:3^N                                   %                      cycle for calc of probabilities d1, d2, ..., dN
    S = BasePres(N, 3, j-1)+1;
    p=1;
    for j1 = 1:N 
      p =p*P(j1, S(j1));
    end;    
Ps(j) = p;
    
end;                                                   %%             /cycle for calc of probabilities d1, d2, ..., dN
% 
for i = 1:3^N                                      %   <1>
if Ps(i) > 0                                          % i -- direction of change (3^N in all!)
Si0 = BasePres(N, 3, i-1)-1;                 % Si0 -- vector of change
for j = 1:N 
  if Si0(j) ==-1
     Anew(j) = max(1, Astate(j)-1);           % Dec
  elseif Si0(j) == 0 
     Anew(j) = Astate(j);
  elseif Si0(j) == 1 
     Anew(j) = min(Astate(j)+1, K);                  % Inc 
end;
 end;
 Anew= Anew-1;
i2 = R2Bas(Anew, K)+1;
MT(i0, i2) =MT(i0, i2) +Ps(i);
end;                                                                %    </1>

 end;
end;
% xlswrite('d:\out.xlsx', MT, 'new-sheet');

% Steady-state vector
PT= MT';
x =ones(KN, 1)/KN;

SSt =PT * x;
while norm(x-SSt, 2)> eps
    x = SSt;
    SSt =PT * x; 
end;


% Correlation matrix
for i1 = 1:N-1
    mxw = ObsTbl(:, i1)' * SSt;
    sigx2 = ((ObsTbl(:, i1)-mxw).^ 2)' * SSt;
        for i2 = i1+1:N
            myw = ObsTbl(:, i2)' * SSt;
            sigy2 = ((ObsTbl(:, i2)-myw).^ 2)' * SSt;
            covxy = (((ObsTbl(:, i1)-mxw) .* (ObsTbl(:, i2)-myw)))' * SSt;
            CorrMtr(i1, i2) = covxy/sqrt(sigx2*sigy2);  
        end;
end;
p=0;
for i = 1:N-1
    for i1 =i+1:N
        p=p+(InitCorrMtr(i, i1)-CorrMtr(i, i1))^2;
    end;
end;

outs = p;


end