function [outs] = CORPEAR2(x)
% x -- vector of relationships'  matric entries
% x = 
% m_11, m_12, m_1N
% m_21, m_22, m_2N
% ................
% m_N1, m_N2, m_NN
% For von Liebig low
global N KN K eps ObsTbl InitCorrMtr SteadyStType Cmatz Cmat PR Anew;
% NR PearCorrMat 

% Ps=zeros(1, 3^N);

ML =zeros(N, N);            %         Relationship matrix
MT=zeros(KN, KN);      
CorrMtr = ones(N, N);
for i1 = 1:N                   %  Fills relatioship matrix
  ML(i1, :) = x(1+(i1-1)*N:i1*N);
end;

P=zeros(N, 3);   
FP = zeros(N, 3);
% Ps := Vector[row](3^N):                                       # # vector of probabilities
for i0 = 1:KN                                                        %%      State's number
Astate =  BasePres(N, K, i0-1)+1;                           %%      Current state

% 
% P =  zeros(N, 3);                                              % Probs for components
 for i = 1:N                                                      %         cycle by components

for j = 1:N                                                        %                     <j>  
  if ML(i, j) < 0    
      if Astate(j) >= Cmatz(i, j) 
        FP(j, 1) = abs(ML(i, j)); FP(j, 2) = 0;   FP(j, 3) = 1-abs(ML(i, j));
      end;
      if (Astate(j) > Cmat(i, j)) && (Astate(j) < Cmatz(i, j)) 
       FP(j, 1) = 0.5; FP(j, 2) = 0; FP(j, 3) = 0.5;    
      end;
      if Astate(j) <= Cmat(i, j) 
        FP(j, 1) = 1-abs(ML(i, j)); FP(j, 2) = 0; FP(j, 3) = abs(ML(i, j));
      end;  
  end;
    if ML(i, j) > 0 
    if Astate(j) >= Cmatz(i, j) 
      FP(j, 1) = 1-ML(i, j);    FP(j, 2) = 0;    FP(j, 3) = ML(i, j);
    end;
    if (Astate(j) > Cmat(i, j)) && (Astate(j) < Cmatz(i, j))
       FP(j, 1) = 0.5; FP(j, 2) = 0; FP(j, 3) = 0.5;    
    end;
    if Astate(j) <= Cmat(i, j) 
      FP(j, 1) = ML(i, j);    FP(j, 2) = 0;  FP(j, 3) = 1-ML(i, j);
    end;
    end;
 end;                                                             %        </j>

PL =0; PC =0; PU =0;
  for j = 1:3^N                                             %             <1>
   Ps = 1;
    for k=1:N
      Ps =Ps * FP(k, PR(j, k)+1);
    end;
    if  Ps > 0 
      if min(PR(j, 1:N)-1) == -1 
        PL = PL+Ps;
      elseif min(PR(j, 1:N)-1) == 0
        PC = PC+Ps;
      elseif min(PR(j, 1:N)-1) == 1 
        PU = PU+Ps;
      end;  
    end;  
  end;                                                            % #             </1>
P(i,1) =PL; P(i,2) =PC; P(i,3) = PU;                                                                 
 end;                                                  %                 /cycle by components
for i = 1:3^N                                               %             <1>
  Ps = 1;
  for k = 1:N
    Ps =Ps* P(k, PR(i, k)+1);
  end;
  
  if  Ps > 0                                                      % # Build Anew 
    TVec = PR(i, 1:N)-1;
      for j = 1:N                                                %         <2>
        if TVec(j) == -1 
          Anew(j) = max(1, Astate(j)-1);
        elseif   TVec(j) == 0 
          Anew(j) = Astate(j);
        elseif   TVec(j) == 1 
          Anew(j) = min(K, Astate(j)+1);
        end;
      end;                                                         %             </2>
      if isequal(Anew, Astate) == 1
        MT(i0, i0) = MT(i0, i0)  + Ps;
      else 
        i2 = R2Bas(Anew-1, K)+1;
        MT(i0, i2) = MT(i0, i2)  + Ps; 
      end;
  end;
end;                                                                %             </1>

end;            %  Cycle by all states of components

csvwrite('d:\matl.dat', MT);

% xlswrite('d:\out.xlsx', MT, 'new-sheet');
% Steady-state vector calculation
% ===================================
if SteadyStType == 1 
% Power method !!
 PT= MT';
 x =ones(KN, 1)/KN;
 SSt =PT * x;
 while norm(x-SSt, 2)> eps
     x = SSt;
     SSt =PT * x; 
     
 end;

 elseif SteadyStType == 2
% Direct method !!
Q = MT;
 RM= sum(MT, 2);
 for i =1:KN
   Q(i,i) = MT(i,i)-RM(i);
 end;
QT = Q';
[L,UP] = lu(QT);
  SSt = zeros(KN, 1);
  SSt(KN) = 1; 
 for i = KN-1:-1:1 
   d1 = 0;
   for j = i+1:KN
     d1 = d1 + UP(i, j)*SSt(j);     
   end;
     SSt(i) = -d1/UP(i, i);
 end;
SSt = SSt/sum(SSt);
end;
% ===================================
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