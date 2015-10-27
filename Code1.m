
KN = 60;
MT = xlsread('c:\MDS\MAPLE\TMP.xlsx', 'S1', 'B69:BI128');

     %       Transfinitisimal generator
t1 = now;
 Q = MT;
 RM= sum(MT, 2);
 for i =1:KN
   Q(i,i) = MT(i,i)-RM(i);
 end;
 
QT = Q';
[L,UP] = lu(QT);
W = zeros(KN, 1);
W(KN) = 1; 
 for i = KN-1:-1:1 
   d1 = 0;
   for j = i+1:KN
     d1 = d1 + UP(i, j)*W(j);     
   end;
     W(i) = -d1/UP(i, i);
 end;
W = W/sum(W);
t2 = now;
disp((t2-t1)*24*3600);
disp(W);
% 
% 
% 
eps = 1/KN/200;
t1 = now;
w = ones(KN, 1)/KN;
  MTT = MT';
  w1 = MTT * w;
 while  norm(w-w1, 2)> eps
% di LinearAlgebra[Norm](w-w1, infinity);
 w = w1;
 w1 = MTT * w;
 end;
t2 = now;
disp((t2-t1)*24*3600);
% t2-t1;
disp(w1);
