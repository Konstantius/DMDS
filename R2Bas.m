function [ O1 ] = R2Bas(x, Ba)
global N;
%   Converts   Ba-base representation of a number in array x to decimal number
s=0;
for i =1:N
  s=s+x(N-i+1)*Ba^(i-1);
end
O1=s;

end

