function [OAr] = BasePres(N0, bas, nb)
%   N0 -- width of output, bas - base, nb - number to be convert (from 0).
Ar = zeros(1,N0);
na=nb;
k=0;
while na >0
  Ar(N0-k)= mod(na, bas);
  na = floor(na/bas);
  k=k+1;
end;
OAr = Ar;

end

