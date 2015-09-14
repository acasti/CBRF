function ind = locate(xx,x)
%-----------------------------------------------------------------------------
% Find the indices in an ordered array xx(1:n) within which x falls.
% Uses a bisection algorithm (see Numerical Recipes, p.111).
%
% USAGE:   ind = locate(xx,x)
%
% INPUT:   xx = ordered array of length n
%           x = real number for which we determine 
%               the pair (xx(j),xx(j+1)) that x falls between.
% OUTPUT:   j = ind = index such that  xx(j) <= x < xx(j+1)
% 
% Note:     ind = 0 is returned if x is out of range and too small
%           ind = n is returned if x is out of range and too large
%
% Written by Alex Casti, FDU Department of Mathematics
% Last modified 13 September 2015
%-----------------------------------------------------------------------------

n = length(xx);
jl = 0;                % Initialize upper and lower limits
ju = n+1;

% Check to see if input value x is out of range
if x < xx(1)
  ind = 0;   % input x is less than every element of array xx
  return
end
if x > xx(end)
  ind = n;   % input x is greater than every element of array xx
  return
end

while ju - jl > 1
 jm = floor((ju+jl)/2);        % Compute a midpoint
 if (xx(n) > xx(1) && x > xx(jm)) || (~(xx(n) > xx(1)) && ~(x > xx(jm)))
  jl = jm;
 else
  ju = jm;
 end
end

ind = jl;
if xx(jl+1) == x
  ind = jl + 1;
end