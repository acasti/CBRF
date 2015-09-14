function M = spikevec2matrix(v,N)
%------------------------------------------------------------------------------------------
% Convert a vector of spike times (i.e. with "trials" concatenated) into a matrix 
%  of spike times with trailing zeros (each row of the matrix is a trial of 
%  spike times with zero padding).
%
% USAGE:     M = spikevec2matrix(v,N)
% INPUT:     v                   * vector of spike times
%            N                   * number of non-null spike times per row 
% OUTPUT:    M                   * spike times (with trailing zeros) in matrix form, with
%
% Comment:   The size of the output matrix M is  length(N) x max(N), since 
%            length(N) indicates the number of trials and max(N) the maximum
%            number of spikes per trial (must be known for the zero padding).
%
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 10 September 2015
%------------------------------------------------------------------------------------------

if (nargin < 2) || any([isempty(v),isempty(N)])
  help spikevec2matrix.m;
  error('Two input arguments required!');
end

% Length of input vector v must equal the number of nonzero spikes indicated by N
if sum(N) ~= length(v)
  error('length(v) must equal the total number of spikes in N vector');
end

% Make v a row vector if it isn't
v = v(:)';

rows = length(N);
cols = max(N);
M = zeros(rows,cols); 
i1 = 1;
for i = 1:rows
  M(i,1:N(i)) = v(i1:i1+N(i)-1);
  i1 = i1 + N(i);
end