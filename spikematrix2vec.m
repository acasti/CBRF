function [v,N] = spikematrix2vec(M)
%------------------------------------------------------------------------------------------
% Convert a matrix of spike times (with trailing zeros) into a vector of spike times with the 
%  trailing zeros removed.  First entries of output vector v consists of non-null M(1,:), the
%  next entries contain the non-null second row M(2,:), etc.
%
% USAGE:     [v,N] = spikematrix2vec(M);
% INPUT:     M              * matrix of spike times (with trailing zeros)
% OUTPUT:    v              * vector of concatenated spike times
%            N              * int vector with number of non-zero events in each row of M
%
% Comment:   In order to reconstruct the original matrix M you must save the vector N
%            that contains the number of non-null entries per row.  Use the associated
%            Matlab script "spikevec2matrix.m"
%
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 10 September 2015
%------------------------------------------------------------------------------------------

% Pre-allocate memory for output v
[rows,cols] = size(M);
v = Inf*ones(1,rows*cols);
N = zeros(1,rows);

i1 = 1;
for i = 1:rows
   N(i) = length(find(M(i,:)));
	 v(i1:i1+N(i)-1) = M(i,1:N(i));
	 i1 = i1 + N(i);
end

% Remove the Infs
v = v( (v<Inf) );