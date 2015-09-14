function S = spiketime_mat2cell(M)
%----------------------------------------------------------------------------------------------------------------
% Convert a matrix of spike times (each row being a different sweep of data, possibly with trailing zeros) 
%   into a cell array of spike times (each cell entry is one trial of spike time data, with no trailing zeros).
%
% USAGE:         S = spiketime_mat2cell(M);
% INPUT:         M            * (N x maxrow) matrix of spike times (zero padded), where 'maxrow' is the maximum
%                                 number of spikes in any sweep.
% OUTPUT:        S            * (1 x N) cell array of spike times with trailing zeros removed
%
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 10 September 2015
%----------------------------------------------------------------------------------------------------------------

if ~isnumeric(M)
  error('Input must be a matrix of spike times (not a cell)...');
end

N = size(M,1);   % Number of sweeps of spike time data
S = cell(1,N);   % Output cell array

% Determine the maximum number of spikes in each sweep
for i = 1:N
  row = M(i,:);
  S{i} = row(row>0);  % Eliminate fake times from zero padding
end

