function numspikes = get_numspikes_each_row(M)
%---------------------------------------------------------------------------------------------
% Get number of spikes per trial for a matrix of spike times (each row is a trial)
% 
% USAGE:           numspikes = get_numspikes_each_row(M);
% INPUT:           M              % m x n matrix of spike times (zero padded)
% OUTPUT:          numspikes      % length m vector of number of spike times per row (trial)
%
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 13 September 2015
%---------------------------------------------------------------------------------------------

if nargin < 1 
  error('Not enough inputs!');
end
if isempty(M)
  numspikes = [];   % return empty vector for empty input
  return
end

rows = size(M,1);           % number of trials
numspikes = zeros(1,rows);  % number of spikes per trial
for i = 1:rows
  numspikes(i) = length(find(M(i,:)));
end
