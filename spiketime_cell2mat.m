function M = spiketime_cell2mat(S)
%----------------------------------------------------------------------------------------------------------------
% Convert a cell array of spike times (i.e. each array entry being a different sweep of data) and convert the cell
%   to a matrix of spike times, with zeros added to the end of rows if necessary.
%
% USAGE:           M = spiketime_cell2mat(S)
% INPUT:           S                    * length N cell array of spike times
% OUTPUT:          M                    * (N x maxrow) matrix of spike times (zero padded), where 'maxrow' is the maximum
%                                          number of spikes in any sweep of the input cell array S.
%
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 10 September 2015
%----------------------------------------------------------------------------------------------------------------

if ~iscell(S)
  error('Input must be a cell array...');
end

N = length(S);   % Number of sweeps of spike time data

% Determine the maximum number of spikes in each sweep
maxspikes = 0;
for i = 1:N
  maxspikes = max(maxspikes,length(S{i}));
end

M = zeros(N,maxspikes);
for i = 1:N
  M(i,:) = [S{i}(:)' , zeros(1,maxspikes-length(S{i}))];
end
