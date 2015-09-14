function [B, numspikes] = get_spike_identification_matrix(M)
%----------------------------------------------------------------------------------------------------
% Returns a matrix B that contains the trial number of a spike (column 1), which spike number it was in the trial
%   (column 2), the acutal spike time (column 3), and the interspike interval (ISI) associated with the spike
%   (column 4).
%
% USAGE:      [B, numspikes] = get_spike_identification_matrix(M);
% INPUT:      M                 * matrix of spike times (zero padded, original times)
% OUTPUT:     B                 * output spike data matrix (see below)
%             numspikes         * (vector) number of spikes in each trial
% Note: The output is a matrix B with the following columns:
%                               column 1 = trial #
%                               column 2 = spike # in trial
%                               column 3 = spike time (original time)
%                               column 4 = ISI associated with spike (1st spike in trial has ISI given by spike time)
%
% Dependencies:   get_isi_each_row.m
% Comments:
%  (1) Fixed problem that can occur if some trials have no spikes.  There was also a problem
%      if the maximum number of spikes in any trial was just 1.  (FIXED)
%
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 11 September 2015
%----------------------------------------------------------------------------------------------------

if ~isnumeric(M)
  M = spiketime_cell2mat(M);   % Input spikes were given as cell array; convert to zero padded matrix. 
end

inrows = size(M,1);
outrows = length(find(M(:)));   % Number of non-null spike times (M is zero padded)
outcols = 4;                          % 4 data values attached to each spike
B = zeros(outrows,outcols);
numspikes = zeros(1,inrows);
% Get interspike intervals (row by row)
ISI = get_isi_each_row(M);
% Cycle through the spike times and define B
i0 = 1;
for i = 1:inrows
  numspikes(i) = length(find(M(i,:)));
  if numspikes(i) > 0   % Don't add any information if a trial has no spikes!
    irange = i0:i0+numspikes(i)-1;
    B(irange,1) = i*ones(numspikes(i),1);      % column 1 has trial number
    B(irange,2) = (1:numspikes(i))';           % column 2 has spike number in trial
    B(irange,3) = M(i,1:numspikes(i))';        % column 3 has actual spike times t (original time)
    B(irange,4) = ISI(i,1:numspikes(i));       % column 4 has ISI values for each spike (original time)
    i0 = irange(end)+1;
  end
end

