function [D,Dmedian] = get_median_spikemetric_across_trials(M1,M2,q)
%--------------------------------------------------------------------------------------------------------------
% Calculate all spike train metric values for sets of 
%   spike responses.  In this version, metric distances are compared for
%   trials at the same level only, and thus returns a cell array
%   D{q(i)}(:) that has the same number of entries as there are trials.
%   This is to be used when comparing "distances" between spike trains arising 
%   from different stimulus parameters, or between model verus experiment, for
%   example.
% 
% Uses Jon Victor's spike time metric.
%
% USAGE:    [D,Dmedian] = get_median_spikemetric_across_trials(M1,M2,q);
% INPUT:    M1{1:N1,:}                       * matrix 1 of spike times (zero padded) 
%           M2{1:N2,:}                       * matrix 2 of spike times (zero padded) 
%           q                                * cost parameter (cost per unit time to translate a spike)
%                                              (q may be a vector of cost parameters)
% OUTPUT:   D{1:length(q)}(min(N1,N2))       * vector of raster line distances (cell array if length(q) > 1)
%           Dmedian(1:length(q))             * median of distance distribution at each q value
%
% Dependencies:  spkd_acmex.m, spiketime_mat2cell.m, get_numspikes_each_row.m
%
% Written by Alex Casti, MSSM, 21 Feb 2006
% Last updated 04 March 2008 (for time-rescaling project)
%--------------------------------------------------------------------------------------------------------------

if nargin < 3
  help get_median_spikemetric_across_trials.m;
  error('Not enough inputs!');
end

% Make sure that both the input matrices have at least 1 non-zero spike time.  If
%  one of the input matrices has no spike times (i.e. all zeros) then return
%  as the cost the number of spikes in the non-null matrix
% First count the number of spikes in each data set
numspikes1 = get_numspikes_each_row(M1);  % gets number of non-null spike times per trial
numspikes2 = get_numspikes_each_row(M2);
N1 = sum(numspikes1);
N2 = sum(numspikes2);
if (N1==0) && (N2==0)
  error('Input to get_median_spikemetric_across_trials.m had no spikes in either data set!')
end
if (N1==0) || (N2==0)
  if N1 == 0
    fprintf('Input M1 had no spikes! Using median spikes per trial in M2 for cost...\n');
    Dmedian = median(numspikes2);  % Median cost is median number of spikes (per trial) in M2
  end
  if N2 == 0
    fprintf('Input M2 had no spikes! Using median spikes per trial in M1 for cost...\n');
    Dmedian = median(numspikes1);  % Median cost is median number of spikes (per trial) in M1
  end
  D = Dmedian;
  return
end

% Convert the zero-padded matrices M1 and M2 into cell arrays
S1 = spiketime_mat2cell(M1);
S2 = spiketime_mat2cell(M2);

N1 = length(S1);
N2 = length(S2);
D = cell(1,length(q));
Dmedian = zeros(1,length(q));
Ncompare = min(N1,N2); 

for nq = 1:length(q)
  disp(sprintf('Calculating D for cost parameter %d of %d:  q = %g',nq,length(q),q(nq)));
  D{nq} = zeros(1,Ncompare);
  for i = 1:Ncompare
    D{nq}(i) = spkd_acmex(S1{i},S2{i},q(nq));
  end
  Dmedian(nq) = median(D{nq});
end

% Do not return a cell array if only 1 value of q is given
if length(q) == 1
  Dtemp = D{1}; 
  D = Dtemp;
end