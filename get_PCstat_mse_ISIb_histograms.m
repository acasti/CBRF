function errdata = get_PCstat_mse_ISIb_histograms(ISIdata,refindices)
%-----------------------------------------------------------------------------------------------
% For any individual data set, obtain a measure of "perfect copy" by pairwise comparison 
%  of time B ISI histograms.  Here we calculate a mean square error between a reference
%  histogram (typically for a moderate rate and small rate derivative) and all the others.
%
% USAGE:    errdata = get_PCstat_mse_ISIb_histograms(ISIdata,refindices)
% INPUT:    ISIdata             * (struct) output of isi_partition.m
%           refindices          * (int vec) length 2 integer vector indicating the reference
% OUTPUT:   errdata             * (struct) output error data (pairwise mse and mean)
% 
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 13 September 2015
%-----------------------------------------------------------------------------------------------

%% Argument check
if nargin < 2
  error('Insufficient number input arguments!');
end

mse = zeros(size(ISIdata.Nhist_normalized));
[numrates,numdrates] = size(ISIdata.Nhist_normalized);

%% Calculate the mean square errors across histograms
refhist = ISIdata.Nhist_normalized{refindices(1),refindices(2)};
for i = 1:numrates
  for j = 1:numdrates
    mse(i,j) = sum( (refhist-ISIdata.Nhist_normalized{i,j}).^2 );
  end
end

%% Define the mean of all mean square departures, excluding the zero value from the mean square departure of the 
%   reference histogram with itself.
mse_mean = sum(mse(:))/(numel(ISIdata.Nhist_normalized)-1);

%% Define another mean mse measure by excluding extremely low and high rates.  Call it mse_mean_exclude_lowhigh.
M = mse(2:end,:);
mse_mean_exclude_lowhigh = mean(M(:))/(numel(M)-1);

%% Define output data structure
errdata.mse = mse;
errdata.mse_mean = mse_mean;
errdata.mse_mean_exclude_lowhigh = mse_mean_exclude_lowhigh;
