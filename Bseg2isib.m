function [ISIb, ISIb_all, rate_mean, drate_mean] = Bseg2isib(ISIdata)
%--------------------------------------------------------------------------------------------
% Just a utility routine to convert the Bseg cell array to a simpler cell array containing
%  only the time B ISI data for the 21 histograms.  The mean rate and mean rate derivative
%  associated with each subset of spikes is also returned.
%
% USAGE:     [ISIb, ISIb_all, rate_mean, drate_mean] = Bseg2isib(ISIdata);
% INPUT:     ISIdata         * (struct) output data structure from the isi_partition.m  routine
% OUTPUT:    ISIb            * (cell) same dimensions as ISIdata.Bseg, but each cell entry
%                                     will have only the time B ISI data.
%            ISIb_all        * (vector) vector of all time B ISIs (unsegregated).
%            rate_mean       * (matrix) matrix of mean rate values among all spikes in 
%                                       each corresponding entry of the ISIb cell array.
%            drate_mean      * (matrix) matrix of mean rate derivative values among all spikes in
%                                       each corresponding entry of the ISIb cell array.
% Comments:  Note the following information contained in each Bseg cell array entry (each entry
%            is a matrix with information about each spike subset segmented by rate and drate)
%      (1) The data field "data.Bseg" contains the array of interval data, for which each cell entry
%          is an 8-column data matrix with the following columns:
%          col 1 = trial #
%          col 2 = spike # in trial
%          col 3 = spike time t (original time)
%          col 4 = spike time B (time B spike time)
%          col 5 = lambda value for spike (rate function, original time)
%          col 6 = d/dt(lambda) for the spike 
%          col 7 = ISI associated with spike (original time)
%          col 8 = ISI associated with spike (time B)
%
% Written by Alex Casti, Cooper Union Math Dept, 20 Jan 2009
% Last updated 28 Jan 2009
%--------------------------------------------------------------------------------------------

%% Argument check
if nargin < 1
  error('Did not supply input to the function!');
end
if ~isstruct(ISIdata)
  error('You did not enter as input the proper ISIdata structure...');
end

%% Obtain the output data of the function
Bseg = ISIdata.Bseg;
ISIb_all = sort(ISIdata.ISIb);   % return a sorted vector of ISIs

[rows,cols] = size(Bseg);
ISIb = cell(rows,cols);
rate_mean = zeros(rows,cols);
drate_mean = zeros(rows,cols);
for i = 1:rows
  for j = 1:cols
    ISIb{i,j} = Bseg{i,j}(:,8);
    rate_mean(i,j) = mean(Bseg{i,j}(:,5));
    drate_mean(i,j) = mean(Bseg{i,j}(:,6));
  end
end