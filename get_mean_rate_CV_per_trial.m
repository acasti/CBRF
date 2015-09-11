function [rate,CV,numspikes] = get_mean_rate_CV_per_trial(M,T)
%------------------------------------------------------------------------------------------------------
% Get mean firing rates over each data trial in a spike matrix (cell) M, and calculate the overall 
%  trial-averaged mean firing rate.
%
% USAGE:      [rate,CV,numspikes] = get_mean_rate_CV_per_trial(M,T);
% INPUT:      M                    * matrix of spike times in each trial (zero padded)
%             T                    * (double) time duration of each trial, spike times in [0,T]
% OUTPUT:     rate.trial           * (double vec) mean firing rate in each trial
%             rate.avg             * (double) mean firing rate averaged across all trials
%             CV.trial             * (double vec) coefficient of variation in each trial
%             CV.all               * (double) coefficient of variation for all data (merged ISIs)
%
% Dependencies:  get_numspikes_each_row.m, spiketime_mat2cell.m, get_isi_each_row.m
%
% Written by Alex Casti, Cooper Union Department of Mathematics, 17 Sept 2008
%------------------------------------------------------------------------------------------------------

if nargin < 2
  error('Not enough inputs!');
end

numtrials = size(M,1);                  % Number of data trials
numspikes = get_numspikes_each_row(M);  % Returns number of spikes in each trial
rate.trial = numspikes/T;               % Mean rate for each trial
rate.avg = mean(rate.trial);            % Mean rate averaged across trials
ISI_temp = get_isi_each_row(M);         % Matrix (with zeros) of ISIs in each trial
ISI = spiketime_mat2cell(ISI_temp);     % Cell array of ISIs for each trial

%% Get Coefficient of Variation (of ISIs) for each trial
CV.trial = zeros(1,numtrials);
for i = 1:numtrials
  CV.trial(i) = std(ISI{i})/mean(ISI{i});
end
ISI_vector = cell2mat(ISI);
CV.all = std(ISI_vector)/mean(ISI_vector);
