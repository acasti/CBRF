function data = get_bigspikematrix(M,Mb,Mlam,Mdlam)
%----------------------------------------------------------------------------------------------------
% Returns a matrix B that keeps track of spike times, spike ordering, which trial that the spike falls,
%  and rate parameter information.
%
% USAGE:       data = get_bigspikematrix(M,Mb,lam,dlam);
% INPUT:       M                 * matrix of spike times (zero padded, original times)
%              Mb                * matrix of spike times (zero padded, new time)
%              Mlam              * matrix of rate values (zero padded, original time)
%              Mdlam             * matrix of rate derivative values (zero padded, original time)
% OUTPUT:      data              * data struct (with B and comments)
% 
% Note: The output data struct includes a matrix B with the following columns:
%                                    B columns
%                               column 1 = trial #
%                               column 2 = spike #
%                               column 3 = spike time (original time)
%                               column 4 = spike time (new time)
%                               column 5 = rate value (original time)
%                               column 6 = rate derivative value (original time)
%                               column 7 = ISI (original time)
%                               column 8 = ISI (new time)
%
% Dependencies:   get_isi_each_row.m
%
% Written by Alex Casti, MSSM, Nov 2007
% Last modified 26 Feb 2008
% Comments:
%  (1) Fixed problem that can occur if some trials have no spikes.  There was also a problem
%      if the maximum number of spikes in any trial was just 1.  (FIXED)
%----------------------------------------------------------------------------------------------------

if nargin < 4
  error('Not enough inputs!');
end

inrows = size(M,1);
outrows = length(find(M(:)));   % Number of non-null spike times (M is zero padded)
outcols = 8;                    % Eight data values attached to each spike
B = zeros(outrows,outcols);
numspikes = zeros(1,inrows);
% Get interspike intervals (row by row)
ISI = get_isi_each_row(M);
ISIb = get_isi_each_row(Mb);
% Cycle through the spike times and define B
i0 = 1;
for i = 1:inrows
  numspikes(i) = length(find(M(i,:)));
  if numspikes(i) > 0  % Don't add any information if a trial has no spikes!
    irange = i0:i0+numspikes(i)-1;
    B(irange,1) = i*ones(numspikes(i),1);    % column 1 has trial number
    B(irange,2) = (1:numspikes(i))';         % column 2 has spike number in trial
    B(irange,3) = M(i,1:numspikes(i))';      % column 3 has actual spike times t (original time)
    B(irange,4) = Mb(i,1:numspikes(i))';     % column 4 has spike times (new time)
    B(irange,5) = Mlam(i,1:numspikes(i))';   % column 5 has lam(t) values for each spike (rate)
    B(irange,6) = Mdlam(i,1:numspikes(i))';  % column 6 has dlam/dt(t) values for each spike (rate derivative)
    B(irange,7) = ISI(i,1:numspikes(i));     % column 7 has ISI values for each spike (original time)
    B(irange,8) = ISIb(i,1:numspikes(i));    % column 8 has ISI values for each spike (new time)
    i0 = irange(end)+1;
  end
end
comments = {'trial','spike#','spike time (orig)','spike time (new)','lam','dlam','ISI (orig)','ISI (new)'};
data.B = B;
data.comments = comments;
data.numspikes = numspikes;
