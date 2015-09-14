function plot_lambar_dlambar(M,lambar,dlambar)
%---------------------------------------------------------------------------------------------------
% Plot lambar (mean lambda between intervals) vs dlambar (mean time derivative of lambar)
%
% USAGE:     plot_lambar_dlambar(M,lambar,dlambar)
% INPUT:     M               * matrix of spike times (so we know number of events; note zero padding)
%            lambar          * matrix of mean rate parameters for each spike
%            dlambar         * matrix of time derivative of mean rates for each spike
% 
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 13 September 2015
%---------------------------------------------------------------------------------------------------

%% Argument check
if nargin < 3
  error('Insufficient input arguments!');
end

%% First get number of non-zero spike times from input spike time matrix M
rows = size(M,1);
N = zeros(1,rows);
for i = 1:rows
  N(i) = length(find(M(i,:)));
end

%% Create concatenated row vectors of lambar and dlambar and make plot
veclen = sum(N);
lambar_vec = Inf*ones(1,veclen);
dlambar_vec = Inf*ones(1,veclen);
i0 = 1;
for i = 1:rows
  i1 = i0 + N(i) - 1;
  lambar_vec(i0:i1) = lambar(i,1:N(i));
  dlambar_vec(i0:i1) = dlambar(i,1:N(i));
  i0 = i1 + 1;
end

figure
plot(lambar_vec,dlambar_vec,'.');