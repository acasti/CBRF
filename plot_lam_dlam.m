function plot_lam_dlam(M,Mlam,Mdlam)
%---------------------------------------------------------------------------------------------------
% Plot lam (rate parameter) vs dlam (rate derivative) at each spike time.
%
% USAGE:     plot_lam_dlam(M,Mlam,Mdlam);
% INPUT:     M               * matrix of spike times (so we know number of events; note zero padding)
%            Mlam            * matrix of rate parameters for each spike
%            Mdlam           * matrix of time derivative of rates for each spike
% 
% Written by Alex Casti, MSSM, November 2007
%---------------------------------------------------------------------------------------------------

% First get number of non-zero spike times from input spike time matrix M
rows = size(M,1);
N = zeros(1,rows);
for i = 1:rows
  N(i) = length(find(M(i,:)));
end

% Create concatenated row vectors of lambar and dlambar
veclen = sum(N);
lam_vec = Inf*ones(1,veclen);
dlam_vec = Inf*ones(1,veclen);
i0 = 1;
for i = 1:rows
  i1 = i0 + N(i) - 1;
  lam_vec(i0:i1) = Mlam(i,1:N(i));
  dlam_vec(i0:i1) = Mdlam(i,1:N(i));
  i0 = i1 + 1;
end

figure
plot(lam_vec,dlam_vec,'.','MarkerSize',1.5);