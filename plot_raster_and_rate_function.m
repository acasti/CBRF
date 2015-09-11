function plot_raster_and_rate_function(M,t,rate,lamgrid,tzoom)
%-----------------------------------------------------------------------------------------------
% Make raster plot of spike times and plot the rate function below it.
%
% USAGE:    plot_raster_and_rate_function(M,t,rate,lamgrid,tzoom);
% INPUT:    M                 * matrix of spike times (row=trial, zero padded)
%           t                 * time grid over which rate is defined
%           rate              * rate function r(t)
%           lamgrid           * (optional) plot horizontal line at these values 
%                                  in rate plot
%           tzoom             * (optional) make a second zoom plot 
%                                 in the time range [tzoom(1),tzoom(2)]
%
% Written by Alex Casti, MSSM 17 Jan 2008
%-----------------------------------------------------------------------------------------------
if nargin < 3
  error('Not enough inputs!');
end
if ~iscell(M)
  spikecell = spiketime_mat2cell(M);
else
  spikecell = M;
end
if nargin < 4 || isempty(lamgrid)
  lamgrid = [];
  plot_lamgrid_lines = false;
else
  plot_lamgrid_lines = true;
end
if nargin < 5 || isempty(tzoom)
  make_zoomplot = false;
else
  make_zoomplot = true;
end
if (isempty(tzoom) || ~exist('tzoom','var')) || length(tzoom)~=2
  tzoom = zeros(1,2);
  fprintf('You have chosen to make a zoom plot but didn''t specify the range!\n');
  tzoom(1) = input('Enter lower bound of zoom plot.  t0 = : ');
  tzoom(2) = input('Enter upper bound of zoom plot.  t1 = : ');
end

figure
t0 = t(1); t1 = t(end);  % First plot full range
subplot(2,1,1)
  rasterplot(spikecell,[t0,t1],[],[],0); hold;
  axis([t0 t1 0 length(spikecell)+1]);
  xlabel(' ');
  ylabel('trial','fontweight','bold','fontsize',12);
  title('Spike Raster Plot','fontweight','bold','fontsize',12);
subplot(2,1,2)
  plot(t,rate,'r-'); hold;
  if plot_lamgrid_lines
    for i = 1:length(lamgrid);
      plot([t0,t1],lamgrid(i)*ones(1,2),'k-');
    end
  end
  xlabel('time (sec)','fontweight','bold','fontsize',12);
  ylabel('Firing Rate (ips)','fontweight','bold','fontsize',12);
  axis([t0 t1 0 1.1*max(rate)]);
  tit = sprintf('Cosine Bell rate function');
  title(tit,'fontweight','bold','fontsize',12);

if make_zoomplot
  figure
  subplot(2,1,1)
    rasterplot(spikecell,[t0,t1],[],[],0); hold;
    axis([tzoom(1) tzoom(2) 0 length(spikecell)+1]);
    xlabel(' ');
    ylabel('trial','fontweight','bold','fontsize',12);
    title('Spike Raster Plot','fontweight','bold','fontsize',12);
  subplot(2,1,2)
    plot(t,rate,'r-'); hold;
    if plot_lamgrid_lines
      for i = 1:length(lamgrid);
        plot([t0,t1],lamgrid(i)*ones(1,2),'k-');
      end
    end
    xlabel('time (sec)','fontweight','bold','fontsize',12);
    ylabel('Firing Rate (ips)','fontweight','bold','fontsize',12);
    axis([tzoom(1) tzoom(2) 0 1.1*max(rate)]);
    tit = sprintf('Cosine Bell rate function');
    title(tit,'fontweight','bold','fontsize',12);
end
