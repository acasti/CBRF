function plot_raster_and_rate_function2(M1,M2,t,rate1,rate2,tzoom,lamgrid)
%-----------------------------------------------------------------------------------------------
% Make raster plot of 2 sets of spike times and plot the rate functions below it.
%
% USAGE:    plot_raster_and_rate_function2(M1,M2,t,rate1,rate2,tzoom,lamgrid);
% INPUT:    M1                 * matrix 1 of spike times (row=trial, zero padded)
%           M2                 * matrix 2 of spike times
%           t                  * time grid over which rate is defined
%           rate1              * first rate function r(t)
%           rate2              * second rate function r(t)
%           tzoom              * (optional) make a second zoom plot 
%                                 in the time range [tzoom(1),tzoom(2)]  (default false)
%           lamgrid            * (optional) plot horizontal line at these values 
%                                  in rate plot (default false)
%
% Written by Alex Casti, MSSM 25 Feb 2008
%-----------------------------------------------------------------------------------------------
if nargin < 5
  error('Not enough inputs!');
end
if ~iscell(M1) 
  spikecell1 = spiketime_mat2cell(M1);
else
  spikecell1 = M1;
end
if ~iscell(M2) 
  spikecell2 = spiketime_mat2cell(M2);
else
  spikecell2 = M2;
end
if nargin < 6 || isempty(tzoom)
  make_zoomplot = false;
else
  make_zoomplot = true;
end
% if (isempty(tzoom) || ~exist('tzoom','var')) || length(tzoom)~=2
%   tzoom = zeros(1,2);
%   fprintf('You have chosen to make a zoom plot but didn''t specify the range!\n');
%   tzoom(1) = input('Enter lower bound of zoom plot.  t0 = : ');
%   tzoom(2) = input('Enter upper bound of zoom plot.  t1 = : ');
% end
if nargin < 7 || isempty(lamgrid)
  lamgrid = [];
  plot_lamgrid_lines = false;
else
  plot_lamgrid_lines = true;
end

figure
t0 = t(1); t1 = t(end);  % First plot full range
subplot(2,1,1)
  rasterplot(spikecell1,[t0,t1],[1 0 0],[],0); hold;  % red dots
  rasterplot(spikecell2,[t0,t1],[0 0 1],[],0);        % blue dots
  axis([t0 t1 0 length(spikecell1)+1]);
  xlabel(' ');
  ylabel('trial','fontweight','bold','fontsize',12);
  title('Spike Raster Plot','fontweight','bold','fontsize',12);
subplot(2,1,2)
  plot(t,rate1,'r-'); hold;
  plot(t,rate2,'b-');
  if plot_lamgrid_lines
    for i = 1:length(lamgrid);
      plot([t0,t1],lamgrid(i)*ones(1,2),'k-');
    end
  end
  xlabel('time (sec)','fontweight','bold','fontsize',12);
  ylabel('Firing Rate (ips)','fontweight','bold','fontsize',12);
  maxrate = max( max(rate1),max(rate2) );
  axis([t0 t1 0 1.1*maxrate]);
  tit = sprintf('Cosine Bell rate function');
  title(tit,'fontweight','bold','fontsize',12);
  legend('rate 1','rate 2');
if make_zoomplot
  figure
  subplot(2,1,1)
    rasterplot(spikecell1,[t0,t1],[1 0 0],[],0); hold;
    rasterplot(spikecell2,[t0,t1],[0 0 1],[],0);
    axis([tzoom(1) tzoom(2) 0 length(spikecell1)+1]);
    xlabel(' ');
    ylabel('trial','fontweight','bold','fontsize',12);
    title('Spike Raster Plot','fontweight','bold','fontsize',12);
  subplot(2,1,2)
    plot(t,rate1,'r-'); hold;
    plot(t,rate2,'b-');
    if plot_lamgrid_lines
      for i = 1:length(lamgrid);
        plot([t0,t1],lamgrid(i)*ones(1,2),'k-');
      end
    end
    xlabel('time (sec)','fontweight','bold','fontsize',12);
    ylabel('Firing Rate (ips)','fontweight','bold','fontsize',12);
    maxrate = max( max(rate1),max(rate2) );
    axis([tzoom(1) tzoom(2) 0 1.1*maxrate]);
    tit = sprintf('Cosine Bell rate function');
    title(tit,'fontweight','bold','fontsize',12);
    legend('rate 1','rate 2');
end
