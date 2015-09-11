function data = isi_partition(M,Mb,Mlam,Mdlam,showplot)
%--------------------------------------------------------------------------------------------------------------------------
% Partition spike times into 21 subsets by their associated rate and rate derivative (time) values.  In the 2D
%  space of rate(t) vs. d/dt(rate(t)), we partition the plane such that there are an equal number of spikes (within
%  1 or 2) in each box.  This segments the spike times into regimes of low firing rate, medium firing rate, high firing rate,
%  and so forth, and also segments them by regimes of slow, medium, and fast rate changes (in time).
%
% USAGE:            data = isi_partition(M,Mb,lam,dlam,showplot);
% INPUT:            M                       * Matrix of spike times in original time (each row is a trial, zero padded)
%                   Mb                      * Matrix of spike times (transformed)
%                   Mlam                    * Matrix of rate values at each spike time (original time) (zero padded)
%                   Mdlam                   * Matrix of time derivative of rate values at each spike time (original time)
%                   showplot                * (logical or int) flag to show plots or not (default true)
% OUTPUT:           data                    * (struct) output data struct; includes 
%
% Dependencies:   get_bigspikematrix.m, get_isi_each_row.m, plot_lam_dlam.m
%
% Comments:
%  (1) The data field "data.Bseg" contains the array of interval data, for which each cell entry
%       is an 8-column data matrix with the following columns:
%      col 1 = trial #
%      col 2 = spike # in trial
%      col 3 = spike time t (original time)
%      col 4 = spike time B (time B spike time)
%      col 5 = lambda value for spike (rate function, original time)
%      col 6 = d/dt(lambda) for the spike 
%      col 7 = ISI associated with spike (original time)
%      col 8 = ISI associated with spike (time B)
%
% Written by Alex Casti, MSSM, Dec 2007
% Last modified March 24, 2008
%--------------------------------------------------------------------------------------------------------------------------

if nargin < 4
  help isi_partition.m;
  error('Not enough inputs!');
end
if (nargin < 5) || isempty(showplot)
  showplot = true;    % show plots by default
end

% Get data matrix (spike times, trial numbers, spike intervals, etc)
spikedata = get_bigspikematrix(M,Mb,Mlam,Mdlam);
B = spikedata.B;      
numspikes = size(B,1);   % total number of spikes in data set

% Sort the rate vector (and permute rows of spike data matrix B accordingly)
lamvec = B(:,5);       % 5th column of B corresponds to rate data
[lamvec,isort] = sort(lamvec);
Bsortlam = B(isort,:);

%-------------------------------------------------------------------------------------------------------------
%                                        Partition the data
%-------------------------------------------------------------------------------------------------------------
% First partition into 7 segments along the x-axis (lam, or rate)
%  and then partition those data sets into 3 segments along the y-axis (dlam, or rate derivative).
% Note: must first partition data using the SORTED lambda vector (rates).
% NOTE: Make number of partitions a variable in future version.
%-------------------------------------------------------------------------------------------------------------

% Partition along x axis (lambda = rate)
numsegx = 7;
nx = floor(length(lamvec)/numsegx);   % number of spikes corresponding to each segment (last segment may have slop)
% Now find the positions of the level lines cutting through the x axis
Blamseg = cell(1,numsegx);     % spike data matrices sorted by lam segment (rate axis)
xgrid = zeros(1,numsegx);
for i = 1:numsegx 
  % Get segment indices for sorted lam values
  if i == numsegx
    ind = 1+(i-1)*nx:numspikes;  % the last segment may have fewer than 'nx' spikes
  else
    ind = 1+(i-1)*nx:i*nx;
  end
  Blamseg{i} = Bsortlam(ind,:);
  xgrid(i) = Blamseg{i}(1,5) - eps;  %  xgrid (vertical lines demarcating x-axis segments)
  % Note: the 5th column of "Blamseg" was chosen because it contains the rate value data for each spike
end

% Partition along y axis (lambda rate derivative)
numsegy = 3;
ygrid = zeros(numsegx,numsegy);  % y-grid changes depending on x value (hence 2 dimensional)
numspikes_seg = zeros(numsegx,numsegy);
Bseg = cell(numsegx,numsegy);     % spike data matrices sorted in both directions
Blamseg_dlamsort = cell(1,numsegx);
for i = 1:numsegx
  % Sort the dlam (lambda rate) values
  dlamvals = Blamseg{i}(:,6);
  [dlamvals,idsort] = sort(dlamvals);
	Blamseg_dlamsort{i} = Blamseg{i}(idsort,:);     % segment of data sorted along y-axis (rate derivative) direction
	%numspikes_seg = length(dlamvals);
	NN = length(dlamvals);  % Number of spikes in the y-direction partition (rate derivative)
	ny = floor(NN/numsegy);
  for j = 1:numsegy
		if j == numsegy
			indj = 1 + (j-1)*ny:NN;
		else
			indj = 1 + (j-1)*ny:j*ny;
		end
		Bseg{i,j} = Blamseg_dlamsort{i}(indj,:);
    numspikes_seg(i,j) = length(Bseg{i,j}); % Number of spikes associated with each segment (partition)
		ygrid(i,j) = Bseg{i,j}(1,6) - eps;  % column 6 of Bseg{i,j} has rate derivative information
  end
end

%-------------------------------------------------------------------------------------------------
%                                             PLOT DATA
%-------------------------------------------------------------------------------------------------

%-----------------------------------------
%      SPIKE TIME SCATTER PLOTS
%-----------------------------------------
% Plot cloud of data points (lam vs dlam) and grid lines (original time)
if showplot
  plot_lam_dlam(M,Mlam,Mdlam); 
    ax = axis(gca); hold on;
    xlabel('\lambda(t)','fontweight','bold','fontsize',13);
    ylabel('d\lambda/dt','fontweight','bold','fontsize',13);
    plot([xgrid; xgrid]',ax(3:4),'r');
    % Plot the hortizontal grid lines one at a time
    for i = 1:numsegx
      if i==numsegx
        xrange = [xgrid(end),ax(2)];
      else
        xrange = xgrid(i:i+1);
      end
      plot(xrange,[ygrid(i,:);ygrid(i,:)],'k');
    end
end
%-----------------------------------------
%         ISI HISTOGRAM PLOTS
%-----------------------------------------
% Plot the interval histograms (new time).  The number of bins used is motivated by a proposal
%  (for exploratory work) mentioned in "Understanding Robust and Exploratory Data Analysis"
%  edited by Hoaglin, Mosteller, and Tukey (p.29).  The "number of events" parameter is set as
%  the mean number of spikes within each segment.
%numbins = ceil(10*log10(mean(numspikes_seg(:))));
%numbins = 32  % fix the number of bins for purposes of commonality across cells, stim, etc
% Define a common ISI axis (new time)
minISI = 0;
maxISI = max(B(:,8));
binwidth = 0.275072878909593;  % fix the bin width
%edges = linspace(minISI,maxISI,numbins+1); 
edges = minISI:binwidth:maxISI;
numbins = length(edges)-1;
dISI = mean(diff(edges));
bincenters = edges(1:end-1)+0.5*dISI;
Nhist = cell(numsegx,numsegy);
Nhist_normalized = cell(numsegx,numsegy);
maxprob = 0;  % Calculate maximum probability in any window for axis scaling
xmin = 0;         % Minimum value on x-axis (ISI) for axis scaling
xmax = maxISI;    % Maximum value on x-axis (ISI) for axis scaling
% Note: Still will generate ISI histogram data even if showplot=false
if showplot, figure; end;
for i = 1:numsegx
  for j = 1:numsegy
    Nhist_temp = histc(Bseg{i,j}(:,8),edges); Nhist{i,j} = Nhist_temp(1:end-1);
    %Nhist_normalized{i,j} = Nhist{i,j}/sum(Nhist{i,j});  % Sum of prob is unity (alternate normalization)
    % Normalize each histogram to have unit area 
    Nhist_normalized{i,j} = (1/dISI)*Nhist{i,j}/sum(Nhist{i,j});  % Integral of prob is unity
    maxprob = max(maxprob,max(Nhist_normalized{i,j}));
    % Order plot so that lower left corner is lowest left box in lam vs dlam plane, etc
    k = i + 2*numsegx - numsegx*(j-1);  % This formula does the trick
    if showplot
      subplot(numsegy,numsegx,k)
        %plot(edges_plot,Nhist_normalized{i,j});
        stairs(bincenters,Nhist_normalized{i,j},'linewidth',2);
        title(sprintf('(%d,%d)',i,j),'fontweight','bold');
        axis([xmin xmax 0 1.5]);
        if j == 1
          xlabel('ISI','fontweight','bold');
        end
        if mod(k-1,numsegx)==0
          %ylabel('norm prob(count)','fontweight','bold');
          ylabel('ISI density','fontweight','bold');
        else 
          set(gca,'Yticklabel','');
        end
    end % end of showplot conditional
  end
end
% Give option to re-scale the axes of the figure
if showplot
  fprintf('maxprob = %g\n',maxprob);
  rescaleY_answer = input('Re-scale y axes using maxprob?  (0) No  (1) Yes : ');
  if rescaleY_answer
    maxprob_plot = maxprob;
    rescale = true;
  else  % Give user option to define own y-axis scale maximum
    rescaleY_answer = input('Re-scale the plot axes using your own vertical axis maximum?  (0) No  or  Enter Value : ');
    if rescaleY_answer
      rescale = true;
      maxprob_plot = rescaleY_answer;
    else 
      maxprob_plot = [];  % In this case each y-axis is scaled by Matlab automatically and is not defined by user
    end
  end
  fprintf('Maximum ISI = %g\n',maxISI);
  rescaleX_answer = input('Re-scale x axis maximum (max ISI is default)? (0) No  or  Enter Value : ');
  % xmin=0 is always assumed (default value is maximum ISI across all histograms)
  if rescaleX_answer
    xmax = rescaleX_answer;
    rescale = true;
  end
  if rescale   % Either x-axis or y=axis re-scale option was chosen
    for i =1:numsegx
      for j = 1:numsegy
         k = i + 2*numsegx - numsegx*(j-1);
         subplot(numsegy,numsegx,k)
           axis([xmin xmax 0 maxprob_plot]);  % No multiplicative factor here
      end
    end
  end
end  % end of showplot conditional

%-----------------------------------------
%         Define some output data
%-----------------------------------------
% Spike segmentation data
data.B = B;   % big data matrix containing spikes, rates, rate derivatives, intervals for orig and new time
data.Blamseg = Blamseg;  % spike times (orig time) segmented along x-axis only (rate axis)
data.Bseg = Bseg;        % spike times (orig time) associated with final 2D segmentation
data.numspikes_seg = numspikes_seg;  % number of spikes in each segment
data.numsegx = numsegx;              % number of segments in x direction (rate axis)
data.numsegy = numsegy;              % number of segments in y direction (rate derivative axis)
data.xgrid = xgrid;                  % segment boundaries (not including x=0) along rate (x) axis
data.ygrid = ygrid;                  % segment boundaries along rate derivative (y) axis
% Binning data for ISI histogram plots
data.numbins = numbins;                    % number of bins for each ISI histogram
data.edges = edges;                        % ISI egdges used to contruct bin counts 
data.bincenters = bincenters;              % bin center values used for histogram plots
data.Nhist = Nhist;                        % counts in each bin
data.Nhist_normalized = Nhist_normalized;  % normalized counts (such that histogram integrates to 1)
data.minISI = minISI;                      % minimum ISI (new time) over all ISIs
data.maxISI = maxISI;                      % maximum ISI (new time) over all ISIs
data.maxprob = maxprob;                    % maximum ISI density probability in any bin (across all partitioned ISI densities)
