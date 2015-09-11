function dataA = get_timeA_cosbells2(M, trange, dt, options)
%---------------------------------------------------------------------------------------------------------------
% Replaces each spike time in a given data set S with a cosine bell function whose width is determined
%   by the maximum of the two adjacent spike intervals.
%
% USAGE:       dataA = get_timeA_cosbells2(M, trange, dt, options);
% INPUT:       M                     * spike data (matrix, with each row a trial of data with zero padding) 
%              trange                * (vector) [t0,t1]  time range over which to determine rate function
%              dt                    * (double) sample time interval (fixed)
%              options               * (struct) optional options structure, which can accept the following fields:
%                     .randflag        (logical) add pseudo-random number to the 1/100th msec decimal place or not
%                     .show_progress   (logical) output progress of calculations to Matlab command window
%                     .method          (string)  'exact' (default) or 'interp' (use quicker interpolation)
% OUTPUT:      dataA                 * (struct)  output data struct with timeA data
% 
% Notes:
%    (1) External functions used:   spiketime_cell2mat.m, locate.m, cumsimps.m, get_pseudo_randnums_for_spiketimes.m
%    (2) All other non-standard Matlab function calls are defined within this routine
%    (3) If input data matrix M is a cell array then it is automatically converted to a zero padded matrix.
%    (4) Fixed problem that can occur if some trials have zero spikes, or if the maximum number of spikes per trial
%        is just 1.
%    (5) Added options structure.  Can now choose to do an exact integration (options.method='exact') of the cosine
%        bells to obtain transformed spike times, or use a faster approximation method (options.method='interp') to
%        do the necessary integrals of the time transformation (accurate to second order in grid spacing).
%
% Written by Alex Casti, MSSM, December 2007
% Last modified 02 Sept 2015
%---------------------------------------------------------------------------------------------------------------

%% Input argument check
if nargin < 3
	error('Not enough inputs!');
end
if ~isnumeric(M)  % Case: M given as cell array; return as zero-padded matrix with each row a trial
	M = spiketime_cell2mat(M);    
end
if any([~isnumeric(trange),length(trange)~=2,trange(end)<=trange(1)])
	error('trange =[t0,t1] must be a vector of length 2 with t0 < t1');
end
if nargin < 4 || isempty(options)
  options = struct;  % define a variable for default struct
end
% Check options structure and set defaults if some fields are undefined
options = check_options(options);

dataA.M = M;                   % input spike times (before any whitening, if randflag = true)

% Define the time grid (constant spacing) on which we define the various functions
t = trange(1):dt:trange(2); 

%---------------------------------------------------------------------------------------------------------------
%                               Spike time data whitening
%---------------------------------------------------------------------------------------------------------------
% Add a pseudo-random number (PRN) to each spike time if instructed to do so.  The idea is to make spike times
%  across different trials unique if they were not collected with sufficient precision to distinguish them.
if options.randflag
	fprintf('Adding pseudo-random numbers to each spike time...\n');
  [Mwhite,prn] = add_PRNs_to_spiketimes(M);
	dataA.Mwhite = Mwhite; dataA.prn = prn;
  M = Mwhite; clear Mwhite;   % redefine M to be the whitened data
end

%---------------------------------------------------------------------------------------------------------------
%                              Cosine Bells definition (and support indices)
%---------------------------------------------------------------------------------------------------------------
% Get cosine bell widths (time A) for each spike in each trial.  
MwidthsA = get_cosbell_widthA(M);     % Gets cosine bell widths for each spike in each trial (returns matrix MwidthsA)

% Define cosine bells and associated time intervals of support (region over which the cosine bell is non-zero) 
%  using the specified grid.
fprintf('Calculating time A cosine bells and time support over the chosen grid...\n');
fprintf('t = [%g,%g] , dt = %g\n',t(1),t(end),mean(diff(t)));
[SUPP, lamA_mean, dlamA_mean, tA, numspikes_trial] = get_ratesA_on_grid(t,M,MwidthsA,options);
% Get lamA (rate A) at each spike time (interpolate on the grid)
fprintf('Calculating rates and timeA at each spike time specifically...\n');
if strcmp(options.method,'exact')
  fprintf('Using exact method for time A transformation...\n');
  [MlamA,MdlamA,MA] = get_lamA_timeA_each_spike_exact(M,SUPP,options);
elseif strcmp(options.method,'interp')
  fprintf('Using interpolation method for time A transformation...\n');
  [MlamA, MdlamA, MA] = get_lamA_timeA_each_spike_interp(M,t,tA,lamA_mean,dlamA_mean,options); 
else
  fprintf('None of the allowed methods were specified in the options structure!\n');
  fprintf('Check the original Matlab function (a default option should have been set!\n');
  return;
end

% Define some output variables
dataA.t = t;        % original time on the grid
dataA.tA = tA;      % new time A values on the original grid t specifically (same length as t)
dataA.dt = dt;      % grid spacing (original time)
dataA.lamA_mean = lamA_mean; % rate function A using cosine bells A
dataA.dlamA_mean = dlamA_mean;  % rate derivative function A
dataA.MA = MA;                       % matrix of time A values associated with each spike
dataA.MlamA = MlamA;              % matrix of rate function A values at each spike time (same dimension as M) 
dataA.MdlamA = MdlamA;           % matrix of rate derivative function A values at each spike time
dataA.MwidthsA = MwidthsA;     % matrix of cosine bell widths A for each spike (same dimension as input M)
dataA.SUPP = SUPP;                 % matrix containing support (on grid) information for each spike
dataA.numspikes_trial = numspikes_trial;    % number of spikes per trial
dataA.options = options;          % options used in time A construction
%********************************************************************************
%                                             Function definitions
%********************************************************************************

%-------------------------------------------------------------------------------------------------------
function [Mwhite,prn] = add_PRNs_to_spiketimes(M)
%-------------------------------------------------------------------------------------------------------
% Add pseudo-random numbers (PRNs) to the 1/100th msec decimal place to each of the spike times
%  in order to whiten the data and eliminate artifacts of duplicate spike times if they were
%  stored to a somewhat low precision (i.e. to a common precision of 1/10th of a millisecond).
% NOTE: spike times assumed in seconds
	seedprime = 1049;  % 4 digit prime will generate 10^4-1 = 9999 non-repeated pseudo-random values
	decplaces = 4;     % adds PRNs after the 4th decimal place (i.e. at 1/100th msec precision)
  offset = (1/1000)*(1/20);  % subtract this offset from the PRNs (some spikes add jitter, some subtract)
	[~, trandscale] = get_pseudo_randnums_for_spiketimes(seedprime,decplaces);
  trandscale = trandscale - offset;   % trandscale is already multiplied by appropriate power of 10
  numspikes = length(find(M(:)));
	% stack the PRN vector so that it equals the length of the spike train
	if length(trandscale) < numspikes
	  prn = repmat(trandscale,1,ceil(numspikes/length(trandscale)));
  	prn = prn(1:numspikes);         % The vector of pseudo-random numbers added 
  else 
    prn = trandscale(1:numspikes)';  % generated more PRNs than spikes, so just truncate vector
  end
  [v,numspikes_trial] = spikematrix2vec(M); % Get vector v of non-zero spike times in M
  v = v + prn;                                  % Whiten the spike data
  Mwhite = spikevec2matrix(v,numspikes_trial);  % Get the whitened data back in matrix form (numspikes_trial is crucial!)

%-------------------------------------------------------------------------------------------------------
function MwidthsA = get_cosbell_widthA(M)
%-------------------------------------------------------------------------------------------------------
% Get cosine bell widths (time A) for the spike train.  Each width is the maximum of the two spike intervals
% adjacent to each spike in each trial.
MwidthsA = zeros(size(M));  % If no spikes in trial then all "bell widths" are simply zero
for i = 1:size(MwidthsA,1)
	temp = M(i,:); spikes = temp(temp>0);   % get non-zero spikes in each row
  if length(spikes) == 1  % Just 1 spike in trial
    MwidthsA(i,1) = spikes(1);
  end
  if length(spikes) > 1
    ISI = diff(spikes);
    %MwidthsA(i,1) = ISI(1);    % first width is the lone interval associated with first spike
    MwidthsA(i,1) = max(ISI(1),M(i,1));   % this version is consistent with Prasad's old code
    MwidthsA(i,2:length(spikes)) = max([ISI; [ISI(2:end),-Inf]]); 
    % The -Inf above ensures that the last spike's width is associated with the lone ISI
  end
end

%-------------------------------------------------------------------------------------------------------
function c = cosbell(tval,t0,width)
%-------------------------------------------------------------------------------------------------------
% Cosine bell function.  This version accepts vector input for all input parameters, which is useful when we 
%  do summations over support spikes associated with a fixed spike time "tval"
  c = (1./(2*width)).*(1+cos(pi*(tval-t0)./width));  
  
	
%-------------------------------------------------------------------------------------------------------
function dc = dcosbell(tvec,t0,width)
%-------------------------------------------------------------------------------------------------------
% Derivative of cosine bell function (also see comments for cosbell function)
	dc = -(pi/2).*((1./width).^2).*sin(pi*(tvec-t0)./width);

%-------------------------------------------------------------------------------------------------------
function intvals = integral_over_interval(a,b,t0,width)
%-------------------------------------------------------------------------------------------------------
% Integral of rate function between times [a,b] 
	intvals = 0.5*((b-a)./width) + (0.5/pi)*sin(pi*(b-t0)./width) - (0.5/pi)*sin(pi*(a-t0)./width);

%-------------------------------------------------------------------------------------------------------
function [SUPP, lamA_mean, dlamA_mean, tA, numspikes_trial] = ...
    get_ratesA_on_grid(t,M,MwidthsA,options)
%-------------------------------------------------------------------------------------------------------
% Get cosine bells for each spike and define the time vector indices of support for each over the chosen
%  time grid.  This is NOT the function that provides the exact calculation at each spike time.  This defines
%  the rate and rate derivative over the grid points only; any spike time likely does not fall on the grid precisely.
% This function also performs the mapping from the original time t  to the new time A  on the grid.
	
  dt = mean(diff(t));
  [~,numspikes_trial] = spikematrix2vec(M);  % need number of spikes per trial (N) to set cell array size
  numspikes = sum(numspikes_trial);          % total number of spikes across all trials
  numtrials = size(M,1);                     % number of trials (repeats)
  lamA_mean = zeros(1,length(t));            % mean lamA (averaged across trials)
  dlamA_mean = zeros(1,length(t));           % mean lamA (averaged across trials)
  SUPP = zeros(numspikes,6);   % matrix of each spike's support information [i j is0 is1 tspk widthA]
	                             % i = trial#
                               % j = spike# in trial
															 % is0 = leftmost support grid index
                               % is1 = rightmost support grid index
															 % tspk = actual spike time
                               % widthA = cosine bell width for spike
	for i = 1:numtrials
    if options.show_progress
      if mod(i,10)==0, fprintf('Time A rate function: processing trial %d of %d\n',i,numtrials); end;
    end
    for j = 1:numspikes_trial(i)
      spikenumber = sum(numspikes_trial(1:i-1)) + j;
      tspk = M(i,j);      % jth spike time (original time) in trial i
      w = MwidthsA(i,j);  % width of associated cosine bell
      tmin = tspk - w; tmax = tspk + w;
      % Check for case where cosine bell terminates within a single grid spacing (better to choose smaller dt!)
      if 2*w < dt   % note that w really spaces 1/2 the support of the cosine bell
        fprintf('For spike %d of trial %d the cosine bell support is smaller than one grid interval!\n',j,i);
        error('Choose a finer time mesh!');
			else
				% Get time indices of support for cosine bell
        is0 = max(locate(t,tmin),1);  % Use max(*,1) in case cosine bell drops to zero off to left of grid 
        is1 = locate(t,tmax);
      end
      if t(is0) < tmin
        is0 = is0 + 1;  % 'locate' function may return index one too far to the left
      end
      SUPP(spikenumber,:) = [i j is0 is1 tspk w];
      % Error check: make sure that all spikes fall on the grid somewhere (modify code later to handle arbitrary grid subsets)
      if isempty(is0) || isempty(is1)
        fprintf('Bad spike processing for spike %d of trial %d\n',j,i);
        error('Some spikes were out of range of the specified time grid');
      end
      indS = is0:is1;   % range of time vector indices which the cosine bell for this spike spans
      lamA_mean(indS) = lamA_mean(indS) + cosbell(t(indS),tspk,w);
      dlamA_mean(indS) = dlamA_mean(indS) + dcosbell(t(indS),tspk,w);
    end
  end
  % Now divide lamA_mean by the number of trials to get the true mean (ensemble average, effectively)
  lamA_mean = lamA_mean/numtrials;
  dlamA_mean = dlamA_mean/numtrials;
	% Now map the original times  't'  on grid into time A  'tA'  values (just the cumumlative integral of the rate function)
	tA = cumsimps(t,lamA_mean);  % 4th-order accurate Simpson rule
	
%-------------------------------------------------------------------------------------------------------
function [MlamA, MdlamA, MA] = get_lamA_timeA_each_spike_exact(M,SUPP,options)
%-------------------------------------------------------------------------------------------------------
% Use exact integrations to get transformed "time A" spike times
  MlamA = zeros(size(M));
  MdlamA = zeros(size(M));
	MA = zeros(size(M));
	[~,numspikes_trial] = spikematrix2vec(M);  % just need number of spikes per trial 
	numspikes = sum(numspikes_trial);
  numtrials = size(M,1);       % number of trials (repeats)
	% Sort the SUPP matrix in order of increasing spike time (SUPP keeps track of trial and spike# in trial)
	spikes_all = SUPP(:,5);
	[spikes_all_sort, ispiketime_sort] = sort(spikes_all);
	SUPPSORT = SUPP(ispiketime_sort,:);  % SUPP matrix sorted in order of increasing spike time
  wvals_all_sort = SUPPSORT(:,6);
  tsupp_left = spikes_all_sort - wvals_all_sort;   % tspk - widthA  for each spike
  tsupp_right = spikes_all_sort + wvals_all_sort;  % tspk + widthA for each spike
  
  for nspike = 1:numspikes
		i = SUPPSORT(nspike,1);   % associated trial# and spike# in trial for the spike
		j = SUPPSORT(nspike,2);
		if options.show_progress
			if mod(nspike,1000) == 0
				fprintf('Processing spike %d of %d\n',nspike,numspikes);
			end
		end
    %tspk = SUPPSORT(nspike,5);   % column 5 has spike times
    tspk = spikes_all_sort(nspike);
		isupp = get_support_indices_single_spike(tspk,tsupp_left,tsupp_right);
		tcentervals = SUPPSORT(isupp,5);           % spike times associated with support of tspk
		wvals = SUPPSORT(isupp,6);                 % widths of spikes associated with support of tspk
		c = cosbell(tspk,tcentervals,wvals);
		dc = dcosbell(tspk,tcentervals,wvals);
		MlamA(i,j) = sum(c)/numtrials;
		MdlamA(i,j) = sum(dc)/numtrials;
		% Calculate spikes in new time A, keeping spike identity within a trial as we step from one time to the next.
		if nspike == 1
			tspk_previous = 0;
			tA_previous = 0;
      [avals,bvals,wvals,tcentervals] = ...
        get_support_ranges_for_integral(tspk_previous,tspk,tsupp_left,tsupp_right,spikes_all_sort,wvals_all_sort);
			intvalsA = integral_over_interval(avals, bvals, tcentervals, wvals);
			tA = tA_previous + sum(intvalsA)/numtrials;
			MA(i,j) = tA;
		end
		if nspike > 1
      [avals,bvals,wvals,tcentervals] = ...
        get_support_ranges_for_integral(tspk_previous,tspk,tsupp_left,tsupp_right,spikes_all_sort,wvals_all_sort);
			intvalsA = integral_over_interval(avals, bvals, tcentervals, wvals);
			tA = tA_previous + sum(intvalsA)/numtrials;
			MA(i,j) = tA;
		end
		tspk_previous = tspk;
		tA_previous = tA;
  end
  
%-------------------------------------------------------------------------------------------------------
	function isupp = get_support_indices_single_spike(tspk,tsupp_left,tsupp_right)
%-------------------------------------------------------------------------------------------------------
% Given a spike time "tspk" return the indices (rows) of SUPP (matrix with support information over the grid, 
% which is not necessary for the exact calculations, incidentally) that correspond to spikes with support 
% including the time "tspk".
  vsum = sign(tspk - tsupp_left) + sign(tsupp_right - tspk);
  isupp = find(vsum > 0); 		% any entries with vsum > 0 constitude the support indices of tspk
    
%-------------------------------------------------------------------------------------------------------
	function [avals,bvals,wvals,tcentervals] = ...
      get_support_ranges_for_integral(tspk1,tspk2,tsupp_left,tsupp_right,S,W)
%-------------------------------------------------------------------------------------------------------
  % any entries with vsum = 2 constitude the support indices of tspk
  vsum1 = sign(tspk1 - tsupp_left) + sign(tsupp_right - tspk1);
  isupp1 = find(vsum1 > 0); 	% support indices for lower bound of integration
  vsum2 = sign(tspk2 - tsupp_left) + sign(tsupp_right - tspk2);
  isupp2 = find(vsum2 > 0);   % support indices for upper bound of integration
  isupp_all = union(isupp1,isupp2);  % all cos bells with support in the integration interval
  %**************************************************************************************
  % Determine appropriate integration range for each cosine bell with support in interval
  %**************************************************************************************
  spikes_supp_all = S(isupp_all);   % spike times associated with cos bell support in range
  widths_supp_all = W(isupp_all);   % all cosine bell widths associated with support spikes
  % Get "left" and "right" cosine bells (those with spike times <= tspk1 and > tspk1)
  check_left_right = spikes_supp_all - tspk1;
  ileft = find(check_left_right <= 0);        % indices associated with "left cosine bells"
  iright = setxor(ileft,1:length(isupp_all)); % "right bell" indices given by setxor operation
  spikes_left = spikes_supp_all(ileft);
  widths_left = widths_supp_all(ileft);
  spikes_right = spikes_supp_all(iright);
  widths_right = widths_supp_all(iright);
  aleft = tspk1*ones(length(ileft),1);
  bleft = min([spikes_left+widths_left , tspk2*ones(length(ileft),1) ],[],2);
  aright = max([spikes_right-widths_right , tspk1*ones(length(iright),1)],[],2);
  bright = tspk2*ones(length(iright),1);
  avals = [aleft;aright];
  bvals = [bleft;bright];
  wvals = [widths_left;widths_right];
  tcentervals = [spikes_left;spikes_right];

%-------------------------------------------------------------------------------------------------------
function [MlamA, MdlamA, MA] = get_lamA_timeA_each_spike_interp(M,t,tA,lamA_mean,dlamA_mean,options)
%-------------------------------------------------------------------------------------------------------
% Use approximation scheme to get transformed "time A" spike times (valid to 3rd order in dt) and 
%  rates and rate derivatives at the spike times.
  MlamA = zeros(size(M));
  MdlamA = zeros(size(M));
	MA = zeros(size(M));
	[~,numspikes_trial] = spikematrix2vec(M);  % just need number of spikes per trial
  numtrials = size(M,1);       % number of trials (repeats)
  % Get rate and d/dt(rate) values at the spike times using cubic Hermite interpolation
  %   First need cubic-Hermite interpolation structure for lamA_mean (rate function) and 
  %   dlamA_mean (rate function derivative) over the entire grid.
  pchip_structure_lamA = pchip(t,lamA_mean);
  pchip_structure_dlamA = pchip(t,dlamA_mean);
  for i = 1:numtrials  % cycle through trials
    if options.show_progress
      if mod(i,10) == 0
        fprintf('Processing trial %d of %d\n',i,numtrials);
      end
    end
    spikes_trial = M(i,1:numspikes_trial(i));
    MlamA(i,1:numspikes_trial(i)) = ppval(pchip_structure_lamA,spikes_trial);
    MdlamA(i,1:numspikes_trial(i)) = ppval(pchip_structure_dlamA,spikes_trial);
    for j = 1:numspikes_trial(i) % cycle through spikes to get interpolated "time A" values
      tspk = spikes_trial(j);
      ispk = locate(t,tspk);  % find time on grid just to the left of actual spike time
      tgrid = t(ispk);
      tgridA = tA(ispk);
      MA(i,j) = tgridA + lamA_mean(ispk)*(tspk-tgrid) + 0.5*dlamA_mean(ispk)*(tspk-tgrid)^2;
    end
  end

%-------------------------------------------------------------------------------------------------------
	function options_out = check_options(options_in)
%-------------------------------------------------------------------------------------------------------    
% Check the options structure and set defaults if necessary
options_out = options_in;
if ~isfield(options_in,'randflag')
  options_out.randflag = false;    % do not whiten data by default
end
if ~isfield(options_in,'show_progress')
  options_out.show_progress = true;    % show calculation progress data by default
end
if ~isfield(options_in,'method')
  options_out.method = 'exact';       % do exact integration of cosine bells by default
end