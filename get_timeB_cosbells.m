function dataB = get_timeB_cosbells(MA, tA, bwidth, show_progress)
%---------------------------------------------------------------------------------------------------------------
% Replaces each spike time in a given data set S with a cosine bell function whose width is determined
%   by a parameter "bwidth" that for each spike marks a time "bwidth" spikes to the right and "bwidth" spikes
%   to the left, and then takes the average of these two distances for the "width" of each cosine bell.
%
% USAGE:      dataB = get_timeB_cosbells(MA, tA, bwidth);
% INPUT:      MA                 * spike data (matrix, with each row a trial of data with zero padding; time A)
%             tA                 * (vector) time A values (mapped from original input time t grid)
%             bwidth             * (int) cosine bell width parameter for time B (# spikes to move each direction)
%             show_progress      * (logical) if true then output progress of computations to Matlab command window
% OUTPUT:     dataB              * (struct)  output data struct
% 
% Notes:
%    (1) External functions used:   spiketime_cell2mat.m, locate.m, cumsimps.m, get_spike_identification_matrix.m
%    (2) All other non-standard Matlab function calls are defined within this routine
%    (3) If input data matrix M is a cell array then it is automatically converted to a zero padded matrix.
%    (4) Fixed problem that can occur if some trials have no spikes.  There was also a problem
%      if the maximum number of spikes in any trial was just 1.  (FIXED)
%    (5) I verified that this fixed version of the function is 100% consistent with the older
%        version for data that doesn't suffer a paucity of spikes as mentioned above in (4).
%
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 13 September 2015
%---------------------------------------------------------------------------------------------------------------

% Error checking
if nargin < 3
	error('Not enough inputs!');
end
if (nargin < 4) || isempty(show_progress)
	show_progress = false;       % (default) do not show progress report of computation as it progresses
end
if ~isnumeric(MA)  % Case: M given as cell array; return as zero-padded matrix with each row a trial
	MA = spiketime_cell2mat(MA);    
end
dtA = mean(diff(tA));

%---------------------------------------------------------------------------------------------------------------
%                              Cosine Bells definition (and support indices)
%---------------------------------------------------------------------------------------------------------------
% Get cosine bell widths (time B) for each spike in each trial.  
MwidthsB = get_cosbell_widthB(MA,bwidth);     

% Define cosine bells and associated time intervals of support (region over which the cosine bell is non-zero) 
%  using the specified grid.
fprintf('Calculating time B cosine bells and time support over the input grid of time A times...\n');
fprintf('tA = [%g,%g] , dt = %g\n',tA(1),tA(end),dtA);
[SUPP, lamB_mean, dlamB_mean, tB, numspikes_trial] = get_ratesB_on_grid(tA,MA,MwidthsB,show_progress);
% Get lamB (rate B) at each spike time (interpolate on the grid)
fprintf('Calculating rates and timeB at each spike time specifically...\n');
[MlamB, MdlamB, MB] = get_lamB_timeB_each_spike(MA,SUPP,show_progress);

% Define some output variables
dataB.MA = MA;                          
dataB.tA = tA;         % original time (time A) on the grid
dataB.tB = tB;         % new time B on the grid (same length as tA)
dataB.dtA = dtA;      % grid spacing (original time)
dataB.lamB_mean = lamB_mean;     % rate function B using cosine bells B
dataB.dlamB_mean = dlamB_mean;  % rate derivative function B
dataB.MB = MB;                       % matrix of time B values associated with each spike
dataB.MlamB = MlamB;              % matrix of rate function B values at each spike time (same dimension as MA) 
dataB.MdlamB = MdlamB;           % matrix of rate derivative function B values at each spike time
dataB.MwidthsB = MwidthsB;     % matrix of cosine bell widths B for each spike (same dimension as input MA)
dataB.SUPP = SUPP;                 % matrix containing support (on grid) information for each spike
dataB.numspikes_trial = numspikes_trial;    % number of spikes per trial
dataB.bwidth = bwidth;

%********************************************************************************
%                                             Function definitions
%********************************************************************************

%-------------------------------------------------------------------------------------------------------
function MwidthsB = get_cosbell_widthB(MA,bwidth)
%-------------------------------------------------------------------------------------------------------
% Get cosine bell widths for the spike train.  Each width is the average of the distances in time "bwidth" 
%  spikes to the left and "bwidth" spikes to the right.  For endpoint effects, if there are not "bwidth" 
%  spikes to the left, for example, then the leftmost spike time on the grid is used to get the 
%  "left width," and the number of spikes missing (num_missing) on the left are added in the 
%  right-going search, so that the "right width" is given by the spike time 'bwidth+num_missing" spikes
%  to the right of the spike in question.  The "left width" and the "right width" are then 
%  averaged as usual.

MwidthsB = zeros(size(MA));
[BIGMA, numspikes_trial] = get_spike_identification_matrix(MA);
numspikes = sum(numspikes_trial);
% NOTE: For BIGMA: col1 = trial# i, col2 = spike# in trial j, col3 = spike time, col4 = ISI associated with spike
spikevec = BIGMA(:,3);
% need to sort the spikes first (note that we have effectively collapsed them onto a single time axis for time B)
[~,isort] = sort(spikevec);
BIGMA_sort = BIGMA(isort,:);
for nspike = 1:numspikes
  i = BIGMA_sort(nspike,1);       % trial# is first column of BIGMA_sort matrix
	j = BIGMA_sort(nspike,2);       % spike# in trial is second column of BIGMA_sort matrix
	tspk = BIGMA_sort(nspike,3);
	nleft = nspike - bwidth;       
	nright = nspike + bwidth;
	% make sure spike indices don't spill beyond minimum spike time
	if nleft < 1
		tleft = BIGMA_sort(1,3);   % just use the left-most spike time
    num_missing = bwidth - nspike + 1;
    nright = nright + num_missing;  % add the missing # of spikes to left to # spikes to the right
    % Note: right-edge spike time will automatically be adjusted below, so no need to set it here
  else
		tleft = BIGMA_sort(nleft,3);
  end
  % make sure spike indices don't spill beyond maximum spike time
	if nright > numspikes
		tright = BIGMA_sort(end,3);  % just use the right-most spike time
    num_missing = nspike + bwidth - numspikes;
    nleft = nleft - num_missing; % add the missing # of spikes to right to # spikes to the left
    tleft = BIGMA_sort(nleft,3); % re-adjust left-most spike time
	else
		tright = BIGMA_sort(nright,3);
  end		
  % NOTE: The above code will fail if 'bwidth' is so large that we spill off the spike time grid
  %       in both directions.  
  
	dist_left = tspk - tleft;
  dist_right = tright - tspk;
	dist_avg = (1/2)*(dist_left + dist_right);
	MwidthsB(i,j) = dist_avg;   % set time B cosine bell widths to the average of these two distances
	if dist_avg <= 0    % error check
		fprintf('In trial %d  spike %d  the bell width = %g\n',i,j,dist_avg);
		fprintf('The input time grid must not have contained all the spikes.  Check for errors in time A calculation...\n');
		error('Exiting...');
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
% Integral of rate function between spikes 
	intvals = 0.5*((b-a)./width) + (0.5/pi)*sin(pi*(b-t0)./width) - (0.5/pi)*sin(pi*(a-t0)./width);

%-------------------------------------------------------------------------------------------------------
function [SUPP, lamB_mean, dlamB_mean, tB, numspikes_trial] = ...
		get_ratesB_on_grid(tA,MA,MwidthsB,show_progress)
%-------------------------------------------------------------------------------------------------------
% Get cosine bells for each spike and define the time vector indices of support for each over the chosen
%  time grid.  This is NOT the function that provides the exact calculation at each spike time.  This defines
%  the rate and rate derivative over the grid points only; any spike time likely does not fall on the grid precisely.
% This function also performs the mapping from the time tA  to the new time B  on the grid.
  dtA = mean(diff(tA));
  [~,numspikes_trial] = spikematrix2vec(MA);  % need number of spikes per trial (N) to set cell array size
  numspikes = sum(numspikes_trial);
  numtrials = size(MA,1);       % number of trials (repeats)
  lamB_mean = zeros(1,length(tA));    % mean lamB (averaged across trials)
  dlamB_mean = zeros(1,length(tA));   % mean lamB (averaged across trials)
  SUPP = zeros(numspikes,6);   % matrix of support information [i j is0 is1 tspk widthA]
	                                        % i = trial#, j = spike# in trial
																					% is0 = leftmost support grid index, is1 = rightmost support grid index
																					% tspk = actual spike time, widthA = cosine bell width for spike
  %show_progress = true;
	for i = 1:numtrials
    if show_progress
      if mod(i,10)==0, fprintf('Processing trial %d of %d\n',i,numtrials); end;
    end
    for j = 1:numspikes_trial(i)
      spikenumber = sum(numspikes_trial(1:i-1)) + j;
      tspk = MA(i,j);      % jth spike time (time A) in trial i
      w = MwidthsB(i,j);  % width of associated cosine bell
      tmin = tspk - w; tmax = tspk + w;
      % Check for case where cosine bell terminates within a single grid spacing (better to choose smaller dt!)
      if 2*w < dtA   % note that w really spaces 1/2 the support of the cosine bell
        fprintf('For spike %d of trial %d the cosine bell support is smaller than one grid interval!\n',j,i);
				save MwidthsB.mat tA MwidthsB;
        error('Choose a finer time mesh!');
			else
				% Get time indices of support for cosine bell
        is0 = max(locate(tA,tmin),1);  % Use max(*,1) in case cosine bell drops to zero off to left of grid 
        is1 = locate(tA,tmax);
      end
      if tA(is0) < tmin
        is0 = is0 + 1;  % 'locate' function may return index one too far to the left
      end
      SUPP(spikenumber,:) = [i j is0 is1 tspk w];
      % Error check: make sure that all spikes fall on the grid somewhere (modify code later to handle arbitrary grid subsets)
      if isempty(is0) || isempty(is1)
        fprintf('Bad spike processing for spike %d of trial %d\n',j,i);
        error('Some spikes were out of range of the specified time grid');
      end
      indS = is0:is1;
      lamB_mean(indS) = lamB_mean(indS) + cosbell(tA(indS),tspk,w);
      dlamB_mean(indS) = dlamB_mean(indS) + dcosbell(tA(indS),tspk,w);
    end
  end
  % Now divide lamA_mean by the number of trials to get the true mean (ensemble average, effectively)
  lamB_mean = lamB_mean/numtrials;
  dlamB_mean = dlamB_mean/numtrials;
	% Now map the original times  'tA'  on grid into time B  'tB'  values (just the cumumlative integral of the rate function)
	tB = cumsimps(tA,lamB_mean);  % 4th-order accurate Simpson rule
	
%-------------------------------------------------------------------------------------------------------
function [MlamB, MdlamB, MB] = get_lamB_timeB_each_spike(MA,SUPP,show_progress)
%-------------------------------------------------------------------------------------------------------
  MlamB = zeros(size(MA));
  MdlamB = zeros(size(MA));
	MB = zeros(size(MA));
	[~,numspikes_trial] = spikematrix2vec(MA);  % just need number of spikes per trial 
	numspikes = sum(numspikes_trial);
  numtrials = size(MA,1);       % number of trials (repeats)
	% Sort the SUPP matrix in order of increasing spike time (SUPP keeps track of trial and spike# in trial)
	spikes_all = SUPP(:,5);
	[spikes_all_sort, ispiketime_sort] = sort(spikes_all);
	SUPPSORT = SUPP(ispiketime_sort,:);  % SUPP matrix sorted in order of increasing spike time
  wvals_all_sort = SUPPSORT(:,6);
  tsupp_left = spikes_all_sort - wvals_all_sort;   % tspk - widthA  for each spike
  tsupp_right = spikes_all_sort + wvals_all_sort;  % tspk + widthA for each spike
  
	%show_progress = true;
  for nspike = 1:numspikes
		i = SUPPSORT(nspike,1);   % associated trial# and spike# in trial for the spike
		j = SUPPSORT(nspike,2);
		if show_progress
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
		MlamB(i,j) = sum(c)/numtrials;
		MdlamB(i,j) = sum(dc)/numtrials;
		% Calculate spikes in new time B, keeping spike identity within a trial as we step from one time to the next.
		if nspike == 1
			tspk_previous = 0;
			tB_previous = 0;
      [avals,bvals,wvals,tcentervals] = ...
        get_support_ranges_for_integral(tspk_previous,tspk,tsupp_left,tsupp_right,spikes_all_sort,wvals_all_sort);
			intvalsB = integral_over_interval(avals, bvals, tcentervals, wvals);
			tB = tB_previous + sum(intvalsB)/numtrials;
			MB(i,j) = tB;
		end
		if nspike > 1
      [avals,bvals,wvals,tcentervals] = ...
        get_support_ranges_for_integral(tspk_previous,tspk,tsupp_left,tsupp_right,spikes_all_sort,wvals_all_sort);
			intvalsB = integral_over_interval(avals, bvals, tcentervals, wvals);
			tB = tB_previous + sum(intvalsB)/numtrials;
			MB(i,j) = tB;
		end
% 		if tspk_previous - tspk == 0
% 			fprintf('nspike = %d , (i,j) = (%d,%d)  had duplicate spike time!\n',nspike,i,j);
% 			fprintf('tspk_old = %g , tspk = %g\n',tspk_previous,tspk);
% 			return
% 		end
		tspk_previous = tspk;
		tB_previous = tB;
  end
  
%-------------------------------------------------------------------------------------------------------
	function isupp = get_support_indices_single_spike(tspk,tsupp_left,tsupp_right)
%-------------------------------------------------------------------------------------------------------
% Given a spike time "tspk" return the indices (rows) of SUPP (matrix with support information over the grid, 
% which is not necessary for the exact calculations, incidentally) that correspond to spikes with support 
% including the time "tspk".
	  vsum = sign(tspk - tsupp_left) + sign(tsupp_right - tspk);
    isupp = find(vsum == 2); 		% any entries with vsum = 2 constitude the support indices of tspk

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
