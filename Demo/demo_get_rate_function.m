%% demo_get_rate_function.m
%------------------------------------------------------------------------------------------------------------------------
%  Demo that constructs the rate function from repeat trials of laboratory spike times
%  Data comes from retinal ganglion cell (X-OFF) responses driven by a noisy flashing spot for
%   128 repeated trials (frozen noise stimulus), each of 8 second duration.  The retinal event
%   times come from the recording of "S-potentials" in the cat lateral geniculate nucleus (LGN).
%
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 13 September 2015
%------------------------------------------------------------------------------------------------------------------------

path(path,'..');  % add path to data and files one directory up

%% Load the demo data (spike times are in seconds)
load data_example.mat;

%% Construct the time A data
dt = .001/20;      % Define rate function over a grid spacing of 1/20 msec
T = 8;             % Trial duration (8 seconds in each repeat trial for the demo data)
trange = [0,T];    % Time interval over which to construct the rate function
randflag = false;  % If 'true' then random numbers are added to the 4th decimal place
                   %   to deal with spike times that may be coincident across trials, which
                   %   can cause problems if the spike times were not finely sampled. (not needed for demo data)
dataA = get_timeA_cosbells(Mspikes, trange, dt, randflag);

%% Time B transformation 
% Choose a range of 'b' parameters for the time B transformation and choose the value that
% makes the time B event times closest to unit-rate homogeneous Poisson. For this example
% data set the optimal value is  b = 18.  Generally, one will need to test a range of values
% by trial and error, although one could write a code that brackets a minimum in the mean square
% error criterion.

% This will take a few minutes to run on a modern desktop.  There are a total of 23,133 spikes
%  in the example data set.
%-----------------------------------------------------------------------------------------------
% NOTE: future releases will include a proper optimization scheme over integer b values
%-----------------------------------------------------------------------------------------------

bvalues = 14:22;       % in the present code each 'b' value must be an integer
show_progress = true;  % display progress to the Matlab command window
[mse,cv] = get_poisson_err_many_bwidths(dataA.MA,dataA.tA,bvalues,show_progress);
% The variable 'mse' is the mean square departure from the expected Poisson order statistics
%  for each value of b.  The variable 'cv' is the associated coefficient of variation for each b.
%  Choose the b value that minimizes the 'mse' vector.  In general, this value will correspond
%  very closely to cv = 1, though rarely exactly.

%% Plot 'mse' versus 'b' to see if we have found a minimum
plotyy_mse_cv(bvalues,mse,cv);
% Define the optimal 'b' value to be that which minimizes 'mse'
[msemin,imin] = min(mse); 
b = bvalues(imin); 
title(sprintf('Optimal "time B" parameter:  b = %d  (CV = %g)',b,cv(imin)),'fontweight','bold','fontsize',10);
fprintf('\nPoisson order statistic error minimized for b = %d\n',b);
fprintf('CV = %g  for b = %d\n',cv(imin),b);

%% Get time B data associated with the optimal parameter b
fprintf('\nGetting time B data for b = %d\n',b);
dataB = get_timeB_cosbells(dataA.MA,dataA.tA, b);


%% Get the rate function r(t) and dr/dt in terms of the original lab time t
rate = dataA.lamA_mean .* dataB.lamB_mean;      % Rate function (see documentation)
drate = dataA.dlamA_mean .* dataB.dlamB_mean;   % Time derivative of rate function
Mrate = dataA.MlamA .* dataB.MlamB;             % Matrix of rate values at each spike time
Mdrate = dataA.MdlamA .* dataB.MdlamB;          % Matrix of rate derivative values for each spike time
% Note: The rows of the matrices Mrate and Mdrate are trials, and each column is an associated
%       event within the trial.  The matrices are zero padded to make them rectangular (each trial
%       can have a different number of spikes).