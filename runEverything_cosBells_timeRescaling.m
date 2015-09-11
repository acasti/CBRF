function runEverything_cosBells_timeRescaling(bwidths,T,dt,dataDirectory,outFile,outDirectory,spikeDataFile,...
  CELLNUM,CELLDATE,CELLTYPE,STIMSIZE,checkOverwriteFlag)
%%----------------------------------------------------------------------------------------------------------------------------
% Get the rate function using the cosine bells and time-rescaling theorem procedure.  This function streamlines
%   all the basic data preparation steps and function calls.
%
% USAGE:   runEverything_cosBells_timeRescaling(bwidths,T,dt,dataDirectory,outFile,outDirectory,spikeDataFile,...
%              CELLNUM,CELLDATE,CELLTYPE,STIMSIZE,checkOverwriteFlag)
% INPUT:   bwidths                  * (int vec) integer vector of "time B" parameter values to cycle through
%          T                        * (double) repeat trial time domain assumed to be [0,T]
%          dt                       * (double) time grid precision (original time)
%          dataDirectory            * (string) path to data directory (root level only required)
%          outFile                  * (string) name of output file, i.e. "dataAB_LGN_cell_date.mat"
%          outDirectory             * (string) output directory path
%          spikeDataFile            * (string) specific name of data file with spike times 
%          CELLNUM                  * (string) cell number 
%          CELLDATE                 * (string) date of recording
%          CELLTYPE.AREA            * (string struct) recording area (e.g. 'LGN')
%                  .CELL                              cell type (e.g. 'XOFF')
%          STIMSIZE                 * (double) stimulus size (i.e. spot size), if relevant
%          checkOverwriteFlag       * (logical) if true then check to make sure you're not overwriting existing data
%
%  Written by Alex Casti, FDU Mathematics
%  Last updated 08/31/2015
%----------------------------------------------------------------------------------------------------------------------------

%% Specify paths to the data and other information for output data file
addpath(genpath(dataDirectory));
load(spikeDataFile,'Srepeat2');    % loads cell array of spike times
outFileDirectory = strcat(outDirectory,'\',outFile);
if checkOverwriteFlag  % be careful not to overwrite existing data file
  if exist(outFileDirectory,'file')   
    overWriteFlag = input('Output data file exists!  Overwrite? (0) No : ');
    if ~overWriteFlag
      fprint('Stopping program.  Choose a different output file name.');
    end
  end
end
fprintf('%s  cell %s  size %g  (%s  %s)\n',CELLDATE,CELLNUM,STIMSIZE,CELLTYPE.AREA,CELLTYPE.CELL);

%% Specify some parameters required of the cosine bell rate function algorithm
trange = [0,T];     % time interval over which the rate function is constructed
randflag = false;   % this parameter indicates whether a small random number should be added to each spike time
M = spiketime_cell2mat(Srepeat2);   % matrix version of spike times in each trial (zero padded)

%% Time A transformation
dataA = get_timeA_cosbells(M, trange, dt, randflag);
% optionsA.randflag = false;
% optionsA.show_progress = true;
% optionsA.method = 'interp';
% dataA = get_timeA_cosbells2(M, trange, dt, optionsA);
%% Time B transformation, determine b value that gives a rate function that best satisfies the time rescaling theorem
%bwidths = 25:50;   % in the present code each 'b' value must be an integer
show_progress = true;
[mse,cv] = get_poisson_err_many_bwidths(dataA.MA,dataA.tA,bwidths,show_progress);

%% Define the optimal 'b' value to be that which minimizes 'mse'
[msemin,imin] = min(mse); 
bwidth = bwidths(imin);        % optimal b value (usign time rescaling theorem)
plotyy_mse_cv(bwidths,mse,cv);
title(sprintf('%s  cell %s  size %g  (%s  %s), b = %g\n',CELLDATE,CELLNUM,STIMSIZE,CELLTYPE.AREA,CELLTYPE.CELL,bwidth));
fprintf('\nPoisson order statistic error minimized for b = %d  (msemin = %g)\n',bwidth,msemin);
fprintf('CV = %g  for b = %d\n',cv(imin),bwidth);

%% Get time B data associated with the optimal parameter b
fprintf('\nGetting time B data for b = %d\n',bwidth);
dataB = get_timeB_cosbells(dataA.MA,dataA.tA, bwidth);


%% Get the rate function r(t) and dr/dt in terms of the original lab time t
lamvec = dataA.lamA_mean .* dataB.lamB_mean;      %#ok<NASGU> % Rate function (see documentation)
dlamvec = dataA.dlamA_mean .* dataB.dlamB_mean;   %#ok<NASGU> % Time derivative of rate function
Mlam = dataA.MlamA .* dataB.MlamB;             % Matrix of rate values for each spike
Mdlam = dataA.MdlamA .* dataB.MdlamB;          % Matrix of rate derivative values for each spike
% Note: The rows of the matrices Mrate and Mdrate are trials, and each column is an associated
%       event within the trial.  The matrices are zero padded to make them rectangular (each trial
%       can have a different number of spikes).

%% Some other useful information to put into data file
[MeanRateData,~,~] = get_mean_rate_CV_per_trial(M,T);
rate = MeanRateData.avg;      %#ok<NASGU> % mean firing rate across all trials
showPlotISIdata = true;
ISIdata = isi_partition2(M,dataB.MB,Mlam,Mdlam,showPlotISIdata);
refindices = [2 4];     % reference ISIb (time B) ISI histogram data for pairwise comparison
errdata = get_ISIb_KS_stat(ISIdata,refindices);  %#ok<NASGU> % pairwise ISIb histogram "error" (difference) data, KS stat

%% Save the data file for later use
save(outFileDirectory,'bwidth','bwidths','CELLNUM','cv','dataA','dataB','CELLDATE','dlamvec','dt','errdata','ISIdata','lamvec',...
  'M','Mlam','Mdlam','mse','msemin','cv','STIMSIZE','Srepeat2','T','CELLTYPE');