%% Run many instances of runEverything_cosBells_timeRescaling.m  to get rate functions for many cells and conditions.
dataDirectory = 'C:\Alex_Casti\Experimental_Data\Casti_LGN_Spotential_data_share';

%% Specific cell data 
outDirectory = 'C:\Alex_Casti\Matlab\My_Matlab_Functions\spikes\time_rescaling_data\Data_Experiment\2007_02_27_cell2_XOFF';
CELLNUM = '2';
CELLDATE = '20070227';
CELLTYPE.CELL = 'XOFF';
T = 8;           % repeat trial duration
dt = .001/20;    % time resolution for rate function
checkOverwriteFlag = false;  % check whether data is being overwritten (requires user reply)
bwidths_default = 5:30;

outFile = 'dataAB_RET_cell_02_size80_27feb2007_2.mat';
spikeDataFile = 'sortspikes_spotentials_cell_02_vanHateren_size80.mat';
CELLTYPE.AREA = 'RET';
STIMSIZE = 80;
bwidths = 10:16;
runEverything_cosBells_timeRescaling(bwidths,T,dt,dataDirectory,outFile,outDirectory,spikeDataFile,...
  CELLNUM,CELLDATE,CELLTYPE,STIMSIZE,checkOverwriteFlag);

% outFile = 'dataAB_LGN_cell_02_size80_27feb2007_2.mat';
% spikeDataFile = 'sortspikes_lgn_cell_02_vanHateren_size80.mat';
% CELLTYPE.AREA = 'LGN';
% STIMSIZE = 80;
% bwidths = 6:12;
% runEverything_cosBells_timeRescaling(bwidths,T,dt,dataDirectory,outFile,outDirectory,spikeDataFile,...
%   CELLNUM,CELLDATE,CELLTYPE,STIMSIZE,checkOverwriteFlag);



