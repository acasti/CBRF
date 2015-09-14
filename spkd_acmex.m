% function d = spkd_acmex(tli,tlj,cost)
%------------------------------------------------------------------------------------------------------------- 
% Calculates the "spike time" distance between two spike trains (Victor & Purpura 1996) for a single cost
%  parameter.  This is a mex-file (C) implementation of the spkd.m Matlab m-file.
%
% USAGE:  d = spkd_acmex(tli,tlj,cost)
% INPUT:  tli: vector of spike times for first spike train 
%         tlj: vector of spike times for second spike train 
%         cost: cost per unit time to move a spike 
% OUTPUT: d: spike time distance at specified cost "cost"
% 
% Copyright (c) 1999 by Daniel Reich and Jonathan Victor. (comments changed by A.Casti) 
% Translated to C from Matlab by Alex Casti, MSSM, 24 jan 2006. 
%-------------------------------------------------------------------------------------------------------------