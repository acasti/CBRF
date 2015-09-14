function [mse,cv] = get_poisson_err_many_bwidths(MA,tA,bwidths,show_progress)
%----------------------------------------------------------------------------------------------------------------
% Run the time B (cosine bell) transformation with many "bwidth" values
%  (used for the time B cosine bells) and return the mean square error 
%  between the generated spike times (time B) and a unit rate Poisson process.
%
% USAGE:     [mse,cv] = get_poisson_err_many_bwidths(MA,tA,bwidths,show_progress);
% INPUT:     MA                     * matrix of time A spike times
%            tA                     * time A values (on chosen grid; needed to get max tB)
%            bwidths                * vector of "bwidth" values for time B transformation
%            show_progress          * (logical) if true then output progress of computations to Matlab command window
% OUTPUT:    mse                    * vector of mean square errors at each input bwidth
%            cv                     * coefficient of variation for each value of bwidth
% Dependencies:  get_timeB_cosbells.m, get_err_poisson_order_stat.m (and others)
%
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 13 September 2015
%----------------------------------------------------------------------------------------------------------------

%% Argument check
if nargin < 3
  error('Not enough inputs!');
end
if (nargin < 4) || isempty(show_progress)
	show_progress = false;       % (default) do not show progress report of each time B computation as it progresses
end

%% Compute the mse and cv values for all elements of the "bwidths" vector
N = length(bwidths);
mse = zeros(1,N);
cv = zeros(1,N);
for i = 1:N
  fprintf('Processing bwidth = %d  (%d of %d)\n',bwidths(i),i,N);
  dataB = get_timeB_cosbells(MA, tA, bwidths(i),show_progress);
  [mse(i),cv(i),~] = get_err_poisson_order_stat(dataB.MB); 
  fprintf('mse(%d) = %g\n',i,mse(i));
	fprintf('cv(%d) = %g\n',i,cv(i));
  clear dataB;  % clear some variables just for the sake of hygiene
end