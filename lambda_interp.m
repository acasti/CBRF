function [t,lam,tinterp,laminterp] = lambda_interp(M,Mlam,trange,dt)
%------------------------------------------------------------------------------------------------------
% Interpolates the rate function over a given time range "trange" with sample interval dt.  Uses
%   cubic Hermite polynomial interpolation (pchip).
%
% USAGE:      [tvec,lamvec,tinterp,laminterp] = lambda_interp(M,Mlam,trange,dt);
% INPUT:      M(:,:)                               * (matrix) spike times (with zero padding)
%             Mlam(:,:)                            * (matrix) matrix of rate function vales at spike times
%             trange                               * (vector) 2 element vector; time range of interpolation [t0,t1]
%             dt                                   * (double) sample interval for interpolation
% OUTPUT:     t                                    * (vector) unique spike times from M collapsed onto single axis
%             lam                                  * (vector) rate function at the unique spike times
%             tinterp                              * (vector) time values at which rate function was interpolated
%             laminterp                            * (vector) interpolated rate function at times "tinterp"
%
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 13 September 2015
%------------------------------------------------------------------------------------------------------

%% Argument check
if nargin < 4
	error('Not enough inputs');
end

%% Preliminaries
tinterp = trange(1):dt:trange(2);  % time interval for interpolation

% Collapse the spike times and rate values into a single vector
tspikes = M(:);
lambda = Mlam(:);
% Eliminate zeros (originally for matrix padding purposes) in spike train
inonzero = find(tspikes);
tspikes = tspikes(inonzero);
lambda = lambda(inonzero);  % Corresponding rates for non-zero spike times
[tspikes,isort] = sort(tspikes);    % Sort the spike times
lambda = lambda(isort);             % Correspondingly sort the rate values

% It is sometimes the case that there are identical spike times, which will cause problems in the interpolation scheme.  
% Remove these repeated spike time values for interpolation purposes (valid to do here)
[t,iunique] = unique(tspikes);   % 'unique' automatically sorts the result, but tspikes is already sorted so it's OK
lam = lambda(iunique);            % Get corresponding rate values associated with unique spike times

%% Do cubic Hermite polynomial interpolation on specified grid, using known spike times and rates
fprintf('Calculating spline interpolation...\n'); 
laminterp = interp1(t,lam,tinterp,'pchip','extrap');

