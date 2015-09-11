function [R,varargout] = get_pseudo_randnums_for_spiketimes(seedprime,decplaces)
%---------------------------------------------------------------------------------------------------------
% Generate pseudo random numbers (PRNs) to add to each spike time (at a certain decimal place of precision).
%
% USAGE:       [R,varargout] = get_pseudo_randnums_for_spiketimes(seedprime,scalefac)
% INPUT:        seedprime                       * A prime number with N digits
%               decplaces                       * (int) number of decimal places to the right that PRNs will be added
%                                                      (e.g.  4  will give numbers at 1/100th msec decimal
%                                                       poistion; default value of this argument is 4)
% OUTPUT:     R                                 * (int vec) non-repeating integer vector of length 10^N - 1
%             varargout{1}                      * (double vec) scaled version of R suitable for adding to spike times (sec)
%
% Method:  Start with a prime number and add it to itself repeatedly.  Keep just the rightmost last N digits
%              (where N is the number of digits of the seed prime) of each result to get a unique integer sequence
%              of pseudo-random numbers (PRNs) that has a period of 10^N.
%
% Comments:
%    (1) Default 'decplaces' will be chosen automatically such that R values (scaled) can be added to spike times
%         (assumed in seconds), and so that the data whitening takes place at a precision 1/100th msec (i.e. a
%         spike time like 3.3419 sec  would have N pseudo-random integers placed after the '9').
% 
% Written by Alex Casti, MSSM January 2008
% Last updated 07 Jan 2008
%---------------------------------------------------------------------------------------------------------

% Argument checking
if ~isprime(seedprime)
	error('Must input a prime number for seeding!');
end
if (nargin < 2) || isempty(decplaces)
	decplaces = 4;   % Default to 4 decimal places (meaning that random numbers will be scaled to 1/100 msec level)
end

num_extraout = nargout - 1;   % number of extra output arguments specified
numdigits = length(num2str(seedprime));   % number of digits in prime number
nexp = numdigits + decplaces;    
scalefac = 10^-nexp;   % Multiply the PRNs by 'scalefac' to get numbers that can directly be added
										   %  to spike times for data whitening chosen or default 'decplaces' position

R = zeros(10^numdigits - 1,1);    % Number of pseudo-random numbers to generate
z = seedprime*(1:length(R))';      % Ignore the final component that will result in all zeros
																						 % Note: it's important that z is a column vector here
zst = int2str(z);          % Turns integer column vector into a column vector of strings
zst_chop = zst(:,end-numdigits+1:end);
R = str2num(zst_chop);  %#ok<ST2NM>

if num_extraout == 1
	varargout{1} = R*scalefac;
end

