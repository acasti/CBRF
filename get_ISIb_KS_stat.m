function errdata = get_ISIb_KS_stat(ISIdata,refindices)
%----------------------------------------------------------------------------------------------------------------------------
% Uses the Kolmogorov-Smirnov (two-sided) test to assess similarity of time B CDF data 
%  for the various data subsets.  The "error data" here simply corresponds to the maximum
%  difference in the 2 CDFs being compared.  Difference measures in the CDF are obtained 
%  only for all time B CDFs relative to a fixed reference CDF (user specified "refindices"),
%  so this function does *not* calculate all possible pairwise differences.
%
% USAGE:    errdata = get_ISIb_KS_stat(ISIdata,refindices)
% INPUT:    ISIdata          * (struct) data struct output of "isi_partition.m"
%           refindices       * (length 2 int vector) reference indices for KS calculation
%                               specifies reference histogram against which to calculate KS stat
% OUTPUT:   errdata          * (struct) output data struct 
%
% Written by Alex Casti, FDU Department of Mathematics
% Last updated 10 September 2015
%----------------------------------------------------------------------------------------------------------------------------

[cols,rows] = size(ISIdata.Bseg);
ISIb = cell(rows,cols);
% Extract the ISI time b interval data
for i = 1:rows
	for j = 1:cols
		ISIb{i,j} = ISIdata.Bseg{j,i}(:,8);
	end
end

% Get the KS statistic for each pair of histograms and keep track of a running average
ksstat = zeros(rows,cols);
h = zeros(rows,cols);
pvalue = zeros(rows,cols);

tail = 'unequal';   % two-sided KS test (tests for maximum differences in 2 CDFs)
alpha = 0.05;       % significance level
for i = 1:rows
	for j = 1:cols
    [h(i,j),pvalue(i,j),ksstat(i,j)] = kstest2(ISIb{refindices(1),refindices(2)},ISIb{i,j},alpha,tail);
	end
end
ksb_avg = sum(ksstat(:))/(numel(ksstat)-1);  % Don't count the obligatory 0 in the average KS value
fprintf('KS average = %g\n',ksb_avg);

% Define output data struct
errdata.ksb_avg = ksb_avg;
errdata.ISIb = ISIb;
errdata.rows = rows;
errdata.cols = cols;
errdata.ksstat = ksstat;
errdata.h = h;
errdata.pvalue = pvalue;
errdata.refindices = refindices;