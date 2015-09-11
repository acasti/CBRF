function ISI = get_isi_each_row(M)
%---------------------------------------------------------------------------------------------
% Get interspike intervals for a matrix of spike times (each row is a trial)
% 
% USAGE:           ISI = get_isi_each_row(M);
% INPUT:           M                     * m x n matrix of spike times
% OUTPUT:          ISI                   * m x n matrix of spike intervals (first ISI in each 
%                                          row is set to the spike time; no predecessor)
% 
% Note:  For each trial the ISI associated with the first spike is set to the spike time, which means
%          it is assumed that the beginning of each trial is marked at time 0.
% Written by Alex Casti, MSSM, Nov 2007
% Last modified 18feb2008
%---------------------------------------------------------------------------------------------

if nargin < 1 
  error('Not enough inputs!');
end
if isempty(M)
  ISI = [];   % return empty ISI matrix for empty input
  return
end

ISI = zeros(size(M));
rows = size(M,1);

for i = 1:rows
  cols = length(find(M(i,:)));   
  for j = 1:cols
    if j == 1
      ISI(i,j) = M(i,j);
    else
      ISI(i,j) = M(i,j)-M(i,j-1);
    end
  end
end