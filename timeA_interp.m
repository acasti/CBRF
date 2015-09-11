function MAI = timeA_interp(M,t,tA)
%--------------------------------------------------------------------------------------------------------
% Given spike times in M, a grid of time values t, and the associated "time A" time values tA on the grid,
%  map the original spike times into time A.
%
%  USAGE:    MAI = timeA_interp(M,t,tA)
%  INPUT:     M
%                t
%                tA
% 
% Written by Alex Casti, MSSM 12 Jan 2008
%--------------------------------------------------------------------------------------------------------
numtrials = size(M,1);
MAI = zeros(size(M));
[v,numspikes_trial] = spikematrix2vec(M);  % just need number of spikes per trial 

Nt = length(t);
for i = 1:numtrials
	for j = 1:numspikes_trial(i)
		tspk = M(i,j);
		ispk = locate(t,tspk);   % locate spike on t grid
		if ispk ~= Nt
			%MAI(i,j) = tA(ispk);
   		MAI(i,j) = 0.5*(tA(ispk)+tA(ispk+1));
		else
		  MAI(i,j) = tA(ispk);
		end
	end
end
