function h = rasterplot(data,win,rgb,spikeheight,newfig)
%----------------------------------------------------------------------------------------------------------------------------
% Make a raster plot.  This code will show blank lines in their proper order for trials in which
%   there are no events.
%
% USAGE:     h = rasterplot(data,win,rgb,spikeheight,newfig);
% INPUT:     data              * (cell array) spiketimes (sec) for raster plot
%            win               * (length-2 vector) for plot windowing (time axis)
%            rgb               * (length-3 RGB vector) raster dot color (default BLACK)
%            spikeheight       * (double) line size for plotting (default .75)
%            newfig            * (integer) flag to create new figure window (default FALSE)
% OUTPUT:    h                 * plot handle
%----------------------------------------------------------------------------------------------------------------------------

if nargin < 2
  error('Not enough inputs!');
end
if nargin < 3 || isempty(rgb)
  rgb = [0,0,0]; % Default black
end
if nargin < 4 || isempty(spikeheight)
  spikeheight = .75;  % Default "dot" size
end
if nargin < 5 || isempty(newfig)
  newfig = 0;   % Do not open new figure by default
end

pre = win(1);
post = win(2);

j=length(data);
w=length(data{1});

for i = 2:j
  w = max(w,length(data{i}));
end

sr = j;
sc = w;

y = (1:sr)';
Y = repmat(y,1,sc);

spikes = zeros(sr,sc);
for i =1:j
  spikes(i,1:length(data{i}))=data{i};
end

N = sr*sc;
SpikeVec(1:N) = spikes(1:N);
TrialVec(1:N) = Y(1:N);
blanks = find(isnan(SpikeVec));
SpikeVec(blanks) = [];
TrialVec(blanks) = [];

s = SpikeVec;
y0 = TrialVec;
y1 = TrialVec-spikeheight;
if newfig, figure;end;
h = line([s;s],[y0;y1],'color',rgb);
axis([pre post 0 sr]);
% axis off;
hold on;
line([0 0],[0 sr],'color','k');
hold off;

xlabel('time (sec)');
ylabel('trial');
