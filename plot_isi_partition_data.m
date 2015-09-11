function plot_isi_partition_data(ISIdata,xmax,ymax,suppressIJ)
%-----------------------------------------
%         ISI HISTOGRAM PLOTS
%-----------------------------------------
B = ISIdata.B;
ISIb = ISIdata.ISIb;
Bseg = ISIdata.Bseg;
numsegx = ISIdata.numsegx;
numsegy = ISIdata.numsegy;
meanrate_seg = ISIdata.meanrate_seg;

%numbins = ceil(10*log10(mean(numspikes_seg(:))));
%numbins = 32  % fix the number of bins for purposes of commonality across cells, stim, etc
% Define a common ISI axis (new time)
minISI = 0;
maxISI = max(B(:,8));
binwidth = 0.275072878909593;  % fix the bin width
%edges = linspace(minISI,maxISI,numbins+1); 
edges = minISI:binwidth:maxISI;
dISI = mean(diff(edges));
bincenters = edges(1:end-1)+0.5*dISI;
Nhist = cell(numsegx,numsegy);
Nhist_normalized = cell(numsegx,numsegy);
maxprob = 0;  % Calculate maximum probability in any window for axis scaling

%% Calculate grand ISI histogram (time B) for the entire data set
Nhist_all_temp = histc(ISIb,edges); Nhist_all = Nhist_all_temp(1:end-1);
Nhist_all_normalized = (1/dISI)*Nhist_all/sum(Nhist_all);   %#ok<NASGU>

%% Calculate the histogram data
for i = 1:numsegx
  for j = 1:numsegy
    Nhist_temp = histc(Bseg{i,j}(:,8),edges); Nhist{i,j} = Nhist_temp(1:end-1);
    %Nhist_normalized{i,j} = Nhist{i,j}/sum(Nhist{i,j});  % Sum of prob is unity (alternate normalization)
    % Normalize each histogram to have unit area 
    Nhist_normalized{i,j} = (1/dISI)*Nhist{i,j}/sum(Nhist{i,j});  % Integral of prob is unity
    maxprob = max(maxprob,max(Nhist_normalized{i,j}));
  end
end

%% Plot the time B ISI histogram data
newfigure = true;
xmin = 0;         % Minimum value on x-axis (ISI) for axis scaling
ymin = 0;         % Minimum value on y-axis (ISI density) for axis scaling
if (nargin < 2) || isempty(xmax)
  xmax = maxISI;    % Maximum value on x-axis (ISI) for axis scaling
end
if (nargin < 3) || isempty(ymax)
  ymax = maxprob;
end
if nargin < 4 || isempty(suppressIJ)
  suppressIJ = false;
end
if newfigure, figure; end;
for i = 1:numsegx
  for j = 1:numsegy
     % Order plot so that lower left corner is lowest left box in lam vs dlam plane, etc
      k = i + 2*numsegx - numsegx*(j-1);  % This formula does the trick
      subplot(numsegy,numsegx,k)
      %plot(edges_plot,Nhist_normalized{i,j});
      stairs(bincenters,Nhist_normalized{i,j},'b','LineWidth',2); hold('on');
      %stairs(bincenters,Nhist_all_normalized,'r','LineWidth',1);
      if ~suppressIJ
        title(sprintf('(%d,%d)',i,j),'fontweight','bold');
      end
      axis([xmin xmax ymin ymax]);
      ratestr = sprintf('r = %4.1f',meanrate_seg(i,j));
      fsz = 10;
      xrate = .3*xmax;
      text(xrate,.9*ymax,ratestr,'fontweight','bold','fontsize',fsz);      
      if (i==1) && (j==1)
        xlabel('ISI','fontweight','bold');
      else
        set(gca,'Xticklabel','');
      end
      %if mod(k-1,numsegx)==0
      if (i==1) && (j==1)
        ylabel('ISI density','fontweight','bold');
      else 
        set(gca,'Yticklabel','');
      end
 end
end

