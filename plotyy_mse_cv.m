function plotyy_mse_cv(bwidths,mse,cv)
%------------------------------------------------------------------------------------
% Plot mse (mean square error of Poisson statistic for all time B values) and cv 
%  (coefficient of variation) for all values of the time B parameter "bwidths".  
%
% USAGE:    plotyy_mse_cv(bwidths,mse,cv);
% INPUT:    bwidths         * vector of time B "bandwidth" values
%           mse             * vector of mean square errors in Poisson test for time B
%           cv              * vector of coefficient of variation values (time B)
%------------------------------------------------------------------------------------

%% Argument check
if nargin < 3
  error('Insufficient input arguments!');
end

%% Make the plot
figure
[AX,H1,H2] = plotyy(bwidths,mse,bwidths,cv,@semilogx);
hold(AX(1), 'on');
hold(AX(2), 'on');
% Plot a line at CV = 1
semilogx(bwidths, ones(1,length(bwidths+2)), 'k', 'Parent', AX(2)); 
%plot(AX(2),[min(bwidths) max(bwidths)],ones(1,2),'k-','LineWidth',1.5);
set(H1,'Marker','*');
set(H2,'Marker','*');
xlabel('b   (cosine bell time B parameter)','fontweight','bold');
set(get(AX(1),'Ylabel'),'String','MSE (Poisson arrival error)','fontweight','bold');
set(get(AX(2),'Ylabel'),'String','CV','fontweight','bold');
