function plotHysteresis(x,y,varargin)
% Very simple function that plots the hysteresis
% It assumes that, if more that one set of data is to be provided, they are
% organised so that the columns of y are in the same number as the columns
% of x

for i = 1:size(y,1)
    quiver(x(1:end-1),y(i,1:end-1),diff(x),diff(y(i,:)),0,varargin{:},'MaxHeadSize',0.03);
end

return