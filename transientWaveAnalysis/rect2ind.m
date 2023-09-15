function [x,y] = rect2ind(rect)
%RECT2IND convert rect vector to matrix index vectors
% [x,y] = rect2ind(rect) converts a rect = [left top width height] vector
% to index vectors x and y (column and row indices, respectively), taking
% into account that rect spedifies the location and size with respect to
% the edge of pixels. 
%
% See also IMRECT, IMCROP

left = rect(1);
top = rect(2);
width = rect(3);
height = rect(4);

x = round( left + 0.5 ):round( left + width - 0.5 );
y = round( top + 0.5 ):round( top + height - 0.5 );
