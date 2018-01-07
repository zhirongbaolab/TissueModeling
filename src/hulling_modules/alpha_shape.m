% ALPHA SHAPE

% ----------------------------------------------
% USAGE
% 
%
%
% ----------------------------------------------


% ----------------------------------------------
% MATLAB DOCS
% Shp = alphaShape(x,y,z,a)
% creates a bounding area or volume 
% that envelops a set of 2-D or 3-D points. 
% You can manipulate the alphaShape object to 
% tighten or loosen the fit around the points to 
% create a nonconvex region. You also can add 
% or remove points or suppress holes or regions.
%
% x, y, z - positions in 3 dimensions
% a - alpha radius, controls the convexity and
%  how tight the shape should fit around the pts
% ----------------------------------------------

function [_] = compute_alpha_shape(pt_cloud)