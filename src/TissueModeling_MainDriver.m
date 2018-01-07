% TissueModeling_MainDriver.m
% Zhirong Bao Lab, Sloan-Kettering Institute
% Author: Braden Katzman
% Created On: January 4, 2017
%
% This is the main driver for the tissue modeling pipeline for the pharynx, muscle and hypoderm
% in the C. elegans embryo.
%



% PART 1 - LOAD THE NUC DATA

% first, access the configurations field for the nuc data:
%     first time
%     last time
%     path to nuc files

% open DESIGN QUESTION - how to represent the tissue specific nuc to isolate this from the whole embryo of nucs
% - .txt file with names (this can make it more generic)
pharyx_nuc = []
% muscle_nuc = []
% hypoderm_nuc = []


% PART 2 - INTERPOLATE TO LANDMARKS, EXPAND POINT CLOUD (OPTIONAL)

% PART 3 - COMPUTE AN ALPHA SHAPE

% PART 4 - COMPUTE A VORONOI TESSELLATION (not on first iteration of development)

% PART 5 - INTERSECT THE ALPHA SHAPE WITH THE VORONOI TESSELLATION

