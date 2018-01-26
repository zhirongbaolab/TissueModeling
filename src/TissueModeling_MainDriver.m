% TissueModeling_MainDriver.m
% Zhirong Bao Lab, Sloan-Kettering Institute
% Author: Braden Katzman
% Created On: January 4, 2017
%
% This is the main driver for the tissue modeling pipeline for the pharynx, muscle and hypoderm
% in the C. elegans embryo.
%

% PART 0 - LOAD THE CONFIG INFO
load_dataset_config.m

% ------------------------ THIS IS ONLY NECESSARY ONCE PER EMRBYO --------------------
% ----- make this into separate tissue_data.m object ----------


% PART 1 - LOAD THE NUC DATA
loadEmbryo.m


% PART 2 - call the respective modeling modules
pharynx_modeling_driver.m


% PART 3 - LOAD TISSUE DATA (tissue specific cells - names, time, diameter, name - at all timepoints)
% nuc data format is 5xNxM
% - 5 items represent x,y,z,diam,name
% - N represents the maximum amount of cells present at a time point for the tissue 
%   (in time points with less than this maximum, blank cells will be indicated by -1s)
% - M represents the number of frames of cells to be modeled  

%/
% dl_muscle_nuc = []
% dr_muscle_nuc = []
% vl_muscle_nuc = []
% vr_muscle_nuc = []
% hypoderm_nuc = []
/%

% PART 4 - SAVE THIS INFO TO FILE


