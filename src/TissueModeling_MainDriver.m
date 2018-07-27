% TissueModeling_MainDriver.m
% Zhirong Bao Lab, Sloan-Kettering Institute
% Author: Braden Katzman
% Created On: January 4, 2018
%
% This is the main driver for the tissue modeling pipeline for the pharynx, muscle and hypoderm
% in the C. elegans embryo.
%

% PART 0 - LOAD THE CONFIG INFO
% config_info: (nx6) containing dataset name, object name, nuclei resource
% location, tissue resource location, start time, end time
num_tissues = 7;
num_fields = 6;
config_info_file_path = 'C:\Users\katzmanb\Desktop\TissueModeling/data/configurations/Tissues_Config.csv';
config_info_file_header_str = 'Dataset Name,Object Name,Nuclei Resource Location,Tissue Resource Location,Start Time,End Time';
config_info = cell(num_tissues, num_fields);
fid = fopen(config_info_file_path);

if fid < 0
    error(['could not open file: ' config_info_file_path]);
end

% iterate until end of file
it = 1;
while ~feof(fid)
    % retrieve the line
    line = fgetl(fid);
    if strcmp(line, config_info_file_header_str)
        line = fgetl(fid);
    end
    
    % format the line into a cell array
    C = textscan(line, '%s%s%s%s%d%d', 'Delimiter', ',');
    
    % add the cell array to the configuration info cell matrix
    config_info(it, :) = C;
    
    % increase the iterator
    it=it+1;
end

%PART 1 - LOAD THE NUC DATA
t = 't';
nuc_str = '-nuclei';
pharynx_first_file = strcat(config_info{1, 3}, t, num2str(config_info{1, 5}), nuc_str);
[embinfo, errors] = loadEmbryo_unzipped(config_info{1, 3}, config_info{8, 5},  config_info{8, 6});

% PART 2 - call the respective modeling modules
% TODO make switch case that looks at a command line arg and determines which thing to model 
%pharynx_modeling_module(embinfo, config_info);
%hypoderm_modeling_module(embinfo);
embryo_modeling_module(embinfo);
% non_outgrowth_modeling_module(embinfo);

% PART 3 - LOAD TISSUE DATA (tissue specific cells - names, time, diameter, name - at all timepoints)
% nuc data format is 5xNxM
% - 5 items represent x,y,z,diam,name
% - N represents the maximum amount of cells present at a time point for the tissue 
%   (in time points with less than this maximum, blank cells will be indicated by -1s)
% - M represents the number of frames of cells to be modeled  


