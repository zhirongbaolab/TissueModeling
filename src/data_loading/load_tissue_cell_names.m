% load_tissue_cell_names.m

function [tissue_cell_names] = load_tissue_cell_names(filename)

tissue_cell_names = []

% open the csv file, make sure that it exists
%READ FILE
fid = fopen(filename);

if fid < 0
     error(['could not open file: ' filename]);  
end


% read the file data into a text array with the given format
C = textscan(fid, '%s%s%d%d', 'delimiter', ',');

% close the file
fclose(fid);

% skip over the header load_tissue_cell_names

% check if line is valid

% add name to list


end