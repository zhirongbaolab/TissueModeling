% module called with all nuc data


% PART 1 - LOAD TISSUE DATA (tissue specific cells - names, time, diameter, name - at all timepoints)

% nuc data format is 5xNxM
% - 5 items represent x,y,z,diam,name
% - N represents the maximum amount of cells present at a time point for the tissue 
%   (in time points with less than this maximum, blank cells will be indicated by -1s)
% - M represents the number of frames of cells to be modeled  
pharynx_proper_nuc_data = []
arcade_nuc_data = []

% PART 2 - INTERPOLATE TO LANDMARKS (OPTIONAL, IF MODELING ARCADE CELLS)

% PART 3 - TEMPORALLY SMOOTH THE POINT CLOUD IN 3 FRAME INCREMENTS

% PART 4 - EXPAND POINT CLOUD

% PART 5 - COMPUTE AN ALPHA SHAPE ON A POINT CLOUD

% PART 6 - GENERATE AN OBJ FILE


 % *************** AFTER FIRST PASS **********************

% PART 7 - COMPUTE A VORONOI TESSELLATION ON A POINT CLOUD

% PART 8 - GENERATE AN OBJ FILE

% PART 9 - INTERSECT THE ALPHA SHAPE WITH THE VORONOI TESSELLATION
