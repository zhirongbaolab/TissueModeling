function [] = pharynx_modeling_module(embinfo, config_info)


% module called with all nuc data


% PART 1 - LOAD TISSUE DATA (tissue specific cells - names, time, diameter, name - at all timepoints)
% pharynx_nuc_data format is 5xNxM
% - 5 items represent x,y,z,diam,name
% - N represents the maximum amount of cells present at a time point for the tissue - this is 83
%       (in time points with less than this maximum, blank cells will be indicated by empty cells)
% - M represents the number of frames of cells to be modeled  
pharynx_nuc_data = cell(83, 5, (config_info{1, 6} - config_info{1, 5})); % row, col, frame (i.e. cell, data fields, time)

% load the pharynx lineage names from the file
pharynx_cell_names = {};
pharynx_cell_names_file_path = config_info{1, 4}{:};
fid = fopen(pharynx_cell_names_file_path);
if fid < 0
    error(['could not open file: ' pharynx_cell_names_file_path]);
end

while ~feof(fid)
    line = fgetl(fid);
    lineage_name = line(1:strfind(line, ' (')-1);
    pharynx_cell_names{end+1} = lineage_name; % consider preallocating this size and then shrinking it after loading
end

% iterate over the embryo data
% TODO --> modularize this for reuse later
for t=config_info{1, 5}:config_info{1, 6}
    it = 1;
    % access the nuclei info at the time point
    cell_data_atT = embinfo(t).celldata;
    cell_names = embinfo(t).cellnames;

    % make sure the arrays are parallel in the number of rows
    if size(cell_data_atT, 1) == size(cell_names, 1)
        % iterate over all the cell names
        for i=1:size(cell_names, 1)
            name = cell_names{i, 1};
            
            % NEEDS TO BE ENCODED AS EXPLICIT RULE IN CONFIG FILE
            if t >= 336 && (strcmp(name, 'ABaraaapaaa') || strcmp(name, 'ABalpaapaaa'))
               continue; 
            end
           
            % check if the lineage name is part of the pharynx by checking the
            % pharynx names file
            idx = find(strcmp(pharynx_cell_names, name), 1);
            if ~isempty(idx)
                % fprintf('%s, %d\n',name, idx);
                % use the index of this cell name in the embryo info 
                % to access cell info in parallel celldata list
                cell_data = cell_data_atT(i, :);
                x = cell_data(4);
                y = cell_data(5);
                z = cell_data(6);
                diam = cell_data(7);
                C = {x,y,z,diam,name};
                
                % add the coordinates, diameter, and name to the
                % pharynx_nuc_data matrix
                pharynx_nuc_data(it, :, (t - config_info{1, 5} + 1)) = C; % entry row, data fields, time frame
                it = it + 1;
            end
        end
    end
end

% PART 2 - INTERPOLATE TO LANDMARKS (OPTIONAL, IF MODELING ARCADE CELLS)

% PART 3 - TEMPORALLY SMOOTH THE POINT CLOUDS IN 3 FRAME INCREMENTS
temporally_smoothed_pc = conv_temp_smoothing(pharynx_nuc_data, 2);
% visualize_pt_cloud(temporally_smoothed_pc, 'Pharynx nuc temporally smoothed');

% PART 4 - EXPAND POINT CLOUDS
temporally_smoothed_spherically_expanded_pc = spherical_expansion(temporally_smoothed_pc);

% PART 5 - COMPUTE ALPHA SHAPES ON POINT CLOUDS
shapes = alpha_shape(temporally_smoothed_spherically_expanded_pc);

% PART 6 - GENERATE OBJ FILES
output_path = 'C:\Users\katzmanb\Desktop\TissueModeling\data\output\pharynx\';
pharynx_str = 'Pharynx_t';
obj_ext_str = '.obj';
offset = 1;
emb_time_it = config_info{1, 5};
for i=1:size(shapes, 1)
    % format the filename
    filename = sprintf('%s%s%s%s', output_path,...
        pharynx_str, num2str(emb_time_it - 1), obj_ext_str);
    
    saveObjFile(filename,...
    shapes{i, 1}{2}, shapes{i, 1}{1});
    
    emb_time_it = emb_time_it + 1;
end
   % config_info{1, 6}

 % *************** AFTER FIRST PASS **********************

% PART 7 - COMPUTE A VORONOI TESSELLATION ON A POINT CLOUD

% PART 8 - GENERATE AN OBJ FILE

% PART 9 - INTERSECT THE ALPHA SHAPE WITH THE VORONOI TESSELLATION

end
