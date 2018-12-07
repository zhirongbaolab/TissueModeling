function [] = non_outgrowth_modeling_module(embinfo)

output_path = 'C:\Users\katzmanb\Desktop\TissueModeling\data\output\non_outgrowth_neurons\';
t_ext = '_t';
obj_ext_str = '.obj';
offset = 1;
start_time = 259;

non_outgrowth_config_filename = 'C:\Users\katzmanb\Desktop\TissueModeling\data\configurations\tissue_cells\non_outgrowth_cells.csv';
num_non_outgrowth_cells = 35;
non_outgrowth_config_info = cell(num_non_outgrowth_cells, 2); % 2 for terminal and lineage name

fid = fopen(non_outgrowth_config_filename);
if fid < 0
    error(['could not open file: ' non_outgrowth_config_filename]);
end

i = 1;
while ~feof(fid)
    line = fgetl(fid);
    tokens = regexp(line, ',', 'split');
    
    if size(tokens, 2) == 2
       non_outgrowth_config_info(i, :) = tokens; 
    end
    
    i = i + 1;
end

non_outgrowth_nuc_data = cell(35, 5, 100); % 35 cells, 5 fields of data (x,y,z,diam,name), 100 time points (260-360)

% load the data
for t=260:360
    itr = 1;
    
    cell_data_atT = embinfo(t).celldata;
    cell_names = embinfo(t).cellnames;
    
    for i=1:size(cell_names, 1)
        name = cell_names{i, 1};
        
        idx = find(strcmp(non_outgrowth_config_info(:, 2), name), 1);
        if ~isempty(idx)
            cell_data = cell_data_atT(i, :);
            x = cell_data(4);
            y = cell_data(5);
            z = cell_data(6);
            diam = cell_data(7);
            C = {x,y,z,diam,name};
            
            non_outgrowth_nuc_data(itr, :, (t - 260) + offset) = C;
            itr = itr + 1;
        end
    end
end

% iterate over all time points
for i=1:101
    
    % iterate over all cells at this time point
    for j=1:size(non_outgrowth_nuc_data, 1)
        if isempty(non_outgrowth_nuc_data{j, 1, i})
            break;
        end
        
        % access data
        x = non_outgrowth_nuc_data{j, 1, i};
        y = non_outgrowth_nuc_data{j, 2, i};
        z = non_outgrowth_nuc_data{j, 3, i};
        diam = non_outgrowth_nuc_data{j, 4, i};
        name = non_outgrowth_nuc_data{j, 5, i};
        
        % create a sphere - keep the parameterization low
        [x_coords_sphere, y_coords_sphere, z_coords_sphere] = sphere(5);
        x_coords_sphere_vec = x_coords_sphere(:);
        y_coords_sphere_vec = y_coords_sphere(:);
        z_coords_sphere_vec = z_coords_sphere(:);
        
        % apply transforms: radius expand, translate, scale
        for k=1:size(x_coords_sphere_vec, 1)
            % apply radius of cell
            q = 1; % homogenous coordinates
            xyz = [x_coords_sphere_vec(k)*(diam/2.0); y_coords_sphere_vec(k)*(diam/2.0); z_coords_sphere_vec(k)*(diam/2.0); q];
            
            % create a scale transform and apply
            scale_transform = [1 0 0 0; 0 1.4 0 0; 0 0 1 0; 0 0 0 1];
            xyz = xyz' * scale_transform;
            
            % center at nuc position
            xyz(1) = xyz(1) + x;
            xyz(2) = xyz(2) + y;
            xyz(3) = xyz(3) + z;
            
            x_coords_sphere_vec(k) = xyz(1);
            y_coords_sphere_vec(k) = xyz(2);
            z_coords_sphere_vec(k) = xyz(3);
        end
        
        x_coords_sphere = [x_coords_sphere_vec(1:6)';...
            x_coords_sphere_vec(7:12)';...
            x_coords_sphere_vec(13:18)';...
            x_coords_sphere_vec(19:24)';...
            x_coords_sphere_vec(25:30)';...
            x_coords_sphere_vec(31:36)'];
 
        y_coords_sphere = [y_coords_sphere_vec(1:6)';...
            y_coords_sphere_vec(7:12)';...
            y_coords_sphere_vec(13:18)';...
            y_coords_sphere_vec(19:24)';...
            y_coords_sphere_vec(25:30)';...
            y_coords_sphere_vec(31:36)'];
        
        z_coords_sphere = [z_coords_sphere_vec(1:6)';...
            z_coords_sphere_vec(7:12)';...
            z_coords_sphere_vec(13:18)';...
            z_coords_sphere_vec(19:24)';...
            z_coords_sphere_vec(25:30)';...
            z_coords_sphere_vec(31:36)'];        
        s = mesh(x_coords_sphere, y_coords_sphere, z_coords_sphere);
        
        [faces, vertices, ~] = surf2patch(s, 'triangles');
        
        % flip the faces, surf2patch has the opposite winding/normal
        % paradigm as JavaFX
        faces_flipped = fliplr(faces);
        
        % save the sphere as an obj
            % format the filename
            
        % first find the terminal name for this lineage name
        idx = find(strcmp(non_outgrowth_config_info(:, 2), name), 1);
        if isempty(idx)
            error(['problem with ' name]);
        end
        
        filename = sprintf('%s%s%s%s%s', output_path, name,...
             t_ext, num2str((start_time + i) - offset), obj_ext_str);
    
        saveObjFile(filename,...
            vertices, faces_flipped);
    end
end


end

