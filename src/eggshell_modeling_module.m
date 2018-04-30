function [] = eggshell_modeling_module(embinfo)

output_path = 'C:\Users\katzmanb\Desktop\TissueModeling/data/output/eggshell/';
eggshell_str = 'Eggshell_t';
obj_ext_str = '.obj';
offset = 1;

% read the eggshell configuation file: 
% start time, end time, resource location, dimension y (rows), dimensions x (columns), comments
eggshell_config_info_filename = 'C:\Users\katzmanb\Desktop\TissueModeling/data/configurations/tissue_cells/eggshell/eggshell_config.csv';
eggshell_config_info = cell (15, 6);
fid = fopen(eggshell_config_info_filename);
if fid < 0
    error(['could not open file: ' eggshell_config_info_filename]);
end

% get the first line out of the way
line = fgetl(fid);
i = 1;
while ~feof(fid)
    line = fgetl(fid);
    tokens = regexp(line, ',', 'split');
    
    if size(tokens, 2) == 6
       eggshell_config_info(i, :) = tokens;
    end
   i = i + 1;    
end

% iterate over the time window segments
for i = 1:2
    % read the names matrix [nxmx1] for the current window
    % rows, columns, names
    names_mat = cell(str2num(eggshell_config_info{i, 4}), str2num(eggshell_config_info{i, 5}));
    names_list = {}; % for easier access later
    fid = fopen(eggshell_config_info{i, 3});
    if fid < 0
        error(['could not open file: ' eggshell_config_info{i, 3}]);
    end
    
    k = 1;
    j = 1;
    while ~feof(fid)
        line = fgetl(fid);
        
        tokens = regexp(line, ',', 'split');
        
        if size(tokens, 2) == str2num(eggshell_config_info{i, 5})
          if k <= str2num(eggshell_config_info{i, 4})
             names_mat(k, :) = tokens; 
             k = k + 1;
             
             for token = tokens
                % need to check if there is a semicolon in the token -->
                % indicates multiple names
                lineage_name = token{:};
                if isempty(strfind(lineage_name, ';'))
                   % just one name, extract lineage name
                   if ~isempty(strfind(lineage_name, '('))
                       lineage_name = lineage_name(1:strfind(lineage_name, ' (')-1);
                   end
                       
                   % add it to list
                   if ~isempty(lineage_name)
                       names_list{j} = lineage_name;
                       j = j + 1;           
                   end
                else
                    % multiple names, separate them
                    mult_names = regexp(token, ';', 'split');
                    
                    % extract lineage names
                    for ln = mult_names{1,1}
                        lineage_name = ln{:};
                        
                        if ~isempty(strfind(lineage_name, '('))
                            lineage_name = lineage_name(1:strfind(lineage_name, ' (')-1);
                        end
                            
                        % make sure it's not a number
                        [~, status] = str2num(lineage_name);
                        
                        % add it to list
                        if ~isempty(lineage_name) && status == 0
                            names_list{j} = lineage_name;
                            j = j + 1;
                        end
                   end
                end
             end
          end
        end
    end
    
    % use the names to access the model data over the entire window
    eggshell_nuc_data = cell((str2num(eggshell_config_info{i, 4}) * str2num(eggshell_config_info{i, 5})),...
        5,...
        (str2num(eggshell_config_info{i, 2}) - str2num(eggshell_config_info{i, 1}))); 
    % row, col, frame (i.e. cells - rowsxcols in names_mat, data fields - 5, time - end_time-start_time)
    
    % iterate over model data to fill out eggshell_nuc_data
    for t = str2num(eggshell_config_info{i, 1}):str2num(eggshell_config_info{i, 2})
        it = 1;
        % access the nuclei info at the time point
        cell_data_atT = embinfo(t).celldata;
        cell_names = embinfo(t).cellnames;

        % make sure the arrays are parallel in the number of rows
        if size(cell_data_atT, 1) == size(cell_names, 1)
            % iterate over all the cell names
            for g=1:size(cell_names, 1)
                name = cell_names{g, 1};

                % check if the lineage name is part of the eggshell in this
                % configuration
                idx = find(strcmp(names_list, name), 1);
                if ~isempty(idx)
                    % fprintf('%s, %d\n',name, idx);
                    % use the index of this cell name in the embryo info 
                    % to access cell info in parallel celldata list
                    cell_data = cell_data_atT(g, :);
                    x = cell_data(4);
                    y = cell_data(5);
                    z = cell_data(6);
                    diam = cell_data(7);
                    C = {x,y,z,diam,name};

                    % add the coordinates, diameter, and name to the
                    % pharynx_nuc_data matrix
                    eggshell_nuc_data(it, :, (t - str2num(eggshell_config_info{i, 1}) + 1)) = C; % entry row, data fields, time frame
                    it = it + 1;
                end
            end
        end
    end
    
    
    % temporally smooth the model data in this window
    temporally_smoothed_eggshell_nuc_data = conv_temp_smoothing(eggshell_nuc_data, 2);
    
    % for each frame in the current window
    for t = str2num(eggshell_config_info{i, 1}):str2num(eggshell_config_info{i, 2})
        % parse for syntax triggering computation:
        % 1. Name --> use position from smoothed model
        % 2. Name; Name --> use midpoint between positions ""
        % 3. Name; Name; Double --> use percentage difference b/w positions ""
        offset_time = str2num(eggshell_config_info{i, 1});
        empty_vert_entry_it = -1;
        positions_mat = cell(str2num(eggshell_config_info{i, 4}), str2num(eggshell_config_info{i, 5}), 3);        
        eggshell_names = temporally_smoothed_eggshell_nuc_data(:, 5, ((t - offset_time) + 1));
        for y=1:str2num(eggshell_config_info{i, 4})
            for x=1:str2num(eggshell_config_info{i, 5})
                str = names_mat{y, x};
                
                % check if there are multiple cell names
                if ~isempty(str) && isempty(strfind(str, ';'))
                    % CASE 1 - add the position data
                    lineage_name = str;
                    if ~isempty(strfind(lineage_name, '('))
                        lineage_name = lineage_name(1:strfind(lineage_name, ' (')-1);
                    end
                    
                    idx = find(strcmp(eggshell_names, lineage_name), 1);
                    
                    if isempty(idx)
                        error(['problems finding idx for ' lineage_name]);
                    end
                    
                    % gather the positional data for this nuc
                    x1 = temporally_smoothed_eggshell_nuc_data{idx, 1, (t - offset_time) + 1};
                    y1 = temporally_smoothed_eggshell_nuc_data{idx, 2, (t - offset_time) + 1};
                    z1 = temporally_smoothed_eggshell_nuc_data{idx, 3, (t - offset_time) + 1};
                    
                    % add the data to the positions_mat
                    positions_mat{y, x, 1} = x1;
                    positions_mat{y, x, 2} = y1;
                    positions_mat{y, x, 3} = z1;
            
                elseif ~isempty(str) && ~isempty(strfind(str, ';'))
                    % gather the positional data for the two nucs
                    r = strfind(str, ';');
                    
                    lineage_name_1 = str(1:r(1, 1)-1);
                    if ~isempty(strfind(lineage_name_1, '('))
                        lineage_name_1 = lineage_name_1(1:strfind(lineage_name_1, ' (')-1);
                    end
                        
                    lineage_name_2 = str(r(1,1)+1:r(1,2)-1);
                    if ~isempty(strfind(lineage_name_2, '('))
                        lineage_name_2 = lineage_name_2(1:strfind(lineage_name_2, ' (')-1);
                    end
                        
                    idx_1 = find(strcmp(eggshell_names, lineage_name_1), 1);
                    if isempty(idx_1)
                        error(['problems finding idx for ' lineage_name_1]);
                    end
                        
                    idx_2 = find(strcmp(eggshell_names, lineage_name_2), 1);
                    if isempty(idx_2)
                        error(['problems finding idx for ' lineage_name_2]);
                    end
                        
                    % gather the positional data for the nucs
                    x1 = temporally_smoothed_eggshell_nuc_data{idx_1, 1, (t - offset_time) + 1};
                    y1 = temporally_smoothed_eggshell_nuc_data{idx_1, 2, (t - offset_time) + 1};
                    z1 = temporally_smoothed_eggshell_nuc_data{idx_1, 3, (t - offset_time) + 1};
                        
                    x2 = temporally_smoothed_eggshell_nuc_data{idx_2, 1, (t - offset_time) + 1};
                    y2 = temporally_smoothed_eggshell_nuc_data{idx_2, 2, (t - offset_time) + 1};
                    z2 = temporally_smoothed_eggshell_nuc_data{idx_2, 3, (t - offset_time) + 1};
                    
                    % if there is just one semicolon
                    if size(r, 2) == 1
                        % CASE 2 - add the midpoint between the two positions
                        
                        
                        % compute the midpoint
                        x_mid = (x1 + x2) / 2;
                        y_mid = (y1 + y2) / 2;
                        x_mid = (z1 + z2) / 2;
                        
                        % add the data to the positions mat
                        positions_mat{y, x, 1} = x_mid;
                        positions_mat{y, x, 2} = y_mid;
                        positions_mat{y, x, 3} = z_mid;
                        
                    elseif size(r, 2) == 2
                        % CASE 3 - compute specified distance between points
                        
                        % get the distance percentage value
                        dist_perc = str2num(str(r(1,2)+1:end));
                        
                        % calculate the distance between two points
                        x_sqr = (x2 - x1).^2;
                        y_sqr = (y2 - y1).^2;
                        z_sqr = (z2 - z1).^2;
                        sum = x_sqr + y_sqr + z_sqr;
                        d = sqrt(sum);
                        
                        % compute the desired point
                        x_perc_dist = ((1 - dist_perc)*x1) + (dist_perc*x2);
                        y_perc_dist = ((1 - dist_perc)*y1) + (dist_perc*y2);
                        z_perc_dist = ((1 - dist_perc)*z1) + (dist_perc*z2);
                        
                        % add the data to the positions mat
                        positions_mat{y, x, 1} = x_perc_dist;
                        positions_mat{y, x, 2} = y_perc_dist;
                        positions_mat{y, x, 3} = z_perc_dist;
                    end
                    
                else
                        % add empty identifiers to the positions mat
                        positions_mat{y, x, 1} = empty_vert_entry_it;
                        positions_mat{y, x, 2} = empty_vert_entry_it;
                        positions_mat{y, x, 3} = empty_vert_entry_it;
                        empty_vert_entry_it = empty_vert_entry_it - 1;
                end    
            end
        end

        % connect up the positions using the deterministic algorithm as
        % follows:
        
        % first, make the vertex list and initialize the faces list
        vertices = zeros((size(positions_mat, 1) * size(positions_mat, 2)), 3);
        it = 1;
        for y=1:size(positions_mat, 1)
            for x=1:size(positions_mat, 2)
               vertices(it, 1) = positions_mat{y, x, 1};
               vertices(it, 2) = positions_mat{y, x, 2};
               vertices(it, 3) = positions_mat{y, x, 3};
               it = it + 1;
            end
        end
        faces = zeros(((size(positions_mat, 1) * size(positions_mat, 2)) * 5), 3);
        
        % counters for faces list
        it = 1;
        for y=1:size(positions_mat, 1)-1
            for x=1:size(positions_mat, 2)-1
               if positions_mat{y,x,1} <= -1
                   continue;
               end
               
               % idx of the current vertex in the vertices list
               v_idx_curr = ((y-1)*size(positions_mat, 2)) + x;
               
               % clear the iterators
               y_itr = 0;
               x_itr = 0;
               
               % look right of current position for non-empty cell
               for x_itr=x+1:size(positions_mat, 2)
                   if positions_mat{y, x_itr, 1} > -1
                  
                       % find the idx in the vertices list that corresponds to
                       % the position to the right of curr
                       v_idx_right = ((y-1)*size(positions_mat, 2)) + x_itr;
                       break;
                   end
               end
                              
               % look down from current position for non-empty cell
               for y_itr=y+1:size(positions_mat, 1)
                   if positions_mat{y_itr, x, 1} > -1

                       % find the idx in the vertices list that corresponds to
                       % the position down from curr
                       v_idx_down = ((y_itr-1)*size(positions_mat, 2)) + x;
                       break;
                   end
               end
               
               
               % look down and right from current position for non-empty cell
               found_down_right = 0;
               for x_itr=x+1:size(positions_mat, 2)
                  if positions_mat{y_itr, x_itr, 1} > -1
                   
                       % find the idx in the vertices list that corresponds to
                       % the position down and right from curr
                       v_idx_down_right = ((y_itr-1)*size(positions_mat, 2)) + x_itr;
                       found_down_right = 1;
                       break;
                   end 
               end
               
               
               % make sure there's actually something to look up from
               if found_down_right == 0
                   % the rest of the row was empty, so we're looking at the
                   % special case: 
                   % _____
                   % | / |
                   % | \ |
                   % -----
                   
                   % right now, should have *current*, *down*, *right*
                   % need to get *down_down*, *down_down_right*
                   y_itr = y_itr + 1; % move down a row
                   x_itr = x; % move x back to the current column
                   if positions_mat{y_itr, x_itr, 1} > 1
                      % found down_down
                      v_idx_down_down = ((y_itr-1)*size(positions_mat, 2)) + x_itr;
                      
                      % now grab down_down_right
                      x_itr = x_itr + 1;
                      if positions_mat{y_itr, x_itr, 1} > 1
                         % found down_down_right
                         v_idx_down_down_right = ((y_itr-1)*size(positions_mat, 2)) + x_itr;
                         
                         
                         % make 3 special case triangles
                         faces(it, 1) = v_idx_curr;
                         faces(it, 2) = v_idx_down;
                         faces(it, 3) = v_idx_right;
                         
                         it = it + 1;
                         
                         faces(it, 1) = v_idx_right;
                         faces(it, 2) = v_idx_down;
                         faces(it, 3) = v_idx_down_down_right;
                         
                         it = it + 1;
                         
                         faces(it, 1) = v_idx_down;
                         faces(it, 2) = v_idx_down_down;
                         faces(it, 3) = v_idx_down_down_right;
                         
                         it = it + 1;
                      end
                   end
                   
               else
               % look up from the position down_right from curr for
               % non-empty cell
               y_itr = y_itr - 1;
               if positions_mat{y_itr, x_itr, 1} > -1    
                  % find the idx in the vertices list that correpsonds
                  % to the position up from the down_right position
                  % from curr
                  v_idx_up_from_down_right = ((y_itr-1)*size(positions_mat, 2)) + x_itr;
                      
                  % check if v_idx_right is same as v_idx_up_from_down_right
                  % if so, this is the base case square, split it into two
                  % triangles. BASE CASE
                  if v_idx_right == v_idx_up_from_down_right
                      faces(it, 1) = v_idx_curr;
                      faces(it, 2) = v_idx_down;
                      faces(it, 3) = v_idx_right;
                      
                      it = it + 1;
                      
                      faces(it, 1) = v_idx_down;
                      faces(it, 2) = v_idx_down_right;
                      faces(it, 3) = v_idx_right;
                      
                      it = it + 1;
 
                   else
                      % in this situation, there is a valid position. first,
                      % look one more above and see if it is the v_idx_right, and
                      % then just skip over this one for the extended base case
                      y_itr = y_itr - 1;
                      if y_itr > 0 && positions_mat{y_itr, x_itr, 1} > -1
                          v_idx_up_from_down_right = ((y_itr-1)*size(positions_mat, 2)) + x_itr;
                       
                          % check if v_idx_right is same as v_idx_up_from_down_right
                          % if so, this is the extended base case square, split it into two
                          % triangles. EXTENDED BASE CASE
                          if v_idx_right == v_idx_up_from_down_right
                              faces(it, 1) = v_idx_curr;
                              faces(it, 2) = v_idx_down;
                              faces(it, 3) = v_idx_right;
                      
                              it = it + 1;
                      
                              faces(it, 1) = v_idx_down;
                              faces(it, 2) = v_idx_down_right;
                              faces(it, 3) = v_idx_right;
                      
                              it = it + 1;
                          else
                              v_idx_up_from_down_right = ((y_itr)*size(positions_mat, 2)) + x_itr;
                              % if that doesn't work and it's not the same as v_idx_right,
                              % this is a special case and we need three triangles: |/_\|
                              % SPECIAL CASE
                      
                              % makes 3 special case triangles
                              faces(it, 1) = v_idx_curr;
                              faces(it, 2) = v_idx_down;
                              faces(it, 3) = v_idx_right;
                       
                              it = it + 1;
                    
                              faces(it, 1) = v_idx_down;
                              faces(it, 2) = v_idx_right;
                              faces(it, 3) = v_idx_down_right;
                       
                              it = it + 1;
                       
                              faces(it, 1) = v_idx_right;
                              faces(it, 2) = v_idx_down_right;
                              faces(it, 3) = v_idx_up_from_down_right;
                       
                              it = it + 1;
                              break;
                          end
                      end
                  end
               else
                   % in this situation, there is no valid position above,
                   % first, look if there is a valid position one more
                   % space above, and see if it's in the same row as curr,
                   % if so, this is just an extension of the base case
                   % where an empty row is in the middle of the square
                   y_itr = y_itr - 1;
                   if y_itr > 0 && positions_mat{y_itr, x_itr, 1} > -1
                      v_idx_up_from_down_right = ((y_itr-1)*size(positions_mat, 2)) + x_itr;
                       
                      % check if v_idx_right is same as v_idx_up_from_down_right
                      % if so, this is the extended base case square, split it into two
                      % triangles. EXTENDED BASE CASE
                      if v_idx_right == v_idx_up_from_down_right
                          faces(it, 1) = v_idx_curr;
                          faces(it, 2) = v_idx_down;
                          faces(it, 3) = v_idx_right;
                      
                          it = it + 1;
                      
                          faces(it, 1) = v_idx_down;
                          faces(it, 2) = v_idx_down_right;
                          faces(it, 3) = v_idx_right;
                      
                          it = it + 1;
                      end
                   else
                       y_itr = y_itr - 1;
                       % try one more up
                       if y_itr > 0 && positions_mat{y_itr, x_itr, 1} > -1
                          % check if v_idx_right is same as v_idx_up_from_down_right
                          % if so, this is the extended base case square, split it into two
                          % triangles. EXTENDED EXTENDED BASE CASE
                          if v_idx_right == v_idx_up_from_down_right
                              faces(it, 1) = v_idx_curr;
                              faces(it, 2) = v_idx_down;
                              faces(it, 3) = v_idx_right;

                              it = it + 1;

                              faces(it, 1) = v_idx_down;
                              faces(it, 2) = v_idx_down_right;
                              faces(it, 3) = v_idx_right;

                              it = it + 1;
                          end
                       else
                           % if that didn't work,look right from down_right, and if valid, make 3
                           % triangles: |_\/_|
                           % SPECIAL CASE
                           y_itr = y_itr + 2; % move back down
                           x_itr = x_itr + 1; % move one space to the right
                           if positions_mat{y_itr, x_itr, 1} > -1

                               v_idx_right_from_down_right = ((y_itr-1)*size(positions_mat, 2)) + x_itr;

                               % makes 3 special case triangles
                               faces(it, 1) = v_idx_curr;
                               faces(it, 2) = v_idx_down;
                               faces(it, 3) = v_idx_down_right;

                               it = it + 1;

                               faces(it, 1) = v_idx_curr;
                               faces(it, 2) = v_idx_down_right;
                               faces(it, 3) = v_idx_right;

                               it = it + 1;

                               faces(it, 1) = v_idx_down_right;
                               faces(it, 2) = v_idx_right_from_down_right;
                               faces(it, 3) = v_idx_right;

                               it = it + 1;
                           end
                       end
                   end
               end
               end
            end 
        end
        
        % remove empty rows from faces list
        faces = faces(any(faces, 2),:);
        
        
        % 'remove' redundant vertices by updating all face ptrs to 
        % redundant vertices to their first instance
        
        % indices of redundant rows
        [~, ind] = unique(vertices(:,:), 'rows');
        sort(ind);
        
        % iterate over vertices
        for itr = 1:size(vertices, 1)
           % check if row vertex is repeat by checking if index is in unique indices list
           ismem_flag = ismember(itr, ind);
           if ~ismem_flag
               % find the first occurence of this vertex
               [~, instances] = ismember(vertices, vertices(itr, :), 'rows');
               first_instance_idx = find(instances, 1, 'first');
               
               % change all ptrs at itr in the face list to first_instance_idx
               faces(faces==itr)=first_instance_idx;
           end
        end
        
        
        % ---- compute per vertex normals ----
        TR = triangulation(faces, vertices);
        u_vert_norms = vertexNormal(TR);
        
        % extrude each vertex in the direction of its normal
        % overcompensate because the smoothing operation will
        % significantly dilute the shape
        translation_amount = 15;
        extruded_verts = zeros(size(vertices));
        for vert_it=1:size(vertices, 1)
            if vertices(vert_it, 1) < 0
                extruded_verts(vert_it, :) = vertices(vert_it, :);
                continue;
            end
            
            xpos = vertices(vert_it, 1);
            ypos = vertices(vert_it, 2);
            zpos = vertices(vert_it, 3);
            
            xdir = u_vert_norms(vert_it, 1);
            ydir = u_vert_norms(vert_it, 2);
            zdir = u_vert_norms(vert_it, 3);
            
            x_extrude = xpos + (translation_amount*xdir);
            y_extrude = ypos + (translation_amount*ydir);
            z_extrude = zpos + (translation_amount*zdir);
            
            extruded_verts(vert_it, 1) = x_extrude;
            extruded_verts(vert_it, 2) = y_extrude;
            extruded_verts(vert_it, 3) = z_extrude;
        end
        
        % smooth the mesh
%         num_smoothing_iterations = 1;
%         smoothed_v = extruded_verts;
%         for p=1:num_smoothing_iterations
%             smoothed_v = lpflow_trismooth(smoothed_v, faces);
%         end
%         fprintf('smoothed mesh\n');
%         toc;
        

        
        % ensure correct poly winding
        [corrected_faces, ~] = unifyMeshNormals(faces, extruded_verts, 'alignTo', 'out');
        
        
        % save to obj file with time offset
        filename = sprintf('%s%s%s%s', output_path, eggshell_str, num2str(t - offset), obj_ext_str);
        saveObjFile(filename, extruded_verts, corrected_faces);
    end
end
        
end