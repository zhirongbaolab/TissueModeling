function [temporally_smoothed_pt_cloud] = conv_temp_smoothing(original_pt_cloud)
fprintf('performing convolutional temporal smoothing using 3x1 moving average kernel');
% helper vars
offset = 1;
x_idx = 1;
y_idx = 2;
z_idx = 3;
diam_idx = 4;
name_idx = 5;

temporally_smoothed_pt_cloud = cell(size(original_pt_cloud));

% use a 3x1 box filter to perform a convolution & temporally smooth the data
 for t=1:size(original_pt_cloud, 3)
    
    % all 3 window positions valid
    if t - offset > 0 && t + offset <= size(original_pt_cloud, 3)
        % iterate over the pts
        for i=1:size(original_pt_cloud, 1)
            
            % reset vars
            x_1 = 0;
            x_2 = 0;
            x_3 = 0;
            y_1 = 0;
            y_2 = 0;
            y_3 = 0;
            z_1 = 0;
            z_2 = 0;
            z_3 = 0;
            
            % access the cell name and its positional data at this time
            name = original_pt_cloud{i, name_idx, t};
            if isempty(name)
               break; 
            end
            
            just_born = 0;
            divided = 0;
            
            x_2 = original_pt_cloud{i, x_idx, t};
            y_2 = original_pt_cloud{i, y_idx, t};
            z_2 = original_pt_cloud{i, z_idx, t};
            
            % check if the cell was alive at the previous time frame
            prev_names = original_pt_cloud(1:size(original_pt_cloud, 1), name_idx, t-1);
            prev_names = prev_names(~cellfun(@isempty, prev_names)); % remove empty cells
            prev_idx = find(strcmp(prev_names, name), 1);
            if ~isempty(prev_idx)
                % save previous position data
                x_1 = original_pt_cloud{prev_idx, x_idx, t-1};
                y_1 = original_pt_cloud{prev_idx, y_idx, t-1};
                z_1 = original_pt_cloud{prev_idx, z_idx, t-1};
            else
                fprintf('%s just born at t=%d\n', name, t);
                just_born = 1;
            end
            
            % check if the cell is alive at the following time frame
            following_names = original_pt_cloud(1:size(original_pt_cloud, 1), name_idx, t+1);
            following_names = following_names(~cellfun(@isempty, following_names)); % remove empty cells
            following_idx = find(strcmp(following_names, name), 1);
            if ~isempty(following_idx)
                % save previous position data
                x_3 = original_pt_cloud{following_idx, x_idx, t+1};
                y_3 = original_pt_cloud{following_idx, y_idx, t+1};
                z_3 = original_pt_cloud{following_idx, z_idx, t+1};
            else
                fprintf('%s divided at t=%d\n', name, t+1);
                divided = 1;
            end
            
            % average the positions x,y,z
            avg_x = 0;
            avg_y = 0;
            avy_z = 0;
            if ~just_born && ~divided % cell alive across all 3 frames
                avg_x = (x_1 + x_2 + x_3) / 3.;
                avg_y = (y_1 + y_2 + y_3) / 3.;
                avg_z = (z_1 + z_2 + z_3) / 3.;
            elseif just_born
                avg_x = (x_2 + x_3) / 2.;
                avg_y = (y_2 + y_3) / 2.;
                avg_z = (z_2 + z_3) / 2.;
            elseif divided % this won't ever happen, but could be useful in another context
                avg_x = (x_1 + x_2) / 2.;
                avg_y = (y_1 + y_2) / 2.;
                avg_z = (z_1 + z_2) / 2.;
            end
            
            % add mean values to temp smoothed matrix at (i, [x,y,z], t)
           C = {avg_x,avg_y, avg_z,original_pt_cloud{i, diam_idx, t},name};
           temporally_smoothed_pt_cloud(i, :, t) = C;
            
        end
        
    end
    
    % first position invalid, last position valid --> front of list
    if t - offset == 0 && t + offset <= size(original_pt_cloud, 3)
        % iterate over the pts
        for i=1:size(original_pt_cloud, 1)
            
            % reset vars
            x_1 = 0;
            x_2 = 0;
            y_1 = 0;
            y_2 = 0;
            z_1 = 0;
            z_2 = 0;
            
            % access the cell name and its positional data at this time
            name = original_pt_cloud{i, name_idx, t};
            if isempty(name)
               break; 
            end
            
            divided = 0;
            
            x_1 = original_pt_cloud{i, x_idx, t};
            y_1 = original_pt_cloud{i, y_idx, t};
            z_1 = original_pt_cloud{i, z_idx, t};
            
            % check if the cell is alive at the following time frame
            following_names = original_pt_cloud(1:size(original_pt_cloud, 1), name_idx, t+1);
            following_names = following_names(~cellfun(@isempty, following_names)); % remove empty cells
            following_idx = find(strcmp(following_names, name), 1);
            if ~isempty(following_idx)
                % save previous position data
                x_2 = original_pt_cloud{following_idx, x_idx, t+1};
                y_2 = original_pt_cloud{following_idx, y_idx, t+1};
                z_2 = original_pt_cloud{following_idx, z_idx, t+1};
            else
                fprintf('%s divided at t=%d\n', name, t+1);
                divided = 1;
            end
            
            % average the positions x,y,z
            avg_x = 0;
            avg_y = 0;
            avg_z = 0;
            if ~divided % cell alive across both frames
                avg_x = (x_1 + x_2) / 2.;
                avg_y = (y_1 + y_2) / 2.;
                avg_z = (z_1 + z_2) / 2.;
            elseif divided % this won't ever happen, but could be useful in another context
                avg_x = x_1;
                avg_y = y_1;
                avg_z = z_1;
            end
            
            % add mean values to temp smoothed matrix at (i, [x,y,z], t)
           C = {avg_x,avg_y, avg_z,original_pt_cloud{i, diam_idx, t},name};
           temporally_smoothed_pt_cloud(i, :, t) = C;
            
        end
        
    end
    
    % first position valid, last position invalid --> back of list
    if t - offset > 0 && t + offset > size(original_pt_cloud, 3)
        % iterate over the pts
        for i=1:size(original_pt_cloud, 1)
            
            % reset vars
            x_1 = 0;
            x_2 = 0;
            y_1 = 0;
            y_2 = 0;
            z_1 = 0;
            z_2 = 0;
            
            % access the cell name and its positional data at this time
            name = original_pt_cloud{i, name_idx, t};
            if isempty(name)
               break; 
            end
            
            just_born = 0;
            
            x_2 = original_pt_cloud{i, x_idx, t};
            y_2 = original_pt_cloud{i, y_idx, t};
            z_2 = original_pt_cloud{i, z_idx, t};
            

            % check if the cell was alive at the previous time frame
            prev_names = original_pt_cloud(1:size(original_pt_cloud, 1), name_idx, t-1);
            prev_names = prev_names(~cellfun(@isempty, prev_names)); % remove empty cells
            prev_idx = find(strcmp(prev_names, name), 1);
            if ~isempty(prev_idx)
                % save previous position data
                x_1 = original_pt_cloud{prev_idx, x_idx, t-1};
                y_1 = original_pt_cloud{prev_idx, y_idx, t-1};
                z_1 = original_pt_cloud{prev_idx, z_idx, t-1};
            else
                fprintf('%s just born at t=%d\n', name, t);
                just_born = 1;
            end
            
            % average the positions x,y,z
            avg_x = 0;
            avg_y = 0;
            avg_z = 0;
            if ~just_born % cell alive across both frames
                avg_x = (x_1 + x_2) / 2.;
                avg_y = (y_1 + y_2) / 2.;
                avg_z = (z_1 + z_2) / 2.;
            elseif just_born % this won't ever happen, but could be useful in another context
                avg_x = x_1;
                avg_y = y_1;
                avg_z = z_1;
            end
            
            % add mean values to temp smoothed matrix at (i, [x,y,z], t)
           C = {avg_x,avg_y, avg_z,original_pt_cloud{i, diam_idx, t},name};
           temporally_smoothed_pt_cloud(i, :, t) = C;
        end   
    end
 end

 
end