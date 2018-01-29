function [temporally_smoothed_pt_cloud] = conv_temp_smoothing(original_pt_cloud)

% helper vars
offset = 1;
x_idx = 1;
y_idx = 2;
z_idx = 3;
name_idx = 5;

temporally_smoothed_pt_cloud = cell(size(original_pt_cloud));

% use a 3x1 box filter to perform a convolution 
% and temporally smooth the data
 for t=1:size(original_pt_cloud, 3)
    
    % all 3 window positions valid
    if t - offset > 0 && t + offset <= size(original_pt_cloud, 3)
        % iterate over the pts
        for i=1:size(original_pt_cloud, 1)
            % access the cell name and its positional data at this time
            name = original_pt_cloud{i, name_idx, t};
            if isempty(name)
               break; 
            end
            
            % TODO -> catch the cases where a cell has just been born or
            % dies *within* the window
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
                fprintf('%s just born at t=%d', name, t);
                just_born = 1;
            end
            
            % check if the cell is alive at the following time frame
            following_names = original_pt_cloud(1:size(original_pt_cloud, 1), name_idx, t+1);
            following_names = prev_names(~cellfun(@isempty, prev_names)); % remove empty cells
            following_idx = find(strcmp(prev_names, name), 1);
            if ~isempty(prev_idx)
                % save previous position data
                x_3 = original_pt_cloud{following_idx, x_idx, t+1};
                y_3 = original_pt_cloud{following_idx, y_idx, t+1};
                z_3 = original_pt_cloud{following_idx, z_idx, t+1};
            else
                fprintf('%s dibvided at t=%d', name, t+1);
                divided = 1;
            end
            
            % check if cell was just born or divided in this window
            
            % average the positions x,y,z
            
            % add mean values to temp smoothed matrix at (i, [x,y,z], t)
        end
        
    end
    
    % first position invalid, last position valid --> front of list
    
    % first position valid, last position invalid --> back of list
    
 end

 
end