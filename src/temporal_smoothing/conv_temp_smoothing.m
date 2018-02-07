function [temporally_smoothed_pt_cloud] = conv_temp_smoothing(original_pt_cloud, offset)
fprintf('performing convolutional temporal smoothing using %dx1 moving average kernel\n', ((offset*2) + 1));
tic;

% helper vars
x_idx = 1;
y_idx = 2;
z_idx = 3;
diam_idx = 4;
name_idx = 5;

temporally_smoothed_pt_cloud = cell(size(original_pt_cloud));

% use an offsetx1 box filter to perform a convolution & temporally smooth the data
 for t=1:size(original_pt_cloud, 3)
        
    window_start = t;
    window_end = t;
    
    if t - offset > 0 && t + offset <= size(original_pt_cloud, 3)
        % all  window positions valid
        window_start = t - offset;
        window_end = t + offset;
    elseif t - offset <= 0 && t + offset <= size(original_pt_cloud, 3)
        % t - offset reaches invalid time points 0 or less
        window_start = 1;
        window_end = t + offset;
    elseif t - offset > 0 && t + offset > size(original_pt_cloud, 3)
        % t + offset reaches invalid time points past movie frames
        window_start = t - offset;
        window_end = size(original_pt_cloud, 3);
    end
        
    % iterate over the pts
    for i=1:size(original_pt_cloud, 1)

        % reset vars
        x = zeros((window_end - window_start) + 1, 1);
        y = zeros((window_end - window_start) + 1, 1);
        z = zeros((window_end - window_start) + 1, 1);

        % access the cell name and its positional data at this time
        name = original_pt_cloud{i, name_idx, t};
        if isempty(name)
           break; 
        end

        % iterate over the full window
        it = 1;
        num_valid_frames = 0;
        for k=window_start:window_end
            % check if the cell was alive at current time frame in the
            % window
            names_at_k = original_pt_cloud(1:size(original_pt_cloud, 1), name_idx, k);
            names_at_k = names_at_k(~cellfun(@isempty, names_at_k)); % remove empty cells
            idx = find(strcmp(names_at_k, name), 1);
            if ~isempty(idx)
                % save previous position data
                x(it) = original_pt_cloud{idx, x_idx, k};
                y(it) = original_pt_cloud{idx, y_idx, k};
                z(it) = original_pt_cloud{idx, z_idx, k};
                num_valid_frames = num_valid_frames + 1;
            else
                if k < t
                    fprintf('%s just born at t=%d\n', name, k);
                elseif k > t 
                    fprintf('%s divided at t=%d\n', name, k);
                    break;
                end
            end

            it = it + 1;
        end

        % average the positions x,y,z
        avg_x = sum(x) / num_valid_frames;
        avg_y = sum(y) / num_valid_frames;
        avg_z = sum(z) / num_valid_frames;

        % add mean values to temp smoothed matrix at (i, [x,y,z], t)
       C = {avg_x,avg_y, avg_z,original_pt_cloud{i, diam_idx, t},name};
       temporally_smoothed_pt_cloud(i, :, t) = C;
    end    
 end
 fprintf('convolutional temporal smoothing complete\n');
 toc;
 fprintf('\n');
end