function [] = visualize_pt_cloud(pc_data, msg)
% visualize point cloud data

% helper vars
x_idx = 1;
y_idx = 2;
z_idx = 3;

figure;
for t=1:size(pc_data, 3)
    x_coords = pc_data(1:size(pc_data, 1), x_idx, t);
    x_coords = x_coords(~cellfun(@isempty, x_coords));
    
    y_coords = pc_data(1:size(pc_data, 1), y_idx, t);
    y_coords = y_coords(~cellfun(@isempty, y_coords));
    
    z_coords = pc_data(1:size(pc_data, 1), z_idx, t);
    z_coords = z_coords(~cellfun(@isempty, z_coords));
    
    pcshow([cell2mat(x_coords),cell2mat(y_coords),cell2mat(z_coords)], 'MarkerSize', 200);
    title([msg, ' at t=', num2str(t)]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');   
end

end