function [spherically_expanded_pc] = spherical_expansion(pt_cloud_data)
% expand the point cloud by treating nuclei as the center of spheres
% and expanding them by a given radius
fprintf('expanding pt cloud by spherical interpolating each nuc\n');
tic;

% helper vars
offset = 1;
x_idx = 1;
y_idx = 2;
z_idx = 3;
diam_idx = 4;
name_idx = 5;
num_coords = 3;
membrane_offset = 2.81;
m5R_membrane_offset = 1.2;
m5R_lineage_name = 'ABaraapppap';

% this is a parallel data structure to the original pt_cloud, with the
% addition of every one nucleus in the original representing a 20x3 matrix
% in this data structure (3-x,y,z coords / 20 is max diameter expected
% which corresponds to number of coordinates generate by the spherical
% expansion. Dimenions = 3x20xmax_num_cellsxnum_frames
spherically_expanded_pc = cell(500, num_coords, size(pt_cloud_data, 1), size(pt_cloud_data, 3));

% iterate over all frames
for t=1:size(pt_cloud_data, 3)
    % iterate over all cells
    for i=1:size(pt_cloud_data, 1)
       % check if the cell array here is empty
       if isempty(pt_cloud_data{i, 1, t})
           break;
       end
       
       % access the cell data
       x = pt_cloud_data{i, x_idx, t};
       y = pt_cloud_data{i, y_idx, t};
       z = pt_cloud_data{i, z_idx, t};
       diam = pt_cloud_data{i, diam_idx, t};
       name = pt_cloud_data{i, name_idx, t};
       
       % generate a sphere with the given diameter as the radius (so that
       % it extends beyond the area occupied by the nucleus so as to better
       % approximate the membrane
       [x_coords_sphere, y_coords_sphere, z_coords_sphere] = sphere(5);
       x_coords_sphere = x_coords_sphere(:);
       y_coords_sphere = y_coords_sphere(:);
       z_coords_sphere = z_coords_sphere(:);
       
       % offset all generic x,y,z coordinates of the sphere with the nuc
       % coords so that the nucleus is the center of this sphere
       xyz_coords_sphere_aligned = cell(size(x_coords_sphere, 1), num_coords);
       
       if strcmp(name, m5R_lineage_name)
          for k=1:size(x_coords_sphere, 1)
            xyz_coords_sphere_aligned{k, x_idx} = (x_coords_sphere(k)*((diam/2.0) + m5R_membrane_offset)) + x;
            xyz_coords_sphere_aligned{k, y_idx} = (y_coords_sphere(k)*((diam/2.0) + m5R_membrane_offset)) + y;
            xyz_coords_sphere_aligned{k, z_idx} = (z_coords_sphere(k)*((diam/2.0) + m5R_membrane_offset)) + z;
          end
       else
          for k=1:size(x_coords_sphere, 1)
            xyz_coords_sphere_aligned{k, x_idx} = (x_coords_sphere(k)*((diam/2.0) + membrane_offset)) + x;
            xyz_coords_sphere_aligned{k, y_idx} = (y_coords_sphere(k)*((diam/2.0) + membrane_offset)) + y;
            xyz_coords_sphere_aligned{k, z_idx} = (z_coords_sphere(k)*((diam/2.0) + membrane_offset)) + z;
          end
       end
       
       
       % add these coordinates to the matrix corresponding to the current
       % nuc in the data structure spherically_expanded_pc(:, :, i, t)
       spherically_expanded_pc(1:size(xyz_coords_sphere_aligned, 1), 1:num_coords,...
           i, t) = xyz_coords_sphere_aligned;
    end
end

fprintf('spherical expansion complete\n');
toc;
fprintf('\n');
end

