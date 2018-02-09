function [shapes] = compute_alpha_shape(pt_cloud_data)
% ALPHA SHAPE

% ----------------------------------------------
% MATLAB DOCS
% Shp = alphaShape(x,y,z,a)
% creates a bounding area or volume 
% that envelops a set of 2-D or 3-D points. 
% You can manipulate the alphaShape object to 
% tighten or loosen the fit around the points to 
% create a nonconvex region. You also can add 
% or remove points or suppress holes or regions.
%
% x, y, z - positions in 3 dimensions
% a - alpha radius, controls the convexity and
%  how tight the shape should fit around the pts
% ----------------------------------------------

% helper vars
x_idx = 1;
y_idx = 2;
z_idx = 3;
num_coords = 3;

shapes = cell(size(pt_cloud_data, 4), 1);

figure(1);
figure(2);
figure(3);
for t=1:size(pt_cloud_data, 4)
   % format the spherical data so that all points are combined: dimensions
   % of pt_cloud are num_pointsXnum_coordsXnum_nucsXnum_frames
   xyz_coords = cell(size(pt_cloud_data, 1) * size(pt_cloud_data, 3), num_coords);
   it = 1;
   for i=1:size(pt_cloud_data, 3)
      for k=1:size(pt_cloud_data, 1)
         x = pt_cloud_data(k, x_idx, i, t);
         y = pt_cloud_data(k, y_idx, i, t);
         z = pt_cloud_data(k, z_idx, i, t);
         
         xyz_coords(it, x_idx) = x;
         xyz_coords(it, y_idx) = y;
         xyz_coords(it, z_idx) = z;
         it = it + 1;
      end
   end
   
   xyz_coords = cell2mat(xyz_coords);
   xyz_coords = unique(xyz_coords, 'rows');

   shp = alphaShape(xyz_coords);
   pc = criticalAlpha(shp,'one-region');
   shp.Alpha = pc*15;
   
   % plot alpha shape
   figure(1);
   clf(figure(1));
   p = plot(shp);
   
   disp(t);
   pause(.1);
   
    % plot alpha shape wire frame
%     figure(2);
%     h = plot(shp);
%     set(h, 'FaceColor', 'none');
%     scatter3(xyz_coords(:, x_idx), xyz_coords(:, y_idx), xyz_coords(:, z_idx));

    % use the Loop subdivision algorithm to subdivide the surface mesh,
    % creating more uniform parameterization (due to the spherical
    % expansion algorith, the mesh is highly parameterized at this point,
    % which high poly count at spheres and large polys among spheres
    [uniform_v, uniform_f] = LoopSubdivisionLimited(p.Vertices, p.Faces, 2);
    
    % smooth the triangulated mesh to remove local variations
    fv_struct.faces = uniform_f;
    fv_struct.vertices = uniform_v;
    fv_struct_smoothed = smoothpatch(fv_struct, 0, 1, 1, 1);
    
    % plot reduced uniformly parameterized shape
    figure(2);
    clf(figure(2));
    patch('Faces', fv_struct_smoothed.faces,...
        'Vertices', fv_struct_smoothed.vertices,...
        'FaceColor', 'red')

    [uniform_smoothed_reduced_f, uniform_smoothed_reduced_v] = ...
        reducepatch(fv_struct_smoothed.faces,...
        fv_struct_smoothed.vertices, .1);
    
    % make faces and vertices into cell object
    C = {uniform_smoothed_reduced_f, uniform_smoothed_reduced_v};
    
    % plot uniform reduced face shape
    figure(3);
    clf(figure(3));
    patch('Faces', uniform_smoothed_reduced_f, ...
        'Vertices', uniform_smoothed_reduced_v,...
        'FaceColor', 'red')
    
    % add shape to cell array
    shapes{t, 1} = C;
   
end
       
end