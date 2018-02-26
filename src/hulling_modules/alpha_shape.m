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
figure(4);
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

   tic;
   shp = alphaShape(xyz_coords);
   pc = criticalAlpha(shp,'one-region');
   shp.Alpha = pc*15;
   fprintf('initial alpha shape computed at time %d\n', t);
   toc;
   
   pause(.1);
   
   % plot the initial alpha shape
   figure(1);
   clf(figure(1));
   p = plot(shp);
   
    % plot alpha shape wire frame
%     figure(2);
%     h = plot(shp);
%     set(h, 'FaceColor', 'none');
%     scatter3(xyz_coords(:, x_idx), xyz_coords(:, y_idx), xyz_coords(:, z_idx));

    % use the Loop subdivision algorithm to subdivide the surface mesh,
    % creating more uniform parameterization (due to the spherical
    % expansion algorith, the mesh is highly parameterized at this point,
    % which high poly count at spheres and large polys among spheres
    tic;
    [uniform_v, uniform_f] = LoopSubdivisionLimited(p.Vertices, p.Faces, 2);
    fprintf('subdivided shape\n');
    toc;
    
    % plot the subdivided (hopefully uniform) shape
    figure(2);
    clf(figure(2));
    patch('Faces', uniform_f,...
        'Vertices', uniform_v,...
        'FaceColor', 'red')
    
    % smooth the more uniform triangulated mesh to remove local variations
    tic;
    num_smoothing_iterations = 3;
    smoothed_uniform_v = uniform_v;
    for i=1:num_smoothing_iterations
        smoothed_uniform_v = lpflow_trismooth(smoothed_uniform_v, uniform_f);
    end
    fprintf('smoothed the mesh\n');
    toc;
    
    % plot the uniform, smoothed mesh
    figure(3);
    clf(figure(3));
    patch('Faces', uniform_f,...
        'Vertices', smoothed_uniform_v,...
        'FaceColor', 'red')
    
    % reduce the uniform, smoothed mesh to a lower poly count shape
%    tic;
%    [uniform_smoothed_reduced_f, uniform_smoothed_reduced_v] = ...
%        reducepatch(uniform_f,...
%        smoothed_uniform_v, .5);
%    fprintf('reduced the poly count');
%    toc;
    
    % validate and correct (if necessary) the winding order of the faces
    uniform_f = correct_poly_winding(...
        uniform_f, smoothed_uniform_v);
    
    % correct any faces that may have inward normals
    [corrected_faces, ~] = unifyMeshNormals(uniform_f, ...
        smoothed_uniform_v, 'alignTo', 'out');
    
    % make faces and vertices into cell object
    C = {corrected_faces, smoothed_uniform_v};
    
    % plot uniform, smoothed, reduced shape
    figure(4);
    clf(figure(4));
    patch('Faces', corrected_faces, ...
        'Vertices', smoothed_uniform_v,...
        'FaceColor', 'red')
    
    % add shape to cell array
    shapes{t, 1} = C;
   
end
       
end