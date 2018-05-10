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

shapes = cell(size(pt_cloud_data, 3), 1);

figure(1);
figure(2);
figure(3);
figure(4);
for t=1:size(pt_cloud_data, 3)
   % format the spherical data so that all points are combined: dimensions
   % of pt_cloud are num_pointsXnum_coordsXnum_nucsXnum_frames
   xyz_coords = cell(size(pt_cloud_data, 1), num_coords);
   it = 1;
      for k=1:size(pt_cloud_data, 1)
         x = pt_cloud_data(k, x_idx, t);
         y = pt_cloud_data(k, y_idx, t);
         z = pt_cloud_data(k, z_idx, t);
         
         xyz_coords(it, x_idx) = x;
         xyz_coords(it, y_idx) = y;
         xyz_coords(it, z_idx) = z;
         it = it + 1;
      end
   
   xyz_coords = cell2mat(xyz_coords);
   xyz_coords = unique(xyz_coords, 'rows');

   tic;
   shp = alphaShape(xyz_coords);
   pc = criticalAlpha(shp,'one-region');
   if t <= 10
      shp.Alpha = pc*40;
   elseif t > 10 && t <= 12
       shp.Alpha = pc*18;
   elseif t > 12
       shp.Alpha = pc*8;
   end
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


    % extrude the shape
    vertices = p.Vertices;
    faces = p.Faces;
    TR = triangulation(faces, vertices);
    u_vert_norms = vertexNormal(TR);
    translation_amount = 13;
%     if t == 1
%         translation_amount = 15;
%     elseif t == 2
%       translation_amount = 13;
%     elseif t == 3
%       translation_amount = 12;
%     elseif t == 4
%       translation_amount = 10;
%     elseif t == 5
%       translation_amount = 9;
%     elseif t == 6
%       translation_amount = 8;
%     elseif t == 7
%       translation_amount = 6;  
%     elseif t == 8
%       translation_amount = 5;
%     elseif t == 9
%         translation_amount = 4;
%     elseif t == 10
%        translation_amount = 3;
%     elseif t == 11
%         translation_amount = 2;
%     elseif t >= 12
%        translation_amount = 1;
%     end
   
    
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

    % use the Loop subdivision algorithm to subdivide the surface mesh,
    % creating more uniform parameterization (due to the spherical
    % expansion algorith, the mesh is highly parameterized at this point,
    % which high poly count at spheres and large polys among spheres
    tic;
    [uniform_v, uniform_f] = LoopSubdivisionLimited(extruded_verts, faces, 1);
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
    num_smoothing_iterations = 1;
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
   tic;
   [uniform_smoothed_reduced_f, uniform_smoothed_reduced_v] = ...
       reducepatch(uniform_f,...
       smoothed_uniform_v, .7);
   fprintf('reduced the poly count');
   toc;
    
    % validate and correct (if necessary) the winding order of the faces
%     uniform_f = correct_poly_winding(...
%         uniform_f, smoothed_uniform_v);
    
    % correct any faces that may have inward normals
    [corrected_faces, ~] = unifyMeshNormals(uniform_smoothed_reduced_f, ...
        uniform_smoothed_reduced_v, 'alignTo', 'out');
    
    % make faces and vertices into cell object
    C = {corrected_faces, uniform_smoothed_reduced_v};
    
    % plot uniform, smoothed, reduced shape
    figure(4);
    clf(figure(4));
    patch('Faces', corrected_faces, ...
        'Vertices', uniform_smoothed_reduced_v,...
        'FaceColor', 'red')
    
    % add shape to cell array
    shapes{t, 1} = C;
end
       
end