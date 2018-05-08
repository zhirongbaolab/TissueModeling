% mesh_shape_metamorphosis.m
% Zhirong Bao Lab, Sloan-Kettering Institute
% Author: Braden Katzman
% Created On: May 8, 2018

% shape metaphorsis between meshes over X number of frames. works by
% establishing dense correspondence and interpolating between frames

% vars
start_frame_mesh_filename = 'C:\Users\katzmanb\Desktop\TissueModeling\data\output\embryo\embryo_t265.obj';
end_frame_mesh_filename = 'C:\Users\katzmanb\Desktop\TissueModeling\data\output\embryo\embryo_t289.obj';

start_time = 265;
end_time = 289;
output_path = 'C:\Users\katzmanb\Desktop\TissueModeling\data\output\embryo\';
embryo_str = 'embryo_t';
obj_ext_str = '.obj';


verts_1 = {};
faces_1 = {};
fid = fopen(start_frame_mesh_filename);
if fid < 0
    error(['could no open file: ' start_frame_mesh_filename]);
end
while ~feof(fid)
    line = fgetl(fid);
    
    tokens = regexp(line, ' ', 'split');
    if size(tokens, 2) == 4
       % make sure vertex line
       if strcmp('v', tokens{1, 1})
           v = {str2num(tokens{1, 2}), str2num(tokens{1, 3}), str2num(tokens{1, 4})};
           verts_1{end+1} = v;
       elseif strcmp('f', tokens{1, 1})
           f = {str2num(tokens{1, 2}), str2num(tokens{1, 3}), str2num(tokens{1, 4})};
           faces_1{end+1} = f;
       end
    else
        h = 0;
    end
end

verts_1 = verts_1';
faces_1 = faces_1';

verts_2 = {};
faces_2 = {};
fid = fopen(end_frame_mesh_filename);
if fid < 0
    error(['could no open file: ' start_frame_mesh_filename]);
end
while ~feof(fid)
    line = fgetl(fid);
    
    tokens = regexp(line, ' ', 'split');
    if size(tokens, 2) == 4
       % make sure vertex line
       if strcmp('v', tokens{1, 1})
           v = {str2num(tokens{1, 2}), str2num(tokens{1, 3}), str2num(tokens{1, 4})};
           verts_2{end+1} = v;
       elseif strcmp('f', tokens{1, 1})
           f = {str2num(tokens{1, 2}), str2num(tokens{1, 3}), str2num(tokens{1, 4})};
           faces_2{end+1} = f;
       end
    else
        h = 0;
    end
end

verts_2 = verts_2';
faces_2 = faces_2';

% create matrix of size verts_1 in y and verts_2 in x
vert_distance_map = zeros(size(verts_1, 1), size(verts_2, 1));

% iterate over verts_1 (because it is larger) and compute its distance to
% every vertex in verts_2
for verts_1_itr=1:size(verts_1, 1)
    x1 = verts_1{verts_1_itr, 1}{1};
    y1 = verts_1{verts_1_itr, 1}{2};
    z1 = verts_1{verts_1_itr, 1}{3};
    for verts_2_itr=1:size(verts_2, 1)
       x2 = verts_2{verts_2_itr, 1}{1};
       y2 = verts_2{verts_2_itr, 1}{2};
       z2 = verts_2{verts_2_itr, 1}{3};
       
       % compute distance
       points_mat = [x1 y1 z1; x2 y2 z2];
       d = pdist(points_mat);
        
       % insert distance calc in map
       vert_distance_map(verts_1_itr, verts_2_itr) = d;
    end
end

% iterate over verts_2, and find the index of the nearest neighbor in verts_1
% of every vert in verts_2
neighbor_indices = zeros(size(verts_2, 1), 1); % parallel in size to record every neighbor
used_vertices = zeros(size(verts_1, 1), 1);
for itr=1:size(verts_2, 1)
   % find the idx of the smallest value in this column of vert_distance_map
   [~, min_dist_idx] = min(vert_distance_map(:, itr));
   
   % make sure the vertex is unused
   if ~used_vertices(min_dist_idx, 1)
       % mark that idx used
       used_vertices(min_dist_idx, 1) = 1;
       
       % record the idx of the neighbor
       neighbor_indices(itr, 1) = min_dist_idx;
   else
       sorted_col = sort(vert_distance_map(:, itr));
       sorted_col_itr = 2; % since we know the first one just didn't work
       found = 0;
       while ~found
          k = find(vert_distance_map(:, itr)==sorted_col(sorted_col_itr));
          if ~used_vertices(k, 1)
              used_vertices(k, 1) = 1;
              neighbor_indices(itr) = k;
              found = 1;
          else
              sorted_col_itr = sorted_col_itr + 1;
          end
       end
   end
end

% iterate over the range between the first frame mesh and the last
% interpolate the values between the two meshes and create an morphed mesh
% at each frame
for t=start_time+1:end_time-1
   morphed_mesh_vertices = zeros(size(verts_2, 1), 3);
   
   for itr=1:size(verts_2, 1)
       % neighbor vert
       x1 = verts_1{neighbor_indices(itr), 1}{1};
       y1 = verts_1{neighbor_indices(itr), 1}{2};
       z1 = verts_1{neighbor_indices(itr), 1}{3};
       
       x2 = verts_2{itr, 1}{1};
       y2 = verts_2{itr, 1}{2};
       z2 = verts_2{itr, 1}{3};
       
       % compute interpolated vertex
       dist_perc = (t - start_time) / (end_time - start_time);
       dist = vert_distance_map(neighbor_indices(itr), itr);
       x_morphed = ((1 - dist_perc)*x1) + (dist_perc*x2);
       y_morphed = ((1 - dist_perc)*y1) + (dist_perc*y2);
       z_morphed = ((1 - dist_perc)*z1) + (dist_perc*z2);
       
       % add to vertices list
       morphed_mesh_vertices(itr, 1) = x_morphed;
       morphed_mesh_vertices(itr, 2) = y_morphed;
       morphed_mesh_vertices(itr, 3) = z_morphed;   
   end
   
   % overwrite .obj file
   
   % convert faces cell to matrix
   faces_2_mat = zeros(size(faces_2, 1), 3);
   for p=1:size(faces_2, 1)
      faces_2_mat(p, :) = [faces_2{p, :}{1} faces_2{p, :}{2} faces_2{p, :}{3}];  
   end
   
   
   filename = sprintf('%s%s%s%s', output_path, embryo_str, num2str(t), obj_ext_str);
   saveObjFile(filename, morphed_mesh_vertices, faces_2_mat);
end