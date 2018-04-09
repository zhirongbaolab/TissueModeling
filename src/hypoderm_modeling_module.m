function [] = hypoderm_modeling_module(embinfo)

% read the hypoderm configuation file: 
% start time, end time, resource location, dimension y (rows), dimensions x (columns), comments
hypoderm_config_info_filename = 'C:\Users\katzmanb\Desktop\TissueModeling\data\configurations\tissue_cells\hypoderm\hypoderm_config.csv';
hypoderm_config_info = cell (12, 6);
fid = fopen(hypoderm_config_info_filename);
if fid < 0
    error(['could not open file: ' hypoderm_config_info_filename]);
end

% get the first line out of the way
line = fgetl(fid);
i = 1;
while ~feof(fid)
    line = fgetl(fid);
    tokens = regexp(line, ',', 'split');
    
    if size(tokens, 2) == 6
       hypoderm_config_info(i, :) = tokens;
    end
   i = i + 1;    
end

% iterate over the time window segments
for i = 1:12
    % read the names matrix [nxmx1] for the current window
    % rows, columns, names
    names_mat = cell(hypoderm_config_info{i, 4}, hypoderm_config_info{i, 5});
    fid = fopen(hypoderm_config_info{i, 3});
    if fid < 0
        error(['could not open file: ' hypoderm_config_info{i, 3}]);
    end
    
    while ~feof(fid)
        line = fgetl(fid);
        
    end
    
    % use the names to access the model data over the entire window
    
    % temporally smooth the model data in this window
    
    % for each frame in the current window
    for t = hypoderm_config_info{i, 1}:hypoderm_config_info{i, 2}
        % parse for syntax triggering computation:
        % 1. Name --> use position from smoothed model
        % 2. Name; Name --> use midpoint between positions ""
        % 3. Name; Name; Double --> use percentage difference b/w positions ""

        % create a positions matrix [nxmx3] --> columns, rows, position in 3D
        % iterate over all cells:
    

        % call the connectivity module to create a mesh
        
        % subdivide the mesh
        
        % smooth the mesh
        
        % save to obj file
    end
end
        
end

