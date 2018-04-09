function [] = hypoderm_modeling_module(embinfo, config_info)

% read the hypoderm configuation file - access the time windows and
% resource locations



% iterate over the time window segments

    % read the names matrix [nxmx1] --> columns, rows, name
    
    % for each frame in the current window:

        % parse for syntax triggering computation:
        % 1. Name --> use position
        % 2. Name; Name --> use midpoint between positions
        % 3. Name; Name; Double --> use percentage difference between positions

        % create a positions matrix [nxmx3] --> columns, rows, position in 3D
        % iterate over all cells:
    

        % call the connectivity module to create a mesh
        
        % subdivide the mesh
        
        % smooth the mesh
        
        % save to obj file
        
end

