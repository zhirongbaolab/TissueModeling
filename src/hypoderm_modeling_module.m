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
    names_mat = cell(str2num(hypoderm_config_info{i, 4}), str2num(hypoderm_config_info{i, 5}));
    names_list = {}; % for easier access later
    fid = fopen(hypoderm_config_info{i, 3});
    if fid < 0
        error(['could not open file: ' hypoderm_config_info{i, 3}]);
    end
    
    k = 1;
    j = 1;
    while ~feof(fid)
        line = fgetl(fid);
        
        tokens = regexp(line, ',', 'split');
        
        if size(tokens, 2) == str2num(hypoderm_config_info{i, 5})
          if k <= str2num(hypoderm_config_info{i, 4})
             names_mat(k, :) = tokens; 
             k = k + 1;
             
             for token = tokens
                % need to check if there is a semicolon in the token -->
                % indicates multiple names
                lineage_name = token{:};
                if isempty(strfind(lineage_name, ';'))
                   % just one name, extract lineage name
                   if ~isempty(strfind(lineage_name, '('))
                       lineage_name = lineage_name(1:strfind(lineage_name, ' (')-1);
                   end
                       
                   % add it to list
                   if ~isempty(lineage_name)
                       names_list{j} = lineage_name;
                       j = j + 1;           
                   end
                else
                    % multiple names, separate them
                    mult_names = regexp(token, ';', 'split');
                    
                    % extract lineage names
                    for ln = mult_names{1,1}
                        lineage_name = ln{:};
                        
                        if ~isempty(strfind(lineage_name, '('))
                            lineage_name = lineage_name(1:strfind(lineage_name, ' (')-1);
                        end
                            
                        % make sure it's not a number
                        [~, status] = str2num(lineage_name);
                        
                        % add it to list
                        if ~isempty(lineage_name) && status == 0
                            names_list{j} = lineage_name;
                            j = j + 1;
                        end
                   end
                end
             end
          end
        end
    end
    
    % use the names to access the model data over the entire window
    hypoderm_nuc_data = cell((str2num(hypoderm_config_info{i, 4}) * str2num(hypoderm_config_info{i, 5})),...
        5,...
        (str2num(hypoderm_config_info{i, 2}) - str2num(hypoderm_config_info{i, 1}))); 
    % row, col, frame (i.e. cells - rowsxcols in names_mat, data fields - 5, time - end_time-start_time)
    
    % iterate over model data to fill out hypoderm_nuc_data
    for t = str2num(hypoderm_config_info{i, 1}):str2num(hypoderm_config_info{i, 2})
        it = 1;
        % access the nuclei info at the time point
        cell_data_atT = embinfo(t).celldata;
        cell_names = embinfo(t).cellnames;

        % make sure the arrays are parallel in the number of rows
        if size(cell_data_atT, 1) == size(cell_names, 1)
            % iterate over all the cell names
            for g=1:size(cell_names, 1)
                name = cell_names{g, 1};

                % check if the lineage name is part of the hypoderm in this
                % configuration
                idx = find(strcmp(names_list, name), 1);
                if ~isempty(idx)
                    % fprintf('%s, %d\n',name, idx);
                    % use the index of this cell name in the embryo info 
                    % to access cell info in parallel celldata list
                    cell_data = cell_data_atT(g, :);
                    x = cell_data(4);
                    y = cell_data(5);
                    z = cell_data(6);
                    diam = cell_data(7);
                    C = {x,y,z,diam,name};

                    % add the coordinates, diameter, and name to the
                    % pharynx_nuc_data matrix
                    hypoderm_nuc_data(it, :, (t - str2num(hypoderm_config_info{i, 1}) + 1)) = C; % entry row, data fields, time frame
                    it = it + 1;
                end
            end
        end
    end
    
    
    % temporally smooth the model data in this window
    temporally_smoothed_hypoderm_nuc_data = conv_temp_smoothing(hypoderm_nuc_data, 2);
    
    % for each frame in the current window
    for t = str2num(hypoderm_config_info{i, 1}):str2num(hypoderm_config_info{i, 2})
        % parse for syntax triggering computation:
        % 1. Name --> use position from smoothed model
        % 2. Name; Name --> use midpoint between positions ""
        % 3. Name; Name; Double --> use percentage difference b/w positions ""
        positions_mat = cell(str2num(hypoderm_config_info{i, 4}), str2num(hypoderm_config_info{i, 5}), 3);        
        hyp_names = temporally_smoothed_hypoderm_nuc_data(:, 5, (t - str2num(hypoderm_config_info{i, 1}) + 1));
        for y=1:str2num(hypoderm_config_info{i, 4})
            for x=1:str2num(hypoderm_config_info{i, 5})
                str = names_mat{y, x};
                
                % check if there are multiple cell names
                if ~isempty(str) && isempty(strfind(str, ';'))
                    % CASE 1
                    lineage_name = str;
                    if ~isempty(strfind(lineage_name, '('))
                        lineage_name = lineage_name(1:strfind(lineage_name, ' (')-1);
                    end
                    
                    idx = find(strcmp(hyp_names, lineage_name), 1);
                    
                    if isempty(idx)
                        error(['problems finding idx for ' lineage_name]);
                    end
                    
                    % query the positional data for this nuc
                    
                    % add the data to the positions_mat
                        
                elseif ~isempty(str) && ~isempty(strfind(str, ';'))
                    r = strfind(str, ';');
                    % if there is just one semicolon
                    if size(r, 2) == 1
                        % CASE 2
                        lineage_name_1 = str(1:r{1, 1});
                        lineage_name_1 = lineage_name_1(1:strfind(lineage_name_1, ' (')-1);
                        
                        lineage
                        
                        idx = find(strcmp(hyp_names, lineage_name_1), 1);
                        
                    elseif size(r, 2) == 2
                        % CASE 3  
                        
                        
                    end
                    
                end
                
            end
        end
        
    

        % call the connectivity module to create a mesh
        
        % subdivide the mesh
        
        % smooth the mesh
        
        % save to obj file with time offset
    end
end
        
end

