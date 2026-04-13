function config = parse_params_file(filePath)

%here I only get the basic data info, that will help me read the .dat files

  % Initialize an empty structure to hold the configuration
    config = struct();
    
    % Open the file for reading
    fid = fopen(filePath, 'r');
    if fid == -1
        error('Cannot open file: %s', filePath);
    end
    
    inDataSection = false;
    
    % Read the file line by line
    while ~feof(fid)
        line = strtrim(fgetl(fid));  % Read a line and trim whitespace
        
        % Check for section header
%         if strcmp(line, '[data]')
%             inDataSection = true;
%             continue;
%         end
%         
%         % If we are in the data section, parse key-value pairs
%         if inDataSection
            if startsWith(line, '###') || isempty(line)  % Ignore comments and empty lines
                continue;
            end
            
            % Split the line into key and value
            parts = strsplit(line, '=');
            if numel(parts) == 2
                key = strtrim(parts{1});
                value = strtrim(parts{2});
                
                % Remove comments from value
                if contains(value, '#')
                    parts2=strsplit(value,'#');
                   value = strtrim(parts2{1});
                end
                
                % Store the key-value pair in the structure
                config.(genvarname(key)) = value;  % Use genvarname to ensure valid field names
            end
%         end
    end
    
    % Close the file
    fclose(fid);
end
