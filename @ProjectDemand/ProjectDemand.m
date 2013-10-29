classdef ProjectDemand < handle
    %ProjectDemand stores and outputs the projected GPCD for the
    %weather-integrated water model
    %   Detailed explanation goes here
    
    properties (GetAccess=public, SetAccess=private)
        gpcdMatrix
        years
        models
    end
    
    methods
        function ReadGPCD( obj, file )
            fid = fopen(file, 'r');
            tline = fgetl(fid);
            
            % Split header
            A(1,:) = regexp(tline, '\,', 'split');
            
            % Parse and read rest of file
            ctr = 1;
            while(~feof(fid))
                if ischar(tline)
                    ctr = ctr + 1;
                    tline = fgetl(fid);
                    A(ctr,:) = regexp(tline, '\,', 'split');
                else
                    break;
                end
            end
            fclose(fid);
            
            % Split into model names, years, and GPCD demands
            obj.gpcdMatrix = cell2mat(cellfun(@str2num, A(2:end,2:end), ...
                'UniformOutput', false));
            obj.years = cell2mat(cellfun(@str2num, A(2:end,1), ...
                'UniformOutput', false));
            obj.models = A(1,2:end);
        end
        
        
        function gpcd = GPCD( obj, model, period )
            gpcd = obj.gpcdMatrix(period, model);
        end
    end
    
end

