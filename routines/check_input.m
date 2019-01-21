function [ec, se]= check_input( lat_data0, lon_data0, orient_data0, method , depth_data, quality , plate , pb_dist, pbe , regime);
% check_input
%   This function evaluates whether the crucial input data is available and
%   returns messages to find the errors.
ec = 0;
se = 0;
for i = 1:length(lat_data0)
    % Hard errors
    if lat_data0(i) > 90 || lat_data0(i) < -90 || isnan(lat_data0(i))
        error = ['Fatal! Check latitude in line ',num2str(i)];
        disp(error)
        ec = ec + 1;
    end
    
    if lon_data0(i) > 180 || lon_data0(i) < -180 || isnan(lon_data0(i))
        error = ['Fatal! Check longitude in line ',num2str(i)];
        disp(error)
        ec = ec + 1;
    end
    
    if orient_data0(i) > 180 || orient_data0(i) < 0 || isnan(orient_data0(i))
        if orient_data0(i) ~= 999
            error = ['Fatal! Check Azimuth in line ',num2str(i)];
            disp(error)
            ec = ec + 1;
        end
    end
    
    %Soft errors
    if ischar(method{i})
    else
        error = ['Check method in line ',num2str(i)];
        disp(error)
        se = se + 1;
    end
    
    if depth_data(i) < 0 || depth_data(i) > 40 || isnan(depth_data(i))
        error = ['Check depth in line ',num2str(i)];
        disp(error)
        se = se + 1;
    end
    
    if isequal(quality{i},'A') || isequal(quality{i},'B') || isequal(quality{i},'C') || ...
            isequal(quality{i},'D') || isequal(quality{i},'E')
    else
        error = ['Fatal! Check quality in line ',num2str(i)];
        disp(error)
        ec = ec + 1;
    end
    
    if isequal(regime{i},'NF') || isequal(regime{i},'NS') || isequal(regime{i},'SS') || ...
            isequal(regime{i},'TS') || isequal(regime{i},'TF') || isequal(regime{i},'U')
    else
        error = ['Check regime in line ',num2str(i)];
        disp(error)
        se = se + 1;
    end
    
    if ischar(plate{i})
    else
        error = ['No plate assigned in line ',num2str(i)];
        disp(error)
        se = se + 1;
    end
    
    if isnan(pb_dist(i)) || pb_dist(i) < 0
        error = ['Check distance to plate boundary in line ',num2str(i)];
        disp(error)
        se = se + 1;
    end
    
    if ischar(pbe{i})
    else
        error = ['Fatal! PBE/No PBE not specified in line ',num2str(i)];
        disp(error)
        ec = ec + 1;
    end
    
end

end

