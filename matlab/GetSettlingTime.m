function [ settling_time ] = GetSettlingTime( time, signals, targetPoint, eps_r, eps_t )

    settling_time = NaN(size(signals, 2), 1);
    for signal_id = 1:size(signals, 2)
        signal = signals(:, signal_id);
        for i = 1:length(time)
            if (prod(abs(signal(i:end)) < eps_r))
                settling_time(signal_id) = time(i);
                break;
            end;
        end 
    end
    
    % Gets the settling time of the position, involving both x and y
    settling_time([1, 3]) = NaN;
    for i = 1:length(time)
        if (prod( (signals(i:end, 1)-targetPoint(1)).^2 + (signals(i:end, 3)-targetPoint(2)).^2 < eps_t^2) )
           settling_time([1, 3]) = time(i); 
           break;
        end
    end

end

