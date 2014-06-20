% Returns samples as predicted by model
function [y] = model_samples(t1, sample_times, T, args, prev_data, t, c, sigma)
% t should be absolute time, in [t1, t1+T)

y = zeros(args.L, 1);

if strcmp(args.method_name, 'ker') || strcmp(args.method_name, 'RMSE') || strcmp(args.method_name, 'RMSE_cal')

    try
        % previous spikes
        for l = 1 : args.L
            for k = 1 : length(prev_data.c)
                y(l) = y(l) + prev_data.c(k) * args.h(sample_times(l) - t1 - (prev_data.t(k) - T), args);
                %prev_contr = prev_data.c(k) * args.h(sample_times(l) - t1 - (prev_data.t(k) - T), args)
            end
        end
    catch
        %caught = 1
    end

    % these spikes
    for l = 1 : args.L
        for k = 1 : args.K
            y(l) = y(l) + c(k) * args.h(sample_times(l) - t(k), args);
            %this_contr = c(k) * args.h(sample_times(l) - t(k), args)
        end
    end

    % correction for sigma
    y = y + sigma / sqrt(2 * pi);
    
elseif strcmp(args.method_name, 'gibbs')
    
    for l = 1 : args.L
        for k = 1 : args.K
            y(l) = y(l) + c(k) * gaussian(sample_times(l) - t(k), args);
        end
    end
    
elseif strcmp(args.method_name, 'ML')
    
    try
        % previous spikes
        for l = 1 : args.L
            for k = 1 : length(prev_data.c)
                y(l) = y(l) + prev_data.c(k) * args.h(sample_times(l) - t1 - (prev_data.t(k) - T), args);
                %prev_contr = prev_data.c(k) * args.h(sample_times(l) - t1 - (prev_data.t(k) - T), args)
            end
        end
    catch
        %caught = 1
    end

    % these spikes
    for l = 1 : args.L
        for k = 1 : args.K
            y(l) = y(l) + c(k) * args.h(sample_times(l) - t(k), args);
            %this_contr = c(k) * args.h(sample_times(l) - t(k), args)
        end
    end

    % correction for sigma
    y = y + sigma;  % here 'sigma' is m_w
    
end