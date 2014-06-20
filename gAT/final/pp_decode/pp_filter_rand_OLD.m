function [RMSE_px, RMSE_py, RMSE_vx, RMSE_vy] = pp_filter_rand(data_set, start_time, end_time, num_reaches, spike_times, plot_flag, method_names, delay, b)
% Runs point process filter (random walk prior) taking spike times (constructed using FRI, CS, etc.) as input
% Simultaneously processes multiple sets of spike times, each is denoted as a 'method' (FRI or CS are 2 possible methods for instance)
% Calls Code/Reaching/rand_walk_filter.m

% Inputs:
% start_time, end_time - start/end time of data in which to look for reaches
% num_reaches - number of reaching movements to process
% spike_times(method_num, neuron_num).data = list of spike times
% plot_flag - plots results if set to 1
% method_names - column cell array of method names (i.e. 'FRI', 'CS'), used for plotting
% delay(method_num, neuron_num) = delay (sec)
% b(method_num, neuron_num, :) = b_0, b_1, b_2

% Outputs:
% RMSE_px(method_num, reach_num) = root mean squared error in x position
% RMSE_py(method_num, reach_num) = root mean squared error in y position

% Parameters
dt = get_dt_info();
A = [1, 0, dt, 0; 0, 1, 0, dt; 0, 0, 1, 0; 0, 0, 0, 1];
W_0 = eye(4) * 1;  % Covariance for first step of filter
Q = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1] * 5e-1;  % TODO: Modify random walk update covariance (for decoding) here

% Data
[d, filename] = data_file(data_set, 0, 0); % 'matthew20080721';
load(filename, '-mat');  % Loads time, px, py, vx, vy, etc. (data regarding monkey's movement)

% Loading b parameters
% See Code/Reaching/Hybrid Test/estimate_b.m
% Data saved upon call to test_rand_hybrid.m, for instance, if estimate_b.m is set to save rather than load
b_filename = [root_dir() 'pp_decode/b_temp.mat'];  % TODO: TEMP
load(b_filename);

% Apply delay
C = size(spike_times, 2);  % Number of neurons
M = size(spike_times, 1);  % Number of methods
for c = 1 : C,
    for m = 1 : M,
        spike_times(m, c).data = spike_times(m, c).data + delay(m, c);
    end
end

dt = (time(length(time)) - time(1)) / (length(time) - 1);
used_ind = find(time >= start_time & time < end_time);  % Indices of used data

% Prune data to used indices
time = time(used_ind);
px = px(used_ind);
py = py(used_ind);
vx = vx(used_ind);
vy = vy(used_ind);
cio = cio(used_ind);
cof = cof(used_ind);
tio = tio(used_ind);
tof = tof(used_ind);

% Segment data
% Selection process for start/end of reach can be modified here
%start_ind = find(cio & [diff(cof); 0] == -1);  % Center removed but hand still inside it
start_ind = find(diff(cio) == -1) - 70;  % Leave center - correct offset ~ -70
%end_ind = find([diff(tof); 0] < 0);  % When target disappears, acknowledging arrival
end_ind = find(diff(tio) == 1);  % As soon as hand gets to target

% Convert spike times to spike vector
% spike_vec(method_num, neuron_num, timestep) = 0 or 1 (indicating presence or absence of spike)
for m = 1 : M,
    spike_vec(m, :, :) = spike_times_to_vec(time, spike_times(m, :));  % TODO: allow multiple spikes per bin?
end

for reach_num = 1 : num_reaches,

    try

        % Indices of this reach
        reach_start = start_ind(reach_num);
        ends = find(end_ind > reach_start);
        reach_end = end_ind(ends(1));
        end_ind(ends(1)) = 0;  % Don't use same end again
        starts = find(start_ind < reach_end);
        reach_start = start_ind(starts(length(starts))) - 1;
        reach_ind = reach_start : reach_end;

    catch

        disp('Out of reaches');
        break;

    end

    % Num timesteps
    K = length(reach_ind) - 1;
    
    % Starting state
    start = [px(reach_ind(1)); py(reach_ind(1)); 0; 0];  % Start at a standstill at correct position

    % Decoded trajectories
    % x_decoded(method_num, :, timestep) = state vector [px;py;vx;vy]
    x_decoded = zeros(M, 4, K+1);

    % Point process filter (random walk prior)
    for m = 1 : M,
        [px_r, py_r, vx_r, vy_r] = rand_walk_filter(K, C, dt, start, W_0, A, Q, squeeze(spike_vec(m, :, reach_ind)), squeeze(b(m, :, 1)), squeeze(b(m, :, 2)), squeeze(b(m, :, 3)));
        x_decoded(m, 1, :) = px_r;
        x_decoded(m, 2, :) = py_r;
        x_decoded(m, 3, :) = vx_r;
        x_decoded(m, 4, :) = vy_r;
        % Calculate RMSE
        RMSE_px(m, reach_num) = RMSE(px(reach_ind), px_r);
        RMSE_py(m, reach_num) = RMSE(py(reach_ind), py_r);
        RMSE_vx(m, reach_num) = RMSE(vx(reach_ind), vx_r);
        RMSE_vy(m, reach_num) = RMSE(vy(reach_ind), vy_r);
    end

    if plot_flag == 1,

        % Plot results
        fig = figure('NumberTitle','off','Units','normalized','Position',[1/8 1/8 3/4 3/4]);
        subplot(2, 3, [1 4]);
        hold on;

        % Dimensions for plotting task (from reach_animate.m)
        cntr = [-31   4];  % Center point
        xlim = [-46 -16];  % Plot bounds
        ylim = [-11 19];   % Plots bounds
        axis equal, set(gca,'XLim',xlim,'YLim',ylim);  % Axes
        targsz = 0.85; % target size (cm)
        cx = [cntr(1)-targsz cntr(1)+targsz cntr(1)+targsz cntr(1)-targsz];
        cy = [cntr(2)-targsz cntr(2)-targsz cntr(2)+targsz cntr(2)+targsz];
        ctar = patch(cx,cy,'k','EdgeColor','k','FaceColor','none');  % Center target
        %ptar = patch(cx,cy,'k','EdgeColor','k','FaceColor','none');  % Peripheral target
        ang = pi/4*tof(reach_start)-pi/8;
        target_dist = 10;
        %set(ptar,'XData',cx+target_dist*cos(ang),'YData',cy+10*sin(ang),'Visible','on');

        num_targets = 8;
        target_ang = pi*2/num_targets*(1:num_targets)-pi/num_targets;
        target_num = tof(reach_start);
        for j = 1 : num_targets,
        color = [0.7, 0.7, 0.7];
        if target_num == j,
            color = [0, 0, 0];
        end
            ang = target_ang(j);
            ptar = patch(cx + target_dist * cos(ang),cy + target_dist * sin(ang), 'k', 'EdgeColor', color, 'FaceColor', 'none');
        end
        xlabel('x (cm)');
        ylabel('y (cm)');

        % Plotting colors
        real_color = 'c';
        method_color = ['b'; 'r'; 'g'; 'k'];

        time_vec = time(reach_ind);

        % Plot 2D
        fig = plot(px(reach_ind), py(reach_ind), real_color);
        h = [fig];
        for m = 1 : M,
            fig = plot(squeeze(x_decoded(m, 1, :)), squeeze(x_decoded(m, 2, :)), method_color(m));
            h = [h, fig];
        end
        xlabel('x (cm)');
        ylabel('y (cm)');
        leg = legend(h, ['real'; method_names]);
        set(leg, 'Location', 'SouthOutside');
        title(['Reach ' int2str(reach_num)]);

        % Plot px
        subplot(2, 3, 2);
        hold on;
        plot(time_vec, px(reach_ind), real_color);
        for m = 1 : M,
            plot(time_vec, squeeze(x_decoded(m, 1, :)), method_color(m));
        end
        xlabel('time (sec)');
        ylabel('x position (cm)');

        % Plot py
        subplot(2, 3, 3);
        hold on;
        plot(time_vec, py(reach_ind), real_color);
        for m = 1 : M,
            plot(time_vec, squeeze(x_decoded(m, 2, :)), method_color(m));
        end
        xlabel('time (sec)');
        ylabel('y position (cm)');

        % Plot vx
        subplot(2, 3, 5);
        hold on;
        plot(time_vec, vx(reach_ind), real_color);
        for m = 1 : M,
            plot(time_vec, squeeze(x_decoded(m, 3, :)), method_color(m));
        end
        xlabel('time (sec)');
        ylabel('x velocity (cm/s)');

        % Plot vy
        subplot(2, 3, 6);
        hold on;
        plot(time_vec, py(reach_ind), real_color);
        for m = 1 : M,
            plot(time_vec, squeeze(x_decoded(m, 4, :)), method_color(m));
        end
        xlabel('time (sec)');
        ylabel('y velocity (cm/s)');

        %saveas(fig, ['/home/alex/Desktop/decodes/reach' int2str(reach_num) '.pdf'], 'pdf');

    end

end
