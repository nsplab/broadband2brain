function [RMSE_px, RMSE_py, RMSE_vx, RMSE_vy, RMSE_p, RMSE_v RMSE_still] = pp_filter_rand(data_set, start_time, end_time, num_reaches, spike_times, plot_flag, method_names, delay, b)
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
% RMSE_still(reach_num) = root mean squared error in position for no motion

% Turn on/off trajectory analysis
traj_test = 0;

% Parameters
dt = get_dt_info();
A = [1, 0, dt, 0; 0, 1, 0, dt; 0, 0, 1, 0; 0, 0, 0, 1];
W_0 = eye(4) * 1;  % Covariance for first step of filter
Q = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1] * 5e-1;  % TODO: Modify random walk update covariance (for decoding) here

% Data
[d, filename] = data_file(data_set, 0, 0); % 'matthew20080721';
load(filename, '-mat');  % Loads time, px, py, vx, vy, etc. (data regarding monkey's movement)

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

if traj_test == 1

    % Plot trajectories
    fig_vec = zeros(8, 1);
    p_fig_vec = zeros(8, 1);
    vx_fig_vec = zeros(8, 1);
    vy_fig_vec = zeros(8, 1);
    for i = 1 : 8
        fig_vec(i) = figure;
        p_fig_vec(i) = subplot(1, 3, 1);
        hold on;
        xlim = [-46 -16];  % Plot bounds
        ylim = [-11 19];   % Plots bounds
        axis equal, set(gca,'XLim',xlim,'YLim',ylim);  % Axes
        title({['target ' int2str(i)], 'position'});
        vx_fig_vec(i) = subplot(1, 3, 2);
        hold on;
        title('vx');
        vy_fig_vec(i) = subplot(1, 3, 3);
        hold on;
        title('vy');
    end

    % Vectors for average reach time
    time_sum = zeros(8, 1);
    num_reach_vec = zeros(8, 1);

    target_nums = zeros(num_reaches, 1);

end

% Convert spike times to spike vector
% spike_vec(method_num, neuron_num, timestep) = 0 or 1 (indicating presence or absence of spike)
for m = 1 : M,
    spike_vec(m, :, :) = spike_times_to_vec(time, spike_times(m, :));  % TODO: allow multiple spikes per bin?
end

RMSE_still = [];

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

        reach_start_inds(reach_num) = reach_start;

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
        RMSE_p(m, reach_num) = RMS(sqrt((px(reach_ind)-px_r).^2 + (py(reach_ind)-py_r).^2));
        RMSE_v(m, reach_num) = RMS(sqrt((vx(reach_ind)-vx_r).^2 + (vy(reach_ind)-vy_r).^2));
    end
    
    RMSE_still = [RMSE_still RMS(sqrt((px(reach_ind)-px(reach_ind(1))).^2 + (py(reach_ind)-py(reach_ind(1))).^2))];
    
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
        ctar = patch(cx,cy,'k','EdgeColor','k','FaceColor','none', 'LineWidth', 2);  % Center target
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
            ptar = patch(cx + target_dist * cos(ang),cy + target_dist * sin(ang), 'k', 'EdgeColor', color, 'FaceColor', 'none', 'LineWidth', 2);
        end
        xlabel('x (cm)', 'FontSize', 12);
        ylabel('y (cm)', 'FontSize', 12);

        % Plotting colors
        real_color = [0 0 0]; %'c';
        method_color = [1 0 0 ; 0 0 1];%['b'; 'r'; 'g'; 'k'];

        time_vec = time(reach_ind);
        time_vec = time_vec - time_vec(1);  % start reach at time = 0

        % Plot 2D
        fig = plot_line(px(reach_ind), py(reach_ind), real_color);
        h = [fig];
        for m = 1 : M,
            fig = plot_line(squeeze(x_decoded(m, 1, :)), squeeze(x_decoded(m, 2, :)), method_color(m,:));
            h = [h, fig];
        end
        box on;
        xlabel('x (cm)', 'FontSize', 12);
        ylabel('y (cm)', 'FontSize', 12);
        %leg = legend(h, ['true trajectory'; method_names]);
        leg = legend(h, 'true trajectory', 'reconstructed (real spikes)', 'reconstructed (gibbs)');
        set(leg, 'Location', 'SouthOutside', 'FontSize', 12);
        %title({['Reach ' int2str(reach_num)], ['RMSE p real: ' num2str(RMSE_p(1, reach_num)) ', reconstructed: ' num2str(RMSE_p(2, reach_num))]});
        title(['RMSE real: ' num2str(RMSE_p(1, reach_num)) ', gibbs: ' num2str(RMSE_p(2, reach_num))], 'FontSize', 12);
        axis square;
        
        % Plot px
        subplot(2, 3, 2);
        hold on;
        plot_line(time_vec, px(reach_ind), real_color);
        for m = 1 : M,
            plot_line(time_vec, squeeze(x_decoded(m, 1, :)), method_color(m,:));
        end
        xlabel('time (sec)', 'FontSize', 12);
        ylabel('x position (cm)', 'FontSize', 12);
        axis square;

        % Plot py
        subplot(2, 3, 3);
        hold on;
        plot_line(time_vec, py(reach_ind), real_color);
        for m = 1 : M,
            plot_line(time_vec, squeeze(x_decoded(m, 2, :)), method_color(m,:));
        end
        xlabel('time (sec)', 'FontSize', 12);
        ylabel('y position (cm)', 'FontSize', 12);
        axis square

        % Plot vx
        subplot(2, 3, 5);
        hold on;
        plot_line(time_vec, vx(reach_ind), real_color);
        for m = 1 : M,
            plot_line(time_vec, squeeze(x_decoded(m, 3, :)), method_color(m,:));
        end
        xlabel('time (sec)', 'FontSize', 12);
        ylabel('x velocity (cm/s)', 'FontSize', 12);
        axis square;

        % Plot vy
        subplot(2, 3, 6);
        hold on;
        plot_line(time_vec, vy(reach_ind), real_color);
        for m = 1 : M,
            plot_line(time_vec, squeeze(x_decoded(m, 4, :)), method_color(m,:));
        end
        xlabel('time (sec)', 'FontSize', 12);
        ylabel('y velocity (cm/s)', 'FontSize', 12);
        axis square;

        %saveas(fig, ['/home/alex/Desktop/decodes/reach' int2str(reach_num) '.pdf'], 'pdf');

        %close;  % Close figure

        if traj_test == 1

            % Plot trajectories
            plot(p_fig_vec(target_num), px(reach_ind), py(reach_ind));
            plot(vx_fig_vec(target_num), time_vec-time_vec(1), vx(reach_ind));
            plot(vy_fig_vec(target_num), time_vec-time_vec(1), vy(reach_ind));
            plot(vx_fig_vec(target_num), time_vec(end)-time_vec(1), vx(reach_ind(end)), 'd');
            plot(vy_fig_vec(target_num), time_vec(end)-time_vec(1), vy(reach_ind(end)), 'd');
            time_sum(target_num) = time_sum(target_num) + time_vec(end) - time_vec(1);
            num_reach_vec(target_num) = num_reach_vec(target_num) + 1;
            target_nums(reach_num) = target_num;

        end

    end

end

if traj_test == 1

    % Calc average reach times
    av_reach_time = time_sum ./ num_reach_vec;
    av_reach_ind = av_reach_time ./ get_dt_info();

    for i = 1 : 8

        count = 0;
        av_x = zeros(av_reach_ind(i), 1);
        av_y = zeros(av_reach_ind(i), 1);
        av_vx = zeros(av_reach_ind(i), 1);
        av_vy = zeros(av_reach_ind(i), 1);

        for reach_num = 1 : num_reaches,

            reach_start = reach_start_inds(reach_num);
            reach_ind = reach_start : reach_start + av_reach_ind(i) - 1;

            % Plot average trajectory
            if target_nums(reach_num) == i
                count = count + 1;
                av_x = av_x + px(reach_ind);
                av_y = av_y + py(reach_ind);
                av_vx = av_vx + vx(reach_ind);
                av_vy = av_vy + vy(reach_ind);
            end
        end

        av_x = av_x / count;
        av_y = av_y / count;
        av_vx = av_vx / count;
        av_vy = av_vy / count;

        time_vec = linspace(0, av_reach_time(i), length(reach_ind));
        plot(p_fig_vec(i), av_x, av_y, 'r');
        plot(vx_fig_vec(i), time_vec, av_vx, 'r');
        plot(vy_fig_vec(i), time_vec, av_vy, 'r');

    end

    % Save figures
    for i = 1 : 8
        saveas(fig_vec(i), ['/home/alex/Desktop/test_traj/target' int2str(i) '.pdf'], 'pdf');
    end
end

%close all;

end
