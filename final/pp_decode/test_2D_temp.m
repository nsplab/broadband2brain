% Creates 2D plots of real vs reconstructed movement using real data and random walk prior

clear all;
close all;

% Parameters
% Start/end times (seconds) for training data and usable data
training_start = 77.554890000000000; %400;  % NOTE: training data should come before data
training_end = 1.8351e+03; %1.830306457500000e+03; %2000;
used_start = 77.554890000000000;
used_end = 1.8351e+03 + 100;
num_reaches = 15;
% Decoding
W_0 = eye(4) * 1;  % Covariance for first step of filter
Q = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1] * 1;  % Random walk update covariance (for decoding)   TODO
% Data
filename = 'matthew20080721';
load([filename '.mat'],'-mat');
delay = get_delay();  % Found using delay_est.m run on time

% Apply delay
for i = 1 : length(spk),
    spk(i).data = spk(i).data + delay(i);
end

% TESTING: only use some neurons    TODO
neurons_to_use = [6, 7, 10, 9, 11];
%spk = spk(neurons_to_use);  % Comment this line to use all neurons

C = length(spk);  % Number of neurons

dt = (time(length(time)) - time(1)) / (length(time) - 1);
used_ind = find(time >= used_start & time < used_end);  % Assumes training data comes before data

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

% Indices of training data
training_ind = find(time >= training_start & time < training_end);

% Reach times   TODO: selection process for start/end of reach can be modified here
cio(training_ind) = 0;  % No reaches in training data
tio(training_ind) = 0;  % No reaches in training data
%start_ind = find(cio & [diff(cof); 0] == -1);  % Center removed but hand still inside it
start_ind = find(diff(cio) == -1) - 70;  % Leave center     TODO: correct offset ~ -70, originally used -20
%end_ind = find([diff(tof); 0] < 0);  % When target disappears, acknowledging arrival
end_ind = find(diff(tio) == 1);  % As soon as hand gets to target

% Get spike vector
spike_vec = spike_times_to_vec(time, spk);

% Estimate b parameters
[b_0, b_1, b_2] = estimate_b(C, dt, vx(training_ind), vy(training_ind), spike_vec(:, training_ind));

for reach_num = 1 : num_reaches,

    reach_num

    % Indices of this reach
    reach_start = start_ind(reach_num);
    ends = find(end_ind > reach_start);
    reach_end = end_ind(ends(1));
    end_ind(ends(1)) = 0;  % Don't use same end again
    starts = find(start_ind < reach_end);
    reach_start = start_ind(starts(length(starts))) - 1;
    data_ind = reach_start : reach_end;

    time(reach_start)
    time(reach_end)

    % Num timesteps
    K = length(time(data_ind)) - 1;

    % Correct starting state  TODO
    %start = [px(data_ind(1)); py(data_ind(1)); vx(data_ind(1)); vy(data_ind(1))];
    start = [px(data_ind(1)); py(data_ind(1)); 0; 0];  % Start at a standstill at correct position

    % Point process filter (random walk prior)
    A = [1, 0, dt, 0; 0, 1, 0, dt; 0, 0, 1, 0; 0, 0, 0, 1];
    [px_r, py_r, vx_r, vy_r] = rand_walk_filter(K, C, dt, start, W_0, A, Q, spike_vec(:, data_ind), b_0, b_1, b_2);

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
    ptar = patch(cx,cy,'k','EdgeColor','k','FaceColor','none');  % Peripheral target
    ang = pi/4*tof(reach_start)-pi/8;
    set(ptar,'XData',cx+10*cos(ang),'YData',cy+10*sin(ang),'Visible','on');
    xlabel('x (cm)');
    ylabel('y (cm)');

    % Plot 2D
    p1 = plot(px(data_ind), py(data_ind), 'c');
    p2 = plot(px_r, py_r, 'b');
    legend([p1, p2], 'real', 'reconstructed');
    title(['Reach ' int2str(reach_num)]);

    % Plot px
    subplot(2, 3, 2);
    hold on;
    plot(time(data_ind), px(data_ind), 'c');
    plot(time(data_ind), px_r, 'b');
    xlabel('time (sec)');
    ylabel('x position (cm)');

    % Plot py
    subplot(2, 3, 3);
    hold on;
    plot(time(data_ind), py(data_ind), 'c');
    plot(time(data_ind), py_r, 'b');
    xlabel('time (sec)');
    ylabel('y position (cm)');

    % Plot vx
    subplot(2, 3, 5);
    hold on;
    plot(time(data_ind), vx(data_ind), 'c');
    plot(time(data_ind), vx_r, 'b');
    xlabel('time (sec)');
    ylabel('x velocity (cm/s)');

    % Plot vy
    subplot(2, 3, 6);
    hold on;
    plot(time(data_ind), py(data_ind), 'c');
    plot(time(data_ind), py_r, 'b');
    xlabel('time (sec)');
    ylabel('y velocity (cm/s)');

    saveas(fig, ['/home/alex/Documents/SURF/Reach plots/2D/reach' int2str(reach_num) '.pdf'], 'pdf');

end
