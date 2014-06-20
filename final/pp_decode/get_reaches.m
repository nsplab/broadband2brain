function [start_times, end_times] = get_reaches(data_set, start_time, end_time)
% Returns start and end times for center-out reaching movements in a given time period

[d, info] = data_file(data_set, 0, 0);
load(info);

used_ind = find(time >= start_time & time < end_time);  % Assumes training data comes before data

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

%start_ind = find(cio & [diff(cof); 0] == -1);  % Center removed but hand still inside it
start_ind = find(diff(cio) == -1) - 70;  % Leave center     TODO: correct offset ~ -70, originally used -20
%end_ind = find([diff(tof); 0] < 0);  % When target disappears, acknowledging arrival
end_ind = find(diff(tio) == 1);  % As soon as hand gets to target

% Potential start and end times
start_times = time(start_ind);
end_times = time(end_ind);

% Select matched start and end times
start_ind = zeros(size(start_times));
end_ind = zeros(size(end_times));

i = length(end_ind);  % Index in end_times
j = length(start_ind);  % Index in start_times
while i >= 1
    while start_times(j) >= end_times(i) && j >= 0
        j = j - 1;
    end
    if j == 0
        break;
    end
    end_ind(i) = 1;
    start_ind(j) = 1;
    i = i - 1;
end

start_times = start_times(find(start_ind));
end_times = end_times(find(end_ind));

end
