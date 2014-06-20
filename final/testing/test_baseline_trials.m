function [base_end, ff_start, ff_end, wash_start, wash_end] = test_baseline_trials(data_set)
%monkey = 'nemo'; % monkey ('nemo' or 'matthew')
%session = '20080721'; % session (yyyymmdd)
%load([root_dir() 'data/' monkey session '.mat']);
%data_set = 2;
[d, file] = data_file(data_set);
load(file);
indres = find(tio==1 & [0;diff(tof)<=-1]); % indices of rewarded delivery
%tbend = time(indres(160)); % end of the baseline trials (s)

total_trials = length(indres);

base_end_ind = 160;
ff_end_ind = 320;
wash_end_ind = 480;

if total_trials < 480
    base_end_ind = round(total_trials / 3);
    ff_end_ind = round(total_trials * 2/3);
    wash_end_ind = total_trials;
end

base_end = time(indres(base_end_ind));
ff_start = base_end;
ff_end = time(indres(ff_end_ind));
wash_start = ff_end;
wash_end = time(indres(wash_end_ind));

%ff_base = sum(ff(find(time < base_end)))
%ff_ff = mean(ff(find(time >= ff_start & time < ff_end)))
%ff_wash = sum(ff(find(time >= wash_start & time < wash_end)))