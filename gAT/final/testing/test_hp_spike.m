% Look at impulse response of various initial hp filters

dt = get_dt();

time = 0 : dt : 0.1;
dat = zeros(size(time));

% Add spike
spk_ind = round(length(dat)/3);
for i = spk_ind : spk_ind %+ round(0.00025/dt)
    dat(i) = 1;
end

[A, B] = butter_hp(1000, 2);
%[A, B] = butter_bp(300, 6000, 2);
args.A_filter = A;
args.B_filter = B;
hp_handle = @filter_generic;  % Causal
%hp_handle = @filtfilt_generic;  % Noncausal

% Apply highpass
dat_hp = hp_handle(dat, args);

figure;
hold on;
plot(time, dat, 'c');
plot(time, dat_hp, 'b');

dat_integral = sum(dat)
dat_hp_integral = sum(abs(dat_hp))
