% Runs LFP analysis
% Plots rayleigh p, KS p, preferred phase, for low/high as a function of T
% methods: real, AT, aFRI
% Saves results in .mat files for use by (for instance) LFP_SNR
function [] = run_LFP_2(data_set, elec, compute)

% Refractory period
Delta = 0.0011;

if compute

    % Data segment
    [start_time, end_time] = get_total_time(data_set);
    start_time = start_time + 1;  % From reconstruct_and_save

    % Real spikes
    real_spk = real_spikes(data_set, elec, start_time, end_time, 1);
    %num_real_spk = length(real_spk)
    %real_start = min(real_spk)
    %real_end = max(real_spk)
    %real_dur = real_end-real_start
    name = struct();
    name.title = ['channel ' int2str(elec) ', real spikes'];
    name.data_set = data_set;
    name.T = 0;
    name.method_name = 'real';
    [plo_real phi_real pdlo_real pdhi_real distrlo_real distrhi_real datlo_real dathi_real cutoff pdglm_real] = nlx_spikephase_2(elec, real_spk, name);
    cutoff
    
    T_vec = [0.001:0.002:0.01 0.015:0.005:0.03 0.04:0.01:0.10];  % Sampling periods (from uni_params)
    plo_AT = zeros(size(T_vec));
    phi_AT = zeros(size(T_vec));
    plo_aFRI = zeros(size(T_vec));
    phi_aFRI = zeros(size(T_vec));
    plo_TD = zeros(size(T_vec));
    phi_TD = zeros(size(T_vec));
    pdlo_AT = zeros(size(T_vec));
    pdhi_AT = zeros(size(T_vec));
    pdlo_aFRI = zeros(size(T_vec));
    pdhi_aFRI = zeros(size(T_vec));
    pdlo_TD = zeros(size(T_vec));
    pdhi_TD = zeros(size(T_vec));
    distrlo_AT = zeros(length(T_vec), length(distrlo_real));
    distrhi_AT = zeros(length(T_vec), length(distrlo_real));
    distrlo_aFRI = zeros(length(T_vec), length(distrlo_real));
    distrhi_aFRI = zeros(length(T_vec), length(distrlo_real));
    distrlo_TD = zeros(length(T_vec), length(distrlo_real));
    distrhi_TD = zeros(length(T_vec), length(distrlo_real));
    datlo_AT = cell(1, length(T_vec));
    dathi_AT = cell(1, length(T_vec));
    datlo_aFRI = cell(1, length(T_vec));
    dathi_aFRI = cell(1, length(T_vec));
    datlo_TD = cell(1, length(T_vec));
    dathi_TD = cell(1, length(T_vec));
    pdglm_AT = zeros(size(T_vec));
    pdglm_aFRI = zeros(size(T_vec));
    pdglm_TD = zeros(size(T_vec));

    % AT
    method = 'AT';
    method_name = 'AT';

    for i = 1 : length(T_vec)
        fprintf(sprintf('AT:\n\telec: %d\n\ti: %d / %d\n', elec, i, length(T_vec)));
        paramID = i + 9;
        
        % TESTING: only do a particular T
        %if i ~= 5
        %    continue;
        %end

        filename = [root_dir() 'reconstruct/saved_spikes/spk_' method '_' int2str(data_set) '_' int2str(paramID) '.mat'];
        load(filename);  % Loads t, c
        t = t{elec}';
        c = c{elec}';
        ind = find(c > 0);
        t = t(ind);
        c = c(ind);
        %num_AT = length(t)
        %AT_start = min(t)
        %AT_end = max(t)
        %AT_dur = AT_end-AT_start
        
        % Refractory period Delta
        t = impose_delta(t, Delta);
        
        name = struct();
        name.title = ['channel ' int2str(elec) ', ' method_name ', T = ' num2str(T_vec(i))];
        name.T = T_vec(i);
        name.data_set = data_set;
        name.method_name = 'AT';
        [plo phi pdlo pdhi distrlo distrhi datlo dathi cutoff pd_glm] = nlx_spikephase_2(elec, t, name, cutoff);
        plo_AT(i) = plo;
        phi_AT(i) = phi;
        pdlo_AT(i) = pdlo;
        pdhi_AT(i) = pdhi;
        distrlo_AT(i,:) = distrlo;
        distrhi_AT(i,:) = distrhi;
        datlo_AT{i} = datlo;
        dathi_AT{i} = dathi;
        pdglm_AT(i) = pd_glm;
    end

    % aFRI
    method = 'analog';
    method_name = 'aFRI';

    for i = 1 : length(T_vec)
        fprintf(sprintf('analog:\n\telec: %d\n\ti: %d / %d\n', elec, i, length(T_vec)));
		paramID = i + 9;
        
        % TESTING: only do a particular T
        %if i ~= 5
        %    continue;
        %end

        filename = [root_dir() 'reconstruct/saved_spikes/spk_' method '_' int2str(data_set) '_' int2str(paramID) '.mat'];
        load(filename);  % Loads t, c
        t = t{elec}';
        c = c{elec}';
        ind = find(c > 0);
        t = t(ind);
        c = c(ind);
        %num_F = length(t)
        %F_start = min(t)
        %F_end = max(t)
        %F_dur = F_end-F_start
        
        % Refractory period Delta
        t = impose_delta(t, Delta);
        
        name.title = ['channel ' int2str(elec) ', ' method_name ', T = ' num2str(T_vec(i))];
        name.T = T_vec(i);
        name.data_set = data_set;
        name.method_name = 'gAT';
        [plo phi pdlo pdhi distrlo distrhi datlo dathi cutoff pd_glm] = nlx_spikephase_2(elec, t, name, cutoff);
        plo_aFRI(i) = plo;
        phi_aFRI(i) = phi;
        pdlo_aFRI(i) = pdlo;
        pdhi_aFRI(i) = pdhi;
        distrlo_aFRI(i,:) = distrlo;
        distrhi_aFRI(i,:) = distrhi;
        datlo_aFRI{i} = datlo;
        dathi_aFRI{i} = dathi;
        pdglm_aFRI(i) = pd_glm;
    end
    
    % TD
    method = 'twodelta';
    method_name = 'twodelta';

    for i = 1 : length(T_vec)
        fprintf(sprintf('twodelta:\n\telec: %d\n\ti: %d / %d\n', elec, i, length(T_vec)));
        paramID = i + 9;
        
        % TESTING: only do a particular T
        %if i ~= 5
        %    continue;
        %end

        filename = [root_dir() 'reconstruct/saved_spikes/spk_' method '_' int2str(data_set) '_' int2str(paramID) '.mat'];
        load(filename);  % Loads t, c
        t = t{elec}';
        c = c{elec}';
        ind = find(c > 0);
        t = t(ind);
        c = c(ind);
        %num_TD = length(t)
        %TD_start = min(t)
        %TD_end = max(t)
        %TD_dur = TD_end-TD_start
        
        % Refractory period Delta
        t = impose_delta(t, Delta);
        
        name = struct();
        name.title = ['channel ' int2str(elec) ', ' method_name ', T = ' num2str(T_vec(i))];
        name.T = T_vec(i);
        name.data_set = data_set;
        name.method_name = 'TD';
        [plo phi pdlo pdhi distrlo distrhi datlo dathi cutoff pd_glm] = nlx_spikephase_2(elec, t, name, cutoff);
        plo_TD(i) = plo;
        phi_TD(i) = phi;
        pdlo_TD(i) = pdlo;
        pdhi_TD(i) = pdhi;
        distrlo_TD(i,:) = distrlo;
        distrhi_TD(i,:) = distrhi;
        datlo_TD{i} = datlo;
        dathi_TD{i} = dathi;
        pdglm_TD(i) = pd_glm;
    end


    %TESTING
    %%{

    % KS test
    KSlo_AT = zeros(size(T_vec));
    KShi_AT = zeros(size(T_vec));
    KSlo_aFRI = zeros(size(T_vec));
    KShi_aFRI = zeros(size(T_vec));
    KSlo_TD = zeros(size(T_vec));
    KShi_TD = zeros(size(T_vec));
    for i = 1 : length(T_vec)
        [h p] = kstest2(datlo_AT{i}, datlo_real);
        KSlo_AT(i) = p;
        [h p] = kstest2(dathi_AT{i}, dathi_real);
        KShi_AT(i) = p;
        [h p] = kstest2(datlo_aFRI{i}, datlo_real);
        KSlo_aFRI(i) = p;
        [h p] = kstest2(dathi_aFRI{i}, dathi_real);
        KShi_aFRI(i) = p;
        [h p] = kstest2(datlo_TD{i}, datlo_real);
        KSlo_TD(i) = p;
        [h p] = kstest2(dathi_TD{i}, dathi_real);
        KShi_TD(i) = p;
    end

    save([root_dir() '../LFP/LFP_' int2str(data_set) '_' int2str(elec)]);
    %%}

else

    load([root_dir() '../LFP/LFP_' int2str(data_set) '_' int2str(elec)]);

end

% TESTING
%%{

% log base 10
%%{
plo_real = log(plo_real)/log(10);
phi_real = log(phi_real)/log(10);
plo_AT = log(plo_AT)/log(10);
phi_AT = log(phi_AT)/log(10);
plo_aFRI = log(plo_aFRI)/log(10);
phi_aFRI = log(phi_aFRI)/log(10);
plo_TD = log(plo_TD)/log(10);
phi_TD = log(phi_TD)/log(10);
%KSlo_AT = log(KSlo_AT)/log(10);
%KShi_AT = log(KShi_AT)/log(10);
%KSlo_aFRI = log(KSlo_aFRI)/log(10);
%KShi_aFRI = log(KShi_aFRI)/log(10);
%KSlo_TD = log(KSlo_TD)/log(10);
%KShi_TD = log(KShi_TD)/log(10);
%%}

% Plot

% rayleigh p-value
fig = figure;
hold on;
h1 = plot_errorbar(T_vec, plo_real*ones(size(T_vec)), 0, 0, lighten([0 0 0],1), 1);
h2 = plot_errorbar(T_vec, phi_real*ones(size(T_vec)), 0, 0, [0 0 0], 1);
h3 = plot_errorbar(T_vec, plo_AT, 0, 0, lighten([1 0 0],1), 1);
h4 = plot_errorbar(T_vec, phi_AT, 0, 0, [1 0 0], 1);
h5 = plot_errorbar(T_vec, plo_aFRI, 0, 0, lighten([0 0 1],1), 1);
h6 = plot_errorbar(T_vec, phi_aFRI, 0, 0, [0 0 1], 1);
h7 = plot_errorbar(T_vec, plo_TD, 0, 0, lighten([0 1 0],1), 1);
h8 = plot_errorbar(T_vec, phi_TD, 0, 0, [0 1 0], 1);
leg = legend([h1 h2 h3 h4 h5 h6 h7 h8], 'real (low amp)', 'real (hi amp)', 'AT (low amp)', 'AT (hi amp)', 'aFRI (low amp)', 'aFRI (hi amp)', 'twodelta (low amp)', 'twodelta (high amp)');
set(leg, 'Location', 'SouthEast');
xlabel('Sampling Period (T)');
ylabel('log p-Value (Rayleigh Test of Uniformity)');
title(['channel ' int2str(elec)]);
print(fig, [root_dir() '../paper2plots/LFP/summary/rayleigh_' int2str(elec)], '-dpdf', '-r0');

% KS p-value (compared to real)
fig = figure;
hold on;
h3 = plot_errorbar(T_vec, KSlo_AT, 0, 0, lighten([1 0 0],1), 1);
h4 = plot_errorbar(T_vec, KShi_AT, 0, 0, [1 0 0], 1);
h5 = plot_errorbar(T_vec, KSlo_aFRI, 0, 0, lighten([0 0 1],1), 1);
h6 = plot_errorbar(T_vec, KShi_aFRI, 0, 0, [0 0 1], 1);
h7 = plot_errorbar(T_vec, KSlo_TD, 0, 0, lighten([0 1 0],1), 1);
h8 = plot_errorbar(T_vec, KShi_TD, 0, 0, [0 1 0], 1);
leg = legend([h3 h4 h5 h6 h7 h8], 'AT (low amp)', 'AT (hi amp)', 'aFRI (low amp)', 'aFRI (hi amp)', 'twodelta (low amp)', 'twodelta (hi amp)');
set(leg, 'Location', 'SouthWest', 'FontSize', 28);
xlabel('Sampling Period (T)', 'FontSize', 28);
ylabel('p-Value (KS Test of Comparison)', 'FontSize', 28);
title(['channel ' int2str(elec)]);
print(fig, [root_dir() '../paper2plots/LFP/summary/KS_' int2str(elec)], '-dpdf', '-r0');

% preferred phase
fig = figure;
hold on;
%h1 = plot_errorbar(T_vec, pdlo_real*ones(size(T_vec)), 0, 0, lighten([0 0 0],1), 1);
h2 = plot_errorbar(T_vec, pdhi_real*ones(size(T_vec)), 0, 0, [0 0 0], 1);
%h3 = plot_errorbar(T_vec, pdlo_AT, 0, 0, lighten([1 0 0],1), 1);
h4 = plot_errorbar(T_vec, pdhi_AT, 0, 0, [1 0 0], 1);
%h5 = plot_errorbar(T_vec, pdlo_aFRI, 0, 0, lighten([0 0 1],1), 1);
h6 = plot_errorbar(T_vec, pdhi_aFRI, 0, 0, [0 0 1], 1);
%h7 = plot_errorbar(T_vec, pdlo_TD, 0, 0, lighten([0 1 0],1), 1);
h8 = plot_errorbar(T_vec, pdhi_TD, 0, 0, [0 1 0], 1);
leg = legend([h2 h4 h6 h8], 'True', 'AT', 'gAT', 'twodelta');
set(leg, 'Location', 'NorthWest', 'FontSize', 28);
xlabel('Sampling Period (T)', 'FontSize', 28);
ylabel('Preferred Phase', 'FontSize', 28);
title(['channel ' int2str(elec)]);
print(fig, [root_dir() '../paper2plots/LFP/summary/pref_' int2str(elec)], '-dpdf', '-r0');

%%}
