function thres = choose_threshold(data_set, elec, id, load_flag, t, c, t_real, tolerance, res)
% Chooses the optimal threshold for pruning spikes given by t, c to match t_real with fewest errors (false positives/negatives)
% Saves threshold to file and can load it later

% Inputs:
% data_set
% elec - single electrode
% id - identifier string, allows for different sets of parameters etc. (name of object saved to file)
% load_flag - 0 or 1, if threshold should be loaded (if threshold is loaded, other parameters not necessary)
% t - vector of spike times
% c - vector of spike amplitudes
% t_real - real spike times
% tolerance - see analyze/compare_spikes.m
% res - resolution of different threshold values tried

filename = [root_dir() 'estimate/thresholds.mat'];
S = load(filename);

if load_flag == 1
    thres = eval(['S.' id '(data_set, elec);']);
else

    % Compute threshold
    %thres_to_try = 0 : res : max(c)+0.1;
    thres_to_try = linspace(0, max(c), 100);  % TESTING: ignore res, fixed number of bins
    
    fpos = zeros(size(thres_to_try));
    fneg = zeros(size(thres_to_try));

    for i = 1 : length(thres_to_try)
        thres = thres_to_try(i);
        ind = find(c >= thres);
        [tp, fp, fn] = compare_spikes(t(ind), t_real, tolerance);
        fpos(i) = length(fp);
        fneg(i) = length(fn);
    end

    [val, in] = min(fpos + fneg);
    thres = thres_to_try(in);

    % TESTING: plot spikes
    %{
    figure;
    hold on;
    stem(t_real, ones(size(t_real))*max(c), 'c');
    stem(t, c, 'b');    
    plot([min(t), max(t)], [thres, thres], 'g--');
    legend('real', 'FRI', 'threshold');
    title(['computing best threshold, electrode ' int2str(elec)]);
    %}

    % Save threshold
    eval(['S.' id '(data_set, elec) = thres; ' id ' = S. ' id ';']);
    save(filename, id, '-append');

end

end
