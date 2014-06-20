% This script is used to analyze twodelta for figure 2.

clear;
close all;

load([root_dir() 'paper2/mat/fig2']);


default_colors = {[0 0 1], [0 0.5 0], [1 0 0], [0 0.75 0.75], [0.75 0 0.75], [0.75 0.75 0], [0.25 0.25 0.25]};




% Computing the confusion matrices
assert(all(size(num_rec_spikes_Flist) == size(num_real_spikeslist)), 'Different number of spikes...');
assert(all(size(num_rec_spikes_TDlist) == size(num_real_spikeslist)), 'Different number of spikes...');

countTD = cell(size(unique(Tlist), 1), 1);
for index = 1:size(unique(Tlist), 1)
    times = unique(Tlist);
    time = times(index);
    fprintf('time: %f\n', time);
    countTD{index} = zeros(max(num_rec_spikes_TDlist) + 1, max(num_real_spikeslist) + 1);
    for i = 1:size(num_real_spikeslist, 1)
        if (Tlist(i) == time)
            countTD{index}(num_rec_spikes_TDlist(i) + 1, num_real_spikeslist(i) + 1) = ....
              countTD{index}(num_rec_spikes_TDlist(i) + 1, num_real_spikeslist(i) + 1) + 1;
        end
    end
    countTD{index}
end

fprintf('All:\n');
%countF = zeros(max(num_rec_spikes_Flist) + 1, max(num_real_spikeslist) + 1);
%countTD = zeros(max(num_rec_spikes_TDlist) + 1, max(num_real_spikeslist) + 1);
%for i = 1:size(num_real_spikeslist, 1)
%    countF(num_rec_spikes_Flist(i) + 1, num_real_spikeslist(i) + 1) = ....
%      countF(num_rec_spikes_Flist(i) + 1, num_real_spikeslist(i) + 1) + 1;
%    countTD(num_rec_spikes_TDlist(i) + 1, num_real_spikeslist(i) + 1) = ....
%      countTD(num_rec_spikes_TDlist(i) + 1, num_real_spikeslist(i) + 1) + 1;
%end
%countF
%countTD
countTD_all = zeros(size(countTD{1}));
for i = 1:size(countTD, 1)
    countTD_all = countTD_all + countTD{i};
end
countTD_all

num_true = cell(max(num_real_spikeslist) + 1, 1);
for i = 1:(max(num_real_spikeslist) + 1)
    num_true{i} = zeros(1, size(unique(Tlist), 1));
    for j = 1:size(unique(Tlist), 1)
        num_true{i}(j) = sum(countTD{j}(:, i));
    end
end

confusion = cell(max(num_rec_spikes_TDlist) + 1, max(num_real_spikeslist) + 1);
for i = 1:(max(num_rec_spikes_TDlist) + 1)
    for j = 1:(max(num_real_spikeslist) + 1)
        confusion{i, j} = zeros(1, size(unique(Tlist), 1));
        for k = 1:size(unique(Tlist), 1)
            confusion{i, j}(k) = countTD{k}(i, j);
        end
    end
end

figure;
plot(1000 * times, confusion{1, 2}, 1000 * times, confusion{2, 2}, 1000 * times, confusion{3, 2});
legend('1 -> 0', '1 -> 1', '1 -> 2')
xlabel('Interval (ms)');
ylabel('Count');
figure;
plot(1000 * times, confusion{1, 2}./num_true{2}, 1000 * times, confusion{2, 2}./num_true{2}, 1000 * times, confusion{3, 2}./num_true{2});
legend('1 -> 0', '1 -> 1', '1 -> 2')
xlabel('Interval (ms)');
ylabel('Fraction');

figure;
plot(1000 * times, confusion{1, 3}, 1000 * times, confusion{2, 3}, 1000 * times, confusion{3, 3});
legend('2 -> 0', '2 -> 1', '2 -> 2')
xlabel('Interval (ms)');
ylabel('Count');
figure;
plot(1000 * times, confusion{1, 3}./num_true{3}, 1000 * times, confusion{2, 3}./num_true{3}, 1000 * times, confusion{3, 3}./num_true{3});
legend('2 -> 0', '2 -> 1', '2 -> 2')
xlabel('Interval (ms)');
ylabel('Fraction');

figure;
plot(1000 * times, num_true{2} - confusion{2, 2}, 1000 * times, num_true{3} - confusion{3, 3});
legend('1', '2');
xlabel('Interval (ms)');
ylabel('Number of Errors');

figure;
plot(1000 * times, confusion{3, 2}, 1000 * times, confusion{2, 3});
legend('1 -> 2', '2 -> 1');
xlabel('Interval (ms)');
ylabel('Count');

% Plot the number of n-spike intervals as a function of interval width
max_spikes = max(num_real_spikeslist);
times = unique(Tlist);
num_intervals = size(times, 1);
prob = zeros(max_spikes + 1, num_intervals);
for i = 1:(max_spikes+1)
    for j = 1:num_intervals
        prob(i, j) = sum(num_real_spikeslist == (i-1) & Tlist == times(j)) / sum(Tlist == times(j));
    end
end
f = figure;
hold on;
h = [];
for i = 1:(max_spikes+1)
    linestyle = '-';
    if (i > size(default_colors, 2))
        linestyle = '--';
    end
    h(i) = plot(1000 * times, prob(i, :)', 'Color', lighten(default_colors{mod(i-1, size(default_colors, 2))+1}), 'LineWidth', 3, 'LineStyle', linestyle);
    plot(1000 * times, prob(i, :)', 'o', 'Color', default_colors{mod(i-1, size(default_colors, 2))+1}, 'LineWidth', 3);
end
xlabel('Interval Length (ms)', 'FontSize', 28);
%ylabel('Fraction', 'FontSize', 28);
ylabel('Fraction of Intervals', 'FontSize', 28);
title('Fraction of Intervals with n Spikes', 'FontSize', 28);
size(h)
legend(h, '0', '1', '2', '3', '4', '5', '6', '7', '8', '9');
set(gca, 'FontSize', 24);
saveas(f, [root_dir() '../paper2plots/final/analyze/n_spike_interval'], 'epsc');
saveas(f, [root_dir() '../paper2plots/final/analyze/n_spike_interval'], 'fig');

f = figure;
hold on;
h = [];
cum_prob = cumsum(prob, 1);
for i = 1:(max_spikes+1)
    linestyle = '-';
    if (i > size(default_colors, 2))
        linestyle = '--';
    end
    h(i) = plot(1000 * times, cum_prob(i, :)', 'Color', lighten(default_colors{mod(i-1, size(default_colors, 2))+1}), 'LineWidth', 3, 'LineStyle', linestyle);
    plot(1000 * times, cum_prob(i, :)', 'o', 'Color', default_colors{mod(i-1, size(default_colors, 2))+1}, 'LineWidth', 3);
end
xlabel('Interval Length (ms)', 'FontSize', 28);
ylabel('Fraction', 'FontSize', 28);
title('Fraction of Intervals with \geq n Spikes', 'FontSize', 28);
legend(h, '0', '1', '2', '3', '4', '5', '6', '7', '8', '9');
set(gca, 'FontSize', 24);
saveas(f, [root_dir() '../paper2plots/final/analyze/n_spike_interval_cumulative'], 'epsc');
saveas(f, [root_dir() '../paper2plots/final/analyze/n_spike_interval_cumulative'], 'fig');



% Attempting to find a good threshold between 1 / 2 spikes.
%   - Real 0-spike intervals count as 1 spike
%   - > 2 real spikes will count as 2 spikes
%   - Intervals with 0 predicted spikes are ignored

% Remove intervals with first successive integral = 0 (zero predicted spikes).
hasSpikes = (ylist(:, 1) ~= 0);
Tlist = Tlist(hasSpikes, :);
ylist = ylist(hasSpikes, :);
yguess3list = yguess3list(hasSpikes, :);
yguess4list = yguess4list(hasSpikes, :);
num_rec_spikes_TDlist = num_rec_spikes_TDlist(hasSpikes, :);
num_real_spikeslist = num_real_spikeslist(hasSpikes, :);
num_rec_spikes_Flist = num_rec_spikes_Flist(hasSpikes, :);

% NOTE: It may make sense to move this *after* the linear regression is run.
%       Doing so would cause the linear regression to be heavily pulled by >2
%       spikes / 0 spikes, which might make sense.
%       Empirically, this seems to only sometimes improve the results.
% Less than 1 spikes (0 spikes) can only be treated as 1 spike
num_real_spikeslist(num_real_spikeslist < 1) = 1;
% Greater than 2 spikes can only be treated as 2 spike
num_real_spikeslist(num_real_spikeslist > 2) = 2;

% Compute the difference
diff3list = ylist(:, 3)-yguess3list;
diff4list = ylist(:, 4)-yguess4list;

% Compute the relative error
ratio3list = (abs(diff3list) ./ ylist(:, 3));
ratio4list = (abs(diff4list) ./ ylist(:, 4));

% Squared error (not used at the moment)
se3list = (diff3list.^2 ./ ylist(:, 3));
se4list = (diff4list.^2 ./ ylist(:, 4));

numIntervals = sum(hasSpikes);
num1spike = sum(num_real_spikeslist == 1);
num2spike = sum(num_real_spikeslist == 2);
prob1 = num1spike / numIntervals;
prob2 = num2spike / numIntervals;

% Equivalent to original twodelta implementation
%pred_num_spike(:) = 1;
%pred_num_spike(ratio3list > 0.09 & ratio4list > 0.24) = 2;

% Generate matrix
%   - rows: values for different intervals
%   - columns: different parameter to be used

%% Optimized version of original twodelta implementation
%%G = [ratio3list ratio4list];
%% Much more complex fit
%G = [Tlist (T.*ratio3list) (T.*ratio4list) (ratio3list.*ratio4list) ones(sum(hasSpikes), 1)];
%
%% Run linear regression to get approximate solution
%%m = (G' * G) \ G' * num_real_spikeslist;
%%pred_num_spike = G * m;
%%pred_num_spike(pred_num_spike < 1.5) = 1;
%%pred_num_spike(pred_num_spike >= 1.5) = 2;
%
%% Weighed linear regression
%G(num_real_spikeslist == 1, :) = G(num_real_spikeslist == 1, :) * sqrt(sqrt(prob1 * prob2) / prob1);
%G(num_real_spikeslist == 2, :) = G(num_real_spikeslist == 2, :) * sqrt(sqrt(prob1 * prob2) / prob2);
%desired = num_real_spikeslist;
%desired(num_real_spikeslist == 1) = -sqrt(sqrt(prob1 * prob2) / prob1);
%desired(num_real_spikeslist == 2) = sqrt(sqrt(prob1 * prob2) / prob2);
%an
%pred_num_spike(pred_num_spike >= 0) = 2;
%pred_num_spike(pred_num_spike < 0) = 1;
%
%solution = m;
%solution
%
%% Absolute correctness
%%best = sum(pred_num_spike == num_real_spikeslist);
%%total = sum(hasSpikes);
%%fprintf(sprintf('Accuracy: %d/%d = %f\n', best, total, best/total));
%
%% Weighted correctness (1 spike vs. 2 spike equalized)
%p1correct = sum(num_real_spikeslist == 1 & pred_num_spike == num_real_spikeslist) / num1spike;
%p2correct = sum(num_real_spikeslist == 2 & pred_num_spike == num_real_spikeslist) / num2spike;
%best = (p1correct + p2correct) / 2;
%fprintf(sprintf('Accuracy: %f\n\t1: %f\n\t2: %f\n', best, p1correct, p2correct));
%
%range = -1:0.5:1;
%[i1 i2 i3 i4 i5 i6] = ndgrid(range);
%increment = [i1(:) i2(:) i3(:) i4(:) i5(:)];
%for i = 1:size(increment, 1)
%    if (mod(i, 5000) == 0)
%        fprintf(sprintf('%d/%d\n = %f\n', i, size(increment, 1), i / size(increment, 1)));
%    end
%
%    m1 = m + increment(i, :)';
%
%    % Linear Regression
%    %pred_num_spike = G * m1;
%    %pred_num_spike(pred_num_spike < 1.5) = 1;
%    %pred_num_spike(pred_num_spike >= 1.5) = 2;
%
%    % Weighted linear regression
%    pred_num_spike = G * m1;
%    pred_num_spike(pred_num_spike >= 0) = 2;
%    pred_num_spike(pred_num_spike < 0) = 1;
%    
%    % Absolute correctness
%    %correct = sum(pred_num_spike == num_real_spikeslist);
%
%    % Weighted correctness (1 spike vs. 2 spike equalized)
%    p1correct = sum(num_real_spikeslist == 1 & pred_num_spike == num_real_spikeslist) / num1spike;
%    p2correct = sum(num_real_spikeslist == 2 & pred_num_spike == num_real_spikeslist) / num2spike;
%    correct = (p1correct + p2correct) / 2;
%
%    if (correct > best)
%        best = correct;
%        solution = m1;
%        solution
%
%        % Absolute correctness
%        %fprintf(sprintf('Accuracy: %d/%d = %f\n', correct, total, correct/total));
%
%        % Weighted correctness (1 spike vs. 2 spike equalized)
%        fprintf(sprintf('Accuracy: %f\n\t1: %f\n\t2: %f\n', correct, p1correct, p2correct));
%    end
%end







% Selecting the best threshold
%total = sum(hasSpikes);
%
%%G = [ratio3list ratio4list (Tlist.*ratio3list) (Tlist.*ratio4list) (ratio3list.*ratio4list)];
%G = [(Tlist.*ratio4list) (Tlist.^(1.2).*ratio4list) (Tlist.^(1.4).*ratio4list) (Tlist.^(1.6).*ratio4list) (Tlist.^(1.8).*ratio4list) (Tlist.^2.*ratio4list)];
%weight = ones(size(Tlist));
%for i = sort(unique(Tlist))'
%    weight(Tlist == i) = 1 / sum(Tlist == i);
%end
%weight = weight / mean(weight);
%weight = weight * 0;
%weight(Tlist == 0.1) = 1;
%num_tests = size(G, 2);
%threshold = size(G);
%for i = 1:size(G, 2)
%    best = -1;
%    for j = linspace(min(G(:, i)), max(G(:, i)), 100)
%        pred(G(:, i) <= j, :) = 1;
%        pred(G(:, i) > j, :) = 2;
%        correct = (pred == num_real_spikeslist)' * weight;
%        if (correct > best)
%            best = correct;
%            threshold(i) = j;
%        end
%    end
%
%    vote(:, i) = (G(:, i) > threshold(i));
%    fprintf(sprintf('Accuracy %d: %d/%d = %f\n', i, best, total, best/total));
%end
%
%votes = sum(vote, 2);
%pred(votes < num_tests / 2) = 1;
%pred(votes >= num_tests / 2) = 2;
%%pred = (ratio3list > 0.09 & ratio4list > 0.24) + 1;
%correct = (pred == num_real_spikeslist)' * weight;
%
%minimum = min(G);
%minimum
%maximum = max(G);
%maximum
%threshold
%fprintf(sprintf('Accuracy: %d/%d = %f\n', correct, total, correct/total));







%% Redefining without absolute value
%ratio3list = (diff3list ./ ylist(:, 3));
%ratio4list = (abs(diff4list) ./ ylist(:, 4));
%
%
%
%total = sum(hasSpikes);
%
%G = [ratio3list ratio4list (Tlist.*ratio3list) (Tlist.*ratio4list) (ratio3list.*ratio4list)];
%G = [(Tlist.*ratio4list) (Tlist.^(1.2).*ratio4list) (Tlist.^(1.4).*ratio4list) (Tlist.^(1.6).*ratio4list) (Tlist.^(1.8).*ratio4list) (Tlist.^2.*ratio4list)];
%weight = ones(size(Tlist));
%%for i = sort(unique(Tlist))'
%%    weight(Tlist == i) = 1 / sum(Tlist == i);
%%end
%weight = weight / mean(weight);
%num_tests = size(G, 2);
%lower = size(G);
%upper = size(G);
%for i = 1:size(G, 2)
%    best = -1;
%    for j = linspace(min(G(:, i)), max(G(:, i)), 100)
%        for k = linspace(min(G(:, i)), max(G(:, i)), 100)
%            if (k >= j)
%                break;
%            end
%
%            onespike = ((G(:, i) <= j) & (G(:, i) >= k));
%            pred(onespike, :) = 1;
%            pred(~onespike, :) = 2;
%            correct = (pred == num_real_spikeslist)' * weight;
%            if (correct > best)
%                best = correct;
%                lower(i) = j;
%                upper(i) = k;
%            end
%        end
%    end
%
%    vote(:, i) = ((G(:, i) > lower(i)) & (G(:, i) < upper(i)));
%    fprintf(sprintf('Accuracy %d: %d/%d = %f\n', i, best, total, best/total));
%end
%
%votes = sum(vote, 2);
%pred(votes < num_tests / 2) = 1;
%pred(votes >= num_tests / 2) = 2;
%pred = (ratio3list > 0.09 & ratio4list > 0.24) + 1;
%correct = (pred == num_real_spikeslist)' * weight;
%
%minimum = min(G);
%minimum
%maximum = max(G);
%maximum
%lower
%upper
%fprintf(sprintf('Accuracy: %d/%d = %f\n', correct, total, correct/total));




Gsave = [ratio3list ratio4list (ratio3list.*ratio4list)];
num_tests = size(Gsave, 2);
times = unique(Tlist);
threshold = zeros(num_tests, size(times, 1));
number_of_tests = zeros(1, size(times, 1));
accuracy = zeros(num_tests, size(times, 1));
overall_accuracy = zeros(1, size(times, 1));
for index = 1:size(unique(Tlist), 1)
    clear pred vote;
    time = times(index);
    fprintf('time: %f\n', time);
    G = Gsave((Tlist == time), :);
    total = size(G, 1);
    number_of_tests(index) = total;
    for i = 1:size(G, 2)
        best = -1;
        for j = linspace(min(G(:, i)), max(G(:, i)), 100)
            pred(G(:, i) <= j, :) = 1;
            pred(G(:, i) > j, :) = 2;
            correct = sum(pred == num_real_spikeslist(Tlist == time));
            if (correct > best)
                best = correct;
                threshold(i, index) = j;
            end
        end
    
        vote(:, i) = (G(:, i) > threshold(i, index));
        fprintf(sprintf('Accuracy %d: %d/%d = %f\n', i, best, total, best/total));
        accuracy(i, index) = best/total;
    end
    
    votes = sum(vote, 2);
    pred(votes < num_tests / 2) = 1;
    pred(votes >= num_tests / 2) = 2;
    correct = sum(pred == num_real_spikeslist(Tlist == time));
    
    minimum = min(G);
    minimum
    maximum = max(G);
    maximum
    threshold
    fprintf(sprintf('Accuracy: %d/%d = %f\n', correct, total, correct/total));
    overall_accuracy(index) = correct/total;
end

number_of_tests
threshold
accuracy
overall_accuracy



Gsave = [ratio3list ratio4list];
num_tests = size(Gsave, 2);
times = unique(Tlist);
threshold = zeros(num_tests, size(times, 1));
number_of_tests = zeros(1, size(times, 1));
accuracy = zeros(1, size(times, 1));
overall_accuracy = zeros(1, size(times, 1));
for index = 1:size(unique(Tlist), 1)
    clear pred vote;
    time = times(index);
    fprintf('time: %f\n', time);
    G = Gsave((Tlist == time), :);
    total = size(G, 1);
    number_of_tests(index) = total;
    best = -1;
    for i = linspace(min(G(:, 1)), max(G(:, 1)), 100)
        for j = linspace(min(G(:, 2)), max(G(:, 2)), 100)
	    pred = ones(size(G, 1), 1);
            pred((G(:, 1) > i | G(:, 2) > j), :) = 2;
            correct = sum(pred == num_real_spikeslist(Tlist == time));
            if (correct > best)
                best = correct;
                threshold(1, index) = i;
                threshold(2, index) = j;
            end
        end
    
    end

    fprintf(sprintf('Accuracy %d: %d/%d = %f\n', i, best, total, best/total));
    accuracy(1, index) = best/total;
    
    minimum = min(G);
    minimum
    maximum = max(G);
    maximum
    threshold
end

number_of_tests
threshold
accuracy
